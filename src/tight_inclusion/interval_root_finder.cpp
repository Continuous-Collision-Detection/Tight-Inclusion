// A root finder using interval arithmetic.
#include <tight_inclusion/interval_root_finder.hpp>

#include <stack>
#include <iostream>
#include <queue>
#include <fstream>

#include <tight_inclusion/timer.hpp>
#include <tight_inclusion/avx.hpp>
// #define COMPARE_WITH_RATIONAL

namespace ticcd {
    Scalar time_predicates = 0, time_width = 0, time_bisect = 0,
           time_eval_origin_1D = 0, time_eval_origin_tuv = 0,
           time_vertex_solving = 0;

    template <typename T, int N>
    inline void min_max_array(const std::array<T, N> &arr, T &min, T &max)
    {
        static_assert(N > 0, "no min/max of empty array");
        min = arr[0];
        max = arr[0];
        for (int i = 1; i < N; i++) {
            if (min > arr[i]) {
                min = arr[i];
            }
            if (max < arr[i]) {
                max = arr[i];
            }
        }
    }

    // ** this version can return the true x or y or z tolerance of the co-domain **
    // eps is the interval [-eps,eps] we need to check
    // if [-eps,eps] overlap, return true
    // bbox_in_eps tell us if the box is totally in eps box
    // ms is the minimum seperation
    template <bool check_vf>
    bool evaluate_bbox_one_dimension_vector(
        std::array<Scalar, 8> &t_up,
        std::array<Scalar, 8> &t_dw,
        std::array<Scalar, 8> &u_up,
        std::array<Scalar, 8> &u_dw,
        std::array<Scalar, 8> &v_up,
        std::array<Scalar, 8> &v_dw,
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
        const int dim,
        const Scalar eps,
        bool &bbox_in_eps,
        const Scalar ms = 0,
        Scalar *tol = nullptr)
    {
        TIGHT_INCLUSION_SCOPED_TIMER(time_vertex_solving);

        std::array<Scalar, 8> vs;
        if constexpr (check_vf) {
            vs = function_vf(
                a0s[dim], a1s[dim], b0s[dim], b1s[dim], a0e[dim], a1e[dim],
                b0e[dim], b1e[dim], t_up, t_dw, u_up, u_dw, v_up, v_dw);
        } else {
            vs = function_ee(
                a0s[dim], a1s[dim], b0s[dim], b1s[dim], a0e[dim], a1e[dim],
                b0e[dim], b1e[dim], t_up, t_dw, u_up, u_dw, v_up, v_dw);
        }

        Scalar minv, maxv;
        min_max_array<Scalar, 8>(vs, minv, maxv);

        if (tol != nullptr) {
            *tol = maxv - minv; // this is the real tolerance
        }

        bbox_in_eps = false;

        const Scalar eps_and_ms = eps + ms;

        if (minv > eps_and_ms || maxv < -eps_and_ms) {
            return false;
        }

        if (minv >= -eps_and_ms && maxv <= eps_and_ms) {
            bbox_in_eps = true;
        }

        return true;
    }

    // ** this version can return the true tolerance of the co-domain **
    // give the result of if the hex overlaps the input eps box around origin
    // use vectorized hex-vertex-solving function for acceleration
    // box_in_eps shows if this hex is totally inside box. if so, no need to do further bisection
    template <bool check_vf>
    bool origin_in_function_bounding_box_vector(
        const Interval3 &paras,
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
        const Array3 &eps,
        bool &box_in_eps,
        const Scalar ms = 0,
        Array3 *tolerance = nullptr)
    {
        box_in_eps = false;

        std::array<Scalar, 8> t_up, t_dw, u_up, u_dw, v_up, v_dw;
        {
            TIGHT_INCLUSION_SCOPED_TIMER(time_eval_origin_tuv);
            convert_tuv_to_array(paras, t_up, t_dw, u_up, u_dw, v_up, v_dw);
        }

        bool box_in[3];
        for (int i = 0; i < 3; i++) {
            TIGHT_INCLUSION_SCOPED_TIMER(time_eval_origin_1D);
            Scalar *tol = tolerance == nullptr ? nullptr : &((*tolerance)[i]);
            if (!evaluate_bbox_one_dimension_vector<check_vf>(
                    t_up, t_dw, u_up, u_dw, v_up, v_dw, a0s, a1s, b0s, b1s, a0e,
                    a1e, b0e, b1e, i, eps[i], box_in[i], ms, tol)) {
                return false;
            }
        }

        if (box_in[0] && box_in[1] && box_in[2]) {
            box_in_eps = true;
        }

        return true;
    }

    // find the largest width/tol dimension that is greater than its tolerance
    int find_next_split(const Array3 &widths, const Array3 &tols)
    {
        // assert((widths > tols).any());
        Array3 tmp =
            (widths > tols)
                .select(
                    widths / tols, -std::numeric_limits<Scalar>::infinity());
        int max_index;
        tmp.maxCoeff(&max_index);
        return max_index;
    }

    bool split_and_push(
        const Interval3 &tuv,
        int split_i,
        std::function<void(const Interval3 &)> push,
        bool check_vf,
        Scalar t_upper_bound = 1)
    {
        std::pair<Interval, Interval> halves = tuv[split_i].bisect();
        if (halves.first.lower >= halves.first.upper
            || halves.second.lower >= halves.second.upper) {
            std::cerr << "OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"
                      << std::endl;
            return true;
        }

        Interval3 tmp = tuv;

        if (split_i == 0) {
            if (t_upper_bound == 1
                || halves.second.overlaps(0, t_upper_bound)) {
                tmp[split_i] = halves.second;
                push(tmp);
            }
            if (t_upper_bound == 1 || halves.first.overlaps(0, t_upper_bound)) {
                tmp[split_i] = halves.first;
                push(tmp);
            }
        } else if (!check_vf) { // edge uv
            tmp[split_i] = halves.second;
            push(tmp);
            tmp[split_i] = halves.first;
            push(tmp);
        } else {
            assert(check_vf && split_i != 0);
            // u + v ≤ 1
            if (split_i == 1) {
                const Interval &v = tuv[2];
                if (NumCCD::is_sum_leq_1(halves.second.lower, v.lower)) {
                    tmp[split_i] = halves.second;
                    push(tmp);
                }
                if (NumCCD::is_sum_leq_1(halves.first.lower, v.lower)) {
                    tmp[split_i] = halves.first;
                    push(tmp);
                }
            } else if (split_i == 2) {
                const Interval &u = tuv[1];
                if (NumCCD::is_sum_leq_1(u.lower, halves.second.lower)) {
                    tmp[split_i] = halves.second;
                    push(tmp);
                }
                if (NumCCD::is_sum_leq_1(u.lower, halves.first.lower)) {
                    tmp[split_i] = halves.first;
                    push(tmp);
                }
            }
        }
        return false; // no overflow
    }

    // this version cannot give the impact time at t=1, although this collision can
    // be detected at t=0 of the next time step, but still may cause problems in
    // line-search based physical simulation
    template <bool check_vf>
    bool interval_root_finder_DFS(
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
        const Array3 &tol,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {
        auto cmp_time = [](const Interval3 &i1, const Interval3 &i2) {
            return i1[0].lower >= i2[0].lower;
        };

        //build interval set [0,1]x[0,1]x[0,1]
        const Interval zero_to_one = Interval(NumCCD(0, 0), NumCCD(1, 0));
        Interval3 iset = {{zero_to_one, zero_to_one, zero_to_one}};

        // Stack of intervals and the last split dimension
        // std::stack<std::pair<Interval3,int>> istack;
        std::priority_queue<
            Interval3, std::vector<Interval3>, decltype(cmp_time)>
            istack(cmp_time);
        istack.emplace(iset);

        Array3 err_and_ms = err + ms;

        int refine = 0;

        toi = std::numeric_limits<Scalar>::infinity();
        NumCCD TOI(1, 0);

        bool collision = false;
        int rnbr = 0;
        while (!istack.empty()) {
            Interval3 current = istack.top();
            istack.pop();

            // if(rnbr>0&&less_than( current[0].first,TOI)){
            //     std::cout<<"not the first"<<std::endl;
            //     // continue;
            // }

            //TOI should always be no larger than current
            if (current[0].lower >= TOI) {
                continue;
            }

            refine++;

            bool zero_in, box_in;
            {
                TIGHT_INCLUSION_SCOPED_TIMER(time_predicates);
                zero_in = origin_in_function_bounding_box_vector<check_vf>(
                    current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, err_and_ms,
                    box_in);
            }

            // #ifdef TIGHT_INCLUSION_USE_GMP // this is defined in the begining of this file
            // zero_in = origin_in_function_bounding_box_rational<check_vf>(
            //     current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
            // #endif

            if (!zero_in) {
                continue;
            }

            Array3 widths;
            {
                TIGHT_INCLUSION_SCOPED_TIMER(time_width);
                widths = width(current);
            }

            if (box_in || (widths <= tol).all()) {
                TOI = current[0].lower;
                collision = true;
                rnbr++;
                // continue;
                toi = TOI.value();
                return true;
            }

            // find the next dimension to split
            int split_i = find_next_split(widths, tol);

            bool overflowed = split_and_push(
                current, split_i,
                [&](const Interval3 &i) { istack.emplace(i); }, check_vf);
            if (overflowed) {
                std::cerr << "OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"
                          << std::endl;
                return true;
            }
        }
        if (collision)
            toi = TOI.value();
        return collision;
    }

    // when check_t_overlap = false, check [0,1]x[0,1]x[0,1]; otherwise, check [0, t_max]x[0,1]x[0,1]
    template <bool check_vf>
    bool interval_root_finder_BFS(
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
        const Interval3 &iset,
        const Array3 &tol,
        const Scalar co_domain_tolerance,
        const Array3 &err,
        const Scalar ms,
        const Scalar max_time,
        const int max_itr,
        Scalar &toi,
        Scalar &output_tolerance)
    {
        long queue_size = 0;
        // if max_itr <0, output_tolerance= co_domain_tolerance;
        // else, output_tolearancewill be the precision after iteration time > max_itr
        output_tolerance = co_domain_tolerance;

        // this is used to catch the tolerance for each level
        Scalar temp_output_tolerance = co_domain_tolerance;
        // return time1 >= time2
        // auto time_cmp = [](const std::pair<Interval3, int> &i1,
        //                    const std::pair<Interval3, int> &i2) {
        //     return i1.first[0].lower >= i2.first[0].lower;
        // };

        // check the tree level by level instead of going deep
        // (if level 1 != level 2, return level 1 >= level 2; else, return time1 >= time2)
        auto horiz_cmp = [](const std::pair<Interval3, int> &i1,
                            const std::pair<Interval3, int> &i2) {
            if (i1.second != i2.second) {
                return i1.second >= i2.second;
            } else {
                return i1.first[0].lower > i2.first[0].lower;
            }
        };

        // Stack of intervals and the last split dimension
        // std::stack<std::pair<Interval3,int>> istack;
        const auto &cmp = horiz_cmp;
        std::priority_queue<
            std::pair<Interval3, int>, std::vector<std::pair<Interval3, int>>,
            decltype(cmp)>
            istack(cmp);
        istack.emplace(iset, -1);

        // current intervals
        Interval3 current;
        int refine = 0;
        Scalar impact_ratio = 1;

        toi = std::numeric_limits<Scalar>::infinity(); //set toi as infinate
        // temp_toi is to catch the toi of each level
        Scalar temp_toi = toi;
        // set TOI to 4. this is to record the impact time of this level
        NumCCD TOI(4, 0);
        // this is to record the element that already small enough or contained in eps-box
        NumCCD TOI_SKIP = TOI;
        bool use_skip = false; // this is to record if TOI_SKIP is used.
        int rnbr = 0;
        int current_level = -2; // in the begining, current_level != level
        int box_in_level = -2;  // this checks if all the boxes before this
        // level < tolerance. only true, we can return when we find one overlaps eps box and smaller than tolerance or eps-box
        bool this_level_less_tol = true;
        bool find_level_root = false;
        Scalar t_upper_bound = max_time; // 2*tol make it more conservative
        while (!istack.empty()) {
#ifdef CHECK_QUEUE_SIZE
            if (istack.size() > queue_size) {
                queue_size = istack.size();
            }
#endif
#ifdef TIGHT_INCLUSION_LIMIT_QUEUE_SIZE
            if (istack.size() > MAX_QSIZE) {
                return true;
            }
#endif

            current = istack.top().first;
            int level = istack.top().second;
            istack.pop();

            // if this box is later than TOI_SKIP in time, we can skip this one.
            // TOI_SKIP is only updated when the box is small enough or totally contained in eps-box
            if (current[0].lower >= TOI_SKIP) {
                continue;
            }
            // before check a new level, set this_level_less_tol=true
            if (box_in_level != level) {
                box_in_level = level;
                this_level_less_tol = true;
            }

            refine++;
            bool zero_in, box_in;
            Array3 true_tol;
            {
                TIGHT_INCLUSION_SCOPED_TIMER(time_predicates);
                // #ifdef TIGHT_INCLUSION_USE_GMP // this is defined in the begining of this file
                // Array3 ms_3d = Array3::Constant(ms);
                // zero_in = origin_in_function_bounding_box_rational_return_tolerance<
                //     check_vf>(
                //     current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, ms_3d, box_in,
                //     true_tol);
                // #else
                zero_in = origin_in_function_bounding_box_vector<check_vf>(
                    current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, err,
                    box_in, ms, &true_tol);
                // #endif
            }

            if (!zero_in)
                continue;

            Array3 widths;
            {
                TIGHT_INCLUSION_SCOPED_TIMER(time_width);
                widths = width(current);
            }

            bool tol_condition = (true_tol <= co_domain_tolerance).all();

            // Condition 1, stopping condition on t, u and v is satisfied. this is useless now since we have condition 2
            bool condition1 = (widths <= tol).all();

            // Condition 2, zero_in = true, box inside eps-box and in this level,
            // no box whose zero_in is true but box size larger than tolerance, can return
            bool condition2 = box_in && this_level_less_tol;
            if (!tol_condition) {
                this_level_less_tol = false;
                // this level has at least one box whose size > tolerance, thus we
                // cannot directly return if find one box whose size < tolerance or box-in
                // TODO: think about it. maybe we can return even if this value is false, so we can terminate earlier.
            }

            // Condition 3, in this level, we find a box that zero-in and size < tolerance.
            // and no other boxes whose zero-in is true in this level before this one is larger than tolerance, can return
            bool condition3 = this_level_less_tol;
            if (condition1 || condition2 || condition3) {
                TOI = current[0].lower;
                rnbr++;
                // continue;
                toi = TOI.value() * impact_ratio;
                // we don't need to compare with TOI_SKIP because we already
                // continue when t >= TOI_SKIP
                return true;
            }

            if (max_itr > 0) { // if max_itr <= 0 ⟹ unlimited iterations
                if (current_level != level) {
                    // output_tolerance=current_tolerance;
                    // current_tolerance=0;
                    current_level = level;
                    find_level_root = false;
                }
                // current_tolerance=std::max(
                // std::max(std::max(current_tolerance,true_tol[0]),true_tol[1]),true_tol[2]
                // );
                if (!find_level_root) {
                    TOI = current[0].lower;
                    // collision=true;
                    rnbr++;
                    // continue;
                    temp_toi = TOI.value() * impact_ratio;

                    // if the real tolerance is larger than input, use the real one;
                    // if the real tolerance is smaller than input, use input
                    temp_output_tolerance = std::max(
                        {true_tol[0], true_tol[1], true_tol[2],
                         co_domain_tolerance});
                    // this ensures always find the earlist root
                    find_level_root = true;
                }
                if (refine > max_itr) {
                    toi = temp_toi;
                    output_tolerance = temp_output_tolerance;

                    // std::cout<<"return from refine"<<std::endl;
                    return true;
                }
                // get the time of impact down here
            }

            // if this box is small enough, or inside of eps-box, then just continue,
            // but we need to record the collision time
            if (tol_condition || box_in) {
                if (current[0].lower < TOI_SKIP) {
                    TOI_SKIP = current[0].lower;
                }
                use_skip = true;
                continue;
            }

            // find the next dimension to split
            int split_i = find_next_split(widths, tol);

            bool overflow = split_and_push(
                current, split_i,
                [&](const Interval3 &i) { istack.emplace(i, level + 1); },
                check_vf, t_upper_bound);
            if (overflow) {
                std::cout << "OVERFLOW HAPPENS WHEN SPLITTING INTERVALS"
                          << std::endl;
                return true;
            }
        }

        if (use_skip) {
            toi = TOI_SKIP.value() * impact_ratio;
            return true;
        }

        return false;
    }

    template <bool check_vf>
    bool interval_root_finder_BFS(
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
        const Array3 &tol,
        const Scalar co_domain_tolerance,
        const Array3 &err,
        const Scalar ms,
        const Scalar max_time,
        const int max_itr,
        Scalar &toi,
        Scalar &output_tolerance)
    {
        // build interval set [0,t_max]x[0,1]x[0,1]
        const Interval zero_to_one = Interval(NumCCD(0, 0), NumCCD(1, 0));
        Interval3 iset = {{
            // Interval(NumCCD(0, 0), NumCCD(max_time)),
            zero_to_one,
            zero_to_one,
            zero_to_one,
        }};

        return interval_root_finder_BFS<check_vf>(
            a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, iset, tol,
            co_domain_tolerance, err, ms, max_time, max_itr, toi,
            output_tolerance);
    }

    void print_times()
    {
        // clang-format off
        std::cout << "origin predicates, " << time_predicates << "\n"
                  << "width, " << time_width << "\n"
                  << "bisect, " << time_bisect << "\n"
                  << "origin part1(evaluate 1 dimension), " << time_eval_origin_1D << "\n"
                  << "origin part2(convert tuv), " << time_eval_origin_tuv << "\n"
                  << "time of call the vertex solving function, " << time_vertex_solving << std::endl;
        // clang-format on
    }

    Array3 get_numerical_error(
        const std::vector<Vector3> &vertices,
        const bool check_vf,
        const bool using_minimum_separation)
    {
        Scalar eefilter;
        Scalar vffilter;
        if (!using_minimum_separation) {
#ifdef TIGHT_INCLUSION_DOUBLE
            eefilter = 6.217248937900877e-15;
            vffilter = 6.661338147750939e-15;
#else
            eefilter = 3.337861e-06;
            vffilter = 3.576279e-06;
#endif
        } else // using minimum separation
        {
#ifdef TIGHT_INCLUSION_DOUBLE
            eefilter = 7.105427357601002e-15;
            vffilter = 7.549516567451064e-15;
#else
            eefilter = 3.814698e-06;
            vffilter = 4.053116e-06;
#endif
        }

        Vector3 max = vertices[0].cwiseAbs();
        for (int i = 0; i < vertices.size(); i++) {
            max = max.cwiseMax(vertices[i].cwiseAbs());
        }
        Vector3 delta = max.cwiseMin(1);
        Scalar filter = check_vf ? vffilter : eefilter;
        return filter * delta.array().pow(3);
    }

    //////////////////////////////////////////////////////////////////////////
    // Template instantiation
    //////////////////////////////////////////////////////////////////////////
    bool edge_edge_interval_root_finder_DFS(
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
        const Array3 &tol,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {
        return interval_root_finder_DFS<false>(
            a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, tol, err, ms, toi);
    }

    bool vertex_face_interval_root_finder_DFS(
        const Vector3 &vertex_start,
        const Vector3 &face_vertex0_start,
        const Vector3 &face_vertex1_start,
        const Vector3 &face_vertex2_start,
        const Vector3 &vertex_end,
        const Vector3 &face_vertex0_end,
        const Vector3 &face_vertex1_end,
        const Vector3 &face_vertex2_end,
        const Array3 &tol,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {
        return interval_root_finder_DFS<true>(
            vertex_start, face_vertex0_start, face_vertex1_start,
            face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
            face_vertex2_end, tol, err, ms, toi);
    }

    bool edge_edge_interval_root_finder_BFS(
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
        const Array3 &tol,
        const Scalar co_domain_tolerance,
        // this is the maximum error on each axis when calculating the vertices, err, aka, filter
        const Array3 &err,
        const Scalar ms,
        const Scalar max_time,
        const int max_itr,
        Scalar &toi,
        Scalar &output_tolerance)
    {
        return interval_root_finder_BFS<false>(
            a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, tol, co_domain_tolerance,
            err, ms, max_time, max_itr, toi, output_tolerance);
    }

    bool vertex_face_interval_root_finder_BFS(
        const Vector3 &vertex_start,
        const Vector3 &face_vertex0_start,
        const Vector3 &face_vertex1_start,
        const Vector3 &face_vertex2_start,
        const Vector3 &vertex_end,
        const Vector3 &face_vertex0_end,
        const Vector3 &face_vertex1_end,
        const Vector3 &face_vertex2_end,
        const Array3 &tol,
        const Scalar co_domain_tolerance,
        // this is the maximum error on each axis when calculating the vertices, err, aka, filter
        const Array3 &err,
        const Scalar ms,
        const Scalar max_time,
        const int max_itr,
        Scalar &toi,
        Scalar &output_tolerance)
    {
        return interval_root_finder_BFS<true>(
            vertex_start, face_vertex0_start, face_vertex1_start,
            face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
            face_vertex2_end, tol, co_domain_tolerance, err, ms, max_time,
            max_itr, toi, output_tolerance);
    }

} // namespace ticcd
