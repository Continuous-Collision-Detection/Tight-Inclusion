// A root finder using interval arithmetic.
#include <tight_inclusion/rational/interval_root_finder.hpp>

#include <stack>

// #define COMPARE_WITH_RATIONAL

namespace ticcd::rational {

    Array3r width(const std::array<std::array<Rational, 2>, 3> &x)
    {
        Array3r w;
        for (int i = 0; i < 3; i++) {
            Rational sub = x[i][1] - x[i][0];
            w[i] = sub >= 0 ? sub : -sub;
            assert(w[i] >= 0);
        }
        return w;
    }

    Vector3r function_f_ee_rational(
        const Rational &t,
        const Rational &u,
        const Rational &v,
        const Vector3 &ea0_t0d,
        const Vector3 &ea1_t0d,
        const Vector3 &eb0_t0d,
        const Vector3 &eb1_t0d,
        const Vector3 &ea0_t1d,
        const Vector3 &ea1_t1d,
        const Vector3 &eb0_t1d,
        const Vector3 &eb1_t1d)
    {
        const Vector3r ea0_t0 = ea0_t0d.cast<Rational>();
        const Vector3r ea1_t0 = ea1_t0d.cast<Rational>();
        const Vector3r eb0_t0 = eb0_t0d.cast<Rational>();
        const Vector3r eb1_t0 = eb1_t0d.cast<Rational>();
        const Vector3r ea0_t1 = ea0_t1d.cast<Rational>();
        const Vector3r ea1_t1 = ea1_t1d.cast<Rational>();
        const Vector3r eb0_t1 = eb0_t1d.cast<Rational>();
        const Vector3r eb1_t1 = eb1_t1d.cast<Rational>();

        const Vector3r ea0 = (ea0_t1 - ea0_t0) * t + ea0_t0;
        const Vector3r ea1 = (ea1_t1 - ea1_t0) * t + ea1_t0;
        const Vector3r va = (ea1 - ea0) * u + ea0;

        const Vector3r eb0 = (eb0_t1 - eb0_t0) * t + eb0_t0;
        const Vector3r eb1 = (eb1_t1 - eb1_t0) * t + eb1_t0;
        const Vector3r vb = (eb1 - eb0) * v + eb0;

        return vb - va;
    }

    Vector3r function_f_ee_rational(
        const NumCCD &tpara,
        const NumCCD &upara,
        const NumCCD &vpara,
        const Vector3 &ea0_t0d,
        const Vector3 &ea1_t0d,
        const Vector3 &eb0_t0d,
        const Vector3 &eb1_t0d,
        const Vector3 &ea0_t1d,
        const Vector3 &ea1_t1d,
        const Vector3 &eb0_t1d,
        const Vector3 &eb1_t1d)
    {
        return function_f_ee_rational(
            Rational(double(tpara.numerator))
                / Rational(double(tpara.denominator())),
            Rational(double(upara.numerator))
                / Rational(double(upara.denominator())),
            Rational(double(vpara.numerator))
                / Rational(double(vpara.denominator())),
            ea0_t0d, ea1_t0d, eb0_t0d, eb1_t0d, ea0_t1d, ea1_t1d, eb0_t1d,
            eb1_t1d);
    }

    Vector3r function_f_vf_rational(
        const Rational &t,
        const Rational &u,
        const Rational &v,
        const Vector3 &v_t0d,
        const Vector3 &f0_t0d,
        const Vector3 &f1_t0d,
        const Vector3 &f2_t0d,
        const Vector3 &v_t1d,
        const Vector3 &f0_t1d,
        const Vector3 &f1_t1d,
        const Vector3 &f2_t1d)
    {
        const Vector3r v_t0 = v_t0d.cast<Rational>();
        const Vector3r f0_t0 = f0_t0d.cast<Rational>();
        const Vector3r f1_t0 = f1_t0d.cast<Rational>();
        const Vector3r f2_t0 = f2_t0d.cast<Rational>();
        const Vector3r v_t1 = v_t1d.cast<Rational>();
        const Vector3r f0_t1 = f0_t1d.cast<Rational>();
        const Vector3r f1_t1 = f1_t1d.cast<Rational>();
        const Vector3r f2_t1 = f2_t1d.cast<Rational>();

        const Vector3r va = (v_t1 - v_t0) * t + v_t0;

        const Vector3r f0 = (f0_t1 - f0_t0) * t + f0_t0;
        const Vector3r f1 = (f1_t1 - f1_t0) * t + f1_t0;
        const Vector3r f2 = (f2_t1 - f2_t0) * t + f2_t0;
        const Vector3r vb = (f1 - f0) * u + (f2 - f0) * v + f0;
        return va - vb;
    }

    Vector3r function_f_vf_rational(
        const NumCCD &tpara,
        const NumCCD &upara,
        const NumCCD &vpara,
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1)
    {
        return function_f_ee_rational(
            Rational(double(tpara.numerator))
                / Rational(double(tpara.denominator())),
            Rational(double(upara.numerator))
                / Rational(double(upara.denominator())),
            Rational(double(vpara.numerator))
                / Rational(double(vpara.denominator())),
            v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1);
    }

    template <typename T, bool is_vertex_face>
    bool origin_in_function_bounding_box_rational(
        const std::array<std::array<T, 2>, 3> &paras,
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Array3 &box,
        bool &box_in_eps,
        Array3 *tolerance = nullptr)
    {
        std::array<T, 2> t, u, v;
        t = paras[0];
        u = paras[1];
        v = paras[2];

        Eigen::Matrix<T, 8, 3, Eigen::ColMajor | Eigen::DontAlign> pts;
        int c = 0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    if constexpr (!is_vertex_face) {
                        pts.row(c) = function_f_ee_rational(
                            t[i], u[j], v[k], a_t0, b_t0, c_t0, d_t0, a_t1,
                            b_t1, c_t1, d_t1);
                    } else {
                        pts.row(c) = function_f_vf_rational(
                            t[i], u[j], v[k], a_t0, b_t0, c_t0, d_t0, a_t1,
                            b_t1, c_t1, d_t1);
                    }
                    c++;
                }
            }
        }

        Eigen::Array<T, 3, 1> minv = pts.colwise().minCoeff();
        Eigen::Array<T, 3, 1> maxv = pts.colwise().maxCoeff();

        if (tolerance != nullptr) {
            *tolerance = (maxv - minv).template cast<Scalar>();
        }
        auto boxR = box.cast<Rational>();
        box_in_eps = (minv >= -boxR).all() && (maxv <= boxR).all();
        return (minv <= boxR).all() && (maxv >= -boxR).all();
    }

    std::pair<RationalInterval, RationalInterval>
    bisect(const RationalInterval &inter)
    {
        Rational mid = (inter[0] + inter[1]) * Rational(0.5);
        return std::make_pair<RationalInterval, RationalInterval>(
            {{inter[0], mid}}, {{mid, inter[1]}});
    }

    template <bool is_vertex_face>
    bool interval_root_finder(
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Array3 &tol,
        const Array3 &err,
        const Scalar ms,
        std::array<RationalInterval, 3> &root)
    {
        RationalInterval interval01 = {{Rational(0), Rational(1)}};
        std::array<RationalInterval, 3> paracube;
        paracube[0] = interval01;
        paracube[1] = interval01;
        paracube[2] = interval01;

        // Stack of intervals and the last split dimension
        std::stack<std::pair<std::array<RationalInterval, 3>, int>> istack;
        istack.emplace(paracube, -1);

        // current intervals
        std::array<RationalInterval, 3> current;
        Array3 err_and_ms;
        err_and_ms[0] = err[0] + ms;
        err_and_ms[1] = err[1] + ms;
        err_and_ms[2] = err[2] + ms;

        Rational max_rational = Rational(std::numeric_limits<Scalar>::max());
        const Eigen::Array<Rational, 3, 1> tolR(
            std::isfinite(tol[0]) ? Rational(tol[0]) : max_rational,
            std::isfinite(tol[1]) ? Rational(tol[1]) : max_rational,
            std::isfinite(tol[2]) ? Rational(tol[2]) : max_rational);

        while (!istack.empty()) {
            current = istack.top().first;
            int last_split = istack.top().second;
            istack.pop();

            bool box_in_eps;
            bool zero_in = origin_in_function_bounding_box_rational<
                Rational, is_vertex_face>(
                current, a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                /*box=*/Array3::Zero(), box_in_eps);

            if (!zero_in)
                continue;
            Array3r widths = width(current);
            if ((widths <= tolR.array()).all()) {
                root = current;
                return true;
            }

            // Bisect the next dimension that is greater than its tolerance
            int split_i;
            for (int i = 1; i <= 3; i++) {
                split_i = (last_split + i) % 3;
                if (widths[split_i] > tolR(split_i)) {
                    break;
                }
            }
            std::pair<RationalInterval, RationalInterval> halves =
                bisect(current[split_i]);

            if (is_vertex_face) {
                if (split_i == 1) {
                    if (halves.first[0] + current[2][0] <= 1) {
                        current[split_i] = halves.first;
                        istack.emplace(current, split_i);
                    }
                    if (halves.second[0] + current[2][0] <= 1) {
                        current[split_i] = halves.second;
                        istack.emplace(current, split_i);
                    }
                }

                if (split_i == 2) {

                    if (halves.first[0] + current[1][0] <= 1) {
                        current[split_i] = halves.first;
                        istack.emplace(current, split_i);
                    }
                    if (halves.second[0] + current[1][0] <= 1) {
                        current[split_i] = halves.second;
                        istack.emplace(current, split_i);
                    }
                }
                if (split_i == 0) {
                    current[split_i] = halves.second;
                    istack.emplace(current, split_i);
                    current[split_i] = halves.first;
                    istack.emplace(current, split_i);
                }
            } else {
                current[split_i] = halves.second;
                istack.emplace(current, split_i);
                current[split_i] = halves.first;
                istack.emplace(current, split_i);
            }
        }
        return false;
    }

    bool edge_edge_interval_root_finder(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        const Array3 &tol,
        const Array3 &err,
        const Scalar ms,
        std::array<RationalInterval, 3> &root)
    {
        return interval_root_finder<false>(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, tol,
            err, ms, root);
    }

    bool vertex_face_interval_root_finder(
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1,
        const Array3 &tol,
        const Array3 &err,
        const Scalar ms,
        std::array<RationalInterval, 3> &root)
    {
        return interval_root_finder<true>(
            v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, tol, err, ms,
            root);
    }
} // namespace ticcd::rational
