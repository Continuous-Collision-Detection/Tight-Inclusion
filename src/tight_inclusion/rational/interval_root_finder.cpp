// A root finder using interval arithmetic.
#include <tight_inclusion/rational/interval_root_finder.hpp>

#include <stack>
#include <iostream>

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
        const NumCCD &tpara,
        const NumCCD &upara,
        const NumCCD &vpara,
        const Vector3 &a0sd,
        const Vector3 &a1sd,
        const Vector3 &b0sd,
        const Vector3 &b1sd,
        const Vector3 &a0ed,
        const Vector3 &a1ed,
        const Vector3 &b0ed,
        const Vector3 &b1ed)
    {
        Rational t = Rational(double(tpara.numerator))
                     / Rational(double(tpara.denominator()));
        Rational u = Rational(double(upara.numerator))
                     / Rational(double(upara.denominator()));
        Rational v = Rational(double(vpara.numerator))
                     / Rational(double(vpara.denominator()));

        Vector3r a0s = a0sd.cast<Rational>(), a1s = a1sd.cast<Rational>(),
                 b0s = b0sd.cast<Rational>(), b1s = b1sd.cast<Rational>(),
                 a0e = a0ed.cast<Rational>(), a1e = a1ed.cast<Rational>(),
                 b0e = b0ed.cast<Rational>(), b1e = b1ed.cast<Rational>();

        Vector3r edge0_vertex0 = (a0e - a0s) * t + a0s;
        Vector3r edge0_vertex1 = (a1e - a1s) * t + a1s;
        Vector3r edge0_vertex =
            (edge0_vertex1 - edge0_vertex0) * u + edge0_vertex0;

        Vector3r edge1_vertex0 = (b0e - b0s) * t + b0s;
        Vector3r edge1_vertex1 = (b1e - b1s) * t + b1s;
        Vector3r edge1_vertex =
            (edge1_vertex1 - edge1_vertex0) * v + edge1_vertex0;

        return edge1_vertex - edge0_vertex;
    }

    Vector3r function_f_ee_rational(
        const Rational &tpara,
        const Rational &upara,
        const Rational &vpara,
        const Vector3 &a0sd,
        const Vector3 &a1sd,
        const Vector3 &b0sd,
        const Vector3 &b1sd,
        const Vector3 &a0ed,
        const Vector3 &a1ed,
        const Vector3 &b0ed,
        const Vector3 &b1ed)
    {
        Vector3r a0s = a0sd.cast<Rational>(), a1s = a1sd.cast<Rational>(),
                 b0s = b0sd.cast<Rational>(), b1s = b1sd.cast<Rational>(),
                 a0e = a0ed.cast<Rational>(), a1e = a1ed.cast<Rational>(),
                 b0e = b0ed.cast<Rational>(), b1e = b1ed.cast<Rational>();

        Vector3r edge0_vertex0 = (a0e - a0s) * tpara + a0s;
        Vector3r edge0_vertex1 = (a1e - a1s) * tpara + a1s;
        Vector3r edge0_vertex =
            (edge0_vertex1 - edge0_vertex0) * upara + edge0_vertex0;

        Vector3r edge1_vertex0 = (b0e - b0s) * tpara + b0s;
        Vector3r edge1_vertex1 = (b1e - b1s) * tpara + b1s;
        Vector3r edge1_vertex =
            (edge1_vertex1 - edge1_vertex0) * vpara + edge1_vertex0;

        return edge1_vertex - edge0_vertex;
    }

    Vector3r function_f_vf_rational(
        const NumCCD &tpara,
        const NumCCD &upara,
        const NumCCD &vpara,
        const Vector3 &a0sd,
        const Vector3 &a1sd,
        const Vector3 &b0sd,
        const Vector3 &b1sd,
        const Vector3 &a0ed,
        const Vector3 &a1ed,
        const Vector3 &b0ed,
        const Vector3 &b1ed)
    {
        Rational t = Rational(double(tpara.numerator))
                     / Rational(double(tpara.denominator()));
        Rational u = Rational(double(upara.numerator))
                     / Rational(double(upara.denominator()));
        Rational v = Rational(double(vpara.numerator))
                     / Rational(double(vpara.denominator()));

        Vector3r vs = a0sd.cast<Rational>(), t0s = a1sd.cast<Rational>(),
                 t1s = b0sd.cast<Rational>(), t2s = b1sd.cast<Rational>(),
                 ve = a0ed.cast<Rational>(), t0e = a1ed.cast<Rational>(),
                 t1e = b0ed.cast<Rational>(), t2e = b1ed.cast<Rational>();

        Vector3r vert = (ve - vs) * t + vs;

        Vector3r t0 = (t0e - t0s) * t + t0s;
        Vector3r t1 = (t1e - t1s) * t + t1s;
        Vector3r t2 = (t2e - t2s) * t + t2s;
        Vector3r p = (t1 - t0) * u + (t2 - t0) * v + t0;
        return vert - p;
    }

    Vector3r function_f_vf_rational(
        const Rational &tpara,
        const Rational &upara,
        const Rational &vpara,
        const Vector3 &a0sd,
        const Vector3 &a1sd,
        const Vector3 &b0sd,
        const Vector3 &b1sd,
        const Vector3 &a0ed,
        const Vector3 &a1ed,
        const Vector3 &b0ed,
        const Vector3 &b1ed)
    {
        Vector3r vs = a0sd.cast<Rational>(), t0s = a1sd.cast<Rational>(),
                 t1s = b0sd.cast<Rational>(), t2s = b1sd.cast<Rational>(),
                 ve = a0ed.cast<Rational>(), t0e = a1ed.cast<Rational>(),
                 t1e = b0ed.cast<Rational>(), t2e = b1ed.cast<Rational>();

        Vector3r v = (ve - vs) * tpara + vs;

        Vector3r t0 = (t0e - t0s) * tpara + t0s;
        Vector3r t1 = (t1e - t1s) * tpara + t1s;
        Vector3r t2 = (t2e - t2s) * tpara + t2s;
        Vector3r p = (t1 - t0) * upara + (t2 - t0) * vpara + t0;
        return v - p;
    }

    template <typename T, bool check_vf>
    bool origin_in_function_bounding_box_rational(
        const std::array<std::array<T, 2>, 3> &paras,
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
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
                    if constexpr (!check_vf) {
                        pts.row(c) = function_f_ee_rational(
                            t[i], u[j], v[k], a0s, a1s, b0s, b1s, a0e, a1e, b0e,
                            b1e);
                    } else {
                        pts.row(c) = function_f_vf_rational(
                            t[i], u[j], v[k], a0s, a1s, b0s, b1s, a0e, a1e, b0e,
                            b1e);
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

    template <bool check_vf>
    bool interval_root_finder(
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
            bool zero_in =
                origin_in_function_bounding_box_rational<Rational, check_vf>(
                    current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e,
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

            if (check_vf) {
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
        std::array<RationalInterval, 3> &root)
    {
        return interval_root_finder<false>(
            a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, tol, err, ms, root);
    }

    bool vertex_face_interval_root_finder(
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
        std::array<RationalInterval, 3> &root)
    {
        return interval_root_finder<true>(
            vertex_start, face_vertex0_start, face_vertex1_start,
            face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
            face_vertex2_end, tol, err, ms, root);
    }
} // namespace ticcd::rational
