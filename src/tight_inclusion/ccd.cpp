#include "ccd.hpp"

#include <tight_inclusion/interval_root_finder.hpp>
#include <tight_inclusion/timer.hpp>

#include <vector>

namespace ticcd {

#ifdef TIGHT_INCLUSION_USE_MAX_ABS_TOL
    static constexpr Scalar CCD_MAX_TIME_TOL = 1e-3;
    static constexpr Scalar CCD_MAX_COORD_TOL = 1e-2;
#else
    static constexpr Scalar CCD_MAX_TIME_TOL =
        std::numeric_limits<double>::infinity();
    static constexpr Scalar CCD_MAX_COORD_TOL =
        std::numeric_limits<double>::infinity();
#endif

    inline std::array<Vector3, 2> bbd_4_pts(
        const Vector3 &p0,
        const Vector3 &p1,
        const Vector3 &p2,
        const Vector3 &p3)
    {
        return {
            {p0.cwiseMin(p1).cwiseMin(p2).cwiseMin(p3),
             p0.cwiseMax(p1).cwiseMax(p2).cwiseMax(p3)}};
    }

    // calculate maximum x, y and z diff
    Scalar get_max_axis_diff(
        const std::array<Vector3, 2> &b1, const std::array<Vector3, 2> &b2)
    {
        return std::max({
            (b1[1] - b1[0]).maxCoeff(),
            (b2[1] - b2[0]).maxCoeff(),
            (b2[0] - b1[1]).cwiseAbs().maxCoeff(),
            (b1[0] - b2[1]).cwiseAbs().maxCoeff(),
        });
    }

    inline Scalar max_linf_4(
        const Vector3 &p1,
        const Vector3 &p2,
        const Vector3 &p3,
        const Vector3 &p4,
        const Vector3 &p1e,
        const Vector3 &p2e,
        const Vector3 &p3e,
        const Vector3 &p4e)
    {
        return std::max(
            {(p1e - p1).lpNorm<Eigen::Infinity>(),
             (p2e - p2).lpNorm<Eigen::Infinity>(),
             (p3e - p3).lpNorm<Eigen::Infinity>(),
             (p4e - p4).lpNorm<Eigen::Infinity>()});
    }

    /// @brief Clamp a/b to [-∞, max_val]
    /// @param a numerator
    /// @param b denominator
    /// @param max_val
    /// @return a/b if b != 0, max_val if b == 0
    inline Scalar
    clamp_div(const Scalar a, const Scalar b, const Scalar max_val)
    {
        if (b == 0) {
            return max_val;
        } else {
            return std::min(a / b, max_val);
        }
    }

    Array3 compute_vertex_face_tolerances(
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1,
        const Scalar distance_tolerance)
    {
        const Vector3 p000 = v_t0 - f0_t0;
        const Vector3 p001 = v_t0 - f2_t0;
        const Vector3 p011 = v_t0 - (f1_t0 + f2_t0 - f0_t0);
        const Vector3 p010 = v_t0 - f1_t0;
        const Vector3 p100 = v_t1 - f0_t1;
        const Vector3 p101 = v_t1 - f2_t1;
        const Vector3 p111 = v_t1 - (f1_t1 + f2_t1 - f0_t1);
        const Vector3 p110 = v_t1 - f1_t1;

        const Scalar dl =
            3 * max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110);
        const Scalar edge0_length =
            3 * max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011);
        const Scalar edge1_length =
            3 * max_linf_4(p000, p100, p110, p010, p001, p101, p111, p011);

        return Array3(
            clamp_div(distance_tolerance, dl, CCD_MAX_TIME_TOL),
            clamp_div(distance_tolerance, edge0_length, CCD_MAX_COORD_TOL),
            clamp_div(distance_tolerance, edge1_length, CCD_MAX_COORD_TOL));
    }

    Array3 compute_edge_edge_tolerances(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        const Scalar distance_tolerance)
    {

        const Vector3 p000 = ea0_t0 - eb0_t0;
        const Vector3 p001 = ea0_t0 - eb1_t0;
        const Vector3 p010 = ea1_t0 - eb0_t0;
        const Vector3 p011 = ea1_t0 - eb1_t0;
        const Vector3 p100 = ea0_t1 - eb0_t1;
        const Vector3 p101 = ea0_t1 - eb1_t1;
        const Vector3 p110 = ea1_t1 - eb0_t1;
        const Vector3 p111 = ea1_t1 - eb1_t1;

        const Scalar dl =
            3 * max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110);
        const Scalar edge0_length =
            3 * max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011);
        const Scalar edge1_length =
            3 * max_linf_4(p000, p100, p110, p010, p001, p101, p111, p011);

        return Array3(
            clamp_div(distance_tolerance, dl, CCD_MAX_TIME_TOL),
            clamp_div(distance_tolerance, edge0_length, CCD_MAX_COORD_TOL),
            clamp_div(distance_tolerance, edge1_length, CCD_MAX_COORD_TOL));
    }

    template <bool is_vertex_face>
    bool
    CCD(const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Array3 &err_in,
        const Scalar ms_in,
        Scalar &toi,
        const Scalar tolerance_in,
        const Scalar t_max_in,
        const long max_itr,
        Scalar &output_tolerance,
        bool no_zero_toi,
        const CCDRootFindingMethod ccd_method)
    {
        constexpr int MAX_NO_ZERO_TOI_ITER = std::numeric_limits<int>::max();
        // unsigned so can be larger than MAX_NO_ZERO_TOI_ITER
        unsigned int no_zero_toi_iter = 0;

        bool is_impacting, tmp_is_impacting;

        // Mutable copies for no_zero_toi
        Scalar t_max = t_max_in;
        Scalar tolerance = tolerance_in;
        Scalar ms = ms_in;

        Array3 tol;
        if constexpr (is_vertex_face) {
            tol = compute_vertex_face_tolerances(
                a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1, tolerance_in);
        } else {
            tol = compute_edge_edge_tolerances(
                a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1, tolerance_in);
        }

        //////////////////////////////////////////////////////////
        // this should be the error of the whole mesh
        Array3 err;
        // if error[0] < 0, means we need to calculate error here
        if (err_in[0] < 0) {
            err = get_numerical_error(
                std::vector<Vector3>{
                    {a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1}},
                is_vertex_face, ms > 0);
        } else {
            err = err_in;
        }
        //////////////////////////////////////////////////////////

        do {
            switch (ccd_method) {
            case CCDRootFindingMethod::DEPTH_FIRST_SEARCH:
                // no handling for zero toi
                return interval_root_finder_DFS<is_vertex_face>(
                    a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1, tol, err,
                    ms, toi);
            case CCDRootFindingMethod::BREADTH_FIRST_SEARCH:
                assert(t_max >= 0 && t_max <= 1);
                tmp_is_impacting = interval_root_finder_BFS<is_vertex_face>(
                    a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1, tol,
                    tolerance, err, ms, t_max, max_itr, toi, output_tolerance);
                break;
            }
            assert(!tmp_is_impacting || toi >= 0);

            if (t_max == t_max_in) {
                // This will be the final output because we might need to
                // perform CCD again if the toi is zero. In which case we will
                // use a smaller t_max for more time resolution.
                is_impacting = tmp_is_impacting;
            } else {
                toi = tmp_is_impacting ? toi : t_max;
            }

            // This modification is for CCD-filtered line-search (e.g., IPC)
            // strategies for dealing with toi = 0:
            // 1. shrink t_max (when reaches max_itr),
            // 2. shrink tolerance (when not reach max_itr and tolerance is big) or
            // ms (when tolerance is too small comparing with ms)
            if (tmp_is_impacting && toi == 0 && no_zero_toi) {
                if (output_tolerance > tolerance) {
                    // reaches max_itr, so shrink t_max to return a more accurate result to reach target tolerance.
                    t_max *= 0.9;
                } else if (10 * tolerance < ms) {
                    ms *= 0.5; // ms is too large, shrink it
                } else {
                    tolerance *= 0.5; // tolerance is too large, shrink it

                    if constexpr (is_vertex_face) {
                        tol = compute_vertex_face_tolerances(
                            a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                            tolerance);
                    } else {
                        tol = compute_edge_edge_tolerances(
                            a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                            tolerance);
                    }
                }
            }

            // Only perform a second iteration if toi == 0.
            // WARNING: This option assumes the initial distance is not zero.
        } while (no_zero_toi && ++no_zero_toi_iter < MAX_NO_ZERO_TOI_ITER
                 && tmp_is_impacting && toi == 0);
        assert(!no_zero_toi || !is_impacting || toi != 0);

        return is_impacting;
    }

    bool edgeEdgeCCD(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        const Array3 &err_in,
        const Scalar ms_in,
        Scalar &toi,
        const Scalar tolerance_in,
        const Scalar t_max_in,
        const long max_itr,
        Scalar &output_tolerance,
        bool no_zero_toi,
        const CCDRootFindingMethod ccd_method)
    {
        return CCD</*is_vertex_face=*/false>(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
            err_in, ms_in, toi, tolerance_in, t_max_in, max_itr,
            output_tolerance, no_zero_toi, ccd_method);
    }

    bool vertexFaceCCD(
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1,
        const Array3 &err_in,
        const Scalar ms_in,
        Scalar &toi,
        const Scalar tolerance_in,
        const Scalar t_max_in,
        const long max_itr,
        Scalar &output_tolerance,
        bool no_zero_toi,
        const CCDRootFindingMethod ccd_method)
    {
        return CCD</*is_vertex_face=*/true>(
            v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, err_in, ms_in,
            toi, tolerance_in, t_max_in, max_itr, output_tolerance, no_zero_toi,
            ccd_method);
    }

#ifdef TIGHT_INCLUSION_FLOAT_WITH_DOUBLE_INPUT
    // these function are designed to test the performance of floating point vertion but with double inputs
    bool edgeEdgeCCD(
        const Eigen::Vector3d &ea0_t0,
        const Eigen::Vector3d &ea1_t0,
        const Eigen::Vector3d &eb0_t0,
        const Eigen::Vector3d &eb1_t0,
        const Eigen::Vector3d &ea0_t1,
        const Eigen::Vector3d &ea1_t1,
        const Eigen::Vector3d &eb0_t1,
        const Eigen::Vector3d &eb1_t1,
        const Eigen::Array3d &err,
        const double ms,
        double &toi,
        const double tolerance,
        const double t_max,
        const long max_itr,
        double &output_tolerance,
        bool no_zero_toi,
        const CCDRootFindingMethod ccd_method)
    {
        Scalar _toi = toi;
        Scalar _output_tolerance = output_tolerance;

        const bool result = edgeEdgeCCD(
            ea0_t0.cast<Scalar>(), ea1_t0.cast<Scalar>(), eb0_t0.cast<Scalar>(),
            eb1_t0.cast<Scalar>(), ea0_t1.cast<Scalar>(), ea1_t1.cast<Scalar>(),
            eb0_t1.cast<Scalar>(), eb1_t1.cast<Scalar>(), err.cast<Scalar>(),
            static_cast<Scalar>(ms), _toi, static_cast<Scalar>(tolerance),
            static_cast<Scalar>(t_max), max_itr, _output_tolerance, no_zero_toi,
            ccd_method);

        toi = _toi;
        output_tolerance = _output_tolerance;

        return result;
    }

    bool vertexFaceCCD(
        const Eigen::Vector3d &v_t0,
        const Eigen::Vector3d &f0_t0,
        const Eigen::Vector3d &f1_t0,
        const Eigen::Vector3d &f2_t0,
        const Eigen::Vector3d &v_t1,
        const Eigen::Vector3d &f0_t1,
        const Eigen::Vector3d &f1_t1,
        const Eigen::Vector3d &f2_t1,
        const Eigen::Array3d &err,
        const double ms,
        double &toi,
        const double tolerance,
        const double t_max,
        const long max_itr,
        double &output_tolerance,
        bool no_zero_toi,
        const CCDRootFindingMethod ccd_method)
    {
        Scalar _toi = toi;
        Scalar _output_tolerance = output_tolerance;

        const bool result = vertexFaceCCD(
            v_t0.cast<Scalar>(), f0_t0.cast<Scalar>(), f1_t0.cast<Scalar>(),
            f2_t0.cast<Scalar>(), v_t1.cast<Scalar>(), f0_t1.cast<Scalar>(),
            f1_t1.cast<Scalar>(), f2_t1.cast<Scalar>(), err.cast<Scalar>(),
            static_cast<Scalar>(ms), _toi, static_cast<Scalar>(tolerance),
            static_cast<Scalar>(t_max), max_itr, _output_tolerance, no_zero_toi,
            ccd_method);

        toi = _toi;
        output_tolerance = _output_tolerance;

        return result;
    }
#endif
} // namespace ticcd
