#include <iostream>
#include <iomanip>
#include <vector>

#include <tight_inclusion/ccd.hpp>
#include <tight_inclusion/interval_root_finder.hpp>
#include <tight_inclusion/timer.hpp>

namespace ticcd {

// #define TIGHT_INCLUSION_MAX_ABS_TOL
#ifdef TIGHT_INCLUSION_MAX_ABS_TOL
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

    Array3 compute_face_vertex_tolerances(
        const Vector3 &vs,
        const Vector3 &f0s,
        const Vector3 &f1s,
        const Vector3 &f2s,
        const Vector3 &ve,
        const Vector3 &f0e,
        const Vector3 &f1e,
        const Vector3 &f2e,
        const Scalar distance_tolerance)
    {
        const Vector3 p000 = vs - f0s;
        const Vector3 p001 = vs - f2s;
        const Vector3 p011 = vs - (f1s + f2s - f0s);
        const Vector3 p010 = vs - f1s;
        const Vector3 p100 = ve - f0e;
        const Vector3 p101 = ve - f2e;
        const Vector3 p111 = ve - (f1e + f2e - f0e);
        const Vector3 p110 = ve - f1e;

        Scalar dl =
            3 * max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110);
        Scalar edge0_length =
            3 * max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011);
        Scalar edge1_length =
            3 * max_linf_4(p000, p100, p110, p010, p001, p101, p111, p011);

        return Array3(
            clamp_div(distance_tolerance, dl, CCD_MAX_TIME_TOL),
            clamp_div(distance_tolerance, edge0_length, CCD_MAX_COORD_TOL),
            clamp_div(distance_tolerance, edge1_length, CCD_MAX_COORD_TOL));
    }

    Array3 compute_edge_edge_tolerances(
        const Vector3 &edge0_vertex0_start,
        const Vector3 &edge0_vertex1_start,
        const Vector3 &edge1_vertex0_start,
        const Vector3 &edge1_vertex1_start,
        const Vector3 &edge0_vertex0_end,
        const Vector3 &edge0_vertex1_end,
        const Vector3 &edge1_vertex0_end,
        const Vector3 &edge1_vertex1_end,
        const Scalar distance_tolerance)
    {

        const Vector3 p000 = edge0_vertex0_start - edge1_vertex0_start;
        const Vector3 p001 = edge0_vertex0_start - edge1_vertex1_start;
        const Vector3 p011 = edge0_vertex1_start - edge1_vertex1_start;
        const Vector3 p010 = edge0_vertex1_start - edge1_vertex0_start;
        const Vector3 p100 = edge0_vertex0_end - edge1_vertex0_end;
        const Vector3 p101 = edge0_vertex0_end - edge1_vertex1_end;
        const Vector3 p111 = edge0_vertex1_end - edge1_vertex1_end;
        const Vector3 p110 = edge0_vertex1_end - edge1_vertex0_end;

        Scalar dl =
            3 * max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110);
        Scalar edge0_length =
            3 * max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011);
        Scalar edge1_length =
            3 * max_linf_4(p000, p100, p110, p010, p001, p101, p111, p011);

        return Array3(
            clamp_div(distance_tolerance, dl, CCD_MAX_TIME_TOL),
            clamp_div(distance_tolerance, edge0_length, CCD_MAX_COORD_TOL),
            clamp_div(distance_tolerance, edge1_length, CCD_MAX_COORD_TOL));
    }

    bool edgeEdgeCCD(
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
        const Array3 &err_in,
        const Scalar ms_in,
        Scalar &toi,
        const Scalar tolerance_in,
        const Scalar t_max_in,
        const int max_itr,
        Scalar &output_tolerance,
        bool no_zero_toi,
        const CCDRootFindingMethod ccd_method)
    {
        const int MAX_NO_ZERO_TOI_ITER = std::numeric_limits<int>::max();
        // unsigned so can be larger than MAX_NO_ZERO_TOI_ITER
        unsigned int no_zero_toi_iter = 0;

        bool is_impacting, tmp_is_impacting;

        // Mutable copies for no_zero_toi
        Scalar t_max = t_max_in;
        Scalar tolerance = tolerance_in;
        Scalar ms = ms_in;

        Array3 tol = compute_edge_edge_tolerances(
            a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, tolerance_in);

        //////////////////////////////////////////////////////////
        // this should be the error of the whole mesh
        Array3 err;
        // if error[0]<0, means we need to calculate error here
        if (err_in[0] < 0) {
            std::vector<Vector3> vlist = {
                {a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e}};
            bool use_ms = ms > 0;
            err = get_numerical_error(vlist, false, use_ms);
        } else {
            err = err_in;
        }
        //////////////////////////////////////////////////////////

        do {
            switch (ccd_method) {
            case CCDRootFindingMethod::DEPTH_FIRST_SEARCH:
                // no handling for zero toi
                return edge_edge_interval_root_finder_DFS(
                    a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, tol, err, ms, toi);
            case CCDRootFindingMethod::BREADTH_FIRST_SEARCH:
                assert(t_max >= 0 && t_max <= 1);
                tmp_is_impacting = edge_edge_interval_root_finder_BFS(
                    a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, tol, tolerance, err,
                    ms, t_max, max_itr, toi, output_tolerance);
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

                    tol = compute_edge_edge_tolerances(
                        a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, tolerance);
                }
            }

            // Only perform a second iteration if toi == 0.
            // WARNING: This option assumes the initial distance is not zero.
        } while (no_zero_toi && ++no_zero_toi_iter < MAX_NO_ZERO_TOI_ITER
                 && tmp_is_impacting && toi == 0);
        assert(!no_zero_toi || !is_impacting || toi != 0);

        return is_impacting;
    }

    bool vertexFaceCCD(
        const Vector3 &vertex_start,
        const Vector3 &face_vertex0_start,
        const Vector3 &face_vertex1_start,
        const Vector3 &face_vertex2_start,
        const Vector3 &vertex_end,
        const Vector3 &face_vertex0_end,
        const Vector3 &face_vertex1_end,
        const Vector3 &face_vertex2_end,
        const Array3 &err_in,
        const Scalar ms_in,
        Scalar &toi,
        const Scalar tolerance_in,
        const Scalar t_max_in,
        const int max_itr,
        Scalar &output_tolerance,
        bool no_zero_toi,
        const CCDRootFindingMethod ccd_method)
    {
        const int MAX_NO_ZERO_TOI_ITER = std::numeric_limits<int>::max();
        // unsigned so can be larger than MAX_NO_ZERO_TOI_ITER
        unsigned int no_zero_toi_iter = 0;

        bool is_impacting, tmp_is_impacting;

        // Mutable copies for no_zero_toi
        Scalar t_max = t_max_in;
        Scalar tolerance = tolerance_in;
        Scalar ms = ms_in;

        Array3 tol = compute_face_vertex_tolerances(
            vertex_start, face_vertex0_start, face_vertex1_start,
            face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
            face_vertex2_end, tolerance);

        //////////////////////////////////////////////////////////
        // this is the error of the whole mesh
        Array3 err;
        // if error[0]<0, means we need to calculate error here
        if (err_in[0] < 0) {
            std::vector<Vector3> vlist = {
                {vertex_start, face_vertex0_start, face_vertex1_start,
                 face_vertex2_start, vertex_end, face_vertex0_end,
                 face_vertex1_end, face_vertex2_end}};
            bool use_ms = ms > 0;
            err = get_numerical_error(vlist, true, use_ms);
        } else {
            err = err_in;
        }
        //////////////////////////////////////////////////////////

        do {
            switch (ccd_method) {
            case CCDRootFindingMethod::DEPTH_FIRST_SEARCH:
                // no handling for zero toi
                return vertex_face_interval_root_finder_DFS(
                    vertex_start, face_vertex0_start, face_vertex1_start,
                    face_vertex2_start, vertex_end, face_vertex0_end,
                    face_vertex1_end, face_vertex2_end, tol, err, ms, toi);
            case CCDRootFindingMethod::BREADTH_FIRST_SEARCH:
                assert(t_max >= 0 && t_max <= 1);
                tmp_is_impacting = vertex_face_interval_root_finder_BFS(
                    vertex_start, face_vertex0_start, face_vertex1_start,
                    face_vertex2_start, vertex_end, face_vertex0_end,
                    face_vertex1_end, face_vertex2_end, tol, tolerance, err, ms,
                    t_max, max_itr, toi, output_tolerance);
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

                    // recompute this
                    tol = compute_face_vertex_tolerances(
                        vertex_start, face_vertex0_start, face_vertex1_start,
                        face_vertex2_start, vertex_end, face_vertex0_end,
                        face_vertex1_end, face_vertex2_end, tolerance);
                }
            }

            // Only perform a second iteration if toi == 0.
            // WARNING: This option assumes the initial distance is not zero.
        } while (no_zero_toi && ++no_zero_toi_iter < MAX_NO_ZERO_TOI_ITER
                 && tmp_is_impacting && toi == 0);
        assert(!no_zero_toi || !is_impacting || toi != 0);

        return is_impacting;
    }

#ifdef TIGHT_INCLUSION_FWDI
    // these function are designed to test the performance of floating point vertion but with double inputs
    bool edgeEdgeCCD(
        const Eigen::Vector3d &a0_start,
        const Eigen::Vector3d &a1_start,
        const Eigen::Vector3d &b0_start,
        const Eigen::Vector3d &b1_start,
        const Eigen::Vector3d &a0_end,
        const Eigen::Vector3d &a1_end,
        const Eigen::Vector3d &b0_end,
        const Eigen::Vector3d &b1_end,
        const std::array<double, 3> &err,
        const double ms,
        double &toi,
        const double tolerance,
        const double t_max,
        const int max_itr,
        double &output_tolerance,
        bool no_zero_toi,
        const CCDRootFindingMethod ccd_method)
    {
        Vector3 fa0_start = a0_start.cast<Scalar>();
        Vector3 fa1_start = a1_start.cast<Scalar>();
        Vector3 fb0_start = b0_start.cast<Scalar>();
        Vector3 fb1_start = b1_start.cast<Scalar>();
        Vector3 fa0_end = a0_end.cast<Scalar>();
        Vector3 fa1_end = a1_end.cast<Scalar>();
        Vector3 fb0_end = b0_end.cast<Scalar>();
        Vector3 fb1_end = b1_end.cast<Scalar>();

        Array3 ferr = {{
            Scalar(err[0]),
            Scalar(err[1]),
            Scalar(err[2]),
        }};

        Scalar fms = ms;
        Scalar ftoi = toi;
        Scalar ftolerance = tolerance;
        Scalar ft_max = t_max;
        Scalar fouttol = output_tolerance;

        bool result = edgeEdgeCCD(
            fa0_start, fa1_start, fb0_start, fb1_start, fa0_end, fa1_end,
            fb0_end, fb1_end, ferr, fms, ftoi, ftolerance, ft_max, max_itr,
            fouttol, no_zero_toi, ccd_method);

        toi = ftoi;
        output_tolerance = fouttol;

        return result;
    }

    bool vertexFaceCCD(
        const Eigen::Vector3d &vertex_start,
        const Eigen::Vector3d &face_vertex0_start,
        const Eigen::Vector3d &face_vertex1_start,
        const Eigen::Vector3d &face_vertex2_start,
        const Eigen::Vector3d &vertex_end,
        const Eigen::Vector3d &face_vertex0_end,
        const Eigen::Vector3d &face_vertex1_end,
        const Eigen::Vector3d &face_vertex2_end,
        const std::array<double, 3> &err,
        const double ms,
        double &toi,
        const double tolerance,
        const double t_max,
        const int max_itr,
        double &output_tolerance,
        bool no_zero_toi,
        const CCDRootFindingMethod ccd_method)
    {
        Vector3 fvertex_start = vertex_start.cast<Scalar>();
        Vector3 fface_vertex0_start = face_vertex0_start.cast<Scalar>();
        Vector3 fface_vertex1_start = face_vertex1_start.cast<Scalar>();
        Vector3 fface_vertex2_start = face_vertex2_start.cast<Scalar>();
        Vector3 fvertex_end = vertex_end.cast<Scalar>();
        Vector3 fface_vertex0_end = face_vertex0_end.cast<Scalar>();
        Vector3 fface_vertex1_end = face_vertex1_end.cast<Scalar>();
        Vector3 fface_vertex2_end = face_vertex2_end.cast<Scalar>();

        Array3 ferr = {{
            Scalar(err[0]),
            Scalar(err[1]),
            Scalar(err[2]),
        }};

        Scalar fms = ms;
        Scalar ftoi = toi;
        Scalar ftolerance = tolerance;
        Scalar ft_max = t_max;
        Scalar fouttol = output_tolerance;

        bool result = vertexFaceCCD(
            fvertex_start, fface_vertex0_start, fface_vertex1_start,
            fface_vertex2_start, fvertex_end, fface_vertex0_end,
            fface_vertex1_end, fface_vertex2_end, ferr, fms, ftoi, ftolerance,
            ft_max, max_itr, fouttol, no_zero_toi, ccd_method);

        toi = ftoi;
        output_tolerance = fouttol;

        return result;
    }
#endif
} // namespace ticcd
