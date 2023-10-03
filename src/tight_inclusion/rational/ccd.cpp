#include <iomanip>
#include <vector>

#include <tight_inclusion/ccd.hpp>
#include <tight_inclusion/rational/ccd.hpp>
#include <tight_inclusion/interval_root_finder.hpp>
#include <tight_inclusion/rational/interval_root_finder.hpp>

namespace ticcd::rational {

    bool edgeEdgeCCD(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {

        Array3 tol = compute_edge_edge_tolerances(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
            DEFAULT_CCD_DISTANCE_TOL);

        //////////////////////////////////////////////////////////
        // TODO this should be the error of the whole mesh
        std::vector<Vector3> vlist = {
            {ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1}};

        bool use_ms = ms > 0;
        Array3 auto_err = get_numerical_error(vlist, false, use_ms);
        //////////////////////////////////////////////////////////

        std::array<RationalInterval, 3> toi_interval;

        bool is_impacting = edge_edge_interval_root_finder(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, tol,
            auto_err, ms, toi_interval);

        // Return a conservative time-of-impact
        if (is_impacting) {
            toi = toi_interval[0][0];
        }
        // This time of impact is very dangerous for convergence
        // assert(!is_impacting || toi > 0);
        return is_impacting;
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
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {
        Array3 tol = compute_vertex_face_tolerances(
            v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1,
            DEFAULT_CCD_DISTANCE_TOL);

        //////////////////////////////////////////////////////////
        // TODO this should be the error of the whole mesh
        std::vector<Vector3> vlist = {
            {v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1}};

        bool use_ms = ms > 0;
        Array3 auto_err = get_numerical_error(vlist, false, use_ms);
        //////////////////////////////////////////////////////////

        std::array<RationalInterval, 3> toi_interval;

        bool is_impacting = vertex_face_interval_root_finder(
            v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, tol, auto_err,
            ms, toi_interval);

        // Return a conservative time-of-impact
        if (is_impacting) {
            toi = toi_interval[0][0];
        }

        // This time of impact is very dangerous for convergence
        // assert(!is_impacting || toi > 0);
        return is_impacting;
    }

} // namespace ticcd::rational
