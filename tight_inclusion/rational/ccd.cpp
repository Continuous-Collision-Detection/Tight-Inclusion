#include <iostream>
#include <iomanip>
#include <vector>

#include <tight_inclusion/ccd.hpp>
#include <tight_inclusion/rational/ccd.hpp>
#include <tight_inclusion/interval_root_finder.hpp>
#include <tight_inclusion/rational/interval_root_finder.hpp>

namespace ticcd::rational {

    bool edgeEdgeCCD(
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {

        Array3 tol = compute_edge_edge_tolerances(
            a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, DEFAULT_CCD_DISTANCE_TOL);

        //////////////////////////////////////////////////////////
        // TODO this should be the error of the whole mesh
        std::vector<Vector3> vlist = {{a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e}};

        bool use_ms = ms > 0;
        Array3 auto_err = get_numerical_error(vlist, false, use_ms);
        //////////////////////////////////////////////////////////

        std::array<RationalInterval, 3> toi_interval;

        bool is_impacting = edge_edge_interval_root_finder(
            a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, tol, auto_err, ms,
            toi_interval);

        // Return a conservative time-of-impact
        if (is_impacting) {
            toi = toi_interval[0][0].to_double();
        }
        // This time of impact is very dangerous for convergence
        // assert(!is_impacting || toi > 0);
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
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {
        Array3 tol = compute_face_vertex_tolerances(
            vertex_start, face_vertex0_start, face_vertex1_start,
            face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
            face_vertex2_end, DEFAULT_CCD_DISTANCE_TOL);
        // std::cout<<"get tolerance successfully"<<std::endl;
        //////////////////////////////////////////////////////////
        // TODO this should be the error of the whole mesh
        std::vector<Vector3> vlist = {
            {vertex_start, face_vertex0_start, face_vertex1_start,
             face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
             face_vertex2_end}};

        bool use_ms = ms > 0;
        Array3 auto_err = get_numerical_error(vlist, false, use_ms);
        // std::cout<<"get error successfully"<<std::endl;
        //////////////////////////////////////////////////////////

        std::array<RationalInterval, 3> toi_interval;

        bool is_impacting = vertex_face_interval_root_finder(
            vertex_start, face_vertex0_start, face_vertex1_start,
            face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
            face_vertex2_end, tol, auto_err, ms, toi_interval);

        // std::cout<<"get result successfully"<<std::endl;
        // Return a conservative time-of-impact
        if (is_impacting) {
            toi = toi_interval[0][0].to_double();
        }
        // std::cout<<"get time successfully"<<std::endl;
        // This time of impact is very dangerous for convergence
        // assert(!is_impacting || toi > 0);
        return is_impacting;
    }

} // namespace ticcd::rational
