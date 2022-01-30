// A root finder using interval arithmetic.
#pragma once

#include <array>
#include <functional>
#include <vector>
#include <tight_inclusion/types.hpp>
#include <tight_inclusion/interval.hpp>

namespace ticcd {
    // this version cannot give the impact time at t=1, although this collision can
    // be detected at t=0 of the next time step, but still may cause problems in
    // line-search based physical simulation

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
        Scalar &toi);

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
        Scalar &toi);

    // this version cannot give the impact time at t=1.
    // max_itr is a user defined maximum iteration time. if < 0, then
    // it will run until stack empty; otherwise the algorithm will stop when
    // iteration time reaches max_itr, and return a solution precision output_tolerance
    // it uses interval t = [0, max_time] instead of t = [0,1]
    // 0<=max_time <=1
    // tree searching order is horizontal

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
        Scalar &output_tolerance);

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
        Scalar &output_tolerance);

    // calculate the sign of f. dim is the dimension we are evaluating.
    template <typename T>
    T function_f_ee(
        const NumCCD &tpara,
        const NumCCD &upara,
        const NumCCD &vpara,
        const int dim,
        const Vector3 &a0s,
        const Vector3 &a1s,
        const Vector3 &b0s,
        const Vector3 &b1s,
        const Vector3 &a0e,
        const Vector3 &a1e,
        const Vector3 &b0e,
        const Vector3 &b1e);

    template <typename T>
    T function_f_vf(
        const NumCCD &tpara,
        const NumCCD &upara,
        const NumCCD &vpara,
        const int dim,
        const Vector3 &vs,
        const Vector3 &t0s,
        const Vector3 &t1s,
        const Vector3 &t2s,
        const Vector3 &ve,
        const Vector3 &t0e,
        const Vector3 &t1e,
        const Vector3 &t2e);

    void print_times();

    // get the filter of ccd. the inputs are the vertices of the bounding box of the simulation scene
    Array3 get_numerical_error(
        const std::vector<Vector3> &vertices,
        const bool check_vf,
        const bool using_minimum_separation);

} // namespace ticcd
