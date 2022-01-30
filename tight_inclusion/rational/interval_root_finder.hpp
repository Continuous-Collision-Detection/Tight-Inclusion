// A root finder using interval arithmetic.
#pragma once

#include <array>
#include <functional>
#include <vector>
#include <tight_inclusion/types.hpp>
#include <tight_inclusion/interval.hpp>
#include <tight_inclusion/rational/rational.hpp>

namespace ticcd::rational {

    typedef std::array<Rational, 2> RationalInterval;

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
        std::array<RationalInterval, 3> &root);

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
        std::array<RationalInterval, 3> &root);

    Scalar print_time_rational();

} // namespace ticcd::rational
