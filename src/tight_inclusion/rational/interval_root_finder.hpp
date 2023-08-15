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
        std::array<RationalInterval, 3> &root);

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
        std::array<RationalInterval, 3> &root);

    Scalar print_time_rational();

} // namespace ticcd::rational
