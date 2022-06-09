// Time-of-impact computation for rigid bodies with angular trajectories.
#pragma once

#include <tight_inclusion/ccd.hpp>

namespace ticcd::rational {

    // this version is an naive implementation of Tight-Inclusion CCD without optimizations
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
        Scalar &toi);

    // this version is an naive implementation of Tight-Inclusion CCD without optimizations
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
        Scalar &toi);

} // namespace ticcd::rational
