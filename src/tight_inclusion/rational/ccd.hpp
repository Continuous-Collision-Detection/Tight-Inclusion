// Time-of-impact computation for rigid bodies with angular trajectories.
#pragma once

#include <tight_inclusion/ccd.hpp>

namespace ticcd::rational {

    // this version is an naive implementation of Tight-Inclusion CCD without optimizations
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
        Scalar &toi);

    // this version is an naive implementation of Tight-Inclusion CCD without optimizations
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
        Scalar &toi);

} // namespace ticcd::rational
