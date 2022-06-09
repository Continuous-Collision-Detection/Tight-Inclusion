#pragma once

#include <tight_inclusion/interval.hpp>

namespace ticcd {
    void convert_tuv_to_array(
        const Interval3 &itv,
        std::array<Scalar, 8> &t_up,
        std::array<Scalar, 8> &t_dw,
        std::array<Scalar, 8> &u_up,
        std::array<Scalar, 8> &u_dw,
        std::array<Scalar, 8> &v_up,
        std::array<Scalar, 8> &v_dw);

    std::array<Scalar, 8> function_vf(
        const Scalar &vs,
        const Scalar &t0s,
        const Scalar &t1s,
        const Scalar &t2s,
        const Scalar &ve,
        const Scalar &t0e,
        const Scalar &t1e,
        const Scalar &t2e,
        const std::array<Scalar, 8> &t_up,
        const std::array<Scalar, 8> &t_dw,
        const std::array<Scalar, 8> &u_up,
        const std::array<Scalar, 8> &u_dw,
        const std::array<Scalar, 8> &v_up,
        const std::array<Scalar, 8> &v_dw);

    std::array<Scalar, 8> function_ee(
        const Scalar &a0s,
        const Scalar &a1s,
        const Scalar &b0s,
        const Scalar &b1s,
        const Scalar &a0e,
        const Scalar &a1e,
        const Scalar &b0e,
        const Scalar &b1e,
        const std::array<Scalar, 8> &t_up,
        const std::array<Scalar, 8> &t_dw,
        const std::array<Scalar, 8> &u_up,
        const std::array<Scalar, 8> &u_dw,
        const std::array<Scalar, 8> &v_up,
        const std::array<Scalar, 8> &v_dw);
} // namespace ticcd
