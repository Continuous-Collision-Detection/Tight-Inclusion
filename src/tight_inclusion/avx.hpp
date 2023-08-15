#pragma once

#include <tight_inclusion/interval.hpp>

#include <array>

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
        const Scalar &v_t0,
        const Scalar &f0_t0,
        const Scalar &f1_t0,
        const Scalar &f2_t0,
        const Scalar &v_t1,
        const Scalar &f0_t1,
        const Scalar &f1_t1,
        const Scalar &f2_t1,
        const std::array<Scalar, 8> &t_up,
        const std::array<Scalar, 8> &t_dw,
        const std::array<Scalar, 8> &u_up,
        const std::array<Scalar, 8> &u_dw,
        const std::array<Scalar, 8> &v_up,
        const std::array<Scalar, 8> &v_dw);

    std::array<Scalar, 8> function_ee(
        const Scalar &ea0_t0,
        const Scalar &ea1_t0,
        const Scalar &eb0_t0,
        const Scalar &eb1_t0,
        const Scalar &ea0_t1,
        const Scalar &ea1_t1,
        const Scalar &eb0_t1,
        const Scalar &eb1_t1,
        const std::array<Scalar, 8> &t_up,
        const std::array<Scalar, 8> &t_dw,
        const std::array<Scalar, 8> &u_up,
        const std::array<Scalar, 8> &u_dw,
        const std::array<Scalar, 8> &v_up,
        const std::array<Scalar, 8> &v_dw);
} // namespace ticcd
