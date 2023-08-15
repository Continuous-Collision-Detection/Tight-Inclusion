#include "avx.hpp"

//#include <immintrin.h>

namespace ticcd {

    // __m512d function_f_ee_vector(
    //     __m512d ea0_t0,
    //     __m512d ea1_t0,
    //     __m512d eb0_t0,
    //     __m512d eb1_t0,
    //     __m512d ea0_t1,
    //     __m512d ea1_t1,
    //     __m512d eb0_t1,
    //     __m512d eb1_t1,
    //     __m512d t_up,
    //     __m512d t_dw,
    //     __m512d u_up,
    //     __m512d u_dw,
    //     __m512d v_up,
    //     __m512d v_dw)
    // {

    //     __m512d ea0 = _mm512_sub_pd(ea0_t1, ea0_t0);
    //     ea0 = _mm512_mul_pd(ea0, t_up);
    //     ea0 = _mm512_div_pd(ea0, t_dw);
    //     ea0 = _mm512_add_pd(ea0, ea0_t0);

    //     __m512d ea1 = _mm512_sub_pd(ea1_t1, ea1_t0);
    //     ea1 = _mm512_mul_pd(ea1, t_up);
    //     ea1 = _mm512_div_pd(ea1, t_dw);
    //     ea1 = _mm512_add_pd(ea1, ea1_t0);

    //     __m512d eb0 = _mm512_sub_pd(eb0_t1, eb0_t0);
    //     eb0 = _mm512_mul_pd(eb0, t_up);
    //     eb0 = _mm512_div_pd(eb0, t_dw);
    //     eb0 = _mm512_add_pd(eb0, eb0_t0);

    //     __m512d eb1 = _mm512_sub_pd(eb1_t1, eb1_t0);
    //     eb1 = _mm512_mul_pd(eb1, t_up);
    //     eb1 = _mm512_div_pd(eb1, t_dw);
    //     eb1 = _mm512_add_pd(eb1, eb1_t0);

    //     __m512d va = _mm512_sub_pd(ea1, ea0);
    //     va = _mm512_mul_pd(va, u_up);
    //     va = _mm512_div_pd(va, u_dw);
    //     va = _mm512_add_pd(va, ea0);

    //     __m512d vb = _mm512_sub_pd(eb1, eb0);
    //     vb = _mm512_mul_pd(vb, v_up);
    //     vb = _mm512_div_pd(vb, v_dw);
    //     vb = _mm512_add_pd(vb, eb0);

    //     return _mm512_sub_pd(va, vb);
    // }

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
        const std::array<Scalar, 8> &v_dw)
    {
        std::array<Scalar, 8> rst;
        for (int i = 0; i < 8; i++) {
            const Scalar ea0 = (ea0_t1 - ea0_t0) * t_up[i] / t_dw[i] + ea0_t0;
            const Scalar ea1 = (ea1_t1 - ea1_t0) * t_up[i] / t_dw[i] + ea1_t0;
            const Scalar eb0 = (eb0_t1 - eb0_t0) * t_up[i] / t_dw[i] + eb0_t0;
            const Scalar eb1 = (eb1_t1 - eb1_t0) * t_up[i] / t_dw[i] + eb1_t0;

            const Scalar va = (ea1 - ea0) * u_up[i] / u_dw[i] + ea0;
            const Scalar vb = (eb1 - eb0) * v_up[i] / v_dw[i] + eb0;
            rst[i] = va - vb;
        }
        return rst;
    }

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
        const std::array<Scalar, 8> &v_dw)
    {
        std::array<Scalar, 8> rst;
        for (int i = 0; i < 8; i++) {
            const Scalar v = (v_t1 - v_t0) * t_up[i] / t_dw[i] + v_t0;
            const Scalar f0 = (f0_t1 - f0_t0) * t_up[i] / t_dw[i] + f0_t0;
            const Scalar f1 = (f1_t1 - f1_t0) * t_up[i] / t_dw[i] + f1_t0;
            const Scalar f2 = (f2_t1 - f2_t0) * t_up[i] / t_dw[i] + f2_t0;
            const Scalar pt = (f1 - f0) * u_up[i] / u_dw[i]
                              + (f2 - f0) * v_up[i] / v_dw[i] + f0;
            rst[i] = v - pt;
        }
        return rst;
    }

    void convert_tuv_to_array(
        const Interval3 &itv,
        std::array<Scalar, 8> &t_up,
        std::array<Scalar, 8> &t_dw,
        std::array<Scalar, 8> &u_up,
        std::array<Scalar, 8> &u_dw,
        std::array<Scalar, 8> &v_up,
        std::array<Scalar, 8> &v_dw)
    {
        // t order: 0,0,0,0,1,1,1,1
        // u order: 0,0,1,1,0,0,1,1
        // v order: 0,1,0,1,0,1,0,1
        const Scalar t0_up = itv[0].lower.numerator;
        const Scalar t0_dw = itv[0].lower.denominator();
        const Scalar t1_up = itv[0].upper.numerator;
        const Scalar t1_dw = itv[0].upper.denominator();
        const Scalar u0_up = itv[1].lower.numerator;
        const Scalar u0_dw = itv[1].lower.denominator();
        const Scalar u1_up = itv[1].upper.numerator;
        const Scalar u1_dw = itv[1].upper.denominator();
        const Scalar v0_up = itv[2].lower.numerator;
        const Scalar v0_dw = itv[2].lower.denominator();
        const Scalar v1_up = itv[2].upper.numerator;
        const Scalar v1_dw = itv[2].upper.denominator();
        t_up = {{t0_up, t0_up, t0_up, t0_up, t1_up, t1_up, t1_up, t1_up}};
        t_dw = {{t0_dw, t0_dw, t0_dw, t0_dw, t1_dw, t1_dw, t1_dw, t1_dw}};
        u_up = {{u0_up, u0_up, u1_up, u1_up, u0_up, u0_up, u1_up, u1_up}};
        u_dw = {{u0_dw, u0_dw, u1_dw, u1_dw, u0_dw, u0_dw, u1_dw, u1_dw}};
        v_up = {{v0_up, v1_up, v0_up, v1_up, v0_up, v1_up, v0_up, v1_up}};
        v_dw = {{v0_dw, v1_dw, v0_dw, v1_dw, v0_dw, v1_dw, v0_dw, v1_dw}};
    }

    // __m512d function_f_vf_vector(
    //     __m512d v_t0,
    //     __m512d t0_t0,
    //     __m512d t1_t0,
    //     __m512d t2_t0,
    //     __m512d v_t1,
    //     __m512d t0_t1,
    //     __m512d t1_t1,
    //     __m512d t2_t1,
    //     __m512d t_up,
    //     __m512d t_dw,
    //     __m512d u_up,
    //     __m512d u_dw,
    //     __m512d v_up,
    //     __m512d v_dw)
    // {
    //     __m512d v = _mm512_sub_pd(v_t1, v_t0);
    //     v = _mm512_mul_pd(v, t_up);
    //     v = _mm512_div_pd(v, t_dw);
    //     v = _mm512_add_pd(v, v_t0);

    //     __m512d t0 = _mm512_sub_pd(f0_t1, f0_t0);
    //     t0 = _mm512_mul_pd(t0, t_up);
    //     t0 = _mm512_div_pd(t0, t_dw);
    //     t0 = _mm512_add_pd(t0, f0_t0);

    //     __m512d t1 = _mm512_sub_pd(f1_t1, f1_t0);
    //     t1 = _mm512_mul_pd(t1, t_up);
    //     t1 = _mm512_div_pd(t1, t_dw);
    //     t1 = _mm512_add_pd(t1, f1_t0);

    //     __m512d t2 = _mm512_sub_pd(f2_t1, f2_t0);
    //     t2 = _mm512_mul_pd(t2, t_up);
    //     t2 = _mm512_div_pd(t2, t_dw);
    //     t2 = _mm512_add_pd(t2, f2_t0);

    //     __m512d t01 = _mm512_sub_pd(t1, t0);
    //     t01 = _mm512_mul_pd(t01, u_up);
    //     t01 = _mm512_div_pd(t01, u_dw);

    //     __m512d t02 = _mm512_sub_pd(t2, t0);
    //     t02 = _mm512_mul_pd(t02, v_up);
    //     t02 = _mm512_div_pd(t02, v_dw);

    //     __m512d pt = _mm512_add_pd(t01, t02);
    //     pt = _mm512_add_pd(pt, t0);

    //     return _mm512_sub_pd(v, pt);
    // }

    // void convert_to_vector_pts(
    //     const Scalar asd,
    //     const Scalar bsd,
    //     const Scalar csd,
    //     const Scalar dsd,
    //     const Scalar aed,
    //     const Scalar bed,
    //     const Scalar ced,
    //     const Scalar ded,

    //     __m512d &as,
    //     __m512d &bs,
    //     __m512d &cs,
    //     __m512d &ds,
    //     __m512d &ae,
    //     __m512d &be,
    //     __m512d &ce,
    //     __m512d &de

    // )
    // {
    //     as = _mm512_setr_pd(asd, asd, asd, asd, asd, asd, asd, asd);
    //     bs = _mm512_setr_pd(bsd, bsd, bsd, bsd, bsd, bsd, bsd, bsd);
    //     cs = _mm512_setr_pd(csd, csd, csd, csd, csd, csd, csd, csd);
    //     ds = _mm512_setr_pd(dsd, dsd, dsd, dsd, dsd, dsd, dsd, dsd);

    //     ae = _mm512_setr_pd(aed, aed, aed, aed, aed, aed, aed, aed);
    //     be = _mm512_setr_pd(bed, bed, bed, bed, bed, bed, bed, bed);
    //     ce = _mm512_setr_pd(ced, ced, ced, ced, ced, ced, ced, ced);
    //     de = _mm512_setr_pd(ded, ded, ded, ded, ded, ded, ded, ded);
    // }
    // void convert_to_vector_pts_uvt(
    //     const std::array<double, 8> &t_up,
    //     const std::array<double, 8> &t_dw,
    //     const std::array<double, 8> &u_up,
    //     const std::array<double, 8> &u_dw,
    //     const std::array<double, 8> &v_up,
    //     const std::array<double, 8> &v_dw,

    //     __m512d &tu,
    //     __m512d &td,
    //     __m512d &uu,
    //     __m512d &ud,
    //     __m512d &vu,
    //     __m512d &vd)
    // {
    //     tu = _mm512_setr_pd(
    //         t_up[0], t_up[1], t_up[2], t_up[3], t_up[4], t_up[5], t_up[6],
    //         t_up[7]);
    //     td = _mm512_setr_pd(
    //         t_dw[0], t_dw[1], t_dw[2], t_dw[3], t_dw[4], t_dw[5], t_dw[6],
    //         t_dw[7]);

    //     uu = _mm512_setr_pd(
    //         u_up[0], u_up[1], u_up[2], u_up[3], u_up[4], u_up[5], u_up[6],
    //         u_up[7]);
    //     ud = _mm512_setr_pd(
    //         u_dw[0], u_dw[1], u_dw[2], u_dw[3], u_dw[4], u_dw[5], u_dw[6],
    //         u_dw[7]);

    //     vu = _mm512_setr_pd(
    //         v_up[0], v_up[1], v_up[2], v_up[3], v_up[4], v_up[5], v_up[6],
    //         v_up[7]);
    //     vd = _mm512_setr_pd(
    //         v_dw[0], v_dw[1], v_dw[2], v_dw[3], v_dw[4], v_dw[5], v_dw[6],
    //         v_dw[7]);
    // }
    // //(a-b)/c
    // __m512d function_test(__m512d a, __m512d b, __m512d c)
    // {
    //     __m512d v = _mm512_sub_pd(a, b);
    //     v = _mm512_div_pd(v, c);
    //     return v;
    // }

} // namespace ticcd
