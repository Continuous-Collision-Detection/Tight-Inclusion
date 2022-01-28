//#include <immintrin.h>
#include <array>

#include <tight_inclusion/avx.hpp>

namespace ticcd {

    // __m512d function_f_ee_vector(
    //     __m512d a0s,__m512d a1s,__m512d b0s,__m512d b1s,
    //     __m512d a0e,__m512d a1e,__m512d b0e,__m512d b1e,
    //     __m512d t_up,__m512d t_dw,__m512d u_up,__m512d u_dw,__m512d v_up,__m512d v_dw){

    //         __m512d edge0_vertex0=_mm512_sub_pd(a0e,a0s);
    //         edge0_vertex0=_mm512_mul_pd(edge0_vertex0,t_up);
    //         edge0_vertex0=_mm512_div_pd(edge0_vertex0,t_dw);
    //         edge0_vertex0=_mm512_add_pd(edge0_vertex0,a0s);

    //         __m512d edge0_vertex1=_mm512_sub_pd(a1e,a1s);
    //         edge0_vertex1=_mm512_mul_pd(edge0_vertex1,t_up);
    //         edge0_vertex1=_mm512_div_pd(edge0_vertex1,t_dw);
    //         edge0_vertex1=_mm512_add_pd(edge0_vertex1,a1s);

    //         __m512d edge1_vertex0=_mm512_sub_pd(b0e,b0s);
    //         edge1_vertex0=_mm512_mul_pd(edge1_vertex0,t_up);
    //         edge1_vertex0=_mm512_div_pd(edge1_vertex0,t_dw);
    //         edge1_vertex0=_mm512_add_pd(edge1_vertex0,b0s);

    //         __m512d edge1_vertex1=_mm512_sub_pd(b1e,b1s);
    //         edge1_vertex1=_mm512_mul_pd(edge1_vertex1,t_up);
    //         edge1_vertex1=_mm512_div_pd(edge1_vertex1,t_dw);
    //         edge1_vertex1=_mm512_add_pd(edge1_vertex1,b1s);

    //         __m512d edge0_vertex=_mm512_sub_pd(edge0_vertex1,edge0_vertex0);
    //         edge0_vertex=_mm512_mul_pd(edge0_vertex,u_up);
    //         edge0_vertex=_mm512_div_pd(edge0_vertex,u_dw);
    //         edge0_vertex=_mm512_add_pd(edge0_vertex,edge0_vertex0);

    //         __m512d edge1_vertex=_mm512_sub_pd(edge1_vertex1,edge1_vertex0);
    //         edge1_vertex=_mm512_mul_pd(edge1_vertex,v_up);
    //         edge1_vertex=_mm512_div_pd(edge1_vertex,v_dw);
    //         edge1_vertex=_mm512_add_pd(edge1_vertex,edge1_vertex0);

    //         return _mm512_sub_pd(edge1_vertex,edge0_vertex);
    //     }

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
        const std::array<Scalar, 8> &v_dw)
    {
        std::array<Scalar, 8> rst;
        for (int i = 0; i < 8; i++) {
            Scalar edge0_vertex0 = (a0e - a0s) * t_up[i] / t_dw[i] + a0s;
            Scalar edge0_vertex1 = (a1e - a1s) * t_up[i] / t_dw[i] + a1s;
            Scalar edge1_vertex0 = (b0e - b0s) * t_up[i] / t_dw[i] + b0s;
            Scalar edge1_vertex1 = (b1e - b1s) * t_up[i] / t_dw[i] + b1s;

            Scalar edge0_vertex =
                (edge0_vertex1 - edge0_vertex0) * u_up[i] / u_dw[i]
                + edge0_vertex0;
            Scalar edge1_vertex =
                (edge1_vertex1 - edge1_vertex0) * v_up[i] / v_dw[i]
                + edge1_vertex0;
            rst[i] = edge0_vertex - edge1_vertex;
        }
        return rst;
    }

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
        const std::array<Scalar, 8> &v_dw)
    {
        std::array<Scalar, 8> rst;
        for (int i = 0; i < 8; i++) {
            Scalar v = (ve - vs) * t_up[i] / t_dw[i] + vs;
            Scalar t0 = (t0e - t0s) * t_up[i] / t_dw[i] + t0s;
            Scalar t1 = (t1e - t1s) * t_up[i] / t_dw[i] + t1s;
            Scalar t2 = (t2e - t2s) * t_up[i] / t_dw[i] + t2s;
            Scalar pt = (t1 - t0) * u_up[i] / u_dw[i]
                        + (t2 - t0) * v_up[i] / v_dw[i] + t0;
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
        const Scalar t0_up = itv[0].lower.numerator,
                     t0_dw = itv[0].lower.denominator(),
                     t1_up = itv[0].upper.numerator,
                     t1_dw = itv[0].upper.denominator(),
                     u0_up = itv[1].lower.numerator,
                     u0_dw = itv[1].lower.denominator(),
                     u1_up = itv[1].upper.numerator,
                     u1_dw = itv[1].upper.denominator(),
                     v0_up = itv[2].lower.numerator,
                     v0_dw = itv[2].lower.denominator(),
                     v1_up = itv[2].upper.numerator,
                     v1_dw = itv[2].upper.denominator();
        t_up = {{t0_up, t0_up, t0_up, t0_up, t1_up, t1_up, t1_up, t1_up}};
        t_dw = {{t0_dw, t0_dw, t0_dw, t0_dw, t1_dw, t1_dw, t1_dw, t1_dw}};
        u_up = {{u0_up, u0_up, u1_up, u1_up, u0_up, u0_up, u1_up, u1_up}};
        u_dw = {{u0_dw, u0_dw, u1_dw, u1_dw, u0_dw, u0_dw, u1_dw, u1_dw}};
        v_up = {{v0_up, v1_up, v0_up, v1_up, v0_up, v1_up, v0_up, v1_up}};
        v_dw = {{v0_dw, v1_dw, v0_dw, v1_dw, v0_dw, v1_dw, v0_dw, v1_dw}};
    }

    //     __m512d function_f_vf_vector(
    //     __m512d vs,__m512d t0s,__m512d t1s,__m512d t2s,
    //     __m512d ve,__m512d t0e,__m512d t1e,__m512d t2e,
    //     __m512d t_up,__m512d t_dw,__m512d u_up,__m512d u_dw,__m512d v_up,__m512d v_dw){
    //          __m512d v=_mm512_sub_pd(ve,vs);
    //         v=_mm512_mul_pd(v,t_up);
    //         v=_mm512_div_pd(v,t_dw);
    //         v=_mm512_add_pd(v,vs);

    //          __m512d t0=_mm512_sub_pd(t0e,t0s);
    //         t0=_mm512_mul_pd(t0,t_up);
    //         t0=_mm512_div_pd(t0,t_dw);
    //         t0=_mm512_add_pd(t0,t0s);

    //          __m512d t1=_mm512_sub_pd(t1e,t1s);
    //         t1=_mm512_mul_pd(t1,t_up);
    //         t1=_mm512_div_pd(t1,t_dw);
    //         t1=_mm512_add_pd(t1,t1s);

    //          __m512d t2=_mm512_sub_pd(t2e,t2s);
    //         t2=_mm512_mul_pd(t2,t_up);
    //         t2=_mm512_div_pd(t2,t_dw);
    //         t2=_mm512_add_pd(t2,t2s);

    //          __m512d t01=_mm512_sub_pd(t1,t0);
    //         t01=_mm512_mul_pd(t01,u_up);
    //         t01=_mm512_div_pd(t01,u_dw);

    //         __m512d t02=_mm512_sub_pd(t2,t0);
    //         t02=_mm512_mul_pd(t02,v_up);
    //         t02=_mm512_div_pd(t02,v_dw);

    //        __m512d pt=_mm512_add_pd(t01,t02);
    //         pt=_mm512_add_pd(pt,t0);

    //        return _mm512_sub_pd(v,pt);

    // }

    // void convert_to_vector_pts(
    //     const Scalar asd,const Scalar bsd,const Scalar csd,const Scalar dsd,
    //     const Scalar aed,const Scalar bed,const Scalar ced,const Scalar ded,

    //     __m512d& as, __m512d& bs, __m512d& cs, __m512d& ds,
    //     __m512d& ae, __m512d& be, __m512d& ce, __m512d& de

    //     ){
    //         as=_mm512_setr_pd(asd,asd,asd,asd,asd,asd,asd,asd);
    //         bs=_mm512_setr_pd(bsd,bsd,bsd,bsd,bsd,bsd,bsd,bsd);
    //         cs=_mm512_setr_pd(csd,csd,csd,csd,csd,csd,csd,csd);
    //         ds=_mm512_setr_pd(dsd,dsd,dsd,dsd,dsd,dsd,dsd,dsd);

    //         ae=_mm512_setr_pd(aed,aed,aed,aed,aed,aed,aed,aed);
    //         be=_mm512_setr_pd(bed,bed,bed,bed,bed,bed,bed,bed);
    //         ce=_mm512_setr_pd(ced,ced,ced,ced,ced,ced,ced,ced);
    //         de=_mm512_setr_pd(ded,ded,ded,ded,ded,ded,ded,ded);
    //     }
    // void convert_to_vector_pts_uvt(
    //     const std::array<double,8> &t_up, const std::array<double,8> &t_dw,
    //     const std::array<double,8> &u_up, const std::array<double,8> &u_dw,
    //     const std::array<double,8> &v_up, const std::array<double,8> &v_dw,

    //     __m512d& tu, __m512d& td, __m512d& uu, __m512d& ud,__m512d& vu, __m512d& vd
    //     ){
    //         tu=_mm512_setr_pd(t_up[0],t_up[1],t_up[2],t_up[3],t_up[4],t_up[5],t_up[6],t_up[7]);
    //         td=_mm512_setr_pd(t_dw[0],t_dw[1],t_dw[2],t_dw[3],t_dw[4],t_dw[5],t_dw[6],t_dw[7]);

    //         uu=_mm512_setr_pd(u_up[0],u_up[1],u_up[2],u_up[3],u_up[4],u_up[5],u_up[6],u_up[7]);
    //         ud=_mm512_setr_pd(u_dw[0],u_dw[1],u_dw[2],u_dw[3],u_dw[4],u_dw[5],u_dw[6],u_dw[7]);

    //         vu=_mm512_setr_pd(v_up[0],v_up[1],v_up[2],v_up[3],v_up[4],v_up[5],v_up[6],v_up[7]);
    //         vd=_mm512_setr_pd(v_dw[0],v_dw[1],v_dw[2],v_dw[3],v_dw[4],v_dw[5],v_dw[6],v_dw[7]);
    //     }
    // //(a-b)/c
    // __m512d function_test(__m512d a,__m512d b, __m512d c){
    //     __m512d v=_mm512_sub_pd(a,b);
    //     v= _mm512_div_pd(v,c);
    //     return v;
    // }

} // namespace ticcd
