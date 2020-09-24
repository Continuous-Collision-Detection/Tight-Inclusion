// A root finder using interval arithmetic.
#include <tight_inclusion/interval_root_finder.hpp>

#include <stack>
#include <tight_inclusion/igl-Timer.h>
#include <iostream>
#include <tight_inclusion/Rational.hpp>
#include <tight_inclusion/avx.h>
#include <queue>
#include <fstream>
// #define COMPARE_WITH_RATIONAL

// #define DEBUGING
namespace inclusion_ccd
{
    double time20 = 0, time21 = 0, time22 = 0, time23 = 0, time24 = 0, time25 = 0, time_rational = 0;
    int refine = 0;
    int refine_return = 0;

    // convert Numccd to double number
    double Numccd2double(const Numccd &n)
    {
        double r = double(n.first) / power(1, n.second);
        return r;
    }

#ifdef TIGHT_INCLUSION_USE_GMP
    std::array<Rational, 3> width(const std::array<std::pair<Rational, Rational>, 3> &x)
    {
        std::array<Rational, 3> w;

        for (int i = 0; i < 3; i++)
        {
            Rational sub = x[i].first - x[i].second;
            w[i] = sub >= 0 ? sub : -sub;
            assert(w[i] >= 0);
        }
        return w;
    }
#endif

    Eigen::VectorX3d width(const Interval3 &x)
    {
        Eigen::VectorX3d w;
        w.resize(3);
        for (int i = 0; i < 3; i++)
        {
            w[i] = Numccd2double(x[i].second) - Numccd2double(x[i].first);
            assert(w[i] >= 0);
        }
        return w;
    }

    std::array<Eigen::Vector3d, 2> bbd_4_pts(const std::array<Eigen::Vector3d, 4> &pts)
    {
        Eigen::Vector3d min, max;
        min = pts[0];
        max = pts[0];
        for (int i = 1; i < 4; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (min[j] > pts[i][j])
                {
                    min[j] = pts[i][j];
                }
                if (max[j] < pts[i][j])
                {
                    max[j] = pts[i][j];
                }
            }
        }
        std::array<Eigen::Vector3d, 2> rst;
        rst[0] = min;
        rst[1] = max;
        return rst;
    }
    std::array<Eigen::Vector3d, 2> bbd_6_pts(const std::array<Eigen::Vector3d, 6> &pts)
    {
        Eigen::Vector3d min, max;
        min = pts[0];
        max = pts[0];
        for (int i = 1; i < 6; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (min[j] > pts[i][j])
                {
                    min[j] = pts[i][j];
                }
                if (max[j] < pts[i][j])
                {
                    max[j] = pts[i][j];
                }
            }
        }
        std::array<Eigen::Vector3d, 2> rst;
        rst[0] = min;
        rst[1] = max;
        return rst;
    }

    // eps is the interval [-eps,eps] we need to check
    // if [-eps,eps] overlap, return true
    template <typename T>
    bool evaluate_bbox_one_dimension(
        const std::array<Numccd, 2> &t,
        const std::array<Numccd, 2> &u,
        const std::array<Numccd, 2> &v,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const int dimension, T tp, const bool check_vf, const double eps)
    {
#ifdef TIGHT_INCLUSION_USE_TIMER
        igl::Timer timer;
#endif

        /*
{// smart way but no possibility of parallazation
    double eva;
    bool flag0=false, flag1=false;
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){
                if(!check_vf){
#ifdef TIGHT_INCLUSION_USE_TIMER
                    timer.start();
#endif
                    eva=function_f_ee(t[i],u[j],v[k],tp,dimension,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
#ifdef TIGHT_INCLUSION_USE_TIMER
                    timer.stop();
                    time25+=timer.getElapsedTimeInMicroSec();
#endif
                }
                else{
#ifdef TIGHT_INCLUSION_USE_TIMER
                    timer.start();
#endif
                    eva=function_f_vf(t[i],u[j],v[k],tp,dimension,a0s,a1s,b0s,b1s,a0e,a1e,b0e,b1e);
#ifdef TIGHT_INCLUSION_USE_TIMER
                    timer.stop();
                    time25+=timer.getElapsedTimeInMicroSec();
#endif
                }
                if(eva<=eps&&eva>=-eps){
                    return true;
                }
                if(eva<-eps){
                    flag0=true;
                }
                if(eva>eps){
                    flag1=true;
                }
                if(flag0&&flag1){
                    return true;
                }
            }
        }
    }
    if(flag0&&flag1)
    return true;
    return false;
}
*/
        if (check_vf)
        { // test
            std::array<double, 8> vs;
            int count = 0;
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.start();
#endif
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {

                        vs[count] = function_f_vf(t[i], u[j], v[k], tp, dimension, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);

                        count++;
                    }
                }
            }
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.stop();
            time25 += timer.getElapsedTimeInMicroSec();
#endif

            double minv = vs[0], maxv = vs[0];

            for (int i = 1; i < 8; i++)
            {
                if (minv > vs[i])
                {
                    minv = vs[i];
                }
                if (maxv < vs[i])
                {
                    maxv = vs[i];
                }
            }
            if (minv > eps || maxv < -eps)
                return false;
            return true;

        } // test end

        else
        { // test
            std::array<double, 8> vs;
            int count = 0;
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.start();
#endif
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {

                        vs[count] = function_f_ee(t[i], u[j], v[k], tp, dimension, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
                        count++;
                    }
                }
            }
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.stop();
            time25 += timer.getElapsedTimeInMicroSec();
#endif

            double minv = vs[0], maxv = vs[0];

            for (int i = 1; i < 8; i++)
            {
                if (minv > vs[i])
                {
                    minv = vs[i];
                }
                if (maxv < vs[i])
                {
                    maxv = vs[i];
                }
            }
            if (minv > eps || maxv < -eps)
                return false;
            return true;

        } // test end
    }
    // eps is the interval [-eps,eps] we need to check
    // if [-eps,eps] overlap, return true
    // bbox_in_eps tell us if the box is totally in eps box
    bool evaluate_bbox_one_dimension_vector(
        std::array<double, 8> &t_up, std::array<double, 8> &t_dw,
        std::array<double, 8> &u_up, std::array<double, 8> &u_dw,
        std::array<double, 8> &v_up, std::array<double, 8> &v_dw,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const int dimension, const bool check_vf, const double eps, bool &bbox_in_eps)
    {
#ifdef TIGHT_INCLUSION_USE_TIMER
        igl::Timer timer;
#endif
        std::array<double, 8> vs;
        int count = 0;
        bbox_in_eps = false;
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.start();
#endif
        if (check_vf)
        {
            vs = function_vf(a0s[dimension], a1s[dimension], b0s[dimension], b1s[dimension], a0e[dimension], a1e[dimension], b0e[dimension], b1e[dimension], t_up, t_dw,
                             u_up, u_dw, v_up, v_dw);
        }
        else
        {
            vs = function_ee(a0s[dimension], a1s[dimension], b0s[dimension], b1s[dimension], a0e[dimension], a1e[dimension], b0e[dimension], b1e[dimension], t_up, t_dw,
                             u_up, u_dw, v_up, v_dw);
        }
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.stop();
        time25 += timer.getElapsedTimeInMicroSec();

#endif
        double minv = vs[0], maxv = vs[0];

        for (int i = 1; i < 8; i++)
        {
            if (minv > vs[i])
            {
                minv = vs[i];
            }
            if (maxv < vs[i])
            {
                maxv = vs[i];
            }
        }
        if (minv > eps || maxv < -eps)
            return false;
        if (minv >= -eps && maxv <= eps)
        {
            bbox_in_eps = true;
        }
        return true;
    }

    // ** this version can return the true x or y or z tolerance of the co-domain **
    // eps is the interval [-eps,eps] we need to check
    // if [-eps,eps] overlap, return true
    // bbox_in_eps tell us if the box is totally in eps box
    // ms is the minimum seperation
    bool evaluate_bbox_one_dimension_vector_return_tolerance(
        std::array<double, 8> &t_up, std::array<double, 8> &t_dw,
        std::array<double, 8> &u_up, std::array<double, 8> &u_dw,
        std::array<double, 8> &v_up, std::array<double, 8> &v_dw,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const int dimension, const bool check_vf, const double eps, const double ms, bool &bbox_in_eps,
        double &tol)
    {
#ifdef TIGHT_INCLUSION_USE_TIMER
        igl::Timer timer;
#endif
        std::array<double, 8> vs;
        int count = 0;
        bbox_in_eps = false;
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.start();
#endif
        if (check_vf)
        {
            vs = function_vf(a0s[dimension], a1s[dimension], b0s[dimension], b1s[dimension], a0e[dimension], a1e[dimension], b0e[dimension], b1e[dimension], t_up, t_dw,
                             u_up, u_dw, v_up, v_dw);
        }
        else
        {
            vs = function_ee(a0s[dimension], a1s[dimension], b0s[dimension], b1s[dimension], a0e[dimension], a1e[dimension], b0e[dimension], b1e[dimension], t_up, t_dw,
                             u_up, u_dw, v_up, v_dw);
        }
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.stop();
        time25 += timer.getElapsedTimeInMicroSec();

#endif
        double minv = vs[0], maxv = vs[0];

        for (int i = 1; i < 8; i++)
        {
            if (minv > vs[i])
            {
                minv = vs[i];
            }
            if (maxv < vs[i])
            {
                maxv = vs[i];
            }
        }
        tol = maxv - minv; // this is the real tolerance
        if (minv - ms > eps || maxv + ms < -eps)
            return false;
        if (minv + ms >= -eps && maxv - ms <= eps)
        {
            bbox_in_eps = true;
        }
        return true;
    }
    // the bounding boxes generated are t0, t1, u, u1, v0, v1 boxes
    template <typename T>
    void evaluate_tuv_bboxes(
        const std::array<Numccd, 2> &t,
        const std::array<Numccd, 2> &u,
        const std::array<Numccd, 2> &v,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        T tp, const bool check_vf, std::array<std::array<Eigen::Vector3d, 2>, 6> &bboxes)
    {

        int count = 0;
        std::array<Eigen::Vector3d, 8> pts;
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    if (!check_vf)
                    {
                        pts[count][0] =
                            function_f_ee(t[i], u[j], v[k], tp, 0, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);

                        pts[count][1] =
                            function_f_ee(t[i], u[j], v[k], tp, 1, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);

                        pts[count][2] =
                            function_f_ee(t[i], u[j], v[k], tp, 2, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
                    }
                    else
                    {
                        pts[count][0] =
                            function_f_vf(t[i], u[j], v[k], tp, 0, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);

                        pts[count][1] =
                            function_f_vf(t[i], u[j], v[k], tp, 1, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);

                        pts[count][2] =
                            function_f_vf(t[i], u[j], v[k], tp, 2, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
                    }
                    count++;
                }
            }
        }
        // now the parameters of pts are:
        // 000, 001, 010, 011, 100, 101, 110, 111
        std::array<Eigen::Vector3d, 4> bps;
        bps[0] = pts[0];
        bps[1] = pts[1];
        bps[2] = pts[2];
        bps[3] = pts[3];
        bboxes[0] = bbd_4_pts(bps); //t0

        bps[0] = pts[4];
        bps[1] = pts[5];
        bps[2] = pts[6];
        bps[3] = pts[7];
        bboxes[1] = bbd_4_pts(bps); //t1

        bps[0] = pts[0];
        bps[1] = pts[1];
        bps[2] = pts[4];
        bps[3] = pts[5];
        bboxes[2] = bbd_4_pts(bps); //u0

        bps[0] = pts[2];
        bps[1] = pts[3];
        bps[2] = pts[6];
        bps[3] = pts[7];
        bboxes[3] = bbd_4_pts(bps); //u1

        bps[0] = pts[0];
        bps[1] = pts[2];
        bps[2] = pts[4];
        bps[3] = pts[6];
        bboxes[4] = bbd_4_pts(bps); //v0

        bps[0] = pts[1];
        bps[1] = pts[3];
        bps[2] = pts[5];
        bps[3] = pts[7];
        bboxes[5] = bbd_4_pts(bps); //v1
    }

#ifdef TIGHT_INCLUSION_USE_GMP
    Vector3r function_f_ee_Rational(
        const Numccd &tpara, const Numccd &upara, const Numccd &vpara,
        const Eigen::Vector3d &a0sd,
        const Eigen::Vector3d &a1sd,
        const Eigen::Vector3d &b0sd,
        const Eigen::Vector3d &b1sd,
        const Eigen::Vector3d &a0ed,
        const Eigen::Vector3d &a1ed,
        const Eigen::Vector3d &b0ed,
        const Eigen::Vector3d &b1ed)
    {

        long tu = tpara.first;
        int td = tpara.second; // t=tu/(2^td)
        long uu = upara.first;
        int ud = upara.second;
        long vu = vpara.first;
        int vd = vpara.second;
        Vector3r
            a0s(a0sd[0], a0sd[1], a0sd[2]),
            a1s(a1sd[0], a1sd[1], a1sd[2]),
            b0s(b0sd[0], b0sd[1], b0sd[2]),
            b1s(b1sd[0], b1sd[1], b1sd[2]),
            a0e(a0ed[0], a0ed[1], a0ed[2]),
            a1e(a1ed[0], a1ed[1], a1ed[2]),
            b0e(b0ed[0], b0ed[1], b0ed[2]),
            b1e(b1ed[0], b1ed[1], b1ed[2]);
        Vector3r edge0_vertex0 = (a0e - a0s) * tu / power(1, td) + a0s;
        Vector3r edge0_vertex1 = (a1e - a1s) * tu / power(1, td) + a1s;
        Vector3r edge0_vertex = (edge0_vertex1 - edge0_vertex0) * uu / power(1, ud) + edge0_vertex0;

        Vector3r edge1_vertex0 = (b0e - b0s) * tu / power(1, td) + b0s;
        Vector3r edge1_vertex1 = (b1e - b1s) * tu / power(1, td) + b1s;
        Vector3r edge1_vertex = (edge1_vertex1 - edge1_vertex0) * vu / power(1, vd) + edge1_vertex0;

        return edge1_vertex - edge0_vertex;
    }
    Vector3r function_f_ee_Rational(
        const Rational &tpara, const Rational &upara, const Rational &vpara,
        const Eigen::Vector3d &a0sd,
        const Eigen::Vector3d &a1sd,
        const Eigen::Vector3d &b0sd,
        const Eigen::Vector3d &b1sd,
        const Eigen::Vector3d &a0ed,
        const Eigen::Vector3d &a1ed,
        const Eigen::Vector3d &b0ed,
        const Eigen::Vector3d &b1ed)
    {

        Vector3r
            a0s(a0sd[0], a0sd[1], a0sd[2]),
            a1s(a1sd[0], a1sd[1], a1sd[2]),
            b0s(b0sd[0], b0sd[1], b0sd[2]),
            b1s(b1sd[0], b1sd[1], b1sd[2]),
            a0e(a0ed[0], a0ed[1], a0ed[2]),
            a1e(a1ed[0], a1ed[1], a1ed[2]),
            b0e(b0ed[0], b0ed[1], b0ed[2]),
            b1e(b1ed[0], b1ed[1], b1ed[2]);

        // Vector3r las = (1-upara)*a0s+upara*a1s;
        // //std::cout<<"las, "<<las[0]<<", "<<las[1]<<", "<<las[2]<<std::endl;
        // Vector3r lae = (1-upara)*a0e+upara*a1e;
        // Vector3r lbs = (1-vpara)*b0s+vpara*b1s;
        // Vector3r lbe = (1-vpara)*b0e+vpara*b1e;
        // Vector3r lla=(1-tpara)*las + tpara*lae;
        // Vector3r llb=(1-tpara)*lbs + tpara*lbe;

        // //std::cout<<"lae, "<<lae[0]<<", "<<lae[1]<<", "<<lae[2]<<std::endl;
        // //std::cout<<"lbs, "<<lbs[0]<<", "<<lbs[1]<<", "<<lbs[2]<<std::endl;
        // //std::cout<<"lbe, "<<lbe[0]<<", "<<lbe[1]<<", "<<lbe[2]<<std::endl;
        // //std::cout<<"lla, "<<lla[0]<<", "<<lla[1]<<", "<<lla[2]<<std::endl;
        // //std::cout<<"llb, "<<llb[0]<<", "<<llb[1]<<", "<<llb[2]<<std::endl<<std::endl;
        // return lla-llb;

        Vector3r edge0_vertex0 = (a0e - a0s) * tpara + a0s;
        Vector3r edge0_vertex1 = (a1e - a1s) * tpara + a1s;
        Vector3r edge0_vertex = (edge0_vertex1 - edge0_vertex0) * upara + edge0_vertex0;

        Vector3r edge1_vertex0 = (b0e - b0s) * tpara + b0s;
        Vector3r edge1_vertex1 = (b1e - b1s) * tpara + b1s;
        Vector3r edge1_vertex = (edge1_vertex1 - edge1_vertex0) * vpara + edge1_vertex0;

        return edge1_vertex - edge0_vertex;
    }
    Vector3r function_f_vf_Rational(
        const Numccd &tpara, const Numccd &upara, const Numccd &vpara,
        const Eigen::Vector3d &a0sd,
        const Eigen::Vector3d &a1sd,
        const Eigen::Vector3d &b0sd,
        const Eigen::Vector3d &b1sd,
        const Eigen::Vector3d &a0ed,
        const Eigen::Vector3d &a1ed,
        const Eigen::Vector3d &b0ed,
        const Eigen::Vector3d &b1ed)
    {

        long tu = tpara.first;
        int td = tpara.second; // t=tu/(2^td)
        long uu = upara.first;
        int ud = upara.second;
        long vu = vpara.first;
        int vd = vpara.second;
        Vector3r
            vs(a0sd[0], a0sd[1], a0sd[2]),
            t0s(a1sd[0], a1sd[1], a1sd[2]),
            t1s(b0sd[0], b0sd[1], b0sd[2]),
            t2s(b1sd[0], b1sd[1], b1sd[2]),

            ve(a0ed[0], a0ed[1], a0ed[2]),
            t0e(a1ed[0], a1ed[1], a1ed[2]),
            t1e(b0ed[0], b0ed[1], b0ed[2]),
            t2e(b1ed[0], b1ed[1], b1ed[2]);

        Vector3r v = (ve - vs) * tu / power(1, td) + vs;

        Vector3r t0 = (t0e - t0s) * tu / power(1, td) + t0s;
        Vector3r t1 = (t1e - t1s) * tu / power(1, td) + t1s;
        Vector3r t2 = (t2e - t2s) * tu / power(1, td) + t2s;
        Vector3r p = (t1 - t0) * uu / power(1, ud) + (t2 - t0) * vu / power(1, vd) + t0;
        return v - p;
    }
    Vector3r function_f_vf_Rational(
        const Rational &tpara, const Rational &upara, const Rational &vpara,
        const Eigen::Vector3d &a0sd,
        const Eigen::Vector3d &a1sd,
        const Eigen::Vector3d &b0sd,
        const Eigen::Vector3d &b1sd,
        const Eigen::Vector3d &a0ed,
        const Eigen::Vector3d &a1ed,
        const Eigen::Vector3d &b0ed,
        const Eigen::Vector3d &b1ed)
    {

        Vector3r
            vs(a0sd[0], a0sd[1], a0sd[2]),
            t0s(a1sd[0], a1sd[1], a1sd[2]),
            t1s(b0sd[0], b0sd[1], b0sd[2]),
            t2s(b1sd[0], b1sd[1], b1sd[2]),

            ve(a0ed[0], a0ed[1], a0ed[2]),
            t0e(a1ed[0], a1ed[1], a1ed[2]),
            t1e(b0ed[0], b0ed[1], b0ed[2]),
            t2e(b1ed[0], b1ed[1], b1ed[2]);

        Vector3r v = (ve - vs) * tpara + vs;

        Vector3r t0 = (t0e - t0s) * tpara + t0s;
        Vector3r t1 = (t1e - t1s) * tpara + t1s;
        Vector3r t2 = (t2e - t2s) * tpara + t2s;
        Vector3r p = (t1 - t0) * upara + (t2 - t0) * vpara + t0;
        return v - p;
    }
#endif

    bool Origin_in_function_bounding_box_double(
        const Interval3 &paras,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const bool check_vf,
        const std::array<double, 3> &box)
    {
#ifdef TIGHT_INCLUSION_USE_TIMER
        igl::Timer timer;
#endif
        std::array<Numccd, 2> t, u, v;

        t[0] = paras[0].first;
        t[1] = paras[0].second;
        u[0] = paras[1].first;
        u[1] = paras[1].second;
        v[0] = paras[2].first;
        v[1] = paras[2].second;
        //bool zero_0=false, zer0_1=false, zero_2=false;
        double input_type;
        bool ck;
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.start();
#endif
        ck = evaluate_bbox_one_dimension(t, u, v, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, 0, input_type, check_vf, box[0]);
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.stop();
        time23 += timer.getElapsedTimeInMicroSec();
#endif
        if (!ck)
            return false;
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.start();
#endif
        ck = evaluate_bbox_one_dimension(t, u, v, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, 1, input_type, check_vf, box[1]);
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.stop();
        time23 += timer.getElapsedTimeInMicroSec();
#endif
        if (!ck)
            return false;
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.start();
#endif
        ck = evaluate_bbox_one_dimension(t, u, v, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, 2, input_type, check_vf, box[2]);
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.stop();
        time23 += timer.getElapsedTimeInMicroSec();
#endif
        if (!ck)
            return false;
        return true;
    }

    // give the result of if the hex overlaps the input box around Origin
    // use vectorized hex-vertex-solving function for acceleration
    // box_in_eps shows if this hex is totally inside box. if so, no need to do further bisection
    bool Origin_in_function_bounding_box_double_vector(
        const Interval3 &paras,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const bool check_vf,
        const std::array<double, 3> &box, bool &box_in_eps)
    {
#ifdef TIGHT_INCLUSION_USE_TIMER
        igl::Timer timer;
#endif
        box_in_eps = false;
        std::array<double, 8> t_up;
        std::array<double, 8> t_dw;
        std::array<double, 8> u_up;
        std::array<double, 8> u_dw;
        std::array<double, 8> v_up;
        std::array<double, 8> v_dw;
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.start();
#endif
        convert_tuv_to_array(paras, t_up, t_dw, u_up, u_dw, v_up, v_dw);
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.stop();
        time24 += timer.getElapsedTimeInMicroSec();
#endif
        bool ck;
        bool box_in[3];
        for (int i = 0; i < 3; i++)
        {
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.start();
#endif
            ck = evaluate_bbox_one_dimension_vector(t_up, t_dw, u_up, u_dw, v_up, v_dw,
                                                    a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, i, check_vf, box[i], box_in[i]);
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.stop();
            time23 += timer.getElapsedTimeInMicroSec();
#endif
            if (!ck)
                return false;
        }
        if (box_in[0] && box_in[1] && box_in[2])
        {
            box_in_eps = true;
        }
        return true;
    }

    // ** this version can return the true tolerance of the co-domain **
    // give the result of if the hex overlaps the input box around Origin
    // use vectorized hex-vertex-solving function for acceleration
    // box_in_eps shows if this hex is totally inside box. if so, no need to do further bisection
    bool Origin_in_function_bounding_box_double_vector_return_tolerance(
        const Interval3 &paras,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const bool check_vf,
        const std::array<double, 3> &box, const double ms, bool &box_in_eps,
        std::array<double, 3> &tolerance)
    {
#ifdef TIGHT_INCLUSION_USE_TIMER
        igl::Timer timer;
#endif
        box_in_eps = false;
        std::array<double, 8> t_up;
        std::array<double, 8> t_dw;
        std::array<double, 8> u_up;
        std::array<double, 8> u_dw;
        std::array<double, 8> v_up;
        std::array<double, 8> v_dw;
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.start();
#endif
        convert_tuv_to_array(paras, t_up, t_dw, u_up, u_dw, v_up, v_dw);
#ifdef TIGHT_INCLUSION_USE_TIMER
        timer.stop();
        time24 += timer.getElapsedTimeInMicroSec();
#endif
        bool ck;
        bool box_in[3];
        for (int i = 0; i < 3; i++)
        {
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.start();
#endif
            ck = evaluate_bbox_one_dimension_vector_return_tolerance(t_up, t_dw, u_up, u_dw, v_up, v_dw,
                                                                     a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, i, check_vf, box[i], ms, box_in[i], tolerance[i]);
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.stop();
            time23 += timer.getElapsedTimeInMicroSec();
#endif
            if (!ck)
                return false;
        }
        if (box_in[0] && box_in[1] && box_in[2])
        {
            box_in_eps = true;
        }
        return true;
    }
    bool bounding_box_intersection(const Eigen::Vector3d &pmin, const Eigen::Vector3d &pmax,
                                   const Eigen::Vector3d &qmin, const Eigen::Vector3d &qmax)
    {
        if (pmax[0] < qmin[0] || pmax[1] < qmin[1] || pmax[2] < qmin[2])
        {
            return false;
        }
        if (qmax[0] < pmin[0] || qmax[1] < pmin[1] || qmax[2] < pmin[2])
        {
            return false;
        }
        return true;
    }
    bool estimate_tuv_through_bbox(
        const Interval3 &paras,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const bool check_vf,
        const std::array<double, 3> &box)
    {
        //igl::Timer timer;
        std::array<Numccd, 2> t, u, v;

        t[0] = paras[0].first;
        t[1] = paras[0].second;
        u[0] = paras[1].first;
        u[1] = paras[1].second;
        v[0] = paras[2].first;
        v[1] = paras[2].second;
        //bool zero_0=false, zer0_1=false, zero_2=false;
        double input_type;
        std::array<std::array<Eigen::Vector3d, 2>, 6> bboxes;
        // get the bounding boxes
        evaluate_tuv_bboxes(t, u, v, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, input_type, check_vf, bboxes);
        //TODO
        if (!evaluate_bbox_one_dimension(t, u, v, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, 0, input_type, check_vf, box[0]))
            return false;
        if (!evaluate_bbox_one_dimension(t, u, v, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, 1, input_type, check_vf, box[1]))
            return false;
        if (!evaluate_bbox_one_dimension(t, u, v, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, 2, input_type, check_vf, box[2]))
            return false;
        return true;
    }

#ifdef TIGHT_INCLUSION_USE_GMP
    bool Origin_in_function_bounding_box_Rational(
        const Interval3 &paras,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e, const bool check_vf)
    {
        //igl::Timer timer;
        std::array<Numccd, 2> t, u, v;

        t[0] = paras[0].first;
        t[1] = paras[0].second;
        u[0] = paras[1].first;
        u[1] = paras[1].second;
        v[0] = paras[2].first;
        v[1] = paras[2].second;
        //std::cout<<"t, ["<<t[0].first/pow(2,t[0].second)<<","<<t[1].first/pow(2,t[1].second)<<"]"<<std::endl;
        //std::cout<<"u, ["<<u[0].first/pow(2,u[0].second)<<","<<u[1].first/pow(2,u[1].second)<<"]"<<std::endl;
        //std::cout<<"v, ["<<v[0].first/pow(2,v[0].second)<<","<<v[1].first/pow(2,v[1].second)<<"]"<<std::endl<<std::endl;
        Vector3r minv, maxv;
        std::array<Vector3r, 8> pts;
        int c = 0;
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    //std::cout<<"c="<<c<<std::endl;
                    if (!check_vf)
                    {
                        //std::cout<<"ee"<<std::endl;
                        pts[c] = function_f_ee_Rational(t[i], u[j], v[k], a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
                    }
                    else
                    {
                        pts[c] = function_f_vf_Rational(t[i], u[j], v[k], a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
                    }

                    c++;
                }
            }
        }
        minv = pts[0];
        maxv = pts[0];
        for (int i = 0; i < 8; i++)
        {
            if (minv[0] > pts[i][0])
            {
                minv[0] = pts[i][0];
            }
            if (minv[1] > pts[i][1])
            {
                minv[1] = pts[i][1];
            }
            if (minv[2] > pts[i][2])
            {
                minv[2] = pts[i][2];
            }
            if (maxv[0] < pts[i][0])
            {
                maxv[0] = pts[i][0];
            }
            if (maxv[1] < pts[i][1])
            {
                maxv[1] = pts[i][1];
            }
            if (maxv[2] < pts[i][2])
            {
                maxv[2] = pts[i][2];
            }
        }
        if (minv[0] <= 0 && minv[1] <= 0 && minv[2] <= 0)
        {
            if (maxv[0] >= 0 && maxv[1] >= 0 && maxv[2] >= 0)
            {
                return true;
            }
        }
        return false;
    }

    // can be used for minumum seperation. in that case, box will be ms-box instead of ms+err-box
    bool Origin_in_function_bounding_box_Rational_return_tolerance(
        const Interval3 &paras,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const bool check_vf,
        const std::array<double, 3> &box, bool &box_in_eps,
        std::array<double, 3> &tolerance)
    {
        //igl::Timer timer;
        std::array<Numccd, 2> t, u, v;
        box_in_eps = true;
        t[0] = paras[0].first;
        t[1] = paras[0].second;
        u[0] = paras[1].first;
        u[1] = paras[1].second;
        v[0] = paras[2].first;
        v[1] = paras[2].second;
        //std::cout<<"t, ["<<t[0].first/pow(2,t[0].second)<<","<<t[1].first/pow(2,t[1].second)<<"]"<<std::endl;
        //std::cout<<"u, ["<<u[0].first/pow(2,u[0].second)<<","<<u[1].first/pow(2,u[1].second)<<"]"<<std::endl;
        //std::cout<<"v, ["<<v[0].first/pow(2,v[0].second)<<","<<v[1].first/pow(2,v[1].second)<<"]"<<std::endl<<std::endl;
        Vector3r minv, maxv;
        std::array<Vector3r, 8> pts;
        int c = 0;
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    //std::cout<<"c="<<c<<std::endl;
                    if (!check_vf)
                    {
                        //std::cout<<"ee"<<std::endl;
                        pts[c] = function_f_ee_Rational(t[i], u[j], v[k], a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
                    }
                    else
                    {
                        pts[c] = function_f_vf_Rational(t[i], u[j], v[k], a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
                    }

                    c++;
                }
            }
        }
        minv = pts[0];
        maxv = pts[0];
        for (int i = 0; i < 8; i++)
        {
            if (minv[0] > pts[i][0])
            {
                minv[0] = pts[i][0];
            }
            if (minv[1] > pts[i][1])
            {
                minv[1] = pts[i][1];
            }
            if (minv[2] > pts[i][2])
            {
                minv[2] = pts[i][2];
            }
            if (maxv[0] < pts[i][0])
            {
                maxv[0] = pts[i][0];
            }
            if (maxv[1] < pts[i][1])
            {
                maxv[1] = pts[i][1];
            }
            if (maxv[2] < pts[i][2])
            {
                maxv[2] = pts[i][2];
            }
        }
        for (int i = 0; i < 3; i++)
        {
            if (maxv[i] > box[i] || minv[i] < -box[i])
            {
                box_in_eps = false;
            }
            tolerance[i] = (maxv[i] - minv[i]).to_double();
        }
        for (int i = 0; i < 3; i++)
        {
            if (maxv[i] < -box[i] || minv[i] > box[i])
            {
                return false;
            }
        }
        return true;
    }

    bool Origin_in_function_bounding_box_Rational(
        const std::array<std::pair<Rational, Rational>, 3> &paras,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e, const bool check_vf)
    {
        //igl::Timer timer;
        std::array<Rational, 2> t, u, v;

        t[0] = paras[0].first;
        t[1] = paras[0].second;
        u[0] = paras[1].first;
        u[1] = paras[1].second;
        v[0] = paras[2].first;
        v[1] = paras[2].second;
        Vector3r minv, maxv;
        std::array<Vector3r, 8> pts;
        int c = 0;
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {

                    if (!check_vf)
                    {
                        pts[c] = function_f_ee_Rational(t[i], u[j], v[k], a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
                    }
                    else
                    {
                        pts[c] = function_f_vf_Rational(t[i], u[j], v[k], a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e);
                    }

                    c++;
                }
            }
        }
        minv = pts[0];
        maxv = pts[0];
        for (int i = 0; i < 8; i++)
        {
            if (minv[0] > pts[i][0])
            {
                minv[0] = pts[i][0];
            }
            if (minv[1] > pts[i][1])
            {
                minv[1] = pts[i][1];
            }
            if (minv[2] > pts[i][2])
            {
                minv[2] = pts[i][2];
            }
            if (maxv[0] < pts[i][0])
            {
                maxv[0] = pts[i][0];
            }
            if (maxv[1] < pts[i][1])
            {
                maxv[1] = pts[i][1];
            }
            if (maxv[2] < pts[i][2])
            {
                maxv[2] = pts[i][2];
            }
        }
        if (minv[0] <= 0 && minv[1] <= 0 && minv[2] <= 0)
        {
            if (maxv[0] >= 0 && maxv[1] >= 0 && maxv[2] >= 0)
            {
                return true;
            }
        }
        return false;
    }
    
#endif
// return power t. n=result*2^t
    long reduction(const long n, long &result)
    {
        int t = 0;
        int newn = n;
        while (newn % 2 == 0)
        {
            newn = newn / 2;
            t++;
        }
        result = newn;
        return t;
    }

    std::pair<Singleinterval, Singleinterval> bisect(const Singleinterval &inter)
    {
        Numccd low = inter.first;
        Numccd up = inter.second;

        // interval is [k1/pow(2,n1), k2/pow(2,n2)], k1,k2,n1,n2 are all not negative
        long k1 = low.first;
        int n1 = low.second;
        long k2 = up.first;
        int n2 = up.second;

        assert(k1 >= 0 && k2 >= 0 && n1 >= 0 && n2 >= 0);

        std::pair<Singleinterval, Singleinterval> result;
        long k;
        int n;
        int p;
        long r;
        if (n2 == n1)
        {
            p = reduction(k1 + k2, r);
            k = r;
            n = n2 - p + 1;
        }
        if (n2 > n1)
        {
            k = k1 * power(1, n2 - n1) + k2;
            assert(k % 2 == 1);
            n = n2 + 1;
        }
        if (n2 < n1)
        {
            k = k1 + k2 * power(1, n1 - n2);
            assert(k % 2 == 1);
            n = n1 + 1;
        }
        Numccd newnum(k, n);
        Singleinterval i1(low, newnum), i2(newnum, up);
        // std::cout<<"low,"<<Numccd2double(low)<<",up,"<<Numccd2double(up)<<", num, "<<Numccd2double(newnum)<<std::endl;
        // std::cout<<"new, k1, "<<newnum.first<<", n1, "<<newnum.second<<std::endl;
        assert(Numccd2double(newnum) > Numccd2double(low) && Numccd2double(newnum) < Numccd2double(up));
        result.first = i1;
        result.second = i2;
        return result;
    }

#ifdef TIGHT_INCLUSION_USE_GMP
    std::pair<std::pair<Rational, Rational>, std::pair<Rational, Rational>> bisect(const std::pair<Rational, Rational> &inter)
    {
        std::pair<std::pair<Rational, Rational>, std::pair<Rational, Rational>> result;
        std::pair<Rational, Rational> single;
        Rational mid = (inter.first + inter.second) / 2;
        single.first = inter.first;
        single.second = mid;
        result.first = single;
        single.first = mid;
        single.second = inter.second;
        result.second = single;
        return result;
    }

#endif

    bool sum_no_larger_1(const Numccd &num1, const Numccd &num2)
    {
        long k1 = num1.first;
        int n1 = num1.second;
        long k2 = num2.first;
        int n2 = num2.second;
        long k;
        int n;
        if (n1 == n2)
        {
            k = k1 + k2;
            n = n1;
        }
        if (n1 < n2)
        {
            k = power(1, n2 - n1) * k1 + k2;
            n = n2;
        }
        if (n1 > n2)
        {
            k = power(1, n1 - n2) * k2 + k1;
            n = n1;
        }
        assert(k >= 0 && n >= 0);
        if (k > power(1, n))
            return false;
        else
            return true;
    }
    //check if num1<num2
    bool less_than(const Numccd &num1, const Numccd &num2)
    {
        long k1 = num1.first;
        int n1 = num1.second;
        long k2 = num2.first;
        int n2 = num2.second;

        if (n1 < n2)
        {
            k1 = power(1, n2 - n1) * k1;
        }
        if (n1 > n2)
        {
            k2 = power(1, n1 - n2) * k2;
        }
        if (k1 < k2)
            return true;
        return false;
    }
    bool less_than(const Numccd &num1, const double num2)
    {
        double res = Numccd2double(num1);
        if (res < num2)
        {
            return true;
        }
        return false;
    }
    bool interval_overlap_region(const Singleinterval &itv, const double r1, const double r2)
    {
        double b1 = Numccd2double(itv.first);
        double b2 = Numccd2double(itv.second);
        if (b2 < r1 || b1 > r2)
            return false;
        return true;
    }

#ifdef TIGHT_INCLUSION_USE_GMP
    bool sum_no_larger_1_Rational(const Numccd &num1, const Numccd &num2)
    {
        long k1 = num1.first;
        int n1 = num1.second;
        long k2 = num2.first;
        int n2 = num2.second;
        Rational nbr1, nbr2;
        nbr1 = Rational(k1) / Rational(power(1, n1));
        nbr2 = Rational(k2) / Rational(power(1, n2));
        Rational rst = nbr1 + nbr2 - Rational(1);
        if (rst > 0)
            return false;
        else
            return true;
    }
#endif

    // calculate the sign of f. dim is the dimension we are evaluating.
    template <typename T>
    T function_f_ee(
        const Numccd &tpara, const Numccd &upara, const Numccd &vpara, const T &type, const int dim,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e)
    {

        long tu = tpara.first;
        int td = tpara.second; // t=tu/(2^td)
        long uu = upara.first;
        int ud = upara.second;
        long vu = vpara.first;
        int vd = vpara.second;

        T edge0_vertex0 = (T(a0e[dim]) - T(a0s[dim])) * tu / power(1, td) + T(a0s[dim]);
        T edge0_vertex1 = (T(a1e[dim]) - T(a1s[dim])) * tu / power(1, td) + T(a1s[dim]);
        T edge0_vertex = (edge0_vertex1 - edge0_vertex0) * uu / power(1, ud) + edge0_vertex0;

        T edge1_vertex0 = (T(b0e[dim]) - T(b0s[dim])) * tu / power(1, td) + T(b0s[dim]);
        T edge1_vertex1 = (T(b1e[dim]) - T(b1s[dim])) * tu / power(1, td) + T(b1s[dim]);
        T edge1_vertex = (edge1_vertex1 - edge1_vertex0) * vu / power(1, vd) + edge1_vertex0;

        return edge1_vertex - edge0_vertex;
    }

    template <typename T>
    T function_f_vf(
        const Numccd &tpara, const Numccd &upara, const Numccd &vpara, const T &type, const int dim,
        const Eigen::Vector3d &vs,
        const Eigen::Vector3d &t0s,
        const Eigen::Vector3d &t1s,
        const Eigen::Vector3d &t2s,

        const Eigen::Vector3d &ve,
        const Eigen::Vector3d &t0e,
        const Eigen::Vector3d &t1e,
        const Eigen::Vector3d &t2e)
    {

        long tu = tpara.first;
        int td = tpara.second; // t=tu/(2^td)
        long uu = upara.first;
        int ud = upara.second;
        long vu = vpara.first;
        int vd = vpara.second;

        T v = (T(ve[dim]) - T(vs[dim])) * tu / power(1, td) + T(vs[dim]);
        T t0 = (T(t0e[dim]) - T(t0s[dim])) * tu / power(1, td) + T(t0s[dim]);
        T t1 = (T(t1e[dim]) - T(t1s[dim])) * tu / power(1, td) + T(t1s[dim]);
        T t2 = (T(t2e[dim]) - T(t2s[dim])) * tu / power(1, td) + T(t2s[dim]);
        T p = (t1 - t0) * uu / power(1, ud) + (t2 - t0) * vu / power(1, vd) + t0;
        return v - p;
    }

    // this version cannot give the impact time at t=1, although this collision can
    // be detected at t=0 of the next time step, but still may cause problems in
    // line-search based physical simulation
    bool interval_root_finder_double_normalCCD(
        const Eigen::VectorX3d &tol,
        double &toi,
        const bool check_vf,
        const std::array<double, 3> err,
        const double ms,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e)
    {
        auto cmp = [](std::pair<Interval3, int> i1, std::pair<Interval3, int> i2) {
            return !(less_than(i1.first[0].first, i2.first[0].first));
        };
        Numccd low_number;
        low_number.first = 0;
        low_number.second = 0; // low_number=0;
        Numccd up_number;
        up_number.first = 1;
        up_number.second = 0; // up_number=1;
        // initial interval [0,1]
        Singleinterval init_interval;
        init_interval.first = low_number;
        init_interval.second = up_number;
        //build interval set [0,1]x[0,1]x[0,1]
        Interval3 iset;
        iset[0] = init_interval;
        iset[1] = init_interval;
        iset[2] = init_interval;
        // Stack of intervals and the last split dimension
        // std::stack<std::pair<Interval3,int>> istack;
        std::priority_queue<std::pair<Interval3, int>, std::vector<std::pair<Interval3, int>>, decltype(cmp)> istack(cmp);
        istack.emplace(iset, -1);

        // current intervals
        Interval3 current;
        std::array<double, 3> err_and_ms;
        err_and_ms[0] = err[0] + ms;
        err_and_ms[1] = err[1] + ms;
        err_and_ms[2] = err[2] + ms;
        refine = 0;
        toi = std::numeric_limits<double>::infinity();
        Numccd TOI;
        TOI.first = 1;
        TOI.second = 0;
        //std::array<double,3>
        bool collision = false;
        int rnbr = 0;
        while (!istack.empty())
        {
            current = istack.top().first;
            int last_split = istack.top().second;
            istack.pop();
            // if(rnbr>0&&less_than( current[0].first,TOI)){
            //     std::cout<<"not the first"<<std::endl;
            //     // continue;
            // }
            if (!less_than(current[0].first, TOI))
            {
                // std::cout<<"not the first"<<std::endl;
                continue;
            }
            //TOI should always be no larger than current

            // if(Numccd2double(current[0].first)>=Numccd2double(TOI)){
            //     std::cout<<"here wrong, comparing"<<std::endl;
            // }
#ifdef TIGHT_INCLUSION_USE_TIMER
            igl::Timer timer;

            timer.start();
#endif
            refine++;
            bool zero_in;
            bool box_in;
            zero_in = Origin_in_function_bounding_box_double_vector(current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, check_vf, err_and_ms, box_in);

#ifdef TIGHT_INCLUSION_USE_TIMER

            timer.stop();
            time20 += timer.getElapsedTimeInMicroSec();
#endif
#ifdef COMPARE_WITH_RATIONAL // this is defined in the begining of this file

            zero_in = Origin_in_function_bounding_box_Rational(current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, check_vf);

#endif
            if (!zero_in)
                continue;
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.start();
#endif
            Eigen::VectorX3d widths = width(current);
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.stop();
            time21 += timer.getElapsedTimeInMicroSec();
#endif
            if ((widths.array() <= tol.array()).all())
            {
                TOI = current[0].first;
                collision = true;
                rnbr++;
                // continue;
                toi = Numccd2double(TOI);
                return true;
            }
            if (box_in)
            {
                TOI = current[0].first;
                collision = true;
                rnbr++;
                // continue;
                toi = Numccd2double(TOI);
                return true;
            }

            std::array<bool, 3> check;
            Eigen::VectorX3d widthratio;
            widthratio.resize(3);
            check[0] = false;
            check[1] = false;
            check[2] = false;
            for (int i = 0; i < 3; i++)
            {
                widthratio(i) = widths(i) / tol(i);
                if (widths(i) > tol(i))
                    check[i] = true; // means this need to be checked
            }

            int split_i = -1;
            for (int i = 0; i < 3; i++)
            {
                if (check[i])
                {
                    if (check[(i + 1) % 3] && check[(i + 2) % 3])
                    {
                        if (widthratio(i) >= widthratio((i + 1) % 3) && widthratio(i) >= widthratio((i + 2) % 3))
                        {
                            split_i = i;
                            break;
                        }
                    }
                    if (check[(i + 1) % 3] && !check[(i + 2) % 3])
                    {
                        if (widthratio(i) >= widthratio((i + 1) % 3))
                        {
                            split_i = i;
                            break;
                        }
                    }
                    if (!check[(i + 1) % 3] && check[(i + 2) % 3])
                    {
                        if (widthratio(i) >= widthratio((i + 2) % 3))
                        {
                            split_i = i;
                            break;
                        }
                    }
                    if (!check[(i + 1) % 3] && !check[(i + 2) % 3])
                    {

                        split_i = i;
                        break;
                    }
                }
            }
            if (split_i < 0)
            {
                std::cout << "ERROR OCCURRED HERE, DID NOT FIND THE RIGHT DIMENSION TO SPLIT" << std::endl;
            }
            // Bisect the next dimension that is greater than its tolerance
            // int split_i;
            // for (int i = 1; i <= 3; i++) {
            //     split_i = (last_split + i) % 3;
            //     if (widths(split_i) > tol(split_i)) {
            //         break;
            //     }
            // }
            std::pair<Singleinterval, Singleinterval> halves = bisect(current[split_i]);
            if (!less_than(halves.first.first, halves.first.second))
            {
                std::cout << "OVERFLOW HAPPENS WHEN SPLITTING INTERVALS" << std::endl;
                return true;
            }
            if (!less_than(halves.second.first, halves.second.second))
            {
                std::cout << "OVERFLOW HAPPENS WHEN SPLITTING INTERVALS" << std::endl;
                return true;
            }
            if (check_vf)
            {
                //std::cout<<"*** check_vf"<<std::endl;
                if (split_i == 1)
                {
                    // assert(sum_no_larger_1(halves.first.first, current[2].first)==sum_no_larger_1_Rational(halves.first.first, current[2].first));
                    // assert(sum_no_larger_1(halves.second.first, current[2].first)==sum_no_larger_1_Rational(halves.second.first, current[2].first));

                    if (sum_no_larger_1(halves.second.first, current[2].first))
                    {
                        current[split_i] = halves.second;
                        istack.emplace(current, split_i);
                    }
                    if (sum_no_larger_1(halves.first.first, current[2].first))
                    {
                        current[split_i] = halves.first;
                        istack.emplace(current, split_i);
                    }
                }

                if (split_i == 2)
                {
                    //assert(sum_no_larger_1(halves.first.first, current[1].first)==sum_no_larger_1_Rational(halves.first.first, current[1].first));
                    //assert(sum_no_larger_1(halves.second.first, current[1].first)==sum_no_larger_1_Rational(halves.second.first, current[1].first));

                    if (sum_no_larger_1(halves.second.first, current[1].first))
                    {
                        current[split_i] = halves.second;
                        istack.emplace(current, split_i);
                    }
                    if (sum_no_larger_1(halves.first.first, current[1].first))
                    {
                        current[split_i] = halves.first;
                        istack.emplace(current, split_i);
                    }
                }
                if (split_i == 0)
                {
                    current[split_i] = halves.second;
                    istack.emplace(current, split_i);
                    current[split_i] = halves.first;
                    istack.emplace(current, split_i);
                }
            }
            else
            {
                current[split_i] = halves.second;
                istack.emplace(current, split_i);
                current[split_i] = halves.first;
                istack.emplace(current, split_i);
            }
        }
        if (collision)
            toi = Numccd2double(TOI);
        // if(toi==0)std::cout<<"infinate roots, "<<std::endl;
        // if(rnbr>5) std::cout<<"nbr of roots, "<<rnbr<<", time, "<<toi<<std::endl;
        return collision;
        return false;
    }

    void sum_up(const Numccd &nbr1, const Numccd &nbr2, Numccd &result)
    {
        long k1 = nbr1.first;
        long k2 = nbr2.first;
        int n1 = nbr1.second;
        int n2 = nbr2.second;
        long k;
        int n;
        int p;
        long r;
        if (n2 == n1)
        {
            p = reduction(k1 + k2, r);
            k = r;
            n = n2 - p;
        }
        if (n2 > n1)
        {
            k = k1 * power(1, n2 - n1) + k2;
            assert(k % 2 == 1);
            n = n2;
        }
        if (n2 < n1)
        {
            k = k1 + k2 * power(1, n1 - n2);
            assert(k % 2 == 1);
            n = n1;
        }
        result.first = k;
        result.second = n;
        return;
    }
    // find a t value that t<tol TODO need to be fixed
    void t_tol_width(Numccd &x, const double tol)
    {
        x.first = 1;
        x.second = 0;
        while (Numccd2double(x) >= tol)
        {
            x.second += 1;
        }
        while (Numccd2double(x) < tol)
        {
            x.first += 1;
        }
    }

    // when check_t_overlap = false, check [0,1]x[0,1]x[0,1]; otherwise, check [0, t_max]x[0,1]x[0,1]
    bool interval_root_finder_double_horizontal_tree(
        const Eigen::VectorX3d &tol,
        const double co_domain_tolerance,
        const Interval3 &iset,
        const bool check_t_overlap,
        const double max_t, // check interval [0, max_t] when check_t_overlap is set as TRUE
        double &toi,
        const bool check_vf,
        const std::array<double, 3> err,
        const double ms,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const int max_itr,
        double &output_tolerance)
    {

        // if max_itr <0, output_tolerance= co_domain_tolerance;
        // else, output_tolearancewill be the precision after iteration time > max_itr
        output_tolerance = co_domain_tolerance;

        // this is used to catch the tolerance for each level
        double temp_output_tolerance = co_domain_tolerance;
        //return time1 >= time2
        auto time_cmp = [](std::pair<Interval3, int> i1, std::pair<Interval3, int> i2) {
            return !(less_than(i1.first[0].first, i2.first[0].first));
        };

        // check the tree level by level instead of going deep
        //(if level 1 != level 2, return level 1 >= level 2;
        // else, return time1 >= time2)
        auto horiz_cmp = [](std::pair<Interval3, int> i1, std::pair<Interval3, int> i2) {
            if (i1.second != i2.second)
            {
                {
                    return (i1.second >= i2.second);
                }
            }
            else
            {
                return !(less_than(i1.first[0].first, i2.first[0].first));
            }
        };

        // Stack of intervals and the last split dimension
        // std::stack<std::pair<Interval3,int>> istack;
        auto cmp = horiz_cmp;
        std::priority_queue<std::pair<Interval3, int>, std::vector<std::pair<Interval3, int>>, decltype(cmp)> istack(cmp);
        istack.emplace(iset, -1);

        // current intervals
        Interval3 current;
        refine = 0;
        double impact_ratio;

        impact_ratio = 1;

        toi = std::numeric_limits<double>::infinity(); //set toi as infinate
        // temp_toi is to catch the toi of each level
        double temp_toi = toi;
        Numccd TOI;
        TOI.first = 4;
        TOI.second = 0;        // set TOI as 4. this is to record the impact time of this level
        Numccd TOI_SKIP = TOI; // this is to record the element that already small enough or contained in eps-box
        bool use_skip = false; // this is to record if TOI_SKIP is used.
        //std::array<double,3>
        bool collision = false;
        int rnbr = 0;
        int current_level = -2; // in the begining, current_level != level
        int box_in_level = -2;  // this checks if all the boxes before this
        // level < tolerance. only true, we can return when we find one overlaps eps box and smaller than tolerance or eps-box
        bool this_level_less_tol = true;
        bool find_level_root = false;
        // double current_tolerance=std::numeric_limits<double>::infinity(); // set returned tolerance as infinite
        double t_upper_bound = max_t; // 2*tol make it more conservative
        while (!istack.empty())
        {
            current = istack.top().first;
            int level = istack.top().second;
            istack.pop();

            // if this box is later than TOI_SKIP in time, we can skip this one.
            // TOI_SKIP is only updated when the box is small enough or totally contained in eps-box
            if (!less_than(current[0].first, TOI_SKIP))
            {
                continue;
            }
            if (box_in_level != level)
            { // before check a new level, set this_level_less_tol=true
                box_in_level = level;
                this_level_less_tol = true;
            }
#ifdef TIGHT_INCLUSION_USE_TIMER
            igl::Timer timer;

            timer.start();
#endif
            refine++;
            bool zero_in;
            bool box_in;
            std::array<double, 3> true_tol;
#ifdef COMPARE_WITH_RATIONAL // this is defined in the begining of this file
            std::array<double, 3> ms_3d;
            ms_3d[0] = ms;
            ms_3d[1] = ms;
            ms_3d[2] = ms;
            zero_in = Origin_in_function_bounding_box_Rational_return_tolerance(current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, check_vf, err_and_ms, box_in, true_tol);

#else
            zero_in = Origin_in_function_bounding_box_double_vector_return_tolerance(current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, check_vf, err, ms, box_in, true_tol);
#endif
#ifdef TIGHT_INCLUSION_USE_TIMER

            timer.stop();
            time20 += timer.getElapsedTimeInMicroSec();
#endif

            if (!zero_in)
                continue;
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.start();
#endif
            Eigen::VectorX3d widths = width(current);
#ifdef TIGHT_INCLUSION_USE_TIMER
            timer.stop();
            time21 += timer.getElapsedTimeInMicroSec();
#endif

            bool tol_condition = true_tol[0] <= co_domain_tolerance && true_tol[1] <= co_domain_tolerance && true_tol[2] <= co_domain_tolerance;

            // Condition 1, stopping condition on t, u and v is satisfied. this is useless now since we have condition 2
            bool condition1 = (widths.array() <= tol.array()).all();

            // Condition 2, zero_in = true, box inside eps-box and in this level,
            // no box whose zero_in is true but box size larger than tolerance, can return
            bool condition2 = box_in && this_level_less_tol;
            if (!tol_condition)
            {
                this_level_less_tol = false;
                // this level has at least one box whose size > tolerance, thus we
                // cannot directly return if find one box whose size < tolerance or box-in
            }

            // Condition 3, in this level, we find a box that zero-in and size < tolerance.
            // and no other boxes whose zero-in is true in this level before this one is larger than tolerance, can return
            bool condition3 = this_level_less_tol;
            if (condition1 || condition2 || condition3)
            {
                TOI = current[0].first;
                collision = true;
                rnbr++;
                // continue;
                toi = Numccd2double(TOI) * impact_ratio;
                return true;
                // we don't need to compare with TOI_SKIP because we already continue
                // when t>=TOI_SKIP
            }

            if (max_itr > 0)
            { // if max_itr < 0, then stop until stack empty
                if (current_level != level)
                {
                    // output_tolerance=current_tolerance;
                    // current_tolerance=0;
                    current_level = level;
                    find_level_root = false;
                }
                // current_tolerance=std::max(
                // std::max(std::max(current_tolerance,true_tol[0]),true_tol[1]),true_tol[2]
                // );
                if (!find_level_root)
                {
                    TOI = current[0].first;
                    // collision=true;
                    rnbr++;
                    // continue;
                    temp_toi = Numccd2double(TOI) * impact_ratio;

                    // if the real tolerance is larger than input, use the real one;
                    // if the real tolerance is smaller than input, use input
                    temp_output_tolerance = std::max(std::max(std::max(true_tol[0], true_tol[1]), true_tol[2]), co_domain_tolerance);
                    find_level_root = true; // this ensures always find the earlist root
                }
                if (refine > max_itr)
                {
                    toi = temp_toi;
                    output_tolerance = temp_output_tolerance;
                    refine_return++;
                    // std::cout<<"return from refine"<<std::endl;
                    return true;
                }
                // get the time of impact down here
            }

            // if this box is small enough, or inside of eps-box, then just continue,
            // but we need to record the collision time
            if (tol_condition || box_in)
            {
                if (less_than(current[0].first, TOI_SKIP))
                {
                    TOI_SKIP = current[0].first;
                }
                use_skip = true;
                continue;
            }

            std::array<bool, 3> check;
            Eigen::VectorX3d widthratio;
            widthratio.resize(3);
            check[0] = false;
            check[1] = false;
            check[2] = false;
            for (int i = 0; i < 3; i++)
            {
                widthratio(i) = widths(i) / tol(i);
                if (widths(i) > tol(i))
                    check[i] = true; // means this need to be checked
            }

            int split_i = -1;
            for (int i = 0; i < 3; i++)
            {
                if (check[i])
                {
                    if (check[(i + 1) % 3] && check[(i + 2) % 3])
                    {
                        if (widthratio(i) >= widthratio((i + 1) % 3) && widthratio(i) >= widthratio((i + 2) % 3))
                        {
                            split_i = i;
                            break;
                        }
                    }
                    if (check[(i + 1) % 3] && !check[(i + 2) % 3])
                    {
                        if (widthratio(i) >= widthratio((i + 1) % 3))
                        {
                            split_i = i;
                            break;
                        }
                    }
                    if (!check[(i + 1) % 3] && check[(i + 2) % 3])
                    {
                        if (widthratio(i) >= widthratio((i + 2) % 3))
                        {
                            split_i = i;
                            break;
                        }
                    }
                    if (!check[(i + 1) % 3] && !check[(i + 2) % 3])
                    {

                        split_i = i;
                        break;
                    }
                }
            }
            if (split_i < 0)
            {
                std::cout << "ERROR OCCURRED HERE, DID NOT FIND THE RIGHT DIMENSION TO SPLIT" << std::endl;
            }
            // Bisect the next dimension that is greater than its tolerance
            // int split_i;
            // for (int i = 1; i <= 3; i++) {
            //     split_i = (last_split + i) % 3;
            //     if (widths(split_i) > tol(split_i)) {
            //         break;
            //     }
            // }
            std::pair<Singleinterval, Singleinterval> halves = bisect(current[split_i]);
            if (!less_than(halves.first.first, halves.first.second))
            {
                std::cout << "OVERFLOW HAPPENS WHEN SPLITTING INTERVALS" << std::endl;
                return true;
            }
            if (!less_than(halves.second.first, halves.second.second))
            {
                std::cout << "OVERFLOW HAPPENS WHEN SPLITTING INTERVALS" << std::endl;
                return true;
            }
            if (check_vf)
            {
                //std::cout<<"*** check_vf"<<std::endl;
                if (split_i == 1)
                {
                    // assert(sum_no_larger_1(halves.first.first, current[2].first)==sum_no_larger_1_Rational(halves.first.first, current[2].first));
                    // assert(sum_no_larger_1(halves.second.first, current[2].first)==sum_no_larger_1_Rational(halves.second.first, current[2].first));

                    if (sum_no_larger_1(halves.second.first, current[2].first))
                    {
                        current[split_i] = halves.second;
                        istack.emplace(current, level + 1);
                    }
                    if (sum_no_larger_1(halves.first.first, current[2].first))
                    {
                        current[split_i] = halves.first;
                        istack.emplace(current, level + 1);
                    }
                }

                if (split_i == 2)
                {
                    //assert(sum_no_larger_1(halves.first.first, current[1].first)==sum_no_larger_1_Rational(halves.first.first, current[1].first));
                    //assert(sum_no_larger_1(halves.second.first, current[1].first)==sum_no_larger_1_Rational(halves.second.first, current[1].first));

                    if (sum_no_larger_1(halves.second.first, current[1].first))
                    {
                        current[split_i] = halves.second;
                        istack.emplace(current, level + 1);
                    }
                    if (sum_no_larger_1(halves.first.first, current[1].first))
                    {
                        current[split_i] = halves.first;
                        istack.emplace(current, level + 1);
                    }
                }
                if (split_i == 0)
                {
                    if (check_t_overlap)
                    {
                        if (interval_overlap_region(halves.second, 0, t_upper_bound))
                        {
                            current[split_i] = halves.second;
                            istack.emplace(current, level + 1);
                        }
                        if (interval_overlap_region(halves.first, 0, t_upper_bound))
                        {
                            current[split_i] = halves.first;
                            istack.emplace(current, level + 1);
                        }
                    }
                    else
                    {
                        current[split_i] = halves.second;
                        istack.emplace(current, level + 1);
                        current[split_i] = halves.first;
                        istack.emplace(current, level + 1);
                    }
                }
            }
            else
            {
                if (check_t_overlap && split_i == 0)
                {
                    if (interval_overlap_region(halves.second, 0, t_upper_bound))
                    {
                        current[split_i] = halves.second;
                        istack.emplace(current, level + 1);
                    }
                    if (interval_overlap_region(halves.first, 0, t_upper_bound))
                    {
                        current[split_i] = halves.first;
                        istack.emplace(current, level + 1);
                    }
                }
                else
                {
                    current[split_i] = halves.second;
                    istack.emplace(current, level + 1);
                    current[split_i] = halves.first;
                    istack.emplace(current, level + 1);
                }
            }
        }

        if (use_skip)
        {
            toi = Numccd2double(TOI_SKIP) * impact_ratio;
            return true;
        }

        return false;
    }

    bool interval_root_finder_double_horizontal_tree(
        const Eigen::VectorX3d &tol,
        const double co_domain_tolerance,
        double &toi,
        const bool check_vf,
        const std::array<double, 3> err, // this is the maximum error on each axis when calculating the vertices, err, aka, filter
        const double ms,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e,
        const double max_time,
        const int max_itr,
        double &output_tolerance)
    {

        bool check_t_overlap = max_time == 1 ? false : true; // if input max_time = 1, then no need to check overlap

        Numccd low_number;
        low_number.first = 0;
        low_number.second = 0; // low_number=0;
        Numccd up_number;
        up_number.first = 1;
        up_number.second = 0; // up_number=1;
        // initial interval [0,1]
        Singleinterval init_interval;
        init_interval.first = low_number;
        init_interval.second = up_number;
        //build interval set [0,1]x[0,1]x[0,1]
        Interval3 iset;
        iset[0] = init_interval;
        iset[1] = init_interval;
        iset[2] = init_interval;

        bool result = interval_root_finder_double_horizontal_tree(
            tol, co_domain_tolerance, iset, check_t_overlap, max_time, toi, check_vf,
            err, ms, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, max_itr, output_tolerance);
        if (result)
            return true;

        return false;
    }

#ifdef TIGHT_INCLUSION_USE_GMP
    bool interval_root_finder_Rational(
        const Eigen::VectorX3d &tol,

        std::array<std::pair<Rational, Rational>, 3> &final,
        const bool check_vf,
        const std::array<double, 3> err,
        const double ms,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e)
    {

        std::pair<Rational, Rational> interval01;
        interval01.first = 0;
        interval01.second = 1;
        std::array<std::pair<Rational, Rational>, 3> paracube;
        paracube[0] = interval01;
        paracube[1] = interval01;
        paracube[2] = interval01;

        // Stack of intervals and the last split dimension
        std::stack<std::pair<std::array<std::pair<Rational, Rational>, 3>, int>> istack;
        istack.emplace(paracube, -1);

        // current intervals
        std::array<std::pair<Rational, Rational>, 3> current;
        std::array<double, 3> err_and_ms;
        err_and_ms[0] = err[0] + ms;
        err_and_ms[1] = err[1] + ms;
        err_and_ms[2] = err[2] + ms;
        while (!istack.empty())
        {
            current = istack.top().first;
            int last_split = istack.top().second;
            istack.pop();
            igl::Timer timer;

            timer.start();
            bool zero_in = Origin_in_function_bounding_box_Rational(current, a0s, a1s, b0s, b1s, a0e, a1e, b0e, b1e, check_vf);
            timer.stop();
            time_rational += timer.getElapsedTimeInMicroSec();

            if (!zero_in)
                continue;
            timer.start();
            std::array<Rational, 3> widths = width(current);
            timer.stop();
            time21 += timer.getElapsedTimeInMicroSec();
            if (widths[0].to_double() <= tol(0) && widths[1].to_double() <= tol(1) && widths[2].to_double() <= tol(2))
            {
                final = current;
                return true;
            }

            // Bisect the next dimension that is greater than its tolerance
            int split_i;
            for (int i = 1; i <= 3; i++)
            {
                split_i = (last_split + i) % 3;
                if (widths[split_i].to_double() > tol(split_i))
                {
                    break;
                }
            }
            std::pair<std::pair<Rational, Rational>, std::pair<Rational, Rational>> halves = bisect(current[split_i]);

            if (check_vf)
            {
                if (split_i == 1)
                {
                    if (halves.first.first + current[2].first <= 1)
                    {
                        current[split_i] = halves.first;
                        istack.emplace(current, split_i);
                    }
                    if (halves.second.first + current[2].first <= 1)
                    {
                        current[split_i] = halves.second;
                        istack.emplace(current, split_i);
                    }
                }

                if (split_i == 2)
                {

                    if (halves.first.first + current[1].first <= 1)
                    {
                        current[split_i] = halves.first;
                        istack.emplace(current, split_i);
                    }
                    if (halves.second.first + current[1].first <= 1)
                    {
                        current[split_i] = halves.second;
                        istack.emplace(current, split_i);
                    }
                }
                if (split_i == 0)
                {
                    current[split_i] = halves.second;
                    istack.emplace(current, split_i);
                    current[split_i] = halves.first;
                    istack.emplace(current, split_i);
                }
            }
            else
            {
                current[split_i] = halves.second;
                istack.emplace(current, split_i);
                current[split_i] = halves.first;
                istack.emplace(current, split_i);
            }
        }
        return false;
    }
#endif
    void print_time_2()
    {
        std::cout << "how many times return from max_itr, " << refine_return << std::endl;
        std::cout << "origin predicates, " << time20 << std::endl;
        std::cout << "width, " << time21 << std::endl;
        std::cout << "bisect, " << time22 << std::endl;
        std::cout << "origin part1(evaluate 1 dimension), " << time23 << std::endl;
        std::cout << "origin part2(convert tuv), " << time24 << std::endl;
        std::cout << "time of call the vertex solving function, " << time25 << std::endl;
        std::cout << "how many times of interval check for this query, " << refine << std::endl;
    }
    double print_time_rational()
    {
        return time_rational;
    }

    std::array<double, 3> get_numerical_error(const std::vector<Eigen::Vector3d> &vertices,
                                              const bool &check_vf, const bool using_minimum_separation)
    {
        double eefilter;
        double vffilter;
        if (!using_minimum_separation)
        {
            eefilter = 6.217248937900877e-15;
            vffilter = 6.661338147750939e-15;
        }
        else
        {
            eefilter = 7.105427357601002e-15;
            vffilter = 7.549516567451064e-15;
        }

        double xmax = fabs(vertices[0][0]);
        double ymax = fabs(vertices[0][1]);
        double zmax = fabs(vertices[0][2]);
        for (int i = 0; i < vertices.size(); i++)
        {
            if (xmax < fabs(vertices[i][0]))
            {
                xmax = fabs(vertices[i][0]);
            }
            if (ymax < fabs(vertices[i][1]))
            {
                ymax = fabs(vertices[i][1]);
            }
            if (zmax < fabs(vertices[i][2]))
            {
                zmax = fabs(vertices[i][2]);
            }
        }
        double delta_x = xmax > 1 ? xmax : 1;
        double delta_y = ymax > 1 ? ymax : 1;
        double delta_z = zmax > 1 ? zmax : 1;
        std::array<double, 3> result;
        if (!check_vf)
        {
            result[0] = delta_x * delta_x * delta_x * eefilter;
            result[1] = delta_y * delta_y * delta_y * eefilter;
            result[2] = delta_z * delta_z * delta_z * eefilter;
        }
        else
        {
            result[0] = delta_x * delta_x * delta_x * vffilter;
            result[1] = delta_y * delta_y * delta_y * vffilter;
            result[2] = delta_z * delta_z * delta_z * vffilter;
        }
        return result;
    }
} // namespace inclusion_ccd
