#pragma once
#include<tight_inclusion/interval.hpp>
namespace inclusion_ccd{
    void convert_tuv_to_array(const Interval3& itv, 
    std::array<double,8>& t_up,std::array<double,8>&t_dw,
    std::array<double,8>& u_up,std::array<double,8>&u_dw,
    std::array<double,8>& v_up,std::array<double,8>&v_dw);

    std::array<double,8> function_vf(
        const double& vs,
        const double& t0s,
        const double& t1s,const double& t2s,
    const double& ve,const double& t0e,
    const double& t1e,const double& t2e,
    const std::array<double,8>& t_up,const std::array<double,8>& t_dw,
    const std::array<double,8>& u_up,const std::array<double,8>& u_dw,
    const std::array<double,8>& v_up,const std::array<double,8>& v_dw);

    std::array<double,8> function_ee(
    const double& a0s,const double& a1s,const double& b0s,const double& b1s,
    const double& a0e,const double& a1e,const double& b0e,const double& b1e,
    const std::array<double,8>& t_up,const std::array<double,8>& t_dw,
    const std::array<double,8>& u_up,const std::array<double,8>& u_dw,
    const std::array<double,8>& v_up,const std::array<double,8>& v_dw);
    long power(const long a, const int b);
    //class get_
}