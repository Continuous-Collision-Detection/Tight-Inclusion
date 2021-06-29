// An interval object.
#pragma once

#include <Eigen/Core>
#include <tight_inclusion/Types.hpp>
#ifdef TIGHT_INCLUSION_USE_GMP
#include <tight_inclusion/Rational.hpp>
#endif

namespace inclusion_ccd
{

    // typedef std::array<std::pair<int,int>,3> Paraccd;//<k,n> pair present parameters (u,v,t) which are presented as k/(2^n)
    typedef std::pair<long, int>
        Numccd; //<k,n> pair present a number k/pow(2,n)
    typedef std::pair<Numccd, Numccd>
        Singleinterval; // a interval presented by two double numbers
    typedef std::array<Singleinterval, 3> Interval3; // 3 dimesional interval
#ifdef TIGHT_INCLUSION_USE_GMP
    typedef Eigen::Matrix<Rational, 3, 1, Eigen::ColMajor | Eigen::DontAlign>
        Vector3r;
#endif

} // namespace inclusion_ccd
