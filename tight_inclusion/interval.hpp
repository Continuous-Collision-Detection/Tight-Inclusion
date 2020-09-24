// An interval object.
#pragma once

#include <string>

#include <Eigen/Core>
#include <boost/numeric/interval.hpp>
#include <tight_inclusion/Rational.hpp>

namespace inclusion_ccd {

// namespace interval_options {
//     typedef boost::numeric::interval_lib::checking_base<double> CheckingPolicy;
// } // namespace interval_options

// #if defined(__APPLE__)

// // clang-format off
// #warning "Rounding modes seem to be broken with trigonometric functions on macOS, unable to compute exact interval arithmetic!"
// // clang-format on
// typedef boost::numeric::interval<
//     double,
//     boost::numeric::interval_lib::policies<
//         boost::numeric::interval_lib::save_state<
//             boost::numeric::interval_lib::rounded_transc_exact<double>>,
//         interval_options::CheckingPolicy>>
//     Interval;

// #else

// // Use proper rounding arithmetic
// typedef boost::numeric::interval<
//     double,
//     boost::numeric::interval_lib::policies<
//         boost::numeric::interval_lib::save_state<
//             boost::numeric::interval_lib::rounded_transc_std<double>>,
//         interval_options::CheckingPolicy>>
//     Interval;

// #endif

// typedef std::array<std::pair<int,int>,3> Paraccd;//<k,n> pair present parameters (u,v,t) which are presented as k/(2^n) 
typedef std::pair<long,int> Numccd;//<k,n> pair present a number k/pow(2,n)
typedef std::pair<Numccd,Numccd> Singleinterval;// a interval presented by two double numbers
typedef std::array<Singleinterval,3> Interval3; // 3 dimesional interval
typedef Eigen::Matrix<Rational, 3, 1, Eigen::ColMajor | Eigen::DontAlign> Vector3r;


} // namespace intervalccd

namespace Eigen {

typedef Matrix<double, 3, 1> Vector3d;
typedef Matrix<double, Dynamic, 1, ColMajor, 3, 1> VectorX3d;

// template <typename T, int dim, int max_dim = dim>
// using Vector = Matrix<T, dim, 1, ColMajor, max_dim, 1>;

// #if EIGEN_MAJOR_VERSION >= 3
// namespace internal {
//     // template <typename X, typename S, typename P>
//     // struct is_convertible<X, boost::numeric::interval<S, P>> {
//     //     enum { value = is_convertible<X, S>::value };
//     // };

//     // template <typename S, typename P1, typename P2>
//     // struct is_convertible<
//     //     boost::numeric::interval<S, P1>,
//     //     boost::numeric::interval<S, P2>> {
//     //     enum { value = true };
//     // };
// } // namespace internal
// #endif
} // namespace Eigen
