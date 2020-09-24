// A root finder using interval arithmetic.
#pragma once

#include <array>
#include <functional>
#include <vector>

#include <tight_inclusion/interval.hpp>

namespace inclusion_ccd
{

    // this version cannot give the impact time at t=1, although this collision can
    // be detected at t=0 of the next time step, but still may cause problems in
    // line-search based physical simulation
    bool interval_root_finder_double_normalCCD(
        const Eigen::VectorX3d &tol,
        //Eigen::VectorX3I& x,// result interval
        // Interval3& final,
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
        const Eigen::Vector3d &b1e);

    // this version cannot give the impact time at t=1.
    // max_itr is a user defined maximum iteration time. if < 0, then
    // it will run until stack empty; otherwise the algorithm will stop when
    // iteration time reaches max_itr, and return a solution precision output_tolerance
    // it uses interval t = [0, max_time] instead of t = [0,1]
    // 0<=max_time <=1
    // tree searching order is horizontal
    bool interval_root_finder_double_horizontal_tree(
        const Eigen::VectorX3d &tol,
        const double co_domain_tolerance,
        //Eigen::VectorX3I& x,// result interval
        // Interval3& final,
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
        double &output_tolerance);
    // return power t. n=result*2^t
    long reduction(const long n, long &result);
    std::pair<Singleinterval, Singleinterval> bisect(const Singleinterval &inter);
    // calculate the sign of f. dim is the dimension we are evaluating.
    template <typename T>
    T function_f_ee(
        const Numccd &tpara,
        const Numccd &upara,
        const Numccd &vpara,
        const T &type,
        const int dim,
        const Eigen::Vector3d &a0s,
        const Eigen::Vector3d &a1s,
        const Eigen::Vector3d &b0s,
        const Eigen::Vector3d &b1s,
        const Eigen::Vector3d &a0e,
        const Eigen::Vector3d &a1e,
        const Eigen::Vector3d &b0e,
        const Eigen::Vector3d &b1e);
    template <typename T>
    T function_f_vf(
        const Numccd &tpara,
        const Numccd &upara,
        const Numccd &vpara,
        const T &type,
        const int dim,
        const Eigen::Vector3d &vs,
        const Eigen::Vector3d &t0s,
        const Eigen::Vector3d &t1s,
        const Eigen::Vector3d &t2s,

        const Eigen::Vector3d &ve,
        const Eigen::Vector3d &t0e,
        const Eigen::Vector3d &t1e,
        const Eigen::Vector3d &t2e);

#ifdef TIGHT_INCLUSION_USE_GMP
    bool interval_root_finder_Rational(
        const Eigen::VectorX3d &tol,
        // Eigen::VectorX3I& x,// result interval
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
        const Eigen::Vector3d &b1e);
#endif
    void print_time_2();
    double print_time_rational();
    int print_refine();

    // get the filter of ccd. the inputs are the vertices of the bounding box of the simulation scene
    std::array<double, 3> get_numerical_error(
        const std::vector<Eigen::Vector3d> &vertices, const bool &check_vf, const bool using_minimum_separation);
} // namespace inclusion_ccd
