// Time-of-impact computation for rigid bodies with angular trajectories.
#pragma once

#include <tight_inclusion/interval.hpp>
#include <tight_inclusion/types.hpp>
#include <array>

namespace ticcd {
    static constexpr bool DEFAULT_NO_ZERO_TOI = false;
    static constexpr Scalar DEFAULT_CCD_DISTANCE_TOL = 1e-6;

    enum class CCDRootFindingMethod {
        DEPTH_FIRST_SEARCH,
        BREADTH_FIRST_SEARCH,
    };

    /// @brief This function can give you the answer of continous collision detection with minimum
    /// seperation, and the earlist collision time if collision happens.
    ///
    /// @param[in] err The filters calculated using the bounding box of the simulation scene.
    ///                If you are checking a single query without a scene, please set it as {-1,-1,-1}.
    /// @param[in] ms The minimum seperation. should set: ms < max(abs(x),1), ms < max(abs(y),1), ms < max(abs(z),1) of the QUERY (NOT THE SCENE!).
    /// @param[out] toi The earlist time of collision if collision happens. If there is no collision, toi will be infinate.
    /// @param[in] tolerance A user - input solving precision. We suggest to use 1e-6.
    /// @param[in] t_max The upper bound of the time interval [0,t_max] to be checked. 0<=t_max<=1
    /// @param[in] max_itr A user-defined value to terminate the algorithm earlier, and return a result under current
    ///                    precision. please set max_itr either a big number like 1e7, or -1 which means it will not be terminated
    ///                    earlier and the precision will be user-defined precision -- tolerance.
    /// @param[out] output_tolerance The precision under max_itr ( > 0). if max_itr < 0, output_tolerance = tolerance;
    /// @param[in] no_zero_toi Refine further if a zero toi is produced (assumes not initially in contact).
    bool edgeEdgeCCD(
        const Vector3 &a0_start,
        const Vector3 &a1_start,
        const Vector3 &b0_start,
        const Vector3 &b1_start,
        const Vector3 &a0_end,
        const Vector3 &a1_end,
        const Vector3 &b0_end,
        const Vector3 &b1_end,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi,
        const Scalar tolerance,
        const Scalar t_max,
        const int max_itr,
        Scalar &output_tolerance,
        bool no_zero_toi = DEFAULT_NO_ZERO_TOI,
        const CCDRootFindingMethod ccd_method =
            CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    /// This function can give you the answer of continous collision detection with minimum
    /// seperation, and the earlist collision time if collision happens.
    ///
    /// @param[in] err The filters calculated using the bounding box of the simulation scene.
    ///                If you are checking a single query without a scene, please set it as {-1,-1,-1}.
    /// @param[in] ms The minimum seperation. should set: ms < max(abs(x),1), ms < max(abs(y),1), ms < max(abs(z),1) of the QUERY (NOT THE SCENE!).
    /// @param[out] toi The earlist time of collision if collision happens. If there is no collision, toi will be infinate.
    /// @param[in] tolerance A user - input solving precision. We suggest to use 1e-6.
    /// @param[in] t_max The upper bound of the time interval [0,t_max] to be checked. 0<=t_max<=1
    /// @param[in] max_itr A user-defined value to terminate the algorithm earlier, and return a result under current
    ///                    precision. please set max_itr either a big number like 1e7, or -1 which means it will not be terminated
    ///                    earlier and the precision will be user-defined precision -- tolerance.
    /// @param[out] output_tolerance The precision under max_itr ( > 0). if max_itr < 0, output_tolerance = tolerance;
    /// @param[in] no_zero_toi Refine further if a zero toi is produced (assumes not initially in contact).
    bool vertexFaceCCD(
        const Vector3 &vertex_start,
        const Vector3 &face_vertex0_start,
        const Vector3 &face_vertex1_start,
        const Vector3 &face_vertex2_start,
        const Vector3 &vertex_end,
        const Vector3 &face_vertex0_end,
        const Vector3 &face_vertex1_end,
        const Vector3 &face_vertex2_end,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi,
        const Scalar tolerance,
        const Scalar t_max,
        const int max_itr,
        Scalar &output_tolerance,
        bool no_zero_toi = DEFAULT_NO_ZERO_TOI,
        const CCDRootFindingMethod ccd_method =
            CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    Array3 compute_face_vertex_tolerances(
        const Vector3 &vs,
        const Vector3 &f0s,
        const Vector3 &f1s,
        const Vector3 &f2s,
        const Vector3 &ve,
        const Vector3 &f0e,
        const Vector3 &f1e,
        const Vector3 &f2e,
        const Scalar distance_tolerance = DEFAULT_CCD_DISTANCE_TOL);

    Array3 compute_edge_edge_tolerances(
        const Vector3 &edge0_vertex0_start,
        const Vector3 &edge0_vertex1_start,
        const Vector3 &edge1_vertex0_start,
        const Vector3 &edge1_vertex1_start,
        const Vector3 &edge0_vertex0_end,
        const Vector3 &edge0_vertex1_end,
        const Vector3 &edge1_vertex0_end,
        const Vector3 &edge1_vertex1_end,
        const Scalar distance_tolerance = DEFAULT_CCD_DISTANCE_TOL);

    long return_queue_size();

    // these function are designed to test the performance of floating point vertion but with double inputs
#ifdef TIGHT_INCLUSION_FWDI
    bool edgeEdgeCCD(
        const Eigen::Vector3d &a0_start,
        const Eigen::Vector3d &a1_start,
        const Eigen::Vector3d &b0_start,
        const Eigen::Vector3d &b1_start,
        const Eigen::Vector3d &a0_end,
        const Eigen::Vector3d &a1_end,
        const Eigen::Vector3d &b0_end,
        const Eigen::Vector3d &b1_end,
        const std::array<double, 3> &err,
        const double ms,
        double &toi,
        const double tolerance,
        const double t_max,
        const int max_itr,
        double &output_tolerance,
        bool no_zero_toi = DEFAULT_NO_ZERO_TOI,
        const CCDRootFindingMethod ccd_method =
            CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    bool vertexFaceCCD(
        const Eigen::Vector3d &vertex_start,
        const Eigen::Vector3d &face_vertex0_start,
        const Eigen::Vector3d &face_vertex1_start,
        const Eigen::Vector3d &face_vertex2_start,
        const Eigen::Vector3d &vertex_end,
        const Eigen::Vector3d &face_vertex0_end,
        const Eigen::Vector3d &face_vertex1_end,
        const Eigen::Vector3d &face_vertex2_end,
        const std::array<double, 3> &err,
        const double ms,
        double &toi,
        const double tolerance,
        const double t_max,
        const int max_itr,
        double &output_tolerance,
        bool no_zero_toi = DEFAULT_NO_ZERO_TOI,
        const CCDRootFindingMethod ccd_method =
            CCDRootFindingMethod::BREADTH_FIRST_SEARCH);
#endif
} // namespace ticcd
