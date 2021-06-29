// Time-of-impact computation for rigid bodies with angular trajectories.
#pragma once

#include <tight_inclusion/interval.hpp>
#include <tight_inclusion/Types.hpp>
#include <array>

/**
 * @namespace ccd:
 * @brief
 */
namespace inclusion_ccd
{

	// This function can give you the answer of continous collision detection with minimum
	// seperation, and the earlist collision time if collision happens.
	// err is the filters calculated using the bounding box of the simulation scene.
	// If you are checking a single query without a scene, please set it as [-1,-1,-1].
	// ms is the minimum seperation. should set: ms < max(abs(x),1), ms < max(abs(y),1), ms < max(abs(z),1) of the QUERY (NOT THE SCENE!).
	// toi is the earlist time of collision if collision happens. If there is no collision, toi will be infinate.
	// tolerance is a user - input solving precision. we suggest to use 1e-6.
	// t_max is the upper bound of the time interval [0,t_max] to be checked. 0<=t_max<=1
	// max_itr is a user-defined value to terminate the algorithm earlier, and return a result under current
	// precision. please set max_itr either a big number like 1e7, or -1 which means it will not be terminated
	// earlier and the precision will be user-defined precision -- tolerance.
	// output_tolerance is the precision under max_itr ( > 0). if max_itr < 0, output_tolerance = tolerance;
	// CCD_TYPE is a switch to choose root-finding methods.
	// 0 is normal ccd,
	// 1 is ccd with input time interval upper bound, using real tolerance, max_itr and horizontal tree,
	bool edgeEdgeCCD_double(
		const Vector3d &a0_start,
		const Vector3d &a1_start,
		const Vector3d &b0_start,
		const Vector3d &b1_start,
		const Vector3d &a0_end,
		const Vector3d &a1_end,
		const Vector3d &b0_end,
		const Vector3d &b1_end,
		const std::array<Scalar, 3> &err,
		const Scalar ms,
		Scalar &toi,
		const Scalar tolerance,
		const Scalar t_max,
		const int max_itr,
		Scalar &output_tolerance,
		const int CCD_TYPE = 1);

	// This function can give you the answer of continous collision detection with minimum
	// seperation, and the earlist collision time if collision happens.
	// err is the filters calculated using the bounding box of the simulation scene.
	// If you are checking a single query without a scene, please set it as [-1,-1,-1].
	// ms is the minimum seperation. should set: ms < max(abs(x),1), ms < max(abs(y),1), ms < max(abs(z),1) of the QUERY (NOT THE SCENE!).
	// toi is the earlist time of collision if collision happens. If there is no collision, toi will be infinate.
	// tolerance is a user - input solving precision. we suggest to use 1e-6.
	// t_max is the upper bound of the time interval [0,t_max] to be checked. 0<=t_max<=1
	// max_itr is a user-defined value to terminate the algorithm earlier, and return a result under current
	// precision. please set max_itr either a big number like 1e7, or -1 which means it will not be terminated
	// earlier and the precision will be user-defined precision -- tolerance.
	// output_tolerance is the precision under max_itr ( > 0). if max_itr < 0, output_tolerance = tolerance;
	// CCD_TYPE is a switch to choose root-finding methods.
	// 0 is normal ccd,
	// 1 is ccd with input time interval upper bound, using real tolerance, max_itr and horizontal tree,
	bool vertexFaceCCD_double(
		const Vector3d &vertex_start,
		const Vector3d &face_vertex0_start,
		const Vector3d &face_vertex1_start,
		const Vector3d &face_vertex2_start,
		const Vector3d &vertex_end,
		const Vector3d &face_vertex0_end,
		const Vector3d &face_vertex1_end,
		const Vector3d &face_vertex2_end,
		const std::array<Scalar, 3> &err,
		const Scalar ms,
		Scalar &toi,
		const Scalar tolerance,
		const Scalar t_max,
		const int max_itr,
		Scalar &output_tolerance,
		const int CCD_TYPE = 1);
#ifdef TIGHT_INCLUSION_USE_GMP

	// this version is an naive implementation of Tight-Inclusion CCD without optimizations
	bool edgeEdgeCCD_rational(
		const Vector3d &a0s,
		const Vector3d &a1s,
		const Vector3d &b0s,
		const Vector3d &b1s,
		const Vector3d &a0e,
		const Vector3d &a1e,
		const Vector3d &b0e,
		const Vector3d &b1e,
		const std::array<Scalar, 3> &err,
		const Scalar ms,
		Scalar &toi);

	// this version is an naive implementation of Tight-Inclusion CCD without optimizations
	bool vertexFaceCCD_rational(
		const Vector3d &vertex_start,
		const Vector3d &face_vertex0_start,
		const Vector3d &face_vertex1_start,
		const Vector3d &face_vertex2_start,
		const Vector3d &vertex_end,
		const Vector3d &face_vertex0_end,
		const Vector3d &face_vertex1_end,
		const Vector3d &face_vertex2_end,
		const std::array<Scalar, 3> &err,
		const Scalar ms,
		Scalar &toi);

#endif
	bool using_rational_method();

} // namespace inclusion_ccd
