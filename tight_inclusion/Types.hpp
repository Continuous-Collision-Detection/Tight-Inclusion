#pragma once
#include <Eigen/Core>
#define CHECK_QUEUE_SIZE

inline bool ticcd_using_double()
{
#ifdef TIGHT_INCLUSION_DOUBLE
    return true;
#else
    return false;
#endif
}
#ifdef TIGHT_INCLUSION_DOUBLE
typedef double Scalar;
typedef Eigen::Vector3d Vector3d;
#else
typedef float Scalar;
typedef Eigen::Vector3f Vector3d;
#endif
