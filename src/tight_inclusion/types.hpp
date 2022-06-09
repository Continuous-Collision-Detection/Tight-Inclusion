#pragma once

#include <Eigen/Core>

//#define CHECK_QUEUE_SIZE

//#define TI_LIMIT_QUEUE_SIZE
#ifdef TI_LIMIT_QUEUE_SIZE
#define MAX_QSIZE 1000
#endif

namespace ticcd {
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
#else
    typedef float Scalar;
#endif
    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
    typedef Eigen::Array<Scalar, 3, 1> Array3;
} // namespace ticcd
