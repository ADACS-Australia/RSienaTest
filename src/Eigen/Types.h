/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Types.h
 * \brief Typedefs for <tt>bool</tt> and <tt>long double</tt> Eigen types.
 *****************************************************************************/

#ifndef RSIENA_EIGEN_TYPES_H_
#define RSIENA_EIGEN_TYPES_H_

#include <Eigen/Core>

namespace Eigen {

//! Dynamic sized boolean array.
//!
typedef Array<bool, Dynamic, Dynamic> ArrayXb;
//! Dynamic sized boolean vector.
//!
typedef Matrix<bool, Dynamic, 1> VectorXb;
//! Dynamic sized boolean matrix.
//!
typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

//! Dynamic sized long double array.
//!
typedef Array<long double, Dynamic, Dynamic> ArrayXld;
//! Dynamic sized long double vector.
//!
typedef Matrix<long double, Dynamic, 1> VectorXld;
//! Dynamic sized long double matrix.
//!
typedef Matrix<long double, Dynamic, Dynamic> MatrixXld;

}

#endif // RSIENA_EIGEN_TYPES_H_
