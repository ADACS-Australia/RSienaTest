/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file RLegacy.h
 * \brief Methods producing numerically identical results to the R versions.
 *
 * Define the R_LEGACY symbol to enable the legacy mode.
 *
 * What is (should be) different in legacy mode:
 *
 *  - Assignment operators in Eigen are not / might not be equal to assignment
 *    of the result. That is use
 *      \code M = M / s; \endcode
 *    instead of
 *      \code M /= s; \endcode
 *
 *  - Matrix product: Use
 *      \code matrixProduct(M, v); \endcode
 *    instead of
 *      \code M * v; \endcode
 *
 *  - Inverse matrix: Use
 *      \code inverse(M); \endcode
 *    instead of
 *      \code M.inverse(); \endcode
 *
 *  - Most R functions use long double precision internally (eg. colSums, ...).
 *    Use
 *      \code M.cast<long double>().colwise().sum().cast<double>(); \endcode
 *    instead of
 *      \code M.colwise().sum(); \endcode
 *    or
 *      \code
 *      MatrixXld sum(MatrixXld::Zero(n, m));
 *       for (MatrixXd m;;) {
 *        sum = sum + m.cast<long double>();
 *      }
 *      MatrixXd mean = (sum / n).cast<double>();
 *      \endcode
 *    instead of
 *      \code
 *      MatrixXd sum(MatrixXd::Zero(n, m));
 *      for (MatrixXd m;;) {
 *        sum += m;
 *      }
 *      MatrixXd mean = sum / n;
 *      \endcode
 *
 *****************************************************************************/

#ifndef RSIENA_R_LEGACY_H_
#define RSIENA_R_LEGACY_H_

#include <Eigen/Dense>

#ifdef R_LEGACY

namespace siena {

Eigen::MatrixXd inverse(const Eigen::MatrixXd& rA);
Eigen::MatrixXd solve(const Eigen::MatrixXd& rA, const Eigen::MatrixXd& rB,
		double tol);
Eigen::VectorXd matrixProduct(const Eigen::MatrixXd& rX,
		const Eigen::VectorXd& rY);

} // namespace siena

#endif

#endif // RSIENA_R_LEGACY_H_
