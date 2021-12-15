/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Util.h
 * \brief Defines some more or less commonly use methods involving Eigen types.
 *
 * Most methods have a R_LEGACY implementation and a straight forward one.
 *****************************************************************************/

#ifndef EIGEN_UTIL_H_
#define EIGEN_UTIL_H_

#include "Eigen/Types.h"

#include <Eigen/Core>
#include <vector>

namespace Eigen {

/**
 * R_LEGACY aware column wise sum.
 *
 * Use this instead of `.colwise().sum()` if it should be aware of the legacy
 * mode.
 *
 * @param rhs Mapped matrix.
 * @return Column wise sum.
 */
inline VectorXd colwiseSum(const Map<MatrixXd>& rhs) {
#ifdef R_LEGACY
	return rhs.cast<long double>().colwise().sum().cast<double>();
#else
	return rhs.colwise().sum();
#endif
}

/**
 * Helper function to apply fixed parameter to matrices.
 *
 * Rows and columns of fixed parameters are set to 0, diagonal elements to 1.
 *
 * @param rFixed Boolean vector. True if the corresponding parameter is fixed.
 * @param[in,out] rM Reference to the matrix.
 */
inline void fixMatrix(const VectorXb& rFixed, MatrixXd& rM) {
	for (int i = 0; i < rFixed.size(); ++i) {
		if (rFixed[i]) {
			rM.row(i).fill(0);
			rM.col(i).fill(0);
			rM(i, i) = 1;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Variance functions
///////////////////////////////////////////////////////////////////////////////

/**
 * R_LEGACY aware mean.
 *
 * Use this instead of `.mean()` if it should be aware of the legacy mode. In
 * legacy mode this uses a two pass mean like the mean function in R.
 *
 * @param rX A vector.
 * @return Mean value.
 */
inline long double mean(const VectorXd& rX) {
#ifdef R_LEGACY
	// from: R/src/library/stats/src/cov.c
	// 2 passes for better accuracy
	long double sum = 0;
	for (int i = 0; i < rX.size(); ++i) {
		sum += static_cast<long double>(rX[i]);
	}
	long double mean = sum / rX.size();
	sum = 0;
	for (int i = 0; i < rX.size(); ++i) {
		sum += (static_cast<long double>(rX[i]) - mean);
	}
	mean += sum / rX.size();
	return mean;
#else
	return rX.mean();
#endif
}

/**
 * R_LEGACY aware variance.
 *
 * @param rX A vector.
 * @return Variance of `rX`.
 */
inline double variance(const VectorXd& rX) {
#ifdef R_LEGACY
	// from: R/src/library/stats/src/cov.c
	long double xm = mean(rX);
	long double sum = 0;
	long double m;
	for (int i = 0; i < rX.size(); ++i) {
		m = rX[i] - xm;
		sum += m * m;
	}
	sum /= rX.size() - 1;
	return sum;
#else
//	assert(rhs.rows() == 1 || rhs.cols() == 1);
	VectorXd centered = rX.array() - rX.mean();
	return centered.cwiseProduct(centered).sum() / (rX.size() - 1);
#endif
}

/**
 * R_LEGACY aware variance.
 *
 * @param rX A vector.
 * @param rY Another vector.
 * @return Covariance of `rX` and `rY`.
 */
inline double covariance(const VectorXd& rX, const VectorXd& rY) {
#ifdef R_LEGACY
	// from: R/src/library/stats/src/cov.c
	long double xm = mean(rX);
	long double ym = mean(rY);
	long double sum = 0;
	for (int i = 0; i < rX.size(); ++i) {
		sum += (rX[i] - xm) * (rY[i] - ym);
	}
	sum /= rX.size() - 1;
	return sum;
#else
	//	assert(rX.rows() == 1 || rX.cols() == 1);
	//	assert(rX.rows() == rY.rows() || rX.cols() == rY.cols());
	return (rX.array() - rX.mean()).cwiseProduct(rY.array() - rY.mean()).sum()
			/ (rX.rows() - 1);
#endif
}

/**
 * Covariance matrix of a series of samples.
 *
 * @param rSamples Matrix of row wise samples.
 * @return Covariance of `rSamples`.
 */
inline MatrixXd covarianceMatrix(const MatrixXd& rSamples) {
	MatrixXd centered = rSamples.rowwise() - rSamples.colwise().mean();
	return (centered.adjoint() * centered) / (rSamples.rows() - 1);
}

///////////////////////////////////////////////////////////////////////////////
// std::vector of Eigen objects
///////////////////////////////////////////////////////////////////////////////

/**
 * @param rV Vector of matrices.
 * @return Entry wise sum of rV.
 */
template<typename MapOrMatrixXd>
inline MatrixXd vSum(const std::vector<MapOrMatrixXd>& rV) {
//#ifdef R_LEGACY
//	MatrixXld sum = MatrixXd::Zero(rV.rows(), rV.cols());
//	typedef typename std::vector<MapOrMatrixXd>::const_iterator CItr;
//	for (CItr it = rV.begin(); it != rV.end(); ++it) {
//		sum = sum + (*it).cast<long double>();
//	}
//	return sum.cast<double>();
//#else
	MatrixXd sum = MatrixXd::Zero(rV[0].rows(), rV[0].cols());
	typedef typename std::vector<MapOrMatrixXd>::const_iterator CItr;
	for (CItr it = rV.begin(); it != rV.end(); ++it) {
		sum = sum + (*it);
	}
	return sum;
//#endif
}

/**
 * @param rV Vector of matrices.
 * @return Entry wise mean of rV.
 */
template<typename MapOrMatrixXd>
inline MatrixXd vMean(const std::vector<MapOrMatrixXd>& rV) {
//#ifdef R_LEGACY
//	MatrixXld sum(rV.rows(), rV.cols());
//	typedef typename std::vector<MapOrMatrixXd>::const_iterator CItr;
//	for (CItr it = rV.begin(); it != rV.end(); ++it) {
//		sum = sum + (*it).cast<long double>();
//	}
//	return (sum / rV.size()).cast<double>();
//#else
	return vSum(rV) / rV.size();
//#endif
}

/**
 * @param rV Vector of matrices.
 * @return Vector of columns sum vectors.
 */
template<typename MapOrMatrixXd>
inline std::vector<VectorXd> vColwiseSum(const std::vector<MapOrMatrixXd>& rV) {
	std::vector<VectorXd> result(rV.size());
	transform(rV.begin(), rV.end(), result.begin(), colwiseSum);
	return result;
}

} // namespace Eigen

#endif // EIGEN_UTIL_H_
