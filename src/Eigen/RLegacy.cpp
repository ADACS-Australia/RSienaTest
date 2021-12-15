/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file RLegacy.cpp
 * \brief Implements RLegacy.h.
 *****************************************************************************/

#include "RLegacy.h"

#include "logger/Logger.h"

#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#include <limits>

using namespace Eigen;
using namespace siena::logger;

#ifdef R_LEGACY

namespace siena {

/**
 * Matrix product from R/src/main/array.c
 *
 * @param x First matrix.
 * @param nrx Number of rows of `x`.
 * @param ncx Number of cols of `x`.
 * @param y Second matrix.
 * @param nry Number of rows of `y`.
 * @param ncy Number of cols of `y`.
 * @param z Result matrix.
 */
void matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy,
		double *z) {
	int NRX = nrx;
	int NRY = nry;

	if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
		Rboolean have_na = FALSE;
		/* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
		 * The test is only O(n) here.
		 */
		for (int i = 0; i < NRX * ncx; i++) {
			if (ISNAN(x[i])) {
				have_na = TRUE;
				break;
			}
		}
		if (!have_na) {
			for (int i = 0; i < NRY * ncy; i++) {
				if (ISNAN(y[i])) {
					have_na = TRUE;
					break;
				}
			}
		}
		if (have_na) {
			long double sum;
			for (int i = 0; i < nrx; i++) {
				for (int k = 0; k < ncy; k++) {
					sum = 0.0;
					for (int j = 0; j < ncx; j++)
						sum += x[i + j * NRX] * y[j + k * NRY];
					z[i + k * NRX] = (double) sum;
				}
			}
		} else {
			char *transa = "N", *transb = "N"; // 'N' no transpose
			double one = 1.0, zero = 0.0;
			F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one, x, &nrx, y,
					&nry, &zero, z, &nrx);
		}
	} else {
		/* zero-extent operations should return zeroes */
		for (int i = 0; i < NRX * ncy; i++) {
			z[i] = 0;
		}
	}
}

/**
 * Solve the equation A*X=B using Lapack.
 *
 * @param rA Square matrix.
 * @param rB Square matrix of the same size as `rA`.
 * @param tol Numerical tolerance.
 * @return Solution to the equation X.
 */
MatrixXd solve(const MatrixXd& rA, const MatrixXd& rB, double tol) {
	assert(rA.rows() == rA.cols()); // Square
	assert(rA.rows() == rB.rows() && rA.cols() == rB.cols()); // Same size

	int n = rA.rows();
	int* ipiv = new int[n];

	MatrixXd aCopy(rA); // A is modified
	MatrixXd bOrSolved(rB); // Return value
	int state; // Return state

	F77_CALL(dgesv)(&n, &n, aCopy.data(), &n, ipiv, bOrSolved.data(), &n,
			&state);

	if (state < 0) {
		LOGF(Priority::VERBOSE,
				"argument %d of Lapack routine %s had invalid value", -state,
				"dgesv");
	}
	if (state > 0) {
		LOGF(Priority::VERBOSE,
				"Lapack routine %s: system is exactly singular: U[%d,%d] = 0",
				"dgesv", state, state);
	}
	if (tol > 0) {
		double anorm = F77_CALL(dlange)("1", &n, &n, rA.data(), &n,
				(double*) NULL);
		double* work = new double[4 * n];
		double rcond;
		F77_CALL(dgecon)("1", &n, aCopy.data(), &n, &anorm, &rcond, work, ipiv,
				&state);
		if (rcond < tol) {
			LOGF(Priority::VERBOSE,
					"system is computationally singular: reciprocal condition number = %g",
					rcond);
		}
		delete[] work;
	}
	delete[] ipiv;
	return bOrSolved;
}

/**
 * Matrix vector product using and implementation equivalent to the R matrix
 * product.
 *
 * @param rX A matrix.
 * @param rY A vector.
 * @return `rX` times `rY`.
 */
VectorXd matrixProduct(const MatrixXd& rX, const VectorXd& rY) {
	LOG(Priority::VERBOSE, "");
	VectorXd z(rX.rows());
	matprod(const_cast<MatrixXd&>(rX).data(), rX.rows(), rX.cols(),
			const_cast<VectorXd&>(rY).data(), rY.rows(), rY.cols(), z.data());
	return z;
}

/**
 * Matrix inversion equivalent to the one in R.
 *
 * @param rA Matrix to invert.
 * @return Inverted matrix.
 */
MatrixXd inverse(const MatrixXd& rA) {
	return solve(rA, MatrixXd::Identity(rA.rows(), rA.cols()),
			std::numeric_limits<double>::epsilon());
}

} // namespace siena

#endif // R_LEGACY
