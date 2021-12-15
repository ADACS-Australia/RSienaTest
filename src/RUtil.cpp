/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file RUtil.cpp
 * \brief Implements RUtil.h.
 *****************************************************************************/

#include "RUtil.h"

#include <Rinternals.h>
#include <R_ext/Random.h>

#include <cstring>

using namespace std;
using namespace Eigen;

namespace siena {

/**
 * Retrieves the R random number generator state.
 *
 * @param seed Vector to be filled with the random number generator state.
 */
void getSeed(std::vector<int>& seed) {
	// Put the RNG state from private R memory to the R object '.Random.seed'.
	PutRNGstate();
	// Get the value of '.Random.seed' from the R context.
	const SEXP sSeed = findVar(install(".Random.seed"), R_GlobalEnv);
	// That is a integer vector.
	const size_t len = length(sSeed);
	if (len != seed.size()) {
		seed.resize(len, 0);
	}
	// Since vectors have guaranteed memory layout, we can copy directly.
	// memcpy(dst, src, n)
	memcpy(&seed[0], INTEGER(sSeed), len * sizeof(int));
}

/**
 * Sets the R random number generator state.
 *
 * @param seed Vector of integers.
 */
void setSeed(const std::vector<int>& seed) {
	// Create a SEXP holding the seed vector.
	SEXP sSeed = PROTECT(allocVector(INTSXP, seed.size()));
	memcpy(INTEGER(sSeed), &seed[0], seed.size() * sizeof(int));
	// Define the '.Random.seed' object.
	defineVar(install(".Random.seed"), sSeed, R_GlobalEnv);
	// Let R copy the RNG state.
	GetRNGstate();
	UNPROTECT(1); // sSeed
}

/**
 * Helper function creating external pointer SEXPs and registers a finalizer on them.
 *
 * @param cObject External pointer.
 * @param finializer Finalizer function.
 * @return EXTPTRSXP
 */
SEXP wrapPointer(void* cObject, void finalizer(SEXP)) {
	SEXP wrapper = PROTECT(R_MakeExternalPtr(cObject, R_NilValue, R_NilValue));
	if (finalizer) {
		R_RegisterCFinalizerEx(wrapper, finalizer, FALSE);
	}
	UNPROTECT(1);
	return wrapper;
}

/**
 * @param rString string to rify.
 * @return STRSXP copy of `rString`.
 */
SEXP rifyString(const std::string& rString) {
	SEXP sString = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(sString, 0, mkChar(rString.c_str()));
	UNPROTECT(1);
	return sString;
}

/**
 * @param rVector Vector of integers to rify.
 * @return R vector (INTSXP) holding a copy of `rVector`.
 */
SEXP rifyIntVector(const std::vector<int>& rVector) {
	SEXP sVector = PROTECT(allocVector(INTSXP, rVector.size()));
	memcpy(INTEGER(sVector), &rVector[0], rVector.size() * sizeof(int));
	UNPROTECT(1);
	return sVector;
}

/**
 * @param rVector Vector of boolean values to rify.
 * @return R vector (LGLSXP) holding a copy of `rVector`.
 */
SEXP rifyVectorXb(const Eigen::VectorXb& rVector) {
	SEXP sVector = PROTECT(allocVector(LGLSXP, rVector.size()));
	int* data = LOGICAL(sVector);
	for (int i = 0; i < rVector.size(); ++i) {
		data[i] = rVector[i];
	}
	UNPROTECT(1);
	return sVector;
}

VectorXb asVectorXb(SEXP sVector) {
	int* data = LOGICAL(sVector);
	VectorXb vector(length(sVector));
	for (int i = 0; i < vector.size(); ++i) {
		vector[i] = data[i] == 1;
	}
	return vector;
}

/**
 * @param rVector Eigen vector to rify.
 * @return R vector (REALSXP) holding a copy of `rVector`.
 */
SEXP rifyVectorXd(const Eigen::VectorXd& rVector) {
	SEXP sVector = PROTECT(allocVector(REALSXP, rVector.size()));
	double* data = REAL(sVector);
	for (int i = 0; i < rVector.size(); ++i) {
		data[i] = rVector[i];
	}
	UNPROTECT(1);
	return sVector;
}

/**
 * @param rMatrix Eigen matrix to rify.
 * @return R vector (REALSXP) holding a copy of `rMatrix`
 */
SEXP rifyMatrixXd(const Eigen::MatrixXd& rMatrix) {
	SEXP sMatrix = PROTECT(allocVector(REALSXP, rMatrix.size()));
	double* data = REAL(sMatrix);
	int i = 0;
	for (int c = 0; c < rMatrix.cols(); ++c) {
		for (int r = 0; r < rMatrix.rows(); ++r) {
			data[i] = rMatrix(r, c);
			++i;
		}
	}
	SEXP sDim = PROTECT(allocVector(INTSXP, 2));
	int* dim = INTEGER(sDim);
	dim[0] = rMatrix.rows();
	dim[1] = rMatrix.cols();
	setAttrib(sMatrix, R_DimSymbol, sDim);
	UNPROTECT(2); // sMatrix, sDim
	return sMatrix;
}

SEXP rifyVectorMatrixXd(const vector<MatrixXd>& rTensor) {
	// dimensions
	SEXP sDim = PROTECT(allocVector(INTSXP, 3));
	int* dim = INTEGER(sDim);
	dim[0] = rTensor.size();
	dim[1] = rTensor[0].rows();
	dim[2] = rTensor[0].cols();
	// data
	SEXP sMatrix = PROTECT(allocVector(REALSXP, dim[0] * dim[1] * dim[2]));
	double* data = REAL(sMatrix);
	int i = 0;
	for (int c = 0; c < dim[2]; ++c) {
		for (int r = 0; r < dim[1]; ++r) {
			for (vector<MatrixXd>::const_iterator it = rTensor.begin(); it != rTensor.end(); ++it) {
				data[i] = (*it)(r, c);
				++i;
			}
		}
	}
	setAttrib(sMatrix, R_DimSymbol, sDim);
	UNPROTECT(2); // sMatrix, sDim
	return sMatrix;
}

MatrixXd asMatrixXd(SEXP sMatrix) {
	SEXP sDim = getAttrib(sMatrix, R_DimSymbol);
	int* dim = INTEGER(sDim);
	MatrixXd matrix(dim[0], dim[1]);
	double* data = REAL(sMatrix);
	int i = 0;
	for (int c = 0; c < matrix.cols(); ++c) {
		for (int r = 0; r < matrix.rows(); ++r) {
			matrix(r, c) = data[i];
			++i;
		}
	}
	return matrix;
}

/**
 * Sets the i-th element of a list and modifies the names vector.
 *
 * @param sList The list.
 * @param i The numeric index of the element.
 * @param name The name of the element.
 * @param sElemnt The element.
 */
void setNamedListElt(SEXP sList, int i, const std::string& name, SEXP sElmnt) {
	SEXP sNames = PROTECT(getAttrib(sList, R_NamesSymbol));
	SET_STRING_ELT(sNames, i, mkChar(name.c_str()));
	SET_VECTOR_ELT(sList, i, sElmnt);
	setAttrib(sList, R_NamesSymbol, sNames);
	UNPROTECT(1);
}

/**
 * Get a element from a list by its name.
 *
 * @param sList The list.
 * @param name The name of the element.
 * @return The element.
 */
SEXP getNamedListElt(SEXP sList, const std::string& name) {
	const char* cname = name.c_str();
	SEXP sNames = getAttrib(sList, R_NamesSymbol);
	for (int i = 0; i < length(sList); ++i) {
		if (strcmp(CHAR(STRING_ELT(sNames, i)), cname) == 0) {
			return VECTOR_ELT(sList, i);
		}
	}
	return R_NilValue;
}

} // namespace siena
