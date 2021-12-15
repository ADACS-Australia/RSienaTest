/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file RUtil.h
 * \brief Utility functions mostly for converting C++ structures to R.
 *****************************************************************************/

#ifndef RUTIL_H_
#define RUTIL_H_

#include "Eigen/Types.h" // custom Eigen type declarations
#include "RInterface.h"

#include <Eigen/Core>

#include <vector>
#include <string>

namespace siena
{

///////////////////////////////////////////////////////////////////////////////
// RNG state
///////////////////////////////////////////////////////////////////////////////

void getSeed(std::vector<int>& seed);
void setSeed(const std::vector<int>& seed);

///////////////////////////////////////////////////////////////////////////////
// C/R conversion
///////////////////////////////////////////////////////////////////////////////

SEXP wrapPointer(void* cObject, void finalizer(SEXP) = 0/*ptr*/);

SEXP rifyString(const std::string& rString);
SEXP rifyIntVector(const std::vector<int>& rVector);
SEXP rifyVectorXb(const Eigen::VectorXb& rVector);
Eigen::VectorXb asVectorXb(SEXP sVector);
SEXP rifyVectorXd(const Eigen::VectorXd& rVector);
SEXP rifyMatrixXd(const Eigen::MatrixXd& rMatrix);
SEXP rifyVectorMatrixXd(const std::vector<Eigen::MatrixXd>& rTensor);
Eigen::MatrixXd asMatrixXd(SEXP sMatrix);

void setNamedListElt(SEXP sList, int i, const std::string& name, SEXP sElmnt);
SEXP getNamedListElt(SEXP sList, const std::string& name);

} // namespace siena

#endif /* RUTIL_H_ */
