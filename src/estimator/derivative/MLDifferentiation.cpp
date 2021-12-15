/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MLDifferentiation.cpp
 * Implementation of MLDifferentiation.h
 *****************************************************************************/

#include "MLDifferentiation.h"

#include "logger/Logger.h"

#include <vector>

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

/**
 * Constructs the MLDifferentiation object.
 *
 * @param nParameters Number of parameters.
 */
MLDifferentiation::MLDifferentiation(int nParameters) :
		Differentiation(), //
		lSumOfDerivative(nParameters, nParameters), //
		lDerivative() {
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with ML_DERIVATIVE.
 */
vector<ResultType>& MLDifferentiation::needs() {
	static vector<ResultType> needs(2, ML_DERIVATIVE);
	return needs;
}

/**
 * Updates the derivative with the new simulation results.
 *
 * @param rResults Simulation results.
 */
void MLDifferentiation::onResults(const Result& rResults) {
	for (vector<MatrixXd>::const_iterator itDerivs =
			rResults.pDerivatives()->begin();
			itDerivs != rResults.pDerivatives()->end(); ++itDerivs) {
		addDerivatives(*itDerivs);
	}
}

/**
 * Add a derivative approximation to the sum.
 *
 * @param rDerivative Approximated derivative.
 */
void MLDifferentiation::addDerivatives(const MatrixXd& rDerivative) {
	lSumOfDerivative += rDerivative;
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void MLDifferentiation::initializePhase() {
	lSumOfDerivative.fill(0);
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void MLDifferentiation::finalizePhase(int nSimulations) {
	lDerivative = lSumOfDerivative / nSimulations;
}

/**
 * \copydoc Differentiation::simulationsBeforeCheck()
 */
int MLDifferentiation::simulationsBeforeCheck(const int n) const {
	return n;
}

/**
 * \copydoc Differentiation::checkDerivative()
 */
bool MLDifferentiation::checkDerivative() const {
	return true;
}

/**
 * \copydoc Differentiation::rDerivative()
 */
const MatrixXd& MLDifferentiation::rDerivative() const {
	return lDerivative;
}

} // namespace siena
