/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ParameterSum.cpp
 * \brief Implements the ParameterSum class.
 *****************************************************************************/

#include "ParameterSum.h"

using namespace Eigen;

namespace siena {

/**
 * Constructs the ParameterSum object.
 *
 * @param nParameters Number of parameters.
 */
ParameterSum::ParameterSum(const int nParameters) :
		lrThetaSum(nParameters) {
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void ParameterSum::initializePhase() {
	lrThetaSum.fill(0);
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void ParameterSum::finalizePhase(int nSimulations) {
	//
}

/**
 * \copybrief EstimatorListener::pResultListener()
 *
 * @return this
 */
ResultListener* ParameterSum::pResultListener() {
	return this;
}

/**
 * @return The sum of parameters used in simulations since the last
 * initializePhase().
 */
const VectorXd& ParameterSum::rThetaSum() {
	return lrThetaSum;
}

} // namespace siena
