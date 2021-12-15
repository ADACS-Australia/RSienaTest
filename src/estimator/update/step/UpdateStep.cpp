/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file UpdateStep.cpp
 * \brief Implements the UpdateStep class.
 *****************************************************************************/

#include "UpdateStep.h"

using namespace std;
using namespace Eigen;

namespace siena {

/**
 * Constructs a UpdateStep object.
 *
 * @param rSimulation The simulation updated in the step.
 * @param gain Initial gain.
 * @param rNormalization Normalization applied to the step.
 */
UpdateStep::UpdateStep(Simulation& rSimulation, const double gain,
		const StepNormalization& rNormalization) :
		lrSimulation(rSimulation), //
		lNormalization(rNormalization), //
		lGain(gain) {
}

UpdateStep::~UpdateStep() {
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with MEAN_STATISTICS_MINUS_TARGET.
 */
vector<ResultType>& UpdateStep::needs() {
	static vector<ResultType> needs(1, MEAN_STATISTICS_MINUS_TARGET);
	return needs;
}

/**
 * \copybrief ResultListener::onResults()
 *
 * Calls step() with the MEAN_STATISTICS_MINUS_TARGET results.
 *
 * @param rResults Results of the simulation.
 */
void UpdateStep::onResults(const Result& rResults) {
	step(*rResults.pMeanStatisticsMinusTargets());
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void UpdateStep::initializePhase() {
	//
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void UpdateStep::finalizePhase(int nSimulations) {
	//
}

/**
 * \copydoc EstimatorListener::pResultListener()
 */
ResultListener* UpdateStep::pResultListener() {
	return this;
}

/**
 * @param gain New step width gain.
 */
void UpdateStep::setGain(const double gain) {
	lGain = gain;
}

/**
 * @return Step width gain.
 */
double UpdateStep::getGain() const {
	return lGain;
}

} // namespace siena
