/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file UpdateController.cpp
 * \brief Implements the UpdateController class.
 *****************************************************************************/

#include "UpdateController.h"

#include "logger/Logger.h"
#include "autocorrelation/MinMaxACStopCondition.h"
#include "step/UpdateStep.h"

#include <R_ext/Utils.h>
//#include <cmath>

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

///////////////////////////////////////////////////////////////////////////////
// Number of sunphases
///////////////////////////////////////////////////////////////////////////////

//! Number of phase 2 subphases.
//! R: sienaModelCreate.r sienaModelCreate() (line ~23)
//const int UpdateController::N_SUBPHASES_DEFAULT = 4; // provided via R
//! Maximum number of repeats of one phase 2 subphase.
//! R: robmon.r robmon() (line ~27)
const int UpdateController::SUBPHASE_MAX_REPEATS = 4;

///////////////////////////////////////////////////////////////////////////////
// Number of iterations
///////////////////////////////////////////////////////////////////////////////

//! Parameter for calculateN2()
//! R: siena07.r AnnouncePhase() (line ~273)
const double UpdateController::N2_MIN_BASE = 2.52;
//! Parameter for calculateN2()
//! R: siena07.r AnnouncePhase() (line ~266)
const double UpdateController::N2_MIN_EFFECT_OFFSET = 7;
//! Parameter for calculateN2()
//! R: siena07.r AnnouncePhase() (line ~267)
const double UpdateController::N2_MIN_EFFECT_MIN = 5;
//! Parameter for calculateN2()
//! R: siena07.r AnnouncePhase() (line ~274)
const int UpdateController::N2_MAX_OFFSET = 200;

/**
 * Calculates the minimum and maximum number of iterations that can be
 * performed in subphases of phase 2.
 *
 * Similar to: siena07.r AnnouncePhase() (line ~266)
 *
 * @param subphase subphase-th sub-phase (0-based).
 * @return Range (pair of min and max) of valid iteration counts for subphases of phase 2.
 */
std::pair<int, int> UpdateController::calculateNIterations(const int subphase) {
	// R: n2min0 is a double, therefore first multiplication here
	int n2min = max(N2_MIN_EFFECT_MIN,
			(N2_MIN_EFFECT_OFFSET + lrSimulation.nStatistics())
					/ lrSimulation.nSimulations()) * N2_MIN_BASE;
	// Multiply with 2.52 for each subphase
	for (int i = 1; i <= subphase; ++i) {
		n2min *= N2_MIN_BASE; // R: Truncates to integer
	}
	int n2max = n2min + N2_MAX_OFFSET;
	return pair<int, int>(n2min, n2max);
}

///////////////////////////////////////////////////////////////////////////////
// Gain control
///////////////////////////////////////////////////////////////////////////////

//! Initial step size parameter. Common name: firstg
//! R: sienaModelCreate.r sienaModelCreate() (line ~26)
//const double UpdateController::GAIN_INITIAL_DEFAULT = .2;
//! Factor multiplied to the gain after each sub-phase of phase 2. // provided via R
//! Common name: reduceg
//! R: sienaModelCreate.r sienaModelCreate() (line ~26)
//const double UpdateController::GAIN_DECAY_DEFAULT = .5; // provided via R

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructs a new UpdateController object.
 *
 * @param rSimulation The simulation producing the statistics for the update.
 * @param nSubphases Number of subphases to perform.
 * @param gainInit Initial step width gain.
 * @param gainDecay Gain decay factor, applied after each subphase.
 */
UpdateController::UpdateController(Simulation& rSimulation, int nSubphases,
		double gainInit, double gainDecay) :
		Controller(rSimulation), //
		lNSubphases(nSubphases), //
		lGainInitial(gainInit * sqrt(rSimulation.nSimulations())), //
		lGainDecay(gainDecay), //
		lNIterations(0), //
		lpStep(0/*ptr*/), //
		lAutoCorrelator(rSimulation.nStatistics()), //
		lParameterSum(rSimulation.nParameters()) {
}

UpdateController::~UpdateController() {
}

/**
 * Issues simulations till a stop condition is fulfilled.
 *
 * @param rCond The stop condition.
 */
void UpdateController::runSubphase(StopCondition &rCond) {
	LOGS(Priority::INFO)<<"(re)started subphase";
	fireInitializePhase();
	lNIterations = 0;
	while (!rCond.fulfilled(lNIterations)) {
		if (lNIterations % 50 == 0) {
			LOGS(Priority::VERBOSE)<<"iteration: "<<lNIterations;
			R_CheckUserInterrupt();
		}
		lrSimulation.simulate();
		++lNIterations;
	}
	fireFinalizePhase(lNIterations * lrSimulation.nSimulations());
}

/**
 * Repeats the subphase till the MaxACStopCondition is fulfilled or the number
 * of retries exceeds the maximum.
 *
 * @param nIterations Lower and upper iteration bound.
 */
void UpdateController::retryPhase(const pair<int, int>& nIterations) {
	// For saving parameters at the start of the phase, they count to the sum.
	VectorXd thetaFirst(lrSimulation.nParameters());
	// Reset the AC to no fail the first for-loop test.
	lAutoCorrelator.initializePhase();
	MinMaxACStopCondition minMaxCond(lAutoCorrelator, nIterations);
	MaxACStopCondition maxCond(lAutoCorrelator, nIterations);
	// Loop starting at 1 (compare R variable z$repeatsubphase)
	// R: phase2.r proc2subphase (line ~80)
	// Condition for repeating is always the maximum ac condition.
	// R: phase2.r proc2subphase (line ~137)
	for (int lPhase2Repeat = 1;
			lPhase2Repeat < SUBPHASE_MAX_REPEATS
					&& !maxCond.fulfilled(lNIterations); ++lPhase2Repeat) {
		thetaFirst = lrSimulation.rParameters();
		// All but the last retry uses the minimum/maximum condition.
		// R: phase2.r doIterations (line ~402)
		runSubphase(minMaxCond);
		if (lAutoCorrelator.max() >= sqrt(2.0 / (lNIterations + 1))) {
			LOG(Priority::INFO, "positive autocorrelation");
		}
	}
	// If maximum condition is still not satisfied perform one last run.
	// R: phase2.r proc2subphase (line ~137)
	if (!maxCond.fulfilled(lNIterations)) {
		thetaFirst = lrSimulation.rParameters();
		runSubphase(maxCond);
	}
	// Use average theta as new theta.
	LOGS(Priority::VERBOSE)<<"sum of parameters: "
	<<lParameterSum.rThetaSum().transpose();
	lrSimulation.updateParameters(
			(thetaFirst + lParameterSum.rThetaSum()) / (1 + lNIterations));
	LOGS(Priority::VERBOSE)<<"parameters after subphase: "
	<<lrSimulation.rParameters().transpose();
}

/**
 * \copydoc Controller::run()
 */
bool UpdateController::run() {
	addListener(lpStep);
	addListener(&lAutoCorrelator);
	addListener(&lParameterSum);
	// Run lPhase2N2subphases sub phases
	lpStep->setGain(lGainInitial);
	for (int subphase = 0; subphase < lNSubphases; ++subphase) {
		const pair<int, int> n = calculateNIterations(subphase);
		LOGS(Priority::INFO)<<"\n\n=== start subphase ==========="
		<<"\nSubphase: "<<subphase
		<<"\nIterations: ["<<n.first<<".."<<n.second<<"]"
		<<"\n==============================\n";
		retryPhase(n);
		// Reduce gain
		lpStep->setGain(lpStep->getGain() * lGainDecay);
		LOGS(Priority::VERBOSE)<<"reduced gain to: "<<lpStep->getGain();
	}
	removeListener(&lParameterSum);
	removeListener(&lAutoCorrelator);
	removeListener(lpStep);
	return true;
}

/**
 * Sets the update step used to update the model.
 *
 * @param pStep The update step.
 */
void UpdateController::step(UpdateStep* pStep) {
	lpStep = pStep;
}

/**
 * @return Number of subphases.
 */
int UpdateController::nSubphases() {
	return lNSubphases;
}

/**
 * @return Initial gain.
 */
double UpdateController::gainInitial() {
	return lGainInitial;
}

/**
 * @return Gain decay.
 */
double UpdateController::gainDecay() {
	return lGainDecay;
}

} // namespace siena
