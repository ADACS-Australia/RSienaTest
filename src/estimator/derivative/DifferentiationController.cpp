/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file DifferentiationController.cpp
 * \brief Implements the DifferentiationController class.
 *****************************************************************************/

#include "DifferentiationController.h"

#include "logger/Logger.h"
#include "ScoreDeviation.h"
#include "FiniteDifference.h"
#include "MLDifferentiation.h"

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

///////////////////////////////////////////////////////////////////////////////
// Finite differences
///////////////////////////////////////////////////////////////////////////////

//! Maximum number of retries with increased epsilon.
//! R: phase1.r phase1.1() (line ~94)
const int DifferentiationController::FD_REPEATS_FOR_EPSILON = 4;
//! Threshold to determine the needed gains. If the number of non-zero
//! elements is less or equal to this, the faster version is applied.
//! R: phase1.r phase1.1() (line ~81)
const double DifferentiationController::FD_EPSILON_LOWER_THRESHOLD =
		FiniteDifference::DERIVATIVE_CHECK_THRESHOLD / 2;
//! Epsilon multiplicator for basic rate parameters.
//! R: phase1.r phase1.1() (line ~83)
const double DifferentiationController::FD_EPSILON_GAIN_FAST_BASIC_RATE = 3;
//! Epsilon multiplicator for other effects.
//! R: phase1.r phase1.1() (line ~84)
const double DifferentiationController::FD_EPSILON_GAIN_FAST = 10;
//! Epsilon multiplicator for basic rate parameters.
//! R: phase1.r phase1.1() (line ~87)
const double DifferentiationController::FD_EPSILON_GAIN_SLOW_BASIC_RATE = 2;
//! Epsilon multiplicator for other effects.
//! R: phase1.r phase1.1() (line ~88)
const double DifferentiationController::FD_EPSILON_GAIN_SLOW = sqrt(10.0);

///////////////////////////////////////////////////////////////////////////////
// Parameter fixing
///////////////////////////////////////////////////////////////////////////////

//! Threshold of positive single derivatives below which the corresponding
//! parameter gets fixed.
//! R: phase1.r phase1.1() (line ~109)
const int DifferentiationController::FD_FIX_THRESHOLD = 2;

///////////////////////////////////////////////////////////////////////////////
// Score deviation
///////////////////////////////////////////////////////////////////////////////

//! Maximum number of iteration performed in phase 1 when using the score
//! deviation method.
//! R: phase1.r CalculateDerivative() (line ~334)
const int DifferentiationController::SD_MAX_ITERATIONS = 200;
//!
//! R: phase1.r CalculateDerivative() (line ~336)
const double DifferentiationController::SD_ITERATION_GAIN = 2;

///////////////////////////////////////////////////////////////////////////////
// Number of iterations
///////////////////////////////////////////////////////////////////////////////

//! Parameter for calculateN1()
//! R: robmon.r robmon() (line ~153)
const int DifferentiationController::N1_OFFSET = 7;
//! Parameter for calculateN1()
//! R: robmon.r robmon() (line ~153)
const int DifferentiationController::N1_FACTOR = 3;

/**
 * Calculates the preferred number of iterations for phase 1.
 *
 * Similar to: robmon.r robmon() (line ~153)
 *
 * @param minSimulations Minimum number of simulations.
 * @return Number of iterations.
 */
int DifferentiationController::calculateNIterations(
		const int minSimulations) {
	int n1 = N1_OFFSET + N1_FACTOR * lrSimulation.nStatistics();
	n1 = max(n1, minSimulations);
	return ceil((double) n1 / (double) lrSimulation.nSimulations());
}

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructs a DifferentiationController object.
 *
 * @param rSimulation The simulation.
 * @param rDiffType Differentiation type that should be used.
 * @param minSimulations Lower bound of number iterations to be performed.
 */
DifferentiationController::DifferentiationController(Simulation& rSimulation,
		const DifferentiationType& rDiffType, const int minSimulations=0) :
		Controller(rSimulation), //
		lDiffType(rDiffType), //
		lNIterations(calculateNIterations(minSimulations)), //
		lpDifferentiation(0/*ptr*/), //
		lFixed(VectorXb::Constant(rSimulation.nParameters(), false)) {
	if (lDiffType == SCORE_DEVIATION) {
		lpDifferentiation = new ScoreDeviation(lrSimulation.nParameters(),
					lrSimulation.nStatistics(), lrSimulation.nPeriods());
	} else {
		lpDifferentiation = new FiniteDifference(lrSimulation);
	}
}

/**
 * \copydoc Controller::run()
 */
bool DifferentiationController::run() {
	// If this is a ML simulation mlDerivative() is the only choice.
	if (lDiffType == MAXIMUM_LIKELIHOOD) {
		return mlDerivative();
	}
	// Otherwise use score deviation (if requested) with finite differences as
	// backup.
	return (lDiffType == SCORE_DEVIATION && scoreDeviation())
			|| finiteDifferences();
}

/**
 * @return Differentiation instance.
 */
Differentiation* DifferentiationController::pDifferentiation() {
	return lpDifferentiation;
}

/**
 * @return The derivative matrix.
 */
const MatrixXd& DifferentiationController::rDerivative() {
	return lpDifferentiation->rDerivative();
}

/**
 * @return Vector of parameters that are fixed because of the low derivative.
 */
const VectorXb& DifferentiationController::rFixed() {
	return lFixed;
}

/**
 * Error handling for the score deviation differentiation method.
 *
 * This extends the number of iterations performed in the phase till a good
 * derivative estimate has been reached.
 *
 * @return true if it has succeeded.
 */
bool DifferentiationController::scoreDeviation() {
	Differentiation* pDdiff = new ScoreDeviation(lrSimulation.nParameters(),
			lrSimulation.nStatistics(), lrSimulation.nPeriods());
	addListener(pDdiff);
	int n1 = lNIterations;
	int n1Check = ceil(
			pDdiff->simulationsBeforeCheck(n1 * lrSimulation.nSimulations())
					/ lrSimulation.nSimulations());
	LOGS(Priority::VERBOSE)<<"running score deviation in phase 1 with "<<n1
		<<" iterations, checking after "<<n1Check<<" iterations";
	// First try.
	fireInitializePhase();
	runPhase(n1Check);
	fireFinalizePhase(n1Check * lrSimulation.nSimulations());
	// Extend the number of iterations till we have a good result or run out of
	// iterations.
	// R: phase1.r CalculateDerivative() (line ~352)
	while (n1 * SD_ITERATION_GAIN < SD_MAX_ITERATIONS &&
			!pDdiff->checkDerivative()) {
		n1 *= SD_ITERATION_GAIN;
		LOGS(Priority::VERBOSE)<<"extend phase 1 to "<<n1<<" iterations";
		n1Check = ceil(
				pDdiff->simulationsBeforeCheck(n1 * lrSimulation.nSimulations())
				/ lrSimulation.nSimulations());
		// No fireInitializePhase() we want to extend the first run.
		runPhase(n1Check);
		fireFinalizePhase(n1Check * lrSimulation.nSimulations());
	}
	// Do the rest (if any).
	LOGS(Priority::VERBOSE)<<"perform the remaining "<<(n1-n1Check)<<" iterations";
	runPhase(n1 - n1Check);
	fireFinalizePhase(n1 * lrSimulation.nSimulations());
	removeListener(pDdiff);
	if (n1 < SD_MAX_ITERATIONS) {
		// We have got a good approximation with not too many iterations.
		// Make it public.
		delete lpDifferentiation;
		lpDifferentiation = pDdiff;
		return true;
	} else {
		// This took to long. Clean up.
		delete pDdiff;
		pDdiff = 0/*ptr*/;
		return false;
	}
}

/**
 * Error handling for the finite differences method.
 *
 * This runs a first a few simulations and checks if the single derivatives are
 * sufficient to approximate the derivative. If not the epsilon is increased and
 * the phase is restarted.
 *
 * @return true if it has succeeded.
 */
bool DifferentiationController::finiteDifferences() {
	LOGS(Priority::VERBOSE);
	FiniteDifference* pDdiff = new FiniteDifference(lrSimulation);
	addListener(pDdiff);
	const int n1 = lNIterations;
	const int n1Check = ceil(
			pDdiff->simulationsBeforeCheck(n1 * lrSimulation.nSimulations())
					/ lrSimulation.nSimulations());
	// First try.
	fireInitializePhase();
	runPhase(n1Check);
	fireFinalizePhase(n1Check * lrSimulation.nSimulations());
	// Restart the phase till we have a good estimate.
	// R: phase1.r phase1.1() (line ~103), robmon.r (line ~225)
	for (int r = 1; r <= FD_REPEATS_FOR_EPSILON && !pDdiff->checkDerivative();
			++r) {
		LOGS(Priority::VERBOSE)<<"Restart finite differences ("<<r
		<<"/"<<FD_REPEATS_FOR_EPSILON<<").";
		increaseEpsilon(pDdiff);
		fireInitializePhase();
		runPhase(n1Check);
		fireFinalizePhase(n1Check * lrSimulation.nSimulations());
	}
	// It may be that not all parameters have good values and we might need to
	// fix some.
	// R: phase1.r phase1.1() (line ~109)
	lFixed =
			lFixed.array()
					|| (pDdiff->rPositiveCount().colwise().maxCoeff().transpose().array()
							< FD_FIX_THRESHOLD);
	// Finish the rest of the iterations.
	runPhase(n1 - n1Check);
	fireFinalizePhase(n1 * lrSimulation.nSimulations());
	removeListener(pDdiff);
	delete lpDifferentiation;
	lpDifferentiation = pDdiff;
	// If all are fixed it is kind of pointless to run phase 2.
	return !lFixed.all();
}

/**
 * Increases epsilon according to the number of non-zero diagonal elements.
 *
 * Similar to: phase1.r phase1.1() (line ~82)
 *
 * @param pFiniteDifference FiniteDifference differentiation method.
 */
void DifferentiationController::increaseEpsilon(
		FiniteDifference* const pFiniteDifference) {
	VectorXd epsilon = pFiniteDifference->rEpsilon();
	const VectorXi maxPosCount =
			pFiniteDifference->rPositiveCount().colwise().maxCoeff();
	// Increase epsilon in each dimension
	for (int p = 0; p < epsilon.size(); ++p) {
		int count = maxPosCount[p];
		bool basicRate = lrSimulation.rSimulationEffects()[p]
				== &Simulation::RATE_EFFECT;
		if (count <= FD_EPSILON_LOWER_THRESHOLD) {
			// Faster for element with less information
			if (basicRate) {
				epsilon[p] *= FD_EPSILON_GAIN_FAST_BASIC_RATE;
			} else {
				epsilon[p] *= FD_EPSILON_GAIN_FAST;
			}
		} else if (count <= FiniteDifference::DERIVATIVE_CHECK_THRESHOLD) {
			// Slower for element with more information
			if (basicRate) {
				epsilon[p] *= FD_EPSILON_GAIN_SLOW_BASIC_RATE;
			} else {
				epsilon[p] *= FD_EPSILON_GAIN_SLOW;
			}
		}
	}
	pFiniteDifference->epsilon(epsilon);
	LOGS(Priority::VERBOSE)<<"Epsilon increased to: "<<epsilon.transpose();
}

/**
 * Error handling for the maximum likelihood method.
 *
 * @return true if it has succeeded.
 */
bool DifferentiationController::mlDerivative() {
	LOGS(Priority::DEBUG);
	delete lpDifferentiation;
	lpDifferentiation = new MLDifferentiation(lrSimulation.nParameters());
	addListener(lpDifferentiation);
	const int n1 = lNIterations;
//	const int n1Check = ceil(
//			lpDifferentiation->simulationsBeforeCheck(
//					n1 * lrSimulation.nSimulations())
//					/ lrSimulation.nSimulations());
	fireInitializePhase();
	runPhase(n1);
	fireFinalizePhase(n1 * lrSimulation.nSimulations());
	removeListener(lpDifferentiation);
	return true;
}

} // namespace siena
