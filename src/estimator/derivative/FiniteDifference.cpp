/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file FiniteDifference.cpp
 * \brief Defines the FiniteDifference class.
 *****************************************************************************/

#include "FiniteDifference.h"

#include <logger/Logger.h>
#include <logger/Priority.h>
#include <model/EffectInfo.h>
#include <sim/listener/Result.h>

#include <Eigen/Core>
#include <cmath>
#include <limits>
#include <ctime>

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

///////////////////////////////////////////////////////////////////////////////
// Espilon
///////////////////////////////////////////////////////////////////////////////

//! Initial epsilon.
//! R: robmon.r robmon() (line ~145)
const double FiniteDifference::INITIAL_EPSILON = .1;

///////////////////////////////////////////////////////////////////////////////
// Number of iterations
///////////////////////////////////////////////////////////////////////////////

//! Number of simulations before the derivative check.
//! R: phase1.r phase1.1() (line ~30)
const int FiniteDifference::SIMULATIONS_BEFORE_CHECK = 10;

///////////////////////////////////////////////////////////////////////////////
// Check
///////////////////////////////////////////////////////////////////////////////

//! Number of single derivative that have to be non-zero to pass the check.
//! R: phase1.r phase1.1() (line ~71)
const int FiniteDifference::DERIVATIVE_CHECK_THRESHOLD =
		SIMULATIONS_BEFORE_CHECK / 2;

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructs a FiniteDifference object.
 *
 * @param rSim Simulation object.
 */
FiniteDifference::FiniteDifference(Simulation& rSim) :
		Differentiation(), //
		lrSimuatlion(rSim), //
		lTheta(rSim.rParameters()), //
		lSeeds(), //
		lEpsilon(VectorXd::Constant(lTheta.size(), INITIAL_EPSILON)), //
		lParameter(0), //
		lStatistics(), //
		lSumOfDerivatives(rSim.nStatistics(), rSim.nParameters()), //
		lDerivative(rSim.nStatistics(), rSim.nParameters()), //
		lPositiveCount(rSim.nStatistics(), rSim.nParameters()), //
		lPrivateSim(this) {
	// Epsilon parameter for basic rates are scaled (up) according to the rate.
	// R: robmon.r robmon() (line ~147)
	lEpsilon =
			(rSim.rSimulationEffects().array() == &Simulation::RATE_EFFECT).select(
					lEpsilon.cwiseProduct(lTheta), lEpsilon);
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with STATISTICS.
 */
vector<ResultType>& FiniteDifference::needs() {
	static vector<ResultType> needs(1, STATISTICS);
	return needs;
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void FiniteDifference::initializePhase() {
	lSumOfDerivatives.fill(0);
	lPositiveCount.fill(0);
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void FiniteDifference::finalizePhase(int nSimulations) {
	lDerivative = lSumOfDerivatives / nSimulations;
}

/**
 * Performs one additional simulation for each dimension to approximate the
 * derivative.
 *
 * Similar to: phase1.r FiniteDifferences()
 *
 * @param rResults Simulation results.
 */
void FiniteDifference::onResults(const Result& rResults) {
	// Remember 'base' seeds and statistics
	lSeeds = *rResults.pSeeds();
	lStatistics = *rResults.pStatistics();
	// Prepare private simulations
	lTheta = lrSimuatlion.rParameters();
	vector<ResultListener*> publicListeners = lrSimuatlion.clearResultListener();
	lrSimuatlion.addResultListener(&lPrivateSim);

	// For each parameter, request a new simulation. The results are catched only
	// in the private simulation listener, that is via onEpsilonResults().
	LOG(Priority::DEBUG, "request private simulations for finite difference");
	for (lParameter = 0; lParameter < lTheta.size(); ++lParameter) {
		lTheta[lParameter] += lEpsilon[lParameter];
		lrSimuatlion.updateParameters(lTheta);
		lrSimuatlion.simulate(true);
		lTheta[lParameter] -= lEpsilon[lParameter];
	}
	// Reset simulation
	lrSimuatlion.clearResultListener();
	for (vector<ResultListener*>::iterator it = publicListeners.begin();
			it != publicListeners.end(); ++it) {
		lrSimuatlion.addResultListener(*it);
	}
	lrSimuatlion.updateParameters(lTheta);
	LOG(Priority::DEBUG, "again normal simulations");
}

/**
 * Update the derivative with private simulations done for each parameter
 * individually.
 *
 * @param rResults Simulation results.
 */
void FiniteDifference::onEpsilonResults(const Result& rResults) {
	LOG(Priority::DEBUG, "");
	// For each simulation in the chunk.
	for (int i = 0; i < lrSimuatlion.nSimulations(); ++i) {
		// Sum up the differences.
		const VectorXd delta = ((*rResults.pStatistics())[i] - lStatistics[i])
				/ lEpsilon[lParameter];
		lSumOfDerivatives.col(lParameter) += delta;
		// Count positive elements.
		const VectorXi count = (delta.array()
				> std::numeric_limits<double>::epsilon()).cast<int>();
		lPositiveCount.col(lParameter) += count;
	}
}

/**
 * Calculates the number of simulations to be performed before the derivative
 * check.
 *
 * Finite difference performs an early derivative check after
 * SIMULATIONS_BEFORE_CHECK simulations.
 *
 * Similar to: phase1.r phase1.1() (line ~30)
 *
 * Common name: firstNit
 *
 * @param n Total number of simulations to be performed in phase 1.
 * @return Simulations before the derivative check in phase 1.
 */
int FiniteDifference::simulationsBeforeCheck(const int) const {
	return SIMULATIONS_BEFORE_CHECK;
}

/**
 * Checks if the collected single derivatives are enough for a reliable approximation.
 *
 * Similar to: phase1.r FiniteDifferenceChecks() (line ~436)
 *
 * @return true if enough derivative were positive.
 */
bool FiniteDifference::checkDerivative() const {
	LOGS(logger::Priority::VERBOSE)<<"Number of positive derivatives:\n"
	<<lPositiveCount<<"\n";
	return (lPositiveCount.colwise().maxCoeff().array()
			>= DERIVATIVE_CHECK_THRESHOLD).all();
}

/**
 * \copydoc Differentiation::rDerivative()
 */
const MatrixXd& FiniteDifference::rDerivative() const {
	return lDerivative;
}

/**
 * @return Matrix counting how many times a single derivative approximation
 *         was positive.
 */
const MatrixXi& FiniteDifference::rPositiveCount() const {
	return lPositiveCount;
}

/**
 * @param rEpsilon New epsilon.
 */
void FiniteDifference::epsilon(const VectorXd& rEpsilon) {
	lEpsilon = rEpsilon;
}

/**
 * @return Epsilon.
 */
const VectorXd& FiniteDifference::rEpsilon() const {
	return lEpsilon;
}

} // namespace siena
