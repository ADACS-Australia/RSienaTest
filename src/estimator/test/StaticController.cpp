/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file StaticController.cpp
 * \brief Implements the StaticController class.
 *****************************************************************************/

#include "StaticController.h"

#include "logger/Logger.h"

using namespace std;
using namespace siena::logger;

namespace siena {

/**
 * Constructs an StaticController object.
 *
 * @param rSimulation Used simulation.
 * @param nSimulations Number of simulations.
 */
StaticController::StaticController(Simulation& rSimulation,
		const int nSimulations) :
		Controller(rSimulation), //
		lNIterations(ceil(nSimulations / rSimulation.nSimulations())) {
}

/**
 * Destroys an StaticController object.
 */
StaticController::~StaticController() {
	//
}

/**
 * \copydoc Controller::run()
 */
bool StaticController::run() {
	fireInitializePhase();
	runPhase(lNIterations);
	fireFinalizePhase(lNIterations * lrSimulation.nSimulations());
  return true;
}

} // namespace siena
