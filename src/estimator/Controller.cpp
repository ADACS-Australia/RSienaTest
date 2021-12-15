/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Controller.cpp
 * \brief Implements the Controller class.
 *****************************************************************************/

#include "Controller.h"

#include "logger/Logger.h"
#include "estimator/listener/EstimatorListener.h"

#include <R_ext/Utils.h>

using namespace std;
using namespace siena::logger;

namespace siena {

/**
 * Constructs an Controller object.
 *
 * @param rSimulation Used simulation.
 */
Controller::Controller(Simulation& rSimulation) :
		lrSimulation(rSimulation), //
		lListenerPtrs() {
}

/**
 * Destroys an Controller object.
 */
Controller::~Controller() {
}

/**
 * Add the listener if not already present.
 *
 * @param pListener the listener.
 */
void Controller::addListener(EstimatorListener* const pListener) {
	LOGS(Priority::DEBUG);
	bool notFound = find(lListenerPtrs.begin(), lListenerPtrs.end(), pListener)
			== lListenerPtrs.end();
	if (notFound) {
		lListenerPtrs.push_back(pListener);
	}
	if (pListener->pResultListener()) {
		lrSimulation.addResultListener(pListener->pResultListener());
	}
}

/**
 * Remove the listener.
 *
 * @param pListener the listener.
 */
void Controller::removeListener(EstimatorListener* const pListener) {
	LOGS(Priority::DEBUG);
	list<EstimatorListener*>::iterator it = find(lListenerPtrs.begin(),
			lListenerPtrs.end(), pListener);
	if (it != lListenerPtrs.end()) {
		lListenerPtrs.erase(it);
	}
	if (pListener->pResultListener()) {
		lrSimulation.removeResultListener(pListener->pResultListener());
	}
}

/**
 * Fire an initializePhase event.
 */
void Controller::fireInitializePhase() {
	LOG(Priority::DEBUG, "");
	for (list<EstimatorListener*>::iterator it = lListenerPtrs.begin();
			it != lListenerPtrs.end(); ++it) {
		(*it)->initializePhase();
	}
}

/**
 * Fire an finalizePhase event.
 */
void Controller::fireFinalizePhase(int nSimulations) {
	LOG(Priority::DEBUG, "");
	for (list<EstimatorListener*>::iterator it = lListenerPtrs.begin();
			it != lListenerPtrs.end(); ++it) {
		(*it)->finalizePhase(nSimulations);
	}
}

/**
 * Requests simulations.
 *
 * @param nIterations Number of iterations each potentially consisting of
 *        multiple parallel simulations.
 */
void Controller::runPhase(const int nIterations) {
	for (int i = 0; i < nIterations; ++i) {
		if (i % 50 == 0) {
			LOGS(Priority::VERBOSE)<<"iteration: "<<i;
			R_CheckUserInterrupt();
		}
		lrSimulation.simulate();
	}
}

}
 // namespace siena
