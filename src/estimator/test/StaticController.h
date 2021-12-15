/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file StaticController.h
 * \brief Defines the StaticController class.
 *****************************************************************************/

#ifndef RSIENA_STATIC_StaticController_H_
#define RSIENA_STATIC_StaticController_H_

#include "estimator/Controller.h"

namespace siena {

/**
 * A controller that runs a fixed number of iterations.
 *
 * This can be use for requesting a fixed number of simulations for tests
 * (phase 3).
 */
class StaticController: public Controller {
public:
	StaticController(Simulation& rSimulation, const int nSimulations);
	~StaticController();

	bool run();

protected:
	//! Iterations to be performed.
	//!
	const int lNIterations;

private:
	// Don't copy it
	StaticController(const StaticController&);
	StaticController& operator=(const StaticController&);

};

} // namespace siena

#endif // RSIENA_STATIC_StaticController_H_
