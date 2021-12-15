/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Controller.h
 * \brief Defines the Controller class.
 *****************************************************************************/

#ifndef RSIENA_CONTROLLER_H_
#define RSIENA_CONTROLLER_H_

#include "sim/StatisticsSimulation.h"

#include <list>

namespace siena {

// Forward declarations
class EstimatorListener;

/**
 * Controls a phase goal of the estimation process.
 *
 * A Controller hold several components required to accomplish the overall
 * goal of the estimation phase. It requests simulations and restarts / runs
 * additional phases till the goal is accomplished. If this is not possible
 * estimate() returns false.
 *
 * Additional components might be added to listen to the Controller and its
 * Simulation, but they have no influence on the number of simulations or
 * subphases.
 */
class Controller {
public:
	explicit Controller(Simulation& rSimulation);
	virtual ~Controller();

	void addListener(EstimatorListener* const pListener);
	void removeListener(EstimatorListener* const pListener);

	/**
	 * Main routine of the Controller.
	 *
	 * @return True if the process finished successful, false otherwise.
	 */
	virtual bool run() = 0;

protected:
	void fireInitializePhase();
	void fireFinalizePhase(int nSimulations);

	void runPhase(const int nSimulations);

	//! Used simulation.
	//!
	Simulation& lrSimulation;

private:
	// Don't copy it
	Controller(const Controller&);
	Controller& operator=(const Controller&);

	//! Estimator listeners. Informed about the start and end of phases.
	//!
	std::list<EstimatorListener*> lListenerPtrs;

};

} // namespace siena

#endif // RSIENA_CONTROLLER_H_
