/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file EstimatorListener.h
 * \brief Defines the EstimatorListener class.
 *****************************************************************************/

#ifndef RSIENA_ESTIMATOR_LISTENER_H_
#define RSIENA_ESTIMATOR_LISTENER_H_

#include <Eigen/Core>

namespace siena {

// Forward declarations
class ResultListener;

/**
 * A listener catching events form the controller.
 */
class EstimatorListener {
public:
	virtual ~EstimatorListener() = 0;

	/**
	 * Called when the Controller starts a new phase.
	 */
	virtual void initializePhase() = 0;

	/**
	 * Called when the Controller ends a phase.
	 *
	 * Any calculation needed at the end of the phase should be performed.
	 *
	 * If the controller chooses to extend the phases there might still be
	 * simulations coming in after this call, followed by an addition call to
	 * finalizePhase().
	 *
	 * @param nSimulations Number of simulations performed in the last phase.
	 */
	virtual void finalizePhase(int nSimulations) = 0;

	/**
	 * ResultListener to register when added to the controller.
	 *
	 * @return A ResultListener to be registered on the simulation of the
	 *         controller or 0 if none.
	 */
	virtual ResultListener* pResultListener() = 0;

};

} // namespace siena

#endif // RSIENA_ESTIMATOR_LISTENER_H_
