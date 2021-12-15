/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ResultListener.h
 * \brief Defines the ResultListener class.
 *****************************************************************************/

#ifndef RSIENA_SIMULATION_LISTENER_H_
#define RSIENA_SIMULATION_LISTENER_H_

#include "Result.h"
#include "ResultType.h"
#include <vector>

namespace siena {

/**
 * The ResultListener catches the results of successful simulations.
 *
 * Results are distributed after the simulation finished a bunch of parallel
 * simulations.
 */
class ResultListener {
public:
	virtual ~ResultListener() = 0;

	/**
	 * @return Vector of simulation results needed by this listener. Anything
	 * not in the vector can not be assumed to have meaningful value in the
	 * result passed to onResults.
	 */
	virtual std::vector<ResultType>& needs();

	/**
	 * Called on every successful simulation.
	 *
	 * @param rResults Container providing access to the simulation results.
	 */
	virtual void onResults(const Result& rResults) = 0;

};

} // namespace siena

#endif // RSIENA_SIMULATION_LISTENER_H_
