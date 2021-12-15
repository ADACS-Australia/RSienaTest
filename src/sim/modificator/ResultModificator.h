/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ResultModificator.h
 * \brief Defines the ResultModificator class.
 *****************************************************************************/

#ifndef RSIENA_RESULT_MODIFICATOR_H_
#define RSIENA_RESULT_MODIFICATOR_H_

#include "sim/listener/Result.h"
#include "sim/listener/ResultType.h"
#include <vector>

namespace siena {

/**
 * A result listener with write access to the result set.
 *
 * The ResultModificator allows to modify the definition of an existing result
 * type. It is called before any ResultListener fires.
 */
class ResultModificator {
public:
	virtual ~ResultModificator() = 0;

	/**
	 * @return Vector of simulation results needed by this modificator. Anything
	 * not in the vector can not be assumed to have meaningful value in the
	 * result passed to onResults.
	 */
	virtual std::vector<ResultType>& needs();

	/**
	 * Called on every successful simulation.
	 *
	 * @param rResults Container providing access to the simulation results.
	 */
	virtual void onResults(Result& rResults) = 0;

};

} // namespace siena

#endif // RSIENA_RESULT_MODIFICATOR_H_
