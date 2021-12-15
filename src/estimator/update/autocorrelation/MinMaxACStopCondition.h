/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MinMaxACStopCondition.h
 * \brief Defines the MinMaxACStopCondition class.
 *****************************************************************************/

#ifndef RSIENA_MIN_MAX_AC_STOP_CONDITION_H_
#define RSIENA_MIN_MAX_AC_STOP_CONDITION_H_

#include "MaxACStopCondition.h"

namespace siena {

/**
 * A stop condition terminating early if the minimum or maximum
 * autocorrelation is low.
 */
class MinMaxACStopCondition: public MaxACStopCondition {
public:
	// Magic numbers
	static const double AC_MIN_EPS;
	static const int AC_MIN_NITERATIONS;

	MinMaxACStopCondition(AutoCorrelator& rAutoCorrelator,
			const std::pair<int, int>& nIterations);

	bool fulfilled(const int iteration);

private:
	// Don't copy it
	MinMaxACStopCondition(const MinMaxACStopCondition&);
	MinMaxACStopCondition& operator=(const MinMaxACStopCondition&);

};

/**
 * @param iteration Current iteration.
 * @return `true` if one of the following conditions is met.
 *   - `iteration` is equal or greater than the maximum iteration
 *   - `iteration` is equal or greater than the minimum iteration and the
 *     maximum autocorrelation is below AC_MAX_EPS
 *   - `iteration` is equal or greater than AC_MIN_NITERATIONS and the minimum
 *     autocorrelation is below AC_MIN_EPS
 */
inline bool MinMaxACStopCondition::fulfilled(const int iteration) {
	return MaxACStopCondition::fulfilled(iteration)
	// Stop early if the minimum autocorrelation is too low.
			|| (iteration >= AC_MIN_NITERATIONS && lrAC.min() < AC_MIN_EPS);
}

} // namespace siena

#endif // RSIENA_MIN_MAX_AC_STOP_CONDITION_H_
