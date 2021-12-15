/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MaxACStopCondition.h
 * \brief Defines the MaxACStopCondition class.
 *****************************************************************************/

#ifndef RSIENA_MAX_AC_STOP_CONDITION_H_
#define RSIENA_MAX_AC_STOP_CONDITION_H_

#include "estimator/update/StopCondition.h"
#include "AutoCorrelator.h"
#include <utility>

namespace siena {

/**
 * A stop condition terminating early if the maximum autocorrelation is low.
 */
class MaxACStopCondition: public StopCondition {
public:
	// Magic numbers
	static const double AC_MAX_EPS;

	MaxACStopCondition(AutoCorrelator& rAutoCorrelator,
			const std::pair<int, int>& nIterations);
	virtual ~MaxACStopCondition();

	virtual bool fulfilled(const int iteration);

protected:
	//! Autocorrelation of the statistics.
	//!
	AutoCorrelator& lrAC;
	//! Iteration bounds.
	//!
	const std::pair<int, int>& lNIterations;

private:
	// Don't copy it
	MaxACStopCondition(const MaxACStopCondition&);
	MaxACStopCondition& operator=(const MaxACStopCondition&);

};

/**
 * @param iteration Current iteration.
 * @return `true` if one of the following conditions is met.
 *   - `iteration` is equal or greater than the maximum iteration
 *   - `iteration` is equal or greater than the minimum iteration and the
 *     maximum autocorrelation is below AC_MAX_EPS
 */
inline bool MaxACStopCondition::fulfilled(const int iteration) {
	return (iteration >= lNIterations.second)
	// Stop early if the maximum autocorrelation is too low.
			|| (iteration >= lNIterations.first && lrAC.max() < AC_MAX_EPS);
}

} // namespace siena

#endif // RSIENA_MAX_AC_STOP_CONDITION_H_
