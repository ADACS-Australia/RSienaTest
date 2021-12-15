/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ResultType.h
 * \brief Defines the ResultType enum.
 *****************************************************************************/

#ifndef RSIENA_RESULT_TYPE_H_
#define RSIENA_RESULT_TYPE_H_

namespace siena {

/**
 * Enumeration of result types produced by StatisticsSimulation
 * MetropolisHastingsSimulation.
 */
enum ResultType {
	STATISTICS,
	PERIOD_STATISTICS,
	MEAN_STATISTICS_MINUS_TARGET,
	PERIOD_SCORES,
	TIMES,
	ML_DERIVATIVE
};

} // namespace siena

#endif // RSIENA_RESULT_TYPE_H_
