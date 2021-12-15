/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MeanStatisticsCalculator.cpp
 * \brief Implements the MeanStatisticsCalculator class.
 *****************************************************************************/

#include "MeanStatisticsCalculator.h"

using namespace std;

namespace siena {

MeanStatisticsCalculator::MeanStatisticsCalculator() {
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with PERIOD_STATISTICS.
 */
vector<ResultType>& MeanStatisticsCalculator::needs() {
	static vector<ResultType> needs(PERIOD_STATISTICS);
	return needs;
}

} // namespace siena
