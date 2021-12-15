/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file StatisticCollector.cpp
 * \brief Implements the StatisticCollector class.
 *****************************************************************************/

#include "StatisticCollector.h"

#include <algorithm>

using namespace std;
using namespace Eigen;

namespace siena {

/**
 * Constructs the StatisticCollector object.
 */
StatisticCollector::StatisticCollector() :
		lStatisticsData(), //
		lStatistics() {
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with STATISTICS.
 */
vector<ResultType>& StatisticCollector::needs() {
	static vector<ResultType> needs(1, STATISTICS);
	return needs;
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void StatisticCollector::initializePhase() {
	lStatisticsData.clear();
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void StatisticCollector::finalizePhase(int nSimulations) {
	lStatistics = MatrixXdRM::Map(&lStatisticsData[0], nSimulations,
			lStatisticsData.size() / nSimulations);
}

/**
 * \copydoc EstimatorListener::pResultListener()
 */
ResultListener* StatisticCollector::pResultListener() {
	return this;
}

/**
 * @return The simulated statistics.
 */
const StatisticCollector::MatrixXdRM& StatisticCollector::rStatistics() {
	return lStatistics;
}

} // namespace siena
