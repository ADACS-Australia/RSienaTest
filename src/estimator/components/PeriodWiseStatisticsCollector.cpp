/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file PeriodWiseStatisticsCollector.cpp
 * \brief Implements the PeriodWiseStatisticsCollector class.
 *****************************************************************************/

#include "PeriodWiseStatisticsCollector.h"

using namespace std;
using namespace Eigen;

namespace siena {

PeriodWiseStatisticsCollector::PeriodWiseStatisticsCollector() :
		lStatistics() {
}

vector<ResultType>& PeriodWiseStatisticsCollector::needs() {
	static vector<ResultType> needs(1, PERIOD_STATISTICS);
	return needs;
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void PeriodWiseStatisticsCollector::initializePhase() {
	lStatistics.clear();
}

/**
 * \copydoc EstimatorListener::pResultListener()
 */
ResultListener* PeriodWiseStatisticsCollector::pResultListener() {
	return this;
}

const vector<MatrixXd>& PeriodWiseStatisticsCollector::rStatistics() {
	return lStatistics;
}

} // namespace siena
