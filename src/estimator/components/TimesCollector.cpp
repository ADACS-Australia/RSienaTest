/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file TimesCollector.cpp
 * \brief Implements the TimesCollector class.
 *****************************************************************************/

#include "TimesCollector.h"

using namespace std;
using namespace Eigen;

namespace siena {

/**
 * Constructs the TimesCollector object.
 */
TimesCollector::TimesCollector() :
		lTimesData(), //
		lTimes() {
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with TIMES.
 */
vector<ResultType>& TimesCollector::needs() {
	static vector<ResultType> needs(1, TIMES);
	return needs;
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void TimesCollector::initializePhase() {
	lTimesData.clear();
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void TimesCollector::finalizePhase(int nSimulations) {
	lTimes = MatrixXdRM::Map(&lTimesData[0], nSimulations,
			lTimesData.size() / nSimulations);
}

/**
 * \copydoc EstimatorListener::pResultListener()
 */
ResultListener* TimesCollector::pResultListener() {
	return this;
}

/**
 * @return The times collected during simulation.
 */
const TimesCollector::MatrixXdRM& TimesCollector::rTimes() {
	return lTimes;
}

} // namespace siena
