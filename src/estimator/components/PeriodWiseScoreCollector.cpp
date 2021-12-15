/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file PeriodWiseScoreCollector.cpp
 * \brief Implements the PeriodWiseScoreCollector class.
 *****************************************************************************/

#include "PeriodWiseScoreCollector.h"

using namespace std;
using namespace Eigen;

namespace siena {

PeriodWiseScoreCollector::PeriodWiseScoreCollector() :
		lScores() {
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with PERIOD_SCORES.
 */
vector<ResultType>& PeriodWiseScoreCollector::needs() {
	static vector<ResultType> needs(1, PERIOD_SCORES);
	return needs;
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void PeriodWiseScoreCollector::initializePhase() {
	lScores.clear();
}

/**
 * \copydoc EstimatorListener::pResultListener()
 */
ResultListener* PeriodWiseScoreCollector::pResultListener() {
	return this;
}

/**
 * @return The scores collected during simulation.
 */
const vector<MatrixXd>& PeriodWiseScoreCollector::rScores() {
	return lScores;
}

} // namespace siena
