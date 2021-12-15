/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ScoreCollector.cpp
 * \brief Implements the ScoreCollector class.
 *****************************************************************************/

#include "ScoreCollector.h"

// #include <algorithm>

using namespace std;
using namespace Eigen;

namespace siena {

/**
 * Constructs the ScoreCollector object.
 */
ScoreCollector::ScoreCollector() :
		lScoresData(), //
		lScores() {
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with PERIOD_SCORES.
 */
vector<ResultType>& ScoreCollector::needs() {
	static vector<ResultType> needs(1, PERIOD_SCORES);
	return needs;
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void ScoreCollector::initializePhase() {
	lScoresData.clear();
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void ScoreCollector::finalizePhase(int nSimulations) {
	lScores = MatrixXdRM::Map(&lScoresData[0], nSimulations,
			lScoresData.size() / nSimulations);
}

/**
 * \copydoc EstimatorListener::pResultListener()
 */
ResultListener* ScoreCollector::pResultListener() {
	return this;
}

/**
 * @return The scores collected during simulation.
 */
const ScoreCollector::MatrixXdRM& ScoreCollector::rScores() {
	return lScores;
}

} // namespace siena
