/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ScoreDeviation.cpp
 * Implementation of ScoreDeviation.h
 *****************************************************************************/

#include "ScoreDeviation.h"

#include "logger/Logger.h"
#include "Eigen/Types.h"

#include <vector>

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

/**
 * Constructs the ScoreDeviation object.
 *
 * @param nParameters Number of parameters.
 * @param nStatistics Number of statistic effects.
 * @param nPeriods Number of periods.
 */
ScoreDeviation::ScoreDeviation(int nParameters, int nStatistics, int nPeriods) :
		Differentiation(), //
		lDerivative(), //
#ifdef R_LEGACY
				lSumStatsM(nPeriods, nStatistics), //
				lSumScoresM(nPeriods, nParameters)
#else
		lSumStatsM(nPeriods, nStatistics), //
		lSumScoresM(nPeriods, nParameters), //
#endif
		lSumOuter(nParameters, nStatistics) {
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with PERIOD_STATISTICS and PERIOD_SCORES.
 */
vector<ResultType>& ScoreDeviation::needs() {
	static ResultType needsc[] = {PERIOD_STATISTICS, PERIOD_SCORES};
	static vector<ResultType> needs(needsc, needsc+2);
	return needs;
}

/**
 * Updates the derivative with the new simulation results.
 *
 * @param rResults Simulation results.
 */
void ScoreDeviation::onResults(const Result& rResults) {
	LOGS(Priority::DEBUG);
	// We're interested in statistics and scores.
	vector<Map<MatrixXd> >::const_iterator itStats =
			rResults.pPeriodStatistics()->begin();
	vector<Map<MatrixXd> >::const_iterator itScores =
			rResults.pPeriodScores()->begin();
	// Loop over the individual simulations and add it to the derivative.
	while (itStats != rResults.pPeriodStatistics()->end()) {
		addScores(*itStats, *itScores);
		++itStats;
		++itScores;
	}
}

/**
 * Sums up the simulated scores, statistics and the outer product of both.
 *
 * Similar to: phase1.r derivativeFromScoresAndDeviations() (part 1)
 *
 * @param rStatistics Simulated statistics.
 * @param rScores Simulated scores.
 */
void ScoreDeviation::addScores(const Map<MatrixXd>& rStatistics,
		const Map<MatrixXd>& rScores) {
	// Sum over outer products of statistics and score in each period
	for (int m = 0; m < rStatistics.rows(); ++m) {
		lSumOuter += rScores.row(m).transpose() * rStatistics.row(m);
	}
	// Sum over both individually
#ifdef R_LEGACY
	lSumStatsM = lSumStatsM + rStatistics.cast<long double>();
	lSumScoresM = lSumScoresM + rScores.cast<long double>();
#else
	lSumStatsM += rStatistics;
	lSumScoresM += rScores;
#endif
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void ScoreDeviation::initializePhase() {
	lSumOuter.fill(0);
	lSumStatsM.fill(0);
	lSumScoresM.fill(0);
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void ScoreDeviation::finalizePhase(int nSimulations) {
	lDerivative = lSumOuter.transpose();
#ifdef R_LEGACY
	lDerivative = lDerivative / lNSimulations;
	MatrixXld meanStatsM = lSumStatsM / lNSimulations;
	MatrixXld meanScoresM = lSumScoresM / lNSimulations;
	// Subtract the outer products of the sums of statistics and score.
	MatrixXld tmp = MatrixXld::Zero(lrSim.nSimulationEffects(), lrSim.nSimulationEffects());
	for (int m = 0; m < lSumStatsM.rows(); ++m) {
		tmp += (meanStatsM.cast<double>().row(m).transpose()
				* meanScoresM.cast<double>().row(m)).cast<long double>();
	}
	lDerivative -= tmp.cast<double>();
#else
	lDerivative /= nSimulations;
	MatrixXd meanStatsM = lSumStatsM / nSimulations;
	MatrixXd meanScoresM = lSumScoresM / nSimulations;
	// Subtract the outer products of the sums of statistics and score.
	for (int m = 0; m < lSumStatsM.rows(); ++m) {
		lDerivative -= meanStatsM.row(m).transpose() * meanScoresM.row(m);
	}
#endif
}

/**
 * Calculates the number of simulations to be performed before the derivative
 * check.
 *
 * Score deviation performs the derivative check at the end of phase 1.
 *
 * @param n Total number of simulations to be performed in phase 1.
 * @return Total number of simulations to be performed in phase 1.
 */
int ScoreDeviation::simulationsBeforeCheck(const int n) const {
	return n;
}

/**
 * Checks if the collected single derivatives are enough for a reliable
 * approximation.
 *
 * Similar to: phase1.r PositiveDerivativeChecks() (line ~436)
 *
 * @return True if for all parameters there is one derivative that is
 *         positive.
 */
bool ScoreDeviation::checkDerivative() const {
	// Create a boolean matrix with 1 for positive derivatives.
	MatrixXb positive = MatrixXd::Constant(lDerivative.rows(),
			lDerivative.cols(), std::numeric_limits<double>::epsilon()).array()
			< lDerivative.array();
	// It used to be that only the diagonal has to be checked, but since GMM the
	// matrix is not square anymore, so we search for the best element in a
	// column. That is for each estimated parameter there is at least one
	// statistic that has an influence (but not every statistics has to
	// influence a parameter).
	return positive.colwise().any().all();
}

/**
 * Calculates the derivative.
 *
 * Similar to: phase1.r derivativeFromScoresAndDeviations() (part 2)
 *
 * @return Derivative matrix.
 */
const MatrixXd& ScoreDeviation::rDerivative() const {
	return lDerivative;
}

} // namespace siena
