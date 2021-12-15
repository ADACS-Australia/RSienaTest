/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ScoreDeviation.h
 *
 * \class siena::ScoreDeviation
 * \brief Scores and deviation differentiation method.
 *****************************************************************************/

#ifndef RSIENA_SCORE_DEVIATION_H_
#define RSIENA_SCORE_DEVIATION_H_

#include "estimator/derivative/Differentiation.h"

namespace siena {

/**
 * Score deviation differentiation method.
 */
class ScoreDeviation: public Differentiation {
public:
	ScoreDeviation(int nParameters, int nStatistics, int nPeriods);

	// SimulationListener
	std::vector<ResultType>& needs();
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations);

	// ScoreDeviation
	int simulationsBeforeCheck(const int n) const;
	bool checkDerivative() const;
	const Eigen::MatrixXd& rDerivative() const;

protected:
	void addScores(const Eigen::Map<Eigen::MatrixXd>& rStatistics,
			const Eigen::Map<Eigen::MatrixXd>& rScores);
	//! The derivative.
	//!
	Eigen::MatrixXd lDerivative;
	//! Sum of period statistics.
	//!
#ifdef R_LEGACY
	Eigen::MatrixXld lSumStatsM;
#else
	Eigen::MatrixXd lSumStatsM;
#endif
	//! Sum of period scores.
	//!
#ifdef R_LEGACY
	Eigen::MatrixXld lSumScoresM;
#else
	Eigen::MatrixXd lSumScoresM;
#endif
	//! Sum of outer products of statistics and scores.
	//!
	Eigen::MatrixXd lSumOuter;

private:
	// Don't copy it
	ScoreDeviation(const ScoreDeviation&);
	ScoreDeviation& operator=(const ScoreDeviation&);

};

} // namespace siena

#endif // RSIENA_SCORE_DEVIATION_H_
