/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ScoreCollector.h
 * \brief Defines the ScoreCollector class.
 *****************************************************************************/

#ifndef RSIENA_SCORE_COLLECTOR_H_
#define RSIENA_SCORE_COLLECTOR_H_

#include "logger/Logger.h"
#include "sim/listener/ResultListener.h"
#include "estimator/listener/EstimatorListener.h"

#include <Eigen/Core>
#include <vector>

namespace siena {

/**
 * Component collecting all simulated scores in a single matrix.
 */
class ScoreCollector: public ResultListener, public EstimatorListener {
public:

	//! Dynamic sized double matrix with RowMajor (used to map a sequentially
	//! filled std::vector).
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
			Eigen::RowMajor> MatrixXdRM;

	ScoreCollector();

	// SimulationListener
	std::vector<ResultType>& needs();
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations);
	ResultListener* pResultListener();

	// ScoreCollector
	const MatrixXdRM& rScores();

protected:
	void addScores(const Eigen::VectorXd& rScores);
	//! Vector of scores, constantly appended.
	//!
	std::vector<double> lScoresData;
	//! Vector of scores mapped on finalizePhase().
	//!
	MatrixXdRM lScores;

private:
	// Don't copy it
	ScoreCollector(const ScoreCollector&);
	ScoreCollector& operator=(const ScoreCollector&);

};

/**
 * Appends all score vectors of the results.
 *
 * Takes the mean score over all periods.
 *
 * @param rResults The simulation results.
 */
inline void ScoreCollector::onResults(const Result& rResults) {
	// LOGS(logger::Priority::DEBUG);
	for (std::vector<Eigen::Map<Eigen::MatrixXd> >::const_iterator it =
			rResults.pPeriodScores()->begin();
			it != rResults.pPeriodScores()->end(); ++it) {
		addScores((*it).colwise().mean());
	}
}

/**
 * Append a single score vector.
 *
 * @param rScores The vector to append.
 */
inline void ScoreCollector::addScores(const Eigen::VectorXd& rScores) {
 // Appends the Eigen::VectorXd to the std::vector.
	lScoresData.insert(lScoresData.end(), rScores.data(), rScores.data()
			+ rScores.size());
}

} // namespace siena

#endif // RSIENA_SCORE_COLLECTOR_H_
