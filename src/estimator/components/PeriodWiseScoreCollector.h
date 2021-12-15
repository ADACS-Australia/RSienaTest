/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file PeriodWiseScoreCollector.h
 * \brief Defines the PeriodWiseScoreCollector class.
 *****************************************************************************/

#ifndef RSIENA_PERIOD_WISE_SCORE_COLLECTOR_H_
#define RSIENA_PERIOD_WISE_SCORE_COLLECTOR_H_

#include "logger/Logger.h"
#include "sim/listener/ResultListener.h"
#include "estimator/listener/EstimatorListener.h"

#include <Eigen/Core>
#include <vector>

namespace siena {

/**
 * Component collecting all simulated scores in a single matrix.
 */
class PeriodWiseScoreCollector: public ResultListener, public EstimatorListener {
public:

	//! Dynamic sized double matrix with RowMajor (used to map a sequentially
	//! filled std::vector).
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
			Eigen::RowMajor> MatrixXdRM;

	PeriodWiseScoreCollector();

	// SimulationListener
	std::vector<ResultType>& needs();
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations) {};
	ResultListener* pResultListener();

	// PeriodWiseScoreCollector
	const std::vector<Eigen::MatrixXd>& rScores();

protected:
	//! Vector of scores mapped on finalizePhase().
	//!
	std::vector<Eigen::MatrixXd> lScores;

private:
	// Don't copy it
	PeriodWiseScoreCollector(const PeriodWiseScoreCollector&);
	PeriodWiseScoreCollector& operator=(const PeriodWiseScoreCollector&);

};

/**
 * Appends all score vectors of the results.
 *
 * Takes the mean score over all periods.
 *
 * @param rResults The simulation results.
 */
inline void PeriodWiseScoreCollector::onResults(const Result& rResults) {
	// LOGS(logger::Priority::DEBUG);
	for (std::vector<Eigen::Map<Eigen::MatrixXd> >::const_iterator it =
			rResults.pPeriodScores()->begin();
			it != rResults.pPeriodScores()->end(); ++it) {
		lScores.push_back(*it);
	}
}

} // namespace siena

#endif // RSIENA_PERIOD_WISE_SCORE_COLLECTOR_H_
