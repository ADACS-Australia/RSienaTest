/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file StatisticCollector.h
 * \brief Defines the StatisticCollector class.
 *****************************************************************************/

#ifndef RSIENA_STATISTIC_COLLECTOR_H_
#define RSIENA_STATISTIC_COLLECTOR_H_

#include "logger/Logger.h"
#include "sim/listener/ResultListener.h"
#include "estimator/listener/EstimatorListener.h"

#include <Eigen/Core>
#include <vector>

namespace siena {

/**
 * Component collecting all simulated statistics in a single matrix.
 */
class StatisticCollector: public ResultListener, public EstimatorListener {
public:

	//! Dynamic sized double matrix with RowMajor (used to map a sequentially
	//! filled std::vector).
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
			Eigen::RowMajor> MatrixXdRM;

	StatisticCollector();

	// SimulationListener
	std::vector<ResultType>& needs();
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations);
	ResultListener* pResultListener();

	// StatisticCollector
	const MatrixXdRM& rStatistics();

protected:
	void addStatistics(const Eigen::VectorXd& rStatistics);
	//! Number of simulations times number of statistics matrix with all
	//! simulated statistics.
	std::vector<double> lStatisticsData;
	//! Number of simulations times number of statistics matrix with all
	//! simulated statistics.
	MatrixXdRM lStatistics;

private:
	// Don't copy it
	StatisticCollector(const StatisticCollector&);
	StatisticCollector& operator=(const StatisticCollector&);

};

/**
 * Appends all statistics vectors of the results.
 *
 * @param rResults The simulation results.
 */
inline void StatisticCollector::onResults(const Result& rResults) {
	// LOGS(logger::Priority::DEBUG);
	for (std::vector<Eigen::VectorXd>::const_iterator it =
			rResults.pStatistics()->begin();
			it != rResults.pStatistics()->end(); ++it) {
		addStatistics(*it);
	}
}

/**
 * Append the single statistics vector.
 *
 * @param rStatistics The vector to append.
 */
inline void StatisticCollector::addStatistics(
		const Eigen::VectorXd& rStatistics) {
	LOGS(logger::Priority::DEBUG)<<"append statistics "<<rStatistics.transpose();
	// Appends the Eigen vector to the std::vector.
	lStatisticsData.insert(lStatisticsData.end(), rStatistics.data(),
			rStatistics.data() + rStatistics.size());
}

} // namespace siena

#endif // RSIENA_STATISTIC_COLLECTOR_H_
