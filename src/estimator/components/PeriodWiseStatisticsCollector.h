/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file PeriodWiseStatisticsCollector.h
 * \brief Defines the PeriodWiseStatisticsCollector class.
 *****************************************************************************/

#ifndef RSIENA_PERIOD_WISE_STATISTICS_COLLECTOR_H_
#define RSIENA_PERIOD_WISE_STATISTICS_COLLECTOR_H_

#include "logger/Logger.h"
#include "sim/listener/ResultListener.h"
#include "estimator/listener/EstimatorListener.h"

#include <Eigen/Core>
#include <vector>

namespace siena {

class PeriodWiseStatisticsCollector: public ResultListener, public EstimatorListener {
public:

	//! Dynamic sized double matrix with RowMajor (used to map a sequentially
	//! filled std::vector).
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
			Eigen::RowMajor> MatrixXdRM;

	PeriodWiseStatisticsCollector();

	// SimulationListener
	std::vector<ResultType>& needs();
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations) {};
	ResultListener* pResultListener();

	// PeriodWiseStatisticsCollector
	const std::vector<Eigen::MatrixXd>& rStatistics();

protected:
	std::vector<Eigen::MatrixXd> lStatistics;

private:
	// Don't copy it
	PeriodWiseStatisticsCollector(const PeriodWiseStatisticsCollector&);
	PeriodWiseStatisticsCollector& operator=(const PeriodWiseStatisticsCollector&);

};

inline void PeriodWiseStatisticsCollector::onResults(const Result& rResults) {
	// LOGS(logger::Priority::DEBUG);
	for (std::vector<Eigen::Map<Eigen::MatrixXd> >::const_iterator it =
			rResults.pPeriodStatistics()->begin();
			it != rResults.pPeriodStatistics()->end(); ++it) {
		lStatistics.push_back(*it);
	}
}

} // namespace siena

#endif // RSIENA_PERIOD_WISE_STATISTICS_COLLECTOR_H_
