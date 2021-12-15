/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file TimesCollector.h
 * \brief Defines the TimesCollector class.
 *****************************************************************************/

#ifndef RSIENA_TIMES_COLLECTOR_H_
#define RSIENA_TIMES_COLLECTOR_H_

#include "logger/Logger.h"
#include "sim/listener/ResultListener.h"
#include "estimator/listener/EstimatorListener.h"

#include <Eigen/Core>
#include <vector>

namespace siena {

/**
 * Component collecting the times of the conditional simulations used to
 * approximate the rate parameters.
 */
class TimesCollector: public ResultListener, public EstimatorListener {
public:

	//! Dynamic sized double matrix with RowMajor (used to map a sequentially
	//! filled std::vector).
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
			Eigen::RowMajor> MatrixXdRM;

	TimesCollector();

	// SimulationListener
	std::vector<ResultType>& needs();
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations);
	ResultListener* pResultListener();

	// TimesCollector
	const MatrixXdRM& rTimes();

protected:
	void addTimes(const Eigen::VectorXd& rTimes);
	//! Vector of scores, constantly appended.
	//!
	std::vector<double> lTimesData;
	//! Vector of scores mapped on finalizePhase().
	//!
	MatrixXdRM lTimes;

private:
	// Don't copy it
	TimesCollector(const TimesCollector&);
	TimesCollector& operator=(const TimesCollector&);

};

/**
 * Appends all score vectors of the results.
 *
 * Takes the mean score over all periods.
 *
 * @param rResults The simulation results.
 */
inline void TimesCollector::onResults(const Result& rResults) {
	for (std::vector<Eigen::Map<Eigen::VectorXd> >::const_iterator it =
			rResults.pTimes()->begin(); it != rResults.pTimes()->end(); ++it) {
		addTimes(*it);
	}
}

/**
 * Append a single time vector.
 *
 * @param rTimes The time to append.
 */
inline void TimesCollector::addTimes(const Eigen::VectorXd& rTimes) {
	lTimesData.insert(lTimesData.end(), rTimes.data(), rTimes.data()
			+ rTimes.size());
}

} // namespace siena

#endif // RSIENA_TIMES_COLLECTOR_H_
