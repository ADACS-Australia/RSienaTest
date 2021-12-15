/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MeanStatisticsCalculator.h
 * \brief Defines the MeanStatisticsCalculator class.
 *****************************************************************************/

#ifndef RSIENA_MEAN_STATISTICS_MODIFICATOR_H_
#define RSIENA_MEAN_STATISTICS_MODIFICATOR_H_

#include "logger/Logger.h"
#include "sim/modificator/ResultModificator.h"
#include "Eigen/Util.h"

#include <Eigen/Core>
#include <vector>

namespace siena {

/**
 * Fills the mean statistics field of the results.
 */
class MeanStatisticsCalculator: public ResultModificator {
public:
	MeanStatisticsCalculator();

	std::vector<ResultType>& needs();
	void onResults(Result& rResults);

private:
	// Don't copy it
	MeanStatisticsCalculator(const MeanStatisticsCalculator&);
	MeanStatisticsCalculator& operator=(const MeanStatisticsCalculator&);

};

/**
 * Appends all statistics vectors of the results.
 *
 * @param rResults The simulation results.
 */
inline void MeanStatisticsCalculator::onResults(Result& rResults) {
	*rResults.pStatistics() = vColwiseSum(*rResults.pPeriodStatistics());
	*rResults.pMeanStatisticsMinusTargets() = vMean(*rResults.pStatistics()) - *rResults.pTargets();
}

} // namespace siena

#endif // RSIENA_MEAN_STATISTICS_MODIFICATOR_H_
