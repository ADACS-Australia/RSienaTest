/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file DolbyModificator.h
 * \brief Defines the DolbyModificator class.
 *****************************************************************************/

#ifndef RSIENA_DOLBY_MODIFICATOR_H_
#define RSIENA_DOLBY_MODIFICATOR_H_

#include "logger/Logger.h"
#include "sim/modificator/ResultModificator.h"
#include "Eigen/Util.h"

#include <Eigen/Core>
#include <vector>

namespace siena {

/**
 * Implements the 'dolby' option.
 */
class DolbyModificator: public ResultModificator {
public:
	DolbyModificator(const Eigen::MatrixXd& rFullStatistics,
			const Eigen::MatrixXd& rFullScores);

	std::vector<ResultType>& needs();
	void onResults(Result& rResults);

	const Eigen::VectorXd& rRegressionCoefficient();

protected:
	//! The score regression calculated in phase 1.
	//!
	const Eigen::VectorXd lrRegressionCoefficient;

	static Eigen::VectorXd regressionCoefficient(const Eigen::MatrixXd& rFullStatistics,
			const Eigen::MatrixXd& rFullScores);

private:
	// Don't copy it
	DolbyModificator(const DolbyModificator&);
	DolbyModificator& operator=(const DolbyModificator&);

};

/**
 * Appends all statistics vectors of the results.
 *
 * @param rResults The simulation results.
 */
inline void DolbyModificator::onResults(Result& rResults) {
	const std::vector<Eigen::VectorXd>& scores = vColwiseSum(*rResults.pPeriodScores());
	for (std::vector<Eigen::VectorXd>::size_type i = 0; i < scores.size(); ++i) {
		*rResults.pMeanStatisticsMinusTargets() -=
			lrRegressionCoefficient.cwiseProduct(scores[i]) / scores.size();
	}
}

} // namespace siena

#endif // RSIENA_DOLBY_MODIFICATOR_H_
