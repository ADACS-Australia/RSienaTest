/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file DolbyModificator.cpp
 * \brief Implements the DolbyModificator class.
 *****************************************************************************/

#include "DolbyModificator.h"

using namespace std;
using namespace Eigen;

namespace siena {

/**
 * Calculates the score regression coefficient.
 *
 * @param rFullStatistics Matrix with all simulated statistics (rowwise).
 * @param rFullScores Matrix with all simulated scores (rowwise).
 * @return Score regression coefficient.
 */
Eigen::VectorXd DolbyModificator::regressionCoefficient( const
		Eigen::MatrixXd& rFullStatistics, const Eigen::MatrixXd& rFullScores) {
	assert(rFullStatistics.cols() == rFullScores.cols());
	VectorXd reg = VectorXd(rFullStatistics.cols());
	for (int i = 0; i < rFullStatistics.cols(); ++i) {
		reg[i] = covariance(rFullStatistics.col(i), rFullScores.col(i)) / variance(rFullScores.col(i));
	}
	return reg;
}

/**
 * Constructs a DolbyModificator object.
 *
 * @param rFullStatistics Matrix with all simulated statistics (rowwise).
 * @param rFullScores Matrix with all simulated scores (rowwise).
 */
DolbyModificator::DolbyModificator(const Eigen::MatrixXd& rFullStatistics,
		const Eigen::MatrixXd& rFullScores) :
	lrRegressionCoefficient(regressionCoefficient(rFullStatistics,
				rFullScores)) {
}

/**
 * \copybrief ResultListener::needs()
 *
 * @return Vector with PERIOD_SCORES and MEAN_STATISTICS_MINUS_TARGET.
 */
vector<ResultType>& DolbyModificator::needs() {
	static ResultType needsc[] = {PERIOD_SCORES, MEAN_STATISTICS_MINUS_TARGET};
	static vector<ResultType> needs(needsc, needsc+2);
	return needs;
}

/**
 * @return The score regression coefficient.
 */
const Eigen::VectorXd& DolbyModificator::rRegressionCoefficient() {
	return lrRegressionCoefficient;
}

} // namespace siena
