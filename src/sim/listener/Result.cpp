/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Result.cpp
 * \brief Implements the Result class.
 *****************************************************************************/

#include "Result.h"

using namespace std;
using namespace Eigen;

namespace siena {

/**
 * Constructs a Result object.
 *
 * The parameters are pointers to the internal stores of the simulation.
 * Anything not provided by the simulation can be set to 0.
 *
 * @param pTheta Model parameters used.
 * @param pTargets Target vector.
 * @param pStatistics Mean statistics over all periods for each simulation.
 * @param pPeriodStatistics Main statistics.
 * @param pMeanStatisticsMinusTargets Mean over all simulations minus the
 *        target statistics.
 * @param pSeeds RNG state at the begin of each period.
 * @param pPeriodScores Scores per period, returned only from the
 *        EpochSimulation.
 * @param pDerivatives Derivatives, returned only from the MLSimulation.
 * @param pTimes Times of the conditional simulation.
 */
Result::Result(Eigen::VectorXd* pTheta, const Eigen::VectorXd* pTargets,
		const Eigen::MatrixXd* pPeriodTargets,
		vector<VectorXd>* pStatistics,
		vector<Map<MatrixXd> >* pPeriodStatistics,
		VectorXd* pMeanStatisticsMinusTargets, vector<vector<int> >* pSeeds,
		vector<Map<MatrixXd> >* pPeriodScores, vector<MatrixXd>* pDerivatives,
		vector<Map<VectorXd> >* pTimes) :
		lpTheta(pTheta), //
		lpTargets(pTargets), //
		lpPeriodTargets(pPeriodTargets), //
		lpStatistics(pStatistics), //
		lpPeriodStatistics(pPeriodStatistics), //
		lpMeanStatisticsMinusTargets(pMeanStatisticsMinusTargets), //
		lpSeeds(pSeeds), //
		lpPeriodScores(pPeriodScores), //
		lpDerivatives(pDerivatives), //
		lpTimes(pTimes) {
}

/**
 * @return Model parameters.
 */
const VectorXd* Result::pTheta() const {
	return lpTheta;
}

/**
 * @return Targets.
 */
const VectorXd* Result::pTargets() const {
	return lpTargets;
}

/**
 * @return Targets per period.
 */
const MatrixXd* Result::pPeriodTargets() const {
	return lpPeriodTargets;
}

/**
 * @return RNG state at the begin of each period.
 */
const vector<vector<int> >* Result::pSeeds() const {
	return lpSeeds;
}

/**
 * @return Main statistics over all periods.
 */
const vector<VectorXd>* Result::pStatistics() const {
	return lpStatistics;
}

/**
 * @return Main statistics over all periods.
 */
vector<VectorXd>* Result::pStatistics() {
	return lpStatistics;
}

/**
 * @return Main statistics per period.
 */
const vector<Map<MatrixXd> >* Result::pPeriodStatistics() const {
	return lpPeriodStatistics;
}

/**
 * @return Main statistics per period.
 */
vector<Map<MatrixXd> >* Result::pPeriodStatistics() {
	return lpPeriodStatistics;
}

/**
 * @return Mean of the main statistics minus the target statistics.
 */
const VectorXd* Result::pMeanStatisticsMinusTargets() const {
	return lpMeanStatisticsMinusTargets;
}

/**
 * @return Mean of the main statistics minus the target statistics.
 */
VectorXd* Result::pMeanStatisticsMinusTargets() {
	return lpMeanStatisticsMinusTargets;
}

/**
 * @return Scores per period, returned only from the EpochSimulation.
 */
const vector<Map<MatrixXd> >* Result::pPeriodScores() const {
	return lpPeriodScores;
}

/**
 * @return Scores per period, returned only from the EpochSimulation.
 */
vector<Map<MatrixXd> >* Result::pPeriodScores() {
	return lpPeriodScores;
}

/**
 * @return Derivatives, returned only from the MLSimulation,
 */
const vector<MatrixXd>* Result::pDerivatives() const {
	return lpDerivatives;
}

/**
 * @return Derivatives, returned only from the MLSimulation,
 */
vector<MatrixXd>* Result::pDerivatives() {
	return lpDerivatives;
}

/**
 * @return Times of the conditional simulation.
 */
const vector<Map<VectorXd> >* Result::pTimes() const {
	return lpTimes;
}

/**
 * @return Main statistics over all periods.
 */
vector<Map<VectorXd> >* Result::pTimes() {
	return lpTimes;
}

} // namespace siena
