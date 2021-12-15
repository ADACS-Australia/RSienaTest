/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Result.h
 * \brief Defines the Result class.
 *****************************************************************************/

#ifndef RSIENA_RESULT_H_
#define RSIENA_RESULT_H_

#include <Eigen/Core>
#include <vector>

namespace siena {

/**
 * Class providing access to results of a bunch of simulations.
 */
class Result {
public:
	Result(Eigen::VectorXd* pTheta, const Eigen::VectorXd* pTargets,
			const Eigen::MatrixXd* pPeriodTargets,
			std::vector<Eigen::VectorXd>* pStatistics,
			std::vector<Eigen::Map<Eigen::MatrixXd> >* pPeriodStatistics,
			Eigen::VectorXd* pMeanStatisticsMinusTargets,
			std::vector<std::vector<int> >* pSeeds,
			std::vector<Eigen::Map<Eigen::MatrixXd> >* pPeriodScores,
			std::vector<Eigen::MatrixXd>* pDerivatives,
			std::vector<Eigen::Map<Eigen::VectorXd> >* pTimes);

	const Eigen::VectorXd* pTheta() const;
	const Eigen::VectorXd* pTargets() const;
	const std::vector<std::vector<int> >* pSeeds() const;

	const std::vector<Eigen::VectorXd>* pStatistics() const;
	std::vector<Eigen::VectorXd>* pStatistics();

	const Eigen::MatrixXd* pPeriodTargets() const;
	Eigen::MatrixXd* pPeriodTargets();

	const std::vector<Eigen::Map<Eigen::MatrixXd> >* pPeriodStatistics() const;
	std::vector<Eigen::Map<Eigen::MatrixXd> >* pPeriodStatistics();

	const Eigen::VectorXd* pMeanStatisticsMinusTargets() const;
	Eigen::VectorXd* pMeanStatisticsMinusTargets();

	const std::vector<Eigen::Map<Eigen::MatrixXd> >* pPeriodScores() const;
	std::vector<Eigen::Map<Eigen::MatrixXd> >* pPeriodScores();

	const std::vector<Eigen::MatrixXd>* pDerivatives() const;
	std::vector<Eigen::MatrixXd>* pDerivatives();

	const std::vector<Eigen::Map<Eigen::VectorXd> >* pTimes() const;
	std::vector<Eigen::Map<Eigen::VectorXd> >* pTimes();


protected:
	//! Model parameters used.
	//!
	Eigen::VectorXd* lpTheta;
	//! Targets vector, that is statistics calculated on the observations.
	//!
	const Eigen::VectorXd* lpTargets;
	//! Targets vector per period.
	//!
	const Eigen::MatrixXd* lpPeriodTargets;
	//! Mean of the main statistics minus the target statistics.
	//!
	std::vector<Eigen::VectorXd>* lpStatistics;
	//! Main statistics per period.
	//!
	std::vector<Eigen::Map<Eigen::MatrixXd> >* lpPeriodStatistics;
	//! Mean of the main statistics minus the target statistics.
	//!
	Eigen::VectorXd* lpMeanStatisticsMinusTargets;
	//! RNG state at the begin of each period.
	//!
	std::vector<std::vector<int> >* lpSeeds;
	//! Scores per period, returned only from the EpochSimulation.
	//!
	std::vector<Eigen::Map<Eigen::MatrixXd> >* lpPeriodScores;
	//! Derivatives, returned only from the MLSimulation.
	//!
	std::vector<Eigen::MatrixXd>* lpDerivatives;
	//! Times of conditional simulation.
	//!
	std::vector<Eigen::Map<Eigen::VectorXd> >* lpTimes;

};

} // namespace siena

#endif // RSIENA_RESULT_H_
