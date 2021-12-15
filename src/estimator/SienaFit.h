/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file SienaFit.h
 * \brief Defines the SienaFit class.
 *****************************************************************************/

#ifndef RSIENA_SIENA_FIT_H_
#define RSIENA_SIENA_FIT_H_

#include "Eigen/Types.h"
#include <Eigen/Dense>
#include <vector>

namespace siena {

/**
 * C++ equivalent of the sienafit object in R.
 *
 * Contains the result of the estimation and some other performance
 * information like the number of iterations.
 */
class SienaFit {
public:
	SienaFit();
	virtual ~SienaFit();

	const std::vector<int>& rPhase1NSimulations() const;
	void phase1NSimulations(const std::vector<int>& rPhase1NSimulations);

	const std::vector<int>& rPhase2NSimulations() const;
	void phase2NSimulations(const std::vector<int>& rPhase2NSimulations);

	int phase3NSimulations() const;
	void phase3NSimulations(const int phase3NSimulations);

	const Eigen::VectorXb& rFixed() const;
	void fixed(const Eigen::VectorXb& rFixed);

	const Eigen::VectorXd& rTargets() const;
	void targets(const Eigen::VectorXd& rTargets);
	
	const Eigen::MatrixXd& rPeriodWiseTargets() const;
	void periodWiseTargets(const Eigen::MatrixXd& rPeriodWiseTargets);
	
	const Eigen::VectorXd& rTheta() const;
	void theta(const Eigen::VectorXd& rTheta);

	const Eigen::VectorXd& rTRatio() const;
	void tRatio(const Eigen::VectorXd& rTRatio);

	const Eigen::VectorXd& rRate() const;
	void rate(const Eigen::VectorXd& rRate);

	const Eigen::VectorXd& rRateError() const;
	void rateError(const Eigen::VectorXd& rRateError);

	const Eigen::MatrixXd& rBlockStructure() const;
	void blockStructure(const Eigen::MatrixXd& rBlockStructure);

	const Eigen::MatrixXd& rGamma() const;
	void gamma(const Eigen::MatrixXd& rGamma);

	const Eigen::MatrixXd& rPhase1Weight() const;
	void phase1Weight(const Eigen::MatrixXd& rWeight);
	const Eigen::MatrixXd& rPhase3Weight() const;
	void phase3Weight(const Eigen::MatrixXd& rWeight);

	const Eigen::MatrixXd& rPhase1Derivative() const;
	void phase1Derivative(const Eigen::MatrixXd& rDerivative);
	const Eigen::MatrixXd& rPhase3Derivative() const;
	void phase3Derivative(const Eigen::MatrixXd& rDerivative);

	const Eigen::MatrixXd& rStatistics() const;
	void statistics(const Eigen::MatrixXd& rStatistics);
	const std::vector<Eigen::MatrixXd>& rPeriodWiseScores() const;
	void periodWiseScores(const std::vector<Eigen::MatrixXd>& rPeriodWiseScores);
	const std::vector<Eigen::MatrixXd>& rPeriodWiseStatistics() const;
	void periodWiseStatistics(const std::vector<Eigen::MatrixXd>& rPeriodWiseStatistics);

	const Eigen::MatrixXd& rPhase1Covariance() const;
	void phase1Covariance(const Eigen::MatrixXd& rCovariance);
	const Eigen::MatrixXd& rPhase3Covariance() const;
	void phase3Covariance(const Eigen::MatrixXd& rCovariance);

private:
	//! Simulations done in phase 1.
	//!
	std::vector<int> lPhase1NSimulations;
	//! Simulations done in each subphase of phase 2. The length of the vector
	//! is the number of subphase performed.
	std::vector<int> lPhase2NSimulations;
	//! Simulations done in phase 3.
	//!
	int lPhase3NSimulations;

	//! Vector of parameters fixed because their phase3Derivative could not be
	//! approximated.
	Eigen::VectorXb lFixed;
	//! Observed statistics.
	//!
	Eigen::VectorXd lTargets;
	//! Estimated parameter values.
	//!
	Eigen::VectorXd lTheta;
	//! T-ratio of the estimated parameter values.
	//!
	Eigen::VectorXd lTRatio;
	//! Estimated rate parameters when using conditional simulations.
	//!
	Eigen::VectorXd lRate;
	//! Error of estimated rate parameters when using conditional simulations.
	//!
	Eigen::VectorXd lRateError;
	//! Block structure mask for the covariance, 0 if the corresponding entry is
	//! neglected, otherwise 1 (apply = element wise multiplication).
	Eigen::MatrixXd lBlockStructure;
	//! Original phase3Derivative from phase 3.
	//!
	Eigen::MatrixXd lPhase3Gamma;
	//! GMM weights from phase 1.
	//!
	Eigen::MatrixXd lPhase1Weight;
	//! GMM weights from phase 3.
	//!
	Eigen::MatrixXd lPhase3Weight;
	//! Derivative from phase 1.
	//!
	Eigen::MatrixXd lPhase1Derivative;
	//! Derivative from phase 3.
	//!
	Eigen::MatrixXd lPhase3Derivative;
	//! All statistics calculated in phase 3. Needed when restarting with
	//! previous answer.
	Eigen::MatrixXd lPhase3Statistics;
	//! Eigen has a tensor module. But we just need this exported to R.
	Eigen::MatrixXd lPeriodWiseTargets;
	std::vector<Eigen::MatrixXd> lPeriodWiseStatistics;
	std::vector<Eigen::MatrixXd> lPeriodWiseScores;
	//! Covariance matrix of statistics from phase 1.
	//!
	Eigen::MatrixXd lPhase1Covariance;
	//! Covariance matrix of statistics from phase 3.
	//!
	Eigen::MatrixXd lPhase3Covariance;

};

} // namespace siena

#endif // RSIENA_SIENA_FIT_H_
