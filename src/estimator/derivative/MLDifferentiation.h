/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MLDifferentiation.h
 *
 * \class siena::MLDifferentiation
 * \brief Scores and deviation differentiation method.
 *****************************************************************************/

#ifndef RSIENA_MAXIMUM_LIKELIHOOD_DIFFERENTIATION_H_
#define RSIENA_MAXIMUM_LIKELIHOOD_DIFFERENTIATION_H_

#include "estimator/derivative/Differentiation.h"

#include "sim/StatisticsSimulation.h"

namespace siena {

/**
 * Maximum likelihood differentiation method.
 *
 * This sums up the derivatives calculated by the
 * MetropolisHastingsSimulation.
 */
class MLDifferentiation: public Differentiation {
public:
	MLDifferentiation(int nParameters);

	// SimulationListener
	std::vector<ResultType>& needs();
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations);

	// MLDifferentiation
	int simulationsBeforeCheck(const int n) const;
	bool checkDerivative() const;
	const Eigen::MatrixXd& rDerivative() const;

protected:
	//! Summing up the derivatives.
	//!
	Eigen::MatrixXd lSumOfDerivative;
	//! The derivative.
	//!
	Eigen::MatrixXd lDerivative;
	void addDerivatives(const Eigen::MatrixXd& rDerivative);

private:
	// Don't copy it
	MLDifferentiation(const MLDifferentiation&);
	MLDifferentiation& operator=(const MLDifferentiation&);

};

} // namespace siena

#endif // RSIENA_MAXIMUM_LIKELIHOOD_DIFFERENTIATION_H_
