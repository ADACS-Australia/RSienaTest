/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file FiniteDifference.h
 *
 * \class siena::FiniteDifference
 * \brief Finite differences differentiation method.
 *
 * Performs additional simulations with slightly modified parameters to
 * approximate the derivative. The parameter variation can be controlled via the
 * epsilon vector.
 *****************************************************************************/

#ifndef RSIENA_FINITE_DIFFERENCE_H_
#define RSIENA_FINITE_DIFFERENCE_H_

#include "estimator/derivative/Differentiation.h"
#include "sim/StatisticsSimulation.h"

#include <vector>

namespace siena {

class FiniteDifference: public Differentiation {
public:
	// Magic numbers
	static const double INITIAL_EPSILON;
	static const int SIMULATIONS_BEFORE_CHECK;
	static const int DERIVATIVE_CHECK_THRESHOLD;

	explicit FiniteDifference(Simulation& rSim);

	// SimulationListener
	std::vector<ResultType>& needs();
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations);

	// Differentiation
	int simulationsBeforeCheck(const int n) const;
	bool checkDerivative() const;
	const Eigen::MatrixXd& rDerivative() const;

	const Eigen::MatrixXi& rPositiveCount() const;
	void epsilon(const Eigen::VectorXd& rEpsilon);
	const Eigen::VectorXd& rEpsilon() const;

protected:

	void onEpsilonResults(const Result& rResults);

	//! The simulation.
	//!
	Simulation& lrSimuatlion;
	//! Modified model parameter vector.
	//!
	Eigen::VectorXd lTheta;
	//! Seeds of the 'base' simulation.
	//!
	std::vector<std::vector<int> > lSeeds;
	//! Near 0 vector. This is the step width per parameter.
	//!
	Eigen::VectorXd lEpsilon;
	//! Index of current parameter.
	//!
	int lParameter;
	//! 'Base' statistics.
	//!
	std::vector<Eigen::VectorXd> lStatistics;

	//! Summing up the derivatives for each simulation.
	//! Common name: sdf
	Eigen::MatrixXd lSumOfDerivatives;
	//! The derivative.
	//!
	Eigen::MatrixXd lDerivative;
	//! Number of simulations a derivative element was positive.
	//!
	Eigen::MatrixXi lPositiveCount;

private:
	// Don't copy it
	FiniteDifference(const FiniteDifference&);
	FiniteDifference& operator=(const FiniteDifference&);

	/**
	 * Helper class used to request private simulations.
	 */
	class Forward: public ResultListener {
	public:
		Forward(FiniteDifference* pFwd) :
				lpFwd(pFwd) {
		}
		std::vector<ResultType>& needs() {
			return lpFwd->needs();
		}
		void onResults(const Result& rResults) {
			lpFwd->onEpsilonResults(rResults);
		}
	private:
		FiniteDifference* const lpFwd;
	};
	//! Instance used to forward simulations.
	//!
	Forward lPrivateSim;
};

} // namespace siena

#endif // RSIENA_FINITE_DIFFERENCE_H_
