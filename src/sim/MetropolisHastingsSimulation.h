/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MetropolisHastingsSimulation.h
 *
 * \class siena::MetropolisHastingsSimulation
 * \brief Simulation wrapper for the MLSimulation used in maximum likelihood
 *        estimation.
 *****************************************************************************/

#ifndef RSIENA_METROPOLIS_HASTINGS_SIMULATION_H_
#define RSIENA_METROPOLIS_HASTINGS_SIMULATION_H_

#include "sim/StatisticsSimulation.h"

#include "logger/Logger.h"

#include "data/LongitudinalData.h"
#include "model/ml/MLSimulation.h"

#include <vector>
#include <algorithm>

namespace siena {

class MetropolisHastingsSimulation: public Simulation {
public:
	// Magic numbers
	static const double N_RUN_MH_MULTIPLICATOR_DEFAULT;

	MetropolisHastingsSimulation(Model* pModel, Data* pData);
	virtual ~MetropolisHastingsSimulation();

	virtual int nSimulations();

	virtual void simulate();
	virtual void simulate(bool);

	virtual void updateParameters(const Eigen::VectorXd& rTheta);

	const Eigen::VectorXd& rTargets();
	const Eigen::MatrixXd& rPeriodWiseTargets();

protected:
	//! Sum of the scores returned by the MLSimulation.
	//!
	std::vector<Eigen::VectorXd> lScores;
	//! Score targets, that is the 0 vector.
	//!
	const Eigen::VectorXd lScoreTargets;
	//! True if the derivative are requested.
	//!
	bool lNeedsDerivative;
	//! Sum of the derivative return from the MLSimulation.
	//!
	std::vector<Eigen::MatrixXd> lDerivative;
	//! Result vector holding the mean scores over all periods and simulations.
	//! Common name: fra
	Eigen::VectorXd lMeanScoresMinusTargets;
	//! Number of steps to run per period.
	//!
	Eigen::VectorXi lNRunMH;
	//! Result access.
	//!
	Result lResult;

	void simulatePeriod(const int m);
	virtual void needsChanged();

	template<typename Derived>
	void addSingleScores(Eigen::MatrixBase<Derived>& rScores, const int period,
			const MLSimulation& rSim);
	template<typename Derived>
	void addSingleDerivatives(Eigen::MatrixBase<Derived>& rDeriv,
			const int period, const MLSimulation& rSim);

private:
	MetropolisHastingsSimulation(const MetropolisHastingsSimulation&);
	MetropolisHastingsSimulation& operator=(
			const MetropolisHastingsSimulation&);

};

/**
 * Adds the scores of a single period to the score vector.
 *
 * Similar to: siena07internals.cpp getScores() (the scores part)
 * Only difference is the sign change which happens in maxlikec.r maxlikec() (line ~128)
 *
 * @param[out] rScores Score vector. The scores of the current period are added.
 * @param period Current period.
 * @param rSim Reference to the simulation.
 */
template<typename Derived>
void MetropolisHastingsSimulation::addSingleScores(
		Eigen::MatrixBase<Derived>& rScores, const int period,
		const MLSimulation& rSim) {
	int i = 0; // Parameter index

	// For each dependent variable
	const std::vector<LongitudinalData*>& vars = lpData->rDependentVariableData();
	for (std::vector<LongitudinalData*>::size_type v = 0; v < vars.size(); ++v) {
		// Add basic rate effects on the diagonal
		rScores[i + period] -=
				rSim.pVariable(vars[v]->name())->basicRateScore();
		i += vars[v]->observationCount() - 1;
		// Other effects
		const int num = nNonBasicRateEffects(vars[v]->name());
		for (int e = 0; e < num; ++e) {
			rScores[i] -= rSim.score(lSimulationEffects[i]);
			++i;
		}
	} // for vars[v]
}

/**
 * Adds the derivative of a single period to the derivative matrix.
 *
 * Similar to: siena07internals.cpp getScores() (the derivative part);
 *             maxlikec.r reformatDerivs()
 *
 * @param[out] rDeriv Derivative matrix. The derivatives of the current period
 *                    are added.
 * @param period Current period.
 * @param rSim Reference to the simulation.
 */
template<typename Derived>
void MetropolisHastingsSimulation::addSingleDerivatives(
		Eigen::MatrixBase<Derived>& rDeriv, const int period,
		const MLSimulation& rSim) {
	int row = 0;

	// For each dependent variable
	const std::vector<LongitudinalData*>& vars = lpData->rDependentVariableData();
	for (std::vector<LongitudinalData*>::size_type v = 0; v < vars.size(); ++v) {
		// Add basic rate effects on the diagonal
		rDeriv(row + period, row + period) -=
				rSim.pVariable(vars[v]->name())->basicRateDerivative();
		row += vars[v]->observationCount() - 1;

		// The layout in the end looks like this (where '.' is 0 / don't care)
		// X . . .  // basic rates
		// . X . .  // basic rates
		// . . X X  // symmetric block with other effects
		// . . x X  // where the entries below the diagonal are not returned by
		//          // the simulation
		const int num = nNonBasicRateEffects(vars[v]->name());
		for (int i = row; i < row + num; ++i) {
			rDeriv(i, i) -= rSim.derivative(lSimulationEffects[i],
					lSimulationEffects[i]);
			for (int j = 1 + i - row; j < num; ++j) {
				double d = rSim.derivative(lSimulationEffects[i],
						lSimulationEffects[i + j]);
				rDeriv(i, i + j) -= d;
				rDeriv(i + j, i) -= d;
			}
		}
		row += num;
	} // for vars[v]
}

} // namespace siena

#endif // RSIENA_METROPOLIS_HASTINGS_SIMULATION_H_
