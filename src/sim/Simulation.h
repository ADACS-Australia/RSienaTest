/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Simulation.h
 *
 * \class siena::Simulation
 * \brief Base class for simulations.
 *
 * The Simulation objects provide the interface from RSienas Model/Data to
 * RSienaCPP. It therefore has a lot in common with the siena07xxx.cpp files,
 * except that it returns Eigen objects not R objects.
 *****************************************************************************/

#ifndef RSIENA_SIMULATION_H_
#define RSIENA_SIMULATION_H_

#include <Eigen/Dense>
#include <list>
#include <vector>

#include "listener/Result.h"
#include "listener/ResultType.h"

#include "data/Data.h"
#include "data/LongitudinalData.h"
#include "model/Model.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"
#include "model/StatisticCalculator.h"

namespace siena {

// Forward declarations
class ResultListener;
class ResultModificator;

//! Vector of EffectInfo pointers.
//!
typedef const std::vector<EffectInfo*> Effects;
//! Constant iterator of EffectInfo pointers.
//!
typedef std::vector<EffectInfo*>::const_iterator EffectIter;
//! Custom Eigen::Vector holding EffectInfo pointers.
//!
typedef Eigen::Matrix<EffectInfo*, Eigen::Dynamic, 1> VectorXEffectPtr;

class Simulation {
public:

	//! Dummy effect info for basic rate parameters.
	//!
	static EffectInfo RATE_EFFECT;

	Simulation(Model* pModel, Data* pData);
	virtual ~Simulation() = 0;

	void addResultModificator(ResultModificator* const pModificator);
	void removeResultModificator(ResultModificator* const pModificator);

	void addResultListener(ResultListener* const pListener);
	void removeResultListener(ResultListener* const pListener);
	std::vector<ResultListener*> clearResultListener();

	/**
	 * @return Number of simulations executed in parallel.
	 */
	virtual int nSimulations() = 0;

	/**
	 * Performs the simulations.
	 */
	virtual void simulate() = 0;

	/**
	 * Performs the simulation resetting seeds for each period to the last saved
	 * state.
	 *
	 * @param withSeeds True when saved seeds should be used.
	 */
	virtual void simulate(bool withSeeds) = 0;

	Model* pModel();
	Data* pData();
	int nPeriods();

	int nParameters();

	/**
	 * @return reference to (a copy) of the actual parameter vector. Changing
	 * this one affects nothing.
	 */
	const Eigen::VectorXd& rParameters();

	/**
	 * Updates the model parameters.
	 *
	 * @param rTheta New parameter vector.
	 */
	virtual void updateParameters(const Eigen::VectorXd& rTheta) = 0;

	const VectorXEffectPtr& rSimulationEffects();

	int nStatistics();
	const VectorXEffectPtr& rStatisticEffects();

	/**
	 * Returns the primary statistics calculated on the observations.
	 *
	 * @return Target vector.
	 */
	virtual const Eigen::VectorXd& rTargets() = 0;

	/**
	 * @return periodwise Target vector.
	 */
	virtual const Eigen::MatrixXd& rPeriodWiseTargets() = 0; // !!!

protected:
	//! The Model, that is parameters.
	//!
	Model* lpModel;
	//! The Data, that is network/behavior variables.
	//!
	Data* lpData;
	//! Result modifications.
	//!
	std::list<ResultModificator*> lResultModificatorPtrs;
	//! Result listeners.
	//!
	std::vector<ResultListener*> lResultListenerPtrs;
	//! The result set requested by the modificators and listeners. This vector
	//! is maintained automatically when adding removing listeners/modificators.
	std::set<ResultType> lNeeds;

	//! Model parameter vector.
	//!
	Eigen::VectorXd lTheta;
	//! Vector of effects used during the simulation to calculate the
	//! contributions of a step. There are as many simulation effects as
	//! parameters. Defined here to provide a single interface for the
	//! parameters (not 4 separated like in Model).
	VectorXEffectPtr lSimulationEffects;
	//! Vector of effects used to calculate the statistics.
	//!
	VectorXEffectPtr lStatisticEffects;

	/**
	 * Called when listeners/modificators were changed.
	 *
	 * Subclasses might want to change their composition because other results
	 * are required.
	 */
	virtual void needsChanged() = 0;

	void fireResultModificators(Result& rResult);
	void fireResult(const Result& rResult);

	void updateLocalParametersFromModel();
	void updateLocalParameters(const Eigen::VectorXd& rTheta);
	int nNonBasicRateEffects(const std::string& variableName);

	template<typename Derived>
	void getSingleStatistics(Eigen::MatrixBase<Derived>& rStats,
			const int period, const StatisticCalculator& rCalculator);

	/**
	 * Calculates the effect statistics on the observed data.
	 *
	 * @return Matrix with row wise period target.
	 */
	Eigen::MatrixXd calculateTargetStatistics();

private:
	Simulation(const Simulation&);
	Simulation& operator=(const Simulation&);

	void addRateEffects(std::vector<EffectInfo*>& rEffects, LongitudinalData* pVar);
	void addObjectiveEffects(std::vector<EffectInfo*>& rEffects,
			LongitudinalData* pVar);

	VectorXEffectPtr createSimulationEffectVector();
	VectorXEffectPtr createStatisticEffectVector();

};

/**
 * Fills the vector with the statistics for the given period.
 *
 * In the first number of periods coefficients only the one for the current
 * period is non-zero.
 *
 * Similar to: siena07internals.cpp getStatistics() (the statistics part)
 *
 * @param[out] rStats
 * @param period
 * @param rCalculator
 */
template<typename Derived>
void Simulation::getSingleStatistics(Eigen::MatrixBase<Derived>& rStats,
		const int period, const StatisticCalculator& rCalculator) {
	int i = -1; // Statistic index

	// For each dependent variable
	const std::vector<LongitudinalData*>& vars = lpData->rDependentVariableData();
	for (std::vector<LongitudinalData*>::size_type v = 0; v < vars.size(); ++v) {
		BehaviorLongitudinalData* pBehavior = lpData->pBehaviorData(
				vars[v]->name());

		// Basic rate effects
		if (!lpModel->conditional()
				|| vars[v]->name() != lpModel->conditionalDependentVariable()) {
			for (int m = 0; m < vars[v]->observationCount() - 1; ++m) {
				if (m == period) {
					rStats[++i] = rCalculator.distance(vars[v], m);
				} else {
					rStats[++i] = 0;
				}
			}
		}

		// For each rate effect
		Effects& rRate = lpModel->rRateEffects(vars[v]->name());
		for (EffectIter it = rRate.begin(); it != rRate.end(); ++it) {
			rStats[++i] = rCalculator.statistic(*it);
		}

		Effects& rGMM = lpModel->rGMMEffects(vars[v]->name());
		if (!rGMM.empty()) {
			// For each gmm effect
			for (EffectIter it = rGMM.begin(); it != rGMM.end(); ++it) {
				rStats[++i] = rCalculator.statistic(*it);
			}
		} else {
			// For each evaluation effect
			Effects& rEval = lpModel->rEvaluationEffects(vars[v]->name());
			for (EffectIter it = rEval.begin(); it != rEval.end(); ++it) {
				rStats[++i] = rCalculator.statistic(*it);
			}

			// For each endownment effect
			Effects& rEndowment = lpModel->rEndowmentEffects(vars[v]->name());
			for (EffectIter it = rEndowment.begin(); it != rEndowment.end();
					++it) {
				if (pBehavior) {
					rStats[++i] = -rCalculator.statistic(*it);
				} else {
					rStats[++i] = rCalculator.statistic(*it);
				}
			}

			// For each creation effect
			Effects& rCreation = lpModel->rCreationEffects(vars[v]->name());
			for (EffectIter it = rCreation.begin(); it != rCreation.end();
					++it) {
				rStats[++i] = rCalculator.statistic(*it);
			}
		}
	} // for vars[v]
	assert(i + 1 == rStats.size());
}

} // namespace siena

#endif // RSIENA_SIMULATION_H_
