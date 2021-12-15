/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file StatisticsSimulation.h
 * \brief Defines the StatisticsSimulation class.
 *****************************************************************************/

#ifndef RSIENA_STATISTICS_SIMULATION_H_
#define RSIENA_STATISTICS_SIMULATION_H_

#include "sim/Simulation.h"

#include "logger/Logger.h"

#include "data/LongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/variables/NetworkVariable.h"

#include "utils/debugging/StopWatch.h"

namespace siena {

/**
 * Wrapper around the EpochSimulation, providing parallelization using OpenMP
 * threading.
 */
class StatisticsSimulation: public Simulation {
public:
	StatisticsSimulation(Model* pModel, Data* pData,
			const unsigned int nThreads = 1);
	virtual ~StatisticsSimulation();

	bool conditional();

	virtual int nSimulations();

	virtual void simulate();
	virtual void simulate(bool withSeeds);

	virtual void updateParameters(const Eigen::VectorXd& rTheta);

	const Eigen::VectorXd& rTargets();
	const Eigen::MatrixXd& rPeriodWiseTargets(); // !!!

protected:
	//! Number of threads use by OpenMP.
	//!
	const unsigned int lNThreads;
	//! Vector of EpochSimulations (one for each thread).
	//!
	std::vector<EpochSimulation*> lEpochSimulations;
	//! Vector of timer use to take the simulation for each thread individually.
	//!
	std::vector<StopWatch> lSimulationTime;
	//! Seeds for each period.
	//!
	// TODO threads should have there own seeds, that is rewrite
	// src/utils/Random.cpp to be OpenMP aware, with all the modifications to
	// the simulation code needed
	std::vector<std::vector<int> > lSeeds;

	//! Target statistics.
	//!
	const Eigen::VectorXd lTargets;
	//! Vector of target statistics per period.
	//!
	const Eigen::MatrixXd lPeriodWiseTargets;
	//! Mean statistics over periods.
	//!
	std::vector<Eigen::VectorXd> lStatistic;
	//! double vector for fixed layout matrices.
	//!
	std::vector<double> lStatisticsData;
	//! Vector of simulated statistics.
	//!
	std::vector<Eigen::Map<Eigen::MatrixXd> > lStatistics;
	//! double vector for fixed layout matrices.
	//!
	std::vector<double> lScoresData;
	//! Vector of simulated scores.
	//!
	std::vector<Eigen::Map<Eigen::MatrixXd> > lScores;

	//! double vector for fixed layout matrices.
	//!
	std::vector<double> lTimesData;
	//! Vector of times needed for conditional simulation.
	//!
	std::vector<Eigen::Map<Eigen::VectorXd> > lTimes;
	//! Result vector holding the mean statistics over all periods and
	//! simulations. Common name: fra
	Eigen::VectorXd lMeanStatisticsMinusTargets;
	//! Result access.
	//!
	Result lResult;

	virtual void needsChanged();

	/** (see implementation below) */
	template<typename Mapped>
	static void setupMaps(std::vector<Eigen::Map<Mapped> >& rMaps,
			std::vector<double>& rData, const unsigned int n, const int rows,
			const int cols);

	void resizeEpochSimulations(const unsigned int num);
	void updateEpochSimulations();

	void simulatePeriods(bool withSeeds);
	void simulatePeriod(const int thread, const int m);

	template<typename Derived>
	void getSingleScores(Eigen::DenseBase<Derived>& rScores, const int period,
			const EpochSimulation& rEpochSim);

private:
	// Don't copy it
	StatisticsSimulation(const StatisticsSimulation&);
	StatisticsSimulation& operator=(const StatisticsSimulation&);

	double getRateScore(const EffectInfo* effectInfo,
			const EpochSimulation& rEpochSim,
			const DependentVariable* pVariable);
	double getStructuralRateScore(const EffectInfo* effectInfo,
			const DependentVariable* pVariable,
			const NetworkVariable* pInteractionVariable);
	double getDiffusionRateScore(const EffectInfo* effectInfo,
			const EpochSimulation& rEpochSim);
	double getCovariateRateScore(const EffectInfo* effectInfo,
			const EpochSimulation& rEpochSim,
			const DependentVariable* pVariable);

};

/**
 * Helper function to set up a bunch of matrices with a fixed memory layout.
 *
 * Since we use MPI we want to have a fixed memory layout to do scattering /
 * gathering. But we don't know the parameter size at compile time, therefore
 * MatrixXd allocates its data dynamically (that is somewhere random). Here we
 * create Maps on a given memory layout.
 *
 * @param[out] rMaps Reference to a vector of maps of Eigen objects.
 * @param[out] rData Reference to a vector of doubles. This is the underlying
 *                   linear memory layout.
 * @param n Numer of Eigen objects to be hold in the rMaps vector.
 * @param rows Number of rows.
 * @param cols Number of columns.
 */
template<typename Mapped>
void StatisticsSimulation::setupMaps(std::vector<Eigen::Map<Mapped> >& rMaps,
		std::vector<double>& rData, const unsigned int n, const int rows,
		const int cols) {
	if (rMaps.size() == n) {
		// Same size, were already done
		return;
	}
	if (n == 0) {
		// Just clearing vectors
		rMaps.clear();
		rData.clear();
		return;
	}
	// Resize vectors
	const int size = rows * cols;
	rData.resize(n * size);
	// Size 1x1 allows dynamic matrix and vector types
	rMaps.resize(n, Eigen::Map<Mapped>(0/*ptr*/, 1, 1));
	LOGS(logger::Priority::DEBUG)<<"setup mapped vector "<<n<<"x"<<rows<<"x"<<cols;
	// Recreate all maps, since resizing might have relocated the memory we have
	// to update all existing maps.
	for (unsigned int i = 0; i < n; ++i) {
		// Would create a new Map (with the correct address) and then copy
		// it into the vector (new MatrixXd with new memory).
//		rMaps[i] = Map<MatrixXd>(&rData[i * size], rows, cols);

		// We want this. Relocating the Map (created with a 0/*ptr*/) to a new
		// data location. This does not allocate memory since we supply the position.
		// See http://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
		new (&rMaps[i]) Eigen::Map<Mapped>(&rData[i * size], rows, cols);
		assert(&rData[i * size] == rMaps[i].data());

		// Would be the way to go but is c++11
//		rMaps.emplace().
	}
	LOGS(logger::Priority::DEBUG);
}

/**
 * Fills the vector with the scores for the given period.
 *
 * In the first number of periods coefficients only the one for the current
 * period is non-zero.
 *
 * Similar to: siena07internals.cpp getStatistics() (the sore part)
 *
 * @param[out] rScores
 * @param period
 * @param rEpochSim
 */
template<typename Derived>
void StatisticsSimulation::getSingleScores(Eigen::DenseBase<Derived>& rScores,
		const int period, const EpochSimulation& rEpochSim) {
	assert(period == rEpochSim.period());

	int i = -1; // Score index

	// For each dependent variable
	const std::vector<LongitudinalData*>& vars = lpData->rDependentVariableData();
	for (std::vector<LongitudinalData*>::size_type v = 0; v < vars.size(); ++v) {
		const DependentVariable* pVariable = rEpochSim.pVariable(
				vars[v]->name());

		// Basic rate effects
		if (!lpModel->conditional()
				|| vars[v]->name() != lpModel->conditionalDependentVariable()) {
			for (int m = 0; m < vars[v]->observationCount() - 1; ++m) {
				if (m == period) {
					rScores[++i] = pVariable->basicRateScore();
				} else {
					rScores[++i] = 0;
				}
			}
		}

		// For each rate effect
		siena::Effects& rRate = lpModel->rRateEffects(vars[v]->name());
		for (EffectIter it = rRate.begin(); it != rRate.end(); ++it) {
			rScores[++i] = getRateScore(*it, rEpochSim, pVariable);
		}

		// For each evaluation effect
		siena::Effects& rEval = lpModel->rEvaluationEffects(vars[v]->name());
		for (EffectIter it = rEval.begin(); it != rEval.end(); ++it) {
			rScores[++i] = rEpochSim.score(*it);
		}

		// For each endownment effect
		siena::Effects& rEndowment = lpModel->rEndowmentEffects(vars[v]->name());
		for (EffectIter it = rEndowment.begin(); it != rEndowment.end(); ++it) {
			rScores[++i] = rEpochSim.score(*it);
		}

		// For each creation effect
		siena::Effects& rCreation = lpModel->rCreationEffects(vars[v]->name());
		for (EffectIter it = rCreation.begin(); it != rCreation.end(); ++it) {
			rScores[++i] = rEpochSim.score(*it);
		}
	} // for vars[v]
	assert(i + 1 == rScores.size());
}

} // namespace siena

#endif // RSIENA_STATISTICS_SIMULATION_H_
