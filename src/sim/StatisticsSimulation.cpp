/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file StatisticsSimulation.cpp
 * \brief Implements the StatisticsSimulation class.
 *****************************************************************************/

#include "StatisticsSimulation.h"

#include "sim/listener/Result.h"
#include "RUtil.h"
#include "Eigen/Util.h"
#include "logger/Logger.h"

#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"

#include <Rinternals.h>
#include <R_ext/Random.h>
#include "utils/Random.h"

#include "sim/modificator/MeanStatisticsCalculator.h"

#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// #include <cassert>

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

/**
 * Constructs the StatisticsSimulation object.
 *
 * @param pModel Pointer to the Model.
 * @param pData Pointer to the Data.
 * @param nThreads Number of threads used to perform parallel simulations.
 */
StatisticsSimulation::StatisticsSimulation(Model* pModel, Data* pData,
		const unsigned int nThreads) :
		Simulation(pModel, pData), //
		lNThreads(nThreads), //
		lEpochSimulations(), //
		lSeeds(pData->observationCount() - 1, vector<int>()), //
		lTargets(calculateTargetStatistics().colwise().sum()), //
		lPeriodWiseTargets(calculateTargetStatistics()), //
		lStatistic(), //
		lStatisticsData(), //
		lStatistics(), //
		lScoresData(), //
		lScores(), //
		lTimesData(), //
		lTimes(), //
		lMeanStatisticsMinusTargets(nStatistics()), //
		lResult(&lTheta, &lTargets, &lPeriodWiseTargets, &lStatistic, &lStatistics,
				&lMeanStatisticsMinusTargets, &lSeeds, &lScores, 0/*ptr*/, &lTimes) {
	assert(lNThreads > 0);
	resizeEpochSimulations(lNThreads);
	// Setup the mapped matrices for statistics and time. Statistics is always
	// needed.
	LOG(Priority::DEBUG, "setup mapped statistics vector");
	setupMaps(lStatistics, lStatisticsData, lNThreads, nPeriods(),
			nStatistics());
	// Conditional is not changed during phases, so no later change required.
	if (lpModel->conditional()) {
		LOG(Priority::DEBUG, "setup mapped times vector");
		setupMaps(lTimes, lTimesData, lNThreads, nPeriods(), 1);
	}
	// Score might not be needed if they are needsChanged, creates the stores.
	int n = lpModel->needScores() ? lNThreads : 0;
	setupMaps(lScores, lScoresData, n, nPeriods(), nParameters());

	static MeanStatisticsCalculator s;
	addResultModificator(&s);
}

StatisticsSimulation::~StatisticsSimulation() {
	// setupMaps(lStatistics, lStatisticsData, 0, 0, 0);
	// setupMaps(lTimes, lTimesData, 0, 0, 0);
	// setupMaps(lScores, lScoresData, 0, 0, 0);
	lStatistics.clear();
	lTimes.clear();
	lScores.clear();
	resizeEpochSimulations(0);
}

/**
 * \copydoc Simulation::needsChanged()
 */
void StatisticsSimulation::needsChanged() {
	set<ResultType>::iterator it = find(lNeeds.begin(), lNeeds.end(),
			PERIOD_SCORES);
	lpModel->needScores(it != lNeeds.end());
	if (lpModel->needScores()) {
		LOGS(Priority::DEBUG)<<"simulation collects scores";
		// Create stores only the first time.
		if (lScores.empty()) {
			setupMaps(lScores, lScoresData, lNThreads, nPeriods(), nParameters());
		}
	}
}

/**
 * @return True if the simulation is conditional.
 */
bool StatisticsSimulation::conditional() {
	return lpModel->conditional();
}

/**
 * Handles the creation and deletion of EpochSimulation objects.
 *
 * @param num Desired number of EpochSimulation objects.
 */
void StatisticsSimulation::resizeEpochSimulations(const unsigned int num) {
	// Delete superfluous simulations
	if (lEpochSimulations.size() > num) {
		for (vector<EpochSimulation*>::iterator it = lEpochSimulations.begin()
				+ num; it != lEpochSimulations.end(); ++it) {
			delete *it;
		}
		lEpochSimulations.resize(num);
	}
	// Clone new simulations in parallel
#pragma omp parallel num_threads(lNThreads)
	{
#pragma omp for schedule(dynamic, 1)
		for (vector<EpochSimulation*>::size_type i = lEpochSimulations.size();
				i < num; ++i) {
			EpochSimulation* ptr = new EpochSimulation(lpData, lpModel);
#pragma omp critical
			{
				lEpochSimulations.push_back(ptr);
			}
		}
	}
	// Debugging
	lSimulationTime.clear();
	for (unsigned int i = 0; i < num; ++i) {
		lSimulationTime.push_back(StopWatch("simulation time"));
	}
}

/**
 * \copydoc Simulation::nSimulations()
 */
int StatisticsSimulation::nSimulations() {
	return lNThreads;
}

/**
 * \copydoc Simulation::simulate()
 */
void StatisticsSimulation::simulate() {
	return simulate(false);
}

/**
 * \copydoc Simulation::simulate(bool)
 */
void StatisticsSimulation::simulate(bool withSeeds) {
	simulatePeriods(withSeeds);
	LOGS(Priority::DEBUG)<<"simulation finished, distribute results";
	fireResultModificators(lResult);
	fireResult(lResult);
}

/**
 * Sequentially simulate each period.
 *
 * @param withSeeds True if stored seeds should be used to reset the RNG at the
 *        beginning of each period.
 */
void StatisticsSimulation::simulatePeriods(bool withSeeds) {
	GetRNGstate();
	#pragma omp parallel num_threads(lNThreads)
	{
		#pragma omp for schedule(dynamic, 1)
		for (unsigned int thread = 0; thread < lNThreads; ++thread) {
			// Run epoch simulation for each period
			LOGS(Priority::DEBUG) << "simulate with theta: " << rParameters().transpose();
			for (int m = 0; m < lpData->observationCount() - 1; ++m) {
				if (withSeeds) {
					setSeed(lSeeds[m]);
				} else {
					getSeed(lSeeds[m]);
				}
				simulatePeriod(thread, m);
			}
		}
	}
	PutRNGstate();
}

/**
 * Runs the epoch simulation for the period and updates the variables.
 *
 * @param thread Thread number.
 * @param m Period.
 */
void StatisticsSimulation::simulatePeriod(const int thread, const int m) {
	printRNGTrace
	lSimulationTime[thread].start();

	lEpochSimulations[thread]->runEpoch(m);

	State state(lEpochSimulations[thread]);
	StatisticCalculator calculator(lpData, lpModel, &state, m);

	lSimulationTime[thread].stop();

	Map<MatrixXd>::RowXpr statsRow = lStatistics[thread].row(m);
	getSingleStatistics(statsRow, m, calculator);

	if (lpModel->needScores()) {
		Map<MatrixXd>::RowXpr scoreRow = lScores[thread].row(m);
		getSingleScores(scoreRow, m, *lEpochSimulations[thread]);
	}

	if (lpModel->conditional()) {
		lTimes[thread][m] = lEpochSimulations[thread]->time();
	}
}

/**
 * \copydoc Simulation::updateParameters()
 */
void StatisticsSimulation::updateParameters(const Eigen::VectorXd& rTheta) {
// Update the model
	updateLocalParameters(rTheta);
// Propagate changes to the simulations
	updateEpochSimulations();
}

/**
 * Propagate the model parameters to all EpochSimulations.
 */
void StatisticsSimulation::updateEpochSimulations() {
	// Update each simulation
#pragma omp parallel num_threads(lNThreads)
	{
#pragma omp for schedule(dynamic, 1)
		for (unsigned int thread = 0; thread < lNThreads; ++thread) {
			// This is basically EpochSimulation::updateParameters(int) for any
			// period (and therefore saves some calls to
			// DependentVariable::updateEffectParameters()). We do need this since
			// we do not recreate EpochSimulations but store them.
			const vector<DependentVariable*>& vars =
					lEpochSimulations[thread]->rVariables();
			for (vector<DependentVariable*>::size_type i = 0; i < vars.size();
					++i) {
				for (int m = 0; m < lpData->observationCount() - 1; ++m) {
					vars[i]->updateBasicRate(m);
				}
				vars[i]->updateEffectParameters();
			}
		}
	}
}

/**
 * Similar to: siena07setup.cpp getTargets()
 *
 * Common name: @f$ s @f$, targets
 *
 * @return Vector of observed statistics.
 */
const VectorXd& StatisticsSimulation::rTargets() {
	return lTargets;
}

/**
 * Common name: targets2
 *
 * @return Vector of periodwise observed statistics.
 */
const MatrixXd& StatisticsSimulation::rPeriodWiseTargets() { // !!!
	return lPeriodWiseTargets;
}

/**
 * Helper function to emulation something like EpochSimulation::score() for
 * rate effects.
 *
 * Mostly copied from: siena07internals.cpp getStatistics()
 */
double StatisticsSimulation::getRateScore(const EffectInfo* effectInfo,
		const EpochSimulation& rEpochSim, const DependentVariable* pVariable) {
	if (effectInfo->rateType() == "structural") {
		if (effectInfo->interactionName1() == "") {
			return getStructuralRateScore(effectInfo, pVariable,
					dynamic_cast<const NetworkVariable*>(pVariable));
		} else {
			return getStructuralRateScore(effectInfo, pVariable,
					dynamic_cast<const NetworkVariable*>(rEpochSim.pVariable(
							effectInfo->interactionName1())));
		}
	} else if (effectInfo->rateType() == "diffusion") {
		return getDiffusionRateScore(effectInfo, rEpochSim);
	} else {
		return getCovariateRateScore(effectInfo, rEpochSim, pVariable);
	}
}

/**
 * Helper function to emulation something like EpochSimulation::score() for
 * rate effects.
 *
 * Mostly copied from: siena07internals.cpp getStatistics()
 */
double StatisticsSimulation::getStructuralRateScore(
		const EffectInfo* effectInfo, const DependentVariable* pVariable,
		const NetworkVariable* pInteractionVariable) {
	if (effectInfo->effectName() == "outRate") {
		return pVariable->outDegreeScore(pInteractionVariable);
	}
	if (effectInfo->effectName() == "inRate") {
		return pVariable->inDegreeScore(pInteractionVariable);
	}
	if (effectInfo->effectName() == "recipRate") {
		return pVariable->reciprocalDegreeScore(pInteractionVariable);
	}
	if (effectInfo->effectName() == "outRateInv") {
		return pVariable->inverseOutDegreeScore(pInteractionVariable);
	}
	LOGF(Priority::ERROR, "Unexpected rate effect '%s'.",
			effectInfo->effectName().c_str());
	return 0;
}

/**
 * Helper function to emulation something like EpochSimulation::score() for
 * rate effects.
 *
 * Mostly copied from: siena07internals.cpp getStatistics()
 */
double StatisticsSimulation::getDiffusionRateScore(const EffectInfo* effectInfo,
		const EpochSimulation& rEpochSim) {
	if (effectInfo->effectName() == "avExposure"
			|| effectInfo->effectName() == "totExposure"
			|| effectInfo->effectName() == "susceptAvIn"
			|| effectInfo->effectName() == "infectDeg"
			|| effectInfo->effectName() == "infectIn"
			|| effectInfo->effectName() == "infectOut"
			|| effectInfo->effectName() == "susceptAvCovar"
			|| effectInfo->effectName() == "infectCovar") {
		return rEpochSim.score(effectInfo);
	}
	LOGF(Priority::ERROR, "Unexpected rate effect '%s'.",
			effectInfo->effectName().c_str());
	return 0;
}

/**
 * Helper function to emulation something like EpochSimulation::score() for
 * rate effects.
 *
 * Mostly copied from: siena07internals.cpp getStatistics()
 */
double StatisticsSimulation::getCovariateRateScore(const EffectInfo* effectInfo,
		const EpochSimulation& rEpochSim, const DependentVariable* pVariable) {
	if (ConstantCovariate* pConstantCov = lpData->pConstantCovariate(
			effectInfo->interactionName1())) {
		return pVariable->constantCovariateScore(pConstantCov);
	}
	if (ChangingCovariate* pChangingCov = lpData->pChangingCovariate(
			effectInfo->interactionName1())) {
		return pVariable->changingCovariateScore(pChangingCov);
	}
	if (const BehaviorVariable* pBehavior =
			dynamic_cast<const BehaviorVariable*>(rEpochSim.pVariable(
					effectInfo->interactionName1()))) {
		return pVariable->behaviorVariableScore(pBehavior);
	}
	LOGF(Priority::ERROR, "No individual covariate '%s'.",
			effectInfo->interactionName1().c_str());
	return 0;
}

} // namespace siena
