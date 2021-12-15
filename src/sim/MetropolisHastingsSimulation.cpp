/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MetropolisHastingsSimulation.cpp
 * Implementation of MetropolisHastingsSimulation.h
 *****************************************************************************/

#include "MetropolisHastingsSimulation.h"

#include "RUtil.h"
#include "Eigen/Util.h"
#include "logger/Logger.h"

#include "data/LongitudinalData.h"
#include "model/State.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/EffectInfo.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "utils/Random.h"
#include "R_ext/Random.h"

#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// #include <cassert>

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

//! Default multiplicator for the number of steps to be performed by the
//! MLSimulation.
const double MetropolisHastingsSimulation::N_RUN_MH_MULTIPLICATOR_DEFAULT = 5;

/**
 * Constructs the MetropolisHastingsSimulation object.
 *
 * @param pModel Pointer to the Model.
 * @param pData Pointer to the Data.
 */
MetropolisHastingsSimulation::MetropolisHastingsSimulation(Model* pModel,
		Data* pData) :
		Simulation(pModel, pData), //
		lScores(1, VectorXd(nStatistics())), //
		lScoreTargets(VectorXd::Zero(nStatistics())), //
		lNeedsDerivative(false), //
		lDerivative(1, MatrixXd(nParameters(), nParameters())), //
		lMeanScoresMinusTargets(), //
		lNRunMH(nPeriods()), //
		lResult(&lTheta, &lScoreTargets, 0/*ptr*/, &lScores, 0/*ptr*/,
				&lMeanScoresMinusTargets, 0/*ptr*/, 0/*ptr*/, &lDerivative, 0/*ptr*/)
	{
	// Conditional is not possible with maximum likelihood
	// R: sienaModelCreate.r sienaModelCreate() (line ~56)
	assert(!pModel->conditional());
	// Calculate the number of simulation steps per period.
	// R: initializeFran.r initializeFRAN() (line ~504)
	const MatrixXd targets = calculateTargetStatistics();
	const ArrayXb basicRates = rStatisticEffects().array()
			== &Simulation::RATE_EFFECT;
	const VectorXd zero = VectorXd::Zero(nStatistics());
	for (int m = 0; m < nPeriods(); ++m) {
		lNRunMH[m] = N_RUN_MH_MULTIPLICATOR_DEFAULT
				* basicRates.select(targets.row(m).transpose(), zero).sum();
		if (lNRunMH[m] > 100) {
			lNRunMH[m] = round((double) lNRunMH[m] / 100.0) * 100;
		}
	}
}

MetropolisHastingsSimulation::~MetropolisHastingsSimulation() {
}

/**
 * \copydoc Simulation::needsChanged()
 */
void MetropolisHastingsSimulation::needsChanged() {
	set<ResultType>::iterator it = find(lNeeds.begin(), lNeeds.end(),
			ML_DERIVATIVE);
	lNeedsDerivative = it != lNeeds.end();
	LOGS(Priority::DEBUG)<<"simulation collects derivatives = "<<lNeedsDerivative;
}

/**
 * \copydoc Simulation::nSimulations()
 */
int MetropolisHastingsSimulation::nSimulations() {
	return 1;
}

/**
 * \copydoc Simulation::simulate()
 */
void MetropolisHastingsSimulation::simulate() {
	GetRNGstate();

	lScores[0].fill(0);
	lDerivative[0].fill(0);

	LOGS(Priority::DEBUG)<<"Simulate with theta: "<<rParameters().transpose();

	for (int m = 0; m < lpData->observationCount() - 1; ++m) {
		printRNGTrace
		simulatePeriod(m);
	}

	PutRNGstate();

	lMeanScoresMinusTargets = lScores[0];
	fireResult(lResult);
}

/**
 * \copydoc Simulation::simulate(bool)
 */
void MetropolisHastingsSimulation::simulate(bool withSeeds) {
	assert(withSeeds == false); // No finite difference available
	simulate();
}

/**
 * Runs the epoch simulation for the period and updates the variables.
 *
 * Similar to: siena07models.cpp mlPeriod()
 *
 * @param m Period.
 */
void MetropolisHastingsSimulation::simulatePeriod(const int m) {
	MLSimulation* pSim = new MLSimulation(lpData, lpModel);

	// Update simulations
	pSim->simpleRates(lpModel->simpleRates());
	pSim->currentPermutationLength(lpModel->currentPermutationLength(m));

	pSim->missingNetworkProbability(
			static_cast<const Model*>(lpModel)->missingNetworkProbability(m));

	pSim->missingBehaviorProbability(
			static_cast<const Model*>(lpModel)->missingBehaviorProbability(m));

	LOGS(Priority::VERBOSE)<<"\nSimple rates: "<<pSim->simpleRates()
	<<"\nCurrent permutation length: "<<pSim->currentPermutationLength()
	<<"\nMissing network probability: "<<pSim->missingNetworkProbability()
	<<"\nMissing behavior probability: "<<pSim->missingBehaviorProbability();

	pSim->pChain(lpModel->rChainStore(m).back()->copyChain());
	lpModel->needScores(false);
	lpModel->needDerivatives(false);
	lpModel->numberMLSteps(lNRunMH[m]);

	LOGS(Priority::VERBOSE)<<"\nNum steps: "<<lpModel->numberMLSteps();

	pSim->runEpoch(m);
	// Run through current state of chain and calculate scores and derivatives.
	lpModel->needScores(true); // !onlyLoglik (bayes)
	lpModel->needDerivatives(lNeedsDerivative);

	pSim->updateProbabilities(pSim->pChain(), pSim->pChain()->pFirst()->pNext(),
			pSim->pChain()->pLast()->pPrevious());
	// Store chain on Model.
	Chain* pChain = pSim->pChain();
	pChain->createInitialStateDifferences();
	pSim->createEndStateDifferences();
	lpModel->chainStore(*pChain, m);
	lpModel->currentPermutationLength(m, pSim->currentPermutationLength());

	// Add up period results.
	addSingleScores(lScores[0], m, *pSim);
	addSingleDerivatives(lDerivative[0], m, *pSim);

	LOGS(Priority::DEBUG)<<"Scores: "<<lScores[0].transpose();
}

/**
 * \copydoc Simulation::updateParameters()
 */
void MetropolisHastingsSimulation::updateParameters(
		const Eigen::VectorXd& rTheta) {
	updateLocalParameters(rTheta);
}

/**
 * \copydoc Simulation::rTargets()
 */
const VectorXd& MetropolisHastingsSimulation::rTargets() {
	static VectorXd targets(VectorXd::Zero(nParameters()));
	return targets;
}

/**
 * \copydoc Simulation::rPeriodWiseTargets()
 */
const MatrixXd& MetropolisHastingsSimulation::rPeriodWiseTargets() {
	static MatrixXd periodWiseTargets(MatrixXd::Zero(nParameters(),nPeriods()));
	return periodWiseTargets;
}

} // namespace siena
