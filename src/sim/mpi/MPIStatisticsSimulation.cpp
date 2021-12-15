/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MPIStatisticsSimulation.cpp
 * \brief Implements the MPIStatisticsSimulation class.
 *****************************************************************************/

#ifdef MPI2

#include "MPIStatisticsSimulation.h"

#include "logger/Logger.h"
#include "RUtil.h"

#include "data/LongitudinalData.h"
#include "model/State.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "utils/Random.h"
#include "R_ext/Random.h"

// #include <cassert>

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

enum CMD {
	DONE,
	SIMULATE,
	SIMULATE_WITH_SEEDS,
	UPDATE_PARAMETERS,
	NEED_SCORES,
	DONT_NEED_SCORES
};

/**
 * Constructs the MPIStatisticsSimulation object.
 *
 * @param pModel Pointer to the Model.
 * @param pData Pointer to the Data.
 * @param nThreads Number of threads used to perform parallel simulations.
 */
MPIStatisticsSimulation::MPIStatisticsSimulation(Model* pModel, Data* pData,
		const unsigned int nThreads) :
		StatisticsSimulation(pModel, pData, nThreads), //
		// lTheta(nParameters()), //
		lMPIStatisticsData(), //
		lMPIStatistics(), //
		lMPIScoresData(), //
		lMPIScores(), //
		lMPITimesData(), //
		lMPITimes() {
	// Set the pointers of the results to the larger MPI stores.
	lResult = Result(&lTheta, &lTargets, 0/*ptr*/, &lStatistic, &lMPIStatistics,
			&lMeanStatisticsMinusTargets, &lSeeds, &lMPIScores, 0/*ptr*/, &lMPITimes);
	// Setup fixed memory layout mapps for the MPI stores.
	if (MPICommunicator::theCommunicator().isRoot()) {
		LOG(Priority::DEBUG, "setup mapped MPI statistics vector");
		setupMaps(lMPIStatistics, lMPIStatisticsData,
				lNThreads * MPICommunicator::theCommunicator().nNodes(),
				nPeriods(), nStatistics());
		if (pModel->conditional()) {
			LOG(Priority::DEBUG, "setup mapped MPI time vectors");
			setupMaps(lMPITimes, lMPITimesData,
					lNThreads * MPICommunicator::theCommunicator().nNodes(),
					nPeriods(), 1);
		}
	}
	LOG(Priority::DEBUG, "MPI simulation is ready");
}

MPIStatisticsSimulation::~MPIStatisticsSimulation() {
	lMPIStatistics.clear();
	lMPIScores.clear();
	if (MPICommunicator::theCommunicator().isRoot()) {
		int cmd = DONE;
		MPICommunicator::theCommunicator().broadcastInt(&cmd);
	}
}

/**
 * \copydoc StatisticsSimulation::needsChanged()
 */
void MPIStatisticsSimulation::needsChanged() {
	 assert(MPICommunicator::theCommunicator().isRoot());
	// Change the local store.
	StatisticsSimulation::needsChanged();
	// Inform the cluster.
	int cmd = lpModel->needScores() ? NEED_SCORES : DONT_NEED_SCORES;
	MPICommunicator::theCommunicator().broadcastInt(&cmd);
	// And prepare the vector used to gather values from the cluster.
	if (lpModel->needScores()) {
		setupMaps(lMPIScores, lMPIScoresData,
				lNThreads * MPICommunicator::theCommunicator().nNodes(),
				nPeriods(), nStatistics());
	} else {
		lMPIScores.clear();
		lMPIScoresData.clear();
	}
}

/**
 * Start listening to the commands of the MPIStatisticsSimulation on the root
 * node.
 *
 * This is the only method of the MPIStatisticsSimulation that should be used
 * on slaves nodes directly. It blocks and terminates after the root
 * simulation was destroyed.
 */
// Ok this is more like emulation of remote function calls. Some of the
// messages send here would be unnecessary if parallelization were done in the
// estimator. But I think the encapsulation achieved this way is more
// desirable.
void MPIStatisticsSimulation::slaveLoop() {
	LOG(Priority::DEBUG, "");
	assert(!MPICommunicator::theCommunicator().isRoot());
	int cmd;
	bool withSeeds;
	bool notdone = true;
	while (notdone) {
		MPICommunicator::theCommunicator().broadcastInt(&cmd);
		withSeeds = false;
		switch (cmd) {
		case DONE:
			notdone = false;
			break;
		case NEED_SCORES:
			// lpModel->needScores(true);
			lNeeds.insert(PERIOD_SCORES);
			StatisticsSimulation::needsChanged();
			break;
		case DONT_NEED_SCORES:
			// lpModel->needScores(false);
			lNeeds.erase(PERIOD_SCORES);
			StatisticsSimulation::needsChanged();
			break;
		case SIMULATE_WITH_SEEDS:
			withSeeds = true;
			// no break, I know
		case SIMULATE:
			StatisticsSimulation::simulatePeriods(withSeeds);
			MPICommunicator::theCommunicator().gatherDouble(&lStatisticsData);
			if (lpModel->needScores())
				MPICommunicator::theCommunicator().gatherDouble(&lScoresData);
			if (lpModel->conditional())
				MPICommunicator::theCommunicator().gatherDouble(&lTimesData);
			break;
		case UPDATE_PARAMETERS:
			MPICommunicator::theCommunicator().broadcastVectorXd(&lTheta);
			StatisticsSimulation::updateParameters(lTheta);
			break;
		}
	}
}

/**
 * \copydoc StatisticsSimulation::nSimulations()
 */
int MPIStatisticsSimulation::nSimulations() {
	return lNThreads * MPICommunicator::theCommunicator().nNodes();
}

/**
 * \copybrief StatisticsSimulation::simulate(bool)
 *
 * Requests one simulation from all nodes in the cluster.
 * This should only be called from the root node.
 *
 * @param withSeeds True when saved seeds should be used.
 */
void MPIStatisticsSimulation::simulate(bool withSeeds) {
	LOG(Priority::DEBUG, "");
	assert(MPICommunicator::theCommunicator().isRoot());
	// Broadcast to slaves that a new simulation is needed.
	int cmd = withSeeds ? SIMULATE_WITH_SEEDS : SIMULATE;
	MPICommunicator::theCommunicator().broadcastInt(&cmd);
	// Do my own simulation.
	StatisticsSimulation::simulatePeriods(withSeeds);
	// Gather the results.
	MPICommunicator::theCommunicator().gatherDouble(&lStatisticsData,
			&lMPIStatisticsData);
	if (lpModel->needScores())
		MPICommunicator::theCommunicator().gatherDouble(&lScoresData,
				&lMPIScoresData);
	if (lpModel->conditional())
		MPICommunicator::theCommunicator().gatherDouble(&lTimesData,
				&lMPITimesData);
	// Distribute results.
	LOGS(Priority::DEBUG)<<"simulation finished, distribute results";
	fireResultModificators(lResult);
	fireResult(lResult);
}

/**
 * \copybrief StatisticsSimulation::updateParameters()
 *
 * Updates the parameters of all simulations in the cluster.
 * This should only be called from the root node.
 *
 * @param rTheta New parameter vector.
 */
void MPIStatisticsSimulation::updateParameters(const VectorXd& rTheta) {
	LOG(Priority::DEBUG, "");
	assert(MPICommunicator::theCommunicator().isRoot());
	// Broadcast command and parameter vector.
	int cmd = UPDATE_PARAMETERS;
	MPICommunicator::theCommunicator().broadcastInt(&cmd);
	// It's save to cast away constness since in root node.
	MPICommunicator::theCommunicator().broadcastVectorXd(
			const_cast<VectorXd*>(&rTheta));
	// Update roots own local model.
	StatisticsSimulation::updateParameters(rTheta);
}

} // namespace siena

#endif // MPI2
