/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file NSimCollector.cpp
 * \brief Implements the NSimCollector class.
 *****************************************************************************/

#include "NSimCollector.h"

#include <algorithm>

using namespace std;
using namespace Eigen;

using namespace std;

namespace siena {

/**
 * Constructs the NSimCollector object.
 */
NSimCollector::NSimCollector() :
		lNSimulations() {
}

/**
 * \copydoc EstimatorListener::initializePhase()
 */
void NSimCollector::initializePhase() {
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void NSimCollector::finalizePhase(int nSimulations) {
	lNSimulations.push_back(nSimulations);
}

/**
 * \copydoc EstimatorListener::pResultListener()
 */
ResultListener* NSimCollector::pResultListener() {
	return 0; // don't need simulations
}

/**
 * Clears the whole vector of number of simulations.
 */
void NSimCollector::clearNSimulations() {
	lNSimulations.clear();
}

/**
 * @return The number of simulated statistics in each phase.
 */
const vector<int>& NSimCollector::rNSimulations() {
	return lNSimulations;
}

} // namespace siena
