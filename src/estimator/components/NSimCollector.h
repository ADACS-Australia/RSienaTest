/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file NSimCollector.h
 * \brief Defines the NSimCollector class.
 *****************************************************************************/

#ifndef RSIENA_N_SIM_COLLECTOR_H_
#define RSIENA_N_SIM_COLLECTOR_H_

#include "logger/Logger.h"
#include "estimator/listener/EstimatorListener.h"

#include <vector>

namespace siena {

/**
 * Component collecting the number of simulations performed in each phase.
 */
class NSimCollector: public EstimatorListener {
public:
	NSimCollector();

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations);
	ResultListener* pResultListener();

	// NSimCollector
	void clearNSimulations();
	const std::vector<int>& rNSimulations();

protected:
	//! Number of simulations in each phase.
	//!
	std::vector<int> lNSimulations;

private:
	// Don't copy it
	NSimCollector(const NSimCollector&);
	NSimCollector& operator=(const NSimCollector&);

};

} // namespace siena

#endif // RSIENA_N_SIM_COLLECTOR_H_
