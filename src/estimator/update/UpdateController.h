/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file UpdateController.h
 * \brief Defines the UpdateController class.
 *****************************************************************************/

#ifndef RSIENA_UPDATE_CONTROLLER_H_
#define RSIENA_UPDATE_CONTROLLER_H_

#include "estimator/Controller.h"
#include "estimator/update/StopCondition.h"
#include "estimator/update/autocorrelation/AutoCorrelator.h"
#include "estimator/components/ParameterSum.h"

#include <utility> // pair

namespace siena {

// Forward declarations
class UpdateStep;

/**
 * Controller for the incremental improvement of the model (phase 2).
 */
class UpdateController: public Controller {
public:
	// Magic numbers
	static const double N2_MIN_BASE;
	static const double N2_MIN_EFFECT_OFFSET;
	static const double N2_MIN_EFFECT_MIN;
	static const int N2_MAX_OFFSET;
	static const int SUBPHASE_MAX_REPEATS;
	// Default parameter values (provided via R)
//	static const int N_SUBPHASES_DEFAULT;
//	static const double GAIN_INITIAL_DEFAULT;
//	static const double GAIN_DECAY_DEFAULT;

	UpdateController(Simulation& rSimulation, int nSubphases, double gainInit,
			double gainDecay);
	~UpdateController();

	bool run();

	void step(UpdateStep* pStep);

	// Variable parameters
	int nSubphases();
	double gainInitial();
	double gainDecay();

protected:
	// Variable parameters
	//! Number of sub phases to run.
	//!
	const int lNSubphases;
	//! Initial gain parameter for the update step.
	//!
	const double lGainInitial;
	//! Gain decay. After each sub phase the gain parameter is multiplied by
	//! value.
	const double lGainDecay;

	std::pair<int, int> calculateNIterations(const int subphase);
	void runSubphase(StopCondition& rCond);
	void retryPhase(const std::pair<int, int>& nIterations);

	//! Iteration counter.
	//!
	int lNIterations;
	//! The step performed in each iteration.
	//!
	UpdateStep* lpStep;
	//! Autocorrelation of the statistics.
	//!
	AutoCorrelator lAutoCorrelator;
	//! Sum of all the parameters used to simulate in the current wave.
	//!
	ParameterSum lParameterSum;

};

} // namespace siena

#endif // RSIENA_UPDATE_CONTROLLER_H_
