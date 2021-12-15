/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file UpdateStep.h
 * \brief Defines the Step class.
 *****************************************************************************/

#ifndef RSIENA_STEP_H_
#define RSIENA_STEP_H_

#include "sim/listener/ResultListener.h"
#include "sim/StatisticsSimulation.h"
#include "estimator/listener/EstimatorListener.h"
#include "logger/Logger.h"
#include "estimator/update/step/normalization/StepNormalization.h"
// #include "Eigen/Types.h"

namespace siena {

/**
 * Base class for parameter update steps.
 *
 * Provides several protected methods to sanitize a step.
 */
class UpdateStep: public ResultListener, public EstimatorListener {
public:
	UpdateStep(Simulation& rSimulation, const double gain,
			const StepNormalization& normalization);
	virtual ~UpdateStep() = 0;

	// SimulationListener
	virtual std::vector<ResultType>& needs();
	virtual void onResults(const Result& rResults);

	// EstimatorListener
	virtual void initializePhase();
	virtual void finalizePhase(int nSimulations);
	ResultListener* pResultListener();

	/**
	 * Performs the update step.
	 *
	 * @param rStatisticsMinusTargets Mean of simulated statistics minus target
	 *        statistics.
	 */
	virtual void step(const Eigen::VectorXd& rStatisticsMinusTargets) = 0;

	void setGain(const double gain);
	double getGain() const;

protected:
	//! Simulation of which the parameters will be changed.
	//!
	Simulation& lrSimulation;
	//! Normalization applied to the step.
	//!
	const StepNormalization& lNormalization;
	//! Gain used to fine tune the step width.
	//!
	double lGain;

private:
	// Don't copy it
	UpdateStep(const UpdateStep&);
	UpdateStep& operator=(const UpdateStep&);

};

} // namespace siena

#endif // RSIENA_STEP_H_
