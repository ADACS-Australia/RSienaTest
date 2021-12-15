/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ParameterSum.h
 * \brief Defines the ParameterSum class.
 *****************************************************************************/

#ifndef RSIENA_PARAMETER_SUM_H_
#define RSIENA_PARAMETER_SUM_H_

#include "logger/Logger.h"
#include "sim/listener/ResultListener.h"
#include "estimator/listener/EstimatorListener.h"

#include <Eigen/Core>

namespace siena {

/**
 * Component summing up all parameters used for simulation.
 */
class ParameterSum: public ResultListener, public EstimatorListener {
public:
	explicit ParameterSum(const int nParameters);

	// SimulationListener
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations);
	ResultListener* pResultListener();

	// ParameterSum
	const Eigen::VectorXd& rThetaSum();

protected:
	//! Sum of parameters.
	//!
	Eigen::VectorXd lrThetaSum;

private:
	// Don't copy it
	ParameterSum(const ParameterSum&);
	ParameterSum& operator=(const ParameterSum&);

};

/**
 * \copydoc ResultListener::onResults()
 *
 * Update parameter sum.
 */
inline void ParameterSum::onResults(const Result& rResults) {
	lrThetaSum += *rResults.pTheta();
	LOGS(logger::Priority::VERBOSE)<<"parameters: "<<rResults.pTheta()->transpose();
	LOGS(logger::Priority::VERBOSE)<<"sum of parameters: "<<lrThetaSum.transpose();
}

} // namespace siena

#endif // RSIENA_PARAMETER_SUM_H_
