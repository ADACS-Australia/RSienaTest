/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file AutoCorrelator.h
 * \brief Defines the AutoCorrelator class.
 *****************************************************************************/

#ifndef RSIENA_AUTO_CORRELATOR_H_
#define RSIENA_AUTO_CORRELATOR_H_

#include "sim/listener/ResultListener.h"
#include "estimator/listener/EstimatorListener.h"

#include <Eigen/Dense>

namespace siena {

/**
 * Component calculating the autocorrelation of the statistics.
 */
class AutoCorrelator: public ResultListener, public EstimatorListener {
public:
	explicit AutoCorrelator(const int n);

	// SimulationListener
	std::vector<ResultType>& needs();
	void onResults(const Result& rResults);

	// EstimatorListener
	void initializePhase();
	void finalizePhase(int nSimulations);
	ResultListener* pResultListener();

	// AutoCorrelator
	const Eigen::VectorXd& rAutoCorrelation();
	double min();
	double max();

private:
	// Don't copy it
	AutoCorrelator(const AutoCorrelator&);
	AutoCorrelator& operator=(const AutoCorrelator&);

	//! Last element set in even iterations.
	//!
	Eigen::VectorXd lPrevious;
	//! Sum of entry wise products of the current vector with itself.
	//!
	Eigen::VectorXd lProd0;
	//! Sum of entry wise products of the current vector with the previous vector.
	//!
	Eigen::VectorXd lProd1;
	//! lProd1 / lProd0
	//!
	Eigen::VectorXd lAutoCorrelation;
	//! Minimum finite autocorrelation.
	//!
	double lMin;
	//! Maximum finite autocorrelation.
	//!
	double lMax;

	//! Member function pointer to the next update.
	//!
	void (siena::AutoCorrelator::*next)(const Eigen::VectorXd& rV);

	void updateEven(const Eigen::VectorXd& rV);
	void updateOdd(const Eigen::VectorXd& rV);

};

} // namespace siena

#endif // RSIENA_AUTO_CORRELATOR_H_
