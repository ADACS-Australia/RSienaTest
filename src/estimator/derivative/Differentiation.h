/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Differentiation.h
 * \brief Defines the Function class.
 *****************************************************************************/

#ifndef RSIENA_DIFFERENTIATION_H_
#define RSIENA_DIFFERENTIATION_H_

#include "sim/listener/ResultListener.h"
#include "estimator/listener/EstimatorListener.h"

#include <Eigen/Core>

namespace siena {

/**
 * Base class for differentiation methods.
 */
class Differentiation: public ResultListener, public EstimatorListener {
public:
	Differentiation();
	virtual ~Differentiation() = 0;

	ResultListener* pResultListener();

	/**
	 * Calculates the number of simulations to be performed before the derivative
	 * check.
	 *
	 * Similar to: phase1.r phase1.1() (line ~30)
	 *
	 * Common name: firstNit
	 *
	 * @param n Total number of simulations to be performed in phase 1.
	 * @return Simulations before the derivative check in phase 1.
	 */
	virtual int simulationsBeforeCheck(const int n) const = 0;

	/**
	 * Performs the derivative check.
	 *
	 * This should be called after simulationsBeforeCheck() simulations.
	 *
	 * @return True if the derivative approximation is good,
	 */
	virtual bool checkDerivative() const = 0;

	/**
	 * Calculates the derivative.
	 *
	 * @return Derivative matrix.
	 */
	virtual const Eigen::MatrixXd& rDerivative() const = 0;

private:
	// Don't copy it
	Differentiation(const Differentiation&);
	Differentiation& operator=(const Differentiation&);

};

} // namespace siena

#endif // RSIENA_DIFFERENTIATION_H_
