/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file DifferentiationController.h
 * \brief Defines the DifferentiationController class.
 *****************************************************************************/

#ifndef RSIENA_DIFF_CONTROLLER_H_
#define RSIENA_DIFF_CONTROLLER_H_

#include "estimator/Controller.h"
#include "estimator/derivative/DifferentiationType.h"

#include "Eigen/Types.h"

namespace siena {

// Forward declarations
class Differentiation;
class FiniteDifference;

/**
 * Controls the approximation of derivatives.
 *
 * There are three differentiation types and specific error handling for each.
 *   - SCORE_DEVIATION: Increases the number of iterations done till a good
 *     approximations is reached. If not fall back to FINITE_DIFFERENCES.
 *   - FINITE_DIFFERENCES: Increase epsilon and restart the phase if the
 *     derivative was not sufficient. Then fix parameters that still have a to
 *     low derivative.
 *   - MAXIMUM_LIKELIHOOD: Just run it.
 */
class DifferentiationController: public Controller {
public:
	// Magic numbers
	static const int FD_REPEATS_FOR_EPSILON;
	static const double FD_EPSILON_LOWER_THRESHOLD;
	static const double FD_EPSILON_GAIN_FAST_BASIC_RATE;
	static const double FD_EPSILON_GAIN_FAST;
	static const double FD_EPSILON_GAIN_SLOW_BASIC_RATE;
	static const double FD_EPSILON_GAIN_SLOW;
	static const int FD_FIX_THRESHOLD;
	static const int SD_MAX_ITERATIONS;
	static const double SD_ITERATION_GAIN;
	static const int N1_OFFSET;
	static const int N1_FACTOR;

	DifferentiationController(Simulation& rSimulation,
			const DifferentiationType& rDiffType, const int minIterations);

	bool run();

	Differentiation* pDifferentiation();
	const Eigen::MatrixXd& rDerivative();
	const Eigen::VectorXb& rFixed();

protected:
	int calculateNIterations(const int minSimulations);

	bool scoreDeviation();

	bool finiteDifferences();
	void increaseEpsilon(FiniteDifference* const pFiniteDifference);

	bool mlDerivative();

	//! The differentiation type. Initially set to the desired type, but changed
	//! when error handling switches to another method.
	DifferentiationType lDiffType;

private:
	//! Number of iterations to run.
	//!
	const int lNIterations;
	//! The differentiation component.
	//!
	Differentiation* lpDifferentiation;
	//! Boolean vector. True for parameters which were fixed as a last
	//! possibility if no derivative could be approximated for them.
	Eigen::VectorXb lFixed;

};

} // namespace siena

#endif // RSIENA_DIFF_CONTROLLER_H_
