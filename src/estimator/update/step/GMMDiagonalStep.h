/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file GMMDiagonalStep.h
 * \brief Defines the GMMDiagonalStep class.
 *****************************************************************************/

#ifndef RSIENA_GMM_DIAGONAL_STEP_H_
#define RSIENA_GMM_DIAGONAL_STEP_H_

#include "estimator/update/step/UpdateStep.h"

namespace siena {

/**
 * Robbins Monro step.
 *
 * \f[ \theta_{t+1} = \theta_t - a \cdot (B (E_{\theta} S-s)) / diag(D) \f]
 */
class GMMDiagonalStep: public UpdateStep {
public:
	// Magic numbers
	static const double MAX_STEP_WIDTH;

	GMMDiagonalStep(Simulation& rSimulation, const double gain,
			const Eigen::VectorXd& rDerivativeDiagonal,
			const Eigen::MatrixXd& rWeight,
			const StepNormalization& normalization);

	void step(const Eigen::VectorXd& rStatisticsMinusTargets);

protected:
	//! Main diagonal of derivative.
	//!
	const Eigen::MatrixXd lDiagonal;
	//! GMM weighting matrix.
	//!
	const Eigen::MatrixXd lWeight;

private:
	// Don't copy it
	GMMDiagonalStep(const GMMDiagonalStep&);
	GMMDiagonalStep& operator=(const GMMDiagonalStep&);

};

} // namespace siena

#endif // RSIENA_GMM_DIAGONAL_STEP_H_
