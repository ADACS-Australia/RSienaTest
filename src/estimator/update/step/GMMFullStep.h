/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file GMMFullStep.h
 * \brief Defines the GMMFullStep class.
 *****************************************************************************/

#ifndef RSIENA_GMM_FULL_STEP_H_
#define RSIENA_GMM_FULL_STEP_H_

#include "estimator/update/step/UpdateStep.h"

namespace siena {

/**
 * A full Newton-Raphson step.
 *
 * \f[ \theta_{t+1} = \theta_t - a \cdot \hat{D}^{-1} B (E_{\theta} S-s) \f]
 *
 * \( D \) is a square matrix:
 *   - the derivative for non-GMM models.
 *   - \( B\Gamma \) for GMM models.
 * \( B \) is the weighting matrix.
 *   - I for non-GMM models.
 *   - \( \Gamma^{T} W^{-1} \).
 * \(\Gamma\) the (non-square) first order derivative of the statistics.
 * \( W \) is the variance covariance matrix of the statistics.
 */
class GMMFullStep: public UpdateStep {
public:
	// Magic numbers
	static const double MAX_STEP_WIDTH;

	GMMFullStep(Simulation& rSimulation, const double gain,
			const Eigen::MatrixXd& rInverseDerivative,
			const Eigen::MatrixXd& rWeight,
			const StepNormalization& normalization);

	void step(const Eigen::VectorXd& rStatisticsMinusTargets);

	const Eigen::MatrixXd& rWeight();
	const Eigen::MatrixXd& rInverseDerivative();

protected:
	//! GMM weighting matrix.
	//!
	const Eigen::MatrixXd lWeight;
	//! This is a square matrix taking the role of the inverse of the
	//! derivative.
	const Eigen::MatrixXd lInverse;

private:
	// Don't copy it
	GMMFullStep(const GMMFullStep&);
	GMMFullStep& operator=(const GMMFullStep&);

};

} // namespace siena

#endif // RSIENA_GMM_FULL_STEP_H_
