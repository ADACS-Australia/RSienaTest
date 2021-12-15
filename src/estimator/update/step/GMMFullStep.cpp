/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file GMMFullStep.cpp
 * \brief Implements the GMMFullStep class.
 *****************************************************************************/

#include "GMMFullStep.h"

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

///////////////////////////////////////////////////////////////////////////////
// Settings
///////////////////////////////////////////////////////////////////////////////

//! Maximum rate of parameter change (before gain).
//! If one element exceeds this value the whole vector is normalized to this.
const double GMMFullStep::MAX_STEP_WIDTH = 5;

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructs the GMMFullStep object. Use this for GMM models.
 *
 * @param rSimulation The simulation updated in the step.
 * @param gain Initial gain.
 * @param rInverseDerivative \( D^{-1} \) a square matrix.
 * @param rWeight \( B \) a matrix weighting the GMM effects.
 * @param rNormalization Normalization applied to the step.
 */
GMMFullStep::GMMFullStep(siena::Simulation& rSimulation, const double gain,
		const Eigen::MatrixXd& rInverseDerivative, const Eigen::MatrixXd& rWeight,
		const StepNormalization& rNormalization) :
		UpdateStep(rSimulation, gain, rNormalization), //
		lWeight(rWeight), //
		lInverse(rInverseDerivative) {
	LOGS(Priority::VERBOSE)<<"1\n";
	LOGS(Priority::VERBOSE)<<"\nweights:\n"<<lWeight
	<<"\ninverse:\n"<<lInverse<<"\n";
}

/**
 * \copydoc UpdateStep::step()
 */
void GMMFullStep::step(const Eigen::VectorXd& rStatisticsMinusTargets) {
	LOGS(Priority::VERBOSE)<<"step\n";
	VectorXd step = lInverse * (lWeight * rStatisticsMinusTargets);
	lNormalization.truncate(step, rStatisticsMinusTargets);
#ifdef R_LEGACY
	step = step * lGain;
#else
	step *= lGain;
#endif

	VectorXd theta = lrSimulation.rParameters();
	lNormalization.positivize(step, theta);

	theta -= step;
	lrSimulation.updateParameters(theta);

	LOGS(Priority::VERBOSE)<<"\nfra: "<<rStatisticsMinusTargets.transpose()
	<<"\nstep: "<<step.transpose()
	<<"\ntheta: "<<theta.transpose();
}

/**
 * @return The \f[ v' \f] matrix.
 */
const Eigen::MatrixXd& GMMFullStep::rWeight() {
	return lWeight;
}

/**
 * @return The \f[ (v'\Gamma)^{-1} \f] matrix.
 */
const Eigen::MatrixXd& GMMFullStep::rInverseDerivative() {
	return lInverse;
}

} // namespace siena
