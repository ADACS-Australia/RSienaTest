/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file GMMDiagonalStep.cpp
 * Implementation of GMMDiagonalStep.h
 *****************************************************************************/

#include "GMMDiagonalStep.h"

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

///////////////////////////////////////////////////////////////////////////////
// Settings
///////////////////////////////////////////////////////////////////////////////

//! Maximum rate of parameter change (before gain).
//! If one element exceeds this value the whole vector is normalized to this.
const double GMMDiagonalStep::MAX_STEP_WIDTH = 5;

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructs the GMMDiagonalStep object.
 *
 * @param rSimulation The simulation updated in the step.
 * @param gain Initial gain.
 * @param rDerivativeDiagonal Main diagonal of the derivative matrix.
 * @param rWeight \( B \) a matrix weighting the GMM effects.
 * @param rNormalization Normalization applied to the step.
 */
GMMDiagonalStep::GMMDiagonalStep(siena::Simulation& rSimulation, const double gain,
		const Eigen::VectorXd& rDerivativeDiagonal, const Eigen::MatrixXd& rWeight,
		const StepNormalization& rNormalization) :
		UpdateStep(rSimulation, gain, rNormalization), //
		lDiagonal(rDerivativeDiagonal), //
		lWeight(rWeight) {
	LOGS(Priority::VERBOSE)<<"\nStep: 1/["<<lDiagonal.transpose()<<"]\n";
}

/**
 * \copydoc UpdateStep::step()
 */
void GMMDiagonalStep::step(const Eigen::VectorXd& rStatisticsMinusTargets) {
	VectorXd step = (lWeight * rStatisticsMinusTargets).cwiseQuotient(lDiagonal);
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

	LOGS(Priority::DEBUG)<<"\nfra: "<<rStatisticsMinusTargets.transpose()
	<<"\nstep: "<<step.transpose()
	<<"\ntheta: "<<theta.transpose();
}

}
 // namespace siena
