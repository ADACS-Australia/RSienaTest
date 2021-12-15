/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file StepNormalization.cpp
 * \brief Implements the StepNormalization class.
 *****************************************************************************/

#include "StepNormalization.h"

using namespace Eigen;

namespace siena {

///////////////////////////////////////////////////////////////////////////////
// Settings
///////////////////////////////////////////////////////////////////////////////

//! If the step for a basic rate effect exceeds its parameter, it is set to the
//! parameter times this value.
const double StepNormalization::POSITIVIZE_FACTOR = 0.5;

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructs a StepNormalization object.
 *
 * @param rBasicRate `True` for basic rate parameters.
 * @param rFixed `True` for fixed parameters.
 * @param maxRate Maximum step width.
 */
StepNormalization::StepNormalization(const Eigen::VectorXb& rBasicRate, const
		Eigen::VectorXb& rFixed, double maxRate) :
		lBasicRate(rBasicRate), //
		lFixed(rFixed), //
		lMaxRate(maxRate) {
}

StepNormalization::~StepNormalization() {
}

/**
 * If any update exceeds lMaxRate the whole step is normalized such that is new
 * maximum is lMaxRate.
 *
 * @param[in,out] rStep The step vector.
 * @param rDeviation The deviation from the target statistics.
 */
void StepNormalization::truncate(Eigen::VectorXd& rStep,
		const Eigen::VectorXd& rDeviation) const {
	// Normalize to lMaxRate if the step exceeds it
	const double max = rStep.cwiseAbs().maxCoeff();
	if (max > lMaxRate) {
#ifdef R_LEGACY
		rStep = rStep * lMaxRate / max;
#else
		rStep *= lMaxRate / max;
#endif
		LOG(logger::Priority::WARNING, "truncated step");
	}
}

/**
 * If the update for basic rate effects exceeds the parameter value it is scaled
 * down to POSITIVIZE_FACTOR times the parameter value.
 *
 * @param[in,out] rStep the step vector.
 * @param rTheta the current parameter vector.
 */
void StepNormalization::positivize(Eigen::VectorXd& rStep,
		const Eigen::VectorXd& rTheta) const {
	// Positivizing: If the step for a basic rate exceeds the current rate, it
	// is set to the half of the current rate. (phase 1 >= ; phase 2 >)
	const Eigen::VectorXb pos = (rStep.array() > rTheta.array()
			&& lBasicRate.array());
	if (pos.any()) {
		rStep = pos.select(rTheta * POSITIVIZE_FACTOR, rStep);
		LOG(logger::Priority::WARNING, "positivized step");
	}
	rStep = lFixed.select(0, rStep);
}

} // namespace siena
