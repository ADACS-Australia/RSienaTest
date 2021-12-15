/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file SDStepNormalization.cpp
 * \brief Implements the SDStepNormalization class.
 *****************************************************************************/

#include "SDStepNormalization.h"
#include "Eigen/Util.h"

using namespace Eigen;

namespace siena {

/**
 * Constructs a StepNormalization object.
 *
 * @param rBasicRate `True` for basic rate parameters.
 * @param rFixed `True` for fixed parameters.
 * @param maxRate Maximum step width.
 * @param rStatistics Full statistics of phase 1.
 */
SDStepNormalization::SDStepNormalization(const Eigen::VectorXb& rBasicRate,
		const Eigen::VectorXb& rFixed, double maxRate,
		const Eigen::MatrixXd& rStatistics) :
		StepNormalization(rBasicRate, rFixed, maxRate), //
		lSd(rBasicRate.size()) {
	for (int c = 0; c < rStatistics.cols(); ++c) {
		double m = mean(rStatistics.col(c));
		lSd[c] = rStatistics.col(c).cwiseProduct(rStatistics.col(c)).sum()
				/ rStatistics.rows() - m * m;
	}
	lSd = lSd.cwiseSqrt();
}

SDStepNormalization::~SDStepNormalization() {
}

/**
 * If any update exceeds lMaxRate the whole step is normalized such that is new
 * maximum is lMaxRate.
 *
 * @param[in,out] rStep The step vector.
 * @param rDeviation The deviation from the target statistics.
 */
void SDStepNormalization::truncate(Eigen::VectorXd& rStep,
		const Eigen::VectorXd& rDeviation) const {
	const double max = lFixed.select(1,
			rDeviation.cwiseAbs().cwiseQuotient(lSd)).maxCoeff();
	if (max > lMaxRate) {
#ifdef R_LEGACY
		rStep = rStep * lMaxRate / max;
#else
		rStep *= lMaxRate / max;
#endif
		LOG(logger::Priority::WARNING, "truncated step");
	}
}

} // namespace siena
