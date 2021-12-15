/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file SDStepNormalization.h
 * \brief Defines the SDStepNormalization class.
 *****************************************************************************/

#ifndef RSIENA_SD_NORMALIZATION_H_
#define RSIENA_SD_NORMALIZATION_H_

#include "StepNormalization.h"

namespace siena
{

/**
 * Step sanitizing with the max step width based on the statistics of phase 1.
 */
class SDStepNormalization: public StepNormalization {
public:
	SDStepNormalization(const Eigen::VectorXb& rBasicRate,
			const Eigen::VectorXb& rFixed, double maxRate,
			const Eigen::MatrixXd& rStatistics);
	virtual ~SDStepNormalization();

	virtual void truncate(Eigen::VectorXd& rStep,
			const Eigen::VectorXd& rDeviation) const;

protected:
	//!
	//! Common name: sd
	Eigen::VectorXd lSd;

private:
	// Don't copy it
	SDStepNormalization(const StepNormalization&);
	SDStepNormalization& operator=(const SDStepNormalization&);

};

} // namespace siena

#endif // RSIENA_SD_NORMALIZATION_H_
