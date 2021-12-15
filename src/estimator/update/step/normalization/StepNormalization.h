/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file StepNormalization.h
 * \brief Defines the StepNormalization class.
 *****************************************************************************/

#ifndef RSIENA_NORMALIZATION_H_
#define RSIENA_NORMALIZATION_H_

#include "logger/Logger.h"
#include "Eigen/Types.h"

namespace siena {

/**
 * Provides step sanitizing routines.
 */
class StepNormalization {
public:
	// Magic numbers
	static const double POSITIVIZE_FACTOR;

	StepNormalization(const Eigen::VectorXb& rBasicRate,
			const Eigen::VectorXb& rFixed, double maxRate);
	virtual ~StepNormalization();

	virtual void truncate(Eigen::VectorXd& rStep,
			const Eigen::VectorXd& rDeviation) const;
	virtual void positivize(Eigen::VectorXd& rStep,
			const Eigen::VectorXd& rTheta) const;

protected:
	//! `true` for basic rate parameters.
	//!
	const Eigen::VectorXb lBasicRate;
	//! `true` for fixed parameters.
	//!
	const Eigen::VectorXb& lFixed;
	//! Maximum step width.
	//!
	const double lMaxRate;

private:
	// Don't copy it
	StepNormalization(const StepNormalization&);
	StepNormalization& operator=(const StepNormalization&);

};

} // namespace siena

#endif // RSIENA_NORMALIZATION_H_
