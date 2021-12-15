/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file GMMEstimator.h
 * \brief Defines the GMMEstimator class.
 *****************************************************************************/

#ifndef RSIENA_GMM_ESTIMATOR_H_
#define RSIENA_GMM_ESTIMATOR_H_

#include "estimator/SienaFit.h"
#include "estimator/components/NSimCollector.h"
#include "estimator/components/StatisticCollector.h"
#include "estimator/derivative/DifferentiationController.h"
#include "estimator/update/UpdateController.h"
#include "estimator/test/StaticController.h"

#include <Eigen/Core>
#include <list>

namespace siena {

// Forward declarations
class Model;
class Data;
class StepNormalization;
class Differentiation;
class EstimatorListener;

/**
 * Estimator for the Generalized Method of Moments.
 *
 * Starts the 3 controllers for the basic phases and registers further
 * components on them if needed.
 */
class GMMEstimator {
public:
	// Magic numbers
	static const double DERIVATIVE_INV_DIAG_THRESHOLD;
	static const double DERIVATIVE_INV_DIAG_MIN;
	static const double PHASE1_GAIN_MULTIPLICATOR;
	static const int N1_DOLBY_MIN_ITERATIONS;

	GMMEstimator(Simulation& rSimulation, const DifferentiationType& rDiffType,
			const bool dolby, const int nsub, const double gainInitial, const double
			gainDecay, const int n1min, const int n3);
	~GMMEstimator();

	void estimate();
	void estimateFromPrevAns(const SienaFit& prevAns);
	const SienaFit& rFit();

protected:
	//! Simulation used.
	//!
	Simulation& lrSimulation;
	//! Indicates whether the dolby option is used.
	//!
	const bool lDolby;
	//! Component to collect the performed number of simulations.
	//!
	NSimCollector lNSimulations;
	//! Statistics collector used in phase 1 and 3.
	//!
	StatisticCollector lStatisticCollector;
	//! Controller for phase 1.
	//!
	DifferentiationController lDiffController;
	//! Update step for phase 2.
	//!
	StepNormalization* lpNormalization;
	//! GMM weights of phase 1 needed to correct tstats of phase 3.
	//!
	Eigen::MatrixXd lWeights;
	//! Update step for phase 2.
	//!
	UpdateStep* lpStep;
	//! Controller for phase 2.
	//!
	UpdateController lUpdateController;
	//! Controller for phase 3.
	//!
	StaticController lStaticController;
	//! Results.
	//!
	SienaFit lFit;

	bool approximateDerivative();
	bool approximateDerivativeDolby();
	void setupStep(const Eigen::MatrixXd& statistics,
		const Eigen::MatrixXd& derivative, const Eigen::MatrixXd& weight);
	void updateModel();
	void testModel();
	void testConditionalModel();

private:
	// Don't copy it
	GMMEstimator(const GMMEstimator&);
	GMMEstimator& operator=(const GMMEstimator&);

	bool updateDerivativeAndWeight(Eigen::MatrixXd& rDerivative, Eigen::MatrixXd& rWeight,
			const Eigen::MatrixXd& rCovariance);
	bool makeDerivativeInvertible(Eigen::MatrixXd& rDerivative);

	Eigen::MatrixXd createBlockStructure();
	// void makeCovarianceBlockDiagonal(Eigen::MatrixXd& rCovariance);

};

} // namespace siena

#endif // RSIENA_GMM_ESTIMATOR_H_
