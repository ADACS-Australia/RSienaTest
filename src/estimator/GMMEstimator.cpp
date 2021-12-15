/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file GMMEstimator.cpp
 * Implementation of GMMEstimator.h
 *****************************************************************************/

#include "GMMEstimator.h"

#include "Eigen/Types.h"
#include "Eigen/Util.h"
#include "components/ScoreCollector.h"
#include "components/PeriodWiseScoreCollector.h"
#include "components/PeriodWiseStatisticsCollector.h"
#include "components/TimesCollector.h"
#include "derivative/Differentiation.h"
#include "derivative/DifferentiationController.h"
#include "listener/EstimatorListener.h"
#include "logger/Logger.h"
#include "model/EffectInfo.h"
#include "model/effects/EffectFactory.h"
#include "sim/StatisticsSimulation.h"
#include "sim/modificator/DolbyModificator.h"
#include "update/UpdateController.h"
#include "update/step/GMMDiagonalStep.h"
#include "update/step/GMMFullStep.h"
#include "update/step/normalization/SDStepNormalization.h"
#ifdef MPI2
#include "sim/mpi/MPIStatisticsSimulation.h"
#endif

#include <Eigen/Dense>
#include <R_ext/Utils.h>

#include <algorithm>
#include <list>
#include <numeric>

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

///////////////////////////////////////////////////////////////////////////////
// Phase 1 fine tuning
///////////////////////////////////////////////////////////////////////////////

//! If a diagonal element of the derivative is below this threshold before
//! inverting it, it is increased.
//! R: phase1.r (line ~202)
const double GMMEstimator::DERIVATIVE_INV_DIAG_THRESHOLD = 1e-8;
//! When a diagonal element was too low it is set to this value.
//! R: phase1.r (line ~203)
const double GMMEstimator::DERIVATIVE_INV_DIAG_MIN = 1e-3;
//! This times the initial gain of phase 2, it the gain used for the step at
//! the end of phase 1.
//! R: phase1.r phase1.2 (line ~257)
const double GMMEstimator::PHASE1_GAIN_MULTIPLICATOR = .5;
//! Minimum nuber of simulations for phase 1 if dolby is true.
//! R: robmon.r robmon() (line ~159)
const int GMMEstimator::N1_DOLBY_MIN_ITERATIONS = 50;

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructs a GMMEstimator objects.
 *
 * @param rSimulation Simulation.
 * @param rDiffType Prefered differentiation type. Error handling might choose
 *        a different type.
 * @param dolby If `true` the dolby options is used.
 * @param nsub Number of subphases of phase 2.
 * @param gainInitial Initial gain level of the update step.
 * @param gainDecay Gain decay factor applied after each subphase.
 * @param n3 Number of simulations to be performed in phase 3.
 */
GMMEstimator::GMMEstimator(Simulation& rSimulation,
		const DifferentiationType& rDiffType, const bool dolby, const int nsub,
		const double gainInitial, const double gainDecay, const int n1min,
		const int n3) :
		lrSimulation(rSimulation), //
		lDolby(dolby), //
		lNSimulations(), //
		lDiffController(rSimulation, rDiffType,
				n1min <= 0 ? (dolby ? N1_DOLBY_MIN_ITERATIONS : 0) : n1min), //
		lpNormalization(0/*ptr*/), //
		lpStep(0/*ptr*/), //
		lUpdateController(rSimulation, nsub, gainInitial, gainDecay), //
		lStaticController(rSimulation, n3) {
	LOGS(Priority::DEBUG);
}

GMMEstimator::~GMMEstimator() {
	delete lpStep;
	delete lpNormalization;
	lpStep = 0/*ptr*/;
	lpNormalization = 0/*ptr*/;
}

/**
 * Does the estimation.
 */
void GMMEstimator::estimate() {
	LOGS(Priority::DEBUG);
	// Approximate derivative, if it fails jump to phase 3.
	bool doPhase2;
	if (lDolby) {
		doPhase2 = approximateDerivativeDolby();
	} else {
		doPhase2 = approximateDerivative();
	}
	// Phase 2 only if we have a derivative.
	if (doPhase2) {
		updateModel();
	}
	if (lrSimulation.pModel()->conditional()) {
		testConditionalModel();
	} else {
		testModel();
	}
}

/**
 * Estimate but skip the first phase using phase 3 from the previous results.
 */
void GMMEstimator::estimateFromPrevAns(const SienaFit& prevAns) {
	// fake phase 1
	lWeights = prevAns.rPhase3Weight();
	lFit.phase1NSimulations(vector<int>(1, 0));
	lFit.fixed(prevAns.rFixed());
	lNSimulations.clearNSimulations();

	setupStep(prevAns.rStatistics(), prevAns.rPhase3Derivative(), prevAns.rPhase3Weight());
	// rest of the estimation
	updateModel();
	if (lrSimulation.pModel()->conditional()) {
		testConditionalModel();
	} else {
		testModel();
	}
}

/**
 * Approximation of derivatives (phase 1).
 *
 * @return `true` if a good approximation of the phase3Derivative was found.
 */
bool GMMEstimator::approximateDerivative() {
	LOGS(Priority::INFO)<<"\n\n=== start phase 1 ============"
	<<"\nparameters: "<<lrSimulation.rParameters().transpose()
	<<"\ntargets: "<<lrSimulation.rTargets().transpose()
	<<"\n==============================\n";

	lDiffController.addListener(&lStatisticCollector);
	lDiffController.addListener(&lNSimulations);
	bool ok = lDiffController.run();
	lDiffController.removeListener(&lNSimulations);
	lDiffController.removeListener(&lStatisticCollector);

	if (!ok) {
		return false;
	}

	LOGS(Priority::INFO)<<"=== end phase 1 ==============";
	const MatrixXd statistics = lStatisticCollector.rStatistics().rowwise()
	- lrSimulation.rTargets().transpose();
	LOGS(Priority::VERBOSE)<<"statistics:\n"<<statistics<<"\n";
	// Common name: mnfra
	const VectorXd statisticsMean = statistics.colwise().mean();
	LOGS(Priority::VERBOSE)<<"mean of statistics: "<<statisticsMean.transpose();
	// Common name: msf, W
	MatrixXd statisticsCovariance = covarianceMatrix(statistics);
	LOGS(Priority::INFO)<<"covariance of statistics:\n"
	<<statisticsCovariance<<"\n";
	if (lrSimulation.pModel()->gmmModel()) {
		// makeCovarianceBlockDiagonal(statisticsCovariance);
		//statisticsCovariance = statisticsCovariance.array() * createBlockStructure().array();
		statisticsCovariance = statisticsCovariance.array();
		LOGS(Priority::INFO) << "blockwise covariance of statistics:\n"
				<< statisticsCovariance << "\n";
	}

	// Construct square derivative and (non-square) weighting matrix.
	MatrixXd derivative = lDiffController.rDerivative();
	LOGS(Priority::VERBOSE)<<"gamma (original phase3Derivative):\n"<<derivative<<"\n";
	MatrixXd weight;
	if (!updateDerivativeAndWeight(derivative, weight, statisticsCovariance)) {
		// If inverting is not possible give up and jump directly to phase 3.
		return false;
	}
	LOGS(Priority::INFO)<<"Derivative:\n"<<derivative<<"\n";
	LOGS(Priority::INFO)<<"weight:\n"<<weight<<"\n";
	// Correct variance covariance matrix.
	const MatrixXd& inverseDerivative = derivative.inverse();
	statisticsCovariance = weight * statisticsCovariance * weight.transpose();
	LOGS(Priority::VERBOSE)<<"(B \\Sigma B^T):\n"<<statisticsCovariance;
	statisticsCovariance = (inverseDerivative * statisticsCovariance * inverseDerivative.transpose());

	lWeights = weight;

	lFit.phase1NSimulations(lNSimulations.rNSimulations());
	lFit.fixed(lDiffController.rFixed());
	lFit.targets(lrSimulation.rTargets());
	
	lFit.periodWiseTargets(lrSimulation.rPeriodWiseTargets()); // !!!
	
	lFit.phase1Weight(weight);
	lFit.phase1Derivative(derivative);
	lFit.phase1Covariance(statisticsCovariance);

	setupStep(statistics, derivative, weight);

	LOGS(Priority::VERBOSE)<<"parameters after step: "
	<<lrSimulation.rParameters().transpose();
	LOGS(Priority::INFO)<<"==============================\n";

	lNSimulations.clearNSimulations();
//	printMemUsage();
	return true;
}

void GMMEstimator::setupStep(const MatrixXd& statistics,
		const MatrixXd& derivative, const MatrixXd& weight) {
	delete lpStep;
	delete lpNormalization;
	// Perform the Newton-Raphson step.
	StepNormalization n(lrSimulation.rSimulationEffects().array() ==
			&Simulation::RATE_EFFECT, lFit.rFixed(), GMMFullStep::MAX_STEP_WIDTH);
	GMMFullStep step(lrSimulation, lUpdateController.gainInitial() *
			PHASE1_GAIN_MULTIPLICATOR, derivative.inverse(), weight, n);
	step.step(statistics.colwise().mean());

	// Setup the Robbins-Monro step used in phase 2.
	if (lrSimulation.pModel()->gmmModel()) {
		lpNormalization = new StepNormalization(lrSimulation.rSimulationEffects().array() ==
			&Simulation::RATE_EFFECT, lFit.rFixed(), GMMFullStep::MAX_STEP_WIDTH);
	} else {
		lpNormalization = new SDStepNormalization(lrSimulation.rSimulationEffects().array() ==
				&Simulation::RATE_EFFECT, lFit.rFixed(), GMMDiagonalStep::MAX_STEP_WIDTH, statistics);
	}
	lpStep = new GMMDiagonalStep(lrSimulation, lUpdateController.gainInitial(),
			derivative.diagonal(), weight, *lpNormalization);
	lUpdateController.step(lpStep);
}


/**
 * Approximation of derivatives. Additionally collect the score regression
 * needed for the dolby option.
 *
 * @return `true` if a good approximation of the phase3Derivative was found.
 */
bool GMMEstimator::approximateDerivativeDolby() {
	ScoreCollector scoreCollector;
	lDiffController.addListener(&scoreCollector);
	bool ok = approximateDerivative();
	lDiffController.removeListener(&scoreCollector);
	static DolbyModificator dolby(lStatisticCollector.rStatistics(),
			scoreCollector.rScores());
	lrSimulation.addResultModificator(&dolby);

	LOGS(Priority::INFO)<<"\n\n=== end phase 1 =============="
	<<"\nScore regession:\n"<<dolby.rRegressionCoefficient()
	<<"\n==============================\n";
	return ok;
}

/**
 * Incremental update of the model (phase 2).
 */
void GMMEstimator::updateModel() {
	lUpdateController.addListener(&lNSimulations);
	lUpdateController.run();
	lUpdateController.removeListener(&lNSimulations);
	lFit.phase2NSimulations(lNSimulations.rNSimulations());
	lNSimulations.clearNSimulations();
}

/**
 * Test the estimated parameters (phase 3).
 */
void GMMEstimator::testModel() {
	LOGS(Priority::INFO)<<"\n\n=== start phase 3 ============"
	<<"\nParameters: "<<lrSimulation.rParameters().transpose()
	<<"\n==============================\n";

	// Collect period-wise statistics and scores to export to R for GOF
	PeriodWiseStatisticsCollector periodStatistics;
	PeriodWiseScoreCollector periodScores;
	lStaticController.addListener(&periodStatistics);
	lStaticController.addListener(&periodScores);
	// Simulate and collect data.
	lStaticController.addListener(&lStatisticCollector);
	lStaticController.addListener(&lNSimulations);
	lStaticController.addListener(lDiffController.pDifferentiation());
	lStaticController.run();
	lStaticController.removeListener(lDiffController.pDifferentiation());
	lStaticController.removeListener(&lNSimulations);
	lStaticController.removeListener(&lStatisticCollector);
	lStaticController.removeListener(&periodScores);
	lStaticController.removeListener(&periodStatistics);

	LOGS(Priority::INFO)<<"=== end phase 3 ==============";
	LOGS(Priority::VERBOSE)<<"parameters:\n"
	<<lrSimulation.rParameters().transpose();
	// Common name: sf
	const MatrixXd statistics = lStatisticCollector.rStatistics().rowwise()
	- lrSimulation.rTargets().transpose();
	LOGS(Priority::VERBOSE)<<"statistics:\n"<<statistics;
	// Common name: mnfra
	const VectorXd statisticsMean = statistics.colwise().mean();
	LOGS(Priority::VERBOSE)<<"mean of statistics:\n"
	<<statisticsMean.transpose();
	// Common name: msf
	MatrixXd statisticsCovariance = covarianceMatrix(statistics);
	LOGS(Priority::INFO)<<"covariance of statistics:\n"<<statisticsCovariance;
	MatrixXd blockStructure(0, 0);
	if (lrSimulation.pModel()->gmmModel()) {
		// makeCovarianceBlockDiagonal(statisticsCovariance);
		// blockStructure = createBlockStructure();
		//statisticsCovariance = statisticsCovariance.array() * blockStructure.array();
		blockStructure = MatrixXd::Ones(lrSimulation.rStatisticEffects().size(),
			lrSimulation.rStatisticEffects().size());
		statisticsCovariance = statisticsCovariance.array();
		LOGS(Priority::INFO) << "blockwise covariance of statistics:\n"
				<< statisticsCovariance << "\n";
	}
	// t-ratios taking the old covariance and phase3Derivative values
	// Common name: tstat
	const VectorXd tstat = lWeights * statisticsMean.cwiseQuotient(
			statisticsCovariance.diagonal().cwiseSqrt());
	LOGS(Priority::INFO)<<"tstat:\n"<<tstat.transpose();
	// Common name: se
	VectorXd statisticsError = statisticsCovariance.diagonal().cwiseSqrt();
	LOGS(Priority::INFO)<<"error:\n"<<statisticsError.transpose();
	// Update derivative and weight.
	MatrixXd derivative = lDiffController.pDifferentiation()->rDerivative();
	LOGS(Priority::VERBOSE)<<"gamma (original phase3Derivative):\n"<<derivative<<"\n";
	MatrixXd weight; // Common name: B
	updateDerivativeAndWeight(derivative, weight, statisticsCovariance);
	LOGS(Priority::INFO)<<"phase3Derivative:\n"<<derivative<<"\n";
	LOGS(Priority::INFO)<<"phase3weight:\n"<<weight<<"\n";
	// Correct variance covariance matrix.
	const MatrixXd& inverseDerivative = derivative.inverse();
	statisticsCovariance = weight * statisticsCovariance * weight.transpose();
	LOGS(Priority::VERBOSE)<<"(B \\Sigma B^T):\n"<<statisticsCovariance;
	statisticsCovariance = (inverseDerivative * statisticsCovariance * inverseDerivative.transpose());
	LOGS(Priority::INFO)<<"corrected covariance of statistics:\n"<<statisticsCovariance;
	// Common name: se
	statisticsError = statisticsCovariance.diagonal().cwiseSqrt();
	LOGS(Priority::INFO)<<"standard errors:\n"<<statisticsError.transpose();

	// const MatrixXd statisticsCorrelation =
	// (statisticsCovariance.array().rowwise()
	// 		/ statisticsError.array().transpose()).colwise()
	// / statisticsError.array();
	// MatrixXd covcor = statisticsCovariance;
	// covcor.triangularView<StrictlyLower>() =
	// statisticsCorrelation.triangularView<StrictlyLower>();

	LOGS(Priority::INFO)<<"==============================";

	// Load the SienaFit object.
	lFit.phase3NSimulations(lNSimulations.rNSimulations()[0]);
	lFit.theta(lrSimulation.rParameters());
	lFit.tRatio(tstat);
	lFit.blockStructure(blockStructure);
	lFit.gamma(lDiffController.pDifferentiation()->rDerivative());
	lFit.phase3Weight(weight);
	lFit.phase3Derivative(derivative);
	lFit.statistics(statistics);
	lFit.phase3Covariance(statisticsCovariance);
	lFit.periodWiseStatistics(periodStatistics.rStatistics());
	lFit.periodWiseScores(periodScores.rScores());
}

/**
 * Performs the normal test, additionally collects times to approximate the
 * rate for conditional simulations.
 */
void GMMEstimator::testConditionalModel() {
	TimesCollector timesCollector;
	lStaticController.addListener(&timesCollector);
	testModel();
	lStaticController.removeListener(&timesCollector);

	const VectorXd rate = timesCollector.rTimes().colwise().mean();
	VectorXd rateError(rate.rows());
	for (int i = 0; i < rate.size(); ++i) {
		rateError[i] = sqrt(variance(timesCollector.rTimes().col(i)));
	}

	lFit.rate(rate);
	lFit.rateError(rateError);
}

/**
 * Prepare phase3Derivative and GMM weights.
 *
 * @param[in,out] rPhase3Derivative the phase3Derivative.
 * @param[in,out] rPhase3Weight the GMM weighting matrix.
 * @param[in] rCovariance variance covariance matrix.
 */
bool GMMEstimator::updateDerivativeAndWeight(MatrixXd& rDerivative, MatrixXd& rWeight,
		const MatrixXd& rCovariance) {
	if (lrSimulation.pModel()->gmmModel()) {
		rWeight = (rDerivative.transpose() * rCovariance.inverse());
		LOGS(Priority::DEBUG) << "B = Gamma * Cov^-1:\n" << rWeight << "\n";
		rWeight = (rWeight.rowwise().norm().eval().asDiagonal().inverse()*rWeight);// Common name: B
		LOGS(Priority::DEBUG) << "B' = normalized:\n" << rWeight << "\n";
		rDerivative = (rWeight * rDerivative); // Common name: D, B\Gamma
	} else {
		if (!makeDerivativeInvertible(rDerivative)) {
			// If inverting is not possible give up and jump directly to phase 3.
			return false;
		}
		rWeight = MatrixXd::Identity(rDerivative.rows(), rDerivative.rows()); // Common name: B
	}
	return true;
}

/**
 * Tests if the phase3Derivative is invertible, if not tries to make it so by
 * increasing the main diagonal entries.
 *
 * @param[in,out] rPhase3Derivative the phase3Derivative.
 */
bool GMMEstimator::makeDerivativeInvertible(MatrixXd& rDerivative) {
	ArrayXb nearZero = rDerivative.diagonal().array()
			< DERIVATIVE_INV_DIAG_THRESHOLD;
	if (nearZero.any()) {
		LOG(Priority::WARNING, "Increased some near-zero diagonal elements");
		rDerivative.diagonal() = nearZero.select(DERIVATIVE_INV_DIAG_MIN,
				rDerivative.diagonal());
	}
	if (!FullPivLU<MatrixXd>(rDerivative).isInvertible()) {
		LOG(Priority::WARNING, "Derivative is not invertible, added identity");
		rDerivative.diagonal().array() += 1;
		if (!FullPivLU<MatrixXd>(rDerivative).isInvertible()) {
			LOG(Priority::WARNING, "Derivative is still not invertible");
			return false;
		}
	}
	return true;
}

// /**
//  * Sets variances for effects not in the same GMM group to 0.
//  *
//  * @param[in,out] rPhase3Derivative the phase3Derivative.
//  */
// void GMMEstimator::makeCovarianceBlockDiagonal(Eigen::MatrixXd &rCovariance) {
// 	// can't be sure that the order is right, so were goining to do zero out element by element
// 	for (int r = 0; r < lrSimulation.rStatisticEffects().size(); r++) {
// 		const string rgroup = EffectFactory::gmmGroup(lrSimulation.rStatisticEffects()[r]);
// 		for (int c = 0; c < lrSimulation.rStatisticEffects().size(); c++) {
// 			const string cgroup = EffectFactory::gmmGroup(lrSimulation.rStatisticEffects()[c]);
// 			if (!((rgroup == "*" || cgroup == "*") || // if either is star, leave it
// 					(rgroup != "" && rgroup == cgroup) || // if non-empty match, leave it
// 					(r == c))) {                          // if index match, leave it
// 				rCovariance(r, c) = 0;
// 			}
// 		}
// 	}
// }

MatrixXd GMMEstimator::createBlockStructure() {
	MatrixXd structure = MatrixXd::Ones(lrSimulation.rStatisticEffects().size(),
			lrSimulation.rStatisticEffects().size());
	// can't be sure that the order is right, so were goining to do zero out element by element
	for (int r = 0; r < lrSimulation.rStatisticEffects().size(); r++) {
		const string rgroup = EffectFactory::gmmGroup(lrSimulation.rStatisticEffects()[r]);
		for (int c = 0; c < lrSimulation.rStatisticEffects().size(); c++) {
			const string cgroup = EffectFactory::gmmGroup(lrSimulation.rStatisticEffects()[c]);
			if (!((rgroup == "*" || cgroup == "*") || // if either is star, leave it
					(rgroup != "" && rgroup == cgroup) || // if non-empty match, leave it
					(r == c))) {                          // if index match, leave it
				structure(r, c) = 0;
			}
		}
	}
	return structure;
}

/**
 * @return A SienaFit object holding the estimated parameters and other
 *         information.
 */
const SienaFit& GMMEstimator::rFit() {
	return lFit;
}
} // namespace siena
