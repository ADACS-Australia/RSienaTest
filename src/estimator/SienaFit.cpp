/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file SienaFit.cpp
 * \brief Implements the SienaFit class.
 *****************************************************************************/

#include "SienaFit.h"

using namespace std;
using namespace Eigen;

namespace siena {

/**
 * Constructs the SienaFit object.
 */
SienaFit::SienaFit() {
}

SienaFit::~SienaFit() {
}

/**
 * @return Number of simulations done in phase 1.
 */
const vector<int>& SienaFit::rPhase1NSimulations() const {
	return lPhase1NSimulations;
}

/**
 * Sets the number of simulations done in phase 1.
 *
 * @param rPhase1NSimulations Number of simulations.
 */
void SienaFit::phase1NSimulations(const vector<int>& rPhase1NSimulations) {
	lPhase1NSimulations = rPhase1NSimulations;
}

/**
 * @return Number of simulations done in each subphase of phase 2.
 */
const vector<int>& SienaFit::rPhase2NSimulations() const {
	return lPhase2NSimulations;
}

/**
 * Sets the number of simulations done in phase 2.
 *
 * @param rPhase2NSimulations Number of simulations done in each subphase of
 *        phase 2.
 */
void SienaFit::phase2NSimulations(const vector<int>& rPhase2NSimulations) {
	lPhase2NSimulations = rPhase2NSimulations;
}

/**
 * @return Number of simulations done in phase 3.
 */
int SienaFit::phase3NSimulations() const {
	return lPhase3NSimulations;
}

/**
 * Sets the number of simulations done in phase 3.
 *
 * @param phase3NSimulations Number of simulations done in phase 3.
 */
void SienaFit::phase3NSimulations(const int phase3NSimulations) {
	lPhase3NSimulations = phase3NSimulations;
}

/**
 * @return Fixed parameters.
 */
const VectorXb& SienaFit::rFixed() const {
	return lFixed;
}

/**
 * @param rFixed Fixed parameters.
 */
void SienaFit::fixed(const VectorXb& rFixed) {
	lFixed = rFixed;
}

/**
 * @return Estimated parameters.
 */
const VectorXd& SienaFit::rTheta() const {
	return lTheta;
}

/**
 * @param rTheta Estimated parameters.
 */
void SienaFit::theta(const VectorXd& rTheta) {
	lTheta = rTheta;
}

/**
 * @return Target statistics.
 */
const VectorXd& SienaFit::rTargets() const {
	return lTargets;
}

/**
 * @param rTargets Target statistics.
 */
void SienaFit::targets(const VectorXd& rTargets) {
	lTargets = rTargets;
}

/**
 * @return Estimated parameters.
 */
const VectorXd& SienaFit::rTRatio() const {
	return lTRatio;
}

/**
 * @param rTRatio Estimated parameters.
 */
void SienaFit::tRatio(const VectorXd& rTRatio) {
	lTRatio = rTRatio;
}

/**
 * @return Estimated rate parameters.
 */
const VectorXd& SienaFit::rRate() const {
	return lRate;
}

/**
 * @param rRate Estimated rate parameters.
 */
void SienaFit::rate(const VectorXd& rRate) {
	lRate = rRate;
}

/**
 * @return Error of estimated rate parameters.
 */
const VectorXd& SienaFit::rRateError() const {
	return lRateError;
}

/**
 * @param rRateError Error of estimated rate parameters.
 */
void SienaFit::rateError(const VectorXd& rRateError) {
	lRateError = rRateError;
}

const MatrixXd& SienaFit::rBlockStructure() const {
	return lBlockStructure;
}

void SienaFit::blockStructure(const MatrixXd& rBlockStructure) {
	lBlockStructure = rBlockStructure;
}

/**
 * @return Original derivative of statistics.
 */
const MatrixXd& SienaFit::rGamma() const {
	return lPhase3Gamma;
}

/**
 * @param rGamma Original derivative of statistics.
 */
void SienaFit::gamma(const MatrixXd& rGamma) {
	lPhase3Gamma = rGamma;
}

/**
 * @return GMM Weights.
 */
const MatrixXd& SienaFit::rPhase1Weight() const {
	return lPhase1Weight;
}

/**
 * @param rWeight GMM Weights.
 */
void SienaFit::phase1Weight(const MatrixXd& rWeight) {
	lPhase1Weight = rWeight;
}

/**
 * @return GMM Weights.
 */
const MatrixXd& SienaFit::rPhase3Weight() const {
	return lPhase3Weight;
}

/**
 * @param rWeight GMM Weights.
 */
void SienaFit::phase3Weight(const MatrixXd& rWeight) {
	lPhase3Weight = rWeight;
}

/**
 * @return Derivative of statistics.
 */
const MatrixXd& SienaFit::rPhase1Derivative() const {
	return lPhase1Derivative;
}

/**
 * @param rDerivative Derivative of statistics.
 */
void SienaFit::phase1Derivative(const MatrixXd& rDerivative) {
	lPhase1Derivative = rDerivative;
}

/**
 * @return Derivative of statistics.
 */
const MatrixXd& SienaFit::rPhase3Derivative() const {
	return lPhase3Derivative;
}

/**
 * @param rDerivative Derivative of statistics.
 */
void SienaFit::phase3Derivative(const MatrixXd& rDerivative) {
	lPhase3Derivative = rDerivative;
}

/**
 * @return Statistics.
 */
const MatrixXd& SienaFit::rStatistics() const {
	return lPhase3Statistics;
}

/**
 * @param rStatistics statistics.
 */
void SienaFit::statistics(const MatrixXd& rStatistics) {
	lPhase3Statistics = rStatistics;
}

const vector<MatrixXd>& SienaFit::rPeriodWiseStatistics() const {
	return lPeriodWiseStatistics;
}
void SienaFit::periodWiseStatistics(const vector<MatrixXd>& rPeriodWiseStatistics) {
	lPeriodWiseStatistics = rPeriodWiseStatistics;
}

const MatrixXd& SienaFit::rPeriodWiseTargets() const {
	return lPeriodWiseTargets;
}
void SienaFit::periodWiseTargets(const MatrixXd& rPeriodWiseTargets) {
	lPeriodWiseTargets = rPeriodWiseTargets;
}
const vector<MatrixXd>& SienaFit::rPeriodWiseScores() const {
	return lPeriodWiseScores;
}
void SienaFit::periodWiseScores(const vector<MatrixXd>& rPeriodWiseScores) {
	lPeriodWiseScores = rPeriodWiseScores;
}

/**
 * @return Covariance matrix of the statistics.
 */
const MatrixXd& SienaFit::rPhase1Covariance() const {
	return lPhase1Covariance;
}

/**
 * @param rCovariance Covariance matrix of the statistics.
 */
void SienaFit::phase1Covariance(const MatrixXd& rCovariance) {
	lPhase1Covariance = rCovariance;
}

/**
 * @return Covariance matrix of the statistics.
 */
const MatrixXd& SienaFit::rPhase3Covariance() const {
	return lPhase3Covariance;;
}

/**
 * @param rCovariance Covariance matrix of the statistics.
 */
void SienaFit::phase3Covariance(const MatrixXd& rCovariance) {
	lPhase3Covariance = rCovariance;
}

} // namespace siena
