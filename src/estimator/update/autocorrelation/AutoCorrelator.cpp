/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file AutoCorrelator.cpp
 * \brief Implements the AutoCorrelator class.
 *****************************************************************************/

#include "AutoCorrelator.h"

#include <limits>

using namespace std;
using namespace Eigen;

namespace siena {

/**
 * Constructs the AutoCorrelator object.
 *
 * @param nStatistics Size of the vectors.
 */
AutoCorrelator::AutoCorrelator(const int nStatistics) :
		lProd0(nStatistics), //
		lProd1(nStatistics), //
		lAutoCorrelation(nStatistics) {
	initializePhase();
}

/**
 * \copybrief  ResultListener::needs()
 *
 * @return Vector with MEAN_STATISTICS_MINUS_TARGET.
 */
vector<ResultType>& AutoCorrelator::needs() {
	static vector<ResultType> needs(1, MEAN_STATISTICS_MINUS_TARGET);
	return needs;
}

/**
 * \copydoc ResultListener::onResults()
 */
void AutoCorrelator::onResults(const Result& rResults) {
	// Calling the next member function. The value of next toggles between
    // updateEven() and updateOdd(). Using this continuation style saves the
	// even/odd iteration comparison.
	(this->*next)(*rResults.pMeanStatisticsMinusTargets());
}

/**
 * \copybrief EstimatorListener::initializePhase()
 *
 * Reset/start a new sequence.
 */
void AutoCorrelator::initializePhase() {
	lProd0.fill(0);
	lProd1.fill(0);
	// Might query this before the second update is completed.
	lAutoCorrelation.fill(0);
	next = &siena::AutoCorrelator::updateEven;
	// Be sure to be no less than the test eps. Both are tested with '< eps'.
	lMax = std::numeric_limits<double>::max();
	lMin = std::numeric_limits<double>::max();
}

/**
 * \copydoc EstimatorListener::finalizePhase()
 */
void AutoCorrelator::finalizePhase(int) {
	//
}

/**
 * \copydoc EstimatorListener::pResultListener()
 */
ResultListener* AutoCorrelator::pResultListener() {
	return this;
}

/**
 * Stores the vector as previous element in the sequence.
 *
 * @param rV Current vector in the sequence.
 */
void AutoCorrelator::updateEven(const Eigen::VectorXd& rV) {
	lPrevious = rV;
	// Next time catching results call updateOdd().
	next = &siena::AutoCorrelator::updateOdd;
}

/**
 * Updates the auto correlation.
 *
 * Similar to: phase2.r doIterations() (line ~322)
 *
 * @param rV Current vector in the sequence.
 */
void AutoCorrelator::updateOdd(const Eigen::VectorXd& rV) {
#ifdef R_LEGACY
	lProd0 = lProd0 + rV.cwiseProduct(rV);
	lProd1 = lProd1 + rV.cwiseProduct(lPrevious);
#else
	lProd0 += rV.cwiseProduct(rV);
	lProd1 += rV.cwiseProduct(lPrevious);
#endif
	// Division by 0 results in [-]infinity.
	lAutoCorrelation = lProd1.cwiseQuotient(lProd0);

	lMax = std::numeric_limits<double>::min(); // -99
	lMin = std::numeric_limits<double>::max(); // 1

	for (int i = 0; i < lProd0.size(); ++i) {
		if (std::isfinite(lAutoCorrelation[i])) {
			if (lAutoCorrelation[i] > -1) {
				lMin = std::min(lMin, lAutoCorrelation[i]);
			}
			lMax = std::max(lMax, lAutoCorrelation[i]);
		}
	}
	// Next time catching results call updateEven().
	next = &siena::AutoCorrelator::updateEven;
}

/**
 * @return Autocorrelation vector. May contain infinity.
 */
const Eigen::VectorXd& AutoCorrelator::rAutoCorrelation() {
	return lAutoCorrelation;
}

/**
 * @return Minimum finite autocorrelation.
 */
double AutoCorrelator::min() {
	return lMin;
}

/**
 * @return Maximum finite autocorrelation.
 */
double AutoCorrelator::max() {
	return lMax;
}

} // namespace siena
