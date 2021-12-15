/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MinMaxACStopCondition.cpp
 * \brief Implements the MinMaxACStopCondition class.
 *****************************************************************************/

#include "MinMaxACStopCondition.h"

using namespace std;

namespace siena {

///////////////////////////////////////////////////////////////////////////////
// Auto correlation conditions
///////////////////////////////////////////////////////////////////////////////

//!
//! R: phase2.r doIterations() (line ~390)
const double MinMaxACStopCondition::AC_MIN_EPS = -.8;
//!
//! R: phase2.r doIterations() (line ~390)
const int MinMaxACStopCondition::AC_MIN_NITERATIONS = 50;

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructs a MinMaxACStopCondition object.
 *
 * @param rAutoCorrelator The auto correlation.
 * @param nIterations Lower and upper iteration bounds.
 */
MinMaxACStopCondition::MinMaxACStopCondition(siena::AutoCorrelator& rAutoCorrelator,
		const std::pair<int, int>& nIterations) :
		MaxACStopCondition(rAutoCorrelator, nIterations) {
}

} // namespace siena
