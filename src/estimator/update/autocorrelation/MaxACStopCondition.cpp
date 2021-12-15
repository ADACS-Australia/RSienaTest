/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MaxACStopCondition.cpp
 * \brief Implements the MaxACStopCondition class.
 *****************************************************************************/

#include "MaxACStopCondition.h"

using namespace std;

namespace siena {

///////////////////////////////////////////////////////////////////////////////
// Auto correlation conditions
///////////////////////////////////////////////////////////////////////////////

//!
//! R: phase2.r doIterations() (line ~388)
const double MaxACStopCondition::AC_MAX_EPS = 1e-10;

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructs a MaxACStopCondition object.
 *
 * @param rAutoCorrelator The auto correlation.
 * @param nIterations Lower and upper iteration bounds.
 */
MaxACStopCondition::MaxACStopCondition(siena::AutoCorrelator& rAutoCorrelator,
		const std::pair<int, int>& nIterations) :
		lrAC(rAutoCorrelator), //
		lNIterations(nIterations) {
}

MaxACStopCondition::~MaxACStopCondition() {
}

} // namespace siena
