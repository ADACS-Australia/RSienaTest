/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Differentiation.cpp
 * \brief Implements the Differentiation class.
 *****************************************************************************/

#include "Differentiation.h"

namespace siena {

Differentiation::Differentiation() {
}

Differentiation::~Differentiation() {
}

/**
 * \copydoc EstimatorListener::pResultListener()
 */
ResultListener* Differentiation::pResultListener() {
	return this;
}

} // namespace siena
