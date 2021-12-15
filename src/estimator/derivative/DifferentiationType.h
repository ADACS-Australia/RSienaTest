/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file DifferentiationType.h
 * \brief Defines the DifferentiationType enum.
 *****************************************************************************/

#ifndef RSIENA_DIFFERENTIATION_TYPE_H_
#define RSIENA_DIFFERENTIATION_TYPE_H_

namespace siena {

/**
 * Types of differentiations.
 */
enum DifferentiationType {
	FINITE_DIFFERENCE,
	SCORE_DEVIATION,
	MAXIMUM_LIKELIHOOD
};

} // namespace siena

#endif // RSIENA_DIFFERENTIATION_TYPE_H_
