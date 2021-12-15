/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * Triangle counting functions.
 *****************************************************************************/
#ifndef SIMMELIAN_TIE_ITERATOR_H_
#define SIMMELIAN_TIE_ITERATOR_H_

namespace siena {

	// forward declarations
	class OneModeNetwork;

	/**
	 * @param pNetwork the source network.
	 * @returns a pointer to a new network (you take care of delete!). The tie
	 * value is the count of fully reciprocated triangles this tie is involved
	 * in.
	 */
	const OneModeNetwork* newSimmelian(const OneModeNetwork* pNetwork);

	/**
	 * @param pNetwork the directed source network.
	 * @param pAggregate a function pointer.
	 * @returns a pointer to a new network (you take care of delete!). It is
	 * symmetric and constructed by aggregating both directions of a dyad into
	 * one value, which is assign to both directions.
	 */
	const OneModeNetwork* newUndirected(const OneModeNetwork* pNetwork,
			const int& (*pAggregate)(const int&, const int&));

	/**
	 * @param pNetwork the symmetric source network.
	 * @returns a pointer to a new network (you take care of delete!). The tie
	 * value is the count of triangles this tie is involved in.
	 */
	const OneModeNetwork* newIncidentTriads(const OneModeNetwork* pNetwork);

} // namespace siena
#endif // SIMMELIAN_TIE_ITERATOR_H_
