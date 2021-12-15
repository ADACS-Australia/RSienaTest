/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file StopCondition.h
 * \brief Defines the StopCondition class.
 *****************************************************************************/

#ifndef RSIENA_STOP_CONDITION_H_
#define RSIENA_STOP_CONDITION_H_

namespace siena {

/**
 * A generic stop condition.
 */
class StopCondition {
public:
	StopCondition();
	virtual ~StopCondition() = 0;

	/**
	 * @param iteration Current iteration.
	 * @return true if the condition is fulfilled.
	 */
	virtual bool fulfilled(const int iteration) = 0;

private:
	// Don't copy it
	StopCondition(const StopCondition&);
	StopCondition& operator=(const StopCondition&);

};

} // namespace siena

#endif // RSIENA_STOP_CONDITION_H_
