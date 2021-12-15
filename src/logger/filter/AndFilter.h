/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file AndFilter.h
 * \brief Defines and implements the AndFilter class.
 *****************************************************************************/

#ifndef RSIENA_LOGGER_AND_FILTER_H_
#define RSIENA_LOGGER_AND_FILTER_H_

#include "logger/LogEntry.h"

namespace siena {
namespace logger {

/**
 * Concatenation of two filters, letting entries pass only if both filters
 * evaluate to true.
 */
class AndFilter: public Filter {
public:

	/**
	 * Constructs a new AndFilter object.
	 *
	 * @param lhs First filter.
	 * @param rhs Second filter.
	 */
	AndFilter(const Filter& lhs, const Filter& rhs) :
			lrLhs(lhs), //
			lrRhs(rhs) {
	}

	~AndFilter() {
	}

	/**
	 * \copydoc Filter::accept()
	 */
	bool accept(const LogEntry& rEntry) const {
		return lrLhs.accept(rEntry) && lrRhs.accept(rEntry);
	}

private:
	AndFilter(const AndFilter&);
	AndFilter& operator=(const AndFilter&);

	const Filter& lrLhs;
	const Filter& lrRhs;

};

} // namespace logger
} // namespace siena

#endif // RSIENA_LOGGER_AND_FILTER_H_
