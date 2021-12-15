/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file PrefixFilter.h
 * \brief Defines and implements the PrefixFilter class.
 *****************************************************************************/

#ifndef RSIENA_LOGGER_PREFIX_FILTER_H_
#define RSIENA_LOGGER_PREFIX_FILTER_H_

#include "Filter.h"

namespace siena {
namespace logger {

/**
 * Filters log entries by matching the start of the message against a prefix
 * string.
 */
class PrefixFilter: public Filter {
public:

	/**
	 * Constructs a new PrefixFilter object.
	 *
	 * @param rPrefix The message prefix.
	 */
	PrefixFilter(const std::string& rPrefix) :
			lPrefix(rPrefix) {
	}

	virtual ~PrefixFilter() {
	}

	/**
	 * \copydoc Filter::accept()
	 */
	virtual bool accept(const LogEntry& rEntry) const {
		return lPrefix.compare(0, lPrefix.length(), rEntry.message()) == 0;
	}

private:
	PrefixFilter(const PrefixFilter&);
	PrefixFilter& operator=(const PrefixFilter&);

	const std::string lPrefix;

};

} // namespace logger
} // namespace siena

#endif // RSIENA_LOGGER_PREFIX_FILTER_H_
