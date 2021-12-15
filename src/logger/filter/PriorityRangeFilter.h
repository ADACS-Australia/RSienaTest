/*
 * ILogFilter.h
 *
 *  Created on: 02.10.2013
 *      Author: ortmann
 */

#ifndef RSIENA_LOGGER_PRIORITY_RANGE_FILTER__H_
#define RSIENA_LOGGER_PRIORITY_RANGE_FILTER__H_

#include "Filter.h"

namespace siena {
namespace logger {

/**
 * Lets log messenges with a priority inside a given range pass.
 */
class PriorityRangeFilter: public Filter {
public:

	/**
	 * Constructs a PriorityRangeFilter object.
	 *
	 * @param rLow Lower priority bound.
	 * @param rHigh Upper priority bound.
	 */
	PriorityRangeFilter(const Priority& rLow, const Priority& rHigh =
			Priority::FATAL) :
			lrLow(rLow), //
			lrHigh(rHigh) {
	}

	virtual ~PriorityRangeFilter() {
	}

	/**
	 * \copybrief Filter::accept()
	 *
	 * @param rEntry Log entry.
	 * @return True if the priority of the message is greater of equal to the
	 *         lower bound and lower or equal to the upper bound.
	 */
	virtual bool accept(const LogEntry& rEntry) const {
		return lrLow <= rEntry.priority() && rEntry.priority() <= lrHigh;
	}

private:
	PriorityRangeFilter(const PriorityRangeFilter&);
	PriorityRangeFilter& operator=(const PriorityRangeFilter&);

	const Priority& lrLow;
	const Priority& lrHigh;

};

} // namespace logger
} // namespace siena

#endif // RSIENA_LOGGER_PRIORITY_RANGE_FILTER__H_
