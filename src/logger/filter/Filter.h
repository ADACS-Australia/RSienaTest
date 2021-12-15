/*
 * ILogFilter.h
 *
 *  Created on: 02.10.2013
 *      Author: ortmann
 */

#ifndef ILOGFILTER_H_
#define ILOGFILTER_H_

#include "logger/LogEntry.h"

namespace siena {
namespace logger {

/**
 * Log entry filter. The default implementation passes all messages.
 */
class Filter {
public:
	Filter() {
	}

	virtual ~Filter() {
	}

	/**
	 * @param rEntry The log entry to be filtered.
	 * @return `true` if the log entry is accepted.
	 */
	virtual bool accept(const LogEntry& rEntry) const {
		return true;
	}

private:
	Filter(const Filter&);
	Filter& operator=(const Filter&);

};

} // namespace logger
} // namespace siena

#endif // ILOGFILTER_H_
