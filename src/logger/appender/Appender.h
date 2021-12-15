/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Appender.h
 *
 *  Created on: 04.09.2013
 *      Author: ortmann
 *****************************************************************************/

#ifndef ABSTRACTAPPENDER_H_
#define ABSTRACTAPPENDER_H_

#include "../formatter/Formatter.h"
#include "../filter/Filter.h"

#include <vector>

namespace siena {
namespace logger {

class LogEntry;

/**
 * Appender base class.
 *
 * Filters and formats the entries and writes them to the sink.
 */
class Appender {
public:
	virtual ~Appender() {
	}

	/**
	 * Filter, format, write a log entry.
	 *
	 * @param rEntry The entry to append.
	 */
	void appendEntry(const LogEntry& rEntry) {
		if (lrFilter.accept(rEntry)) {
			writeLog(lrFormatter.format(rEntry));
		}
	}

protected:

	/**
	 * Constructs an Appender.
	 *
	 * @param rFilter Log entry filter. Only entries passing the filter are
	 *        formatted and appended.
	 * @param rFormatter Formatter to be applied before appending.
	 */
	Appender(const Filter& rFilter, const Formatter& rFormatter) :
			lrFilter(rFilter), //
			lrFormatter(rFormatter) {
	}

	/**
	 * Writes the message.
	 *
	 * @param msg The formatted message.
	 */
	virtual void writeLog(const std::string& msg) = 0;

private:
	Appender(const Appender&);

	const Filter& lrFilter;
	const Formatter& lrFormatter;

};

} // namespace logger
} // namespace siena

#endif // ABSTRACTAPPENDER_H_
