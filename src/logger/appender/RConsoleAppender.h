/*
 * RConsoleAppender.h
 *
 *  Created on: 04.09.2013
 *      Author: ortmann
 */

#ifndef RCONSOLEAPPENDER_H_
#define RCONSOLEAPPENDER_H_

#include "Appender.h"

/* #include <iostream> */

namespace siena {
namespace logger {

/**
 * Appends log messages to the r console.
 */
class RConsoleAppender: public Appender {
public:

	/**
	 * \copydoc Appender::Appender()
	 */
	RConsoleAppender(const Filter& rFilter, const Formatter& rFormatter) :
			Appender(rFilter, rFormatter) {
	}

	virtual ~RConsoleAppender() {
	}

protected:

	/**
	 * \copydoc Appender::writeLog()
	 */
	void writeLog(const std::string& msg) {
		Rprintf("%s\n", msg.c_str());
	}

};

} // namespace logger
} // namespace siena

#endif /* RCONSOLEAPPENDER_H_ */
