/*
 * ConsoleAppender.h
 *
 *  Created on: 04.09.2013
 *      Author: ortmann
 */

#ifndef CONSOLEAPPENDER_H_
#define CONSOLEAPPENDER_H_

#include "Appender.h"

#include <iostream>

namespace siena {
namespace logger {

/**
 * Appends log messages to `std::cout`.
 */
class ConsoleAppender: public Appender {
public:

	/**
	 * \copydoc Appender::Appender()
	 */
	ConsoleAppender(const Filter& rFilter, const Formatter& rFormatter) :
			Appender(rFilter, rFormatter) {
	}

	virtual ~ConsoleAppender() {
	}

protected:

	/**
	 * \copydoc Appender::writeLog()
	 */
	void writeLog(const std::string& msg) {
		std::cout << msg << std::endl;
	}

};

} // namespace logger
} // namespace siena

#endif /* CONSOLEAPPENDER_H_ */
