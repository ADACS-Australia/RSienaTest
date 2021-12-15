/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file LogEntry.h
 * \brief Defines the LogEntry class.
 *****************************************************************************/

#ifndef RSIENA_LOGGER_LOG_ENTRY_H_
#define RSIENA_LOGGER_LOG_ENTRY_H_

#include "Priority.h"

#include <ctime>
#include <sstream>

namespace siena {
namespace logger {

/**
 * A log message.
 *
 * Formats and sends the message to the log sinks. A LogEntry is created and
 * destroyed every time one of the LOG* macros is called.
 */
class LogEntry {
public:
	LogEntry(const Priority& priority, const std::tm& time,
			const std::string& file, const std::string& function,
			const int line);
	~LogEntry();

	void message(const char* msg);
	void formatMessage(const char* format, ...);
	std::ostringstream& streamMessage();

	/**
	 * @return Priority of the log entry.
	 */
	const Priority& priority() const {
		return lPriority;
	}

	/**
	 * @return Time when the message was dispatched.
	 */
	const std::tm& time() const {
		return lTime;
	}

	/**
	 * @return File in which the message was dispatched.
	 */
	const std::string& file() const {
		return lFile;
	}

	/**
	 * @return Function in which the message was dispatched.
	 */
	const std::string& function() const {
		return lFunction;
	}

	/**
	 * @return Line in which the message was dispatched.
	 */
	int line() const {
		return lLine;
	}

	/**
	 * @return The actual log message.
	 */
	std::string message() const {
		return lStream.str();
	}

private:
	const Priority& lPriority;

	const std::tm lTime;
	const std::string lFile;
	const std::string lFunction;
	const int lLine;

	std::ostringstream lStream;

};

} // namespace logger
} // namespace siena

#endif // RSIENA_LOGGER_LOG_ENTRY_H_
