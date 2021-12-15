/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file LogEntry.cpp
 * \brief Implements the LogEntry class.
 *****************************************************************************/

#include "LogEntry.h"
#include "appender/AppenderPool.h"

#include <ctime>
#include <cstdarg>
#include <fstream>
#include <iostream>
//#include <memory>

using namespace std;

namespace siena {
namespace logger {

//! The maximum length of a printf-like log message.
//!
const int MAX_MESSAGE_LENGTH = 300;

/**
 * Constructs a LogEntry object.
 *
 * @param priority Priority of the log entry.
 * @param time Time the log entry was dispatched.
 * @param file File from where the log entry was dispatched.
 * @param line Line in which the log entry was dispatched.
 * @param function Function in which the log entry was dispatched.
 */
LogEntry::LogEntry(const Priority& priority, const tm& time, const string& file,
		const string& function, const int line) :
		lPriority(priority), //
		lTime(time), //
		lFile(file), //
		lFunction(function), //
		lLine(line), //
		lStream() {
}

/**
 * Submits the log message.
 */
LogEntry::~LogEntry() {
	AppenderPool::informAppenders(*this);
}

/**
 * Logs a simple message.
 *
 * @param msg Simple single string message.
 */
void LogEntry::message(const char* msg) {
	lStream << msg;
}

/**
 * Logs a printf-like message.
 *
 * @param format printf-like format string.
 * @param ... objects passed to printf.
 */
void LogEntry::formatMessage(const char* format, ...) {
	char msg[MAX_MESSAGE_LENGTH];

	// Format the ellipsis parameters
	va_list args;
	va_start(args, format);
	const int n = vsnprintf(msg, sizeof(msg), format, args);
	va_end(args);

	// Handle format error / to long messages
	if (n <= 0) {
		string emsg = "Format error in: ";
		emsg.append(format);
		message(emsg.c_str());
	} else if (n > MAX_MESSAGE_LENGTH) {
		string emsg = msg;
		emsg.append("... Exceeded log message length!");
		message(emsg.c_str());
	} else {
		message(msg);
	}
}

/**
 * Logs a streaming syntax message.
 *
 * @return Reference to the log stream, for streaming syntax logging.
 */
ostringstream& LogEntry::streamMessage() {
	return lStream;
}

} // namespace logger
} // namespace siena
