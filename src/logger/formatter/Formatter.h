#ifndef RSIENA_LOGGER_FORMATTER_H_
#define RSIENA_LOGGER_FORMATTER_H_

#include "logger/LogEntry.h"

#include <string>

namespace siena {
namespace logger {

/**
 * Log formatter building strings out of log entries.
 */
class Formatter {
public:
	virtual ~Formatter() {
	}

	/**
	 * @param rEntry The log entry.
	 * @return Formatted message.
	 */
	virtual std::string format(const LogEntry& rEntry) const = 0;

protected:
	Formatter() {
	}

private:
	Formatter(const Formatter&);
	Formatter& operator=(const Formatter&);

};

} // namespace logger
} // namespace siena

#endif // RSIENA_LOGGER_FORMATTER_H_
