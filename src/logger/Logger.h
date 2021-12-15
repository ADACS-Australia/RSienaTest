#ifndef LOGGER_H_
#define LOGGER_H_

#include "Priority.h"
#include "LogEntry.h"

namespace siena {
namespace logger {

inline tm getLocalTime() {
	time_t t;
	time(&t);
	return *localtime(&t);
}

} // namespace logger
} // namespace siena

#define LOG(priority, msg) \
siena::logger::LogEntry(priority, siena::logger::getLocalTime(), \
__FILE__, __FUNCTION__, __LINE__).message(msg);

#define LOGF(priority, format, ...) \
siena::logger::LogEntry(priority, siena::logger::getLocalTime(), \
__FILE__, __FUNCTION__, __LINE__).formatMessage(format, __VA_ARGS__);

#define LOGS(priority) \
siena::logger::LogEntry(priority, siena::logger::getLocalTime(), \
__FILE__, __FUNCTION__, __LINE__).streamMessage()

#endif /* LOGGER_H_ */
