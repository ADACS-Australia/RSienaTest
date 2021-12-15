/*
 * AppenderPool.h
 *
 *  Created on: 04.09.2013
 *      Author: ortmann
 */

#ifndef APPENDERPOOL_H_
#define APPENDERPOOL_H_

namespace siena {
namespace logger {

class Appender;
class LogEntry;

/**
 * Global static pool of appenders.
 */
class AppenderPool {
public:
	static void addAppender(Appender* appender);
	static void removeAppender(Appender* appender);
	static void informAppenders(const LogEntry& rEntry);

private:
	AppenderPool();

};

} // namespace logger
} // namespace siena

#endif /* APPENDERPOOL_H_ */
