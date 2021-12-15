/*
 * AppenderPool.cpp
 *
 *  Created on: 04.09.2013
 *      Author: ortmann
 */

#include "AppenderPool.h"

#include "Appender.h"

#include <vector>
#include <algorithm> // find

namespace siena {
namespace logger {

class LogEntry;

typedef std::vector<Appender*>::iterator Iter;
static std::vector<Appender*> rAppenders;

/**
 * @param appender The appender to add.
 */
void AppenderPool::addAppender(Appender* appender) {
	if (std::find(rAppenders.begin(), rAppenders.end(), appender)
			== rAppenders.end()) {
		rAppenders.push_back(appender);
	}
}

/**
 * @param appender The appender to remove.
 */
void AppenderPool::removeAppender(Appender* appender) {
	rAppenders.erase(std::find(rAppenders.begin(), rAppenders.end(), appender));
}

/**
 * Sends an log entry to all appenders in the pool.
 *
 * @param rEntry The log entry.
 */
void AppenderPool::informAppenders(const LogEntry& rEntry) {
	for (Iter iter = rAppenders.begin(); iter != rAppenders.end(); ++iter) {
		(*iter)->appendEntry(rEntry);
	}
}

} // namespace logger
} // namespace siena
