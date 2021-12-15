#include "utils/debugging/StopWatch.h"

#include "logger/Logger.h"

using namespace std;
using namespace siena;
using namespace siena::logger;

#ifndef __unix__ // __APPLE__ and _WIN32
// clock_gettime for Mac
#include <sys/time.h>
#define CLOCK_MONOTONIC 0
int clock_gettime(int, struct timespec* t) {
	struct timeval now;
	int r = gettimeofday(&now, 0);
	if (r)
		return r;
	t->tv_sec = now.tv_sec;
	t->tv_nsec = now.tv_usec * 1e3;
	return 0;
}
#endif

/**
 * Constructs a new StopWatch object.
 *
 * @param rName Name to log before the time.
 */
StopWatch::StopWatch(const std::string& rName) :
		lName(rName) {
}

StopWatch::~StopWatch() {
}

/**
 * Copy constructs a StopWatch with the same name.
 *
 * @param r StopWatch to duplicate.
 */
StopWatch::StopWatch(const StopWatch& r) :
		lName(r.lName) {
}

/**
 * Copy assigns the name of a StopWatch.
 *
 * @param r StopWatch to copy the name from.
 * @return *this.
 */
StopWatch& StopWatch::operator=(const StopWatch& r) {
	lName = r.lName;
	return *this;
}

/**
 * Start taking the time.
 */
void StopWatch::start() {
//	lStart = Clock::now();
	clock_gettime(CLOCK_MONOTONIC, &lStart);
}

/**
 * Stop the time and log the duration since the last start as a debug message.
 */
void StopWatch::stop() {
//	LOGS(siena::logger::Priority::DEBUG)<<lName<<": "
//	<<std::chrono::duration_cast < milliseconds > (Clock::now() - lStart).count()
//	<<" ms";
	struct timespec end;
	clock_gettime(CLOCK_MONOTONIC, &end);
	long duration = ((end.tv_sec * 1e3) + (end.tv_nsec / 1e6))
			- ((lStart.tv_sec * 1e3) + (lStart.tv_nsec / 1e6));
	LOGS(siena::logger::Priority::DEBUG)<<lName<<": "<<duration<<" ms";
}
