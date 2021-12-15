#ifndef STOPWATCH_H_
#define STOPWATCH_H_

#include <string>

//#include <chrono>
#include <ctime>

//typedef std::chrono::high_resolution_clock Clock;
//typedef std::chrono::duration<long, std::milli> milliseconds;

/**
 * Simple stop watch for debugging.
 */
class StopWatch {
public:
	StopWatch(const std::string& rName);
	~StopWatch();

	StopWatch(const StopWatch& r);
	StopWatch& operator=(const StopWatch& r);

	void start();
	void stop();

private:
	std::string lName;

//	mutable Clock::time_point lStart;
	struct timespec lStart;
};

#endif // STOPWATCH_H_
