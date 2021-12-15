/*
 * Priority.h
 *
 *  Created on: 04.09.2013
 *      Author: ortmann
 */

#ifndef PRIORITY_H_
#define PRIORITY_H_

#include <string>

namespace siena {
namespace logger {

/**
 * Log levels.
 */
class Priority {
public:
	//! Priority for logging fatal errors.
	//!
	static const Priority FATAL;
	//! Priority for logging errors.
	//!
	static const Priority ERROR;
	//! Priority for logging warnings.
	//!
	static const Priority WARNING;
	//! Priority for normal information messages.
	//!
	static const Priority INFO;
	//! Priority for verbose messages.
	//!
	static const Priority VERBOSE;
	//! Priority for debugging output.
	//!
	static const Priority DEBUG;

	static const Priority& getByName(const std::string& rName,
			const Priority& rDefault);

	/**
	 * Tests whether the current priority is lower.
	 *
	 * @param rhs Second priority.
	 * @return `true` if the priority is lower (level higher) than the `rhs`
	 *         Priority.
	 */
	bool operator<=(const Priority& rhs) const {
		return rLevel >= rhs.rLevel;
	}

	/**
	 * @return Human readable priority name.
	 */
	const std::string& toString() const {
		return rName;
	}

private:
	Priority(int level, std::string name) :
			rLevel(level), //
			rName(name) {
	}

	Priority(const Priority&);
	Priority& operator=(const Priority&);

	const int rLevel;
	const std::string rName;
};

} // namespace logger
} // namespace siena

#endif // PRIORITY_H_
