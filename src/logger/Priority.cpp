/*
 * Priority.cpp
 *
 *  Created on: 04.09.2013
 *      Author: ortmann
 */

#include "Priority.h"

namespace siena {
namespace logger {

const Priority Priority::FATAL = Priority(0, "FATAL");
const Priority Priority::ERROR = Priority(1, "ERROR");
const Priority Priority::WARNING = Priority(2, "WARNING");
const Priority Priority::INFO = Priority(3, "INFO");
const Priority Priority::VERBOSE = Priority(4, "VERBOSE");
const Priority Priority::DEBUG = Priority(5, "DEBUG");

/**
 * Get the priority to give name.
 *
 * @param rName Priority name.
 * @param rDefault Priority to return if none matches the name.
 * @return The priority named `rName` or `rDefault` if no priority matching
 *         `rName` was found.
 */
const Priority& Priority::getByName(const std::string& rName,
		const Priority& rDefault) {
	if (rName == Priority::DEBUG.toString())
		return Priority::DEBUG;
	if (rName == Priority::VERBOSE.toString())
		return Priority::VERBOSE;
	if (rName == Priority::INFO.toString())
		return Priority::INFO;
	if (rName == Priority::WARNING.toString())
		return Priority::WARNING;
	if (rName == Priority::ERROR.toString())
		return Priority::ERROR;
	if (rName == Priority::FATAL.toString())
		return Priority::FATAL;
	return rDefault;
}

} // namespace logger
} // namespace siena
