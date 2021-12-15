/*
 * FileAppender.h
 *
 *  Created on: 04.09.2013
 *      Author: ortmann
 */

#ifndef FILEAPPENDER_H_
#define FILEAPPENDER_H_

#include "Appender.h"

//#include <fstream> // Conflicts with R
#include <iosfwd> // Forward declarations of ofstream

namespace siena {
namespace logger {

/**
 * Log appender writing to a file.
 */
class FileAppender: public Appender {
public:
	FileAppender(const std::string& fName, const Filter& rFilter,
			const Formatter& rFormatter);

	virtual ~FileAppender();

protected:
	void writeLog(const std::string& msg);

private:
	std::ofstream * const rWriter;

};

} // namespace logger
} // namespace siena

#endif /* FILEAPPENDER_H_ */
