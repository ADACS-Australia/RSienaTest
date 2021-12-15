#include "FileAppender.h"

#include <fstream>

namespace siena {
namespace logger {

/**
 * \copybrief Appender::Appender()
 *
 * @param fName The file log entry are written to.
 * @param rFilter Log entry filter. Only entries passing the filter are
 *        formatted and appended.
 * @param rFormatter Formatter to be applied before appending.
 */
FileAppender::FileAppender(const std::string& fName, const Filter& rFilter,
		const Formatter& rFormatter) :
		Appender(rFilter, rFormatter), //
		rWriter(new std::ofstream(fName.c_str())) {
}

FileAppender::~FileAppender() {
	rWriter->flush();
	rWriter->close();
	delete rWriter;
}

/**
 * \copydoc Appender::writeLog()
 */
void FileAppender::writeLog(const std::string& msg) {
	(*rWriter) << msg << std::endl;
}

} // namespace logger
} // namespace siena

