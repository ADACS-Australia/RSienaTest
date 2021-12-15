/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MPIAwareFormatter.cpp
 * \brief Implements the MPIAwareFormatter class.
 *****************************************************************************/

#include "MPIAwareFormatter.h"

#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif
#ifdef MPI2
#include <mpi.h>
#endif
#include <iomanip>

using namespace std;

namespace siena {
namespace logger {

//! Time stamp format.
//!
const char* TIME_FORMAT = "%Y-%m-%d %H:%M:%S";

/**
 * Constructs a new MPIAwareFormatter object.
 *
 * @param location Whether to include the code location (file, function, line)
 *        in the formatted message.
 */
MPIAwareFormatter::MPIAwareFormatter(const bool location) :
		Formatter(), //
		lLocation(location), //
		lRank(0) {
#ifdef MPI2
	MPI_Comm_rank(MPI_COMM_WORLD, &lRank);
#endif
}

MPIAwareFormatter::~MPIAwareFormatter() {
}

/**
 * \copydoc Formatter::format()
 */
string MPIAwareFormatter::format(const LogEntry& rEntry) const {
	strftime(lStamp, sizeof(lStamp), TIME_FORMAT, &rEntry.time());
	ostringstream entry;
	entry << lRank;
#ifdef SUPPORT_OPENMP
	entry << ":" << omp_get_thread_num();
#endif
	entry << " " << setw(7) << rEntry.priority().toString() << setw(0)
		<< " [" << lStamp << "] ";
	if (lLocation) {
		entry << rEntry.file() << "." << rEntry.function() << ":"
			<< rEntry.line();
	}
	entry << " " << rEntry.message();
	return entry.str();
}

} // namespace logger
} // namespace siena
