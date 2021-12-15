/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file OpenMPThreadFilter.h
 * \brief Defines and implements the OpenMPThreadFilter class.
 *****************************************************************************/

#ifndef RSIENA_LOGGER_OPENMP_THREAD_FILTER_H_
#define RSIENA_LOGGER_OPENMP_THREAD_FILTER_H_

#include "Filter.h"

#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

namespace siena {
namespace logger {

/**
 * Filters log entries by their OpenMP thread number.
 */
class OpenMPThreadFilter: public Filter {
public:

	/**
	 * Constructs a new OpenMPThreadFilter object.
	 *
	 * @param number Only entries dispatched in a thread with this number are
	 *        accepted.
	 */
	OpenMPThreadFilter(const int number) :
			lNumber(number) {
	}

	virtual ~OpenMPThreadFilter() {
	}

	/**
	 * \copydoc Filter::accept()
	 */
	virtual bool accept(const LogEntry& rEntry) const {
#ifdef SUPPORT_OPENMP
		return lNumber == omp_get_thread_num();
#else
		return true;
#endif
	}

private:
	OpenMPThreadFilter(const OpenMPThreadFilter&);
	OpenMPThreadFilter& operator=(const OpenMPThreadFilter&);

	const int lNumber;

};

} // namespace logger
} // namespace siena

#endif // RSIENA_LOGGER_OPENMP_THREAD_FILTER_H_
