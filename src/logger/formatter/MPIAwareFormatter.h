/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MPIAwareFormatter.h
 * \brief Defines the MPIAwareFormatter class.
 *****************************************************************************/

#ifndef MPIAWAREFORMATTER_H_
#define MPIAWAREFORMATTER_H_

#include "Formatter.h"

namespace siena {
namespace logger {

/**
 * Log formatter prepending the MPI rank and thread number.
 */
class MPIAwareFormatter: public Formatter {
public:
	explicit MPIAwareFormatter(const bool location);
	virtual ~MPIAwareFormatter();

	virtual std::string format(const LogEntry& rEntry) const;

private:
	MPIAwareFormatter(const MPIAwareFormatter&);
	MPIAwareFormatter& operator=(const MPIAwareFormatter&);

	const bool lLocation;
	int lRank;
	mutable char lStamp[30];

};

} // namespace logger
} // namespace siena

#endif // MPIAWAREFORMATTER_H_
