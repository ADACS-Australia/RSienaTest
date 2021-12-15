/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ResultModificator.cpp
 * \brief Implements the ResultModificator class.
 *****************************************************************************/

#include "ResultModificator.h"

using namespace std;

namespace siena {

ResultModificator::~ResultModificator() {
}

vector<ResultType>& ResultModificator::needs() {
	static vector<ResultType> needs(0);
	return needs;
}

} // namespace siena
