/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file ResultListener.cpp
 * \brief Implements the ResultListener class.
 *****************************************************************************/

#include "ResultListener.h"

using namespace std;

namespace siena {

ResultListener::~ResultListener() {
}

vector<ResultType>& ResultListener::needs() {
	static vector<ResultType> needs(0);
	return needs;
}

} // namespace siena
