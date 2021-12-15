/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file RInterface.h
 * \brief Provides an R interface to the C objects via external pointers.
 *****************************************************************************/

#ifndef RSIENA_RINTERFACE_H_
#define RSIENA_RINTERFACE_H_

#include <Rinternals.h>

extern "C" {

SEXP sienaSetupLogger(SEXP sPriorityNameConsole, SEXP sPriorityNameFile,
		SEXP sBaseName, SEXP sIncludeLocatiologIncludeLocationn, SEXP sNThreads);
SEXP sienaFinalizeLogger();
SEXP sienaLog(SEXP sPriority, SEXP sEntry);

SEXP sienaInitialize();
SEXP sienaFinalize();

SEXP sienaMPIRank();
SEXP sienaMPISize();
SEXP sienaMPIBcastInt(SEXP sInteger);

SEXP sienaEstimateGroup(SEXP sSienaModel, SEXP sModelPtr, SEXP sGroupDataPtr,
		SEXP sNThreads);
SEXP sienaEstimate(SEXP sSienaModel, SEXP sModelPtr, SEXP sDataPtr, SEXP
		sNThreads);
SEXP sienaSlaveSimulation(SEXP sModelPtr, SEXP sDataPtr, SEXP sNThreads);

} // extern "C"

#endif // RSIENA_RINTERFACE_H_
