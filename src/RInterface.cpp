/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file RInterface.cpp
 * \brief Implements RInterface.h.
 *****************************************************************************/

#include "estimator/SienaFit.h"
#include "RInterface.h"

#include <Eigen/Core>

#include "logger/Logger.h"
#include "logger/appender/AppenderPool.h"
#include "logger/appender/RConsoleAppender.h"
#include "logger/appender/FileAppender.h"
#include "logger/filter/AndFilter.h"
#include "logger/filter/PriorityRangeFilter.h"
#include "logger/filter/OpenMPThreadFilter.h"
#include "logger/formatter/MPIAwareFormatter.h"

#include "sim/StatisticsSimulation.h"
#include "sim/MetropolisHastingsSimulation.h"
#include "estimator/GMMEstimator.h"
#include "RUtil.h"

#include "data/Data.h"
#include "model/Model.h"

#include "sim/mpi/MPICommunicator.h"
#ifdef MPI2
#include "sim/mpi/MPIStatisticsSimulation.h"
#include <mpi.h>
#ifdef OPENMPI
#include <dlfcn.h>
#endif
#endif // MPI2

using namespace std;
using namespace siena;
using namespace siena::logger;
using namespace Eigen;

// The following two are defined here because I had some problem with the old
// logging setup on Mac OS X.
//! Filter family assigning each OpenMP thread its own file.
//!
vector<Filter*> logFileThreadFilters;
//! Final file filters, this filter for OpenMP thread and priority.
//!
vector<Filter*> logFileFilters;
//! Log file appenders.
//!
vector<Appender*> logFileAppenders;
MPIAwareFormatter* logFormatter;
Filter* logFilePriorityFilter;
Filter* logConsoleFilter;
Appender* logConsoleAppender;

/**
 * R-ify (construct a R list with reasonably named field, that is copy all
 * fields to R) of an SienaFit object.
 *
 * @param rFit SienaFit object.
 * @return R list containing the information in `rFit`.
 */
static SEXP rifySienaFit(const SienaFit& rFit) {
	int n = 21;
	SEXP sList;
	PROTECT(sList = allocVector(VECSXP, n));
	setAttrib(sList, R_NamesSymbol, allocVector(STRSXP, n));

	setNamedListElt(sList, --n, "n1", rifyIntVector(rFit.rPhase1NSimulations()));
	setNamedListElt(sList, --n, "n2", rifyIntVector(rFit.rPhase2NSimulations()));
	setNamedListElt(sList, --n, "n3", ScalarInteger(rFit.phase3NSimulations()));
	setNamedListElt(sList, --n, "fixed", rifyVectorXb(rFit.rFixed()));
	setNamedListElt(sList, --n, "theta", rifyVectorXd(rFit.rTheta()));
	setNamedListElt(sList, --n, "targets", rifyVectorXd(rFit.rTargets()));
	setNamedListElt(sList, --n, "period_wise_targets", rifyMatrixXd(rFit.rPeriodWiseTargets()));
	setNamedListElt(sList, --n, "tstat", rifyVectorXd(rFit.rTRatio()));
	setNamedListElt(sList, --n, "rate", rifyVectorXd(rFit.rRate()));
	setNamedListElt(sList, --n, "rateError", rifyVectorXd(rFit.rRateError()));
	setNamedListElt(sList, --n, "blockStructure", rifyMatrixXd(rFit.rBlockStructure()));
	setNamedListElt(sList, --n, "gamma", rifyMatrixXd(rFit.rGamma()));
	setNamedListElt(sList, --n, "phase1_gmmweight", rifyMatrixXd(rFit.rPhase1Weight()));
	setNamedListElt(sList, --n, "phase3_gmmweight", rifyMatrixXd(rFit.rPhase3Weight()));
	setNamedListElt(sList, --n, "phase1_derivative", rifyMatrixXd(rFit.rPhase1Derivative()));
	setNamedListElt(sList, --n, "phase3_derivative", rifyMatrixXd(rFit.rPhase3Derivative()));
	setNamedListElt(sList, --n, "phase3_statistics", rifyMatrixXd(rFit.rStatistics()));
	setNamedListElt(sList, --n, "phase3_period_wise_statistics", rifyVectorMatrixXd(rFit.rPeriodWiseStatistics()));
	setNamedListElt(sList, --n, "phase3_period_wise_scores", rifyVectorMatrixXd(rFit.rPeriodWiseScores()));
	setNamedListElt(sList, --n, "phase1_covariance", rifyMatrixXd(rFit.rPhase1Covariance()));
	setNamedListElt(sList, --n, "phase3_covariance", rifyMatrixXd(rFit.rPhase3Covariance()));
	assert(n == 0);

	UNPROTECT(1); // sList
	return sList;
}

extern "C" {

/**
 * Setups logging filters and appenders.
 *
 * @param sPriorityNameConsole Priority threshold for logging to console.
 * @param sPriorityNameFile Priority threshold for logging to files.
 * @param sBaseName Base name of log files. The final file name is
 *        'BASE_NAME-MPI_RANK:THREAD_ID.log'
 * @param sIncludeLocation True if the location (file, method, line) should be
 *        included in the output.
 * @param sNThreads Number of OpenMP threads. Each gets its own file.
 * @return R_NilValue
 */
SEXP sienaSetupLogger(SEXP sPriorityNameConsole, SEXP sPriorityNameFile,
		SEXP sBaseName, SEXP sIncludeLocation, SEXP sNThreads) {
	// clear previous logging setup
	sienaFinalizeLogger();
	// parse arguments
	const string priorityNameConsole(CHAR(STRING_ELT(sPriorityNameConsole, 0)));
	const string priorityNameFile(CHAR(STRING_ELT(sPriorityNameFile, 0)));
	const string baseName(CHAR(STRING_ELT(sBaseName, 0)));
	const bool includeLocation = asLogical(sIncludeLocation);
	int nThreads = asInteger(sNThreads);
	logFormatter = new MPIAwareFormatter(includeLocation);
	// console
	const Priority& priorityConsole = Priority::getByName(priorityNameConsole, Priority::DEBUG);
	logConsoleFilter = new PriorityRangeFilter(priorityConsole);
	logConsoleAppender = new RConsoleAppender(*logConsoleFilter, *logFormatter);
	AppenderPool::addAppender(logConsoleAppender);
	// files
	const Priority& priorityFile = Priority::getByName(priorityNameFile, Priority::DEBUG);
	logFilePriorityFilter = new PriorityRangeFilter(priorityFile);
	logFileThreadFilters = vector<Filter*>(nThreads);
	logFileFilters = vector<Filter*>(nThreads);
	logFileAppenders = vector<Appender*>(nThreads);
	for (int i = 0; i < nThreads; ++i) {
		logFileThreadFilters[i] = new OpenMPThreadFilter(i);
		logFileFilters[i] = new AndFilter(*logFileThreadFilters[i], *logFilePriorityFilter);
		ostringstream s;
		s << baseName << "-" << i << ".log";
		logFileAppenders[i] = new FileAppender(s.str(), *logFileFilters[i], *logFormatter);
		AppenderPool::addAppender(logFileAppenders[i]);
	}
	return R_NilValue;
}

/**
 * Terminate the file logging.
 *
 * @return R_NilValue
 */
SEXP sienaFinalizeLogger() {
	if (logConsoleFilter) {
		AppenderPool::removeAppender(logConsoleAppender);
		delete logConsoleAppender;
		delete logConsoleFilter;
		logConsoleAppender = 0;
		logConsoleFilter = 0;
	}
	for (vector<Appender*>::size_type i = 0; i < logFileAppenders.size(); ++i) {
		AppenderPool::removeAppender(logFileAppenders[i]);
		delete logFileAppenders[i];
		delete logFileFilters[i];
		delete logFileThreadFilters[i];
	}
	logFileThreadFilters.clear();
	logFileFilters.clear();
	logFileAppenders.clear();
	delete logFilePriorityFilter;
	logFilePriorityFilter = 0;
	delete logFormatter;
	logFormatter = 0;
	return R_NilValue;
}

/**
 * Wraps the logger to make it accessible from R.
 *
 * @param sPriority Name of the priority level. If no matching priority found
 *        defaults to 'FATAL'.
 * @param sEntry Log message content.
 * @return R_NilValue
 */
SEXP sienaLog(SEXP sPriority, SEXP sEntry) {
	const string priorityName(CHAR(STRING_ELT(sPriority, 0)));
	const Priority& priority = Priority::getByName(priorityName,
			Priority::FATAL);
	const char* entry(CHAR(STRING_ELT(sEntry, 0)));
	LOG(priority, entry);
	return R_NilValue;
}

/**
 * Initialize the libraries for use in parallelized environment.
 *
 * This has to be called before using MPI or OpenMP.
 *
 * @return R_NilValue
 */
SEXP sienaInitialize() {
	Eigen::initParallel();
	// Initialize MPI, if present
#ifdef MPI2
#ifdef OPENMPI
	// Don't load Open MPI to the private R symbol space.
	// http://www.open-mpi.de/faq/?category=troubleshooting#missing-symbols
	// http://www.open-mpi.org/community/lists/devel/2010/05/7931.php
#ifndef MAC
	if (!dlopen("libmpi.so.0", RTLD_GLOBAL | RTLD_LAZY)
			&& !dlopen("libmpi.so", RTLD_GLOBAL | RTLD_LAZY)) {
		Rprintf("%s\n", dlerror());
		return Rf_ScalarInteger(-1);
	}
#endif // MAC
#endif // OPENMPI
	int argc = 0;
	char** argv;
	MPI_Init(&argc, &argv);
#endif // MPI2
	return R_NilValue;
}

/**
 * Terminate MPI.
 *
 * @return R_NilValue
 */
SEXP sienaFinalize() {
	LOG(Priority::DEBUG, "MPI_Finalize");
#ifdef MPI2
	MPI_Finalize();
#endif
	sienaFinalizeLogger();
	return R_NilValue;
}

/**
 * @return Rank of the current node or 0 without MPI support.
 */
SEXP sienaMPIRank() {
	int rank = 0;
#ifdef MPI2
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	return ScalarInteger(rank);
}

/**
 * @return Number of nodes or 1 without MPI support.
 */
SEXP sienaMPISize() {
	int size = 1;
#ifdef MPI2
	MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
	return ScalarInteger(size);
}

/**
 * R sided MPI broadcasting of integers.
 *
 * @param sInteger Integer to broadcast. Only the root needs to supply a
 *        meaningful value.
 * @return The integer distributed by the root.
 */
SEXP sienaMPIBcastInt(SEXP sInteger) {
	int integer = asInteger(sInteger);
#ifdef MPI2
	MPI_Bcast(&integer, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
	return ScalarInteger(integer);
}

/**
 * Runs the estimator on each data set in the group.
 *
 * @param sSienaModel An R object (usually called `z`) of the class
 *        "sienaAlgorithm".
 * @param sModelPtr A pointer to the C++ model object.
 * @param sGroupDataPtr A pointer to a C++ vector holding pointers to data
 *        objects.
 * @param sNThreads Number of OpenMP threads to use.
 * @return A list with the results of each estimation.
 */
SEXP sienaEstimateGroup(SEXP sSienaModel, SEXP sModelPtr, SEXP sGroupDataPtr,
		SEXP sNThreads) {
	// Type check parameters and convert them to C
	vector<Data*>& pGroupData = *static_cast<vector<Data*>*>(R_ExternalPtrAddr(
			sGroupDataPtr));

	vector<Data*>::size_type n = pGroupData.size();
	SEXP sList;
	PROTECT(sList = allocVector(VECSXP, n));
	// Go through the data sets and append results to the list.
	while (n > 0) {
		--n;
		SET_VECTOR_ELT(sList, n,
				sienaEstimate(sSienaModel, sModelPtr,
						R_MakeExternalPtr(pGroupData[n], R_NilValue,
								R_NilValue), sNThreads));
	}
	UNPROTECT(1); // sList
	return sList;
}

/**
 * Run the estimator.
 *
 * @param sSienaModel An R object (usually called `z`) of the class
 *        "sienaAlgorithm".
 * @param sModelPtr A pointer to the C++ model object.
 * @param sDataPtr A pointer to the C++ data object.
 * @param sNThreads Number of OpenMP threads to use.
 * @return An R lists representation of the `SienaFit` object returned by the
 *         estimator.
 */
SEXP sienaEstimate(SEXP sSienaModel, SEXP sModelPtr, SEXP sDataPtr,
		SEXP sNThreads) {
	// Type check parameters and convert them to C
	Model* pModel = static_cast<Model*>(R_ExternalPtrAddr(sModelPtr));
	Data* pData = static_cast<Data*>(R_ExternalPtrAddr(sDataPtr));
	SEXP x = getNamedListElt(sSienaModel, "x");
	const bool maxLike = asLogical(getNamedListElt(x, "maxlike"));
	const bool finDiff = asLogical(
			getNamedListElt(x, "FinDiff.method"));
	const bool dolby = asLogical(getNamedListElt(x, "dolby"));
#ifdef SUPPORT_OPENMP
	const int nThreads = asInteger(sNThreads);
#else
	const int nThreads = 1;
#endif
	const int nsub = asInteger(getNamedListElt(x, "nsub"));
	const int n1min = asInteger(getNamedListElt(x, "n1min"));
	const int n3 = asInteger(getNamedListElt(x, "n3"));
	const double firstg = asReal(getNamedListElt(x, "firstg"));
	const double reduceg = asReal(getNamedListElt(x, "reduceg"));

#ifdef MPI2
	if (!MPICommunicator::theCommunicator().isRoot()) {
		// Using MPI and the current node is not the root. Start the simulation in
		// slave mode.
		MPIStatisticsSimulation sim(pModel, pData, nThreads);
		sim.slaveLoop();
		return R_NilValue;
	}
#endif
	// Choose simulation and differentiation method.
	Simulation* pSim = 0;
	DifferentiationType diff = finDiff ? FINITE_DIFFERENCE : SCORE_DEVIATION;
	if (maxLike) {
		pSim = new MetropolisHastingsSimulation(pModel, pData);
		diff = MAXIMUM_LIKELIHOOD;
		}
#ifdef MPI2
	else if (MPICommunicator::theCommunicator().nNodes() > 1) {
		// Having more than one node needs the MPI aware simulation.
		pSim = new MPIStatisticsSimulation(pModel, pData, nThreads);
	}
#endif
	else {
		// Single core uses the normal simulation.
		pSim = new StatisticsSimulation(pModel, pData, nThreads);
	}
	GMMEstimator con(*pSim, diff, dolby, nsub, firstg, reduceg, n1min, n3);
	bool haveDfra = asLogical(getNamedListElt(sSienaModel, "haveDfra"));
	if (haveDfra) {
		LOGS(Priority::INFO) << "estimate from previous answer";
		SienaFit s;
		s.fixed(asVectorXb(getNamedListElt(sSienaModel, "fixed")));
		s.phase3Weight(asMatrixXd(getNamedListElt(sSienaModel, "gmmweight")));
		s.phase3Derivative(asMatrixXd(getNamedListElt(sSienaModel, "dfra")));
		s.statistics(asMatrixXd(getNamedListElt(sSienaModel, "sf")));
		con.estimateFromPrevAns(s);
	} else {
		con.estimate();
	}
	const SEXP ret = rifySienaFit(con.rFit());
	delete pSim;
	return ret;
}

/**
 * Runs the MPIStatisticsSimulation in slave mode.
 *
 * This can be used the setup separate R scripts for master and slave.
 * (At the moment this does not bring any benefit since the serialization is
 * in R.)
 *
 * @param sModelPtr A pointer to the C++ model object.
 * @param sDataPtr A pointer to the C++ data object.
 * @param sNThreads Number of OpenMP threads to use.
 * @return R_NilValue
 */
SEXP sienaSlaveSimulation(SEXP sModelPtr, SEXP sDataPtr, SEXP sNThreads) {
#ifdef MPI2
	Model* pModel = static_cast<Model*>(R_ExternalPtrAddr(sModelPtr));
	Data* pData = static_cast<Data*>(R_ExternalPtrAddr(sDataPtr));
	const int nThreads = asInteger(sNThreads);
	MPIStatisticsSimulation sim(pModel, pData, nThreads);
	sim.slaveLoop();
#endif
	return R_NilValue;
}

} // extern "C"
