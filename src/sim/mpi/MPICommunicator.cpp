/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MPICommunicator.cpp
 * Implementation of MPICommunicator.h
 *****************************************************************************/

#include "MPICommunicator.h"

using namespace std;

namespace siena {

/**
 * Returns the singleton communicator.
 *
 * @return Communicator.
 */
MPICommunicator& MPICommunicator::theCommunicator() {
	static MPICommunicator com;
	return com;
}

/**
 * Constructs the MPICommunicator object.
 */
MPICommunicator::MPICommunicator() :
#ifdef MPI2
		lrComm(MPI_COMM_WORLD), //
#endif // MPI2
		lRootRank(0), //
		lRank(0), //
		lCommSize(1) {
#ifdef MPI2
	MPI_Comm_rank(lrComm, &lRank);
	MPI_Comm_size(lrComm, &lCommSize);
#endif // MPI2
}

/**
 * @return Rank (ID) of the current process.
 */
int MPICommunicator::rank() const {
	return lRank;
}

/**
 * @return True if the current process is the root process.
 */
bool MPICommunicator::isRoot() const {
	return lRank == lRootRank;
}

/**
 * @return Number of processes.
 */
int MPICommunicator::nNodes() const {
	return lCommSize;
}

///////////////////////////////////////////////////////////////////////////////
#ifdef MPI2
/**
 * Broadcasts an integer from the root process to all other.
 *
 * @param pInt Pointer to the int to broadcast.
 */
bool MPICommunicator::broadcastInt(int* pInt) const {
	return MPI_SUCCESS == MPI_Bcast(pInt, 1, MPI_INT, lRootRank, lrComm);
}

/**
 * Broadcasts a Vector from the root process to all other.
 *
 * @param pVector Pointer to the Vector to broadcast.
 */
bool MPICommunicator::broadcastVectorXd(Eigen::VectorXd* pVector) const {
	return MPI_SUCCESS
			== MPI_Bcast(pVector->data(), pVector->size(), MPI_DOUBLE,
					lRootRank, lrComm);
}

/**
 * Collects double vectors from all processes and puts them in one double
 * vector in the root process.
 *
 * @param[in] pSend Pointer to the double vector to send.
 * @param[out] pRecv Pointer to the double vector the received values are
 * written to. This should be nullptr on slave node and have a size equal to
 * nNodes() times size of the pSend vector on the root process.
 */
bool MPICommunicator::gatherDouble(const vector<double>* pSend,
		vector<double>* pRecv) const {
	assert(pSend != 0/*ptr*/);              // All nodes contribute send
	assert(!isRoot() || (pRecv != 0/*ptr*/       // Either we are not root...
		&& lCommSize * pSend->size() == pRecv->size())); // ...or have a matching recv
	// We actually need a pointer to the data
	double* recv_data = pRecv ? &(*pRecv)[0] : 0/*ptr*/;
	// Casting away constness of the send buffer should not be a problem
	return MPI_SUCCESS
			== MPI_Gather(const_cast<double*>(&(*pSend)[0]), pSend->size(),
			MPI_DOUBLE, recv_data, pSend->size(), MPI_DOUBLE, lRootRank, lrComm);
}
#endif // MPI2
///////////////////////////////////////////////////////////////////////////////

}// namespace siena
