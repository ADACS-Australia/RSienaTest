/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MPICommunicator.h
 *
 * \class siena::MPICommunicator
 * \brief The MPICommunicator provides some higher level methods to send/receive
 *        objects.
 *****************************************************************************/

#ifndef RSIENA_MPI_COMMUNICATOR_H_
#define RSIENA_MPI_COMMUNICATOR_H_

#include <Eigen/Core>
#include <vector>

#ifdef MPI2
#include <mpi.h>
#endif // MPI2

namespace siena {

class MPICommunicator {
public:
	static MPICommunicator& theCommunicator();

	int rank() const;
	bool isRoot() const;
	int nNodes() const;

#ifdef MPI2
	bool broadcastInt(int* pInt) const;
	bool broadcastVectorXd(Eigen::VectorXd* pVector) const;
	bool gatherDouble(const std::vector<double>* pSend,
			std::vector<double>* pRecv = 0/*ptr*/) const;
#endif // MPI2

private:
	MPICommunicator();
	MPICommunicator(const MPICommunicator&);
	MPICommunicator& operator=(const MPICommunicator&);

#ifdef MPI2
	//! The MPI Comm (that is a group of nodes).
	//! Usually MPI_COMM_WORLD.
	const MPI_Comm lrComm;
#endif // MPI2

	//! The rank of the root node.
	//!
	const int lRootRank;
	//! The rank of the current node.
	//!
	int lRank;
	//! Number of nodes in the lrComm.
	//!
	int lCommSize;

};

} // namespace siena

#endif // RSIENA_MPI_COMMUNICATOR_H_
