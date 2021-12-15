/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file MPIStatisticsSimulation.h
 * \brief Defines the MPIStatisticsSimulation class.
 *****************************************************************************/

#ifdef MPI2

#ifndef RSIENA_MPI_STATISTICS_SIMULATION_H_
#define RSIENA_MPI_STATISTICS_SIMULATION_H_

#include "sim/StatisticsSimulation.h"
#include "sim/mpi/MPICommunicator.h"

namespace siena {

/**
 * Wrapper around the EpochSimulation, providing parallelization using a
 * combination of OpenMP and MPI.
 */
class MPIStatisticsSimulation: public StatisticsSimulation {
public:
	MPIStatisticsSimulation(Model* pModel, Data* pData,
			const unsigned int nThreads);
	virtual ~MPIStatisticsSimulation();

	virtual int nSimulations();

	void slaveLoop();

//	virtual void saveSeeds(const bool save);
	virtual void simulate(bool withSeeds = false);

	virtual void updateParameters(const Eigen::VectorXd& rTheta);

protected:
	//! double vector for fixed layout matrices.
	//!
	std::vector<double> lMPIStatisticsData;
	//! Vector of simulated statistics, collected from all processes.
	//!
	std::vector<Eigen::Map<Eigen::MatrixXd> > lMPIStatistics;

	//! double vector for fixed layout matrices.
	//!
	std::vector<double> lMPIScoresData;
	//! Vector of simulated scores, collected from all processes.
	//!
	std::vector<Eigen::Map<Eigen::MatrixXd> > lMPIScores;

	//! double vector for fixed layout matrices.
	//!
	std::vector<double> lMPITimesData;
	//! Vector of times needed for conditional simulation.
	//!
	std::vector<Eigen::Map<Eigen::VectorXd> > lMPITimes;

	virtual void needsChanged();

private:
	MPIStatisticsSimulation(const MPIStatisticsSimulation&);
	MPIStatisticsSimulation& operator=(const MPIStatisticsSimulation&);

};

} // namespace siena

#endif // RSIENA_MPI_STATISTICS_SIMULATION_H_
#endif // MPI2
