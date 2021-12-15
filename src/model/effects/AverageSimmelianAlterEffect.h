/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageSimmelianAlterEffect.h
 *
 * Description: This file contains the definition of the
 * AverageSimmelianAlterEffect class.
 *****************************************************************************/

#ifndef AVERAGESIMMELIANALTEREFFECT_H_
#define AVERAGESIMMELIANALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"
#include "network/OneModeNetwork.h"

namespace siena
{

/**
 * Average alter effect defined as the product of the ego with the average
 * of its neighbors (with respect to a certain network).
 */
class AverageSimmelianAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	AverageSimmelianAlterEffect(const EffectInfo * pEffectInfo, bool divide);
	~AverageSimmelianAlterEffect();

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoEndowmentStatistic(int ego, const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);

protected:
	double statistic(const Network * pSummationTieNetwork);

private:
	// ldivide indicates whether there will be division by the outdegree
	bool ldivide;
	const OneModeNetwork* pSimmelian;
	void clearSimmelian();
	void updateSimmelian(const OneModeNetwork* pNetwork);
};

}

#endif /*AVERAGESIMMELIANALTEREFFECT_H_*/
