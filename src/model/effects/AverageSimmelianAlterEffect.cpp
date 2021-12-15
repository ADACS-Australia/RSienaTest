/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageSimmelianAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageSimmelianAlterEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageSimmelianAlterEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "network/Simmelian.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
AverageSimmelianAlterEffect::AverageSimmelianAlterEffect(
		const EffectInfo * pEffectInfo, bool divide) :
	NetworkDependentBehaviorEffect(pEffectInfo), //
	pSimmelian(0)
{
	this->ldivide = divide;
	// Indicates whether there will be division by the outdegree of ego
}


AverageSimmelianAlterEffect::~AverageSimmelianAlterEffect() {
	clearSimmelian();
}

/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AverageSimmelianAlterEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	int simmelianOutdegree = 0;
	clearSimmelian();
	updateSimmelian((const OneModeNetwork*) this->pNetwork()); // simulated state
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the effect:
		// s_i(x) = v_i * avg(v_j) over all simmelian neighbors j of i.
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// This is d * avg(v_j), the average being taken over all neighbors
		// of i. This is what is calculated below.
		// if (not divide), instead of avg the total is used.
		
		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			if (pSimmelian->tieValue(actor, iter.actor()) > 0)
			{
				contribution += this->centeredValue(actor);
				simmelianOutdegree ++;
			}
		}			
		contribution *= difference;
		if ((this->ldivide) && (simmelianOutdegree >= 1))
		{
			contribution /= simmelianOutdegree;
		}
	}

	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AverageSimmelianAlterEffect::egoStatistic(int i, double * currentValues)
{
	double statistic = 0;
	clearSimmelian();
	updateSimmelian((const OneModeNetwork*) this->pNetwork()); // simulated state
	const Network * pNetwork = this->pNetwork();
	int simmelianOutdegree = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(i);
		 iter.valid();
		 iter.next())
	{
		if (pSimmelian->tieValue(i, iter.actor()) > 0)
			{
				statistic += currentValues[iter.actor()];
				simmelianOutdegree ++;
			}		
	}

	if (simmelianOutdegree > 0)
	{
		statistic *= currentValues[i];
		if (this->ldivide)
		{
			statistic /= simmelianOutdegree;
		}
	}

	return statistic;
}

/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AverageSimmelianAlterEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;
	int simmelianOutdegree = 0;

	clearSimmelian();
	updateSimmelian((const OneModeNetwork*) this->pNetwork()); // simulated state
	const Network * pNetwork = this->pNetwork();

	if (difference[ego] > 0)
	{
		if (pNetwork->outDegree(ego) > 0)
		{
			double thisStatistic = 0;
			double previousStatistic = 0;

			for (IncidentTieIterator iter = pNetwork->outTies(ego);
				 iter.valid();
				 iter.next())
			{
				double alterValue = currentValues[iter.actor()];
				double alterPreviousValue = currentValues[iter.actor()]
													+ difference[iter.actor()]; 
				if (pSimmelian->tieValue(ego, iter.actor()) > 0)
				{
					simmelianOutdegree ++;
					thisStatistic += alterValue;
					previousStatistic += alterPreviousValue;
				}
			}
			if (simmelianOutdegree >= 1)
			{
				thisStatistic *= currentValues[ego];
				previousStatistic *= (currentValues[ego] + difference[ego]);
				statistic = thisStatistic - previousStatistic;
				if (this->ldivide)
				{
					statistic /= simmelianOutdegree;
				}
			}
		}
	}

	return statistic;
}

// duplications with other simmelians - to be removed:
/**
 * Deletes pSimmelian.
 */
void AverageSimmelianAlterEffect::clearSimmelian() {
	if (pSimmelian != 0) {
		delete pSimmelian;
		pSimmelian = 0;
	}
}

/**
 * Updates pSimmelian with the simmelian strength on pNetwork.
 */
void AverageSimmelianAlterEffect::updateSimmelian(const OneModeNetwork* pNetwork) {
	pSimmelian = newSimmelian((const OneModeNetwork*) pNetwork);
}

}
