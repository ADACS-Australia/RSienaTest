/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *****************************************************************************/
#include "SimmelianEffect.h"

#include "logger/Logger.h"
#include "data/NetworkLongitudinalData.h"
#include "model/effects/NetworkEffect.h"
#include "network/IncidentTieIterator.h"
#include "network/Simmelian.h"

/**
 * This is a very basic (naive/dumb) implementation of a
 * NetworkEffect where the network is filtered to only
 * consider simmelian ties. The simmelian network is reconstructed on each
 * call to preprocessEgo and statistic.
 * Constructed from CovariateSimmelianAlterEffect
 *
 * assumes: simmelian strength cannot be non-0 when summationTieNetwork is 0
 */
namespace siena
{

using namespace logger;

SimmelianEffect::SimmelianEffect(
		const EffectInfo* pEffectInfo):
	NetworkEffect(pEffectInfo), //
	pSimmelian(0)
{
}

SimmelianEffect::~SimmelianEffect() {
	clearSimmelian();
}

/**
 * Initializes contribution collection.
 *
 * Called once before all calls to calculateContribution.
 */
void SimmelianEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);
	clearSimmelian();
	updateSimmelian((const OneModeNetwork*) this->pNetwork()); // simulated state
	// This is now done for all simmelian effects;
	// must be inefficient. How to do better?
}

/**
 * Contribution of one tie flip.
 *
 * Called for every alter.
 */
double SimmelianEffect::calculateContribution(int alter) const {
	// only consider simmelian ties
	if (pSimmelian->tieValue(this->ego(), alter) > 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

/**
 * Initializes statistics collection.
 */
double SimmelianEffect::statistic(const Network * pSummationTieNetwork) {
	// update pSimmelian
	clearSimmelian();
	updateSimmelian((const OneModeNetwork*) pSummationTieNetwork);
	// default
	return NetworkEffect::statistic(pSummationTieNetwork);
}

/**
 * Contribution to the effect statistics.
 *
 * Called only for ties that exist in pSummationTieNetwork.
 */
double SimmelianEffect::tieStatistic(int alter) {
	if (pSimmelian->tieValue(this->ego(), alter) > 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

/**
 * Deletes pSimmelian.
 */
void SimmelianEffect::clearSimmelian() {
	if (pSimmelian != 0) {
		delete pSimmelian;
		pSimmelian = 0;
	}
}

/**
 * Updates pSimmelian with the simmelian strength on pNetwork.
 */
void SimmelianEffect::updateSimmelian(const OneModeNetwork* pNetwork) {
	pSimmelian = newSimmelian((const OneModeNetwork*) pNetwork);
}

}
