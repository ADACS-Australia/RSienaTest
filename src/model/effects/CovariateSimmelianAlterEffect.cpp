/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * This is a very basic (naive/dumb) implementation of a
 * CovariateDependentNetworkEffect where the network is filtered to only
 * consider simmelian ties. The simmelian network is reconstructed on each
 * call to preprocessEgo and statistic.
 *
 * assumes: simmelian strength cannot be non-0 when summationTieNetwork is 0
 *****************************************************************************/

#include "CovariateSimmelianAlterEffect.h"

#include "logger/Logger.h"
#include "data/NetworkLongitudinalData.h"
#include "model/effects/NetworkEffect.h"
#include "network/IncidentTieIterator.h"
#include "network/Simmelian.h"

namespace siena
{

using namespace logger;

CovariateSimmelianAlterEffect::CovariateSimmelianAlterEffect(
		const EffectInfo* pEffectInfo):
	CovariateDependentNetworkEffect(pEffectInfo), //
	pSimmelian(0)
{
}

CovariateSimmelianAlterEffect::~CovariateSimmelianAlterEffect() {
	clearSimmelian();
}

/**
 * Initializes contribution collection.
 *
 * Called once before all calls to calculateContribution.
 */
void CovariateSimmelianAlterEffect::preprocessEgo(int ego)
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
double CovariateSimmelianAlterEffect::calculateContribution(int alter) const {
	// only consider simmelian ties
	if (pSimmelian->tieValue(this->ego(), alter) > 0)
		return this->value(alter);
	return 0;
}

/**
 * Initializes statistics collection.
 */
double CovariateSimmelianAlterEffect::statistic(const Network * pSummationTieNetwork) {
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
double CovariateSimmelianAlterEffect::tieStatistic(int alter) {
	if (!this->missing(alter) && pSimmelian->tieValue(this->ego(), alter) > 0)
		return this->value(alter);
	return 0;
}

/**
 * Deletes pSimmelian.
 */
void CovariateSimmelianAlterEffect::clearSimmelian() {
	if (pSimmelian != 0) {
		delete pSimmelian;
		pSimmelian = 0;
	}
}

/**
 * Updates pSimmelian with the simmelian strength on pNetwork.
 */
void CovariateSimmelianAlterEffect::updateSimmelian(const OneModeNetwork* pNetwork) {
	pSimmelian = newSimmelian((const OneModeNetwork*) pNetwork);
}

}
