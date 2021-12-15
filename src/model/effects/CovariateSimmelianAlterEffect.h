/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *****************************************************************************/
#ifndef COVAR_SIMMELIAN_ALTER_EFFECT_H_
#define COVAR_SIMMELIAN_ALTER_EFFECT_H_

#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "CovariateDependentNetworkEffect.h"

namespace siena {

	class CovariateSimmelianAlterEffect: public CovariateDependentNetworkEffect {
		public:
			CovariateSimmelianAlterEffect(const EffectInfo* pEffectInfo);
			~CovariateSimmelianAlterEffect();

			void preprocessEgo(int ego);
			double calculateContribution(int) const;

		protected:
			double statistic(const Network * pSummationTieNetwork);
			double tieStatistic(int alter);

		private:
			const OneModeNetwork* pSimmelian;
			void clearSimmelian();
			void updateSimmelian(const OneModeNetwork* pNetwork);
	};

} // namespace siena
#endif // COVAR_SIMMELIAN_ALTER_EFFECT_H_

