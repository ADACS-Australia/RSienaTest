/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *****************************************************************************/
#ifndef SIMMELIAN_EFFECT_H_
#define SIMMELIAN_EFFECT_H_

#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "NetworkEffect.h"

namespace siena {

	class SimmelianEffect: public NetworkEffect {
		public:
			SimmelianEffect(const EffectInfo* pEffectInfo);
			~SimmelianEffect();

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
#endif // SIMMELIAN_EFFECT_H_

