#include "ApproxSize.h"

#include "logger/Logger.h"
#include "data/ActorSet.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingDyadicCovariate.h"
#include "network/Network.h"

#include <vector>

using namespace std;
using namespace siena;
using namespace siena::logger;

int approxActorSet(const ActorSet& r) {
	int size = sizeof(r);
	LOGF(Priority::DEBUG, "ActorSet\t%s\t%d b", r.name().c_str(), size);
	return size;
}
int approxNetwork(const Network& r) {
	int size = sizeof(r);
	// lpOutTies
	for (int n = 0; n < r.n(); ++n) {
		size += r.inDegree(n) * 2 * sizeof(int); // each map entry is 2 ints
	}
	// lpInTies
	for (int m = 0; m < r.m(); ++m) {
		size += r.inDegree(m) * 2 * sizeof(int); // each map entry is 2 ints
	}
	LOGF(Priority::DEBUG, "Network\t\t%d b", size);
	return size;
}
int approxNetworkLongitudinalData(const NetworkLongitudinalData& r) {
	int size = sizeof(r);
	// networks
	for (int i = 0; i < r.observationCount(); ++i) {
		size += approxNetwork(*r.pNetwork(i));
		size += approxNetwork(*r.pStructuralTieNetwork(i));
		size += approxNetwork(*r.pMissingTieNetwork(i));
		size += approxNetwork(*r.pNetworkLessMissing(i));
		size += approxNetwork(*r.pNetworkLessMissingStart(i));
	}
	// ldensity
	size += r.observationCount() * sizeof(double);
	LOGF(Priority::DEBUG, "NetworkLongitudinalData\t%s\t%d b", r.name().c_str(),
			size);
	return size;
}
int approxBehaviorLongitudinalData(const BehaviorLongitudinalData& r) {
	int size = sizeof(r);
	// lvalues, lvaluesLessMissings, lvaluesLessMissingStarts
	size += 3 * r.observationCount() * r.n() * sizeof(int);
	// lmissing, lstructural
	size += 2 * r.observationCount() * r.n() * sizeof(bool);
	// lobservedDistributions
	LOGF(Priority::DEBUG, "BehaviorLongitudinalData\t%s\t%d b",
			r.name().c_str(), size);
	return size;
}
int approxLongitudinalData(const LongitudinalData& r) {
	int size = sizeof(r);
	// lupOnly, ldownOnly
	size += 2 * (r.observationCount() - 1) * sizeof(bool);

	if (const NetworkLongitudinalData* a =
			dynamic_cast<const NetworkLongitudinalData*>(&r)) {
		size += approxNetworkLongitudinalData(*a);
	} else if (const BehaviorLongitudinalData* b =
			dynamic_cast<const BehaviorLongitudinalData*>(&r)) {
		size += approxBehaviorLongitudinalData(*b);
	}
	LOGF(Priority::DEBUG, "LongitudinalData\t%s\t%d b", r.name().c_str(), size);
	return size;
}
int approxConstantCovariate(const ConstantCovariate& r) {
	int size = sizeof(r);
	// lvalues
	size += r.pActorSet()->n() * sizeof(double);
	// lmissing
	size += r.pActorSet()->n() * sizeof(double);
	LOGF(Priority::DEBUG, "ConstantCovariate\t%s\t%d b", r.name().c_str(), size);
	return size;
}
int approxChangingCovariate(const ChangingCovariate& r, int obs) {
	int size = sizeof(r);
	// lvalues
	size += obs * r.pActorSet()->n() * sizeof(double);
	// lmissing
	size += obs * r.pActorSet()->n() * sizeof(double);
	LOGF(Priority::DEBUG, "ChangingCovariate\t%s\t%d b", r.name().c_str(), size);
	return size;
}
int approxSize(const Data& r) {
	int size = sizeof(r);
	typedef const vector<const ActorSet*> ActorPV;
	ActorPV& a = r.rActorSets();
	for (ActorPV::const_iterator it = a.begin(); it != a.end(); ++it) {
		size += approxActorSet(**it);
	}
	typedef const vector<LongitudinalData*> DataPV;
	DataPV& b = r.rDependentVariableData();
	LOGF(Priority::DEBUG, "LongitudinalData\t\t%d #", b.size());
	for (DataPV::const_iterator it = b.begin(); it != b.end(); ++it) {
		size += approxLongitudinalData(**it);
	}
	typedef const vector<ConstantCovariate*> CCovPV;
	CCovPV& c = r.rConstantCovariates();
	LOGF(Priority::DEBUG, "ConstantCovariate\t\t%d #", c.size());
	for (CCovPV::const_iterator it = c.begin(); it != c.end(); ++it) {
		size += approxConstantCovariate(**it);
	}
	typedef const vector<ChangingCovariate*> CovPV;
	CovPV& d = r.rChangingCovariates();
	LOGF(Priority::DEBUG, "ChangingCovariate\t\t%d #", d.size());
	for (CovPV::const_iterator it = d.begin(); it != d.end(); ++it) {
		size += approxChangingCovariate(**it, r.observationCount());
	}
	typedef const vector<ConstantDyadicCovariate*> CDCovPV;
	CDCovPV& e = r.rConstantDyadicCovariates();
	LOGF(Priority::DEBUG, "ConstantDyadicCovariate\t\t%d #", e.size());
	typedef const vector<ChangingDyadicCovariate*> DCovPV;
	DCovPV& f = r.rChangingDyadicCovariates();
	LOGF(Priority::DEBUG, "ChangingDyadicCovariate\t\t%d #", f.size());

	LOGF(Priority::DEBUG, "Data\t\t%d b", size);
	LOGF(Priority::DEBUG, "Data\t\t%d mb", size / 1024 / 1024);
	return size;
}
