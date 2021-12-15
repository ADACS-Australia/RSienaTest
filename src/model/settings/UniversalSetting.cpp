/*
 * UniversalSetting.cpp
 *
 *  Created on: 18.06.2014
 *      Author: ortmann
 */

#include "UniversalSetting.h"
#include "../../network/Network.h"
#include "../../network/iterators/IntVecIterator.h"
#include "../../logger/Logger.h"

using namespace siena::logger;

namespace siena {

UniversalSetting::UniversalSetting() :
		GeneralSetting(), //
		rSteps() {
}

UniversalSetting::~UniversalSetting() {
}

void UniversalSetting::initSetting(Network* const lpNetwork) {
	if (!rSteps.empty()) {
		LOGS(Priority::ERROR)<< "setting has not been terminated\n";
		throw "setting has not been terminated";
	}
	rSteps.reserve(lpNetwork->m());
	for (int i = 0; i < lpNetwork->m(); i++) {
		rSteps.push_back(i);
	}
}

void UniversalSetting::terminateSetting(Network* const lpNetwork) {
	rSteps.clear();
}

void UniversalSetting::initSetting() {
	return;
}

ITieIterator* UniversalSetting::getSteps() {
	if (!rSteps.empty()) {
		return new IntVecIterator(rSteps.begin(), rSteps.end());
	}
	LOGS(Priority::ERROR)<< "the setting has not been initialized";
	throw "the setting has not been initialized";
}

int UniversalSetting::getSize() {
	return rSteps.size();
}

} /* namespace siena */
