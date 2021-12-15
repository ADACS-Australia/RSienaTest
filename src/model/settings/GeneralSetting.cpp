/*
 * GeneralSetting.cpp
 *
 *  Created on: 24.07.2014
 *      Author: ortmann
 */

#include "GeneralSetting.h"

#include "../../logger/Logger.h"
#include "../../network/iterators/UnionTieIterator.h"
#include "../../network/iterators/SingleIterator.h"
#include "../../network/iterators/FilteredIterator.h"

using namespace siena::logger;

namespace siena {

GeneralSetting::~GeneralSetting() {
	if (lppermittedIter != 0) {
		delete lppermittedIter;
	}
}

ITieIterator* GeneralSetting::getPermittedSteps() {
	return lppermittedIter->clone();
}

int GeneralSetting::getPermittedSize() {
	return lppermittedIter->size();
}

GeneralSetting::GeneralSetting() :
		Setting(), //
		lppermittedIter(0) {
}

void GeneralSetting::terminateSetting() {
	if (lppermittedIter != 0) {
		delete lppermittedIter;
		lppermittedIter = 0;
	} else {
		LOGS(Priority::ERROR)<< "the setting has not been initialized";
		throw "the setting has not been initialized";
	}
}

void GeneralSetting::initPermittedSteps(const bool* const permitted) {
	if (lppermittedIter == 0) {
	ITieIterator* iter = getSteps();
	lppermittedIter = new FilteredIterator(*iter,permitted);
	delete iter;
	} else {
		LOGS(Priority::ERROR)<<"the setting has not been terminated";
		throw "the setting has not been terminated";
}
}

}
/* namespace siena */
