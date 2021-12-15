/*
 * PrimarySetting.h
 *
 *  Created on: 18.06.2014
 *      Author: ortmann
 */

#ifndef BASICPRIMARYSETTING_H_
#define BASICPRIMARYSETTING_H_

#include "GeneralSetting.h"
#include "network/layers/DistanceTwoLayer.h"
#include "network/layers/PrimaryLayer.h"

namespace siena {

class Network;

class BasicPrimarySetting: public GeneralSetting {
public:
	BasicPrimarySetting();

	virtual ~BasicPrimarySetting();

	/** @copydoc ASetting::initSetting(Network* const lpNetwork) */
	virtual void initSetting(Network* const lpNetwork);
	/** @copydoc ASetting::terminateSetting(Network* const lpNetwork) */
	virtual void terminateSetting(Network* const lpNetwork);

	const Network * pPrimaryNetwork() const;

	ITieIterator* getSteps();

	int getSize();


protected:
	void initSetting(); // called on each initSetting(int ego)
	void terminateSetting();

private:
	Network* lpNetwork; // the underlying dependent

	// old
	// DistanceTwoLayer rDistTwoLayer;
	// ITieIterator* lpiter;

	// new
	PrimaryLayer lPLayer; // the 2-path extension

};

} /* namespace siena */
#endif /* BASICPRIMARYSETTING_H_ */
