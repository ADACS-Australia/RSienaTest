/**************************************************************************//**
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * \file Simulation.cpp
 * Implementation of Simulation.h
 *****************************************************************************/

#include "sim/Simulation.h"

#include "modificator/ResultModificator.h"
#include "listener/ResultListener.h"

#include "logger/Logger.h"
#include "data/LongitudinalData.h"
#include "model/State.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"

using namespace std;
using namespace Eigen;
using namespace siena::logger;

namespace siena {

// Dummy effect info for basic rate parameters.
EffectInfo Simulation::RATE_EFFECT("", "BasicRate", "rate", 0, 0, "", "",
		"dummy");

/**
 * Constructs the Simulation object.
 *
 * @param pModel Pointer to the Model.
 * @param pData Pointer to the Data.
 */
Simulation::Simulation(Model* pModel, Data* pData) :
		lpModel(pModel), //
		lpData(pData), //
		lResultModificatorPtrs(), //
		lResultListenerPtrs(), //
		lNeeds() {
	if (pModel != 0/*ptr*/&& pData != 0/*ptr*/) {
		lSimulationEffects = createSimulationEffectVector();
		lTheta = VectorXd(lSimulationEffects.size());
		lStatisticEffects = createStatisticEffectVector();
		LOGF(Priority::INFO, "%d simulation effects, %d statistic effects",
				lSimulationEffects.size(), lStatisticEffects.size());
		updateLocalParametersFromModel();
	} else {
		LOG(Priority::ERROR,
				"created simulation with null pointer to model or data");
	}
	LOGS(Priority::DEBUG);
}

Simulation::~Simulation() {
}

/**
 * Adds a result modificator to the simulation.
 *
 * @param pModificator The modificator to add.
 */
void Simulation::addResultModificator(ResultModificator* const pModificator) {
	LOGS(Priority::DEBUG);
	list<ResultModificator*>::iterator it = find(lResultModificatorPtrs.begin(),
			lResultModificatorPtrs.end(), pModificator);
	if (it == lResultModificatorPtrs.end()) {
		lResultModificatorPtrs.push_back(pModificator);
		lNeeds.insert(pModificator->needs().begin(),
				pModificator->needs().end());
		needsChanged();
	}
}

/**
 * Removes a result modificator from the simulation.
 *
 * @param pModificator The modificator to remove.
 */
void Simulation::removeResultModificator(
		ResultModificator* const pModificator) {
	LOGS(Priority::DEBUG);
	list<ResultModificator*>::iterator it = find(lResultModificatorPtrs.begin(),
			lResultModificatorPtrs.end(), pModificator);
	if (it != lResultModificatorPtrs.end()) {
		lResultModificatorPtrs.erase(it);
	}
	lNeeds.clear();
	for (it = lResultModificatorPtrs.begin();
			it != lResultModificatorPtrs.end(); ++it) {
		lNeeds.insert((*it)->needs().begin(), (*it)->needs().end());
	}
	needsChanged();
}

/**
 * @param rResult Result access to distribute.
 */
void Simulation::fireResultModificators(Result& rResult) {
	LOGS(Priority::DEBUG)<<"run "<<lResultModificatorPtrs.size()<<" modificators";
	for (list<ResultModificator*>::iterator it = lResultModificatorPtrs.begin();
			it != lResultModificatorPtrs.end(); ++it) {
		(*it)->onResults(rResult);
	}
}

/**
 * Adds a result listener to the simulation.
 *
 * The union of all result sets requested by the listeners is then broad
 * casted when the simulation has finished.
 *
 * @param pListener The listener to add.
 */
void Simulation::addResultListener(ResultListener* const pListener) {
	LOGS(Priority::DEBUG);
	vector<ResultListener*>::iterator it = find(lResultListenerPtrs.begin(),
			lResultListenerPtrs.end(), pListener);
	if (it == lResultListenerPtrs.end()) {
		lResultListenerPtrs.push_back(pListener);
		lNeeds.insert(pListener->needs().begin(), pListener->needs().end());
		needsChanged();
	}
}

/**
 * Removes a result listener from the simulation.
 *
 * @param pListener The listener to remove.
 */
void Simulation::removeResultListener(ResultListener* const pListener) {
	LOGS(Priority::DEBUG);
	vector<ResultListener*>::iterator it = find(lResultListenerPtrs.begin(),
			lResultListenerPtrs.end(), pListener);
	if (it != lResultListenerPtrs.end()) {
		lResultListenerPtrs.erase(it);
	}
	lNeeds.clear();
	for (it = lResultListenerPtrs.begin(); it != lResultListenerPtrs.end();
			++it) {
		lNeeds.insert((*it)->needs().begin(), (*it)->needs().end());
	}
	needsChanged();
}

/**
 * Clear the result listeners and returns a copy of the list with the old
 * listener pointers.
 *
 * @return Copy of the current listeners.
 */
vector<ResultListener*> Simulation::clearResultListener() {
	LOGS(Priority::DEBUG);
	vector<ResultListener*> copy(lResultListenerPtrs);
	lResultListenerPtrs.clear();
	lNeeds.clear();
	needsChanged();
	return copy;
}

/**
 * @param rResult Result access to distribute.
 */
void Simulation::fireResult(const Result& rResult) {
	LOGS(Priority::DEBUG)<<"distribute to "<<lResultListenerPtrs.size()<<" listeners";
	for (vector<ResultListener*>::iterator it = lResultListenerPtrs.begin();
			it != lResultListenerPtrs.end(); ++it) {
		(*it)->onResults(rResult);
	}
}

/**
 * @return Pointer to the Model object.
 */
Model* Simulation::pModel() {
	return lpModel;
}

/**
 * @return Pointer to the Data object.
 */
Data* Simulation::pData() {
	return lpData;
}

/**
 * @return number of periods.
 */
int Simulation::nPeriods() {
	return lpData->observationCount() - 1;
}

/**
 * Return the simulation effects.
 *
 * If you need to know which of the is a basic rate effect compare it with the
 * RATE_EFFECT dummy.
 * \code
 * (rEffects().array() == &Simulation::RATE_EFFECT).select(a, b);
 * \endcode
 *
 * @return Reference to the simulation effects.
 */
const VectorXEffectPtr& Simulation::rSimulationEffects() {
	return lSimulationEffects;
}

/**
 * @return The number of statistics / statistic effects.
 */
int Simulation::nStatistics() {
	return lStatisticEffects.size();
}

/**
 * Return the statistic effects.
 *
 * If you need to know which of the is a basic rate effect compare it with the
 * RATE_EFFECT dummy.
 * \code
 * (rEffects().array() == &Simulation::RATE_EFFECT).select(a, b);
 * \endcode
 *
 * @return Reference to the statistic effects.
 */
const VectorXEffectPtr& Simulation::rStatisticEffects() {
	return lStatisticEffects;
}

/**
 * @return the number of parameters.
 */
int Simulation::nParameters() {
	return lTheta.size();
}

/**
 * Common name: @f$ \theta @f$, theta
 *
 * @return Vector of model parameters. First number-of-period coefficients are
 * the basic rate parameters, after that the non-basic rate, evaluation,
 * endownment, creation effects.
 */
const VectorXd& Simulation::rParameters() {
	return lTheta;
}

/**
 * Fill the parameter vector with the current model parameters.
 */
void Simulation::updateLocalParametersFromModel() {
	int i = 0; // Parameter index

	// For each dependent variable
	const vector<LongitudinalData*>& vars = lpData->rDependentVariableData();
	for (vector<LongitudinalData*>::size_type v = 0; v < vars.size(); ++v) {
		// Basic rate parameters
		if (!lpModel->conditional()
				|| vars[v]->name() != lpModel->conditionalDependentVariable()) {
			for (int e = 0; e < vars[v]->observationCount() - 1; ++e) {
				lTheta[i] = lpModel->basicRateParameter(vars[v], e);
				++i;
			}
		}
		// Other effects
		const int num = nNonBasicRateEffects(vars[v]->name());
		for (int e = 0; e < num; ++e) {
			lTheta[i] = lSimulationEffects[i]->parameter();
			++i;
		}
	} // for vars[v]

	assert(i == lTheta.size());
}

/**
 * Set parameters of the local model.
 *
 * Default implementation for updateParamters.
 *
 * Similar to: siena07internals.cpp updateParameters()
 *
 * @param rTheta Parameter vector.
 */
void Simulation::updateLocalParameters(const VectorXd& rTheta) {
	int i = 0; // Parameter index

	// For each dependent variable
	const vector<LongitudinalData*>& vars = lpData->rDependentVariableData();
	for (vector<LongitudinalData*>::size_type v = 0; v < vars.size(); ++v) {
		// Basic rate parameters
		if (!lpModel->conditional()
				|| vars[v]->name() != lpModel->conditionalDependentVariable()) {
			for (int e = 0; e < vars[v]->observationCount() - 1; ++e) {
				lpModel->basicRateParameter(vars[v], e, rTheta[i]);
				lTheta[i] = rTheta[i];
				++i;
			}
		}
		// Other effects
		const int num = nNonBasicRateEffects(vars[v]->name());
		for (int e = 0; e < num; ++e) {
			lSimulationEffects[i]->parameter(rTheta[i]);
			lTheta[i] = rTheta[i];
			++i;
		}
	} // for vars[v]

	assert(i == rTheta.size());
}

/**
 * Helper function to create the effect vector.
 *
 * The basic rate and other rate effects of the variable are appended to the
 * effect vector.
 *
 * @param[out] rEffects
 * @param pVar
 */
void Simulation::addRateEffects(vector<EffectInfo*>& rEffects,
		LongitudinalData* pVar) {
	if (!lpModel->conditional()
			|| pVar->name() != lpModel->conditionalDependentVariable()) {
		for (int e = 0; e < pVar->observationCount() - 1; ++e) {
			rEffects.push_back(&RATE_EFFECT);
		}
	}
	Effects& rRate = lpModel->rRateEffects(pVar->name());
	rEffects.insert(rEffects.end(), rRate.begin(), rRate.end());
}

/**
 * Helper function to create the effect vector.
 *
 * The effects of the objective function of the variable are appended to the
 * effect vector.
 *
 * @param[out] rEffects
 * @param pVar
 */
void Simulation::addObjectiveEffects(vector<EffectInfo*>& rEffects,
		LongitudinalData* pVar) {
	Effects& rEval = lpModel->rEvaluationEffects(pVar->name());
	rEffects.insert(rEffects.end(), rEval.begin(), rEval.end());
	Effects& rEndowment = lpModel->rEndowmentEffects(pVar->name());
	rEffects.insert(rEffects.end(), rEndowment.begin(), rEndowment.end());
	Effects& rCreation = lpModel->rCreationEffects(pVar->name());
	rEffects.insert(rEffects.end(), rCreation.begin(), rCreation.end());
}

/**
 * Helper function to create the effect vector.
 *
 *  If the model is unconditional the first number-of-period coefficients are
 *  dummy informations, since the basic rates are handled in a slightly
 *  different way.
 *
 *  If the model is a GMM model the evaluation, creation and endowment effects
 *  are replaced by the GMM effects.
 *
 * @return Vector of effect informations.
 */
VectorXEffectPtr Simulation::createSimulationEffectVector() {
	vector<EffectInfo*> effects;
	// For each dependent variable
	const vector<LongitudinalData*>& vars = lpData->rDependentVariableData();
	for (vector<LongitudinalData*>::size_type v = 0; v < vars.size(); ++v) {
		addRateEffects(effects, vars[v]);
		addObjectiveEffects(effects, vars[v]);
	}
	return Map<VectorXEffectPtr>(&effects[0], effects.size(), 1);
}

VectorXEffectPtr Simulation::createStatisticEffectVector() {
	vector<EffectInfo*> effects;
	// For each dependent variable
	const vector<LongitudinalData*>& vars = lpData->rDependentVariableData();
	for (vector<LongitudinalData*>::size_type v = 0; v < vars.size(); ++v) {
		addRateEffects(effects, vars[v]);
		Effects& rGMM = lpModel->rGMMEffects(vars[v]->name());
		if (!rGMM.empty()) {
			effects.insert(effects.end(), rGMM.begin(), rGMM.end());
		} else {
			addObjectiveEffects(effects, vars[v]);
		}
	}
	return Map<VectorXEffectPtr>(&effects[0], effects.size(), 1);
}

/**
 * Helper function to calculate the number of effects that are not basic rate
 * effects.
 *
 * @param variableName Name of the network variable.
 * @return Number of effects that are not basic rate effects.
 */
int Simulation::nNonBasicRateEffects(const string& variableName) {
	return lpModel->rRateEffects(variableName).size()
			+ lpModel->rEvaluationEffects(variableName).size()
			+ lpModel->rEndowmentEffects(variableName).size()
			+ lpModel->rCreationEffects(variableName).size();
}

MatrixXd Simulation::calculateTargetStatistics() {
	const int M = lpData->observationCount() - 1;
	MatrixXd stats(M, lStatisticEffects.size());
	// For each period calculate the target statistics
	for (int m = 0; m < M; ++m) {
		State state(lpData, m + 1);
		StatisticCalculator calculator(lpData, lpModel, &state, m);
		MatrixXd::RowXpr row = stats.row(m);
		getSingleStatistics(row, m, calculator);
	}
	return stats;
}

} // namespace siena
