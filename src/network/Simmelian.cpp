/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *****************************************************************************/
#include "Simmelian.h"

#include "logger/Logger.h"
#include "network/IncidentTieIterator.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"

#include <algorithm> // sort
#include <vector>

namespace siena {
	using namespace logger;

	const OneModeNetwork* newSimmelian(const OneModeNetwork* pNetwork) {
		// assuming binary networks: the minimun is full reciprocity
		const OneModeNetwork* pReciprocal = newUndirected(pNetwork, std::min);
		// simmelian strength: number of fully reciprocated triads
		const OneModeNetwork* pSimmelian = newIncidentTriads(pReciprocal);
		delete pReciprocal;
		return pSimmelian;
	}

	const OneModeNetwork* newUndirected(const OneModeNetwork* pNetwork,
			const int& (*aggregate)(const int&, const int&)) {
		OneModeNetwork* pUndirected = new OneModeNetwork(pNetwork->n(), false);
		// int count = 0;
		// LOGS(Priority::INFO) << "network has " << pNetwork->tieCount() << " ties";
		for (TieIterator tie = pNetwork->ties(); tie.valid(); tie.next()) {
			const int src = tie.ego();
			const int snk = tie.alter();
			const int val = (*aggregate)(tie.value(), pNetwork->tieValue(snk, src));
			if (val != 0) {
				// count += 1;
				pUndirected->setTieValue(src, snk, val);
				pUndirected->setTieValue(snk, src, val);
			}
		}
		// LOGS(Priority::INFO) << "undirected has " << count << " ties";
		return pUndirected;
	}

	/**
	 * Chiba and Nishizeki. Could by faster with bucket sort or Marks stuff.
	 * First run it on RSiena size and see if more speed is needed.
	 */
	const OneModeNetwork* newIncidentTriads(const OneModeNetwork* pNetwork) {
		const int n = pNetwork->n();
		// sort nodes descending by degree
		std::vector<std::pair<int, int> > ascending;
		for (int i = 0; i < n; i++) {
			ascending.push_back(std::make_pair(-pNetwork->outDegree(i), i));
		}
		std::sort(ascending.begin(), ascending.end());
		// output: triad count valued network
		OneModeNetwork* pIncidentTriads = new OneModeNetwork(n, false);
		std::vector<bool> removed(n); // flag nodes as "removed" from the network
		std::vector<bool> marked(n);  // marker
		for (int vi = 0; vi < n; vi++) {
			const int v = ascending[vi].second;
			if (removed[v]) continue;
			// mark adjacent nodes
			std::fill(marked.begin(), marked.end(), false);
			for (IncidentTieIterator ui = pNetwork->outTies(v); ui.valid(); ui.next()) {
				marked[ui.actor()] = true;
			}
			// for each marked (adjacent again, skipping cleared marks)
			for (IncidentTieIterator ui = pNetwork->outTies(v); ui.valid(); ui.next()) {
				const int u = ui.actor();
				// LOGS(Priority::INFO) << "node: " << v << " -> " << u;
				if (removed[u] || !marked[u]) continue;
				for (IncidentTieIterator wi = pNetwork->outTies(u); wi.valid(); wi.next()) {
					const int w = wi.actor();
					if (removed[w] || !marked[w]) continue;
					// LOGS(Priority::INFO) << "triad: " << v << " " << u << " " << w;
					pIncidentTriads->increaseTieValue(v, u, 1);
					pIncidentTriads->increaseTieValue(u, w, 1);
//					pIncidentTriads->increaseTieValue(w, u, 1); // this now is replaced by the following line (TS)
					pIncidentTriads->increaseTieValue(w, v, 1);
					pIncidentTriads->increaseTieValue(u, v, 1);
					pIncidentTriads->increaseTieValue(w, u, 1);
					pIncidentTriads->increaseTieValue(v, w, 1);
					marked[u] = false;
				}
			}
			removed[v] = true;
		}
		return pIncidentTriads;
	}

} // namespace siena
