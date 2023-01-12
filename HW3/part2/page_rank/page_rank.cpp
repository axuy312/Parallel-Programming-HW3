#include "page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"
#include <omp.h>

// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double *solution, double damping, double convergence)
{

  // initialize vertex weights to uniform probability. Double
  // precision scores are used to avoid underflow for large graphs

  int numNodes = num_nodes(g);
  double equal_prob = 1.0 / numNodes;
  for (int i = 0; i < numNodes; ++i)
  {
    solution[i] = equal_prob;
  }

  /*
     For PP students: Implement the page rank algorithm here.  You
     are expected to parallelize the algorithm using openMP.  Your
     solution may need to allocate (and free) temporary arrays.

     Basic page rank pseudocode is provided below to get you started:
  */
	// initialization: see example code above
	double *score_old = new double[numNodes];
	bool converged = false;
	while (!converged) {
		double no_outgoing_sum = 0.0;
		double global_diff = 0.0;
		for (int vi = 0; vi < numNodes; ++vi)
		{
			score_old[vi] = solution[vi];
			if (outgoing_size(g, vi)==0)
				no_outgoing_sum += damping * score_old[vi] / numNodes;
		}
		#pragma omp parallel for
		for (int vi = 0; vi < numNodes; ++vi)
		{
			solution[vi] = 0.0;
			const Vertex* start = incoming_begin(g, vi);
			int inSize = incoming_size(g, vi);
			for (int j = 0; j < inSize; ++j)
			{
				int vj = *(start + j);
				solution[vi] += (score_old[vj] / outgoing_size(g, vj));
			}
			solution[vi] = (damping * solution[vi]) + (1.0 - damping) / numNodes;
			solution[vi] += no_outgoing_sum;
		}
		for (int vi = 0; vi < numNodes; ++vi)
		{
			global_diff += fabs(score_old[vi] - solution[vi]);
		}
		converged = (global_diff < convergence);
		// compute score_new[vi] for all nodes vi:
		//score_new[vi] = sum over all nodes vj reachable from incoming edges
        //                  { score_old[vj] / number of edges leaving vj  }
       //score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

       //score_new[vi] += sum over all nodes v in graph with no outgoing edges
       //                   { damping * score_old[v] / numNodes }

       // compute how much per-node scores have changed
       // quit once algorithm has converged

       //global_diff = sum over all nodes vi { abs(score_new[vi] - score_old[vi]) };
       //converged = (global_diff < convergence)
	}
}
