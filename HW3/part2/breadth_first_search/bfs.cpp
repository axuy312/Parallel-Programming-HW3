#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

#include <algorithm> 

#include "../common/CycleTimer.h"
#include "../common/graph.h"


#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1
//#define VERBOSE3

void vertex_set_clear(vertex_set *list)
{
    list->count = 0;
}

void vertex_set_init(vertex_set *list, int count)
{
    list->max_vertices = count;
    list->vertices = (int *)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set *frontier,
    int *distances,
	int frontier_level)
{
	int next_frontier_step = 0;
	
	#pragma omp parallel for reduction(+: next_frontier_step)
	for (int u = 0; u < g->num_nodes; u++)
    {
		if (frontier->vertices[u] == frontier_level)
		{
			const Vertex* start = outgoing_begin(g, u);
			int outSize = outgoing_size(g, u);

			// attempt to add all neighbors to the new frontier
			for (int iv = 0; iv < outSize; iv++)
			{
				int v = *(start + iv);

				if (distances[v] == NOT_VISITED_MARKER)
				{
					distances[v] = distances[u] + 1;
					frontier->vertices[v] = frontier_level + 1;
					next_frontier_step++;
				}
			}
		}
    }
	
	frontier->count += next_frontier_step;
}

void bottom_up_step(
    Graph g,
    vertex_set *frontier,
    int *distances,
	int frontier_level)
	{
		int next_frontier_step = 0;
		
		#pragma omp parallel for reduction(+: next_frontier_step)
		for(int v = 0; v < g->num_nodes; v++)
		{
			if (distances[v] == NOT_VISITED_MARKER)
			{
				const Vertex* start = incoming_begin(g, v);
				int inSize = incoming_size(g, v);
				for (int iu = 0; iu < inSize; iu++)
				{
					int u = *(start + iu);
					if (frontier->vertices[u] == frontier_level)
					{
						distances[v] = distances[u] + 1;
						frontier->vertices[v] = frontier_level + 1;
						next_frontier_step++;
						break;
					}
				}
			}
		}
		
		frontier->count += next_frontier_step;
	}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution *sol)
{

    vertex_set list1;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set *frontier = &list1;
	int frontier_level = 1;

    // initialize all nodes to NOT_VISITED
	#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
	{
        sol->distances[i] = NOT_VISITED_MARKER;
	}

    // setup frontier with the root node
    frontier->vertices[ROOT_NODE_ID] = frontier_level;
	frontier->count++;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0)
    {

#ifdef VERBOSE1
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(frontier);

        top_down_step(graph, frontier, sol->distances, frontier_level);
		frontier_level++;

#ifdef VERBOSE1
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.10f sec\n", frontier->count, end_time - start_time);
#endif
		
    }
}

void bfs_bottom_up(Graph graph, solution *sol)
{
    // For PP students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.

    vertex_set list1;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set *frontier = &list1;
	int frontier_level = 1;
	

    // initialize all nodes to NOT_VISITED
	#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
	{
        sol->distances[i] = NOT_VISITED_MARKER;
	}

    // setup frontier with the root node
    frontier->vertices[ROOT_NODE_ID] = frontier_level;
	frontier->count++;
    sol->distances[ROOT_NODE_ID] = 0;
	

    while (frontier->count != 0)
    {

#ifdef VERBOSE2
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(frontier);

        bottom_up_step(graph, frontier, sol->distances, frontier_level);
		frontier_level++;

#ifdef VERBOSE2
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.10f sec\n", frontier->count, end_time - start_time);
#endif

    }
}

void bfs_hybrid(Graph graph, solution *sol)
{
    // For PP students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
    // For PP students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.

    vertex_set list1;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set *frontier = &list1;
	int frontier_level = 1;
	

    // initialize all nodes to NOT_VISITED
	#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
	{
        sol->distances[i] = NOT_VISITED_MARKER;
	}

    // setup frontier with the root node
    frontier->vertices[ROOT_NODE_ID] = frontier_level;
	frontier->count++;
    sol->distances[ROOT_NODE_ID] = 0;
	

    while (frontier->count != 0)
    {

#ifdef VERBOSE3
        double start_time = CycleTimer::currentSeconds();
#endif

		if (frontier->count < 7800000)//9770000
		{
			vertex_set_clear(frontier);
			top_down_step(graph, frontier, sol->distances, frontier_level);
		}
		else
		{
			vertex_set_clear(frontier);
			bottom_up_step(graph, frontier, sol->distances, frontier_level);\
		}
		frontier_level++;

#ifdef VERBOSE3
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.10f sec\n", frontier->count, end_time - start_time);
#endif

    }
}
