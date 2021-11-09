#ifndef GRAPH_OPERATIONS_H
#define GRAPH_OPERATIONS_H
#include "graph/graph.h"

class GraphOperations {
public:
	static void getKCore(const Graph *graph, int *core_table);
    static void bfsTraversal(const Graph *graph, VertexID root_vertex, TreeNode *&tree, VertexID *&bfs_order);
};

#endif