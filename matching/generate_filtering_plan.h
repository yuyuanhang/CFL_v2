#ifndef GENERATE_FILTERING_PLAN_H
#define GENERATE_FILTERING_PLAN_H
#include <queue>
#include "graph/graph.h"
#include "configuration/types.h"

class GenerateFilteringPlan {
public:
    static void generateCFLFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                      VertexID *&order, int &level_count, ui *&level_offset);

private:
    static VertexID selectCFLFilterStartVertex(const Graph *data_graph, const Graph *query_graph);

};

#endif