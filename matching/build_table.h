#ifndef BUILD_TABLE_H
#define BUILD_TABLE_H
#include "graph/graph.h"
#include <vector>

class BuildTable {
public:
    static void buildTables(const Graph* data_graph, const Graph* query_graph, ui** candidates, ui* candidates_count,
                            Edges*** edge_matrix);
    static void printTableCardinality(const Graph* query_graph, Edges*** edge_matrix);
    static size_t computeMemoryCostInBytes(const Graph* query_graph, ui* candidates_count, Edges*** edge_matrix);
};

#endif