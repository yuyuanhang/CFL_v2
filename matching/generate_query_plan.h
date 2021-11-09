#ifndef GENERATE_QUERY_PLAN_H
#define GENERATE_QUERY_PLAN_H
#include "graph/graph.h"
#include <vector>

class GenerateQueryPlan {
public:
    static void
    generateCFLQueryPlan(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix,
                             ui *&order, ui *&pivot, TreeNode *tree, ui *bfs_order, ui *candidates_count);

    static void checkQueryPlanCorrectness(const Graph* query_graph, ui* order, ui* pivot);

    static void printSimplifiedQueryPlan(const Graph *query_graph, ui *order);
    
private:
    static void estimatePathEmbeddsingsNum(std::vector<ui> &path, Edges ***edge_matrix,
                                           std::vector<size_t> &estimated_embeddings_num);

    static void generateCorePaths(const Graph* query_graph, TreeNode* tree_node, VertexID cur_vertex, std::vector<ui> &cur_core_path,
                                  std::vector<std::vector<ui>> &core_paths);

    static void generateTreePaths(const Graph* query_graph, TreeNode* tree_node, VertexID cur_vertex,
                                  std::vector<ui> &cur_tree_path, std::vector<std::vector<ui>> &tree_paths);

    static void generateLeaves(const Graph* query_graph, std::vector<ui>& leaves);

    static ui generateNoneTreeEdgesCount(const Graph *query_graph, TreeNode *tree_node, std::vector<ui> &path);
};
#endif