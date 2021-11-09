#include <queue>
#include "generate_filtering_plan.h"
#include "filter_vertices.h"
#include "utility/graph_operations.h"

void GenerateFilteringPlan::generateCFLFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                  VertexID *&order, int &level_count, ui *&level_offset) {
    ui query_vertices_num = query_graph->getVerticesCount();
    // The vertex that has less candidate vertices (with NLF filter) and more neighbors is returned
    VertexID start_vertex = selectCFLFilterStartVertex(data_graph, query_graph);
    // get a BFS tree and a vertex just has one parent although it can have more than one neighbor on the last level.
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, order);

    // an inverted index to record the order of vertex i.
    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    // record the number of vertices on each level.
    // level count means the current level.
    level_count = -1;
    level_offset = new ui[query_vertices_num + 1];

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        if (tree[u].level_ != level_count) {
            level_count += 1;
            level_offset[level_count] = 0;
        }

        level_offset[level_count] += 1;

        ui u_in_nbrs_count;
        const VertexID* u_in_nbrs = query_graph->getVertexInNeighbors(u, u_in_nbrs_count);
        // For each neighbor of current vertex,
        for (ui j = 0; j < u_in_nbrs_count; ++j) {
            VertexID u_in_nbr = u_in_nbrs[j];

            // bn neighbor means that this neighbor is in the front of current vertex in the bfs order.
            // fn neighbor means that this neighbor is behind current vertex in the bfs order.
            // under level neighbor means that this neighbor is behind current vertex in the bfs order and it is on higher level than current vertex.
            if (tree[u].level_ == tree[u_in_nbr].level_) {
                if (order_index[u_in_nbr] < order_index[u]) {
                    tree[u].bn_[tree[u].bn_count_] = u_in_nbr;
                    tree[u].bn_directions_[tree[u].bn_count_] = IN;
                    tree[u].bn_count_++;
                }
                else {
                    tree[u].fn_[tree[u].fn_count_] = u_in_nbr;
                    tree[u].fn_directions_[tree[u].fn_count_] = IN;
                    tree[u].fn_count_++;
                }
            }
            else if (tree[u].level_ > tree[u_in_nbr].level_) {
                tree[u].bn_[tree[u].bn_count_] = u_in_nbr;
                tree[u].bn_directions_[tree[u].bn_count_] = IN;
                tree[u].bn_count_++;
            }
            else {
                tree[u].under_level_[tree[u].under_level_count_] = u_in_nbr;
                tree[u].under_level_directions_[tree[u].under_level_count_] = IN;
                tree[u].under_level_count_++;
            }
        }

        ui u_out_nbrs_count;
        const VertexID* u_out_nbrs = query_graph->getVertexOutNeighbors(u, u_out_nbrs_count);
        // For each neighbor of current vertex,
        for (ui j = 0; j < u_out_nbrs_count; ++j) {
            VertexID u_out_nbr = u_out_nbrs[j];

            // bn neighbor means that this neighbor is in the front of current vertex in the bfs order.
            // fn neighbor means that this neighbor is behind current vertex in the bfs order.
            // under level neighbor means that this neighbor is behind current vertex in the bfs order and it is on higher level than current vertex.
            if (tree[u].level_ == tree[u_out_nbr].level_) {
                if (order_index[u_out_nbr] < order_index[u]) {
                    tree[u].bn_[tree[u].bn_count_] = u_out_nbr;
                    tree[u].bn_directions_[tree[u].bn_count_] = OUT;
                    tree[u].bn_count_++;
                }
                else {
                    tree[u].fn_[tree[u].fn_count_] = u_out_nbr;
                    tree[u].fn_directions_[tree[u].fn_count_] = OUT;
                    tree[u].fn_count_++;
                }
            }
            else if (tree[u].level_ > tree[u_out_nbr].level_) {
                tree[u].bn_[tree[u].bn_count_] = u_out_nbr;
                tree[u].bn_directions_[tree[u].bn_count_] = OUT;
                tree[u].bn_count_++;
            }
            else {
                tree[u].under_level_[tree[u].under_level_count_] = u_out_nbr;
                tree[u].under_level_directions_[tree[u].under_level_count_] = OUT;
                tree[u].under_level_count_++;
            }
        }
    }

    level_count += 1;

    // compute prefix sum.
    ui prev_value = 0;
    for (int i = 1; i <= level_count; ++i) {
        ui temp = level_offset[i];
        level_offset[i] = level_offset[i - 1] + prev_value;
        prev_value = temp;
    }
    level_offset[0] = 0;
}

VertexID GenerateFilteringPlan::selectCFLFilterStartVertex(const Graph *data_graph, const Graph *query_graph) {
    auto rank_compare = [](std::pair<VertexID, double> l, std::pair<VertexID, double> r) {
        return l.second < r.second;
    };

    std::priority_queue<std::pair<VertexID, double>, std::vector<std::pair<VertexID, double>>, decltype(rank_compare)> rank_queue(rank_compare);

    // Compute the ranking.
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;

        // If there is a 2-core, we only consider vertices in 2-core; Otherwise, we consider all vertices.
        if (query_graph->get2CoreSize() == 0 || query_graph->getCoreValue(query_vertex) > 1) {
            LabelID label = query_graph->getVertexLabel(query_vertex);
            ui degree = query_graph->getVertexDegree(query_vertex);
            ui frequency = data_graph->getLabelsFrequency(label);

            // Vertices considered are ranked by label frequency/degree (consider both cardinality and structure).
            double rank = frequency / (double) degree;
            rank_queue.push(std::make_pair(query_vertex, rank));
        }
    }

    // Keep the top-3.
    while (rank_queue.size() > 3) {
        rank_queue.pop();
    }

    VertexID start_vertex = 0;
    double min_score = data_graph->getGraphMaxLabelFrequency() + 1;

    while (!rank_queue.empty()) {
        VertexID query_vertex = rank_queue.top().first;
        ui count;
        // Now, candidate sets of top-3 vertices are filtered further based on NLF (neighbor label frequency).
        FilterVertices::computeCandidateWithNLF(data_graph, query_graph, query_vertex, count);
        // And, the score is computed again.
        double cur_score = count / (double) query_graph->getVertexDegree(query_vertex);

        if (cur_score < min_score) {
            start_vertex = query_vertex;
            min_score = cur_score;
        }
        rank_queue.pop();
    }

    return start_vertex;
}