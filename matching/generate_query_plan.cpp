#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include "generate_query_plan.h"

void GenerateQueryPlan::generateCFLQueryPlan(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix,
                                             ui *&order, ui *&pivot, TreeNode *tree, ui *bfs_order, ui *candidates_count) {
    ui query_vertices_num = query_graph->getVerticesCount();
    VertexID root_vertex = bfs_order[0];
    order = new ui[query_vertices_num];
    pivot = new ui[query_vertices_num];
    std::vector<bool> visited_vertices(query_vertices_num, false);

    std::vector<std::vector<ui>> core_paths;
    std::vector<std::vector<std::vector<ui>>> forests;
    std::vector<ui> leaves;

    // collect vertices with degree 1.
    generateLeaves(query_graph, leaves);
    // If a core structure exist,
    if (query_graph->getCoreValue(root_vertex) > 1) {
        // generate core paths.
        std::vector<ui> temp_core_path;
        generateCorePaths(query_graph, tree, root_vertex, temp_core_path, core_paths);
        // generate tree paths. repeat on all vertices because tree paths can start at any vertex in core structure.
        for (ui i = 0; i < query_vertices_num; ++i) {
            VertexID cur_vertex = i;
            if (query_graph->getCoreValue(cur_vertex) > 1) {
                std::vector<std::vector<ui>> temp_tree_paths;
                std::vector<ui> temp_tree_path;
                generateTreePaths(query_graph, tree, cur_vertex, temp_tree_path, temp_tree_paths);
                if (!temp_tree_paths.empty()) {
                    forests.emplace_back(temp_tree_paths);
                }
            }
        }
    }
    else {
        std::vector<std::vector<ui>> temp_tree_paths;
        std::vector<ui> temp_tree_path;
        // generate tree paths from root vertex.
        generateTreePaths(query_graph, tree, root_vertex, temp_tree_path, temp_tree_paths);
        if (!temp_tree_paths.empty()) {
            forests.emplace_back(temp_tree_paths);
        }
    }

    // Order core paths.
    ui selected_vertices_count = 0;
    order[selected_vertices_count++] = root_vertex;
    visited_vertices[root_vertex] = true;

    if (!core_paths.empty()) {
        std::vector<std::vector<size_t>> paths_embededdings_num;
        std::vector<ui> paths_non_tree_edge_num;
        for (auto& path : core_paths) {
            ui non_tree_edge_num = generateNoneTreeEdgesCount(query_graph, tree, path);
            paths_non_tree_edge_num.push_back(non_tree_edge_num + 1); // plus one so that the denominator is at least 1.

            std::vector<size_t> path_embeddings_num;
            estimatePathEmbeddsingsNum(path, edge_matrix, path_embeddings_num);
            paths_embededdings_num.emplace_back(path_embeddings_num);
        }

        // Select the start path that has less embeddings and more non-tree edges.
        double min_value = std::numeric_limits<double>::max();
        ui selected_path_index = 0;

        for (ui i = 0; i < core_paths.size(); ++i) {
            double cur_value = paths_embededdings_num[i][0] / (double) paths_non_tree_edge_num[i];

            if (cur_value < min_value) {
                min_value = cur_value;
                selected_path_index = i;
            }
        }

        // load the selected path into order sequence.
        for (ui i = 1; i < core_paths[selected_path_index].size(); ++i) {
            order[selected_vertices_count] = core_paths[selected_path_index][i];
            pivot[selected_vertices_count] = core_paths[selected_path_index][i - 1];
            selected_vertices_count += 1;
            visited_vertices[core_paths[selected_path_index][i]] = true;
        }

        core_paths.erase(core_paths.begin() + selected_path_index);
        paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
        paths_non_tree_edge_num.erase(paths_non_tree_edge_num.begin() + selected_path_index);

        while (!core_paths.empty()) {
            min_value = std::numeric_limits<double>::max();
            selected_path_index = 0;

            // For other paths, 
            for (ui i = 0; i < core_paths.size(); ++i) {
                ui path_root_vertex_idx = 0;

                // first select start vertex of this path (the first unvisited vertex)
                for (ui j = 0; j < core_paths[i].size(); ++j) {
                    VertexID cur_vertex = core_paths[i][j];

                    if (visited_vertices[cur_vertex])
                        continue;

                    path_root_vertex_idx = j - 1;
                    break;
                }

                // estimate cardinality (|embeddings| / |candidate set|)
                double cur_value = paths_embededdings_num[i][path_root_vertex_idx] / (double)candidates_count[core_paths[i][path_root_vertex_idx]];
                if (cur_value < min_value) {
                    min_value = cur_value;
                    selected_path_index = i;
                }
            }

            for (ui i = 1; i < core_paths[selected_path_index].size(); ++i) {
                if (visited_vertices[core_paths[selected_path_index][i]])
                    continue;

                order[selected_vertices_count] = core_paths[selected_path_index][i];
                pivot[selected_vertices_count] = core_paths[selected_path_index][i - 1];
                selected_vertices_count += 1;
                visited_vertices[core_paths[selected_path_index][i]] = true;
            }

            core_paths.erase(core_paths.begin() + selected_path_index);
            paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
        }
    }

    // Order tree paths. Similar to core path.
    for (auto& tree_paths : forests) {
        std::vector<std::vector<size_t>> paths_embededdings_num;
        for (auto& path : tree_paths) {
            std::vector<size_t> path_embeddings_num;
            estimatePathEmbeddsingsNum(path, edge_matrix, path_embeddings_num);
            paths_embededdings_num.emplace_back(path_embeddings_num);
        }

        while (!tree_paths.empty()) {
            double min_value = std::numeric_limits<double>::max();
            ui selected_path_index = 0;

            for (ui i = 0; i < tree_paths.size(); ++i) {
                ui path_root_vertex_idx = 0;
                for (ui j = 0; j < tree_paths[i].size(); ++j) {
                    VertexID cur_vertex = tree_paths[i][j];

                    if (visited_vertices[cur_vertex])
                        continue;

                    path_root_vertex_idx = j == 0 ? j : j - 1;
                    break;
                }

                double cur_value = paths_embededdings_num[i][path_root_vertex_idx] / (double)candidates_count[tree_paths[i][path_root_vertex_idx]];
                if (cur_value < min_value) {
                    min_value = cur_value;
                    selected_path_index = i;
                }
            }

            for (ui i = 0; i < tree_paths[selected_path_index].size(); ++i) {
                if (visited_vertices[tree_paths[selected_path_index][i]])
                    continue;

                order[selected_vertices_count] = tree_paths[selected_path_index][i];
                pivot[selected_vertices_count] = tree_paths[selected_path_index][i - 1];
                selected_vertices_count += 1;
                visited_vertices[tree_paths[selected_path_index][i]] = true;
            }

            tree_paths.erase(tree_paths.begin() + selected_path_index);
            paths_embededdings_num.erase(paths_embededdings_num.begin() + selected_path_index);
        }
    }

    // Order the leaves.
    while (!leaves.empty()) {
        double min_value = std::numeric_limits<double>::max();
        ui selected_leaf_index = 0;

        // select the leaf that has the least candidates.
        for (ui i = 0; i < leaves.size(); ++i) {
            VertexID vertex = leaves[i];
            double cur_value = candidates_count[vertex];

            if (cur_value < min_value) {
                min_value = cur_value;
                selected_leaf_index = i;
            }
        }

        if (!visited_vertices[leaves[selected_leaf_index]]) {
            order[selected_vertices_count] = leaves[selected_leaf_index];
            pivot[selected_vertices_count] = tree[leaves[selected_leaf_index]].parent_;
            selected_vertices_count += 1;
            visited_vertices[leaves[selected_leaf_index]] = true;
        }
        leaves.erase(leaves.begin() + selected_leaf_index);
    }
}

void GenerateQueryPlan::estimatePathEmbeddsingsNum(std::vector<ui> &path, Edges ***edge_matrix,
                                                   std::vector<size_t> &estimated_embeddings_num) {
    assert(path.size() > 1);
    std::vector<size_t> parent;
    std::vector<size_t> children;

    estimated_embeddings_num.resize(path.size() - 1);
    Edges& last_edge = *edge_matrix[path[path.size() - 2]][path[path.size() - 1]];
    children.resize(last_edge.vertex_count_);

    // For the last edge in the path, the source vertex has last_edge.vertex_count_ instances.
    // For i-th instance of source vertex, the number of instances of last edge is children[i].
    size_t sum = 0;
    for (ui i = 0; i < last_edge.vertex_count_; ++i) {
        children[i] = last_edge.offset_[i + 1] - last_edge.offset_[i];
        sum += children[i];
    }

    estimated_embeddings_num[path.size() - 2] = sum;

    for (int i = path.size() - 2; i >= 1; --i) {
        ui begin = path[i - 1];
        ui end = path[i];

        Edges& edge = *edge_matrix[begin][end];
        parent.resize(edge.vertex_count_);

        sum = 0;
        // For each instances of parent vertex, 
        for (ui j = 0; j < edge.vertex_count_; ++j) {

            size_t local_sum = 0;
            // sum up the number of instances of child vertex.
            for (ui k = edge.offset_[j]; k < edge.offset_[j + 1]; ++k) {
                ui nbr = edge.edge_[k];
                local_sum += children[nbr];
            }

            parent[j] = local_sum;
            sum += local_sum;
        }

        estimated_embeddings_num[i - 1] = sum;
        // parent vertex becomes child vertex.
        parent.swap(children);
    }
}

ui GenerateQueryPlan::generateNoneTreeEdgesCount(const Graph *query_graph, TreeNode *tree_node, std::vector<ui> &path) {
    ui non_tree_edge_count = query_graph->getVertexDegree(path[0]) - tree_node[path[0]].children_count_;

    for (ui i = 1; i < path.size(); ++i) {
        VertexID vertex = path[i];
        // the number of non-tree edges is degree minus the number of children minus the number of parents (1).
        non_tree_edge_count += query_graph->getVertexDegree(vertex) - tree_node[vertex].children_count_ - 1;
    }

    return non_tree_edge_count;
}

void GenerateQueryPlan::generateCorePaths(const Graph *query_graph, TreeNode *tree_node, VertexID cur_vertex,
                                          std::vector<ui> &cur_core_path, std::vector<std::vector<ui>> &core_paths) {
    TreeNode& node = tree_node[cur_vertex];
    cur_core_path.push_back(cur_vertex);

    bool is_core_leaf = true;
    // find vertices in core structure.
    for (ui i = 0; i < node.children_count_; ++i) {
        VertexID child = node.children_[i];
        if (query_graph->getCoreValue(child) > 1) {
            generateCorePaths(query_graph, tree_node, child, cur_core_path, core_paths);
            is_core_leaf = false;
        }
    }

    // We reach the border of core structure.
    if (is_core_leaf) {
        core_paths.emplace_back(cur_core_path);
    }
    // backtrack to find another potential path.
    cur_core_path.pop_back();
}

void GenerateQueryPlan::generateTreePaths(const Graph *query_graph, TreeNode *tree_node, VertexID cur_vertex,
                                          std::vector<ui> &cur_tree_path, std::vector<std::vector<ui>> &tree_paths) {
    TreeNode& node = tree_node[cur_vertex];
    cur_tree_path.push_back(cur_vertex);

    bool is_tree_leaf = true;
    for (ui i = 0; i < node.children_count_; ++i) {
        VertexID child = node.children_[i];
        if (query_graph->getVertexDegree(child) > 1) {
            generateTreePaths(query_graph, tree_node, child, cur_tree_path, tree_paths);
            is_tree_leaf = false;
        }
    }

    if (is_tree_leaf && cur_tree_path.size() > 1) {
        tree_paths.emplace_back(cur_tree_path);
    }
    cur_tree_path.pop_back();
}

void GenerateQueryPlan::generateLeaves(const Graph *query_graph, std::vector<ui> &leaves) {
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID cur_vertex = i;
        if (query_graph->getVertexDegree(cur_vertex) == 1) {
            leaves.push_back(cur_vertex);
        }
    }
}

void GenerateQueryPlan::checkQueryPlanCorrectness(const Graph *query_graph, ui *order, ui *pivot) {
    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<bool> visited_vertices(query_vertices_num, false);
    // Check whether each query vertex is in the order.
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];
        assert(vertex < query_vertices_num && vertex >= 0);

        visited_vertices[vertex] = true;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID vertex = i;
        assert(visited_vertices[vertex]);
    }

    // Check whether the order is connected.
    std::fill(visited_vertices.begin(), visited_vertices.end(), false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];
        VertexID pivot_vertex = pivot[i];
        assert(query_graph->checkEdgeExistence(vertex, pivot_vertex));
        assert(visited_vertices[pivot_vertex]);
        visited_vertices[vertex] = true;
    }
}

void GenerateQueryPlan::printSimplifiedQueryPlan(const Graph *query_graph, ui *order) {
    ui query_vertices_num = query_graph->getVerticesCount();
    printf("Query Plan: ");
    for (ui i = 0; i < query_vertices_num; ++i) {
        printf("%u ", order[i]);
    }
    printf("\n");
}