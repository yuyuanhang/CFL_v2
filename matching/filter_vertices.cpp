#include "filter_vertices.h"
#define INVALID_VERTEX_ID 100000000

bool
FilterVertices::CFLFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                          ui *&order, TreeNode *&tree) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);
    int level_count;
    ui* level_offset;
    GenerateFilteringPlan::generateCFLFilterPlan(data_graph, query_graph, tree, order, level_count, level_offset);

    VertexID start_vertex = order[0];
    // compute candidate set (with NLF filter) for start vertex.
    computeCandidateWithNLF(data_graph, query_graph, start_vertex, candidates_count[start_vertex], candidates[start_vertex]);

    ui* updated_flag = new ui[data_graph->getVerticesCount()];
    ui* flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);

    // Top-down generation.
    for (int i = 1; i < level_count; ++i) {
        // generate candidate sets for vertices on current level based on order.
        for (ui j = level_offset[i]; j < level_offset[i + 1]; ++j) {
            VertexID query_vertex = order[j];
            TreeNode& node = tree[query_vertex];

            generateCandidates(data_graph, query_graph, query_vertex, node.bn_, node.bn_directions_, node.bn_count_, candidates, candidates_count, flag, updated_flag);
        }

        // prune invalid candidates for vertices on current level based on reverse order.
        for (ui j = level_offset[i + 1] - 1; j >= level_offset[i]; --j) {
            VertexID query_vertex = order[j];
            TreeNode& node = tree[query_vertex];

            if (node.fn_count_ > 0) {
                pruneCandidates(data_graph, query_graph, query_vertex, node.fn_, node.fn_directions_, node.fn_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
    }

    // Bottom-up refinement.
    for (int i = level_count - 2; i >= 0; --i) {
        for (ui j = level_offset[i]; j < level_offset[i + 1]; ++j) {
            VertexID query_vertex = order[j];
            TreeNode& node = tree[query_vertex];

            if (node.under_level_count_ > 0) {
                pruneCandidates(data_graph, query_graph, query_vertex, node.under_level_, node.under_level_directions_, node.under_level_count_, candidates, candidates_count, flag, updated_flag);
            }
        }
    }

    compactCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    delete[] updated_flag;
    delete[] flag;
    return isCandidateSetValid(candidates, candidates_count, query_graph->getVerticesCount());
}

void FilterVertices::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui **&candidates,
                                    ui *&candidates_count) {
    ui query_vertex_num = query_graph->getVerticesCount();
    ui candidates_max_num = data_graph->getGraphMaxLabelFrequency();

    candidates_count = new ui[query_vertex_num];
    memset(candidates_count, 0, sizeof(ui) * query_vertex_num);

    candidates = new ui*[query_vertex_num];

    for (ui i = 0; i < query_vertex_num; ++i) {
        candidates[i] = new ui[candidates_max_num];
    }
}

void
FilterVertices::computeCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                               ui &count, ui *buffer) {
    LabelID label = query_graph->getVertexLabel(query_vertex);
    ui in_degree = query_graph->getVertexInDegree(query_vertex);
    ui out_degree = query_graph->getVertexOutDegree(query_vertex);
#if OPTIMIZED_LABELED_GRAPH == 1
    const std::unordered_map<LabelID, ui>* query_vertex_in_nlf = query_graph->getVertexInNLF(query_vertex);
    const std::unordered_map<LabelID, ui>* query_vertex_out_nlf = query_graph->getVertexOutNLF(query_vertex);
#endif
    // get vertices with the same label.
    ui data_vertex_num;
    const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);
    count = 0;

    // For each vertex with the same label,
    for (ui j = 0; j < data_vertex_num; ++j) {
        ui data_vertex = data_vertices[j];

        // check if it has larger degree than target vertex in query graph.
        if (data_graph->getVertexInDegree(data_vertex) >= in_degree && data_graph->getVertexOutDegree(data_vertex) >= out_degree) {

            // check if it has more neighbors for each label.
#if OPTIMIZED_LABELED_GRAPH == 1
            const std::unordered_map<LabelID, ui>* data_vertex_in_nlf = data_graph->getVertexInNLF(data_vertex);
            const std::unordered_map<LabelID, ui>* data_vertex_out_nlf = data_graph->getVertexOutNLF(data_vertex);

            if (data_vertex_in_nlf->size() >= query_vertex_in_nlf->size() && data_vertex_out_nlf->size() >= query_vertex_out_nlf->size()) {
                bool is_valid = true;

                for (auto element : *query_vertex_in_nlf) {
                    auto iter = data_vertex_in_nlf->find(element.first);
                    if (iter == data_vertex_in_nlf->end() || iter->second < element.second) {
                        is_valid = false;
                        break;
                    }
                }

                if (is_valid) {
                    for (auto element : *query_vertex_out_nlf) {
                        auto iter = data_vertex_out_nlf->find(element.first);
                        if (iter == data_vertex_out_nlf->end() || iter->second < element.second) {
                            is_valid = false;
                            break;
                        }
                    }
                }

                if (is_valid) {
                    if (buffer != NULL) {
                        buffer[count] = data_vertex;
                    }
                    count += 1;
                }
            }
#else
            // if buffer is not null, data vertex need to be recorded.
            if (buffer != NULL) {
                buffer[count] = data_vertex;
            }
            count += 1;
#endif
        }
    }

}

void FilterVertices::generateCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                       VertexID *pivot_vertices, bool *pivot_directions, ui pivot_vertices_count, VertexID **candidates,
                                       ui *candidates_count, ui *flag, ui *updated_flag) {
    LabelID query_vertex_label = query_graph->getVertexLabel(query_vertex);
    ui query_vertex_in_degree = query_graph->getVertexInDegree(query_vertex);
    ui query_vertex_out_degree = query_graph->getVertexOutDegree(query_vertex);
#if OPTIMIZED_LABELED_GRAPH == 1
    const std::unordered_map<LabelID , ui>* query_vertex_in_nlf = query_graph->getVertexInNLF(query_vertex);
    const std::unordered_map<LabelID , ui>* query_vertex_out_nlf = query_graph->getVertexOutNLF(query_vertex);
#endif
    ui count = 0;
    ui updated_flag_count = 0;
    // For each neighbor in the front of current vertex, 
    for (ui i = 0; i < pivot_vertices_count; ++i) {
        VertexID pivot_vertex = pivot_vertices[i];
        bool pivot_direction = pivot_directions[i];

        // For each candidate vertex in the candidate set of current neighbor, 
        for (ui j = 0; j < candidates_count[pivot_vertex]; ++j) {
            VertexID v = candidates[pivot_vertex][j];

            if (v == INVALID_VERTEX_ID)
                continue;

            if (pivot_direction == OUT) {
                ui v_in_nbrs_count;
                const VertexID* v_in_nbrs = data_graph->getVertexInNeighbors(v, v_in_nbrs_count);

                for (ui k = 0; k < v_in_nbrs_count; ++k) {
                    VertexID v_in_nbr = v_in_nbrs[k];
                    LabelID v_in_nbr_label = data_graph->getVertexLabel(v_in_nbr);
                    ui v_in_nbr_in_degree = data_graph->getVertexInDegree(v_in_nbr);
                    ui v_in_nbr_out_degree = data_graph->getVertexOutDegree(v_in_nbr);
    
                    // record its neighbor in data graph that can pass label and degree filter.
                    if (flag[v_in_nbr] == count && v_in_nbr_label == query_vertex_label && v_in_nbr_in_degree >= query_vertex_in_degree && v_in_nbr_out_degree >= query_vertex_out_degree) {
                        flag[v_in_nbr] += 1;
    
                        if (count == 0) {
                            updated_flag[updated_flag_count++] = v_in_nbr;
                        }
                    }
                }
            } else if (pivot_direction == IN) {
                ui v_out_nbrs_count;
                const VertexID* v_out_nbrs = data_graph->getVertexOutNeighbors(v, v_out_nbrs_count);

                for (ui k = 0; k < v_out_nbrs_count; ++k) {
                    VertexID v_out_nbr = v_out_nbrs[k];
                    LabelID v_out_nbr_label = data_graph->getVertexLabel(v_out_nbr);
                    ui v_out_nbr_in_degree = data_graph->getVertexInDegree(v_out_nbr);
                    ui v_out_nbr_out_degree = data_graph->getVertexOutDegree(v_out_nbr);
    
                    // record its neighbor in data graph that can pass label and degree filter.
                    if (flag[v_out_nbr] == count && v_out_nbr_label == query_vertex_label && v_out_nbr_in_degree >= query_vertex_in_degree && v_out_nbr_out_degree >= query_vertex_out_degree) {
                        flag[v_out_nbr] += 1;
    
                        if (count == 0) {
                            updated_flag[updated_flag_count++] = v_out_nbr;
                        }
                    }
                }
            }
        }

        count += 1;
    }

    // For each potential candidate vertex,
    for (ui i = 0; i < updated_flag_count; ++i) {
        VertexID v = updated_flag[i];
        // If it has neighbors in each candidate set of its neighbors in query graph, 
        if (flag[v] == count) {
            // NLF filter.
#if OPTIMIZED_LABELED_GRAPH == 1
            const std::unordered_map<LabelID, ui>* data_vertex_in_nlf = data_graph->getVertexInNLF(v);
            const std::unordered_map<LabelID, ui>* data_vertex_out_nlf = data_graph->getVertexOutNLF(v);

            if (data_vertex_in_nlf->size() >= query_vertex_in_nlf->size() && data_vertex_out_nlf->size() >= query_vertex_out_nlf->size()) {
                bool is_valid = true;

                for (auto element : *query_vertex_in_nlf) {
                    auto iter = data_vertex_in_nlf->find(element.first);
                    if (iter == data_vertex_in_nlf->end() || iter->second < element.second) {
                        is_valid = false;
                        break;
                    }
                }

                if (is_valid) {
                    for (auto element : *query_vertex_out_nlf) {
                        auto iter = data_vertex_out_nlf->find(element.first);
                        if (iter == data_vertex_out_nlf->end() || iter->second < element.second) {
                            is_valid = false;
                            break;
                        }
                    }
                }

                // and pass NLF filter, add it in candidate set.
                if (is_valid) {
                    candidates[query_vertex][candidates_count[query_vertex]++] = v;
                }
            }
#else
            candidates[query_vertex][candidates_count[query_vertex]++] = v;
#endif
        }
    }

    // recover array for next query vertex.
    for (ui i = 0; i < updated_flag_count; ++i) {
        ui v = updated_flag[i];
        flag[v] = 0;
    }
}

void FilterVertices::pruneCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                    VertexID *pivot_vertices, bool *pivot_directions, ui pivot_vertices_count, VertexID **candidates,
                                    ui *candidates_count, ui *flag, ui *updated_flag) {
    LabelID query_vertex_label = query_graph->getVertexLabel(query_vertex);
    ui query_vertex_in_degree = query_graph->getVertexInDegree(query_vertex);
    ui query_vertex_out_degree = query_graph->getVertexOutDegree(query_vertex);

    ui count = 0;
    ui updated_flag_count = 0;
    // For each neighbor that is behind current vertex and on the same level.
    for (ui i = 0; i < pivot_vertices_count; ++i) {
        VertexID pivot_vertex = pivot_vertices[i];
        bool pivot_direction = pivot_directions[i];

        // For each candidate vertex in the candidate set of current neighbor,
        for (ui j = 0; j < candidates_count[pivot_vertex]; ++j) {
            VertexID v = candidates[pivot_vertex][j];

            if (v == INVALID_VERTEX_ID)
                continue;

            if (pivot_direction == OUT) {
                ui v_in_nbrs_count;
                const VertexID* v_in_nbrs = data_graph->getVertexInNeighbors(v, v_in_nbrs_count);
    
                // record intersection set.
                for (ui k = 0; k < v_in_nbrs_count; ++k) {
                    VertexID v_in_nbr = v_in_nbrs[k];
                    LabelID v_in_nbr_label = data_graph->getVertexLabel(v_in_nbr);
                    ui v_in_nbr_in_degree = data_graph->getVertexInDegree(v_in_nbr);
                    ui v_in_nbr_out_degree = data_graph->getVertexOutDegree(v_in_nbr);
    
                    if (flag[v_in_nbr] == count && v_in_nbr_label == query_vertex_label && v_in_nbr_in_degree >= query_vertex_in_degree && v_in_nbr_out_degree >= query_vertex_out_degree) {
                        flag[v_in_nbr] += 1;
    
                        if (count == 0) {
                            updated_flag[updated_flag_count++] = v_in_nbr;
                        }
                    }
                }
            } else if (pivot_direction == IN) {
                ui v_out_nbrs_count;
                const VertexID* v_out_nbrs = data_graph->getVertexOutNeighbors(v, v_out_nbrs_count);
    
                // record intersection set.
                for (ui k = 0; k < v_out_nbrs_count; ++k) {
                    VertexID v_out_nbr = v_out_nbrs[k];
                    LabelID v_out_nbr_label = data_graph->getVertexLabel(v_out_nbr);
                    ui v_out_nbr_in_degree = data_graph->getVertexInDegree(v_out_nbr);
                    ui v_out_nbr_out_degree = data_graph->getVertexOutDegree(v_out_nbr);
    
                    if (flag[v_out_nbr] == count && v_out_nbr_label == query_vertex_label && v_out_nbr_in_degree >= query_vertex_in_degree && v_out_nbr_out_degree >= query_vertex_out_degree) {
                        flag[v_out_nbr] += 1;
    
                        if (count == 0) {
                            updated_flag[updated_flag_count++] = v_out_nbr;
                        }
                    }
                }
            }
        }

        count += 1;
    }

    for (ui i = 0; i < candidates_count[query_vertex]; ++i) {
        ui v = candidates[query_vertex][i];
        if (v == INVALID_VERTEX_ID)
            continue;
        // If current candidate vertex is not in the intersection set, prune it.
        if (flag[v] != count) {
            candidates[query_vertex][i] = INVALID_VERTEX_ID;
        }
    }

    // recover array for next query vertex.
    for (ui i = 0; i < updated_flag_count; ++i) {
        ui v = updated_flag[i];
        flag[v] = 0;
    }
}

void FilterVertices::compactCandidates(ui **&candidates, ui *&candidates_count, ui query_vertex_num) {
    for (ui i = 0; i < query_vertex_num; ++i) {
        VertexID query_vertex = i;
        ui next_position = 0;
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];

            if (data_vertex != INVALID_VERTEX_ID) {
                candidates[query_vertex][next_position++] = data_vertex;
            }
        }

        candidates_count[query_vertex] = next_position;
    }
}

bool FilterVertices::isCandidateSetValid(ui **&candidates, ui *&candidates_count, ui query_vertex_num) {
    for (ui i = 0; i < query_vertex_num; ++i) {
        if (candidates_count[i] == 0)
            return false;
    }
    return true;
}

void FilterVertices::sortCandidates(ui **candidates, ui *candidates_count, ui num) {
    for (ui i = 0; i < num; ++i) {
        std::sort(candidates[i], candidates[i] + candidates_count[i]);
    }
}