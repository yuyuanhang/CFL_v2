#include "build_table.h"

void BuildTable::buildTables(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                             Edges ***edge_matrix) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ui* flag = new ui[data_graph->getVerticesCount()];
    ui* updated_flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);

    for (ui i = 0; i < query_vertices_num; ++i) {
        for (ui j = 0; j < query_vertices_num; ++j) {
            edge_matrix[i][j] = NULL;
        }
    }

    std::vector<VertexID> build_table_order(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        build_table_order[i] = i;
    }

    std::sort(build_table_order.begin(), build_table_order.end(), [query_graph](VertexID l, VertexID r) {
        if (query_graph->getVertexDegree(l) == query_graph->getVertexDegree(r)) {
            return l < r;
        }
        return query_graph->getVertexDegree(l) > query_graph->getVertexDegree(r);
    });

    std::vector<ui> temp_edges(data_graph->getEdgesCount());

    for (auto u : build_table_order) {
        ui u_in_nbrs_count;
        const VertexID* u_in_nbrs = query_graph->getVertexInNeighbors(u, u_in_nbrs_count);

        ui updated_flag_count = 0;

        for (ui i = 0; i < u_in_nbrs_count; ++i) {
            VertexID u_in_nbr = u_in_nbrs[i];

            if (edge_matrix[u][u_in_nbr] != NULL)
                continue;

            if (updated_flag_count == 0) {
                for (ui j = 0; j < candidates_count[u]; ++j) {
                    VertexID v = candidates[u][j];
                    flag[v] = j + 1;
                    updated_flag[updated_flag_count++] = v;
                }
            }

            edge_matrix[u_in_nbr][u] = new Edges;
            edge_matrix[u_in_nbr][u]->vertex_count_ = candidates_count[u_in_nbr];
            edge_matrix[u_in_nbr][u]->offset_ = new ui[candidates_count[u_in_nbr] + 1];
            edge_matrix[u_in_nbr][u]->direction_ = OUT;

            edge_matrix[u][u_in_nbr] = new Edges;
            edge_matrix[u][u_in_nbr]->vertex_count_ = candidates_count[u];
            edge_matrix[u][u_in_nbr]->offset_ = new ui[candidates_count[u] + 1];
            edge_matrix[u][u_in_nbr]->direction_ = IN;
            std::fill(edge_matrix[u][u_in_nbr]->offset_, edge_matrix[u][u_in_nbr]->offset_ + candidates_count[u] + 1, 0);

            ui local_edge_count = 0;
            ui local_max_degree = 0;

            // For each candidate vertex of neighbor, 
            for (ui j = 0; j < candidates_count[u_in_nbr]; ++j) {
                VertexID v = candidates[u_in_nbr][j];
                edge_matrix[u_in_nbr][u]->offset_[j] = local_edge_count;

                ui v_out_nbrs_count;
                const VertexID* v_out_nbrs = data_graph->getVertexOutNeighbors(v, v_out_nbrs_count);

                ui local_degree = 0;

                // check if its neighbor is in candidate set of current neighbor.
                for (ui k = 0; k < v_out_nbrs_count; ++k) {
                    VertexID v_out_nbr = v_out_nbrs[k];

                    if (flag[v_out_nbr] != 0) {
                        // if it is, revert its position in candidate set.
                        ui position = flag[v_out_nbr] - 1;
                        // record this edge for this neighbor.
                        temp_edges[local_edge_count++] = position;
                        // accumulate this edge for current vertex.
                        edge_matrix[u][u_in_nbr]->offset_[position + 1] += 1;
                        // count the number of edge for this neighbor.
                        local_degree += 1;
                    }
                }

                if (local_degree > local_max_degree) {
                    local_max_degree = local_degree;
                }
            }

            edge_matrix[u_in_nbr][u]->offset_[candidates_count[u_in_nbr]] = local_edge_count;
            edge_matrix[u_in_nbr][u]->max_degree_ = local_max_degree;
            edge_matrix[u_in_nbr][u]->edge_count_ = local_edge_count;
            edge_matrix[u_in_nbr][u]->edge_ = new ui[local_edge_count];
            // copy the list recorded above
            std::copy(temp_edges.begin(), temp_edges.begin() + local_edge_count, edge_matrix[u_in_nbr][u]->edge_);

            edge_matrix[u][u_in_nbr]->edge_count_ = local_edge_count;
            edge_matrix[u][u_in_nbr]->edge_ = new ui[local_edge_count];

            // count statistics
            local_max_degree = 0;
            for (ui j = 1; j <= candidates_count[u]; ++j) {
                if (edge_matrix[u][u_in_nbr]->offset_[j] > local_max_degree) {
                    local_max_degree = edge_matrix[u][u_in_nbr]->offset_[j];
                }
                edge_matrix[u][u_in_nbr]->offset_[j] += edge_matrix[u][u_in_nbr]->offset_[j - 1];
            }

            edge_matrix[u][u_in_nbr]->max_degree_ = local_max_degree;

            // For each candidate vertex of neighbor, check if they are connected; and create inverted edge list for vertex u.
            for (ui j = 0; j < candidates_count[u_in_nbr]; ++j) {
                ui begin = j;
                for (ui k = edge_matrix[u_in_nbr][u]->offset_[begin]; k < edge_matrix[u_in_nbr][u]->offset_[begin + 1]; ++k) {
                    ui end = edge_matrix[u_in_nbr][u]->edge_[k];

                    edge_matrix[u][u_in_nbr]->edge_[edge_matrix[u][u_in_nbr]->offset_[end]++] = begin;
                }
            }

            // revert offset array.
            for (ui j = candidates_count[u]; j >= 1; --j) {
                edge_matrix[u][u_in_nbr]->offset_[j] = edge_matrix[u][u_in_nbr]->offset_[j - 1];
            }
            edge_matrix[u][u_in_nbr]->offset_[0] = 0;
        }

        for (ui i = 0; i < updated_flag_count; ++i) {
            VertexID v = updated_flag[i];
            flag[v] = 0;
        }

        ui u_out_nbrs_count;
        const VertexID* u_out_nbrs = query_graph->getVertexOutNeighbors(u, u_out_nbrs_count);

        updated_flag_count = 0;

        for (ui i = 0; i < u_out_nbrs_count; ++i) {
            VertexID u_out_nbr = u_out_nbrs[i];

            if (edge_matrix[u][u_out_nbr] != NULL)
                continue;

            if (updated_flag_count == 0) {
                for (ui j = 0; j < candidates_count[u]; ++j) {
                    VertexID v = candidates[u][j];
                    flag[v] = j + 1;
                    updated_flag[updated_flag_count++] = v;
                }
            }

            edge_matrix[u_out_nbr][u] = new Edges;
            edge_matrix[u_out_nbr][u]->vertex_count_ = candidates_count[u_out_nbr];
            edge_matrix[u_out_nbr][u]->offset_ = new ui[candidates_count[u_out_nbr] + 1];
            edge_matrix[u_out_nbr][u]->direction_ = IN;

            edge_matrix[u][u_out_nbr] = new Edges;
            edge_matrix[u][u_out_nbr]->vertex_count_ = candidates_count[u];
            edge_matrix[u][u_out_nbr]->offset_ = new ui[candidates_count[u] + 1];
            edge_matrix[u][u_out_nbr]->direction_ = OUT;
            std::fill(edge_matrix[u][u_out_nbr]->offset_, edge_matrix[u][u_out_nbr]->offset_ + candidates_count[u] + 1, 0);

            ui local_edge_count = 0;
            ui local_max_degree = 0;

            // For each candidate vertex of neighbor, 
            for (ui j = 0; j < candidates_count[u_out_nbr]; ++j) {
                VertexID v = candidates[u_out_nbr][j];
                edge_matrix[u_out_nbr][u]->offset_[j] = local_edge_count;

                ui v_in_nbrs_count;
                const VertexID* v_in_nbrs = data_graph->getVertexInNeighbors(v, v_in_nbrs_count);

                ui local_degree = 0;

                // check if its neighbor is in candidate set of current neighbor.
                for (ui k = 0; k < v_in_nbrs_count; ++k) {
                    VertexID v_in_nbr = v_in_nbrs[k];

                    if (flag[v_in_nbr] != 0) {
                        // if it is, revert its position in candidate set.
                        ui position = flag[v_in_nbr] - 1;
                        // record this edge for this neighbor.
                        temp_edges[local_edge_count++] = position;
                        // accumulate this edge for current vertex.
                        edge_matrix[u][u_out_nbr]->offset_[position + 1] += 1;
                        // count the number of edge for this neighbor.
                        local_degree += 1;
                    }
                }

                if (local_degree > local_max_degree) {
                    local_max_degree = local_degree;
                }
            }

            edge_matrix[u_out_nbr][u]->offset_[candidates_count[u_out_nbr]] = local_edge_count;
            edge_matrix[u_out_nbr][u]->max_degree_ = local_max_degree;
            edge_matrix[u_out_nbr][u]->edge_count_ = local_edge_count;
            edge_matrix[u_out_nbr][u]->edge_ = new ui[local_edge_count];
            // copy the list recorded above
            std::copy(temp_edges.begin(), temp_edges.begin() + local_edge_count, edge_matrix[u_out_nbr][u]->edge_);

            edge_matrix[u][u_out_nbr]->edge_count_ = local_edge_count;
            edge_matrix[u][u_out_nbr]->edge_ = new ui[local_edge_count];

            // count statistics
            local_max_degree = 0;
            for (ui j = 1; j <= candidates_count[u]; ++j) {
                if (edge_matrix[u][u_out_nbr]->offset_[j] > local_max_degree) {
                    local_max_degree = edge_matrix[u][u_out_nbr]->offset_[j];
                }
                edge_matrix[u][u_out_nbr]->offset_[j] += edge_matrix[u][u_out_nbr]->offset_[j - 1];
            }

            edge_matrix[u][u_out_nbr]->max_degree_ = local_max_degree;

            // For each candidate vertex of neighbor, check if they are connected; and create inverted edge list for vertex u.
            for (ui j = 0; j < candidates_count[u_out_nbr]; ++j) {
                ui begin = j;
                for (ui k = edge_matrix[u_out_nbr][u]->offset_[begin]; k < edge_matrix[u_out_nbr][u]->offset_[begin + 1]; ++k) {
                    ui end = edge_matrix[u_out_nbr][u]->edge_[k];

                    edge_matrix[u][u_out_nbr]->edge_[edge_matrix[u][u_out_nbr]->offset_[end]++] = begin;
                }
            }

            // revert offset array.
            for (ui j = candidates_count[u]; j >= 1; --j) {
                edge_matrix[u][u_out_nbr]->offset_[j] = edge_matrix[u][u_out_nbr]->offset_[j - 1];
            }
            edge_matrix[u][u_out_nbr]->offset_[0] = 0;
        }

        for (ui i = 0; i < updated_flag_count; ++i) {
            VertexID v = updated_flag[i];
            flag[v] = 0;
        }
    }
}

void BuildTable::printTableCardinality(const Graph *query_graph, Edges ***edge_matrix) {
    std::vector<std::pair<std::pair<VertexID, VertexID >, ui>> core_edges;
    std::vector<std::pair<std::pair<VertexID, VertexID >, ui>> tree_edges;
    std::vector<std::pair<std::pair<VertexID, VertexID >, ui>> leaf_edges;

    ui query_vertices_num = query_graph->getVerticesCount();

    double sum = 0;
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID begin_vertex = i;

        for (ui j = i + 1; j < query_vertices_num; ++j) {
            VertexID end_vertex = j;

            if (query_graph->checkEdgeExistence(begin_vertex, end_vertex)) {
                ui cardinality = (*edge_matrix[begin_vertex][end_vertex]).edge_count_;
                sum += cardinality;
                // this edge belongs to core structure.
                if (query_graph->getCoreValue(begin_vertex) > 1 && query_graph->getCoreValue(end_vertex) > 1) {
                    core_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
                // this edge belongs to leaf structure because the degree of one vertex is 1.
                else if (query_graph->getVertexDegree(begin_vertex) == 1 || query_graph->getVertexDegree(end_vertex) == 1) {
                    leaf_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
                // Otherwise, it belongs to forest structure.
                else {
                    tree_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
            }
        }
    }

    printf("Index Info: CoreTable(%zu), TreeTable(%zu), LeafTable(%zu)\n", core_edges.size(), tree_edges.size(), leaf_edges.size());

    for (auto table_info : core_edges) {
        printf("CoreTable %u-%u: %u\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    for (auto table_info : tree_edges) {
        printf("TreeTable %u-%u: %u\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    for (auto table_info : leaf_edges) {
        printf("LeafTable %u-%u: %d\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    printf("Total Cardinality: %.1lf\n", sum);
}

size_t BuildTable::computeMemoryCostInBytes(const Graph *query_graph, ui *candidates_count, Edges ***edge_matrix) {
    size_t memory_cost_in_bytes = 0;
    size_t per_element_size = sizeof(ui);

    ui query_vertices_num = query_graph->getVerticesCount();
    for (ui i = 0; i < query_vertices_num; ++i) {
        memory_cost_in_bytes += candidates_count[i] * per_element_size;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID begin_vertex = i;
        for (ui j = i + 1; j < query_vertices_num; ++j) {
            VertexID end_vertex = j;

            if (query_graph->checkEdgeExistence(begin_vertex, end_vertex)) {
                // edge list and offset array.
                Edges& edge = *edge_matrix[begin_vertex][end_vertex];
                memory_cost_in_bytes += edge.edge_count_ * per_element_size + (edge.vertex_count_ + 1) * per_element_size;

                Edges& reverse_edge = *edge_matrix[end_vertex][begin_vertex];
                memory_cost_in_bytes += reverse_edge.edge_count_ * per_element_size + (reverse_edge.vertex_count_ + 1) * per_element_size;
            }

        }
    }

    return memory_cost_in_bytes;
}