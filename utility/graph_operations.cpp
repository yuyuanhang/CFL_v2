#include <queue>
#include "graph_operations.h"

void GraphOperations::getKCore(const Graph *graph, int *core_table) {
    int vertices_count = graph->getVerticesCount();
    int max_in_degree = graph->getGraphInMaxDegree();
    int max_out_degree = graph->getGraphOutMaxDegree();
    int max_degree = max_in_degree + max_out_degree;

    int* vertices = new int[vertices_count];          // Vertices sorted by degree.
    int* position = new int[vertices_count];          // The position of vertices in vertices array.
    int* degree_bin = new int[max_degree + 1];      // Record the number of vertices the degree of that range from 0 to max_degree.
    int* offset = new int[max_degree + 1];          // The offset in vertices array according to degree.

    std::fill(degree_bin, degree_bin + (max_degree + 1), 0);

    // Record the degree of vertices.
    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexInDegree(i);
        degree += graph->getVertexOutDegree(i);
        core_table[i] = degree;
        degree_bin[degree] += 1;
    }

    // offset[i] = offset[i - 1] + degree_bin[i - 1]; offset[0] = 0;
    int start = 0;
    for (int i = 0; i < max_degree + 1; ++i) {
        offset[i] = start;
        start += degree_bin[i];
    }

    // Here, offset[i] represents the sum of the number of vertices the degree of that is less than i.
    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexInDegree(i);
        degree += graph->getVertexOutDegree(i);
        position[i] = offset[degree];
        vertices[position[i]] = i;
        offset[degree] += 1;
    }

    // Then, offset[i] represents the sum of the number of vertices the degree of that is no greater than i. Hence, we recover it.
    for (int i = max_degree; i > 0; --i) {
        offset[i] = offset[i - 1];
    }
    offset[0] = 0;

    // For each vertex ordered by degree (ascending), 
    for (int i = 0; i < vertices_count; ++i) {
        int v = vertices[i];

        ui in_count;
        const VertexID * in_neighbors = graph->getVertexInNeighbors(v, in_count);

        // we check its all incoming neighbors. 
        for(ui j = 0; j < in_count; ++j) {
            int u = in_neighbors[j];

            // If its neighbor has larger degree, we need to process it; Otherwise, we have processed it.
            // After peeling, the degree of vertices in core_table is the same.
            if (core_table[u] > core_table[v]) {

                // Get the position and vertex which is with the same degree
                // and at the start position of vertices array.
                int cur_degree_u = core_table[u];
                int position_u = position[u];
                int position_w = offset[cur_degree_u];
                int w = vertices[position_w];

                if (u != w) {
                    // Swap u and w.
                    position[u] = position_w;
                    position[w] = position_u;
                    vertices[position_u] = w;
                    vertices[position_w] = u;
                }

                offset[cur_degree_u] += 1;
                core_table[u] -= 1;
            }
        }

        ui out_count;
        const VertexID * out_neighbors = graph->getVertexOutNeighbors(v, out_count);

        // we check its all outgoing neighbors. 
        for(ui j = 0; j < out_count; ++j) {
            int u = out_neighbors[j];

            // If its neighbor has larger degree, we need to process it; Otherwise, we have processed it.
            // After peeling, the degree of vertices in core_table is the same.
            if (core_table[u] > core_table[v]) {

                // Get the position and vertex which is with the same degree
                // and at the start position of vertices array.
                int cur_degree_u = core_table[u];
                int position_u = position[u];
                int position_w = offset[cur_degree_u];
                int w = vertices[position_w];

                if (u != w) {
                    // Swap u and w.
                    position[u] = position_w;
                    position[w] = position_u;
                    vertices[position_u] = w;
                    vertices[position_w] = u;
                }

                offset[cur_degree_u] += 1;
                core_table[u] -= 1;
            }
        }
    }

    delete[] vertices;
    delete[] position;
    delete[] degree_bin;
    delete[] offset;
}

void GraphOperations::bfsTraversal(const Graph *graph, VertexID root_vertex, TreeNode *&tree, VertexID *&bfs_order) {
    ui vertex_num = graph->getVerticesCount();

    std::queue<VertexID> bfs_queue;
    std::vector<bool> visited(vertex_num, false);

    tree = new TreeNode[vertex_num];
    for (ui i = 0; i < vertex_num; ++i) {
        tree[i].initialize(vertex_num);
    }
    bfs_order = new VertexID[vertex_num];

    ui visited_vertex_count = 0;
    // visit root vertex.
    bfs_queue.push(root_vertex);
    visited[root_vertex] = true;
    tree[root_vertex].level_ = 0;
    tree[root_vertex].id_ = root_vertex;

    while(!bfs_queue.empty()) {
        const VertexID u = bfs_queue.front();
        bfs_queue.pop();
        bfs_order[visited_vertex_count++] = u;

        ui u_in_nbrs_count;
        ui u_out_nbrs_count;
        // get neighbors of current vertex.
        const VertexID* u_in_nbrs = graph->getVertexInNeighbors(u, u_in_nbrs_count);
        const VertexID* u_out_nbrs = graph->getVertexOutNeighbors(u, u_out_nbrs_count);

        for (ui i = 0; i < u_in_nbrs_count; ++i) {
            VertexID u_in_nbr = u_in_nbrs[i];

            if (!visited[u_in_nbr]) {
                // record this neighbor.
                bfs_queue.push(u_in_nbr);
                visited[u_in_nbr] = true;
                tree[u_in_nbr].id_ = u_in_nbr;
                // make it a child.
                tree[u_in_nbr].parent_ = u;
                tree[u_in_nbr].level_ = tree[u].level_ + 1;
                tree[u_in_nbr].p_direction_ = OUT;
                tree[u].children_[tree[u].children_count_] = u_in_nbr;
                tree[u].c_directions_[tree[u].children_count_] = IN;
                tree[u].children_count_++;
            }
        }

        for (ui i = 0; i < u_out_nbrs_count; ++i) {
            VertexID u_out_nbr = u_out_nbrs[i];

            if (!visited[u_out_nbr]) {
                // record this neighbor.
                bfs_queue.push(u_out_nbr);
                visited[u_out_nbr] = true;
                tree[u_out_nbr].id_ = u_out_nbr;
                // make it a child.
                tree[u_out_nbr].parent_ = u;
                tree[u_out_nbr].level_ = tree[u].level_ + 1;
                tree[u_out_nbr].p_direction_ = IN;
                tree[u].children_[tree[u].children_count_] = u_out_nbr;
                tree[u].c_directions_[tree[u].children_count_] = OUT;
                tree[u].children_count_++;
            }
        }
    }
}