#include "graph.h"
#include "utility/graph_operations.h"

Graph::Graph(){
    
    vertices_count_ = 0;
    in_edges_count_ = 0;
    out_edges_count_ = 0;
    labels_count_ = 0;
    max_in_degree_ = 0;
    max_out_degree_ = 0;
    max_label_frequency_ = 0;

    in_offsets_ = NULL;
    out_offsets_ = NULL;
    in_neighbors_ = NULL;
    out_neighbors_ = NULL;
    labels_ = NULL;

    reverse_index_offsets_ = NULL;
    reverse_index_ = NULL;
    labels_frequency_.clear();

    core_length_ = 0;
    
#if OPTIMIZED_LABELED_GRAPH == 1
    in_nlf_ = NULL;
    out_nlf_ = NULL;
#endif
}

Graph::~Graph() {
    delete[] in_offsets_;
    delete[] out_offsets_;
    delete[] in_neighbors_;
    delete[] out_neighbors_;
    delete[] labels_;
    delete[] reverse_index_offsets_;
    delete[] reverse_index_;
#if OPTIMIZED_LABELED_GRAPH == 1
    delete[] in_nlf_;
    delete[] out_nlf_;
#endif
}

void Graph::loadGraphFromFile(const std::string &file_path) {
    std::ifstream infile(file_path, std::ios::binary);

    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    // read basic characteristics
    infile.read((char*)&vertices_count_, 4);
    infile.read((char*)&labels_count_, 4);
    labels_count_--;

    //read labels
    reverse_index_offsets_ = new ui[labels_count_ + 1];
    reverse_index_ = new ui[vertices_count_];
    for (ui i = 0; i < vertices_count_; i++) {
        reverse_index_[i] = i;
    }
    infile.read((char*)reverse_index_offsets_, 4 * (labels_count_ + 1));
    labels_ = new LabelID[vertices_count_];
    for (ui i = 0; i < labels_count_; i++) {
        std::fill(labels_ + reverse_index_offsets_[i], labels_ + reverse_index_offsets_[i + 1], i);
    }
    for (ui i = 0; i < labels_count_; i++) {
        labels_frequency_[i] = reverse_index_offsets_[i + 1] - reverse_index_offsets_[i];
    }

    // read in-coming degree
    in_offsets_ = new ui[vertices_count_ + 1];
    in_offsets_[0] = 0;
    infile.read((char*)(in_offsets_ + 1), 4 * vertices_count_);
    for (ui i = 1; i <= vertices_count_; i++) {
        if (in_offsets_[i] > max_in_degree_) {
            max_in_degree_ = in_offsets_[i];
        }
    }
    for (ui i = 1; i <= vertices_count_; i++) {
        in_offsets_[i] = in_offsets_[i - 1] + in_offsets_[i];
    }
    in_edges_count_ = in_offsets_[vertices_count_];

    // read out-going degree
    out_offsets_ = new ui[vertices_count_ + 1];
    out_offsets_[0] = 0;
    infile.read((char*)(out_offsets_ + 1), 4 * vertices_count_);
    for (ui i = 1; i <= vertices_count_; i++) {
        if (out_offsets_[i] > max_out_degree_) {
            max_out_degree_ = out_offsets_[i];
        }
    }
    for (ui i = 1; i <= vertices_count_; i++) {
        out_offsets_[i] = out_offsets_[i - 1] + out_offsets_[i];
    }
    out_edges_count_ = out_offsets_[vertices_count_];

    // read in-coming neighbors
    in_neighbors_ = new ui[in_edges_count_];
    infile.read((char*)in_neighbors_, 4 * in_edges_count_);

    // read out-going neighbors
    out_neighbors_ = new ui[out_edges_count_];
    infile.read((char*)out_neighbors_, 4 * out_edges_count_);

    // extract information
    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    // sort for set intersection
    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(in_neighbors_ + in_offsets_[i], in_neighbors_ + in_offsets_[i + 1]);
    }
    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(out_neighbors_ + out_offsets_[i], out_neighbors_ + out_offsets_[i + 1]);
    }

    infile.close();

#if OPTIMIZED_LABELED_GRAPH == 1
    BuildNLF();
#endif

}

void Graph::loadGraphFromFileCompressed(const std::string &degree_path, const std::string &edge_path,
                                        const std::string &label_path) {
    
}

void Graph::printGraphMetaData() {
    std::cout << "|V|: " << vertices_count_ << ", |E-|(|E+|): " << in_edges_count_ << ", |\u03A3|: " << labels_count_ << std::endl;
    std::cout << "Max Incoming Degree: " << max_in_degree_ << ", Max Outgoing Degree: " << max_out_degree_ << ", Max Label Frequency: " << max_label_frequency_ << std::endl;
}

void Graph::buildCoreTable() {
    core_table_ = new int[vertices_count_];
    GraphOperations::getKCore(this, core_table_);

    for (ui i = 0; i < vertices_count_; ++i) {
        // 2-core
        if (core_table_[i] > 1) {
            core_length_ += 1;
        }
    }
}

void Graph::BuildReverseIndex() {
    reverse_index_ = new ui[vertices_count_];
    reverse_index_offsets_= new ui[labels_count_ + 1];
    reverse_index_offsets_[0] = 0;

    ui total = 0;
    // reverse_index_offsets_[0] = 0;
    // reverse_index_offsets_[1] = 0;
    // reverse_index_offsets_[i + 2] = reverse_index_offsets_[i + 1] + labels_frequency_[i];
    // reverse_index_offsets_[i] stores the start position of vertex with label i - 1;
    for (ui i = 0; i < labels_count_; ++i) {
        reverse_index_offsets_[i + 1] = total;
        total += labels_frequency_[i];
    }

    // write vertices and update reverse_index_offsets_ array so that reverse_index_offsets_[i] stores the start position of vertex with label i;
    for (ui i = 0; i < vertices_count_; ++i) {
        LabelID label = labels_[i];
        reverse_index_[reverse_index_offsets_[label + 1]++] = i;
    }
}

#if OPTIMIZED_LABELED_GRAPH == 1
void Graph::BuildNLF() {
    in_nlf_ = new std::unordered_map<LabelID, ui>[vertices_count_];
    out_nlf_ = new std::unordered_map<LabelID, ui>[vertices_count_];

    for (ui i = 0; i < vertices_count_; ++i) {
        ui in_neighbors_count;
        const VertexID * in_neighbors = getVertexInNeighbors(i, in_neighbors_count);

        for (ui j = 0; j < in_neighbors_count; ++j) {
            VertexID u = in_neighbors[j];
            LabelID label = getVertexLabel(u);
            if (in_nlf_[i].find(label) == in_nlf_[i].end()) {
                in_nlf_[i][label] = 0;
            }

            in_nlf_[i][label] += 1;
        }

        ui out_neighbors_count;
        const VertexID * out_neighbors = getVertexOutNeighbors(i, out_neighbors_count);

        for (ui j = 0; j < out_neighbors_count; ++j) {
            VertexID u = out_neighbors[j];
            LabelID label = getVertexLabel(u);
            if (out_nlf_[i].find(label) == out_nlf_[i].end()) {
                out_nlf_[i][label] = 0;
            }

            out_nlf_[i][label] += 1;
        }
    }
}
#endif