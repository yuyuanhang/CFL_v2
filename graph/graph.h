#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include "configuration/types.h"
#include "configuration/config.h"

class Graph {
private:
    ui vertices_count_;
    ui in_edges_count_;
    ui out_edges_count_;
    ui labels_count_; // for the construction of reverse index (label index).
    ui max_degree_;
    ui max_in_degree_; // for 2-core decomposition
    ui max_out_degree_;
    ui max_label_frequency_; // upper bound of the size of candidate set.

    ui* in_offsets_;
    ui* out_offsets_;
    VertexID * in_neighbors_;
    VertexID * out_neighbors_;
    LabelID* labels_;
    // store vertices based on labels so that we can directly access vertices with the same label.
    ui* reverse_index_offsets_;
    ui* reverse_index_;

    std::unordered_map<LabelID, ui> labels_frequency_;

    //for 2-core decomposition
    int* core_table_;
    ui core_length_;

#if OPTIMIZED_LABELED_GRAPH == 1
    std::unordered_map<LabelID, ui>* in_nlf_;
    std::unordered_map<LabelID, ui>* out_nlf_;
#endif

private:
    void BuildReverseIndex();

#if OPTIMIZED_LABELED_GRAPH == 1
    void BuildNLF();
#endif

public:
    Graph();
    ~Graph();

    void loadGraphFromFile(const std::string& file_path);
    void loadGraphFromFileCompressed(const std::string& degree_path, const std::string& edge_path,
                                     const std::string& label_path);
    void printGraphMetaData();
    void buildCoreTable();

    void printCoreTable() {
        for (ui i = 0; i < vertices_count_; ++i) {
            std::cout << core_table_[i] << std::endl;
        }
    }

    const ui getVerticesCount() const {
        return vertices_count_;
    }

    const ui getEdgesCount() const {
        return in_edges_count_;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const ui getGraphInMaxDegree() const {
        return max_in_degree_;
    }

    const ui getGraphOutMaxDegree() const {
        return max_out_degree_;
    }

    const ui getVertexDegree(const VertexID id) const {
        return getVertexInDegree(id) + getVertexOutDegree(id);
    }

    const ui getVertexInDegree(const VertexID id) const {
        return in_offsets_[id + 1] - in_offsets_[id];
    }

    const ui getVertexOutDegree(const VertexID id) const {
        return out_offsets_[id + 1] - out_offsets_[id];
    }

    const ui get2CoreSize() const {
        return core_length_;
    }

    const ui getCoreValue(const VertexID id) const {
        return core_table_[id];
    }

    const LabelID getVertexLabel(const VertexID id) const {
        return labels_[id];
    }

    const ui getLabelsFrequency(const LabelID label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }

    const std::unordered_map<LabelID, ui>* getVertexNLF(const VertexID id) const {
        return in_nlf_ + id;
    }

    const std::unordered_map<LabelID, ui>* getVertexInNLF(const VertexID id) const {
        return in_nlf_ + id;
    }

    const std::unordered_map<LabelID, ui>* getVertexOutNLF(const VertexID id) const {
        return out_nlf_ + id;
    }

    const ui * getVerticesByLabel(const LabelID id, ui& count) const {
        count = reverse_index_offsets_[id + 1] - reverse_index_offsets_[id];
        return reverse_index_ + reverse_index_offsets_[id];
    }

    const ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = in_offsets_[id + 1] - in_offsets_[id];
        return in_neighbors_ + in_offsets_[id];
    }

    const ui * getVertexInNeighbors(const VertexID id, ui& count) const {
        count = in_offsets_[id + 1] - in_offsets_[id];
        return in_neighbors_ + in_offsets_[id];
    }

    const ui * getVertexOutNeighbors(const VertexID id, ui& count) const {
        count = out_offsets_[id + 1] - out_offsets_[id];
        return out_neighbors_ + out_offsets_[id];
    }

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        ui count = 0;
        const VertexID* in_neighbors = getVertexInNeighbors(v, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (in_neighbors[mid] == u) {
                return true;
            }
            else if (in_neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        count = 0;
        const VertexID* out_neighbors = getVertexOutNeighbors(v, count);

        begin = 0;
        end = count - 1;

        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (out_neighbors[mid] == u) {
                return true;
            }
            else if (out_neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }
};

#endif