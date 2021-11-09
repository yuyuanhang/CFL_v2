#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
#include <stdlib.h>

typedef unsigned int ui;

typedef uint32_t VertexID;
typedef ui LabelID;

#define IN 0
#define OUT 1

class TreeNode {
public:
    VertexID id_;
    VertexID parent_;
    bool p_direction_;
    ui level_;
    ui under_level_count_;
    ui children_count_;
    ui bn_count_;
    ui fn_count_;
    VertexID* under_level_;
    VertexID* children_;
    VertexID* bn_;
    VertexID* fn_;
    size_t estimated_embeddings_num_;
    bool* c_directions_;
    bool* bn_directions_;
    bool* fn_directions_;
    bool* under_level_directions_;
public:
    TreeNode() {
        id_ = 0;
        under_level_ = NULL;
        bn_ = NULL;
        fn_ = NULL;
        children_ = NULL;
        c_directions_ = NULL;
        parent_ = 0;
        level_ = 0;
        under_level_count_ = 0;
        children_count_ = 0;
        bn_count_ = 0;
        fn_count_ = 0;
        estimated_embeddings_num_ = 0;
    }

    ~TreeNode() {
        delete[] under_level_;
        delete[] bn_;
        delete[] fn_;
        delete[] children_;
        delete[] c_directions_;
        delete[] bn_directions_;
        delete[] fn_directions_;
        delete[] under_level_directions_;
    }

    void initialize(const ui size) {
        under_level_ = new VertexID[size];
        bn_ = new VertexID[size];
        fn_ = new VertexID[size];
        children_ = new VertexID[size];
        c_directions_ = new bool[size];
        bn_directions_ = new bool[size];
        fn_directions_ = new bool[size];
        under_level_directions_ = new bool[size];
    }
};

class Edges {
public:
    ui* offset_;
    ui* edge_;
    ui vertex_count_;
    ui edge_count_;
    ui max_degree_;
    bool direction_;
public:
    Edges() {
        offset_ = NULL;
        edge_ = NULL;
        vertex_count_ = 0;
        edge_count_ = 0;
        max_degree_ = 0;
    }

    ~Edges() {
        delete[] offset_;
        delete[] edge_;
    }
};

#endif