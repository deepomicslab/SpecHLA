//
// Created by yyh on 6/27/2019.
//
#include "graph.h"
#include <list>
//---------------------------------------VariantGraph---------------------------------------
VariantGraph::VariantGraph()
        : visited(nullptr), hap_visited(nullptr), filtered(nullptr) {}

VariantGraph::VariantGraph(uint variant_count)
        : variant_count(variant_count), hap_visited(nullptr)
{
    visited = new bool[variant_count];
    filtered = new bool[variant_count];
    hap_visited = new bool[2 * variant_count];
    for (uint i = 0; i < variant_count; i++)
    {
        visited[i] = false;
        filtered[i] = false;
        hap_visited[2 * i] = hap_visited[2 * i + 1] = false;
        graph.emplace_back();
        haplotype_graph.emplace_back();
        haplotype_graph.emplace_back();
    }
}

VariantGraph::~VariantGraph()
{
    delete[] visited;
    visited = nullptr;
    delete[] hap_visited;
    hap_visited = nullptr;
    delete[] filtered;
    filtered = nullptr;
    graph.clear();
    haplotype_graph.clear();
    connected_component.clear();
}

void VariantGraph::load_connected_component()
{
    for (uint i = 0; i < variant_count; i++)
    {
        if (visited[i])         // if a variant is part of connected component, just skip
            continue;
        if (check_disconnected_variant(i))
            copy_uoset2set(i);  // a disjointed variant, can be accessed with starting index (key of map)
        else
            copy_insert_connected_component(i); //connected component, can be accessed with starting index (key of map)
    }
}

void VariantGraph::copy_insert_connected_component(uint idx)
{
    std::set<uint> dest;
    dfs_util(dest, idx);
    connected_component.emplace(idx, dest);
}

void VariantGraph::dfs_util(std::set<uint> &dest, uint idx)
{
    if (visited[idx] || filtered[idx])
        return;
    auto &current_variant = graph[idx];
    visited[idx] = true;
    for (auto i : current_variant)
    {
        if (filtered[i])
            continue;
        dest.insert(i);
    }
    for (auto i : current_variant)
        dfs_util(dest, i);
}

//call this function after find connected component
bool VariantGraph::fully_seperatable(uint idx)
{
    uint s = 2 * idx;
    uint d = 2 * idx + 1;

    std::list<uint> queue;
    hap_visited[s] = true;
    queue.push_back(s);

    while (!queue.empty())
    {
        s = queue.front();
        queue.pop_front();
        for (auto i : haplotype_graph[s])
        {
            if (this->filtered[s/2])
                continue;
            if (i == d)
                return false;
            if (!hap_visited[i])
            {
                hap_visited[i] = true;
                queue.push_back(i);
            }
        }
    }

    return true;
}

void VariantGraph::split_after_filtering(uint idx, std::vector<std::set<uint>> &blks)
{
    std::set<uint> &sub_graph_variant = this->connected_component[idx];
    for (auto &i : sub_graph_variant)
        visited[i] = false;

    for (auto &i : sub_graph_variant)
    {
        if (visited[i] || filtered[i])
            continue;
        blks.emplace_back();
        dfs_util(blks.back(), i);
    }
}