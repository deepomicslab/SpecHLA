//
// Created by yyh on 6/27/2019.
//

#ifndef SPECHAP_GRAPH_H
#define SPECHAP_GRAPH_H

#include <set>
#include "type.h"

//TODO: remove edge and redo clustering
//TODO: change the dfs method into disjointed set
//note that the variant graph and spectral uses matrix idx (insertion order) for calculation purpose

class VariantGraph              //linked list implementation for graph of variant, for checking connected component
{
public:
    uint variant_count;
    bool *visited;
    bool *filtered;
    bool *hap_visited;
    std::vector<std::set<uint>> graph;                //linked list implementation of graph, disconnected snp is represented empty set
    std::vector<std::set<uint>> haplotype_graph;
    std::unordered_map<uint, std::set<uint>> connected_component;

public:
    VariantGraph();
    explicit VariantGraph(uint variant_count);
    ~VariantGraph();
    bool fully_seperatable(uint idx);
    void load_connected_component();
    //zz: true:00, false:01
    inline void add_edge(uint start, uint end, bool zz)
    {
        graph[start].insert(end);
        graph[end].insert(start);
        if (zz)
        {
            haplotype_graph[2*start].insert(2*end);
            haplotype_graph[2*start+1].insert(2*end+1);
        }
        else
        {
            haplotype_graph[2*start].insert(2*end + 1);
            haplotype_graph[2*start+1].insert(2*end);
        }
    }
    inline bool check_disconnected_variant(uint idx)
    {
        return graph[idx].empty();
    }
    //TODO: remove this function
    inline void copy_uoset2set(uint  idx)
    {
        connected_component[idx];
    }
    void copy_insert_connected_component(uint idx);
    void dfs_util(std::set<uint> &dest, uint idx);
    inline void clear()
    {
        delete [] hap_visited;
        hap_visited = nullptr;
        delete []visited;
        visited = nullptr;
        delete []filtered;
        filtered = nullptr;
        haplotype_graph.clear();
        graph.clear();
        connected_component.clear();
    }
    inline void reset(uint variant_count)
    {
        //delete[] visited;
        graph.clear();
        haplotype_graph.clear();
        connected_component.clear();
        this->variant_count = variant_count;
        visited = new bool[variant_count];
        filtered = new bool[variant_count];
        hap_visited = new bool[2 * variant_count];
        for (uint i = 0; i < variant_count; i++)
        {
            visited[i] = false;
            filtered[i] = false;
            hap_visited[2*i] = hap_visited[2*i + 1] = false;
            graph.emplace_back();
            haplotype_graph.emplace_back();
            haplotype_graph.emplace_back();
        }
    }
    // for query with connected components
    inline bool contain(uint idx) { return connected_component.find(idx) != connected_component.end();}
    inline bool disjointedatpos(uint idx) { return connected_component.count(idx) == 0;}
    inline void remove_variant(uint idx)
    {
        filtered[idx] = true;
    }
    void split_after_filtering(uint idx, std::vector<std::set<uint>> &blks);
};



#endif //LAPHAP_GRAPH_H
