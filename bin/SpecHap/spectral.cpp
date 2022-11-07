//
// Created by yonghanyu2 on 8/10/2018.
//

#include "spectral.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include "indexing.h"
#include "phaser.h"
#include "type.h"
#include <limits>
#include <unordered_set>
#include <cmath>
#include <tuple>
#include <ctime>

void print_to_file(const char *filename, GMatrix &mat)
{
    std::ofstream ofs(filename);
    ofs.precision(std::numeric_limits<double>::max_digits10);
    ofs << mat;
    ofs.close();
}

void print_to_file(BarcodeLinkers *barcodes, const char *filename)
{
    std::ofstream ofs(filename);
    for (auto &barcode: barcodes->barcodes)
    {
        ofs << barcode.first << " ";
        for (auto &info : barcode.second.barcode_info)
            ofs << info.first << " " << info.second.first << " ";
        ofs << "\n";
    }
}

void print_map(std::unordered_map<uint, uint> &buff, uint key)
{
    std::cout << buff[key] << std::endl;
}

void print(double ** matrix, int n)
{
    int i, j;
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
            std::cout <<  matrix[i][j] << " ";
        std::cout << "\n";
    }
}

void print_hap(ptr_ResultforSingleVariant &block)
{
    if (!block->variant_phased())
        std::cout << -1 << std::endl;
    if (block->is_REF())
        std::cout << 0 << std::endl;
    else
        std::cout << 1 << std::endl;
}

void print_hap(ptr_PhasedBlock &block, uint idx)
{
    if (!block->results[idx]->variant_phased())
        std::cout << -1 << std::endl;
    if (block->results[idx]->is_REF())
        std::cout << 0 << std::endl;
    else
        std::cout << 1 << std::endl;
}

// initialization
Spectral::Spectral(std::vector<FragmentReader *>& frfrags, std::vector<double>& fr_weights, BEDReader *frbed, double threhold, int coverage, bool use_secondary)
        :use_secondary(use_secondary), raw_graph(nullptr), raw_count(nullptr), threhold(threhold), epsilon(10e-2),
        coverage(coverage), q_aver(0.0),
        q_sum(0.0), barcode_linker(nullptr), barcode_linker_index_set(false), frbed(frbed), phasing_window(nullptr), chromo_phaser(nullptr)
{
    for (auto item : frfrags) {
        frs.push_back(item);
    }
    block_no = 1;
    this->barcode_linker = new BarcodeLinkers(MAX_BARCODE_SPANNING);
    this->region_frag_stats = new RegionFragStats();
    this->barcode_linker->regionFragStats = region_frag_stats;
    this->offset = 0;
    for (auto item : fr_weights) {
        this->fr_weights.push_back(item);
    }
}

Spectral::~Spectral()
{
    delete[] raw_graph;
    raw_graph = nullptr;
    delete[] raw_count;
    raw_count = nullptr;
//    fr = nullptr;
    frs.clear();
    frs.shrink_to_fit();
    frbed = nullptr;
    delete barcode_linker;
    barcode_linker = nullptr;
    delete region_frag_stats;
    region_frag_stats = nullptr;
}

void Spectral::clean()
{
    barcode_linker_index_set = false;
    delete[] raw_graph;
    raw_graph = nullptr;
    delete[] raw_count;
    raw_count = nullptr;
    variant_graph.clear();
    barcode_linker->clear();
    phased_block_starts.clear();
    block_tobe_split.clear();
}

void Spectral::reset()
{
    this->clean();
    this->variant_count = phasing_window->curr_window_block_count();
    this->start_variant_idx = phasing_window->start;
    this->end_variant_idx_intended = phasing_window->intended_end;
    this->end_variant_idx_overlap = phasing_window->end;
    this->n = 2 * variant_count;
    this->q_sum = 0;
    for (auto item : frs) {
        item->set_window_info(start_variant_idx, end_variant_idx_overlap, end_variant_idx_intended);
    }
//    fr->set_window_info(start_variant_idx, end_variant_idx_overlap, end_variant_idx_intended);
    this->variant_graph.reset(variant_count);
    this->raw_graph = new double[n * n];
    this->raw_count = new int[n * n];

    for (unsigned int i = 0; i < n * n; i++)
    {
        this->raw_graph[i] = 0.0;
        this->raw_count[i] = 0;
    }
    //    read graph
    ViewMap weighted_graph(raw_graph, n, n);
    CViewMap count_graph(raw_count, n, n);
    //    used for lst that large than window_size
    set_prev_buff(weighted_graph, count_graph);
    this->frag_buffer.clear();

    for (int i = 0; i < OPERATIONS.size(); i++) {
        double w = this->fr_weights[i];
//        if (i == 1) {
//            w = 0.5;
//        }
        auto item = OPERATIONS[i];
        if (item == MODE_10X)
            read_fragment_10x(i,weighted_graph,count_graph);
        else if (item == MODE_HIC)
            read_fragment_hic(i,weighted_graph,count_graph);
        else if (item == MODE_PE)
            read_fragment(i,weighted_graph,count_graph, w);
        else if (item == MODE_NANOPORE)
            read_fragment_nanopore(i,weighted_graph,count_graph);
        else if (item == MODE_PACBIO)
            read_fragment_pacbio(i,weighted_graph,count_graph, w);
        else if (item == MODE_MATRIX)
            read_fragment_matrix(i,weighted_graph,count_graph, w);
    }
//    cal prob graph
//    if (HAS_TENX){
//        add_snp_edge_barcode(weighted_graph, count_graph);
//    }
    cal_prob_matrix(weighted_graph, count_graph, nullptr, nullptr, nullptr);
}
void Spectral::set_prev_buff(ViewMap& weighted_graph, CViewMap& count_graph) {
    double w = 1;
    for (auto &item : this->frag_buffer) {
        if (item.snps.back().first > this->start_variant_idx) {
            this->add_snp_edge(item, weighted_graph,count_graph,1);
        }
    }
}


GMatrix Spectral::slice_submat(std::set<uint> &variants_mat, GMatrix &adj_mat)
{
    Eigen::ArrayXi index(variants_mat.size() * 2);
    uint count = 0;
    for (auto j = variants_mat.begin(); j != variants_mat.end(); ++j, count++)
    {
        index(2 * count) = 2 * (*j);
        index(2 * count + 1) = 2 * (*j) + 1;
    }
    GMatrix sub_mat = mat_indexing(adj_mat, index, index);
    return sub_mat;
}

GMatrix Spectral::slice_submat(std::set<uint> &variants_mat)
{
    Eigen::Ref<GMatrix> adj_mat = this->adjacency_matrix;
    Eigen::ArrayXi index(variants_mat.size() * 2);
    uint count = 0;
    for (auto j = variants_mat.begin(); j != variants_mat.end(); ++j, count++)
    {
        index(2 * count) = 2 * (*j); 
        index(2 * count + 1) = 2 * (*j) + 1;
    }
    GMatrix sub_mat = mat_indexing(adj_mat, index, index);
    return sub_mat;
}

CMatrix Spectral::slice_submat(std::set<uint> &variants_mat, bool t, CMatrix &adj_mat)
{
    Eigen::ArrayXi index(variants_mat.size() * 2);
    uint count = 0;
    for (auto j = variants_mat.begin(); j != variants_mat.end(); ++j, count++)
    {
        index(2 * count) = 2 * (*j);
        index(2 * count + 1) = 2 * (*j) + 1;
    }
    CMatrix sub_mat = mat_indexing(adj_mat, index, index);
    return sub_mat;
}

CMatrix Spectral::slice_submat(std::set<uint> &variants_mat, bool t)
{
    Eigen::Ref<CMatrix> adj_mat = this->count_matrix;
    Eigen::ArrayXi index(variants_mat.size() * 2);
    uint count = 0;
    for (auto j = variants_mat.begin(); j != variants_mat.end(); ++j, count++)
    {
        index(2 * count) = 2 * (*j);
        index(2 * count + 1) = 2 * (*j) + 1;
    }
    CMatrix sub_mat = mat_indexing(adj_mat, index, index);
    return sub_mat;
}

// read fragment matrix
void Spectral::read_fragment(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph, double w)
{
    auto fr = frs[frIdx];
//    this->frag_buffer.clear();
    Fragment fragment;
//    ViewMap weighed_graph(raw_graph, n, n);
//    CViewMap count_graph(raw_count, n, n);
    while (fr->get_next_pe(fragment))
    {
//        std::cout<<fragment.read_qual<<std::endl;
        add_snp_edge(fragment, weighted_graph, count_graph, w);
        this->frag_buffer.push_back(fragment);
        fragment.reset();
    }
//    cal_prob_matrix(weighted_graph, count_graph, nullptr, nullptr, nullptr);
}

void Spectral::read_fragment_10x(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph)
{
    auto fr = frs[frIdx];
    this->region_frag_stats->clear();
    uint start_pos = 0, end_pos = 0;
    get_current_window_pos(start_pos, end_pos);
    //TODO: make start_pos dynamic to the size of phasing window
    frbed->read_region_frag(this->chromo_phaser->chr_name.c_str(), start_pos, end_pos, this->region_frag_stats);
    //this->barcode_linker->regionFragStats = region_frag_stats;
    Fragment fragment;
//    ViewMap weighed_graph(raw_graph, n, n);
//    CViewMap count_graph(raw_count, n, n);

    while (fr->get_next_tenx(fragment))
    {
        //add_snp_edge(fragment, weighed_graph);
        populate_variant_info(fragment);
        if (fragment.start == 4294967205) {
            int dd = 9;
        }
        if (fragment.barcode != "NULL")
            barcode_linker->insert_barcode(fragment);
            //add_barcode_info(fragment, barcode_linker);
        fragment.reset();
    }
    add_snp_edge_barcode(weighted_graph, count_graph);
//    cal_prob_matrix(weighted_graph, count_graph, nullptr, nullptr, nullptr);
}

//TODO: only unused HiC linker should be stored
void Spectral::read_fragment_hic(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph)
{
    auto fr = frs[frIdx];
    Fragment fragment;
//    this->frag_buffer.clear();
//    CViewMap count_graph(raw_count, n, n);
//    ViewMap weighed_graph(raw_graph, n, n);
    while (fr->get_next_hic(fragment))
    {
        add_snp_edge(fragment, weighted_graph, count_graph,1);
        this->frag_buffer.push_back(fragment);
        if ( fragment.snps[0].first >= phasing_window->prev_window_start)
            if (fragment.insertion_size >= 5000 && fragment.insertion_size <= 40000000)
                this->hic_linker_container.add_HiC_info(fragment);
        fragment.reset();
    }
//    cal_prob_matrix(weighted_graph, count_graph, nullptr, nullptr, nullptr);
}

void Spectral::read_fragment_nanopore(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph)
{
    auto fr = frs[frIdx];
//    this->frag_buffer.clear();
    Fragment fragment;
//    ViewMap weighed_graph(raw_graph, n, n);
//    CViewMap count_graph(raw_count, n, n);
    while (fr->get_next_nanopore(fragment))
    {
        add_snp_edge(fragment, weighted_graph, count_graph,1);
        this->frag_buffer.push_back(fragment);
        fragment.reset();
    }
    //this->q_aver = q_sum / this->frag_buffer.size() / 6;
//    cal_prob_matrix(weighted_graph, count_graph, nullptr, nullptr, nullptr);
}

void Spectral::read_fragment_pacbio(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph, double w)
{
    auto fr = frs[frIdx];
//    this->frag_buffer.clear();
    Fragment fragment;
//    ViewMap weighed_graph(raw_graph, n, n);
//    CViewMap count_graph(raw_count, n, n);
    while (fr->get_next_pacbio(fragment))
    {
        add_snp_edge(fragment, weighted_graph, count_graph,w);
        this->frag_buffer.push_back(fragment);
        fragment.reset();
    }
//    cal_prob_matrix(weighted_graph, count_graph, nullptr, nullptr, nullptr);
}

void Spectral::read_fragment_matrix(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph, double w){
    auto fr = frs[frIdx];
    std::vector<Fragment> frags;
    for (int i = 0; i < 4; i++) {
        Fragment f;
        frags.push_back(f);
    }

    while (fr->get_next_matrix(frags))
    {
        for (int i = 0; i < 4; i++) {
            Fragment &f = frags[i];
            add_snp_edge(f, weighted_graph, count_graph,w);
            this->frag_buffer.push_back(frags[i]);
            f.reset();
        }
    }
}


// aid function
void Spectral::add_snp_edge(Fragment &fragment, ViewMap &weighted_graph, CViewMap &count_graph, double w)
{
    if (fragment.snps.empty())
        return;

    fragment.read_qual == 0 ? fragment.read_qual -= 1 : 0;
    std::map<uint, std::pair<double,double> > block_supporting_likelyhood;
    for (auto &i : fragment.snps) {
        if (!phasing_window->in_range(i.first))
            continue;
        auto a = this->phasing_window->var_idx2mat_idx(i.first);
        if (a < 0 or a >= this->variant_count)
            continue;
        uint blk_idx = this->phasing_window->mat_idx2var_idx(a);
        ptr_PhasedBlock blk = this->phasing_window->blocks[blk_idx];
        if (blk->results.count(i.first) == 0)
            continue;
        if (block_supporting_likelyhood.count(blk_idx) == 0)
            block_supporting_likelyhood[blk_idx] = std::make_pair(0.0, 0.0);

        if ((i.second.first == 0 && blk->results[i.first]->is_REF()) ||
            (i.second.first == 1 && blk->results[i.first]->is_ALT())) {
            if (fragment.getType() == FRAG_MATRIX) {
//                Here, beacause the matrix format with count 00(11) twice
                block_supporting_likelyhood[blk_idx].first += i.second.second / 2;
                block_supporting_likelyhood[blk_idx].second += 0;
            } else {
                block_supporting_likelyhood[blk_idx].first += log10(i.second.second);
                block_supporting_likelyhood[blk_idx].second += log10(1 - i.second.second);
            }
        } else {
            if (fragment.getType() == FRAG_MATRIX) {
                block_supporting_likelyhood[blk_idx].first += 0;
                block_supporting_likelyhood[blk_idx].second += i.second.second / 2;
            } else {
                block_supporting_likelyhood[blk_idx].first += log10(1 - i.second.second);
                block_supporting_likelyhood[blk_idx].second += log10(i.second.second);
            }
        }

    }

    for (auto &i : block_supporting_likelyhood)
    {
        auto a = this->phasing_window->var_idx2mat_idx(i.first);

        for (auto &j : block_supporting_likelyhood)
        {
            auto b = this->phasing_window->var_idx2mat_idx(j.first);
            if (a == b)
                continue;
            double P_H1, P_H2;
            if (fragment.getType() == FRAG_MATRIX) {
                if ((i.second.first != 0 && j.second.first != 0 ) || (i.second.second != 0 && j.second.second != 0)) {
                    P_H1 = (i.second.first + j.second.first + i.second.second + j.second.second) / 2;
                    P_H2 = 0;
                } else {
                    P_H1 = 0;
                    P_H2 = (i.second.first + j.second.second + i.second.second + j.second.first) / 2;
                }
            } else {
                P_H1 = cal_score(i.second.first + j.second.first, i.second.second + j.second.second);
                P_H2 = cal_score(i.second.first + j.second.second, i.second.second + j.second.first);
            }
//            double P_H1 = cal_score(i.second.first + j.second.first, i.second.second + j.second.second);
//            double P_H2 = cal_score(i.second.first + j.second.second, i.second.second + j.second.first);
            if (P_H2 == P_H1)
                continue;
            double score = abs(P_H1 - P_H2);
            int connection;
            P_H1 > P_H2 ? connection = 0 : connection = 1;
            //connection == 0 ? variant_graph.add_edge(a, b, true) : variant_graph.add_edge(a, b, false);

            //weighted_graph(2*a, 2*b + connection) += score;
            //weighted_graph(2*a + 1, 2*b + 1 - connection) += score;
            //count_graph(2 * a, 2 * b + connection) ++;
            //count_graph(2 * a + 1, 2 * b + 1 - connection) ++;

            weighted_graph(2*a, 2*b) += P_H1*w;
            weighted_graph(2*a + 1, 2*b + 1) += P_H1*w;
            weighted_graph(2*a, 2*b + 1) += P_H2*w;
            weighted_graph(2*a + 1, 2*b) += P_H2*w;
        }
    }
}

void Spectral::add_snp_edge_barcode(ViewMap &weighted_graph, CViewMap &count_graph)
{

    for (auto &content : this->barcode_linker->barcodes)
    {
        std::map<uint, std::pair<double,double> > block_supporting_likelyhood;
        Barcode &barcode = content.second;
        for (auto &i : barcode.barcode_info)
        {
            if (!phasing_window->in_range(i.first))
                continue;
            auto a = this->phasing_window->var_idx2mat_idx(i.first);
            if (a < 0 or a >= this->variant_count)
                continue;
            uint blk_idx = this->phasing_window->mat_idx2var_idx(a);
            ptr_PhasedBlock blk = this->phasing_window->blocks[blk_idx];
            if (blk->results.count(i.first) == 0)
                continue;
            if (block_supporting_likelyhood.count(blk_idx) == 0)
                block_supporting_likelyhood[blk_idx] = std::make_pair(0.0, 0.0);

            if ((i.second.first == 0 && blk->results[i.first]->is_REF()) ||
                (i.second.first == 1 && blk->results[i.first]->is_ALT())) {
                block_supporting_likelyhood[blk_idx].first += log10(i.second.second);
                block_supporting_likelyhood[blk_idx].second += log10(1 - i.second.second);
            } else {
                block_supporting_likelyhood[blk_idx].first += log10(1 - i.second.second);
                block_supporting_likelyhood[blk_idx].second += log10(i.second.second);
            }
        }

        for (auto &i : block_supporting_likelyhood)
        {
            auto a = this->phasing_window->var_idx2mat_idx(i.first);

            for (auto &j : block_supporting_likelyhood)
            {
                auto b = this->phasing_window->var_idx2mat_idx(j.first);
                if (a == b)
                    continue;
                double P_H1 = cal_score(i.second.first + j.second.first, i.second.second + j.second.second);
                double P_H2 = cal_score(i.second.first + j.second.second, i.second.second + j.second.first);
                if (P_H2 == P_H1)
                    continue;
                double score = abs(P_H1 - P_H2);
                int connection;
                P_H1 > P_H2 ? connection = 0 : connection = 1;
                //connection == 0 ? variant_graph.add_edge(a, b, true) : variant_graph.add_edge(a, b, false);

                //weighted_graph(2*a, 2*b + connection) += score;
                //weighted_graph(2*a + 1, 2*b + 1 - connection) += score;
                //count_graph(2 * a, 2 * b + connection) ++;
                //count_graph(2 * a + 1, 2 * b + 1 - connection) ++;

                weighted_graph(2*a, 2*b) += P_H1;
                weighted_graph(2*a + 1, 2*b + 1) += P_H1;
                weighted_graph(2*a, 2*b + 1) += P_H2;
                weighted_graph(2*a + 1, 2*b) += P_H2;
            }
        }
    }
    //std::cout << weighted_graph << std::endl;
}

void Spectral::add_snp_edge_hic(ViewMap &weighted_graph, CViewMap &count_graph)
{
    for (auto it = this->hic_linker_container.linker.begin(); it != this->hic_linker_container.linker.end(); it++)
    {
        std::map<uint, std::pair<double,double> > block_supporting_likelyhood;
        HiCLinker &hic_linker = it->second;
        for (auto &i : hic_linker.hic_info)
        {
            if (!phasing_window->in_range(i.first))
                continue;
            auto a = this->phasing_window->var_idx2mat_idx(i.first);
            if (a < 0 or a >= this->variant_count)
                continue;
            uint blk_idx = this->phasing_window->mat_idx2var_idx(a);
            ptr_PhasedBlock blk = this->phasing_window->blocks[blk_idx];
            if (blk->results.count(i.first) == 0)
                continue;
            if (block_supporting_likelyhood.count(blk_idx) == 0)
                block_supporting_likelyhood[blk_idx] = std::make_pair(0.0, 0.0);

            if ((i.second.first == 0 && blk->results[i.first]->is_REF()) ||
                (i.second.first == 1 && blk->results[i.first]->is_ALT())) {
                block_supporting_likelyhood[blk_idx].first += log10(i.second.second);
                block_supporting_likelyhood[blk_idx].second += log10(1 - i.second.second);
            } else {
                block_supporting_likelyhood[blk_idx].first += log10(1 - i.second.second);
                block_supporting_likelyhood[blk_idx].second += log10(i.second.second);
            }

        }
        for (auto &i : block_supporting_likelyhood)
        {
            auto a = this->phasing_window->var_idx2mat_idx(i.first);

            for (auto &j : block_supporting_likelyhood)
            {
                auto b = this->phasing_window->var_idx2mat_idx(j.first);
                if (a == b)
                    continue;
                double P_H1 = cal_score(i.second.first + j.second.first, i.second.second + j.second.second);
                double P_H2 = cal_score(i.second.first + j.second.second, i.second.second + j.second.first);
                if (P_H2 == P_H1)
                    continue;
                double score = abs(P_H1 - P_H2);
                int connection;
                P_H1 > P_H2 ? connection = 0 : connection = 1;

                weighted_graph(2*a, 2*b) += P_H1;
                weighted_graph(2*a + 1, 2*b + 1) += P_H1;
                weighted_graph(2*a, 2*b + 1) += P_H2;
                weighted_graph(2*a + 1, 2*b) += P_H2;
            }
        }
    }
}

void Spectral::cal_prob_matrix(ViewMap &weighted_graph, CViewMap &count_graph, GMatrix *weight, CMatrix *count, VariantGraph *variant_graph)
{
    GMatrix *ptr_adj_mat;
    CMatrix *ptr_count_mat;
    VariantGraph *ptr_var_graph;
    if (weight == nullptr)
    {
        ptr_adj_mat = &this->adjacency_matrix;
        ptr_count_mat = &this->count_matrix;
        ptr_var_graph = &this->variant_graph;
    }
    else
    {
        ptr_adj_mat = weight;
        ptr_count_mat = count;
        ptr_var_graph = variant_graph;
    }
    uint var_count = ptr_var_graph -> variant_count;
    GMatrix &adj_mat = *ptr_adj_mat;
    CMatrix &count_mat = *ptr_count_mat;
    VariantGraph &var_graph = *ptr_var_graph;
    auto tsize = weighted_graph.size();

    adj_mat = weighted_graph; // + GMatrix::Identity(2*var_count, 2*var_count);

    for (int i = 0; i < var_count; i++)
    {
        for (int j = 0; j < var_count; j++)
        {
            if (i == j)
                continue;
            //adj_mat(i, j)  = abs(log10(weighted_graph(i, j)));
            //the connection provides no sufficient information for phasing
            auto tm1 = weighted_graph(2 * i, 2*j);
            auto tm3 = weighted_graph(2*i, 2*j+1);
            double score = weighted_graph(2 * i, 2*j) - weighted_graph(2*i, 2*j+1);
            if (score > 0)
            {
                adj_mat(2*i, 2*j) = adj_mat(2*i + 1, 2 * j + 1) = score;
                adj_mat(2*i, 2*j + 1) = adj_mat(2*i + 1, 2 * j) = 0;
                count_graph(2*i, 2*j) ++;
                count_graph(2*i + 1, 2*j + 1) ++;
                var_graph.add_edge(i, j, true);
                
            }
            else if (score < 0)
            {
                adj_mat(2*i, 2*j) = adj_mat(2*i + 1, 2 * j + 1) = 0;
                adj_mat(2*i, 2*j + 1) = adj_mat(2*i + 1, 2 * j) = -1 * score;
                count_graph(2*i, 2*j + 1) ++;
                count_graph(2*i + 1, 2*j) ++;
                var_graph.add_edge(i, j, false);
                //var_graph.add_edge(i, j, true);
            }

            else
            {
                adj_mat(2*i, 2*j) = adj_mat(2*i + 1, 2 * j + 1) = 0;
                adj_mat(2*i, 2*j + 1) = adj_mat(2*i + 1, 2 * j) = 0;
            }
            //adj_mat(2*i, 2*j) = adj_mat(2*i + 1, 2 * j + 1) = abs(log10(weighted_graph(2*i, 2*j)));
            //adj_mat(2*i, 2*j + 1) = adj_mat(2*i + 1, 2 * j ) = abs(log10(weighted_graph(2*i, 2*j + 1)));
        }
    }

    count_mat = count_graph;

    //if (op == op_mode::TENX)
    //{
    //    for (auto &i : phasing_window->current_window_idxes)
    //        cal_10x_filter(phasing_window->blocks[i]);
    //}
    var_graph.load_connected_component();
}

void Spectral::solver_recursive()
{
    std::unordered_set<uint> met_idx;
    //load_hic_linker();
    Eigen::Ref<GMatrix> adj_mat = adjacency_matrix;
    uint mat_idx;

    //for each block (variant) in phasing window
    for (auto i : phasing_window->current_window_idxes)
    {
        mat_idx = phasing_window->var_idx2mat_idx(i);
        //not met yet
        if (met_idx.find(mat_idx) == met_idx.end())
            met_idx.insert(mat_idx);
        else
            continue;
        // variant inside phasing block
        if (!variant_graph.contain(mat_idx))
            continue;

            // disjointed snp, do nothing
        else if (variant_graph.disjointedatpos(mat_idx))
            //phasing_window->insert_phased_block_starting_idx(i);
            ;

            // connected phased block
        else
        {
            std::set<uint> &variants_mat = variant_graph.connected_component[mat_idx];
            Eigen::ArrayXi index(variants_mat.size() * 2);
            uint count = 0;
            for (auto j = variants_mat.begin(); j != variants_mat.end(); ++j, count++)
            {
                index(2 * count) = 2 * (*j);
                index(2 * count + 1) = 2 * (*j) + 1;
            }
            GMatrix sub_mat = mat_indexing(adj_mat, index, index);
            if (variant_graph.fully_seperatable(mat_idx))
                find_connected_component(sub_mat, variants_mat);
            else
            {
                int block_count = 0;
                std::map<uint, int> subroutine_map;
                std::map<uint, uint> subroutine_blk_start;
                std::map<uint, double> block_qualities;
                call_haplotype(sub_mat, variants_mat, block_count, subroutine_map, subroutine_blk_start, false, block_qualities);
            }
            block_no++;
        }
    }
}

//TODO: change for loop into direct indexing

void Spectral::solver()
{
    std::unordered_set<uint> met_idx;
    reset();
    uint mat_idx;
    //for each block (variant) in phasing window
    auto var_count = this->variant_count;
    for (int k = 0; k < var_count - 1; k++) {
        auto j = k + 1;
//        for (int j = i; j < var_count; j++) {
//            if (i == j)
//                continue;
        //adj_mat(i, j)  = abs(log10(weighted_graph(i, j)));
        //the connection provides no sufficient information for phasing
//        auto break_idx = this->phasing_window->mat_idx2var_idx(j);
//        if(break_idx +451839  == 465951) {
//            auto tmpdebug = 1;
//        }
        auto tsize = this->adjacency_matrix.size();
        auto tm1 = this->adjacency_matrix(2 * k, 2 * j);
        auto tm3 = this->adjacency_matrix(2 * k, 2 * j + 1);
        double score = this->adjacency_matrix(2 * k, 2 * j) - this->adjacency_matrix(2 * k, 2 * j + 1);

        if (CHECK_SCORE) {
            if (((score > 0 && score < 6) || (score < 0 && score > -6))) {
    //                split_phased_blk(i);
    //                auto blk_no = chromo_phaser->variant_to_block_id
                auto break_idx = this->phasing_window->mat_idx2var_idx(j);
                this->setBlkIdx(break_idx);
            }
       }
    }
    for (auto i : phasing_window->current_window_idxes)
    {
        if (i == 561)
            int tmp = 9;
        mat_idx = phasing_window->var_idx2mat_idx(i);
        if (met_idx.find(mat_idx) == met_idx.end())
            met_idx.insert(mat_idx);
        else
            continue;
        // variant inside phasing block
        if (!variant_graph.contain(mat_idx))
            continue;

            // disjointed snp, do nothing
        else if (variant_graph.disjointedatpos(mat_idx))
            //phasing_window->insert_phased_block_starting_idx(i);
            ;

            // connected phased block
        else
        {
            std::set<uint> &variants_mat = variant_graph.connected_component[mat_idx];
            GMatrix sub_mat = this->slice_submat(variants_mat);
            CMatrix sub_count = this->slice_submat(variants_mat, true);
            if (variant_graph.fully_seperatable(mat_idx))
                find_connected_component(sub_count, variants_mat);
            else{
                int block_count = 0;
                std::map<uint, int> subroutine_map;
                std::map<uint, uint> subroutine_blk_start;
                std::map<uint,double> block_qualities;
                call_haplotype(sub_mat, variants_mat, block_count, subroutine_map, subroutine_blk_start, false, block_qualities);


                // merge these block
                if (block_count > 1)
                {
                    //std::cout << block_count << std::endl;
                    solver_subroutine(block_count, subroutine_map, subroutine_blk_start, block_qualities);
                }
            }

            block_no++;
        }
    }
    this->offset += this->chromo_phaser->intended_window_size;

//    TODO, if we need this
//    if (HAS_TENX)
//    {
//        for (auto start_idx: this->phased_block_starts)
//            barcode_aware_filter(start_idx.first);
//    }

//    split_phased_blk(2);

}

void Spectral::add_snp_edge_barcode_subroutine(ViewMap &sub_weighted_graph, CViewMap &sub_count_graph, VariantGraph &sub_variant_graph, std::map<uint, int> &subroutine_map, std::map<uint, uint> & subroutine_blk_start)
{
    for (auto &content : this->barcode_linker->barcodes)
    {
        std::map<uint, std::pair<double,double> > block_supporting_likelyhood;
        Barcode &barcode = content.second;
        for (auto &i : barcode.barcode_info)
        {
            if (!phasing_window->in_range(i.first))
                continue;
            auto a = this->phasing_window->var_idx2mat_idx(i.first);
            if (a < 0 or a >= this->variant_count)
                continue;
            uint blk_idx = this->phasing_window->mat_idx2var_idx(subroutine_blk_start[subroutine_map[a]]);
            ptr_PhasedBlock blk = this->phasing_window->blocks[blk_idx];

            if (blk->results.count(i.first) == 0)
                continue;
            if (block_supporting_likelyhood.count(blk_idx) == 0)
                block_supporting_likelyhood[blk_idx] = std::make_pair(0.0, 0.0);

            if ((i.second.first == 0 && blk->results[i.first]->is_REF()) ||
                (i.second.first == 1 && blk->results[i.first]->is_ALT())) {
                block_supporting_likelyhood[blk_idx].first += log10(i.second.second);
                block_supporting_likelyhood[blk_idx].second += log10(1 - i.second.second);
            } else {
                block_supporting_likelyhood[blk_idx].first += log10(1 - i.second.second);
                block_supporting_likelyhood[blk_idx].second += log10(i.second.second);
            }
        }

        for (auto &i : block_supporting_likelyhood) {
            auto a = this->phasing_window->var_idx2mat_idx(i.first);

            for (auto &j : block_supporting_likelyhood) {
                auto b = this->phasing_window->var_idx2mat_idx(j.first);
                if (a == b)
                    continue;
                double P_H1 = cal_score(i.second.first + j.second.first, i.second.second + j.second.second);
                double P_H2 = cal_score(i.second.first + j.second.second, i.second.second + j.second.first);
                if (P_H2 == P_H1)
                    continue;
                double score = abs(P_H1 - P_H2);
                int connection;
                P_H1 > P_H2 ? connection = 0 : connection = 1;


                sub_weighted_graph(2 * subroutine_map[a], 2 * subroutine_map[b]) += P_H1;
                sub_weighted_graph(2 * subroutine_map[a] + 1, 2 * subroutine_map[b] + 1) += P_H1;
                sub_weighted_graph(2 * subroutine_map[a], 2 * subroutine_map[b] + 1) += P_H2;
                sub_weighted_graph(2 * subroutine_map[a] + 1, 2 * subroutine_map[b]) += P_H2;

            }
        }
    }
}

void Spectral::add_snp_edge_subroutine(ViewMap &sub_weighted_graph, CViewMap &sub_count_graph, VariantGraph &sub_variant_graph, std::map<uint, int> & subroutine_map, std::map<uint, uint> & subroutine_blk_start, std::map<uint,double> &block_qualities)
{

    //fragment.read_qual = ceil(fragment.read_qual / 3);
    for (Fragment & fragment : frag_buffer) {
        if (fragment.snps.empty())
            continue;
        std::map<uint, std::pair<double,double> > block_supporting_likelyhood;
        for (auto &i : fragment.snps) {
            if (!phasing_window->in_range(i.first))
                continue;
            auto a = this->phasing_window->var_idx2mat_idx(i.first);
            if (a < 0 or a >= this->variant_count)
                continue;
            uint blk_idx = this->phasing_window->mat_idx2var_idx(subroutine_blk_start[subroutine_map[a]]);
            ptr_PhasedBlock blk = this->phasing_window->blocks[blk_idx];
            if (blk->results.count(i.first) == 0)
                continue;
            if (block_supporting_likelyhood.count(blk_idx) == 0)
                block_supporting_likelyhood[blk_idx] = std::make_pair(0.0, 0.0);

//            if ((i.second.first == 0 && blk->results[i.first]->is_REF()) ||
//                (i.second.first != 1 && blk->results[i.first]->is_ALT())) {
//                block_supporting_likelyhood[blk_idx].first += log10(i.second.second);
//                block_supporting_likelyhood[blk_idx].second += log10(1 - i.second.second);
//            } else {
//                block_supporting_likelyhood[blk_idx].first += log10(1 - i.second.second);
//                block_supporting_likelyhood[blk_idx].second += log10(i.second.second);
//            }

            if ((i.second.first == 0 && blk->results[i.first]->is_REF()) ||
                (i.second.first == 1 && blk->results[i.first]->is_ALT())) {
                if (fragment.getType() == FRAG_MATRIX) {
//                Here, beacause the matrix format with count 00(11) twice
                    block_supporting_likelyhood[blk_idx].first += i.second.second / 2;
                    block_supporting_likelyhood[blk_idx].second += 0;
                } else {
                    block_supporting_likelyhood[blk_idx].first += log10(i.second.second);
                    block_supporting_likelyhood[blk_idx].second += log10(1 - i.second.second);
                }
            } else {
                if (fragment.getType() == FRAG_MATRIX) {
                    block_supporting_likelyhood[blk_idx].first += 0;
                    block_supporting_likelyhood[blk_idx].second += i.second.second / 2;
                } else {
                    block_supporting_likelyhood[blk_idx].first += log10(1 - i.second.second);
                    block_supporting_likelyhood[blk_idx].second += log10(i.second.second);
                }
            }


        }

        for (auto &i : block_supporting_likelyhood)
        {
            auto a = this->phasing_window->var_idx2mat_idx(i.first);

            for (auto &j : block_supporting_likelyhood)
            {
                auto b = this->phasing_window->var_idx2mat_idx(j.first);
                if (a == b)
                    continue;
                double P_H1, P_H2;
                if (fragment.getType() == FRAG_MATRIX) {
                    if ((i.second.first != 0 && j.second.first != 0 ) || (i.second.second != 0 && j.second.second != 0)) {
                        P_H1 = (i.second.first + j.second.first + i.second.second + j.second.second) / 2;
                        P_H2 = 0;
                    } else {
                        P_H1 = 0;
                        P_H2 = (i.second.first + j.second.second + i.second.second + j.second.first) / 2;
                    }
                } else {
                    P_H1 = cal_score(i.second.first + j.second.first, i.second.second + j.second.second);
                    P_H2 = cal_score(i.second.first + j.second.second, i.second.second + j.second.first);
                }
                if (P_H2 == P_H1)
                    continue;
                double score = abs(P_H1 - P_H2);
                int connection;
                P_H1 > P_H2 ? connection = 0 : connection = 1;
                //connection == 0 ? sub_variant_graph.add_edge(subroutine_map[a], subroutine_map[b], true) : sub_variant_graph.add_edge(subroutine_map[a], subroutine_map[b], false);

                //sub_weighted_graph(2 * subroutine_map[a], 2 * subroutine_map[b] + connection) += score;
                //sub_weighted_graph(2 * subroutine_map[a] + 1, 2 * subroutine_map[b] + 1 - connection) += score;
                //sub_count_graph(2 * subroutine_map[a], 2 * subroutine_map[b] + connection) ++;
                //sub_count_graph(2 * subroutine_map[a] + 1, 2 * subroutine_map[b] + 1 - connection) ++;

                sub_weighted_graph(2 * subroutine_map[a], 2 * subroutine_map[b]) += P_H1;
                sub_weighted_graph(2 * subroutine_map[a] + 1, 2 * subroutine_map[b] + 1) += P_H1;
                sub_weighted_graph(2 * subroutine_map[a], 2 * subroutine_map[b] + 1) += P_H2;
                sub_weighted_graph(2 * subroutine_map[a] + 1, 2 * subroutine_map[b] ) += P_H2;

            }
        }
    }
}

void Spectral::solver_subroutine(int block_count, std::map<uint, int> & subroutine_map, std::map<uint, uint> &subroutine_blk_start, std::map<uint,double> &block_qualities)
{
    VariantGraph sub_variant_graph;
    sub_variant_graph.reset(block_count);
    int N = 2 * block_count;
    int *sub_count = new int [N * N];
    double *sub_weight = new double [N * N];

    for (int i = 0; i < N*N; i++)
    {
        sub_count[i] = 0;
        sub_weight[i] = 0.0;
    }

    ViewMap sub_weighed_graph(sub_weight, N, N);
    CViewMap sub_count_graph(sub_count, N, N);
    GMatrix weight_mat;
    CMatrix count_mat;
    if (HAS_TENX)
        add_snp_edge_barcode_subroutine(sub_weighed_graph, sub_count_graph, sub_variant_graph, subroutine_map, subroutine_blk_start);
    else
        add_snp_edge_subroutine(sub_weighed_graph, sub_count_graph, sub_variant_graph, subroutine_map, subroutine_blk_start, block_qualities);
    cal_prob_matrix(sub_weighed_graph, sub_count_graph, &weight_mat, &count_mat, &sub_variant_graph);

    std::unordered_set<uint> met_idx;
    for (uint i = 0; i < block_count; i++)
    {
        if (met_idx.find(i) == met_idx.end())
            met_idx.insert(i);
        if (!sub_variant_graph.contain(i))
            continue;
        else if (sub_variant_graph.disjointedatpos(i))
            continue;
        else
        {
            std::map<uint, int> new_subroutine_map;
            std::set<uint> &ori_variants_mat = sub_variant_graph.connected_component[i];
            std::set<uint> variants_mat;
            for (auto idx : ori_variants_mat)
                variants_mat.insert(subroutine_blk_start[idx]);
            GMatrix sub_mat = this->slice_submat(ori_variants_mat, weight_mat);
            CMatrix sub_count = this->slice_submat(ori_variants_mat, true, count_mat);
            if (sub_variant_graph.fully_seperatable(i))
                find_connected_component(count_mat, ori_variants_mat, sub_variant_graph, subroutine_blk_start);
            else
            {
                int new_block_count = 0;
                std::map<uint, uint> new_subroutine_blk_start;
                std::map<uint,double> _new_block_qualities;
                call_haplotype(sub_mat, variants_mat, new_block_count, new_subroutine_map, new_subroutine_blk_start, true, _new_block_qualities);
                if (new_block_count > 1)
                {
                    std::map<uint, int> _new_subroutine_map;
                    for (auto &j : subroutine_map)
                    {
                        _new_subroutine_map[j.first] = new_subroutine_map[subroutine_blk_start[j.second]];
                    }
                    solver_subroutine(new_block_count, _new_subroutine_map, new_subroutine_blk_start, _new_block_qualities);;
                }
            }

        }
    }

    delete [] sub_count;
    delete [] sub_weight;
    subroutine_blk_start.clear();
}

template<typename Derived>
void Spectral::find_connected_component_dfs(const Eigen::MatrixBase<Derived> &adj_mat, uint var_idx,
                                            std::unordered_map<uint, bool> &settled, ptr_PhasedBlock starting_block,
                                            bool prev_is_ref, uint prev_var_idx, VariantGraph &sub_variant_graph, std::map<uint, uint> &subroutine_blk_start)
{
    bool curr_is_ref = true;
    if (settled[var_idx])
        return;
        //not settled hap
    else
    {
        settled[var_idx] = true;
        std::set<uint> &nxt_vars = sub_variant_graph.graph[var_idx];

        uint idx = phasing_window->mat_idx2var_idx(subroutine_blk_start[var_idx]);
        ptr_PhasedBlock block_to_merge = phasing_window->blocks[idx];
        if (adj_mat(2 * prev_var_idx, 2 * var_idx) > 0)
        {
            if (block_to_merge->size() == 1)
            {
                block_to_merge->results[idx]->set_hap(prev_is_ref);
                block_to_merge->results[idx]->set_phased();
                curr_is_ref = block_to_merge->results[idx]->is_REF();
            }
            else
            {
                //prev_flipped
                if (!prev_is_ref)
                {
                    block_to_merge->flip();
                    curr_is_ref = false;
                }
            }
        }
        else if (adj_mat(2 * prev_var_idx + 1, 2 * var_idx) > 0)
        {
            if (block_to_merge->size() == 1)
            {
                block_to_merge->results[idx]->set_hap(!prev_is_ref);
                block_to_merge->results[idx]->set_phased();
                curr_is_ref = block_to_merge->results[idx]->is_REF();
            }
            else
            {
                //prev_flipped
                if (prev_is_ref)
                {
                    block_to_merge->flip();
                    curr_is_ref = false;
                }
            }
        }


        for (auto &results : block_to_merge->results)
            phasing_window->mat2variant_index[phasing_window->var_idx2mat_idx(results.first)] = starting_block->start_variant_idx;
        starting_block->join_block_no_overlap(block_to_merge);
        phasing_window->mat2variant_index[subroutine_blk_start[var_idx]] = starting_block->start_variant_idx;
        if (idx != starting_block->start_variant_idx)
            phasing_window->destroy_merged_block(idx);

        for (auto nxt_var : nxt_vars)
        {
            if (this->variant_graph.filtered[nxt_var])
                continue;
            if (!settled[nxt_var] && nxt_var != prev_var_idx)
                find_connected_component_dfs(adj_mat, nxt_var, settled, starting_block, curr_is_ref,
                                             var_idx, sub_variant_graph, subroutine_blk_start);
        }
    }
}

template<typename Derived>
void Spectral::find_connected_component_dfs(const Eigen::MatrixBase<Derived> &adj_mat, uint var_idx,
                                            std::unordered_map<uint, bool> &settled, ptr_PhasedBlock starting_block,
                                            bool prev_is_ref, uint prev_var_idx)
{
    // treat prev_is_ref as prev_not_flapped when calling recursive solver
    // already settled hap
    bool curr_is_ref = true;
    if (settled[var_idx])
        return;
        //not settled hap
    else
    {
        settled[var_idx] = true;
        std::set<uint> &nxt_vars = this->variant_graph.graph[var_idx];

        uint idx = phasing_window->mat_idx2var_idx(var_idx);
        ptr_PhasedBlock block_to_merge = phasing_window->blocks[idx];
        if (adj_mat(2 * prev_var_idx, 2 * var_idx) > 0)
        {
            if (block_to_merge->size() == 1)
            {
                block_to_merge->results[idx]->set_hap(prev_is_ref);
                block_to_merge->results[idx]->set_phased();
                curr_is_ref = block_to_merge->results[idx]->is_REF();
            }
            else
            {
                //prev_flipped
                if (!prev_is_ref)
                {
                    block_to_merge->flip();
                    curr_is_ref = false;
                }
            }
        }
        else if (adj_mat(2 * prev_var_idx + 1, 2 * var_idx) > 0)
        {
            if (block_to_merge->size() == 1)
            {
                block_to_merge->results[idx]->set_hap(!prev_is_ref);
                block_to_merge->results[idx]->set_phased();
                curr_is_ref = block_to_merge->results[idx]->is_REF();
            }
            else
            {
                //prev_flipped
                if (prev_is_ref)
                {
                    block_to_merge->flip();
                    curr_is_ref = false;
                }
            }
        }

        for (auto &results : block_to_merge->results)
            phasing_window->mat2variant_index[phasing_window->var_idx2mat_idx(results.first)] = starting_block->start_variant_idx;

        starting_block->join_block_no_overlap(block_to_merge);
         phasing_window->mat2variant_index[var_idx] = starting_block->start_variant_idx;
        if (idx != starting_block->start_variant_idx)
            phasing_window->destroy_merged_block(idx);

        for (auto nxt_var : nxt_vars)
        {
            if (this->variant_graph.filtered[nxt_var])
                continue;
            if (!settled[nxt_var] && nxt_var != prev_var_idx)
                find_connected_component_dfs(adj_mat, nxt_var, settled, starting_block, curr_is_ref,
                                             var_idx);
        }
    }
}


template<typename Derived>
void Spectral::find_connected_component(const Eigen::MatrixBase<Derived> &adj_mat, const std::set<uint> &variant_idx, VariantGraph &sub_variant_graph, std::map<uint, uint> &subroutine_blk_start)
{
    std::unordered_map<uint, bool> settled;
    for (auto i : variant_idx)
        settled[i] = false;

    auto it = variant_idx.begin();
    settled[*it] = true;
    uint prev_idx, start_idx;
    prev_idx = start_idx = phasing_window->mat_idx2var_idx(subroutine_blk_start[*it]);
    ptr_PhasedBlock phased_blk = phasing_window->blocks[start_idx];
    // not phased yet
    if (phased_blk->size() == 1)
    {
        phased_blk->results[start_idx]->set_phased();
        phased_blk->results[start_idx]->block = phased_blk;
    }
    std::set<uint> &nxt_vars = sub_variant_graph.graph[*it];
    for (auto nxt_var : nxt_vars)
    {
        if (sub_variant_graph.filtered[nxt_var])
            continue;
        if (phasing_window->blocks.count(phasing_window->mat_idx2var_idx(subroutine_blk_start[nxt_var])) > 0)
        {
            ptr_PhasedBlock block = phasing_window->blocks[phasing_window->mat_idx2var_idx(subroutine_blk_start[nxt_var])];
            if (block->size() == 1)
                find_connected_component_dfs(adj_mat, nxt_var, settled, phased_blk, true, *it, sub_variant_graph, subroutine_blk_start);
            else
                find_connected_component_dfs(adj_mat, nxt_var, settled, phased_blk, true, *it, sub_variant_graph, subroutine_blk_start);
        }
    }

    phasing_window->insert_phased_block_starting_idx(start_idx);
}

template<typename Derived>
void Spectral::find_connected_component(const Eigen::MatrixBase<Derived> &adj_mat, const std::set<uint> &variant_idx)
{
    Eigen::Ref<CMatrix> adjmat = count_matrix;
    std::unordered_map<uint, bool> settled;
    for (auto i : variant_idx)
        settled[i] = false;

    auto it = variant_idx.begin();
    settled[*it] = true;
    uint prev_idx, idx, start_idx, i;
    prev_idx = start_idx = phasing_window->mat_idx2var_idx(*it);
    ptr_PhasedBlock phased_blk = phasing_window->blocks[start_idx];
    // not phased yet
    if (phased_blk->size() == 1)
    {
        phased_blk->results[start_idx]->set_phased();
        phased_blk->results[start_idx]->block = phased_blk;
    }
    std::set<uint> &nxt_vars = this->variant_graph.graph[*it];
    for (auto nxt_var : nxt_vars)
    {
        if (this->variant_graph.filtered[nxt_var])
            continue;
        if (phasing_window->blocks.count(phasing_window->mat_idx2var_idx(nxt_var)) > 0)
        {
            ptr_PhasedBlock block = phasing_window->blocks[phasing_window->mat_idx2var_idx(nxt_var)];
            if (block->size() == 1)
                find_connected_component_dfs(adjmat, nxt_var, settled, phased_blk, true, *it);
            else
                find_connected_component_dfs(adjmat, nxt_var, settled, phased_blk, true, *it);
        }
    }

    phasing_window->insert_phased_block_starting_idx(start_idx);
}

//TODO: starting variant should not be filtered
int Spectral::locate_block_valid_start(const Eigen::VectorXd &vec)
{
    int i = 0;
    for (; i < (vec.rows() / 2); i++)
    {
        double diff =abs(abs(vec(2*i)) - abs(vec(2*i + 1)));
        diff = diff / abs(vec(2*i));
        //if (abs(vec(2 * i)) < threhold)
        //    continue;
        if (vec(2 * i) * vec(2*i + 1) < 0 && diff < 0.001)
            break;
    }
    return i;
}

//separate connected component according to Fielder Vector
void Spectral::separate_connected_component(const Eigen::VectorXd &vec, const std::set<uint> &variant_idx_mat)
{
    int valid_start_offset = locate_block_valid_start(vec);
    if (valid_start_offset >= vec.rows() / 2 - 2 && valid_start_offset > 0)
        return;
    auto it = variant_idx_mat.begin();
    std::advance(it, valid_start_offset);

    uint idx, start_idx, i;
    start_idx = phasing_window->mat_idx2var_idx(*it);
    this->phased_block_starts[start_idx] = *variant_idx_mat.begin();
    ptr_PhasedBlock phased_blk = phasing_window->blocks[start_idx];
    std::advance(it, 1);
    ViewMap weighted_graph(raw_graph, n, n);


    if (phased_blk->size() == 1)
    {
        phased_blk->results[start_idx]->block = phased_blk;
        phased_blk->results[start_idx]->set_phased();
    }

    //REF belongs to class 1
    bool ref_is_one = true;

    if (vec(2*valid_start_offset) < vec(2*valid_start_offset + 1))
        ref_is_one = false;

    i = valid_start_offset + 1;

    for (; it != variant_idx_mat.end(); ++it, i++)
    {
        idx = phasing_window->mat_idx2var_idx(*it);
        ptr_PhasedBlock block_to_merge = phasing_window->blocks[idx];
        //shall not use field with zero value
        auto tmp2 = abs(vec(2*i));
        auto tmp3 = abs(vec(2*i+1));
        if (abs(vec(2*i)) < threhold || abs(vec(2*i+1)) < threhold)
        {
            if (block_to_merge->size() == 1)
            {
                block_to_merge->results[idx]->set_filter(filter_type::POOLRESULT);
            }
            this->variant_graph.remove_variant(*it);
            continue;
        }
        //TODO: potential bug here， seperate the if condition for phased block
        if (block_to_merge->results[idx]->get_filter() == filter_type::PASS)
        {
            double diff =abs(abs(vec(2*i)) - abs(vec(2*i + 1)));
            diff = diff / abs(vec(2*i));
            //class 1
            if (vec(2 * i) > vec(2 * i + 1) && vec(2 * i) * vec(2 * i + 1) < 0 && diff < 0.001) // || diff > 0.5))
            {
                //if (su(2*i) >= -6)
                //{
                //    block_to_merge[idx].set_filter( filter_type::CONFILCTINGRESULT );
                //    phasing_window.insert_phased_block_starting_idx(idx);
                //}
                //else
                if (block_to_merge->size() == 1)
                {
                    block_to_merge->results[idx]->set_phased();
                    block_to_merge->results[idx]->set_hap(ref_is_one);
                }
                else
                {
                    if (!ref_is_one)
                        block_to_merge->flip();
                }
                //merge
                for (auto &results : block_to_merge->results)
                    phasing_window->mat2variant_index[phasing_window->var_idx2mat_idx(results.first)] = start_idx;
                phased_blk->join_block_no_overlap(block_to_merge);

                phasing_window->mat2variant_index[*it] = start_idx;
                if (idx != phased_blk->start_variant_idx)
                    phasing_window->destroy_merged_block(idx);

            }
                //class 2
            else if (vec(2 * i) < vec(2 * i + 1) && vec(2 * i) * vec(2 * i + 1) < 0 && diff < 0.001)//|| diff > 0.5))
            {
                //if (su(2*i) >= -6)
                //{
                //    block_to_merge[idx].set_filter( filter_type::CONFILCTINGRESULT );
                //    phasing_window.insert_phased_block_starting_idx(idx);
                //}
                //else
                if (block_to_merge->size() == 1)
                {
                    block_to_merge->results[idx]->set_phased();
                    block_to_merge->results[idx]->set_hap(!ref_is_one);
                }
                else
                {
                    if (ref_is_one)
                        block_to_merge->flip();

                }
                //merge
                for (auto &results : block_to_merge->results)
                    phasing_window->mat2variant_index[phasing_window->var_idx2mat_idx(results.first)] = start_idx;
                phased_blk->join_block_no_overlap(block_to_merge);
                phasing_window->mat2variant_index[*it] = start_idx;
                if (idx != phased_blk->start_variant_idx)
                    phasing_window->destroy_merged_block(idx);

            }
            else
            {
                if (block_to_merge->size() == 1)
                {
                    block_to_merge->results[idx]->set_filter(filter_type::POOLRESULT);
                }
                this->variant_graph.remove_variant(*it);
            }
        }
        else
            this->variant_graph.remove_variant(*it);
    }
}

void Spectral::barcode_aware_filter(uint block_start_idx)
{
    std::set<uint> filtered_idx;
    ptr_PhasedBlock phased_block = phasing_window->blocks[block_start_idx];
    // Notice that after phasing, every barcode provides support for at most one block
    BarcodeLinker::iterator iter;
    this->barcode_linker_index_set ? iter = this->barcode_linker->barcodes.find(this->barcode_index) : iter = this->barcode_linker->barcodes.begin();
    this->barcode_linker_index_set = false;
    for (; iter != this->barcode_linker->barcodes.end(); ++iter)
    {
        if (iter->second.end_var_idx < phased_block->start_variant_idx)
            continue;
        if (iter->second.start_var_idx > *phased_block->variant_idxes.rbegin())
        {
            this->barcode_index = iter->first;
            this->barcode_linker_index_set = true;
            break;
        }
        iter->second.get_barcode_support(phased_block);
    }
    phased_block->filter_inconsistency(filtered_idx);
    if (filtered_idx.size() > 5)
        this->block_tobe_split.insert(block_start_idx);

    for (auto i: filtered_idx)
        this->variant_graph.remove_variant(this->phasing_window->var_idx2mat_idx(i));
}

template<typename Derived>
bool Spectral::cal_fiedler_vec(int nev, const Eigen::MatrixBase<Derived> &adj_mat, GMatrix &vecs, Eigen::VectorXd &vals)
{
    int size = adj_mat.cols();
    const GMatrix &D = adj_mat.colwise().sum().asDiagonal();
    const GMatrix &L = D - adj_mat;
    bool converge = true;

    if (nev >= size - 1)
        nev = size - 1;

    Arpack arpack;
    arpack.compute(L,  nev, "SM");
    if (arpack.info() != Eigen::Success)
        converge = false;
    else
    {
        vecs = arpack.eigenvectors();
        vals = arpack.eigenvalues();
        //std::cout << pow(size, 0.5) * vals(1) << std::endl;
    }

    return coverage;
}

//TODO: check whether we can directly yield the haplotype in this function. urgent! this is a bug
//TODO: refine code for the logics to work same for Hi-C NGS, 10X, Nanopore and pacbio 
void Spectral::call_haplotype(GMatrix &adj_mat, const std::set<uint> &variant_idx_mat, int& block_count, std::map<uint, int> &subroutine_map, std::map<uint, uint> &subroutine_blk_start, bool sub, std::map<uint, double> &block_quality)
{
    if (variant_idx_mat.size() == 1)
        return;
    int nev = 3;
    GMatrix vecs;
    Eigen::VectorXd vals;
    bool converge = cal_fiedler_vec(nev, adj_mat, vecs, vals);
    std::map<uint, uint> reverse_map;
    if (sub)
    {
        int c = 0;
        for (auto i: variant_idx_mat)
            reverse_map[i] = c++;
    }
    // not converge, rare case since we are solving Lrm
    if (!converge)
    {
        for (auto i : variant_idx_mat)
        {
            uint idx = phasing_window->mat_idx2var_idx(i);
            ptr_PhasedBlock phased_block = phasing_window->blocks[idx];
            if (phased_block->size() == 1)
                phased_block->results[idx]->set_filter(filter_type::POOLRESULT);
        }
    }
        //converge
    else
    {
        uint fiedler_idx = 1;
        GMatrix tmp = adj_mat;
        if (vecs.cols() < 1)
            return;
        if (trivial_fiedler(vecs.col(fiedler_idx)))
        {
            std::set<uint> p_var, n_var;
            fiedler_guided_group_count(&p_var, &n_var, vecs.col(fiedler_idx), variant_idx_mat);
            if (p_var.empty() || n_var.empty()) //failed fiedler vector, use secondary,
            {
                fiedler_idx = 2;
                if (trivial_fiedler(vecs.col(fiedler_idx)))
                    ;//std::cout << "fail to solve fielder vector at \n"; //failed we seperate block by doing nothing.
                //else
                //   separate_connected_component(vecs.col(fiedler_idx), variant_idx_mat);
            }
            else  //block to be cut, enter recursive solver
            {
                if (p_var.size() > 1)
                {
                    GMatrix l_submat;
                    std::set<uint> sub_p_var;
                    if (sub)
                    {
                        for (auto j : p_var)
                            sub_p_var.insert(reverse_map[j]);
                        l_submat = this->slice_submat(sub_p_var, adj_mat);
                    }
                    else
                        l_submat = this->slice_submat(p_var);
                    call_haplotype(l_submat, p_var, block_count, subroutine_map, subroutine_blk_start, sub, block_quality);
                }
                if (n_var.size() > 1)
                {
                    GMatrix r_submat;
                    std::set<uint> sub_n_var;
                    if (sub)
                    {
                        for (auto j : n_var)
                            sub_n_var.insert(reverse_map[j]);
                        r_submat = this->slice_submat(sub_n_var, adj_mat);
                    }
                    else
                        r_submat = this->slice_submat(n_var);
                    call_haplotype(r_submat, n_var, block_count, subroutine_map, subroutine_blk_start, sub, block_quality);
                }
            }

        }
        else
        {
            separate_connected_component(vecs.col(fiedler_idx), variant_idx_mat);
            //if (!sub)
            //    poss_phase_error_correction(phasing_window->mat_idx2var_idx(*variant_idx_mat.begin()));
            block_quality[block_count] = log10(vals(1) *  pow(vecs.rows(), 0.5));
            for (auto & idx : variant_idx_mat)
                subroutine_map.emplace(idx, block_count);
            subroutine_blk_start[block_count] =  *variant_idx_mat.begin();
            block_count++;
        }

    }
}

void Spectral::fragment_supported_flipping_score(ptr_PhasedBlock &phased_block, Fragment &fragment, int *supporting_reads_count, double *supporting_weight_count, std::map<uint, std::set<uint>> &connection_map)
{
    int supporting_snp_count = 0;
    int total_snp_count = 0;
    int count = 0;
    std::map<uint, int> idx_map;
    for (auto &i : phased_block->results)
        idx_map[i.first] = count++;

    for (auto &i : fragment.snps)
    {
        if (!phased_block->contain_variant(i.first))
            continue;
        if (connection_map.count(i.first) == 0)
            connection_map[i.first];
        for (auto &j : fragment.snps)
        total_snp_count++;
        if ((i.second.first == 0 && phased_block->results[i.first]->is_REF()) || (i.second.first == 1 && phased_block->results[i.first]->is_ALT()))
            supporting_snp_count ++;
    }
    if (total_snp_count <= 1)
        return;
    double supporting_score = double(supporting_snp_count) / total_snp_count;
    for (auto &i : fragment.snps)
    {
        if (!phased_block->contain_variant(i.first))
            continue;
        if ((i.second.first == 0 && phased_block->results[i.first]->is_REF()) || (i.second.first == 1 && phased_block->results[i.first]->is_ALT()))
        {
            supporting_weight_count[idx_map[i.first]] += 1 - supporting_score;    //  flipping
            //supporting_weight_count[2*idx_map[i.first] + 1] += supporting_score;    // not flipping
            supporting_reads_count[idx_map[i.first]] ++;
        }
        else{
            supporting_weight_count[idx_map[i.first]] += supporting_score;        // flipping
            //supporting_weight_count[2*idx_map[i.first] + 1] += 1 - supporting_score;// not flipping
            supporting_reads_count[idx_map[i.first]] ++;
        }
    }

}

void Spectral::poss_phase_error_correction(uint block_start_idx)
{
    ptr_PhasedBlock phased_blk = this->phasing_window->blocks[block_start_idx];
    int block_snp_count = phased_blk->results.size();
    int *supporting_reads_count = new int [block_snp_count];
    double *supporting_weight_count = new double [block_snp_count]; //flipping/not flipping
    double *supporting_score = new double[block_snp_count];
    double **flip_penalty = new double *[block_snp_count];
    int count = 0;
    std::map<uint, int> idx_map;
    for (auto &i : phased_blk->results)
        idx_map[i.first] = count++;
    std::map<uint, std::set<uint>> connection_map;
    for (int i = 0; i < block_snp_count; i++)
    {
        supporting_reads_count[i] = 0;
        supporting_weight_count[i] = 0.0;
        //supporting_weight_count[2 * i + 1] = 0.0;
        supporting_score[i] = 0.0;
        //supporting_score[2*i + 1] = 0.0;
        flip_penalty[i] = new double[block_snp_count];
    }
    for (int i = 0; i < block_snp_count; i++)
        for (int j = 0; j < block_snp_count; j++)
            flip_penalty[i][j] = 0.0;

    for (Fragment &fragment : frag_buffer)
    {
        if (fragment.end < block_start_idx)
            continue;
        fragment_supported_flipping_score(phased_blk, fragment, supporting_reads_count, supporting_weight_count, connection_map);
    }

    for (int i = 0; i < block_snp_count; i++)
    {
        supporting_weight_count[i] /= supporting_reads_count[i];
        //supporting_weight_count[2 * i + 1] /= supporting_reads_count[i];
    }

    for (auto &i : phased_blk->results)
    {
        i.second->complement_hap();
        for (int j = 0; j < block_snp_count; j++)
            supporting_reads_count[j] = 0;
        for (Fragment &fragment : frag_buffer)
        {
            if (fragment.end < block_start_idx)
                continue;
            fragment_supported_flipping_score(phased_blk, fragment, supporting_reads_count, flip_penalty[idx_map[i.first]], connection_map);
        }
        for (int j = 0; j < block_snp_count; j++)
        {
            flip_penalty[idx_map[i.first]][j] /= supporting_reads_count[j];
            //supporting_weight_count[2 * i + 1] /= supporting_reads_count[i];
        }
        i.second->complement_hap();
    }
    double *penalty = new double[block_snp_count];
    for (int i = 0; i < block_snp_count; i++)
    {
        penalty[i] = 0.0;
        for (int j = 0; j < block_snp_count; j++)
        {
            if (i == j)
                continue;

            penalty[i] += flip_penalty[i][j];
            //no_flipping_weight_count[j] += 1.0 - flip_penalty[i][j];
        }
    }

    delete[] supporting_reads_count;
    delete[] supporting_weight_count;
    delete[] supporting_score;
    delete[] penalty;
    for (int i = 0; i < block_snp_count; i++)
        delete [] flip_penalty[i];
    delete [] flip_penalty;
}

void Spectral::split_phased_blk(uint idx)
{
    ptr_PhasedBlock phased_blk = this->phasing_window->blocks[idx];
    std::vector<std::set<uint>> split_var_idx;
    this->variant_graph.split_after_filtering(phasing_window->var_idx2mat_idx(idx), split_var_idx);
    this->phasing_window->split_after_phasing(idx, split_var_idx);
}

std::unordered_map<uint, std::set<uint>> Spectral::load_hic_poss_info()
{
    std::set<uint>var_idx;
    std::map<uint, uint> var_idx_map;
    for (auto &linker : hic_linker_container.linker)
    {
        for (auto &info : linker.second.hic_info)
        {
            if (info.first >= this->chromo_phaser->results_for_variant.size())
                continue;
            ptr_ResultforSingleVariant variant =  this->chromo_phaser->results_for_variant[info.first];
            if (is_uninitialized(variant->block))
            {
                var_idx.insert(info.first);
                var_idx_map[info.first] = info.first;
            }
            else
            {
                uint start_var_idx = variant->block.lock()->start_variant_idx;
                var_idx.insert(start_var_idx);
                var_idx_map[info.first] = start_var_idx;
            }
        }
    }

    std::unordered_map<uint, uint> var2idx;
    std::unordered_map<uint, uint> idx2var;
    std::unordered_map<uint, uint> start_var2idx;
    int count = 0;
    
    for (auto &i : var_idx)
    {
        idx2var[count] = i;
        start_var2idx[i] = count++;
    }

    count = 0;
    for (auto &i : var_idx_map)
    {
        var2idx[i.first] = start_var2idx[i.second];
    }

    // now use disjoint set to find connected component
    std::vector<uint> parent;
    std::vector<uint> size;
    

    parent.resize(var_idx.size());
    size.resize(var_idx.size());
    for (int i = 0; i < var_idx.size(); i++)
        this->make_set(i, parent, size);

    for (auto &linker : hic_linker_container.linker)
    {
        uint linker_start = var2idx[linker.second.start_var_idx];
        for (auto &info: linker.second.hic_info)
        {
            uint linker_idx = var2idx[info.first];
            if (linker_idx == linker_start)
                continue;
            if (this->find_set(linker_start, parent) != this->find_set(linker_idx, parent))
                this->union_sets(linker_start, linker_idx, parent, size);
        }
    }

    std::set<uint> comps_ids;
    std::unordered_map<uint, std::set<uint>> comps;
    for (int i = 0; i < var_idx.size(); i++)
    {
        if (comps.count(parent[i]) == 0)
        {
            comps[parent[i]] = std::set<uint>();
            comps_ids.insert(parent[i]);
        }
        comps[parent[i]].insert(idx2var[i]);
    }
    return comps;
}


void Spectral::load_hic_linker(int nblock)
{
    // variant_count = curr_window_block_count
    
    // start_variant_idx = starting_block's start_variant_idx
    // end_variant_idx_overlap = last block's start_variant_idx
    this->clean();
    this->variant_count = nblock;
    this->n = 2 * variant_count;
    this->variant_graph.reset(variant_count);
    this->raw_graph = new double[n * n];
    this->raw_count = new int[n * n];

    for (unsigned int i = 0; i < n * n; i++)
    {
        this->raw_graph[i] = 0.0;
        this->raw_count[i] = 0;
    }
    ViewMap weighed_graph(raw_graph, n, n);
    CViewMap count_graph(raw_count, n, n);
    add_snp_edge_hic(weighed_graph, count_graph);
    cal_prob_matrix(weighed_graph, count_graph, nullptr, nullptr, nullptr);
}


void Spectral::hic_poss_solver(int nblock)
{
    load_hic_linker(nblock);
        //main driven code 

    Eigen::Ref<GMatrix> adj_mat = adjacency_matrix;
    uint mat_idx;
        //for each block (variant) in phasing window
    std::unordered_set<uint> met_idx;
    for (auto i : phasing_window->current_window_idxes)
    {
        mat_idx = phasing_window->var_idx2mat_idx(i);
        if (met_idx.find(mat_idx) == met_idx.end())
            met_idx.insert(mat_idx);
        else
            continue;
        // variant inside phasing block
        if (!variant_graph.contain(mat_idx))
            continue;

                // disjointed snp, do nothing
        else if (variant_graph.disjointedatpos(mat_idx))
                //phasing_window->insert_phased_block_starting_idx(i);
            ;

                // connected phased block
        else
        {
            std::set<uint> &variants_mat = variant_graph.connected_component[mat_idx];
            Eigen::ArrayXi index(variants_mat.size() * 2);
            uint count = 0;
            for (auto j = variants_mat.begin(); j != variants_mat.end(); ++j, count++)
            {
                index(2 * count) = 2 * (*j);
                index(2 * count + 1) = 2 * (*j) + 1;
            }
            GMatrix sub_mat = mat_indexing(adj_mat, index, index);
            if (variant_graph.fully_seperatable(mat_idx))
                find_connected_component(sub_mat, variants_mat);
            else
            {
                int block_count = 0;
                std::map<uint, int> subroutine_map;
                std::map<uint, uint> subroutine_blk_start;
                std::map<uint,double> block_qualities;
                call_haplotype(sub_mat, variants_mat, block_count, subroutine_map, subroutine_blk_start, false, block_qualities);
            }
            block_no++;
        }
    }
    delete[] this->raw_graph;
    delete[] this->raw_count;
    this->raw_count = nullptr;
    this->raw_graph = nullptr;
}

const std::set<uint> &Spectral::getBreakIdxs() const {
    return break_idxs;
}

uint Spectral::getOffset() const {
    return offset;
}

void Spectral::setOffset(uint offset) {
    Spectral::offset = offset;
}

