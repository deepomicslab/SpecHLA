//
// Created by yonghanyu2 on 8/10/2018.
//
#ifndef SPECHAP_SPECTRAL_H
#define SPECHAP_SPECTRAL_H

#include <map>
#include <set>
#include "frag_io.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/ArpackSupport>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include "type.h"
#include "tenx_util.h"
#include "hic_util.h"
#include "graph.h"
#include "util.h"
//------------------------------------------------alias----------------------------------------------
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> GMatrix;
typedef Eigen::Map<GMatrix> ViewMap;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> CMatrix;
typedef Eigen::Map<CMatrix> CViewMap;
typedef Eigen::FullPivLU<GMatrix> denseChol;
typedef Eigen::ArpackGeneralizedSelfAdjointEigenSolver<GMatrix , denseChol> Arpack;

//----------------------------------------------------------------------------------------------------
extern int MAX_BARCODE_SPANNING;
extern int OPERATION;
extern std::vector<int> OPERATIONS;
extern bool HAS_HIC;
extern bool HAS_TENX;
extern bool CHECK_SCORE;

inline double cal_score(double a, double b)
{
    double score;
    if (a > b)
    {
        score = a + log10(1 + pow(10, b - a)) - log10(2);
    }
    else
    {
        score = b + log10(1 + pow(10, a - b)) - log10(2);
    }
    return score;
}

//---------------------------------------solver class------------------------------------------------
class Spectral
{
private:
    bool use_secondary;
    double *raw_graph;          //raw data to hold the Graph
    int *raw_count;
    GMatrix adjacency_matrix;
    CMatrix count_matrix;
    uint offset;
public:
    uint getOffset() const;

    void setOffset(uint offset);

private:
    //offset for current phasing window
    double threhold;                //precision, to be calculated from read depth                         //operation mode, PE, 10X or HiC
    uint variant_count{};                 //number of point
    uint start_variant_idx{};
    uint end_variant_idx_intended{};
    uint end_variant_idx_overlap{};
    uint n{};                             //size of array
    double epsilon;
    int coverage;
    int max_barcode_spanning_length{};
    std::vector<Fragment> frag_buffer;
    double q_aver;
    double q_sum;

    bool has_hic{};
    bool has_tenx{};
    std::set<uint> break_idxs;
    std::vector<double> fr_weights;
public:
    const std::set<uint> &getBreakIdxs() const;

private:

    //-----------------pack these into another classes------------------------
    BarcodeLinkers * barcode_linker;
    bool barcode_linker_index_set;
    std::string barcode_index;
    RegionFragStats *region_frag_stats;
    //------------------------------------------------------------------------

    std::map<uint, uint> phased_block_starts;
    std::set<uint> block_tobe_split;
//    FragmentReader *fr;             //Fragment matrix reader
    std::vector<FragmentReader* > frs;
    BEDReader *frbed;
    HiCLinkerContainer hic_linker_container;
    VariantGraph variant_graph;
    uint block_no;                  //start from 1
    ptr_PhasingWindow phasing_window;
    ChromoPhaser *chromo_phaser;


public:
    Spectral(std::vector<FragmentReader *>& frs, std::vector<double>& fr_weights, BEDReader *frbed, double threhold, int coverage, bool use_secondary);
    void clean();
    void reset();
    void set_prev_buff(ViewMap& weighted_graph, CViewMap& count_graph);
    void solver();
    void hic_poss_solver(int nblock);
    std::unordered_map<uint, std::set<uint>> load_hic_poss_info();
    void solver_recursive();
    void solver_subroutine(int block_count, std::map<uint, int> &subroutine_map, std::map<uint, uint>& subroutine_blk_start, std::map<uint,double> &block_qualities);
    inline void set_chromo_phaser(ChromoPhaser *chromo_phaser)
    {
        this->chromo_phaser = chromo_phaser;
        this->barcode_linker->chromo_phaser = chromo_phaser;
        this->phasing_window = chromo_phaser->phased;
    }
    inline void release_chromo_phaser()
    {
        this->clean();
        this->chromo_phaser = nullptr;
        this->phasing_window = nullptr;
        this->barcode_linker->chromo_phaser = nullptr;
        this->hic_linker_container.clear();
    }

    ~Spectral();
    uint setBlkIdx(uint k) {
        this->break_idxs.emplace(k);
    }

private:
    void add_snp_edge(Fragment &fragment, ViewMap &weighted_graph, CViewMap &count_graph, double w);
//    void add_snp_edge_matrix(Fragment &fragment, ViewMap &weighted_graph, CViewMap &count_graph, double w);
    void add_snp_edge_barcode(ViewMap &weighted_graph, CViewMap &count_graph);
    void add_snp_edge_hic(ViewMap &weighted_graph, CViewMap &count_graph);
    void add_snp_edge_subroutine(ViewMap &sub_weighted_graph, CViewMap &sub_count_graph, VariantGraph & sub_variant_graph, std::map<uint, int> & subroutine_map, std::map<uint, uint> & subroutine_blk_start, std::map<uint,double> &block_qualities);
    void add_snp_edge_barcode_subroutine(ViewMap &sub_weighted_graph, CViewMap &sub_count_graph, VariantGraph & sub_variant_graph, std::map<uint, int> & subroutine_map, std::map<uint, uint> & subroutine_blk_start);
    void add_snp_edge_tgs();
    GMatrix slice_submat(std::set<uint> &variants_mat);
    GMatrix slice_submat(std::set<uint> &variants_mat, GMatrix &adj_mat);
    CMatrix slice_submat(std::set<uint> &variants_mat, bool t);
    CMatrix slice_submat(std::set<uint> &variants_mat, bool t, CMatrix &adj_mat);
    void filter_inconsistency();
    void read_fragment_10x(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph);
    void read_fragment(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph, double w);
    void read_fragment_hic(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph);
    void read_fragment_pacbio(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph, double w);
    void read_fragment_matrix(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph, double w);
    void read_fragment_nanopore(int frIdx, ViewMap &weighted_graph, CViewMap &count_graph);
    void poss_phase_error_correction(uint block_start_idx);
    void fragment_supported_flipping_score(ptr_PhasedBlock &phased_block, Fragment & fragment, int *supporting_reads_count, double *supporting_weight_count, std::map<uint, std::set<uint>> &connection_map);
    int locate_block_valid_start(const Eigen::VectorXd &vec);


    template <typename Derived>
    void find_connected_component(const Eigen::MatrixBase<Derived> &adj_mat, const std::set<uint> &variant_id, VariantGraph &sub_variant_graph, std::map<uint, uint> &subroutine_blk_start);
    template <typename Derived>
    void find_connected_component(const Eigen::MatrixBase<Derived> &adj_mat, const std::set<uint> &variant_idx);
    template <typename Derived>
    void find_connected_component_dfs(const Eigen::MatrixBase<Derived> &adj_mat, uint var_idx, std::unordered_map<uint ,bool> &settled, ptr_PhasedBlock starting_block, bool prev_is_ref, uint prev_var_idx);
    template <typename Derived>
    void find_connected_component_dfs(const Eigen::MatrixBase<Derived> &adj_mat, uint var_idx, std::unordered_map<uint ,bool> &settled, ptr_PhasedBlock starting_block, bool prev_is_ref, uint prev_var_idx, VariantGraph &sub_variant_graph, std::map<uint, uint> &subroutine_blk_start);
    void separate_connected_component(const Eigen::VectorXd &vec, const std::set<uint> &variant_idx_mat);

    void call_haplotype(GMatrix &adj_mat, const std::set<uint> &variant_idx_mat, int &block_count, std::map<uint, int> &subroutine_map, std::map<uint, uint> &subroutine_blk_start, bool sub, std::map<uint, double> &block_quality);
    //void call_haplotype_hic_link(const Eigen::MatrixBase<Derived> &adj_mat, HapStruct &phased_blocks, const unsigned int &idx);
    void load_hic_linker(int nblock);
    void cal_prob_matrix(ViewMap &weighted_graph, CViewMap &count_graph, GMatrix *weight, CMatrix *count, VariantGraph *variant_graph);
    void merge_hap_block();
    template <typename Derived>
    bool cal_fiedler_vec(int nev, const Eigen::MatrixBase<Derived> &adj_mat, GMatrix &vecs, Eigen::VectorXd & vals);
    void barcode_aware_filter(uint block_start_idx);
    bool use_vec(const Eigen::VectorXd &first);
    void split_phased_blk(uint index);
    //inline func


    inline int count_0_eigen_val(const Eigen::VectorXd &eigen_vals)
    {
        int count = 0;
        for (int i = 0; i < eigen_vals.size(); i++)
        {
            const double &val = eigen_vals(i);
            if (val < threhold and val > -1 * threhold)
                count++;
        }
        return count;
    }



    inline double get_qual_barcode(const char &qa, const char &qb)
    {
        if(qa > qb)
            return (int(qb)-33) / (-10.0);
        return (int(qa) -33) / (-10.0);
        //return (int(qa) - 66 + int(qb)) / (-10);
    }

    inline bool equal(double first, double second){return abs(first - second) < pow(10, -8);}
    inline void cal_10x_filter(ptr_PhasedBlock &phased_block)
    {
        if (phased_block->size() > 1)
            return;
        uint var_idx = phased_block->start_variant_idx;
        ptr_ResultforSingleVariant result = phased_block->results[var_idx];
        if (result->filtered)
            return;
        if (result->alt_count < 2 || (double(result->alt_count) / (result->alt_count + result->ref_count)) < 0.15)
        {
            result->set_filter(filter_type::TENX_ALLELE_FRENCENCY_FILTER);
            this->variant_graph.remove_variant(this->phasing_window->var_idx2mat_idx(var_idx));
        }
        if (result->qual < 15 || (result->af > 0.5 && result->qual < 50))
        {
            if (result->get_filter() == filter_type::PASS)
                result->set_filter(filter_type::TENX_QUAL_FILTER);
            this->variant_graph.remove_variant(this->phasing_window->var_idx2mat_idx(var_idx));
        }
        if (double(result->ad1) / (result->ad0 + result->ad1) < 0.15)
        {
            result->set_filter(filter_type::TENX_ALLELE_FRENCENCY_FILTER);
            this->variant_graph.remove_variant(this->phasing_window->var_idx2mat_idx(var_idx));
        }
        result->nMMD > 0 ? result->MMD /= result->nMMD : result->MMD = -1;
        if (result->rescued + result->not_rescued > 0)
        {
            if (double (result->rescued) / (result->rescued + result->not_rescued) > 0.1 && (result->MMD == -1 || result->MMD >= 3))
            {
                result->set_filter(filter_type::TENX_RESCUED_MOLECUE_HIGH_DIVERSITY);
                this->variant_graph.remove_variant(this->phasing_window->var_idx2mat_idx(var_idx));
            }
        }
        result->filtered = true;
        /*
        if (result->dp < (this->coverage / 5 * 2) && result->get_filter() == filter_type::PASS)
        {
            result->set_filter(filter_type::LOW_COVERAGE);
            this->variant_graph.remove_variant(this->phasing_window->var_idx2mat_idx(var_idx));
        }
        */
    }
    inline void get_current_window_pos(uint &start_pos, uint &end_pos)
    {
        uint start_ix = this->phasing_window->current_window_idxes[0];
        uint end_ix = this->phasing_window->current_window_idxes.back();
        start_pos = this->chromo_phaser->get_var_pos(start_ix);
        end_pos = this->chromo_phaser->get_var_pos(end_ix);
        const uint offset = 300000;
        if ( start_pos > offset)
            start_pos -= offset;
        else
            start_pos = 1;
        end_pos += 1;
    }

    inline void populate_variant_info(const Fragment &fragment)
    {
        for (auto &i :fragment.snps)
        {
            if (i.first >= this->chromo_phaser->variant_count)
                continue;
            ptr_ResultforSingleVariant snp = this->chromo_phaser->results_for_variant[i.first];
            if (fragment.read_qual <= -3)
            {
                if (i.second.first == 0)
                    snp->ref_count ++;
                else  if (i.second.first == 1)
                    snp->alt_count ++;
            }
            if (fragment.rescued)
                snp->rescued ++;
            else
                snp->not_rescued ++;
            if (i.second.first == 1)
            {
                if (fragment.dm > 0.001)
                {
                    snp->nMMD++;
                    snp->MMD += fragment.dm;
                }

            }
        }

    }

    inline bool trivial_fiedler (const Eigen::VectorXd & vec)
    {
        return vec(0) * vec(1) >= 0 && vec(2) * vec(3) >= 0;
    }

    inline void fiedler_guided_group_count(std::set<uint> *pvar, std::set<uint> *nvar, const Eigen::VectorXd &vec, const std::set<uint> &variant_idx_mat)
    {
        size_t count = 0;

        if (vec(0) < 0)
        {
            auto * tmp = nvar;
            nvar = pvar;
            pvar = tmp;
            tmp = nullptr;
        }

        for (auto i : variant_idx_mat)
        {
            if (vec(2 * count) > 0 && vec(2 * count + 1) > 0)
                pvar->insert(i);
            if (vec(2 * count) < 0 && vec(2 * count + 1) < 0)
                nvar->insert(i);
            count++;
        }
    }

    // function for disjoint set implementation
    inline void make_set(uint v, std::vector<uint> &parent, std::vector<uint> &size) 
    {
        parent[v] = v;
        size[v] = 1;
    }

    inline int find_set(uint v, std::vector<uint> &parent) 
    {
        if(parent[parent[v]] != parent[v])
            parent[v]= find_set(parent[v], parent);
        return parent[v];
        //if (v == parent[v])
        //    return v;
        //return parent[v] = find_set(parent[v], parent);
    }

    void union_sets(uint a, uint b, std::vector<uint> &parent, std::vector<uint> &size) 
    {
        a = find_set(a, parent);
        b = find_set(b, parent);
        if (a != b) 
        {
            if (size[a] < size[b])
                swap(a, b);
            parent[b] = a;
            size[a] += size[b];
        }
    }

    inline void swap(uint &a, uint &b)
    {
        int tmp = a;
        a = b;
        b = tmp;
    }
};


#endif //SPECHAP_SPECTRAL_H
