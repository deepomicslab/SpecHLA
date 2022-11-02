
// Note that the index maintains different meaning in different scope
// Variant idx: idx of variants on chromosome according to its order in VCF file, 0 based
// Variant idx in fragment: order according to VCF file, 1 based
// Block idx: always start from 0
// Block Matrix idx: same as above

#ifndef SPECHAP_TYPE_H
#define SPECHAP_TYPE_H

#include "results.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <string>
#include <cmath>
#include <map>
typedef std::pair<int, double> allele_info;
typedef std::pair<uint, allele_info> snp_info;
typedef std::vector<snp_info> snp_container;

#define FRAG_HIC    0
#define FRAG_10X    1
#define FRAG_NORMAL 2
#define FRAG_MATRIX 3

#define MODE_10X        0
#define MODE_HIC        1
#define MODE_PE         2
#define MODE_PACBIO     3
#define MODE_NANOPORE   4
#define MODE_MATRIX     5
#define MODE_HYBRID     6
class Fragment
{
public:
    Fragment() = default;
    Fragment(const Fragment &rhs);
    ~Fragment() = default;
    uint start, end;
    int type = FRAG_NORMAL;

    int getType() const;

    void setType(int type);

    snp_container snps;         //<variant_idx, allele>
    std::string barcode;
    double read_qual = 0;
    int insertion_size = 0;
    bool rescued = false;
    double dm = 0;
    inline void reset() { this->snps.clear(); }
    inline void insert(snp_info &snp) { snps.push_back(snp); } //potential bug here
    inline void update_start_end()
    {
        start = this->snps[0].first;
        end = this->snps.back().first;
    }
};

typedef std::shared_ptr<ResultforSingleVariant> ptr_ResultforSingleVariant;

template class std::unordered_map<uint, ptr_ResultforSingleVariant>;

class PhasedBlock : public std::enable_shared_from_this<PhasedBlock>
{
public:
    uint start_variant_idx;                                             //variant idx
    std::unordered_map<uint, ptr_ResultforSingleVariant> results;           //key variant idx
    std::set<uint> variant_idxes;
    uint block_id = 0;                                                        //for output

public:///=//public:

    PhasedBlock();
    PhasedBlock(uint start_idx);
    PhasedBlock(uint start_idx, ptr_ResultforSingleVariant result);
    PhasedBlock(const PhasedBlock &rhs);
    PhasedBlock(const PhasedBlock &rhs, std::unordered_map<uint, ptr_ResultforSingleVariant> &results_dup);
    ~PhasedBlock();

    //getter
    inline std::size_t size() { return  variant_idxes.size(); }
    inline uint get_start_position() { return start_variant_idx; }
    //setter
    void flip();                                                          //flip the haplotype of the block


    //utility
    void fragment_supported_score(Fragment &fragment);
    void join_block_no_overlap(std::shared_ptr<PhasedBlock> rhs);                       //merge the content of another non-overlapping block
    void join_block_overlap(std::shared_ptr<PhasedBlock> rhs);                          //merge the content of overlapping block
    inline bool contain_variant(uint idx) { return results.count(idx) != 0; }
    void filter_inconsistency(std::set<uint> &filtered_idx);
    inline void insert_variant(uint variant_idx, ptr_ResultforSingleVariant result)
    {
        results[variant_idx] = result;
        variant_idxes.insert(variant_idx);
    }
    std::string get_hap();
};


typedef std::shared_ptr<PhasedBlock> ptr_PhasedBlock;
template class std::unordered_map<uint, ptr_PhasedBlock>;

//this class should contain infomation about overlap region
class PhasingWindow
{
public:
    std::unordered_map<uint, std::shared_ptr<PhasedBlock> > blocks;       //key stored in the block idxes
    uint overlap_length;                                //blk count of overlap region
    uint intended_window_length;                                 //window length with overlap
    // all the blocks are indexed by starting variant index, listed below
    std::unordered_map<uint, uint> variant2mat_index;
    std::unordered_map<uint, uint> mat2variant_index;
    std::vector<uint> prev_block_idxes;                    //the starting variant index of blocks, by inserting order, used before phasing
    std::map<uint, uint> block_idxes;                      //the starting variant index of blocks, by inserting order, initialized after phasing
    std::vector<uint> current_window_idxes;
    std::map<uint, int> current_window_overlap_region_block_size;
    uint rest_blk_count, phased_blk_count, phasing_block_size;
    uint total_blk_count;
    uint intended_start, start, end, intended_end, prev_window_start;


public:
    PhasingWindow() = default;
    PhasingWindow(uint intended_windows_size, uint overlap_length);      //window_size without overlap
    PhasingWindow(const PhasingWindow &rhs);
    ~PhasingWindow() = default;
    void clean();
    void clear();


    //inline function
    inline void insert_phased_block_starting_idx(uint idx) {block_idxes[idx] = idx;}
    inline void destroy_merged_block(uint idx)                  //call the function to destroy merged block after mergign or joining PhasedBlock
    {
        //if (blocks_not_to_be_deleted.count(idx) == 0)
        //    delete blocks[idx];
        blocks.erase(idx);
        block_idxes.erase(idx);
    }
    void update_phasing_info(int max_barcode_spanning_length);
    void update_phasing_info();
    void update_phasing_info_keep_phased();
    void initialize_for_recursive_solver(uint intended_window_size, uint overlap);

    inline bool in_range(uint var_idx) { return variant2mat_index.count(var_idx) != 0;}
    inline uint mat_idx2var_idx(uint mat_idx) { return mat2variant_index[mat_idx];}          //bug here
    inline uint var_idx2mat_idx(uint var_idx) { return variant2mat_index[var_idx];}

    inline uint curr_window_block_count() { return  uint(current_window_idxes.size());}

    inline void insert_block_initial(uint start_var_idx, uint var_idx, ptr_ResultforSingleVariant result)
    {
        if (blocks.count(start_var_idx) == 0)
        {
            blocks[start_var_idx] = std::make_shared<PhasedBlock> (start_var_idx, result);
            result->block = blocks[start_var_idx];
            prev_block_idxes.push_back(start_var_idx);
            block_idxes[start_var_idx] = start_var_idx;
        }
        else {
            auto block = blocks[start_var_idx];
            block->insert_variant(var_idx, result);
            result->block = blocks[start_var_idx];
            //block_idxes[var_idx] = start_var_idx;
        }    
    }
    
    inline void insert_block_initial_recursive_solver(uint start_var_idx, std::shared_ptr<PhasedBlock> &block, bool overlap)
    {
        blocks[start_var_idx] = block;
        prev_block_idxes.push_back(start_var_idx);
        for (auto i : block->variant_idxes)
            variant2mat_index[i] = prev_block_idxes.size() - 1;

    }


    //00: 0; 01: 1
    inline int check_connection_block_recursive_solver(bool is_ref_a, bool is_ref_b, uint variant_a, uint variant_b)
    {
        std::shared_ptr<PhasedBlock> block_a = blocks[mat_idx2var_idx(var_idx2mat_idx(variant_a))];
        std::shared_ptr<PhasedBlock> block_b = blocks[mat_idx2var_idx(var_idx2mat_idx(variant_b))];
        bool is_ref_a_in_block = block_a->results[variant_a]->is_REF();
        bool is_ref_b_in_block = block_b->results[variant_b]->is_REF();

        return (is_ref_a == is_ref_b) != (is_ref_a_in_block == is_ref_b_in_block);
    }

    //operator
    PhasingWindow& operator=(const PhasingWindow &phasing_window);       //copy
    void split_after_phasing(uint idx, std::vector<std::set<uint>> &set_var_idxes);

};

typedef std::shared_ptr<PhasingWindow> ptr_PhasingWindow;

class ChromoPhaser
{
public:
    //initialized parameter
    uint chr_id;
    std::string chr_name;
    uint variant_count;         //variant no to phase
    uint init_block_count;
    uint window_overlap;
    uint intended_window_size, pre_intended_window_size;    //intended_window_size: recusivly increasing
    uint window_size;          //window size with overlap
    uint no_of_recursion;
    uint max_no_of_recursion;

    //to be calculated afterwards

    //ChromoPhaser is obliged to manage the following two container
    std::vector<std::shared_ptr<ResultforSingleVariant> > results_for_variant;
    std::unordered_map<uint, uint> variant_to_block_id;
    

    ptr_PhasingWindow phased;


public:
    ChromoPhaser(const uint &chr_id, const std::string &chr_name, const uint &window_overlap, const uint &intended_window_size);
    ~ChromoPhaser();
    void construct_phasing_window_initialize();                             //initialize struct
    void construct_phasing_window_r_initialize();
    inline void cal_window_recursive_solver()
    {
        // max 3000
        intended_window_size = pre_intended_window_size * no_of_recursion;
        window_overlap = intended_window_size / 20;
        window_size = intended_window_size + window_overlap;
        //window_no = uint(ceil(phased->curr_window_block_count() / double (window_size) ));
    }
    inline uint get_block_count() { return uint(phased->blocks.size());}
    inline uint get_var_pos(uint idx) { return results_for_variant[idx]->get_pos(); }


};

#endif

