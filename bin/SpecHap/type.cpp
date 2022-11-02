//
// Created by yyh on 12/10/2018.
//
#include "type.h"
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>

Fragment::Fragment(const Fragment &rhs)
{
    this->read_qual = rhs.read_qual;
    this->snps = rhs.snps;
    this->insertion_size = rhs.insertion_size;
    this->rescued = rhs.rescued;
    this->dm = rhs.dm;
    this->start = rhs.start;
    this->end = rhs.end;
}

int Fragment::getType() const {
    return type;
}

void Fragment::setType(int Ftype) {
    Fragment::type = Ftype;
}

//------------------------------- PhasedBlock Class ------------------------------------//

PhasedBlock::PhasedBlock() { std::cout << "You shouldn't be here"; }
PhasedBlock::PhasedBlock(uint start_idx)
{
    this->start_variant_idx = start_idx;
}
PhasedBlock::PhasedBlock(uint variant_start_idx, ptr_ResultforSingleVariant result)
{
    this->start_variant_idx = variant_start_idx;
    results[variant_start_idx] = result;
    variant_idxes.insert(start_variant_idx);
}

//do not use this constructor
PhasedBlock::PhasedBlock(const PhasedBlock &rhs)
{
    start_variant_idx = rhs.start_variant_idx;
    results = rhs.results;                  //bug here, should make duplicate
    variant_idxes = rhs.variant_idxes;
}

PhasedBlock::PhasedBlock(const PhasedBlock &rhs, std::unordered_map<uint, ptr_ResultforSingleVariant> &results_dup)
{
    start_variant_idx = rhs.start_variant_idx;
    variant_idxes = rhs.variant_idxes;
    for (uint var_idx : variant_idxes)
        results[var_idx] = results_dup[var_idx];
}

PhasedBlock::~PhasedBlock()
{
    results.clear();
    variant_idxes.clear();
}

//PhasedBlock::PhasedBlock(const PhasedBlock &rhs) {}

void PhasedBlock::flip()
{
    for (auto it : results)
        it.second->complement_hap();
}

void PhasedBlock::join_block_no_overlap(std::shared_ptr<PhasedBlock> rhs)
{
    for (auto it: rhs->variant_idxes)
    {
        variant_idxes.insert(it);
        results[it] = rhs->results[it];
        results[it]->block = shared_from_this();
    }
    //delete rhs with phasing window
}


//Note not to use the results from overlap region
void PhasedBlock::join_block_overlap(std::shared_ptr<PhasedBlock> rhs)
{
    std::set<uint> intersection;    //variant idx of intersection variant, merge and delete the dup
    std::set<uint> diff;            //variant idx in rhs while not in this block

    std::set_intersection(variant_idxes.begin(), variant_idxes.end(), rhs->variant_idxes.begin(), rhs->variant_idxes.end(), std::inserter(intersection, intersection.begin()));
    std::set_difference(rhs->variant_idxes.begin(), rhs->variant_idxes.end(), variant_idxes.begin(), variant_idxes.end(), std::inserter(diff, diff.begin()));

    uint intersection_length, identical_length;
    intersection_length = uint(intersection.size());
    identical_length = 0;

    if (rhs->size() > 1)
    {
        results[start_variant_idx]->set_phased();
        results[start_variant_idx]->block = shared_from_this();
    }

    //use left hand result
    for (auto it : intersection)
    {
        ptr_ResultforSingleVariant a = this->results[it];
        ptr_ResultforSingleVariant b = rhs->results[it];
        if (a->variant_phased() && b->variant_phased())
        {
            if ((a->is_REF() && b->is_REF()) || (a->is_ALT() && b->is_ALT()))
                identical_length++;
        }
    }
    if (identical_length < intersection_length / 2)
        rhs->flip();


    for (auto it: diff)
    {
        variant_idxes.insert(it);
        results[it] = rhs->results[it];
        results[it]->block = shared_from_this();
    }
    //delete rhs out of this block in the phasing window
}

void PhasedBlock::filter_inconsistency(std::set<uint> &filtered_idx)
{
    for (auto &i : results)
    {
        if (i.second->filter_inconsistency())
            filtered_idx.insert(i.first);
    }
}

void PhasedBlock::fragment_supported_score(Fragment &fragment)
{
    for (auto & i : fragment.snps)
    {
        if (!contain_variant(i.first))
            continue;
    }
}

//------------------------------- PhasedBlock Class End ------------------------------------//


//------------------------------- PhasingWindow Class ------------------------------------//
// intended_window_size: window_size without overlap
PhasingWindow::PhasingWindow(uint intended_windows_size,  uint overlap_length)
{
    intended_window_length = intended_windows_size;
    phased_blk_count = phasing_block_size = 0;
    this->end = 0;

    this->overlap_length = overlap_length;
    this->intended_start = 0;
}

PhasingWindow::PhasingWindow(const PhasingWindow &rhs)
{
    std::cerr << "error copying phasing window" << std::endl;
    exit(255);
}


// create single variant phased block

void PhasingWindow::initialize_for_recursive_solver(uint intended_window_size, uint overlap)
{
    this->intended_window_length = intended_window_size;
    this->overlap_length = overlap;
    intended_start = phased_blk_count = phasing_block_size = 0;
    total_blk_count = rest_blk_count = uint(blocks.size());
    prev_block_idxes.clear();
    current_window_idxes.clear();
    for (auto &i : block_idxes)
        prev_block_idxes.push_back(i.first);
}

void PhasingWindow::update_phasing_info(int max_barcode_spanning_length)
{
    //three parameters to be determined:
    //starting point, ending point, and block count
    auto i = block_idxes.find(prev_block_idxes[intended_start]);
    uint count = 0;
    //the first window
    phasing_block_size = 0;
    variant2mat_index.clear();
    mat2variant_index.clear();
    current_window_idxes.clear();
    current_window_overlap_region_block_size.clear();
    uint temp = intended_start;

    // first window
    if (i == block_idxes.begin())
    {
        phased_blk_count = 0;
        start = intended_start;
        if (phased_blk_count + intended_window_length >= prev_block_idxes.size())
            end = intended_end = prev_block_idxes.back() + 1;
        else
        {
            intended_end = prev_block_idxes[phased_blk_count + intended_window_length];
            if (phased_blk_count + intended_window_length + overlap_length >= prev_block_idxes.size())
                end = prev_block_idxes.back() + 1;
            else
                end = prev_block_idxes[phased_blk_count + intended_window_length + overlap_length];
        }
        intended_start = intended_start + intended_window_length + overlap_length;
    }
    else
    {
        //determine true start
        while (i != block_idxes.begin() && count < overlap_length  )
        {
            i = prev(i);
            if (i->first < intended_start)
            {
                ptr_PhasedBlock &blk = blocks[i->first];
                if (blk->size() == 1)
                    continue;
                else
                    break;
                //if (blk->results[*blk->variant_idxes.rbegin()]->get_pos() + max_barcode_spanning_length < start_blk->get_pos())
                //    break;
            }
            count++;
        }

        start = i->first;
        if (phased_blk_count + intended_window_length - overlap_length >= prev_block_idxes.size())
            end = intended_end = prev_block_idxes.back() + 1;
        else
        {
            intended_end = prev_block_idxes[phased_blk_count + intended_window_length - overlap_length];
            if (phased_blk_count + intended_window_length >= prev_block_idxes.size())
                end = prev_block_idxes.back() + 1;
            else
                end = prev_block_idxes[phased_blk_count + intended_window_length];
        }
        intended_start = intended_start + intended_window_length;
    }
    auto it = i;
    uint threhold;
    long(temp) - int(overlap_length) <= 0 ? threhold = 0 : threhold = temp - overlap_length;

    for (; it != block_idxes.end() && it->first != end; it++)
    {
        int size_count = 0;
        ptr_PhasedBlock blk = blocks[it->first];
        if ((*blk->variant_idxes.rbegin()) < threhold)
            continue;
        if (blk->size() == 1)
        {
            if (blk->results.begin()->second->get_filter() != filter_type::PASS)
                continue;
            variant2mat_index[blk->start_variant_idx] = phasing_block_size;
        }
        else
        {

            std::set<uint>::iterator var_itr = blk->variant_idxes.begin();
            for (uint j = threhold; j < temp; j++)
            {
                uint idx = prev_block_idxes[j];
                if (blk->results.count(idx) == 0)
                    continue;
                else
                {
                    var_itr = blk->variant_idxes.find(idx);
                    break;
                }
            }
            size_count = *blk->variant_idxes.rbegin() - threhold;
            variant2mat_index[blk->start_variant_idx] = phasing_block_size;
            for (; var_itr != blk->variant_idxes.end(); ++var_itr)
            {
                uint idx = *var_itr;
                variant2mat_index[idx] = phasing_block_size;
                if (idx < temp)
                    current_window_overlap_region_block_size[idx] = size_count;
            }
        }
        mat2variant_index[phasing_block_size] = it->first;
        current_window_idxes.push_back(it->first);
        phasing_block_size++;
    }


    phased_blk_count += (end - temp);
    prev_window_start = temp;
    rest_blk_count = total_blk_count - phased_blk_count;
}

void PhasingWindow::update_phasing_info_keep_phased()
{
    auto i = block_idxes.find(end);
    uint count = 0;
    //the first window
    phasing_block_size = 0;
    variant2mat_index.clear();
    mat2variant_index.clear();
    current_window_idxes.clear();
    uint temp = intended_start;
    uint temp_start = 0;
    current_window_overlap_region_block_size.clear();

    // if it is the first interval
    if (i == block_idxes.begin())
    {
        phased_blk_count = 0;
        start = intended_start;
        if (phased_blk_count + intended_window_length >= prev_block_idxes.size())
            end = intended_end = prev_block_idxes.back() + 1;
        else
        {
            intended_end = prev_block_idxes[phased_blk_count + intended_window_length];
            if (phased_blk_count + intended_window_length + overlap_length >= prev_block_idxes.size())
                end = prev_block_idxes.back() + 1;
            else
                end = prev_block_idxes[phased_blk_count + intended_window_length + overlap_length];
        }
        intended_start = intended_start + intended_window_length + overlap_length;
    }
    else {
        
        //determine true start
        temp_start = i ->first;
        while ( count < overlap_length  )
        {
            count++;
            
            if (i == block_idxes.begin())
                break;
            i = prev(i);
        }
        start = i->first;
        
        if (phased_blk_count + intended_window_length - overlap_length >= prev_block_idxes.size())
            end = intended_end = prev_block_idxes.back() + 1;
        else
        {
            intended_end = prev_block_idxes[phased_blk_count + intended_window_length - overlap_length];
            if (phased_blk_count + intended_window_length >= prev_block_idxes.size())
                end = prev_block_idxes.back() + 1;
            else
                end = prev_block_idxes[phased_blk_count + intended_window_length];
        }
        intended_start = intended_start + intended_window_length;
    }
    auto it = i;

    int offset = 0;
    for (; it != block_idxes.end() && it->first != end; it++, phasing_block_size++)
    {
        
        ptr_PhasedBlock blk = blocks[it->first];
        if (blk->size() == 1)
            if (blk->results.begin()->second->get_filter() == filter_type::POOLRESULT)
                continue;
        for (auto idx : blk->variant_idxes)
            variant2mat_index[idx] = phasing_block_size;
        mat2variant_index[phasing_block_size] = it->first;
        current_window_idxes.push_back(it->first);
        if (it->first >= temp_start)
        {
            offset ++; 
            //if (blk->size() > 1)
                //offset = offset - blk->size() + 1; 
        }
    }


    phased_blk_count += offset;
    prev_window_start = temp;
    rest_blk_count = total_blk_count - phased_blk_count;
}

void PhasingWindow::update_phasing_info()
{
    //three parameters to be determined:
    //starting point, ending point, and block count

    //TODO found a bug here, prev_block_idxes might not contain intended start if block read from vcf 
    // The index system works for single snv node instead of block
    auto i = block_idxes.find(prev_block_idxes[intended_start]);
    uint count = 0;
    //the first window
    phasing_block_size = 0;
    variant2mat_index.clear();
    mat2variant_index.clear();
    current_window_idxes.clear();
    uint temp = intended_start;
    current_window_overlap_region_block_size.clear();


    if (i == block_idxes.begin())
    {
        phased_blk_count = 0;
        start = intended_start;
        if (phased_blk_count + intended_window_length >= prev_block_idxes.size())
            end = intended_end = prev_block_idxes.back() + 1;
        else
        {
            intended_end = prev_block_idxes[phased_blk_count + intended_window_length];
            if (phased_blk_count + intended_window_length + overlap_length >= prev_block_idxes.size())
                end = prev_block_idxes.back() + 1;
            else
                end = prev_block_idxes[phased_blk_count + intended_window_length + overlap_length];
        }
        intended_start = intended_start + intended_window_length + overlap_length;
    }
    else {
        //determine true start
        while ( count < overlap_length  )
        {
            count++;

            if (i == block_idxes.begin())
                break;
            i = prev(i);
        }
        start = i->first;
        if (phased_blk_count + intended_window_length - overlap_length >= prev_block_idxes.size())
            end = intended_end = prev_block_idxes.back() + 1;
        else
        {
            intended_end = prev_block_idxes[phased_blk_count + intended_window_length - overlap_length];
            if (phased_blk_count + intended_window_length >= prev_block_idxes.size())
                end = prev_block_idxes.back() + 1;
            else
                end = prev_block_idxes[phased_blk_count + intended_window_length];
        }
        intended_start = intended_start + intended_window_length;
    }
    auto it = i;

    for (; it != block_idxes.end() && it->first != end; it++, phasing_block_size++)
    {
        ptr_PhasedBlock blk = blocks[it->first];
        if (blk->size() == 1)
            if (blk->results.begin()->second->get_filter() == filter_type::POOLRESULT)
                continue;
        for (auto idx : blk->variant_idxes)
            variant2mat_index[idx] = phasing_block_size;
        mat2variant_index[phasing_block_size] = it->first;
        current_window_idxes.push_back(it->first);
    }


    phased_blk_count += (end - temp);
    prev_window_start = temp;
    rest_blk_count = total_blk_count - phased_blk_count;
}

void PhasingWindow::clean()
{
    blocks.clear();
    variant2mat_index.clear();
    prev_block_idxes.clear();
    block_idxes.clear();
}

void PhasingWindow::clear()
{
    variant2mat_index.clear();
    mat2variant_index.clear();
    current_window_idxes.clear();
    current_window_overlap_region_block_size.clear();
}

PhasingWindow& PhasingWindow::operator=(const PhasingWindow &phasing_window)
{
    overlap_length = phasing_window.overlap_length;

    blocks = phasing_window.blocks;
    prev_block_idxes = phasing_window.prev_block_idxes;
    block_idxes = phasing_window.block_idxes;
    return *this;
}

void PhasingWindow::split_after_phasing(uint idx, std::vector<std::set<uint>> &set_var_idxes)
{
    ptr_PhasedBlock block_tobe_split = blocks[idx];
    auto unary_op = [=] (uint i) {return mat2variant_index[i];};
    for (auto &var_idxes_mat : set_var_idxes)
    {
        std::set<uint> var_idxes;
        std::transform(var_idxes_mat.begin(), var_idxes_mat.end(), std::inserter(var_idxes,var_idxes.begin()), unary_op);
        uint start_idx = *var_idxes.begin();

        //disjoint, remove from phased blk
        if (var_idxes.size() == 1)
        {
            block_tobe_split->results[start_idx]->set_unphased();
            block_tobe_split->results[start_idx]->block.reset();
        }

        //be kept in origin blk
        if (start_idx == block_tobe_split->start_variant_idx)
        {
            block_tobe_split->variant_idxes = var_idxes; //copy
            continue;
        }
        ptr_PhasedBlock phased_block = std::make_shared<PhasedBlock>(start_idx);
        phased_block->variant_idxes = var_idxes;
        for (auto var_idx : var_idxes)
        {
            //TODO: make sure block_tobe_split->result[var_idx] always exists
            ptr_ResultforSingleVariant result = block_tobe_split->results[var_idx];
            phased_block->results[var_idx] = result;
            result->update_phased_block(phased_block);
            block_tobe_split->results.erase(var_idx);
            block_tobe_split->variant_idxes.erase(var_idx);
        }

        blocks[start_idx] = phased_block;
        block_idxes[start_idx] = start_idx;
    }
}

//------------------------------- PhasingWindow Class End ------------------------------------//


//------------------------------- ChromoPhaser Class ------------------------------------//

ChromoPhaser::ChromoPhaser(const uint &chr_id, const std::string &chr_name, const uint &window_overlap, const uint &intended_window_size)
         : chr_id(chr_id), chr_name(chr_name), window_overlap(window_overlap), intended_window_size(intended_window_size)
{
    pre_intended_window_size = intended_window_size;
    phased = std::make_shared<PhasingWindow>(intended_window_size,  window_overlap);
    no_of_recursion = 1;
    max_no_of_recursion = 3;            //for control of window size
    window_size = window_overlap + intended_window_size;

}

ChromoPhaser::~ChromoPhaser()
{
    //for (auto i : results_for_variant)
    //    delete i;
    //for (auto i : overlap_region_result_dup)
    //    delete i.second;
    //for (auto i : overlap_region_block_dup)
    //    delete i.second;

    results_for_variant.clear();
}

void ChromoPhaser::construct_phasing_window_initialize()
{
    for (uint i = 0; i < variant_count; i++)
        phased->insert_block_initial(this->variant_to_block_id[i], i, results_for_variant[i]);
    phased->total_blk_count = phased->rest_blk_count = phased->blocks.size();
}

void ChromoPhaser::construct_phasing_window_r_initialize()
{
    if (no_of_recursion < max_no_of_recursion)
        no_of_recursion++;
    cal_window_recursive_solver();
    phased->initialize_for_recursive_solver(intended_window_size, window_overlap);
}

/*
void ChromoPhaser::update_overlap_region_block_idx()
{
    for (uint i = 0; i < work.size(); i++)
    {
        PhasingWindow &phasing_window = work[i];
        PhasedBlock *prev_block_ptr = nullptr;
        for (uint j = phasing_window.window_end_idx; j < phasing_window.window_end_overlap; j++)
        {
            ResultforSingleVariant *resultforSingleVariant = results_for_variant[j];

            if (resultforSingleVariant->block == nullptr) //disjointed
                phasing_window.overlap_region_block_idxes.push_back(j);
            else
            {
                if (prev_block_ptr == resultforSingleVariant->block)
                    continue;
                else
                    phasing_window.overlap_region_block_idxes.push_back(resultforSingleVariant->block->start_variant_idx);
            }
            prev_block_ptr = resultforSingleVariant->block;
        }
    }
}
*/


//------------------------------- ChromoPhaser Class End ------------------------------------//