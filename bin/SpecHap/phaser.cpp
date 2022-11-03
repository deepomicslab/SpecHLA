//
// Created by yyh on 3/4/2019.
//

#include "phaser.h"
#include "htslib/vcf.h"
#include "type.h"
#include "util.h"
#include <iostream>
#include <cstdlib>
#include <unordered_map>

// TODO clarify between variant count and block count

Phaser::Phaser(const std::string &fnvcf, const std::string &fnout, std::vector<std::string> & fnfrags, const std::string &fnbed, std::vector<double> & fr_weights)
{
    frbed = nullptr;
    for (int i = 0; i < fnfrags.size(); i++) {
        auto item = fnfrags[i];
        frfrags.push_back(new FragmentReader(item.data()));
        if(OPERATIONS[i] == MODE_10X) {
            frbed = new BEDReader(fnbed.c_str());
        }
    }
    frvcf = new VCFReader(fnvcf.data());
    fwvcf = new VCFWriter(frvcf->header, fnout.data());
//    frfrag = new FragmentReader(fnfrag.data());
//    frbed = nullptr;
    coverage = 30;  //deprecated
//
//    if (HAS_TENX)
//        frbed = new BEDReader(fnbed.data());

    bool use_secondary = false;
    threshold = 1e-5;
//    threshold = 0;

    spectral = new Spectral(frfrags, fr_weights, frbed, threshold, coverage, use_secondary);
}

Phaser::~Phaser()
{
    delete frvcf;
    delete fwvcf;
//    delete frfrag;
    for (auto & item : frfrags) {
        delete item;
        item = nullptr;
    }
    frfrags.clear();
    frfrags.shrink_to_fit();
    delete spectral;
    if (frbed != nullptr)
        delete frbed;
}


int Phaser::load_contig_records(ChromoPhaser *chromo_phaser)
{
    int status = 0;
    while (true)
    {
        ptr_ResultforSingleVariant result = std::make_shared<ResultforSingleVariant>();
        int status = this->frvcf->get_next_record_contig(result, false);
        if (status < 0)
            break;
//        else if (status > 0)
//            continue;
        chromo_phaser->results_for_variant.push_back(result);
    }


    chromo_phaser->variant_count = chromo_phaser->results_for_variant.size();
    for (int i = 0; i < chromo_phaser->variant_count; i++)
    {
        chromo_phaser->variant_to_block_id[i] = i;
    }
    chromo_phaser->init_block_count = chromo_phaser->variant_count;
    return status;
}


int Phaser::load_contig_blocks(ChromoPhaser *chromo_phaser)
{
    int status = 0;
    while (true)
    {
        ptr_ResultforSingleVariant result = std::make_shared<ResultforSingleVariant>();
        int status = this->frvcf->get_next_record_contig(result, true);
        if (status < 0) //eof, new chromosome 
            break;
//        else if (status > 0) //homo
//            continue;
        chromo_phaser->results_for_variant.push_back(result);
    }

    chromo_phaser->variant_count = chromo_phaser->results_for_variant.size();
    
    
    std::unordered_map<uint, uint> ps2block_ids;
    uint block_count = 0;
    for (int i = 0; i < chromo_phaser->variant_count; i++)
    {
        auto result = chromo_phaser->results_for_variant[i];
        uint ps = result->ps;
        if (ps == 0) //not phased 
        {
            chromo_phaser->variant_to_block_id[i] = i;
            block_count++;
        }
        else {      //phased
            if (ps2block_ids.count(ps) == 0)
            {   // not met before
                ps2block_ids[ps] = i;
                chromo_phaser->variant_to_block_id[i] = ps2block_ids[ps];
                block_count++;
            }
            else {
                chromo_phaser->variant_to_block_id[i] = ps2block_ids[ps];
            }
        }
    }
    chromo_phaser->init_block_count = block_count;
    return status;
}


void Phaser::phasing()
{
    uint prev_variant_count = 0;

    for (uint rid = 0; rid < frvcf->contigs_count; rid++)
    {
        if (frvcf->jump_to_contig(rid) != 0)
            break;

        ChromoPhaser *chromo_phaser = new ChromoPhaser(rid, frvcf->contigs[rid], WINDOW_OVERLAP, WINDOW_SIZE);
        std::string mess = "phasing haplotype for " + std::string(frvcf->contigs[rid]);
        logging(std::cerr, mess);
        if (KEEP_PS)
            load_contig_blocks(chromo_phaser);
        else
            load_contig_records(chromo_phaser);
        chromo_phaser->construct_phasing_window_initialize();
        for (auto item : frfrags) {
            item->set_prev_chr_var_count(prev_variant_count);
        }
//        frfrag->set_prev_chr_var_count(prev_variant_count);
        spectral->set_chromo_phaser(chromo_phaser);
        if (!contigs.empty() && (std::find(contigs.begin(), contigs.end(), frvcf->contigs[rid]) == contigs.end())) {
            std::string mess = " skipp " + std::string(frvcf->contigs[rid]);
            logging(std::cerr, mess);
//            continue;
            fwvcf->write_nxt_contigs(frvcf->contigs[rid].data(), chromo_phaser, *frvcf,spectral->getBreakIdxs());
//            prev_variant_count += chromo_phaser->variant_count;
//            spectral->release_chromo_phaser();
//            delete chromo_phaser;
        } else {
            phasing_by_chrom(chromo_phaser->variant_count, chromo_phaser);
            fwvcf->write_nxt_contigs(frvcf->contigs[rid].data(), chromo_phaser, *frvcf, spectral->getBreakIdxs());
        }

        //write vcf
        prev_variant_count += chromo_phaser->variant_count;
        spectral->release_chromo_phaser();
        delete chromo_phaser;
    }
}


void Phaser::phasing_by_chrom(uint var_count, ChromoPhaser *chromo_phaser)
{
    for (auto item : frfrags) {
        item->set_curr_chr_var_count(var_count);
    }
//    frfrag->set_curr_chr_var_count(var_count);
    this->spectral->setOffset(0);
    while (chromo_phaser->phased->rest_blk_count > 0)
    {
        if (chromo_phaser->phased->rest_blk_count > chromo_phaser->init_block_count)
            break;
        if (HAS_TENX)
            chromo_phaser->phased->update_phasing_info(MAX_BARCODE_SPANNING);
        else
            {
                if (KEEP_PS)
                    chromo_phaser->phased->update_phasing_info_keep_phased();
                else 
                    chromo_phaser->phased->update_phasing_info();
            }
        spectral->solver();
    }
        //std::cout << chromo_phaser->phased->phased_blk_count << std::endl;
    if (HAS_HIC)
    {
        phase_HiC_poss(chromo_phaser);
    }
}


void Phaser::phase_HiC_poss(ChromoPhaser *chromo_phaser)
{
    std::unordered_map<uint, std::set<uint>> connected_comps = spectral->load_hic_poss_info();
    for (auto i : connected_comps)
    {
        std::set<uint> &connected_comp = i.second;
        int nblocks = connected_comp.size();
        if (nblocks == 1)
            continue;

        int count = 0;
    
        //split into window again
        int HiC_poss_block = WINDOW_SIZE, overlap = WINDOW_OVERLAP;
        if (nblocks > HiC_poss_block + overlap)
        {
            this->phase_HiC_recursive(chromo_phaser, connected_comp);  
        }
        //direct phase
        else 
        {
            chromo_phaser->phased->clear();
            count = 0;
            //update indexing scheme
            for (auto blk_start_id: connected_comp)
            {
                chromo_phaser->phased->current_window_idxes.push_back(blk_start_id);
                chromo_phaser->phased->mat2variant_index[count] = blk_start_id;
                ptr_ResultforSingleVariant variant =  chromo_phaser->results_for_variant[blk_start_id];
                if (is_uninitialized(variant->block))
                {
                    chromo_phaser->phased->variant2mat_index[blk_start_id] = count;
                }
                else
                {
                    auto blk = variant->block.lock();
                    for (auto _var_id : blk->variant_idxes)
                    {
                        chromo_phaser->phased->variant2mat_index[_var_id] = count;
                    }   
                }
                count ++;
            }

            spectral->hic_poss_solver(nblocks);        
        }
    }
    
}

void Phaser::phase_HiC_recursive(ChromoPhaser *chromo_phaser, std::set<uint> &connected_comp)
{
    int count = 0;
    std::unordered_map<uint, uint> var2id;
    std::unordered_map<uint, uint> id2var;
    std::vector<uint> prev_block_idxes;
    std::map<uint, uint> block_idxes;
    int HiC_poss_block = WINDOW_SIZE, overlap = WINDOW_OVERLAP;
    int n_recursion = 0;

    for (auto blk_start_id: connected_comp)
    {
        var2id[blk_start_id] = count;
        id2var[count] = blk_start_id;    
        block_idxes[count] = count; 
        prev_block_idxes.push_back(count);
        count++; 
    }
    int prev_blk_count = block_idxes.size();
    while (n_recursion < RECURSIVE_LIMIT || block_idxes.size() != prev_blk_count)
    {
        //now for each recurssion

        int nblocks = prev_block_idxes.size();

        int phased = 0, rest = nblocks;
        uint start = 0, end = 0;
        uint intend_start = 0;
        uint tt = 0;
        int current_window_size = 0;

        while (rest > 0)
        {
            chromo_phaser->phased->clear();
            auto i = block_idxes.find(end);
                tt = intend_start;
            uint count = 0;
            //the first window
            current_window_size = 0;

            if (i == block_idxes.begin())
            {
                phased = 0;

                if (phased + HiC_poss_block >= prev_block_idxes.size())
                    end = prev_block_idxes.back() + 1;
                else
                {
                    if (HiC_poss_block + phased + overlap >= prev_block_idxes.size())
                        end = prev_block_idxes.back() + 1;
                    else
                        end = prev_block_idxes[phased + HiC_poss_block + overlap];
                }
                    intend_start = intend_start + HiC_poss_block + overlap;
            }
            else {
            //determine true start

                while ( count < overlap )
                {
                    count++;

                    if (i == block_idxes.begin())
                        break;
                    i = prev(i);
                }
                start = i->first;
                if (phased + HiC_poss_block - overlap >= prev_block_idxes.size())
                    end  = prev_block_idxes.back() + 1;
                else
                {
                    if (phased + HiC_poss_block >= prev_block_idxes.size())
                        end = prev_block_idxes.back() + 1;
                    else
                        end = prev_block_idxes[HiC_poss_block + phased];
                }
                intend_start = intend_start + HiC_poss_block ;
            }

            auto it = i;

            for (; it != block_idxes.end() && it->first != end; it++, current_window_size++)
            {

                uint blk_idx = id2var[it->first];
                ptr_PhasedBlock blk = chromo_phaser->phased->blocks[blk_idx];
                if (blk->size() == 1)
                    if (blk->results.begin()->second->get_filter() == filter_type::POOLRESULT)
                            continue;
                for (auto idx : blk->variant_idxes)
                    chromo_phaser->phased->variant2mat_index[idx] = current_window_size;
                chromo_phaser->phased->mat2variant_index[current_window_size] = blk_idx;
                chromo_phaser->phased->current_window_idxes.push_back(blk_idx);
            }

            phased += (end - tt);
            rest = nblocks - phased;

            //now do poss phasing 
            spectral->hic_poss_solver(chromo_phaser->phased->current_window_idxes.size());

                    //now update the index
            for (auto it : chromo_phaser->phased->current_window_idxes)
            {
                        //its been phased! update index accordingly
                if (chromo_phaser->phased->block_idxes.count(it) == 0)
                {
                    block_idxes.erase(var2id[it]);
                }
            }
        }

        //update index after the recurssion
        std::unordered_map<uint, uint> _var2id;
        std::unordered_map<uint, uint> _id2var;
        std::vector<uint> _prev_block_idxes;
        std::map<uint, uint> _block_idxes;

        count = 0;
        for (auto i : block_idxes)
        {
            _block_idxes[count] = count; 
            _prev_block_idxes.push_back(count);
            uint var = id2var[i.first];
            _var2id[var] = count;
            _id2var[count] = var;
            count ++;
        }
        prev_blk_count = block_idxes.size();
        prev_block_idxes = _prev_block_idxes;
        block_idxes = _block_idxes;
        var2id = _var2id;
        id2var = _id2var;

        n_recursion++;
    }
}

void Phaser::set_contigs(std::string &s) {
    size_t pos = 0;
    std::string token;
    if ((pos = s.find(',')) == std::string::npos){
        contigs.push_back(s);
        return;
    }
    while ((pos = s.find(',')) != std::string::npos) {
        token = s.substr(0, pos);
        contigs.push_back(token);
        s.erase(0, pos + 1);
    }
    contigs.push_back(s);
}

