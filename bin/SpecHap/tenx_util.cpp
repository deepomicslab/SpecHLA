//
// Created by yyh on 6/25/2019.
//
#include "tenx_util.h"
#include <climits>

//TODO: rearranging the method for looping through barcode
Barcode::Barcode()
{
    start_var_idx = UINT_MAX;
    end_var_idx = 0;
}

Barcode::Barcode(Fragment &fragment)
{
    start_var_idx = UINT_MAX;
    end_var_idx = 0;

    for (auto i : fragment.snps)
    {
        const uint &idx = i.first;
        if (idx > end_var_idx)
            end_var_idx = idx;
        if (idx < start_var_idx)
            start_var_idx = idx;
        if (bad_snp.find(idx) == bad_snp.end())
        {
            var_idx.insert(idx);
            barcode_info.insert(i);
        }
    }
}

void Barcode::insert_info(Fragment &fragment)
{
    if (fragment.read_qual > -3)
        return;
    for (auto i : fragment.snps)
    {
        const uint &idx = i.first;
        if (bad_snp.find(idx) != bad_snp.end())
            continue;
        if (idx > end_var_idx)
            end_var_idx = idx;
        if (idx < start_var_idx)
            start_var_idx = idx;

        std::pair<std::multimap<uint, std::pair<int, double>>::iterator, std::multimap<uint, std::pair<int,double>>::iterator> iter = barcode_info.equal_range(idx);
        //not found
        if (iter.first == iter.second)
        {
            var_idx.insert(idx);
            barcode_info.insert(i);
        }
        else
        {
            for (auto it = iter.first; it != iter.second; ++it)
            {
                //bad snp
                if (i.second.first != it->second.first)
                {
                    bad_snp.insert(idx);
                    var_idx.erase(var_idx.find(idx));
                    barcode_info.erase(iter.first, iter.second);
                    break;
                }
            }
        }
    }
}

void Barcode::get_barcode_support(ptr_PhasedBlock phased_block)
{
    // notice that phased snvs belong to only one block
    // barcode support at most one block
    if (this->end_var_idx < phased_block->start_variant_idx || this->start_var_idx > *phased_block->variant_idxes.rbegin())
        return;
    std::set<uint> supported_idx; //the snv in the phased block that this barcode supported
    double count_same_hap = 0;
    double count_same_hap_c = 0;
    for (auto variant_idx: var_idx)
    {
        if (phased_block->contain_variant(variant_idx))
        {
            supported_idx.insert(variant_idx);
            auto hap_qual = this->barcode_info.find(variant_idx);
            double qual = (int(hap_qual->second.second) - 33) / 10.0;
            if (hap_qual->second.first == phased_block->results[variant_idx]->get_hap())
            {
                count_same_hap += qual;
                count_same_hap_c -= qual;
            }
            else
            {
                count_same_hap -= qual;
                count_same_hap_c += qual;
            }
        }
    }
    if (supported_idx.empty() || supported_idx.size() == 1)
        return;
    count_same_hap /= supported_idx.size();
    count_same_hap_c /= supported_idx.size();
    double h0p, h1p;
    h0p = double(count_same_hap) / supported_idx.size();
    h1p = 1 - h0p;

    if (count_same_hap == count_same_hap_c)
        return;

    if (count_same_hap > count_same_hap_c && count_same_hap <= 3)
        return;
    if (count_same_hap_c > count_same_hap && count_same_hap_c <= 3)
        return;

    int l_hap;
    for (auto i: supported_idx)
    {
        l_hap = this->barcode_info.find(i)->second.first;
        ptr_ResultforSingleVariant result = phased_block->results[i];
        //global haplotype 0
        if (count_same_hap > count_same_hap_c)
        {
            //local haplotype 0
            if (l_hap == result->get_hap())
                result->L0_G0++;
            else
                result->L1_G0++;
        }
            //global haplotype 1
        else
        {
            //local haplotype 0
            if (l_hap == result->get_hap())
                result->L0_G1++;
            else
                result->L1_G1++;
        }
    }
}

BarcodeLinkers::BarcodeLinkers(ChromoPhaser *chromo_phaser, uint max_barcode_spanning):
chromo_phaser(chromo_phaser), max_barcode_spanning(max_barcode_spanning), regionFragStats(nullptr)
{}

void BarcodeLinkers::insert_barcode(Fragment &fragment)
{
    uint fragment_start = fragment.snps[0].first;
    //locate the range for fragment
    int idx = this->regionFragStats->get_index_for_barcode_variant(fragment.barcode, chromo_phaser->get_var_pos(fragment_start));
    //range not found, this cannot be used
    if (idx <= 0)
        return;

    //first time see this barcode in window

    std::string barcode_key = fragment.barcode + std::string("_") + std::to_string(idx);

    //not met
    if (barcodes.count(barcode_key) == 0)
        barcodes[barcode_key] = Barcode(fragment);

    else
        barcodes[barcode_key].insert_info(fragment);

}

BarcodeLinkers::BarcodeLinkers(uint max_barcode_spanning)
:max_barcode_spanning(max_barcode_spanning), chromo_phaser(nullptr), regionFragStats(nullptr)
{}
