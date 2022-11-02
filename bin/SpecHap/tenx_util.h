//
// Created by yyh on 6/25/2019.
//

#ifndef SPECHAP_TENX_UTIL_H
#define SPECHAP_TENX_UTIL_H

#include "type.h"
#include <set>
#include <map>
#include <vector>
#include "frag_io.h"

//TODO: change barcode_info into unordered_map, for multiplicity of snv, sum up their log quality
class Barcode
{
public:
    uint start_var_idx;
    uint end_var_idx;
    std::set<uint> var_idx;
    std::set<uint> bad_snp;
    std::multimap<uint, std::pair<int, double>> barcode_info;                         //var_idx and haplotype

public:
    Barcode();
    explicit Barcode(Fragment &fragment);
    void insert_info(Fragment &fragment);
    void get_barcode_support(ptr_PhasedBlock phased_block);
};


//TODO: change BarcodeLinker into class
typedef std::unordered_map<std::string, Barcode> BarcodeLinker;

class BarcodeLinkers
{
public:
    std::unordered_map<std::string, Barcode> barcodes;
    std::unordered_map<std::string, int> barcode_index;
    ChromoPhaser *chromo_phaser;
    RegionFragStats *regionFragStats;
    uint max_barcode_spanning;

public:
    explicit BarcodeLinkers(uint max_barcode_spanning);
    BarcodeLinkers(ChromoPhaser *chromo_phaser, uint max_barcode_spanning);
    void insert_barcode(Fragment &fragment);

    inline void clear()
    {
        barcodes.clear();
        barcode_index.clear();
    }
};

#endif //SPECHAP_TENX_UTIL_H
