//
// Created by yyh on 6/25/2019.
//

#ifndef SPECHAP_HIC_UTIL_H
#define SPECHAP_HIC_UTIL_H
#include <set>
#include "type.h"

class HiCLinker
{
public:
    uint start_var_idx;
    uint end_var_idx;
    std::set<uint> var_idx;
    std::set<std::pair<uint, std::pair<int, double>>> hic_info;
    double read_qual;

public:
    explicit HiCLinker(Fragment &fragment);
};

class HiCLinkerContainer
{
public:
    std:: multimap<uint, HiCLinker> linker;

public:
    HiCLinkerContainer() = default;
    inline void add_HiC_info(Fragment &fragment)
    {
        uint start_idx = fragment.snps[0].first;
        linker.emplace(start_idx, fragment);
    }
    inline void clear()
    {
        this->linker.clear();
    }
};

#endif //SPEC_HIC_UTIL_H
