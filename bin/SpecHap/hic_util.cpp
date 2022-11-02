//
// Created by yyh on 6/25/2019.
//

#include "hic_util.h"

HiCLinker::HiCLinker(Fragment &fragment)
{
    start_var_idx = fragment.snps[0].first;
    end_var_idx = fragment.snps.back().first;
    for (auto i : fragment.snps)
    {
        const uint idx = i.first;
        var_idx.insert(idx);
        hic_info.insert(i);
    }
    read_qual = fragment.read_qual;
}