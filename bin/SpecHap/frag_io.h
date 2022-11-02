//
// Created by yonghanyu2 on 8/10/2018.
//

#ifndef PHASER_FRAG_IO_H
#define PHASER_FRAG_IO_H
#include <htslib/vcf.h>
#include <fstream>
#include <cstring>
#include <string>
#include <htslib/tbx.h>
#include "type.h"
#include <cmath>
#include "util.h"

extern int BASE_OFFSET;
extern bool NEW_FORMAT; 
extern bool HYBRID;

template < class ContainerT >
void tokenize(const std::string &str, ContainerT &tokens, const std::string &delimiters=" ", bool trimEmpty=false);

class FragmentReader
{
private:
    uint32_t curr_chr_var_count;
    uint32_t prev_chr_var_count;
    uint32_t intended_window_end;
    uint32_t window_start, window_end;              //note that value corresponding with the index in the fragment file
    std::fstream::streampos nxt_window_start;
    bool nxt_window_set;
    std::fstream frag_file;
    std::vector<std::string> buffer;

private:
    inline std::streampos tell() { return this->frag_file.tellg(); }
    inline void seek(std::streampos pos) { this->frag_file.seekg(pos); }
    inline void get_chr(std::string &name) {}
    inline double cal_base_qual(char & qual) {
        auto v = 1 - pow(0.1, double(qual - BASE_OFFSET) / 10);
        return  v;
    }

public:
    explicit FragmentReader(const char *file_name);
    ~FragmentReader();
    inline void set_prev_chr_var_count(uint32_t prev_chr_var_count) {this->prev_chr_var_count = prev_chr_var_count;}
    inline void set_curr_chr_var_count(uint32_t curr_chr_var_count) {this->curr_chr_var_count = curr_chr_var_count;}
    inline void set_window_info(uint32_t s, uint32_t e, uint32_t intended_e)
    {
        this->window_start = s + prev_chr_var_count;
        this->window_end = e + prev_chr_var_count;
        this->intended_window_end = intended_e + prev_chr_var_count;
        this->nxt_window_set = false;
    }
    bool get_next_pe(Fragment &fragment);
    bool get_next_tenx(Fragment &fragment);
    bool get_next_hic(Fragment &fragment);
    bool get_next_pacbio(Fragment &fragment);
    bool get_next_nanopore(Fragment &fragment);
    bool get_next_matrix(std::vector<Fragment>& frags);
    bool get_next_hybrid(Fragment &fragment);
    
};

struct FragStat
{
    uint32_t start;
    uint32_t end;
    std::string barcode;
};

struct Range
{
    uint start;
    uint end;
    Range(uint start, uint end);
    bool operator==(const Range &rhs) const;
};

struct RangeHasher
{
    std::size_t operator() (const Range &key) const;
};


class RegionFragStats
{
public:
    std::unordered_map<std::string, std::unordered_map<Range, int, RangeHasher>> frags_range;
    std::unordered_map<std::string, int> count;
    RegionFragStats() = default;
    ~RegionFragStats() = default;
    inline void clear()
    {
        this->frags_range.clear();
        this->count.clear();
    }
    inline void insert(FragStat *fragStat)
    {
        if (frags_range.count(fragStat->barcode) == 0)
        {
            frags_range[fragStat->barcode];
            count[fragStat->barcode] = 1;
            auto &temp = frags_range[fragStat->barcode];
            temp[Range(fragStat->start, fragStat->end)] = 1;
        }
        else
        {
            auto &temp = frags_range[fragStat->barcode];
            temp[Range(fragStat->start, fragStat->end)] = ++count[fragStat->barcode];
        }
    }
    inline int get_index_for_barcode_variant(std::string &barcode, uint start_pos)
    {
        auto &ranges = this->frags_range[barcode];
        for (auto &info : ranges)
        {
            const Range &range = info.first;
            if (start_pos > range.start && start_pos < range.end)
                return info.second;
        }
        return 0;
    }
};

class BEDReader
{
private:
    tbx_t *tbx;
    htsFile *inf;
    hts_itr_t *iter;

public:
    kstring_t tmp;

public:
    BEDReader(const char *in);
    ~BEDReader();
    int jump_to_region(const char *chromo, uint begin, uint end);
    int read_region_frag(const char *chromo, uint begin, uint end, RegionFragStats *region_frag_stats);
    inline int get_next_record(FragStat *stat)
    {
        int ret = tbx_itr_next(inf, tbx, iter, &tmp);
        if (ret < 0)
            return ret;

        std::string buffer = std::string(tmp.s);
        std::vector<std::string> token;
        tokenize(buffer, token, "\t");
        if (token.size() <= 2)
            return -1;
        stat->start = std::stoi(token[1]);
        stat->end = std::stoi(token[2]);
        stat->barcode = token[3];
        return 0;
    }
};



#endif //PHASER_FRAG_IO_H
