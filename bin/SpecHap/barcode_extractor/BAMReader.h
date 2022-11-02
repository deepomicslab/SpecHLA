//
// Created by yyh on 8/7/2019.
//

#ifndef BARCODE_EXTRACTOR_BAMREADER_H
#define BARCODE_EXTRACTOR_BAMREADER_H
#include "htslib/sam.h"
#include <string>
class AlignedReads
{
public:
    int read_len;
    std::string barcode;
    uint32_t pos, mpos;
    int flag;
    int mapping_quality;
    int tid, mtid;
    int span;
    uint32_t end_pos;

public:
    ~AlignedReads() = default;
    explicit AlignedReads(bam1_t *record);
    inline bool mapped()
    {
        return ((flag & BAM_FUNMAP) == 0) && ((flag & BAM_FSUPPLEMENTARY) == 0) && ((flag & BAM_FSUPPLEMENTARY) == 0);
    }
    inline bool useful()
    {
        return (!this->barcode.empty()) && mapped();
    }
};

class BAMReader
{
public:
    samFile *inf;
    bam_hdr_t *header;
    hts_idx_t *index;
    hts_itr_t *iter;
    bam1_t *buffer;
    int curr_tid;
public:
    explicit BAMReader(const char *in_name);
    ~BAMReader();
    int get_next_contig_read();
    int seek_to_contig(const char *chromosome);
    int read_into_struct(AlignedReads *&aligned_reads);
};


#endif //BARCODE_EXTRACTOR_BAMREADER_H
