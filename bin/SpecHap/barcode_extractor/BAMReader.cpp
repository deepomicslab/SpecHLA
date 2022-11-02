//
// Created by yyh on 8/7/2019.
//

#include <cstdlib>
#include <iostream>
#include "BAMReader.h"

BAMReader::BAMReader(const char *in_name)
: inf(nullptr), header(nullptr), index(nullptr), iter(nullptr), buffer(nullptr)
{
    inf = sam_open(in_name, "r");
    if (inf == nullptr)
    {
        std::cerr << "extractFrags: Error reading bam file, check whether it exists." << std::endl;
        exit(1);
    }
    index = bam_index_load(in_name);
    if (index == nullptr)
    {
        std::cerr << "extractFrags: Error loading bam index, check whether it exists." << std::endl;
        exit(1);
    }
    header = sam_hdr_read(inf);
    if (header == nullptr)
    {
        std::cerr << "extractFrags: Error loading bam header, check whether it is corrupted." << std::endl;
        exit(1);
    }
    buffer = bam_init1();
}

BAMReader::~BAMReader()
{
    sam_close(inf);
    bam_hdr_destroy(header);
    bam_itr_destroy(iter);
    hts_idx_destroy(index);
    bam_destroy1(buffer);
}

int BAMReader::seek_to_contig(const char *chromosome)
{
    curr_tid = bam_name2id(header, chromosome);
    bam_itr_destroy(iter);
    iter = nullptr;
    iter = bam_itr_querys(index, header, chromosome);
    if (iter == nullptr)
        return -1;
    return 0;
}

int BAMReader::read_into_struct(AlignedReads *&aligned_reads)
{
    if (this->iter == nullptr)
        return -1;
    int ret = bam_itr_next(inf, iter, buffer);
    if (buffer->core.tid != curr_tid)
        return -1;
    if (ret <0)
        return ret;
    if ((buffer->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)))
        aligned_reads = nullptr;
    else
        aligned_reads = new AlignedReads(buffer);
    return 0;
}

AlignedReads::AlignedReads(bam1_t *record)
{
    span = 0;
    const bam1_core_t &data = record->core;
    this->read_len = data.l_qseq;
    this->pos = data.pos;
    this->flag = data.flag;
    this->tid = data.tid;
    this->mtid = data.mtid;
    this->mapping_quality = data.qual;
    const char *bar = (char*) bam_aux_get(record, "BX");

    bar == nullptr?barcode = "":this->barcode = std::string(bar);

    uint32_t *cigar = bam_get_cigar(record);
    int op, opl;
    for (auto i = 0; i < data.n_cigar; i++)
    {
        op = bam_cigar_op(cigar[i]);
        opl = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
        {
            this->span += opl;
        }
        else if (op == BAM_CDEL)
        {
            this->span += opl;
        }
        else if (op == BAM_CINS);
        else if (op == BAM_CREF_SKIP)
            this->span += opl;
        else if (op == BAM_CSOFT_CLIP);
        else if (op == BAM_CHARD_CLIP);
        else
            ;
    }

    end_pos = pos + span;
}