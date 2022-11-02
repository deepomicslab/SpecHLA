//
// Created by yyh on 8/7/2019.
//

#include "Barcode.h"

Barcode::Barcode(const Barcode &rhs)
{
    this->barcode = rhs.barcode;
    this->start_pos = rhs.start_pos;
    this->end_pos = rhs.end_pos;
    this->last_read_pos = rhs.last_read_pos;
    this->last_q30_read_end_pos = rhs.last_q30_read_end_pos;
    reads_count = rhs.reads_count;
    confident_reads_count = rhs.reads_count;
}
//involve this constructor when read_qual > q30
Barcode::Barcode(const std::string &barcode, const int &start_pos, const int &end_pos, const int &last_read_pos)
{
    this->barcode = barcode;
    this->start_pos = start_pos;
    this->end_pos = end_pos;
    this->last_q30_read_end_pos = end_pos;
    this->last_read_pos = last_read_pos;
    this->reads_count = 1;
    this->confident_reads_count = 0;
}

int Barcode::extend(AlignedReads *read)
{
    if (read->pos > this->last_read_pos)
        this->last_read_pos = read->pos;
    if (read->end_pos > this->end_pos)
        this->end_pos = read->end_pos;
    if (read->mapping_quality >= 30 && read->end_pos > this->last_q30_read_end_pos)
        this->last_q30_read_end_pos = read->end_pos;
    reads_count++;
    return 0;
}