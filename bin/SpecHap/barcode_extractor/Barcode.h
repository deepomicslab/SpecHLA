//
// Created by yyh on 8/7/2019.
//

#ifndef BARCODE_EXTRACTOR_BARCODE_H
#define BARCODE_EXTRACTOR_BARCODE_H


#include <string>
#include <unordered_map>
#include "BAMReader.h"

class Barcode
{
public:
    std::string barcode;
    int reads_count;
    int confident_reads_count;
    uint32_t start_pos;
    uint32_t end_pos;
    uint32_t last_read_pos;
    uint32_t last_q30_read_end_pos;

public:
    Barcode(const std::string &barcode, const int &start_pos, const int &end_pos, const int &last_read_pos);
    Barcode(const Barcode &rhs);
    ~Barcode() = default;
    int extend(AlignedReads *read);
};




#endif //BARCODE_EXTRACTOR_BARCODE_H
