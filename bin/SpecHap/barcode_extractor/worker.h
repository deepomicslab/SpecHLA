//
// Created by yyh on 8/8/2019.
//

#ifndef BARCODE_EXTRACTOR_WORKER_H
#define BARCODE_EXTRACTOR_WORKER_H

#include <map>
#include <string>
#include <fstream>
#include "Barcode.h"
#include "BAMReader.h"

class worker
{
public:
    std::ofstream ofs;
    std::multimap<std::string, Barcode> frags;
    BAMReader *bam_reader;
    ~worker();
    worker(const char *infile, const char *out);
    int extract();
    int extract_one_target(const char *target);
    int extract_read_barcode(AlignedReads *read);
    void write_frags(const char * target);
};


#endif //BARCODE_EXTRACTOR_WORKER_H
