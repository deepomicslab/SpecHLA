//
// Created by yyh on 8/8/2019.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include "worker.h"

worker::worker(const char *infile, const char *out)
{
    bam_reader = new BAMReader(infile);
    ofs.open(out, std::ostream::out);
}

worker::~worker()
{
    delete bam_reader;
    ofs.close();
    frags.clear();
}



int worker::extract_one_target(const char *target)
{
    if (bam_reader->seek_to_contig(target) != 0)
        return -1;

    while (true)
    {
        AlignedReads *alignedReads = nullptr;
        int ret = bam_reader->read_into_struct(alignedReads);
        if (ret != 0)   //target end
            break;
        if (alignedReads == nullptr)    // reads not valid;
            continue;
        if (alignedReads->barcode.empty())
        {
            delete alignedReads;
            continue;
        }
        this->extract_read_barcode(alignedReads);
        delete alignedReads;
    }
    write_frags(target);
    this->frags.clear();
    return 0;
}

int worker::extract()
{
    for (int i = 0; i < bam_reader->header->n_targets; i++)
    {
        const char *target = bam_reader->header->target_name[i];
        if (bam_reader->seek_to_contig(target) != 0)
            continue;
        std::cout << "Processing reads on Chromosomo " << target << std::endl;
        extract_one_target(target);
    }
    return 0;
}

int worker::extract_read_barcode(AlignedReads *read)
{
    std::string barcode = read->barcode;
    auto range = this->frags.equal_range(barcode);
    if (range.first == range.second)    // not encountered
    {
        if (read->mapping_quality >= 30)
            frags.insert(std::make_pair(barcode, Barcode(read->barcode, read->pos, read->end_pos, read->pos)));
    }
    else                                // already met, check distant
    {
        Barcode &frag = (--range.second)->second;
        if (frag.last_read_pos + 60000 < read->pos) // exceed max range, create new frag
        {
            if (read->mapping_quality >= 30)
                frags.insert(std::make_pair(barcode, Barcode(read->barcode, read->pos, read->end_pos, read->pos)));
        }
        else
        {
            frag.extend(read);
        }
    }
    return 0;
}

void worker::write_frags(const char *target)
{
    typedef std::pair<std::string, Barcode> val_pair;
    std::vector<val_pair> temp;
    for (auto &i: frags)
        temp.push_back(i);
    std::sort(temp.begin(), temp.end(), [=](val_pair &a, val_pair &b)
    {
        return a.second.start_pos < b.second.start_pos;
    });

    for (auto &i: temp)
    {
        Barcode &frag = i.second;
        if (frag.reads_count == 1)
            continue;
        ofs << target << "\t";
        ofs << frag.start_pos << "\t" << frag.last_q30_read_end_pos << "\t" << frag.barcode.substr(1) << "\n";
    }
}