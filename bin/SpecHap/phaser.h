//
// Created by yyh on 3/4/2019.
//

#ifndef SPECHAP_PHASER_H
#define SPECHAP_PHASER_H

#include "type.h"
#include "vcf_io.h"
#include "frag_io.h"
#include "optionparser.h"
#include "util.h"

extern int WINDOW_SIZE;
extern int WINDOW_OVERLAP;
extern int MAX_BARCODE_SPANNING;
extern bool HYBRID;
extern bool KEEP_PS;
extern int MAX_HIC_INSERTION;
extern int OPERATION;
extern int RECURSIVE_LIMIT;

class Phaser
{
public:
    Phaser() = default;
    explicit Phaser(const std::string & fnvcf, const std::string & fnout, std::vector<std::string> & fnfrags, const std::string &fnbed, std::vector<double> & fr_weights);
    ~Phaser();
    void phasing();
    void set_contigs(std::string& contigs);

private:
//    contigs to phase
    std::vector<std::string> contigs;
    void sort_frag_file(std::string file_name);
    double threshold;
    VCFReader *frvcf;
    VCFWriter *fwvcf;
    std::vector<FragmentReader* > frfrags;
//    FragmentReader *frfrag;
    BEDReader *frbed;
    Spectral *spectral;
    int coverage;
    void phasing_by_chrom(uint var_count, ChromoPhaser *chromo_phaser);
    void phase_HiC_recursive(ChromoPhaser *chromo_phaser, std::set<uint> &connected_comp);
    void phase_HiC_poss(ChromoPhaser *chromo_phaser);
    void update_HiC_phasing_window();
    int load_contig_records(ChromoPhaser *chromo_phaser);
    int load_contig_blocks(ChromoPhaser *chromo_phaser);
};

#endif //SPECHAP_PHASER_H
