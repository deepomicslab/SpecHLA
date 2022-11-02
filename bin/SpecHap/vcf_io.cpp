//
// Created by yonghanyu2 on 10/10/2018.
//

#include "vcf_io.h"
#include <cstring>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

VCFReader::VCFReader(const char *filename)
: iter(nullptr), tbx_index(nullptr)
{
    size_t l = strlen(filename);
    if (l < 1)
    {
        std::cerr << "Error: Check input vcf file" << std::endl;
        exit(-1);
    }
    this->filename = new char[l + 1];
    strcpy(this->filename, filename);
    vcf_file = bcf_open(filename, "r");
    if (vcf_file == nullptr)
    {
        std::cerr << "Error: Fail to load VCF file" << std::endl;
    }
    tbx_index = tbx_index_load(filename);
    if (tbx_index == nullptr)
    {
        std::cerr<< "Error: Fail to load tabix file, check whether it exists." << std::endl;
        exit(-1);
    }
    header = bcf_hdr_read(vcf_file);
    if (header == nullptr)
    {
        std::cerr << "Error: corrupted vcf header." << std::endl;
        exit(-1);
    }
    tmp = {0, 0, nullptr};
    get_contigs();
    buffer = bcf_init();
    //count_contigs();
    //reset();
}

VCFReader::~VCFReader()
{
    delete []filename;
    bcf_close(vcf_file);
    bcf_hdr_destroy(header);
    contigs.clear();
    var_count.clear();
    hts_itr_destroy(iter);
    tbx_destroy(tbx_index);
    free(tmp.s);
    bcf_destroy(buffer);
}

void VCFReader::read_into_struct(bcf1_t *record, ptr_ResultforSingleVariant &result)
{
    record->qual == NAN ? result->qual = 40: result->qual = record->qual;
    result->pos = record->pos / 1000;
    int *dps = nullptr; int dp_n;
    bcf_get_format_int32(header, record, "DP", &dps, &dp_n);
    result->dp = dps[0];
    free(dps);
}

template<typename T> struct map_init_helper
{
    T& data;
    explicit map_init_helper(T& d) : data(d) {}
    map_init_helper& operator() (typename T::key_type const& key, typename T::mapped_type const& value)
    {
        data[key] = value;
        return *this;
    }
};

template<typename T> map_init_helper<T> map_init(T& item)
{
    return map_init_helper<T>(item);
}

VCFWriter::VCFWriter(const bcf_hdr_t *hdr, const char *outfile)
{

    fp = bcf_open(outfile, "w");
    map_init(this->FilterN) (filter_type::NOINFO, "INFO_NOT_ENOUGH") (filter_type::POOLRESULT, "POOL_SPECTRAL_RESULT") (filter_type::TENXINCONSISTENCY, "10X_PHASING_INCONSISTENCY") (filter_type::PASS, "PASS");
    map_init(this->FormatN) (GT, "GT") (PS, "PS");

    header = bcf_hdr_dup(hdr);
    header_init();
    ngt = bcf_hdr_nsamples(hdr);

//    gt = new int[2*ngt];
    gt = (int*)malloc(ngt*2*sizeof(int));
}

void VCFWriter::header_init()
{
    char GT[] = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    char PS[] = "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing Block No.\">";
    char FP[] = "##FILTER=<ID=PASS,Description=\"All filters passed\">";
    char FI[] = "##FILTER=<ID=INFO_NOT_ENOUGH,Description=\"Provided Information is not enough for phasing\">";
    char FPR[] = "##FILTER=<ID=POOL_SPECTRAL_RESULT,Description=\"Phasing results with spetral graph theory is poor to be used\">";
    char FP1[] = "##FILTER=<ID=10X_PHASING_INCONSISTENCY,Description=\"Phasing results are inconsistent\">";
    char FP2[] = "##FILTER=<ID=WINDOW_RESULT_INCONSISTENCY,Description=\"Phasing results are inconsistent\">";
    char FP3[] = "##FILTER=<ID=10X_ALLELE_FREQUENCY_FILTER,Description=\"Phasing results are inconsistent\">";
    char FP4[] = "##FILTER=<ID=LOW_COVERAGE,Description=\"Low sequencing coverage\">";
    char FP5[] = "##FILTER=<ID=TENX_QUAL_FILTER,Description=\"Low qual\">";
    char FP6[] = "##FILTER=<ID=TENX_RESCUED_MOLECUE_HIGH_DIVERSITY,Description=\"Filter by mapq rescue and mean molecule divergence\">";

    std::string line_;
    kstring_t kstring = {0, 0, nullptr};
    bcf_hdr_format(header, 0, &kstring);
    char *htxt = new char[kstring.l + 5000];
    char *hdr_txt = ks_release(&kstring);
    std::stringstream ss(hdr_txt);
    int str_l;
    bool Format = false; bool Filter = false;

    std::getline(ss, line_, '\n');
    strcpy(htxt, line_.data()); strcat(htxt, "\n");

    while(std::getline(ss, line_, '\n'))
    {
        if (line_[0] != '#')
            break;
        bcf_hrec_t *hrec = bcf_hdr_parse_line(header, line_.c_str(), &str_l);
        if (hrec == nullptr)
        {
            strcat(htxt, line_.c_str());
            break;
        }
        else
        {
            if (strcmp(hrec->key, "FILTER") == 0)
            {
                if (!Filter)
                {
                    Filter = true;
                    strcat(htxt, FP);
                    strcat(htxt, "\n");
                    strcat(htxt, FI);
                    strcat(htxt, "\n");
                    strcat(htxt, FPR);
                    strcat(htxt, "\n");
                    strcat(htxt, FP1);
                    strcat(htxt, "\n");
                    strcat(htxt, FP2);
                    strcat(htxt, "\n");
                    strcat(htxt, FP3);
                    strcat(htxt, "\n");
                    strcat(htxt, FP4);
                    strcat(htxt, "\n");
                    strcat(htxt, FP5);
                    strcat(htxt, "\n");
                    strcat(htxt, FP6);
                    strcat(htxt, "\n");
                }
                if ((hrec->nkeys != 0) && (strcmp(hrec->vals[0], "PASS") == 0))
                {
                    bcf_hrec_destroy(hrec);
                    continue;
                }
                strcat(htxt, line_.c_str());
                strcat(htxt, "\n");
            }
                // format
            else if (strcmp(hrec->key, "FORMAT") == 0)
            {
                if (!Format)
                {
                    Format = true;
                    strcat(htxt, GT);
                    strcat(htxt, "\n");
                    strcat(htxt, PS);
                    strcat(htxt, "\n");
                }
                if (hrec->nkeys != 0 && strcmp(hrec->vals[0], "GT") == 0)
                {
                    bcf_hrec_destroy(hrec);
                    continue;
                }
                if (hrec->nkeys != 0 && strcmp(hrec->vals[0], "PS") == 0)
                {
                    bcf_hrec_destroy(hrec);
                    continue;
                }
                strcat(htxt, line_.c_str());
                strcat(htxt, "\n");
            }
            else
            {
                strcat(htxt, line_.c_str());
                strcat(htxt, "\n");
            }
        }
        bcf_hrec_destroy(hrec);
    }
    bcf_hdr_destroy(header);
    header = bcf_hdr_init("w");
    bcf_hdr_parse(header, htxt);
    int rr = bcf_hdr_write(fp, header);

    delete []htxt;
    free(hdr_txt);
}

VCFWriter::~VCFWriter()
{
    bcf_hdr_destroy(header);
    bcf_close(fp);
    delete [] gt;
}

void VCFWriter::write_nxt_record(bcf1_t *record, ptr_ResultforSingleVariant resultforSingleVariant, const unsigned int blk_no)
{
    if (record->pos == 126330351) {
        int tmp = 33;
    }
    const char * filter = filter_map.find(resultforSingleVariant->get_filter())->second.data();
    int ref;
    int alt;
    int32_t *gt_arr = NULL, ngt_arr = 0;
    bcf_get_genotypes(header, record, &gt_arr, &ngt_arr);
    int32_t *ptr = gt_arr;
    int allele0 = bcf_gt_allele(ptr[2*SIDX]);
    int allele1 = bcf_gt_allele(ptr[2*SIDX+1]);
    int temp = bcf_alleles2gt(allele0, allele1);;
    bcf_gt2alleles(temp, &ref, &alt);
//    ref = allele0;
//    alt = allele1;
    for(int i = 0; i < ngt_arr; i++) {
        gt[i] = gt_arr[i];
    }
    if (resultforSingleVariant->is_REF())
    {
        if (resultforSingleVariant->variant_phased() && blk_no != 0)
        {
            gt[2*SIDX] = bcf_gt_phased(ref);
            gt[2*SIDX+1] = bcf_gt_phased(alt);
        }
        else
        {
            gt[2*SIDX] = bcf_gt_unphased(ref);
            gt[2*SIDX+1] = bcf_gt_unphased(alt);
        }

    }
    else
    {
        if (resultforSingleVariant->variant_phased())
        {
            gt[2*SIDX] = bcf_gt_phased(alt);
            gt[2*SIDX+1] = bcf_gt_phased(ref);
        }
        else
        {
            gt[2*SIDX] = bcf_gt_unphased(alt);
            gt[2*SIDX+1] = bcf_gt_unphased(ref);
        }
    }

    int k = bcf_hdr_id2int(this->header, BCF_DT_ID, filter);
    bcf_update_filter(this->header, record, &k, 1);
    bcf_update_genotypes(this->header, record, gt, 2*ngt);
    
    if (resultforSingleVariant->variant_phased()) {
        int *ps_s = nullptr;
        int nps = 0, ps =0;
        nps = bcf_get_format_int32(this->header ,record, "PS", &ps_s, &ps);
        if (nps <=0 ) {
            ps_s = (int*)malloc(ngt*2*sizeof(int));
            for (int i = 0; i < ngt; ++i) {
                ps_s[i] = 0;
            }
        }else {
            for (int i = 0; i < ngt; ++i) {
                ps_s[i] = (ps_s[i] >= INT_MAX || ps_s[i] < 0) ? 0:ps_s[i];
            }
        }
        ps_s[SIDX] = blk_no;
        bcf_update_format_int32(this->header, record, "PS", ps_s, ngt);
        if (ps_s != nullptr)
        free(ps_s);
    } else {
        int *ps_s = nullptr;
        int nps = 0, ps =0;
        nps = bcf_get_format_int32(this->header ,record, "PS", &ps_s, &ps);
        if (nps <= 0) {
            ps_s = (int*)malloc(ngt*2*sizeof(int));
            for (int i = 0; i < ngt; ++i) {
                ps_s[i] = 0;
            }
        } else {
            for (int i = 0; i < ngt; ++i) {
                ps_s[i] = (ps_s[i] >= INT_MAX || ps_s[i] < 0) ? 0:ps_s[i];
            }
        }
        bcf_update_format_int32(this->header, record, "PS", ps_s, ngt);
        if (ps_s != nullptr)
            free(ps_s);
    }

    bcf1_t *record_w = bcf_dup(record);
    //bcf_unpack(record_w, BCF_UN_ALL);
    int rr = bcf_write(this->fp, this->header, record_w);
    bcf_destroy(record_w);
    free(gt_arr);
}



void VCFWriter::write_nxt_contigs(const char *contig, ChromoPhaser *chromo_phaser, VCFReader &frvcf, const std::set<uint> &break_idx)
{
    bcf1_t *record = bcf_init();
    uint blk_count = 0;
    int gap_count = 0;
    std::unordered_map<ptr_PhasedBlock, uint > encountered_phased_block;
    frvcf.jump_to_contig(frvcf.curr_contig);
    std::unordered_map<uint, uint> idx2pos;
    std::unordered_map<uint, uint> blk_no2_nblk_no;
    for (uint idx = 0; idx < chromo_phaser->variant_count; idx++)
    {
        frvcf.get_next_record(record);
        idx2pos[idx] = record->pos + 1;
        ptr_ResultforSingleVariant resultforSingleVariant = chromo_phaser->results_for_variant[idx];
        bcf_translate(this->header, frvcf.header, record);
        if (record->pos == 4565) {
            auto tmp = 44;
        }
        if (!is_uninitialized(resultforSingleVariant->block) && resultforSingleVariant->variant_phased() && resultforSingleVariant->get_filter() == filter_type::PASS)
        {
            if (gap_count >= 30)
            {
                blk_count++;
            }
            ptr_PhasedBlock block = resultforSingleVariant->block.lock();
            if (break_idx.find(idx) != break_idx.end()) {
                blk_no2_nblk_no.emplace(idx2pos[block->start_variant_idx], record->pos + 1);
            }
            if (block->results.size() == 1) {
                int tmp = 2;
            }
            if (block.get() == nullptr)
            {
                write_nxt_record(record, resultforSingleVariant, 0);
                continue;
            }
            //already met
            int mm = encountered_phased_block.count(block);
            if (encountered_phased_block.count(block) != 0)
            {
                uint blk_no = encountered_phased_block[block];
                uint cur_blk_no = idx2pos[block->start_variant_idx];
                if (blk_no2_nblk_no.find(cur_blk_no) != blk_no2_nblk_no.end() && record->pos +1 >= blk_no2_nblk_no[cur_blk_no]) {
                    auto new_blk_no = blk_no2_nblk_no[cur_blk_no];
                    write_nxt_record(record, resultforSingleVariant, new_blk_no);
                } else {
                    write_nxt_record(record, resultforSingleVariant, cur_blk_no);
                }
            }
            else
            {
                encountered_phased_block.emplace(block, ++blk_count);
                if(block->results.size() == 1) {
                    write_nxt_record(record, resultforSingleVariant, 0);
                } else{
                    if (blk_no2_nblk_no.find(idx2pos[block->start_variant_idx]) != blk_no2_nblk_no.end() && record->pos >= blk_no2_nblk_no[idx2pos[block->start_variant_idx]]) {
                        auto new_blk_no = blk_no2_nblk_no[idx2pos[block->start_variant_idx]];
                        write_nxt_record(record, resultforSingleVariant, new_blk_no);
                    } else {
                        write_nxt_record(record, resultforSingleVariant, idx2pos[block->start_variant_idx]);
                    }
                }
            }

            gap_count = 0;
        }
        else
        {
            gap_count++;
            write_nxt_record(record, resultforSingleVariant, 0);
        }
    }
    encountered_phased_block.clear();
    bcf_destroy(record);
}

