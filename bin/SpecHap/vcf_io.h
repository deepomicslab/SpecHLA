//
// Created by yonghanyu2 on 10/10/2018.
//

#ifndef SPECHAP_VCF_IO_H
#define SPECHAP_VCF_IO_H

#include "htslib/vcf.h"
#include <map>
#include "spectral.h"
#include <vector>
#include <string>
#include "htslib/tbx.h"
#include "results.h"
extern int SIDX;
typedef std::map<filter_type , std::string> filter_map_t;
//const int DOT_PS = 2147483648;

enum format_n{GT, PS};
class VCFReader
{
public:
    htsFile *vcf_file;
    bcf_hdr_t *header;
    int curr_contig, curr_bcf_contig;
    char *filename;
    int contigs_count;
    std::vector<std::string> contigs;   //rid to contig name appeared in header
    std::map<std::string, int> bcf_contig_map;
    std::map<uint, uint> var_count;      //rid to contig_variant_count appeared in vcf file
    hts_itr_t *iter;
    tbx_t *tbx_index;
    kstring_t tmp;
    
    bcf1_t *buffer;

public:
    explicit VCFReader(const char *filename);
    ~VCFReader();

    inline void reset()
    {
        bcf_close(vcf_file);
        bcf_hdr_destroy(header);
        vcf_file = bcf_open(filename, "r");
        header = bcf_hdr_read(vcf_file);
        curr_contig = 0;
        hts_itr_destroy(iter);
        tbx_destroy(tbx_index);
    }

    inline const void get_contigs()
    {
        const char **temp = tbx_seqnames(tbx_index, &contigs_count);
        for (int i = 0; i < contigs_count; i++)
        {
            contigs.emplace_back(std::string(temp[i]));
            temp[i] = nullptr;
        }
        free(temp);

        int bcf_contig_count;
        temp = bcf_hdr_seqnames(header, &bcf_contig_count);
        for (int i = 0; i < bcf_contig_count; i++)
        {
            bcf_contig_map[std::string(temp[i])] = i;
            temp[i] = nullptr;
        }
        free(temp);
    }

    /*
    inline int get_AF0(float *af0, bcf1_t *record)
    {
        float **AF = nullptr;
        int32_t l_AF;
        bcf_get_info_float(header, record, "AF", &AF, &l_AF);
    }*/

    inline int get_next_record(bcf1_t *record)
    {
        if (iter == nullptr)
            return -1;
        int status = tbx_itr_next(vcf_file, tbx_index, iter, &tmp);
        status = vcf_parse1(&tmp, header, record);
        if (status != 0) //eof or corruption in file
            return status;
        if (record->rid != this->curr_bcf_contig)
            return -1;
        bcf_unpack(record, BCF_UN_ALL);
        return status;
    }
    void read_into_struct(bcf1_t *record, ptr_ResultforSingleVariant &result);
    inline int get_next_record_contig(ptr_ResultforSingleVariant &result, bool use_gt)
    {
        //bcf_empty(record);
        //no such contig or iter not initialized
        if (iter == nullptr)
            return -1;   
        //bcf1_t *r = bcf_init();
        int status = tbx_itr_next(vcf_file, tbx_index, iter, &tmp);
        if (status < 0)
            return  status;
        status = vcf_parse1(&tmp, header, buffer);
        if (status != 0) //eof or corruption in file
            return -1;
        //new contig
        if (buffer->rid != this->curr_bcf_contig)
            return -1;
        //bcf_subset_format(header, record);
        bcf_unpack(buffer, BCF_UN_ALL);

        int dp = 0, ndp;
        int af = 0, naf;
        int ad = 0, nad;
        int ps = 0, nps;
        int *dps = nullptr, *ads = nullptr;
        int *ps_s = nullptr;                  //PS
        float *afs = nullptr;
        ndp = bcf_get_info_int32(header, buffer, "DP", &dps, &dp);
        naf = bcf_get_info_float(header, buffer, "AF", &afs, &af);
        nad = bcf_get_format_int32(header ,buffer, "AD", &ads, &ad);
        nps = bcf_get_format_int32(header ,buffer, "PS", &ps_s, &ps);
        buffer->qual == NAN ? result->qual = 40 : result->qual = buffer->qual;
        result->pos = buffer->pos;
        if (buffer->pos == 126330351) {
            int tmp = 33;
        }
        dps == nullptr ? result->dp = 30 : result->dp = dps[SIDX];
        afs == nullptr ? result->af = 0.5 : result->af = afs[SIDX];
        ps_s == nullptr ? result->ps = 0: result->ps = ps_s[SIDX];
        if (ads != nullptr)
        {
            if (nad % 2 == 0)
            {
                result->ad0 = ads[2*SIDX];
                result->ad1 = ads[2*SIDX+1];
            }
            if (nad % 3 == 0)
            {
                result->ad0 = ads[2*SIDX+1];
                result->ad1 = ads[2*SIDX+2];
            }
        }

        if (use_gt)
        {   
            int32_t *gt_arr = nullptr, ngt_arr = 0;
            bcf_get_genotypes(header, buffer, &gt_arr, &ngt_arr);
            int32_t *ptr = gt_arr;
            int allele0 = bcf_gt_allele(ptr[2*SIDX]);
            int allele1 = bcf_gt_allele(ptr[2*SIDX+1]);
            

            if (allele1 == allele0) //homozy
            {
                free(gt_arr);
                return 1;
            }
            if (ps_s != nullptr)
            {
                if (allele0 == 2 || allele1 == 2) {
                    if (allele0 == 0 || allele0 == 1)
                        result->set_hap_0();
                    else
                        result->set_hap_1();
                } else {
                    if (allele0 == 0)
                        result->set_hap_0();
                    else if (allele0 == 1)
                        result->set_hap_1();
                }
                result->set_phased();
                //std::cout<< allele0 << ";" << allele1 << std::endl;
            }
            free(gt_arr);
            
        }
        free(dps); free(ads); free(ps_s); free(afs);

        //bcf_destroy(r);
        return 0;
    }

    inline int jump_to_contig(int contig_id)
    {
        curr_contig = contig_id;
        tbx_itr_destroy(iter);
        iter = nullptr;
        iter = tbx_itr_querys(tbx_index, contigs[contig_id].c_str());
        if (iter == nullptr) // no such contig
            return -1;
        curr_bcf_contig = bcf_contig_map[contigs[curr_contig]];
        return 0; //success
    }
};

class VCFWriter
{
private:
    bcf_hdr_t *header;
    htsFile *fp;
    std::map<filter_type, const char*> FilterN;
    std::map<format_n, const char*> FormatN;
    filter_map_t filter_map =
        {
                {filter_type::NOINFO, "INFO_NOT_ENOUGH"},
                {filter_type::POOLRESULT, "POOL_SPECTRAL_RESULT"},
                {filter_type::TENXINCONSISTENCY, "10X_PHASING_INCONSISTENCY"},
                {filter_type::PASS, "PASS"},
                {filter_type::CONFILCTINGRESULT, "WINDOW_RESULT_INCONSISTENCY"},
                {filter_type::TENX_ALLELE_FRENCENCY_FILTER, "10X_ALLELE_FREQUENCY_FILTER"},
                {filter_type::LOW_COVERAGE, "LOW_COVERAGE"},
                {filter_type::TENX_QUAL_FILTER, "TENX_QUAL_FILTER"},
                {filter_type ::TENX_RESCUED_MOLECUE_HIGH_DIVERSITY, "TENX_RESCUED_MOLECUE_HIGH_DIVERSITY"}
        };
    int ngt;
    int *gt;

private:
    void header_init();

public:
    VCFWriter(const bcf_hdr_t *hdr, const char *outfile);
    ~VCFWriter();
    void write_nxt_record(bcf1_t *record, ptr_ResultforSingleVariant resultforSingleVariant, const unsigned int blk_no);
    void write_nxt_contigs(const char *contig, ChromoPhaser *chromo_phaser, VCFReader &frvcf, const std::set<uint> &break_idx);
};


#endif //SPECHAP_VCF_IO_H
