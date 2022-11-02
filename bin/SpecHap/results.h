//
// Created by yyh on 1/7/2019.
//

#ifndef SPECHAP_RESULTS_H
#define SPECHAP_RESULTS_H
#include <bitset>
#include <memory>
#include <htslib/vcf.h>

enum class filter_type
        {
            PASS,
            NOINFO,
            POOLRESULT,
            LOW_COVERAGE,
            TENXINCONSISTENCY,
            CONFILCTINGRESULT,
            TENX_ALLELE_FRENCENCY_FILTER,
            TENX_QUAL_FILTER,
            TENX_RESCUED_MOLECUE_HIGH_DIVERSITY
        };



typedef std::bitset<1> shap;
typedef unsigned int uint;

template <typename T>
bool is_uninitialized(std::weak_ptr<T> const& weak)
{
    using wt = std::weak_ptr<T>;
    return !weak.owner_before(wt{}) && !wt{}.owner_before(weak);
}


class PhasedBlock;

class ResultforSingleVariant
{
private:
    bool phased;
    filter_type filter;
    shap a;                                     //Ref: 0, Alt: 1

public:
    std::weak_ptr<PhasedBlock> block;
    int ref_count, alt_count;
    int ad0, ad1;
    double qual;
    int dp;
    double af;
    uint pos;
    uint ps;                                    //phased block no;
    int L0_G0, L0_G1, L1_G1, L1_G0;            //L0: local hap 0, L1: local hap 1, G1: global hap 1, G0: global hap 0
    double MMD;
    int rescued, not_rescued;
    int nMMD;
    bool filtered;
    
public:
    //constructor
    ResultforSingleVariant();
    ResultforSingleVariant(bcf1_t *record);
    ResultforSingleVariant(const double &qual, const uint &pos, const int &dp);
    ResultforSingleVariant(const ResultforSingleVariant &rhs);
    ~ResultforSingleVariant() = default;


    //getter
    inline bool variant_phased() const { return phased; }
    inline filter_type get_filter() const { return filter; }
    inline bool is_REF() const { return  !a.test(0); }
    inline bool is_ALT() const { return  a.test(0); }
    inline uint get_pos() const {return pos;}

    //setter
    inline void complement_hap() {a.flip();}
    inline void set_hap_0() {a.reset();}
    inline void set_hap_1() {a.set();}
    inline void update_phased_block(std::shared_ptr<PhasedBlock> phasedBlock)
    {
        this->block.reset();
        this->block = phasedBlock;
    }
    inline void set_hap(bool is_ref)
    {
        if (is_ref)
            a.reset();
        else
            a.set();
    }
    inline int get_hap() {if (is_REF()) return 0; else return 1;}
    inline void set_phased() {phased = true;}
    inline void set_unphased() {phased = false;}
    inline void set_filter(filter_type filter) {this->filter = filter;}
    bool filter_inconsistency();

};

#endif //SPECHAP_RESULTS_H