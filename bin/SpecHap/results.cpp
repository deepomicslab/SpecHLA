//
// Created by yyh on 1/7/2019.
//

#include "results.h"
#include <cmath>

//------------------------------- Result Class ------------------------------------//

ResultforSingleVariant::ResultforSingleVariant()
        :phased(false), filter(filter_type ::PASS), ref_count(0), alt_count(0),
         L0_G0(0), L0_G1(0), L1_G1(0), L1_G0(0), ad0(1), ad1(1), rescued(0), not_rescued(0), MMD(0),
         nMMD(0), filtered(false)
{
    
}

ResultforSingleVariant::ResultforSingleVariant(bcf1_t *record)
        :phased(false), filter(filter_type ::PASS), ref_count(0), alt_count(0), pos(record->pos / 1000),
        L0_G0(0), L0_G1(0), L1_G1(0), L1_G0(0), ad0(1), ad1(1), rescued(0), not_rescued(0), MMD(0),nMMD(0)
        , filtered(false)
{
    if (record->qual == NAN)
        this->qual = 40;
}

ResultforSingleVariant::ResultforSingleVariant(const double &qual, const uint &pos, const int &dp)
        : phased(false), filter(filter_type ::PASS), ref_count(0), alt_count(0), pos(pos / 1000), dp(dp),
        L0_G0(0), L0_G1(0), L1_G1(0), L1_G0(0), ad0(1), ad1(1), rescued(0), not_rescued(0), MMD(0), nMMD(0)
        , filtered(false)
{   a.reset();
    if (qual == NAN)
        this->qual = 40;
    if (dp == NAN)
        this->filter = filter_type::LOW_COVERAGE;
}

ResultforSingleVariant::ResultforSingleVariant(const ResultforSingleVariant &rhs)
        : phased(rhs.phased), filter(rhs.filter), ref_count(rhs.ref_count), alt_count(alt_count), qual(rhs.qual),
        L0_G0(rhs.L0_G0), L0_G1(rhs.L0_G1), L1_G0(rhs.L1_G0), L1_G1(rhs.L1_G1), pos(rhs.pos), dp(rhs.dp), ad0(rhs.ad0), ad1(rhs.ad1),
          rescued(rhs.rescued), not_rescued(rhs.not_rescued), MMD(rhs.MMD), nMMD(rhs.nMMD), filtered(rhs.filtered)
{
    if (rhs.a.test(0))
        a.reset();
    else
        a.set();

}

bool ResultforSingleVariant::filter_inconsistency()
{
    if (!phased)
        return false;

    double score, p = 0.001;
    int good = L0_G0 + L1_G1;
    int bad = L0_G1 + L1_G0;
    int total = bad + good;
    if (total == 0)
        return false;
    score = good * log(1 - p) + bad * log(p) - total * log(0.5);
    //filter
    if (score < 0)
    {
        this->phased = false;
        this->filter = filter_type::TENXINCONSISTENCY;
        this->set_hap_0();
        this->block.reset();
        return true;
    }
    return false;
}

//------------------------------- Result Class End ------------------------------------//

