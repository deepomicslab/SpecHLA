/* functions for comparing an aligned sequence read to the set of variants to identify alleles and haplotype-informative reads */

#include <fstream>
#include <map>
#include "parsebamread.h"
#include "readvariant.h"

// only for comparing read to SNP (or block substitution with multiple SNPs that don't change length of haplotype)

int compare_read_SNP(struct alignedread* read, VARIANT* varlist, int ss, int start, int l1, int l2, FRAGMENT* fragment) {
    int offset = varlist[ss].position - start - l2;
    int j = 0;
    char match = '0';       // match -> 0 -> REF | 1 -> ALT | - -> not match
    for (j = 0; j < strlen(varlist[ss].allele1); j++) {
        if (read->sequence[l1 + offset + j] != varlist[ss].allele1[j]) {    // SNV -> gt0 ?
            match = '-';
            break;
        }
    }
    if (match == '-') {
        match = '1';
        for (j = 0; j < strlen(varlist[ss].allele2); j++) {
            if (read->sequence[l1 + offset + j] != varlist[ss].allele2[j]) {  // SNV -> gt1 ?
                match = '-';
                break;
            }
        }
    }

    if (match != '-' && read->quality[l1 + offset] - QVoffset >= MINQ && varlist[ss].type == 0) {
        fragment->alist[fragment->variants].varid = ss;
        fragment->alist[fragment->variants].allele = match;
        //assign base quality to be minimum of base quality and mapping quality
        if (read->mquality < (int) read->quality[l1 + offset] - QVoffset) fragment->alist[fragment->variants].qv = (char) (read->mquality + QVoffset);
        else fragment->alist[fragment->variants].qv = read->quality[l1 + offset];
        fragment->variants++;
        varlist[ss].depth++;
        if (match == '0') {
            if ((read->flag & 16) == 16) varlist[ss].A1 += 1 << 16;
            else varlist[ss].A1 += 1;
        } else // match == '1'
        {
            if ((read->flag & 16) == 16) varlist[ss].A2 += 1 << 16;
            else varlist[ss].A2 += 1;
        }
        if (TRI_ALLELIC == 1 && varlist[ss].heterozygous == '2') fprintf(stderr, "comparing read %s to non-ref-het SNP %d:%d %s %s | allele = %c\n", read->readid, ss + 1, varlist[ss].position, varlist[ss].allele1, varlist[ss].allele2, match);
        return 1;
    }
    return 0;
    //printf("allele %s %d %c/%c %c %c %d\t",varlist[ss].chrom,varlist[ss].position,varlist[ss].allele1,varlist[ss].allele2,read->sequence[l1+offset],read->quality[l1+offset],offset);
}


int compare_read_BND(struct alignedread *read, VARIANT *varlist, int ss, int start, int l1, int l2, FRAGMENT *fragment)
{
    // I dont think we will ever be here
    return 0;
}

//CL is cigar offset that is being evaluated 4S 60M 10I 5M  then CL = 1 for 60M, cigarlist is copy of bam->cigar

int compare_read_INDEL(struct alignedread* read, VARIANT* varlist, int ss, int start, int l1, int l2, int clength, FRAGMENT* fragment, int CL, REFLIST* reflist) {
    int pflag = 0;
    int shift = calculate_rightshift(varlist, ss, reflist);
    if (pflag) fprintf(stdout, "rightshift %d \n", shift);
    int offset = varlist[ss].position - start - l2;
    int j = 0, k = 0;
    int a1 = 0, a2 = 0;
    int s1 = 0, s2 = 0;
    int score1 = 0, score2 = 0, score3 = 0, score4 = 0;
    int b1 = 0, b2 = 0, b3 = 0, b4 = 0;
    int partialflag = 0;

    if (pflag) {
        for (j = 0; j < read->cigs; j++) fprintf(stdout, "%d%c", read->cigarlist[j] >> 4, read->cigarlist[j]&0xf);
        fprintf(stdout, " %d%c %d/%d ", read->cigarlist[CL] >> 4, read->cigarlist[CL]&0xf, CL, read->cigs);
        fprintf(stdout, "%s reference-indel-read %d-%d %s:%d:%s:%s %s\n", read->readid, start + l2, start + l2 + clength, varlist[ss].chrom, varlist[ss].position - 1, varlist[ss].allele1, varlist[ss].allele2, read->sequence);
    }

    a1 = strlen(varlist[ss].allele1);
    a2 = strlen(varlist[ss].allele2);
    for (j = 1; j < a1 && j + offset - 1 < clength; j++) {
        if (varlist[ss].allele1[j] != read->sequence[l1 + j + offset - 1] && read->sequence[l1 + j + offset - 1] != 'N') score1++;
        if (read->sequence[l1 + j + offset - 1] != 'N') b1++;
    }
    for (j = 1; j < a2 && j + offset - 1 < clength; j++) {
        if (varlist[ss].allele2[j] != read->sequence[l1 + j + offset - 1] && read->sequence[l1 + j + offset - 1] != 'N') score2++;
        if (read->sequence[l1 + j + offset - 1] != 'N') b2++;
    }

    // need to limit these loops to reduce computation time for long reads....
    //for (j=offset-1+a1;j<clength;j++) { if (reflist->sequences[reflist->current][start+l2-1+j] != read->sequence[l1+j])score1++; }
    //for (j=offset-1+a2;j<clength;j++) { if (reflist->sequences[reflist->current][start+l2-1+j+a1-a2] != read->sequence[l1+j])score2++; }

    j = offset - 1 + a1;
    while (j < clength || j - a1 + a2 < clength) {
        if (j < clength) {
            if (start + l2 - 1 + j < reflist->lengths[reflist->current] && reflist->sequences[reflist->current][start + l2 - 1 + j] != read->sequence[l1 + j] && read->sequence[l1 + j] != 'N')score1++;
            if (read->sequence[l1 + j] != 'N') b1++;
        }
        if (j - a1 + a2 < clength) {
            if (start + l2 - 1 + j < reflist->lengths[reflist->current] && reflist->sequences[reflist->current][start + l2 - 1 + j] != read->sequence[l1 + j - a1 + a2] && read->sequence[l1 + j - a1 + a2] != 'N') score2++;
            if (read->sequence[l1 + j - a1 + a2] != 'N') b2++;
        }
        if ((score1 > score2 + 1 || score2 > score1 + 1) && b1 >= 10 && b2 >= 10) break;
        j++;
    }
    k = j;
    //fprintf(stdout,"offset %d cigar %d \n",offset,clength);
    if (pflag) {
        if (offset - 1 >= 0) {
            for (j = offset - 1; j < clength; j++) fprintf(stdout, "%c", read->sequence[l1 + j]);
        } else {
            fprintf(stdout, "_");
            for (j = offset; j < clength; j++) fprintf(stdout, "%c", read->sequence[l1 + j]);
        }
        //fprintf(stdout," read\n");
        fprintf(stdout, "%s", varlist[ss].allele1);
        for (j = offset - 1 + a1; j < k; j++) fprintf(stdout, "%c", reflist->sequences[reflist->current][start + l2 - 1 + j] + 32);
        fprintf(stdout, " allele1 %d/%d\t", score1, b1);
        fprintf(stdout, "%s", varlist[ss].allele2);
        for (j = offset - 1 + a2; j < k; j++) fprintf(stdout, "%c", reflist->sequences[reflist->current][start + l2 - 1 + j + a1 - a2] + 32);
        fprintf(stdout, " allele2 %d/%d\t", score2, b2);
    }

    // calculate score3 where we assume that the base-base alignment from right-end of read is correct
    // last base in read l1+read->cigarlist[i]-1 aligns to position start+l2+clength-2 on read
    // varlist[ss].position + length of deletion /length of insertion
    //s1 = varlist[ss].position + a1; s2 = l1+s1-start-l2-1;

    /* score to the left of indel with reference allele */
    s1 = start - 1 + l2 + clength - 1;
    s2 = l1 + clength - 1;
    s1 = varlist[ss].position - 2 + a1 - 1;
    s2 -= start - 1 + l2 + clength - 1;
    s2 += varlist[ss].position - 2 + a1 - 1;

    if (s2 < clength) {
        b3 = 0;
        while (s2 >= 0 && b3 < 10 && s1 >= 0) {
            if (read->sequence[s2] != reflist->sequences[reflist->current][s1] && read->sequence[s2] != 'N') score3++;
            if (read->sequence[s2] != 'N') b3++;
            //fprintf(stdout,"%d-%d %c-%c \n",s2,s1,read->sequence[s2],reflist->sequences[reflist->current][s1]);
            s2--;
            s1--;
        }
        if (pflag) fprintf(stdout, "%d %d scoreref %d| %d\t", s1, s2, score3, b3);
    }
    /* score to the left of indel with reference allele */

    s1 = start - 1 + l2 + clength - 1;
    s2 = l1 + clength - 1;
    s1 = varlist[ss].position - 2 + a1 - 1;
    s2 -= start - 1 + l2 + clength - 1;
    s2 += varlist[ss].position - 2 + a1 - 1;
    s1 -= a1;
    if (s2 < clength) {
        b4 = 0;
        for (j = a2 - 1; j >= 0; j--) {
            if (s2 < 0) {
                partialflag = 1;
                break;
            } // should be flagged so that if the read doesn't cover the full insertion to the left, it is not marked as allele 1
            if (read->sequence[s2] != varlist[ss].allele2[j] && read->sequence[s2] != 'N') score4++;
            if (read->sequence[s2] != 'N')b4++;
            //fprintf(stdout,"%d-%d %c-%c \n",s2,j,read->sequence[s2],varlist[ss].allele2[j]);
            s2--;
        }
        while (s2 >= 0 && b4 < 10 && s1 >= 0) {
            if (read->sequence[s2] != reflist->sequences[reflist->current][s1] && read->sequence[s2] != 'N') score4++;
            //fprintf(stdout,"%d-%d %c-%c \n",s2,s1,read->sequence[s2],reflist->sequences[reflist->current][s1]);
            s2--;
            s1--;
            if (read->sequence[s2] != 'N') b4++;
        }
        //while (s1 >= varlist[ss].position+a1) { fprintf(stdout,"%c-%c ",reflist->sequences[reflist->current][s1-1],read->sequence[s2-1]); s1--; s2--; }
        if (pflag) fprintf(stdout, "%d %d scorealt %d| %d ", s1, s2, score4, b4);
    }

    // if difference between scores ==1, assign low base quality to allele, if both scores > 0, assign LBQ
    // determine the read supports ref-allele, var-allele or ambiguous
    int allele = -1;
    int quality = 20; // base quality
    int op;
    if (score1 < score2) // reference allele matches the read better than variant allele going left -> right
    {
        //if (b3 == 0 && b4 ==0) allele = 0;
        if (score3 < score4) allele = 0; // should this be changed to strictly less than
        else if (score4 == 0 && score3 >= 2) allele = 1;
        // additional filter for long cigars
        op = read->cigarlist[0]&0xf;
        if (allele == 1 && CL != 0 && (CL != 1 || op != BAM_CSOFT_CLIP)) allele = -2;
    } else if (score2 < score1) {
        op = read->cigarlist[read->cigs - 1]&0xf;
        if (score3 < score4 && (score1 - score2 >= 2 || score2 == 0)) allele = 1;
        if (allele == 1 && CL != read->cigs - 1 && (CL != read->cigs - 2 || op != BAM_CSOFT_CLIP)) allele = -2;
        if (score1 - score2 == 1) quality = 13;
        //else if (score2 ==0) quality = 30; else quality = 20;
    } else if (score1 == score2) // equal match score going from left to right
    {
        op = read->cigarlist[read->cigs - 1]&0xf;
        if (score3 <= score4) allele = -1;
        else if (score3 - score4 >= 2 && score4 < 2) allele = 1;
        if (allele == 1 && CL != read->cigs - 1 && (CL != read->cigs - 2 || op != BAM_CSOFT_CLIP)) allele = -2;
        //if (score4 ==0) quality = 30; else quality = 20;
    }

    if (partialflag == 1 && allele == 1) allele = -1;
    if (pflag) fprintf(stdout, "allele %d\n\n", allele);
    if (allele == 0 || allele == 1) {
        fragment->alist[fragment->variants].varid = ss;
        fragment->alist[fragment->variants].qv = (char) (quality + QVoffset);
        //if (read->quality[l1+offset] < fragment->alist[fragment->variants].qv) fragment->alist[fragment->variants].qv = read->quality[l1+offset];
        if (allele == 0) fragment->alist[fragment->variants].allele = '0';
        else if (allele == 1) fragment->alist[fragment->variants].allele = '1';
        varlist[ss].depth++;
        if (allele == 0) {
            if ((read->flag & 16) == 16) varlist[ss].A1 += 1 << 16;
            else varlist[ss].A1 += 1;
        } else if (allele == 1) {
            if ((read->flag & 16) == 16) varlist[ss].A2 += 1 << 16;
            else varlist[ss].A2 += 1;
        }
        fragment->variants++;
    }
    return 0;
}
std::map<int, float> parse_freq(std::string freq_file){
    std::ifstream infile(freq_file);
    std::map<int, float> result;
    int pos;
    float freq;
    while (infile >> pos >> freq) {
        result[pos] = freq;
    }
    return result;
}

int reads_in_sv_region(VARIANT* varlist, int* prev_bnd_pos, alignedread* read){
    char* chrom = read->chrom; int start = read->position, end = read->position + read->span;
    auto v = varlist[*prev_bnd_pos];
    if (((start > v.bnd_pos + 20 && start < v.bnd_mate_pos - 20) || (end  < v.bnd_mate_pos -20 && end > v.bnd_pos + 20)) && strcmp(v.chrom, chrom) == 0 && v.bnd_mate_pos != 0) {
        if (v.bnd_type == BND_DEL) {
            return  1;
        } else if (v.bnd_type == BND_INV) {
            //reads pair cannot all inside the invertion region.
            if (abs(read->IS) >= 300 and abs(read->IS) <= 500 && (std::min(read->position, read->mateposition) < v.bnd_pos - 20 || std::max(read->position, read->mateposition) > v.bnd_mate_pos + 20)) {
                return 1;
            }
        } else if (v.bnd_type == BND_DUP) {
//            for dup, we will process those position according depth.
            return 2;
        }
//        return v.bnd_type;
    }
    return -1;

//    if (((start > v.bnd_pos && start < v.bnd_mate_pos) || (end  < v.bnd_mate_pos && end > v.bnd_pos)) && strcmp(v.chrom, chrom) == 0 && v.bnd_mate_pos != 0) {
//        return v.bnd_type;
//    } else{
//        return -1;
//    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// this function needs some additional checks july 18 2012
// find the variants that are covered by the read and determine the alleles at each of those variants
// TODO: if SV is 1/1, we need to extract the linkage between pair-ended data
// TODO: if SV is 0/1, we only connect the bnd with snp. caution: the direction matters 
// TODO: test for spanning read, i.e. split alignment: 1) we dont need to consider the BND direction for them 2) notice that the bnd is covered through hardclip/softclip. 3) notice that the bnd alignment is not exactly matching
// TODO: test for discordant support reads, i.e. abnormal insertion size. 1) need to consider the direction? 2) need to use the pair-ended info to locate breake-end
int extract_variants_read(struct alignedread* read, HASHTABLE* ht, CHROMVARS* chromvars, VARIANT* varlist, int paired, FRAGMENT* fragment, int chrom, REFLIST* reflist, int * prev_bnd_pos, bool is_found) {
    if (strcmp(read->readid,"chr22_mrecom_10560865_10561181_0:0:0_0:0:0_b71a7") == 0){
        auto tmppp = 3;
    }
    std::set<int> bnd_sses; // reads cover bnds
    int start = read->position;
    int end = start + read->span;
    int ss = 0, firstvar = 0, j = 0, ov = 0, i = 0;     // ss -> variant id; j -> hash index, ov -> no of variants covered by the reads
    j = (int) (start / BSIZE);
    if (j >= chromvars[chrom].blocks) return 0; // another BUG april29 2011 found here

    ss = chromvars[chrom].intervalmap[j];   // if reads out of hashing range, return, else find hashed variant id

    if (ss < 0 || ss >= VARIANTS) return 0;

    // check if ss is less than first variant for the chromosome 'chrom', if so assign it to the first variant
    if (ss < chromvars[chrom].first) ss = chromvars[chrom].first;

    if (varlist[ss].position <= end) {
        while (ss < VARIANTS && varlist[ss].position < start && ss <= chromvars[chrom].last) ss++;  //locate the first variant on reads
        firstvar = ss;
        while (ss < VARIANTS && varlist[ss].position <= end && ss <= chromvars[chrom].last) {       //locate the last variant on reads
            //printf("variant %s %d %c %c \n",varlist[ss].chrom,varlist[ss].position,varlist[ss].allele1,varlist[ss].allele2);
            ov++;
            ss++;
        }
    }

    //fprintf(stderr,"chrom %d variants %d ss %d first %d-%d span %d-%d\n",chrom,VARIANTS,ss,chromvars[chrom].first,chromvars[chrom].last,start,end);
    //fprintf(stderr,"ov %d %d\n",ov,firstvar);
    if ((paired == 0 && ov < 2 && SINGLEREADS == 0) || (paired == 0 && ov < 1 && SINGLEREADS == 1) || (paired == 1 && ov < 1)) return 0; //no variant on reads, return
    ss = firstvar; // use variable firstvar to store first variant that overlaps this read

    int l1 = 0, l2 = 0; // l1 is advance on read, l2 is advance on reference genome
    int op = 0, ol = 0; // op == operation; ol == operation length
    bool support_ref_bnd_reads = false;
//    if pos in bnd region
    auto region_tag = reads_in_sv_region(varlist, prev_bnd_pos, read);
    int prev_dup_pos = 0;
    if((region_tag == 1 and !is_found)) {
        support_ref_bnd_reads = true;
        bnd_sses.insert(*prev_bnd_pos);
    } else if (region_tag == 2) prev_dup_pos = *prev_bnd_pos;
    ss = firstvar;
    for (i = 0; i < read->cigs; i++) {          //iter through base with CIGAR
        //fprintf(stdout,"%c %d \t",(char)read->cigarlist[i+1],read->cigarlist[i]);
        while (varlist[ss].position < start + l2 && ss <= chromvars[chrom].last) ss++;
        op = read->cigarlist[i]&0xf;
        ol = read->cigarlist[i] >> 4;
//        if(ss == 281850) {
//            int m = 000;
//        }

        if (op == BAM_CMATCH) {
            while (ss <= chromvars[chrom].last && varlist[ss].position >= start + l2 && varlist[ss].position < start + l2 + ol) {
                if(varlist[ss].bnd == 1 && reads_in_sv_region(varlist, &ss, read) != -1 && varlist[ss].heterozygous == '1') {
                    bnd_sses.insert(ss);
                    support_ref_bnd_reads = true;
                    *prev_bnd_pos = ss;
                }
                // function call
                if (varlist[ss].heterozygous == '1' && varlist[ss].type == 0 && varlist[ss].bnd == 0)
                    compare_read_SNP(read, varlist, ss, start, l1, l2, fragment);
                else if (varlist[ss].heterozygous == '2' && varlist[ss].type == 0 && varlist[ss].bnd == 0) {
                    compare_read_SNP(read, varlist, ss, start, l1, l2, fragment);
                } else if (varlist[ss].heterozygous == '1' && varlist[ss].type != 0 && varlist[ss].position < start + l2 + ol - 1 && reflist->current >= 0 && PARSEINDELS == 1 && varlist[ss].bnd == 0) {
                    compare_read_INDEL(read, varlist, ss, start, l1, l2, ol, fragment, i, reflist);
                }
                ss++;
            }

            l1 += ol;
            l2 += ol;
        } else if (op == BAM_CINS) {
            fragment->is_all_m = false;
            support_ref_bnd_reads = false;
            if (varlist[ss].heterozygous == '1' && varlist[ss].position == start + l2 && varlist[ss].type == ol && PARSEINDELS == 1 && ss <= chromvars[chrom].last && varlist[ss].bnd == 0) {
                if (IFLAG) fprintf(stdout, "%s INSERTION %d %s:%d:%s:%s\n", read->readid, start + l2, varlist[ss].chrom, varlist[ss].position, varlist[ss].RA, varlist[ss].AA);
                fragment->alist[fragment->variants].varid = ss;
                fragment->alist[fragment->variants].allele = '1';
                fragment->alist[fragment->variants].qv = read->quality[l1];
                fragment->variants++;
                varlist[ss].depth++;
                if ((read->flag & 16) == 16) varlist[ss].A2 += 1 << 16;
                else varlist[ss].A2 += 1;
                ss++;
            } else if(varlist[ss].bnd == 1) {
                *prev_bnd_pos = ss;
                ss++;
            }
            l1 += ol;
        } else if (op == BAM_CDEL) {
            fragment->is_all_m = false;
            support_ref_bnd_reads = false;
            if (varlist[ss].heterozygous == '1' && varlist[ss].position == start + l2 && varlist[ss].type == -1 * ol  && PARSEINDELS == 1 && ss <= chromvars[chrom].last && varlist[ss].bnd == 0) {
                if (IFLAG) fprintf(stdout, "%s DELETION %d %s:%d:%s:%s\n", read->readid, start + l2, varlist[ss].chrom, varlist[ss].position, varlist[ss].RA, varlist[ss].AA);
                fragment->alist[fragment->variants].varid = ss;
                fragment->alist[fragment->variants].allele = '1';
                fragment->alist[fragment->variants].qv = read->quality[l1];
                varlist[ss].depth++;
                if ((read->flag & 16) == 16) varlist[ss].A2 += 1 << 16;
                else varlist[ss].A2 += 1;
                fragment->variants++;
                ss++;
            } else if(varlist[ss].bnd == 1) {
//                support_ref_bnd_reads = true;
                *prev_bnd_pos = ss;
                ss++;
            }
            l2 += ol;
        } else if (op == BAM_CREF_SKIP){
            fragment->is_all_m = false;
            support_ref_bnd_reads = false;
            l2 += ol;
        }
        else if ( op == BAM_CSOFT_CLIP) // primary alignment
        {
            fragment->is_all_m = false;
            support_ref_bnd_reads = false;

            //fixme sv position is not a precise location, a range +-5,not used , reads provided
            if (PARSEBND && SUPPORT_READS_TAG == nullptr)
            {
                if (varlist[ss].heterozygous == '1' && (varlist[ss].position <= start + l2 + BND_RANGE + 1 && varlist[ss].position >= start + l2 - BND_RANGE - 1) && varlist[ss].bnd == 1 && ss <= chromvars[chrom].last && ol >= 10 )
                {
//                    fragment->alist[fragment->variants].varid = ss;
//                    fragment->alist[fragment->variants].allele = '1';
//                    fragment->alist[fragment->variants].qv = read->quality[l1];
//                    fragment->variants++;
//                    varlist[ss].depth++;
//                    if ((read->flag & 16) == 16) varlist[ss].A2 += 1 << 16;
//                    else varlist[ss].A2 += 1;
                    ss++;
                }
            } else {
                if (varlist[ss].heterozygous == '1' && varlist[ss].position == start + l2 && varlist[ss].bnd == 1 && ss <= chromvars[chrom].last)
                {
//                    fragment->alist[fragment->variants].varid = ss;
//                    fragment->alist[fragment->variants].allele = '1';
//                    fragment->alist[fragment->variants].qv = read->quality[l1];
//                    fragment->variants++;
//                    varlist[ss].depth++;
//                    if ((read->flag & 16) == 16) varlist[ss].A2 += 1 << 16;
//                    else varlist[ss].A2 += 1;
                    ss++;
                }
            }
            l1 += ol;
        }
        else if (op == BAM_CHARD_CLIP) // secondary alignment ,notice that this is quite uncommon for reads with leght less than 150bp
        {
            fragment->is_all_m = false;
            support_ref_bnd_reads = false;

            if (PARSEBND && SUPPORT_READS_TAG == nullptr)
            {
                if (varlist[ss].heterozygous == '1' && (varlist[ss].position <= start + l2 + BND_RANGE + 1 && varlist[ss].position >= start + l2 - BND_RANGE - 1) && varlist[ss].bnd == 1 && ss <= chromvars[chrom].last && ol >= 10 )
                {
//                    fragment->alist[fragment->variants].varid = ss;
//                    fragment->alist[fragment->variants].allele = '1';
//                    fragment->alist[fragment->variants].qv = read->quality[l1];
//                    fragment->variants++;
//                    varlist[ss].depth++;
//                    if ((read->flag & 16) == 16) varlist[ss].A2 += 1 << 16;
//                    else varlist[ss].A2 += 1;
                    ss++;
                }
            } else {
                if (varlist[ss].heterozygous == '1' && varlist[ss].position == start + l2 && varlist[ss].bnd == 1 && ss <= chromvars[chrom].last)
                {
//                    fragment->alist[fragment->variants].varid = ss;
//                    fragment->alist[fragment->variants].allele = '1';
//                    fragment->alist[fragment->variants].qv = read->quality[l1];
//                    fragment->variants++;
//                    varlist[ss].depth++;
//                    if ((read->flag & 16) == 16) varlist[ss].A2 += 1 << 16;
//                    else varlist[ss].A2 += 1;
                    ss++;
                }
            }
            l2 += 0;
        }
    }

if (SUPPORT_READS_TAG == nullptr && PARSEBND && DATA_TYPE == 1)
{//this works for het SV BND only

    // start - 300
    //TODO: add direction check
    if (read->IS > MAX_IS) {
        support_ref_bnd_reads = false;

        int start_pos = start - 300;
            int sss = firstvar;
            while (varlist[sss].position >= start_pos && sss >= chromvars[chrom].first) {
                if (varlist[sss].bnd ==1 && varlist[sss].heterozygous == '1' ) {
                    if (varlist[sss].bnd_type != BNDTYPE_PAIRED || varlist[sss].bnd_type != BND_INS) //we only consider paired at this moment
                        break;
                    if (abs(varlist[sss].bnd_pair_distance - read->IS) < 1500)
                    {
                        fragment->alist[fragment->variants].varid = sss;
                        fragment->alist[fragment->variants].allele = '1';
                        fragment->alist[fragment->variants ].qv = read->quality[l1];
                        fragment->variants++;
                        varlist[sss].depth++;
                        if ((read->flag & 16) == 16) varlist[sss].A2 += 1 << 16;
                        else varlist[sss].A2 += 1;
                        break;
                    }
                    else {
                        break;
                    }
                }
                sss--;
            }
        }
    //TODO: add direction check
    //end + 300
        if (read->IS > MAX_IS) {
            support_ref_bnd_reads = false;

            int end_pos = end + 300; // positive abnormal IS, search sv bnd after reads alignment end 300 bp
            while (varlist[ss].position <= end_pos && ss <= chromvars[chrom].last) {
                if (varlist[ss].bnd ==1 && varlist[ss].heterozygous == '1' ) {
                    if (varlist[ss].bnd_type != BNDTYPE_PAIRED) //we only consider paired at this moment
                        break;
                    if (abs(varlist[ss].bnd_pair_distance - read->IS) < 1500)
                    {
                        fragment->alist[fragment->variants].varid = ss;
                        fragment->alist[fragment->variants].allele = '1';
                        fragment->alist[fragment->variants ].qv = read->quality[l1];
                        fragment->variants++;
                        varlist[ss].depth++;
                        if ((read->flag & 16) == 16) varlist[ss].A2 += 1 << 16;
                        else varlist[ss].A2 += 1;
                        break;
                    }
                    else {
                        break;
                    }
                }
                ss++;
            }
        }
    }
//if
if (300 > abs(read->IS) || abs(read->IS) > 500) {
    support_ref_bnd_reads = false;
}

if (region_tag == 2) {
    for (int k = 0; k < fragment->variants; k++) {
        if (fragment->alist[k].allele == '0') {
            if (varlist[prev_dup_pos].snp0_dup_region->find(fragment->alist[k].varid) != varlist[prev_dup_pos].snp0_dup_region->end()) {
                (*(varlist[prev_dup_pos].snp0_dup_region))[fragment->alist[k].varid] = (*(varlist[prev_dup_pos].snp0_dup_region))[fragment->alist[k].varid] + 1;
            } else {
                (*(varlist[prev_dup_pos].snp0_dup_region))[fragment->alist[k].varid] = 1;
                (*(varlist[prev_dup_pos].snp1_dup_region))[fragment->alist[k].varid] = 0;
            }
        } else {
            if (varlist[prev_dup_pos].snp1_dup_region->find(fragment->alist[k].varid) != varlist[prev_dup_pos].snp1_dup_region->end()) {
                (*(varlist[prev_dup_pos].snp1_dup_region))[fragment->alist[k].varid] = (*(varlist[prev_dup_pos].snp1_dup_region))[fragment->alist[k].varid] + 1;
            } else {

                (*(varlist[prev_dup_pos].snp1_dup_region))[fragment->alist[k].varid] = 1;
                (*(varlist[prev_dup_pos].snp0_dup_region))[fragment->alist[k].varid] = 0;
            }
        }
    }
}

if (PARSEBND && support_ref_bnd_reads) {
    for (auto cbnd_ss : bnd_sses) {
        fragment->alist[fragment->variants].varid = cbnd_ss;
        fragment->alist[fragment->variants].allele = '0';
        fragment->alist[fragment->variants ].qv = 60;
        fragment->alist[fragment->variants].is_bnd = true;
        fragment->variants++;
        varlist[cbnd_ss].depth++;
        if ((read->flag & 16) == 16) varlist[cbnd_ss].A2 += 1 << 16;
        else varlist[cbnd_ss].A2 += 1;
    }
}
//if(PARSEBND && SUPPORT_READS_TAG == nullptr) {
//    if(support_ref_bnd_reads) {
//        fragment->alist[fragment->variants].varid = prev_bnd_pos;
//        fragment->alist[fragment->variants].allele = '0';
//        fragment->alist[fragment->variants ].qv = read->quality[l1];
//        fragment->variants++;
//        varlist[prev_bnd_pos].depth++;
//        if ((read->flag & 16) == 16) varlist[prev_bnd_pos].A2 += 1 << 16;
//        else varlist[prev_bnd_pos].A2 += 1;
//    }
//}
     return 0;
    //	printf("read %s %d %s %d %d %d %s XM %d parsed %d %d %d vars %d\n",read->readid,read->flag,read->chrom,read->position,read->mquality,read->IS,read->cigar,read->XM,chrom,start,end,ov);
}

//int copy_fragment(FRAGMENT* fnew, FRAGMENT* fragment, struct alignedread* read) {
//    int i = 0, sl;
//    fnew->variants = fragment->variants;
//    fnew->paired = 1;
//    fnew->alist = (allele*) malloc(sizeof (allele) * fragment->variants);
//    fnew->absIS   = fragment->absIS;
//    for (i = 0; i < fragment->variants; i++) {
//        fnew->alist[i].varid = fragment->alist[i].varid;
//        fnew->alist[i].allele = fragment->alist[i].allele;
//        fnew->alist[i].qv = fragment->alist[i].qv;
//    }
//}


// add a fragment to the flist whose mate is yet to be seen

int add_fragment(FRAGMENT* flist, FRAGMENT* fragment, struct alignedread* read, int fragments) {
    int i = 0, sl = 0;
    flist[fragments].variants = fragment->variants;
    flist[fragments].paired = 1;
    flist[fragments].absIS = fragment->absIS;
    flist[fragments].read_qual = fragment->read_qual;
    flist[fragments].bnd_reads = fragment->bnd_reads;
    flist[fragments].is_all_m = fragment->is_all_m;

    // april 10 2012 change made for Poplar data to not use matepos...
    if (read->IS > 0) flist[fragments].matepos = read->mateposition;
        //if (read->IS > 0) flist[fragments].matepos = read->position + read->IS;
    else flist[fragments].matepos = read->position;

    flist[fragments].alist = (allele*) malloc(sizeof (allele) * fragment->variants);
    for (i = 0; i < fragment->variants; i++) {
        flist[fragments].alist[i].varid = fragment->alist[i].varid;
        flist[fragments].alist[i].allele = fragment->alist[i].allele;
        flist[fragments].alist[i].qv = fragment->alist[i].qv;
        flist[fragments].alist[i].is_bnd = fragment->alist[i].is_bnd;
    }
    sl = strlen(read->readid);
    flist[fragments].id = (char*) malloc(sl + 1);
    for (i = 0; i < sl; i++) flist[fragments].id[i] = read->readid[i];
    flist[fragments].id[i] = '\0';
    if (read->barcode == NULL){
        flist[fragments].barcode = NULL;
    }else{
        sl = strlen(read->barcode);
        flist[fragments].barcode = (char*) malloc(sl + 1);
        for (i = 0; i < sl; i++) flist[fragments].barcode[i] = read->barcode[i];
        flist[fragments].barcode[i] = '\0';
    }
    flist[fragments].rescued = read->rescued;
    flist[fragments].dm = read->dm;

    return 0;
}
