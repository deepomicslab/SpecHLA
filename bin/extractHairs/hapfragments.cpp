#include "hapfragments.h"

int MAX_BQ=93; // base quality cannot exceed this, QV=60

// when --vcf-phased specified, filter fragment by phasing info

int filter_by_phasing_info(FRAGMENT* fragment, VARIANT* varlist){
    int len = fragment->variants;
    int phase_set_count = 0;
    //h0 > left side count //h1 > right side count // fragment score for each phase_set // side: vcf gt x/x, left 0, right 1
    HASHTABLE h0, h1;
    h0.htsize = len;
    h1.htsize = len;
    init_hashtable(&h0);
    init_hashtable(&h1);
//    int* score = (int*) malloc(sizeof(int) * len);
    char** phase_set_keys = (char**) malloc(sizeof(char*) * len);
    int phase_set = 0;
    VARIANT* var = NULL;
    for (int i = 0; i < len; i++) {
        var = &varlist[fragment->alist[i].varid];
        phase_set = var->phase_set;
        if(phase_set == 0) continue;
//        phase_set to phase_set_str
        int length = snprintf( NULL, 0, "%d", phase_set );
        char* phase_set_str = (char *)malloc( sizeof(char)*(length + 1));
        snprintf( phase_set_str, length + 1, "%d", phase_set );
//        specific side left or right
        int v = -1;
        if (var->genotype[0] == fragment->alist[i].allele) {
            v = getindex(&h0, phase_set_str);
            if (v == -1) insert_keyvalue(&h0, phase_set_str, length, 1);
            else update_value(&h0, phase_set_str, v+1);
        }else {
            v = getindex(&h1, phase_set_str);
            if (v == -1) insert_keyvalue(&h1, phase_set_str, length, 1);
            else update_value(&h1, phase_set_str, v+1);
        }
        int found = 0;
        for (int k = 0; k < phase_set_count; k++) {
            if(strcmp(phase_set_keys[k], phase_set_str) == 0) {found = 1;break;};
        }
        if (found == 0) {
            phase_set_keys[phase_set_count] = (char *)malloc( sizeof(char)*(length + 1));
            phase_set_keys[phase_set_count][length] = '\0';
            strcpy(phase_set_keys[phase_set_count++], phase_set_str);
        }
        free(phase_set_str);
    }
    if (phase_set_count == 0) return 1;
// score
    char* p_s_k;
    float c_0, c_1, score, score_gt_65 = 0;
    for (int j = 0; j < phase_set_count; j++) {
        p_s_k = phase_set_keys[j];
        if (p_s_k == NULL) continue;
        c_0 = getindex(&h0, p_s_k);
        c_1 = getindex(&h1, p_s_k);
        score = c_0 / (c_1 + c_0);
        if (score >= 0.65) score_gt_65++;
        free(phase_set_keys[j]);
    }
    free(phase_set_keys);
    if(score_gt_65/phase_set_count > 0.65) return 1;
    return 0;
}

// sort such that mate pairs are together and reverse sorted by starting position of second read in a mate-piar

int compare_fragments(const void *a, const void *b) {
    FRAGMENT* f1 = (FRAGMENT*) a;
    FRAGMENT* f2 = (FRAGMENT*) b;
    if (f1->matepos == f2->matepos) return strcmp(f1->id, f2->id);
    else return f2->matepos - f1->matepos;
}

int compare_alleles(const void *a, const void *b) {
    if (((allele*) a)->varid == ((allele*) b)->varid) return ((allele*) a)->allele - ((allele*) b)->allele;
    else return ((allele*) a)->varid - ((allele*) b)->varid;
}

int filter_ref_bnd(FRAGMENT* fragment) {
    int count = fragment->variants;
    allele* tmp_alleles = (allele*) malloc(sizeof (allele)*fragment->variants+10);
    for (int i = 0; i < fragment->variants;i++) {
        tmp_alleles[i].varid = fragment->alist[i].varid;
        tmp_alleles[i].allele= fragment->alist[i].allele;
        tmp_alleles[i].qv = fragment->alist[i].qv;
        tmp_alleles[i].is_bnd = fragment->alist[i].is_bnd;
    }
    fragment->variants = 0;
    int j = 0;
    int k  = 0;
    for (volatile int k = 0; k < count;k++) {
        if (tmp_alleles[k].is_bnd) {
            if (!fragment->is_all_m) {
                continue;
            }
        }
        fragment->alist[j].varid = tmp_alleles[k].varid;
        fragment->alist[j].allele = tmp_alleles[k].allele;
        fragment->alist[j].qv = tmp_alleles[k].qv;
        fragment->alist[j].is_bnd = tmp_alleles[k].is_bnd;
        fragment->variants++;
        j++;
    }
}

int print_fragment(FRAGMENT* fragment, VARIANT* varlist, FILE* outfile)  {
    if (strcmp(fragment->id, "SRR5114981.518061_MP") == 0){
        auto tmpp = 33;
    }
    if (SUPPORT_READS_TAG == NULL) {
        filter_ref_bnd(fragment);
//        if (fragment->support_reads < SUPPORT_READS) return 0;
    }

    if (fragment->variants < 2 && DATA_TYPE != 2) return 0;
    if (strcmp(fragment->id, "D00360:95:H2YWMBCXX:2:2214:11806:83148") == 0) {
        int tmp = 33;
    }
    if (PRINT_FRAGMENTS == 0) return 0;

    if (VCF_PHASED) {
        int is_print = filter_by_phasing_info(fragment, varlist);
        if (is_print == 0) return 0;
    }
    int i = 0;
    sort_framgment(fragment);
    /*
       fprintf(stdout,"HAIR %s %d \t",fragment->id,fragment->variants);
       j = fragment->alist[i].varid;
       for (i=0;i<fragment->variants;i++) fprintf(stdout,"%d %s %d %c/%c %c %c \t",j,varlist[j].chrom,varlist[j].position,varlist[j].allele1,varlist[j].allele2,fragment->alist[i].allele,fragment->alist[i].qv);
       fprintf(stdout,"\n");
     */
    // varid is the index of the variant in the list of variants (this should actually  be the index in the VCF file)
    // fragment is printed using 1-based coordinate system instead of 0-based since this is encoded in HapCUT

    fragment->blocks = 1;

//    if bnd and not print,


    for (i = 0; i < fragment->variants - 1; i++) {
        if (fragment->alist[i + 1].varid - fragment->alist[i].varid != 1) fragment->blocks++;
    }
    fprintf(outfile, "%d %s", fragment->blocks,fragment->id);
    //fprintf(outfile, "%d %c:%s", fragment->blocks, fragment->strand,fragment->id);

    //new format prints col 3 as data type (0 for normal, 1 for HiC) and col 4 as mate 2 index
    if (DATA_TYPE == 2){
        if(fragment->barcode != NULL)
//            if(fragment->barcode[0] == '\0')
            fprintf(outfile, " 2 %s -1", fragment->barcode);
        else
            fprintf(outfile, " 2 NULL -1");
    }else if (NEW_FORMAT)
        fprintf(outfile, " %d -1 -1", DATA_TYPE);

    //for (i=0;i<fragment->variants;i++) fprintf(stdout,"%c",fragment->alist[i].qv);
    // varid is printed with offset of 1 rather than 0 since that is encoded in the Hapcut program
    if (PRINT_COMPACT ==1)
    {

        fprintf(outfile, " %d %c", fragment->alist[0].varid + 1, fragment->alist[0].allele);
        for (i = 1; i < fragment->variants; i++) {
            if (fragment->alist[i].varid - fragment->alist[i - 1].varid == 1) fprintf(outfile, "%c", fragment->alist[i].allele);
            else fprintf(outfile, " %d %c", fragment->alist[i].varid + 1, fragment->alist[i].allele);
        }
        fprintf(outfile, " ");
        for (i = 0; i < fragment->variants; i++) {
            fprintf(outfile, "%c", fragment->alist[i].qv);
        }
        fprintf(outfile, " ");
        fprintf(outfile, "%d", fragment->read_qual);
        if (DATA_TYPE == 2)
        {
            if (fragment->rescued != 0 && fragment->rescued != 1)
                fragment->rescued = 0;
            fprintf(outfile, " %d %f", fragment->rescued, fragment->dm);
        }
    }
    else // individual variant format, for debugging
    {
        //fprintf(outfile,":%c",fragment->strand);
        for (i = 0; i < fragment->variants; i++) fprintf(outfile, " %d:%c:%d",fragment->alist[i].varid+1,fragment->alist[i].allele,(char)fragment->alist[i].qv-33);
    }
    if (fragment->bnd_reads == 1) {
        fprintf(outfile, " SV");
    } else if (fragment->bnd_reads == 2) {
        fprintf(outfile, " SV_REF");
    }
    fprintf(outfile, "\n");

    return 0;
}

// make sure they are in the correct order, i+1 could be < i

int print_matepair(FRAGMENT* f1, FRAGMENT* f2, VARIANT* varlist, FILE* outfile) {
//    sort_framgment(f1);
//    sort_framgment(f2);
    if (PRINT_FRAGMENTS == 0) return 0;
    int i = 0, j = 0;
    int f2_size = f2->variants;
    for (i = 0; i < f2->variants; i++) {
        bool ee = false;
        for (j = 0; j < f1->variants; j++) {
            if (f1->alist[j].varid == f2->alist[i].varid) {
                ee = true;
                break;
            }
        }
        if (ee)
            f2_size--;
    }
//    if (VCF_PHASED) {
    FRAGMENT* f = (FRAGMENT*)malloc(sizeof(FRAGMENT));
    f->id = (char*) malloc(strlen(f1->id) + 4);
    if (DATA_TYPE == 2 && f1->barcode != nullptr) {
        f->barcode = (char*) malloc(strlen(f1->barcode) + 1);
        strcpy(f->barcode, f1->barcode);
    } else {
        f->barcode = nullptr;
    }
    f->read_qual = (f1->read_qual + f2->read_qual) / 2;
    strcpy(f->id, f1->id);
    strcat(f->id,"_MP");
    if (strcmp(f1->id, "D00360:94:H2YT5BCXX:1:1109:5815:14301") == 0) {
        int tmp = 33;
    }
    f->alist = (allele*) malloc(sizeof (allele) * (f1->variants + f2_size + 1));
    f->variants = f1->variants + f2_size;
    for (i = 0; i < f1->variants; i++) {
        f->alist[i] = f1->alist[i];
    }
    int pos = 0 ;
    for (i = 0; i < f2->variants; i++) {
        bool ee = false;
        for (j = 0; j < f1->variants; j++) {
            if (f1->alist[j].varid == f2->alist[i].varid) {
                ee = true;
                break;
            }
        }
        if (!ee) {
            f->alist[pos+f1->variants] = f2->alist[i];
            pos++;
        }
    }
    f->bnd_reads = f1->bnd_reads;
    f->is_all_m = f1->is_all_m && f2->is_all_m && f1->bnd_reads == 0 && f2->bnd_reads == 0;
    print_fragment(f, varlist,outfile);
    free(f);
//        int is_print = filter_by_phasing_info(f, varlist);
//        free(f);
//        if (is_print == 0) return 0;
//    }
//    f1->blocks = 1;
//    for (i = 0; i < f1->variants - 1; i++) {
//        if (f1->alist[i + 1].varid - f1->alist[i].varid != 1) f1->blocks++;
//    }
//    if (f2->alist[0].varid - f1->alist[f1->variants - 1].varid != 1) f1->blocks++;
//    for (i = 0; i < f2->variants - 1; i++) {
//        if (f2->alist[i + 1].varid - f2->alist[i].varid != 1) f1->blocks++;
//    }
//
//    fprintf(outfile, "%d %s_MP", f1->blocks, f1->id);
//    //fprintf(outfile,"%d %s_%s_MP",f1->blocks,varlist[f1->alist[0].varid].chrom,f1->id);
//    //	for (i=0;i<f1->variants;i++) fprintf(outfile,"%c",f1->alist[i].qv);
//    //	for (i=0;i<f2->variants;i++) fprintf(outfile,"%c",f2->alist[i].qv);
//
//    //new format prints col 3 as data type (0 for normal, 1 for HiC) and col 4 as mate 2 index
//    if (DATA_TYPE == 2){
//        if(f1->barcode != NULL)
//            fprintf(outfile, " 2 %s -1", f1->barcode);
//        else
//            fprintf(outfile, " 2 NULL -1");
//    }
//    else if (NEW_FORMAT)
//        fprintf(outfile, " %d %d %d", DATA_TYPE, f2->alist[0].varid+1, f1->absIS);
////    std::vector<int> f1_pos;
//    // varid is printed with offset of 1 rather than 0 since that is encoded in the Hapcut program
//    fprintf(outfile, " %d %c", f1->alist[0].varid + 1, f1->alist[0].allele);
//    for (i = 1; i < f1->variants; i++) {
//        if (f1->alist[i].varid - f1->alist[i - 1].varid == 1) fprintf(outfile, "%c", f1->alist[i].allele);
//        else fprintf(outfile, " %d %c", f1->alist[i].varid + 1, f1->alist[i].allele);
//    }
//    if (f2->alist[0].varid - f1->alist[f1->variants - 1].varid == 1) fprintf(outfile, "%c", f2->alist[0].allele);
//    else fprintf(outfile, " %d %c", f2->alist[0].varid + 1, f2->alist[0].allele);
//    for (i = 1; i < f2->variants; i++) {
//        if (f2->alist[i].varid - f2->alist[i - 1].varid == 1) fprintf(outfile, "%c", f2->alist[i].allele);
//        else fprintf(outfile, " %d %c", f2->alist[i].varid + 1, f2->alist[i].allele);
//    }
//    fprintf(outfile, " ");
//    for (i = 0; i < f1->variants; i++) fprintf(outfile, "%c", f1->alist[i].qv);
//    for (i = 0; i < f2->variants; i++) fprintf(outfile, "%c", f2->alist[i].qv);
//	fprintf(outfile, " ");
//	fprintf(outfile, "%d", f1->read_qual);
//    if (DATA_TYPE == 2)
//    {
//        if (f1->rescued != 0 && f1->rescued != 1)
//            f1->rescued = 0;
//        fprintf(outfile, " %d %f", f1->rescued, f1->dm);
//    }
//    fprintf(outfile, "\n");

    //	fprintf(outfile," type:");
    //	for (i=0;i<f1->variants;i++) fprintf(outfile,"%d:%d,",varlist[f1->alist[i].varid].type,varlist[f1->alist[i].varid].position);
    //	for (i=0;i<f2->variants;i++) fprintf(outfile,"%d:%d,",varlist[f2->alist[i].varid].type,varlist[f2->alist[i].varid].position);
    /*
       fprintf(outfile,"mated frag %s matepos %d vars %d \t",f1->id,f1->matepos,f1->variants);
       for (j=0;j<f1->variants;j++)
       {
       k = f1->alist[j].varid;
       fprintf(outfile,"%d %s %d %c/%c %c %c \t",k,varlist[k].chrom,varlist[k].position,varlist[k].allele1,varlist[k].allele2,f1->alist[j].allele,f1->alist[j].qv);
       }
       fprintf(outfile,"\n");
       fprintf(outfile,"mated frag %s matepos %d vars %d \t",f2->id,f2->matepos,f2->variants);
       for (j=0;j<f2->variants;j++)
       {
       k = f2->alist[j].varid;
       fprintf(stdout,"%d %s %d %c/%c %c %c \t",k,varlist[k].chrom,varlist[k].position,varlist[k].allele1,varlist[k].allele2,f2->alist[j].allele,f2->alist[j].qv);
       }
       fprintf(stdout,"\n");
     */
    return 0;
}

int print_mate_bnd_fragment(std::unordered_map<std::string , std::pair<int, int>> & BNDs, FILE* outfile){
    for (auto bnd: BNDs) {
        if (bnd.second.first == 0 || bnd.second.second == 0)
            continue;
//        1 SRR3760936.8700_MP 16 000 77A 60
// 2 chr14_2_106757572_106757738_0_1_0_0_0:0:0_0:0:0_a30_MP 2 AAGCTACTCGAAACTC-1 -1 145078 0 145081 0 ID 60 0 0.317771
        if (NEW_FORMAT) {
            for (int i = 0; i < 3; i++){
                fprintf(outfile, "2 BND 2 NNNNNNNNNNNN-1 -1 %d 1 %d 1 ID 60 0 0.31",bnd.second.first, bnd.second.second);
                fprintf(outfile,"\n");
            }
        } else {
            for (int i = 0; i < 3; i++){
                fprintf(outfile, "2 BND %d 1 %d 1 77A 60",bnd.second.first, bnd.second.second);
                fprintf(outfile,"\n");
            }
        }
    }
    return 1;
}

int print_dup_region_snp(VARIANT* variant, FILE* outfile, int idx){
    auto v = variant[idx];
    if (v.bnd_type != BND_DUP)return 0;
    int diff = 0;
    int minv = 0;
    int snp1 = 0;
    if (v.snp0_dup_region->empty() || v.snp0_dup_region == nullptr)
        return 0;
    for(auto item : *(v.snp0_dup_region)) {
        snp1 = (*(v.snp1_dup_region))[item.first];
        diff = abs(item.second - (*(v.snp1_dup_region))[item.first]);
        minv = std::min(item.second,snp1);
        if (diff >= ((minv * 3) / 4)) {
            if (item.second > snp1) {
                if (NEW_FORMAT) {
                    for (int i = 0; i < 3; i++){
                        fprintf(outfile, "2 BND 2 NNNNNNNNNNNN-1 -1 %d 1 %d 0 ID 60 0 0.31",idx + 1, item.first + 1);
                        fprintf(outfile,"\n");
                    }
                } else {
                    for (int i = 0; i < 3; i++){
                        fprintf(outfile, "2 BND %d 1 %d 0 77A 60",idx + 1, item.first + 1);
                        fprintf(outfile,"\n");
                    }
                }
            } else {
                if (NEW_FORMAT) {
                    for (int i = 0; i < 3; i++){
                        fprintf(outfile, "2 BND 2 NNNNNNNNNNNN-1 -1 %d 1 %d 1 ID 60 0 0.31",idx + 1, item.first + 1);
                        fprintf(outfile,"\n");
                    }
                } else {
                    for (int i = 0; i < 3; i++){
                        fprintf(outfile, "2 BND %d 1 %d 1 77A 60",idx + 1, item.first + 1);
                        fprintf(outfile,"\n");
                    }
                }
            }
        }
    }
}

// sort the fragment list by 'mate-position or position of 2nd read' so that reads that are from the same DNA fragment are together
// also takes care of overlapping paired-end reads to avoid duplicates in fragments
void clean_fragmentlist(FRAGMENT* flist, int* fragments, VARIANT* varlist, int currchrom, int currpos, int prevchrom) {
    int i = 0, j = 0, k = 0, first = 0, sl = 0, bl = 0;
    FRAGMENT fragment;
    fragment.variants = 0;
    fragment.alist = (allele*) malloc(sizeof (allele)*16184);
    if (*fragments > 1) qsort(flist, *fragments, sizeof (FRAGMENT), compare_fragments);
    // sort such that mate pairs are together and reverse sorted by starting position of second read in a mate-piar
    //for (i=0;i<*fragments;i++) fprintf(stdout,"frag %s %d vars %d \n",flist[i].id,flist[i].alist[0].varid,flist[i].variants);
    if (currchrom == prevchrom) // need to ignore the top of the fragment list
    {
        first = 0;
        while (flist[first].matepos >= currpos && first < *fragments) first++;
    }
    //	fprintf(stdout,"cleaning the fragment list: current chrom %d %d first %d fragments %d\n",currchrom,currpos,first,*fragments);

    if (*fragments > 1) // bug fixed jan 13 2012, when there is only one fragment, we don't need to check if it is part of mate-pair
    {
        // serious bug fixed here: mate-pairs being examined twice April 5 2012
        // check this code for corrrectness: mate-pairs will be adjacent to each other.
        i = first;
        while (i < (*fragments) - 1) {
            if (strcmp(flist[i].id, flist[i + 1].id) == 0) // mate pair with both ends having at least one variant
            {
                fragment.bnd_reads = flist[i].bnd_reads;
                fragment.is_all_m = flist[i].is_all_m;
                //fprintf(stdout,"mate-pair %s %s %s\n",flist[i].id);
                if (flist[i].alist[flist[i].variants - 1].varid <= flist[i + 1].alist[0].varid) print_matepair(&flist[i], &flist[i + 1], varlist, fragment_file);
                else if (flist[i + 1].alist[flist[i + 1].variants - 1].varid < flist[i].alist[0].varid) print_matepair(&flist[i + 1], &flist[i], varlist, fragment_file);
                else if (flist[i].variants + flist[i + 1].variants >= 2) {
                    j = 0;
                    k = 0;
                    for (int t = 0; t < fragment.variants; t++) {
                        fragment.alist[t].is_bnd = false;
                    }
                    fragment.variants = 0;
                    while (j < flist[i].variants || k < flist[i + 1].variants) {
                        if (j >= flist[i].variants) {
                            fragment.alist[fragment.variants].varid = flist[i + 1].alist[k].varid;
                            fragment.alist[fragment.variants].allele = flist[i + 1].alist[k].allele;
                            fragment.alist[fragment.variants].qv = flist[i + 1].alist[k].qv;
                            fragment.alist[fragment.variants].is_bnd = flist[i + 1].alist[k].is_bnd;
                            fragment.variants++;
                            k++;
                            continue;
                        }
                        if (k >= flist[i + 1].variants) {
                            fragment.alist[fragment.variants].varid = flist[i].alist[j].varid;
                            fragment.alist[fragment.variants].allele = flist[i].alist[j].allele;
                            fragment.alist[fragment.variants].qv = flist[i].alist[j].qv;
                            fragment.alist[fragment.variants].is_bnd = flist[i + 1].alist[k].is_bnd;
                            fragment.variants++;
                            j++;
                            continue;
                        }

                        if (flist[i].alist[j].varid < flist[i + 1].alist[k].varid) {
                            fragment.alist[fragment.variants].varid = flist[i].alist[j].varid;
                            fragment.alist[fragment.variants].allele = flist[i].alist[j].allele;
                            fragment.alist[fragment.variants].qv = flist[i].alist[j].qv;
                            fragment.alist[fragment.variants].is_bnd = flist[i + 1].alist[k].is_bnd;

                            fragment.variants++;
                            j++;
                        } else if (flist[i].alist[j].varid > flist[i + 1].alist[k].varid) {
                            fragment.alist[fragment.variants].varid = flist[i + 1].alist[k].varid;
                            fragment.alist[fragment.variants].allele = flist[i + 1].alist[k].allele;
                            fragment.alist[fragment.variants].qv = flist[i + 1].alist[k].qv;
                            fragment.alist[fragment.variants].is_bnd = flist[i + 1].alist[k].is_bnd;

                            fragment.variants++;
                            k++;
                        } else if (flist[i].alist[j].allele == flist[i + 1].alist[k].allele) // consistent
                        {
                            fragment.alist[fragment.variants].varid = flist[i].alist[j].varid;
                            fragment.alist[fragment.variants].allele = flist[i].alist[j].allele;
                            fragment.alist[fragment.variants].qv = flist[i].alist[j].qv;
                            fragment.alist[fragment.variants].is_bnd = flist[i + 1].alist[k].is_bnd;

                            if (flist[i + 1].alist[k].qv > flist[i].alist[j].qv) fragment.alist[fragment.variants].qv = flist[i + 1].alist[k].qv;
                            fragment.variants++;
                            j++;
                            k++;
                        } else {
                            j++;
                            k++;
                        }
                    }
                    if (fragment.variants >= 2 || (SINGLEREADS ==1 && fragment.variants >=1)) {
                        sl = strlen(flist[i].id);
                        fragment.id = (char*) malloc(sl + 1);
                        for (j = 0; j < sl; j++) fragment.id[j] = flist[i].id[j];
                        fragment.id[j] = '\0';
                        fragment.read_qual = flist[i].read_qual;
                        if (flist[i].barcode == NULL){
                            fragment.barcode = NULL;
                        }else{
                            bl = strlen(flist[i].barcode);
                            fragment.barcode = (char*) malloc(bl + 1);
                            for (j = 0; j < bl; j++) fragment.barcode[j] = flist[i].barcode[j];
                            fragment.barcode[j] = '\0';
                        }

                        //for (j=0;j<flist[i].variants;j++) fprintf(stdout,"%d ",flist[i].alist[j].varid); fprintf(stdout,"| ");
                        //for (j=0;j<flist[i+1].variants;j++) fprintf(stdout,"%d ",flist[i+1].alist[j].varid);
                        //fprintf(stdout,"order of variants not correct %s \t",flist[i].id);
                        print_fragment(&fragment, varlist, fragment_file);
                        free(fragment.id);
                    }

                }
                else if (flist[i].variants+flist[i+1].variants ==2 && SINGLEREADS ==1 && flist[i].variants >= 1)print_fragment(&flist[i],varlist,fragment_file); // added 05/31/2017 for OPE

                //else if (flist[i].variants ==1 && flist[i+1].variants >1) print_fragment(&flist[i+1],varlist);
                //else if (flist[i].variants > 1 && flist[i+1].variants ==1) print_fragment(&flist[i],varlist);
                // april 27 2012 these PE reads were being ignored until now
                i += 2;
                // what about overlapping paired-end reads.... reads..... ???? jan 13 2012,
            } else if (flist[i].variants >= 2 || (SINGLEREADS == 1 && flist[i].variants >=1)) {
                print_fragment(&flist[i], varlist, fragment_file);
                i++;
            } else i++;

        }
        // last read examined if it is not paired
        if (i < *fragments) {
            if (flist[i].variants >= 2 || (SINGLEREADS == 1 && flist[i].variants >=1)) print_fragment(&flist[i], varlist, fragment_file);
        }
    } else // only one fragment in fraglist single end
    {
        if (flist[first].variants >= 2 || (SINGLEREADS == 1 && flist[i].variants >=1)) print_fragment(&flist[first], varlist, fragment_file);
    }

    // free the fragments starting from first....
    if (*fragments > 0)// check added jan 13 2012
    {
        for (i = first; i<*fragments; i++) {
            free(flist[i].id);
            free(flist[i].alist);
        }
    }
    (*fragments) = first;
    free(fragment.alist);
}

void sort_framgment(FRAGMENT* fragment) {
    int i,j;
    int n = fragment->variants;
    allele t;
    for (i = 0; i < n; i++) {

        for (j = i + 1; j < n; j++) {

            if ((fragment->alist+ j)->varid < (fragment->alist+ i)->varid) {
                t = *(fragment->alist+ i);
                *(fragment->alist + i) = *(fragment->alist + j);
                *(fragment->alist + j) = t;
            }
        }
    }
}
