"""
Link the unlinked blocks by mapping them to the database
regard each sub-haplotype in the phase block as a node
construct the linkage graph to link blocks
"""

import numpy as np
from scipy.sparse.linalg import eigsh 
import sys
import os

def get_fiedler_vec(mat):
    D = np.diag(np.sum(mat, axis= 0))
    L = np.matrix(D - mat)
    vals, vecs = eigsh(L, k=2, which='SM')
    fiedler_vec = vecs[:,[1]]
    return fiedler_vec

def constr_graph_freq(freqs): # construct the linkage graph with the allele frequencies at each locus.
    locus_num = len(freqs)
    mat = np.zeros((2*locus_num, 2*locus_num))
    for i in range(locus_num):
        for j in range(i+1, locus_num):
            same = freqs[i] *  freqs[j] + (1-freqs[i]) *  (1-freqs[j])
            diff = freqs[i] *(1-freqs[j]) + (1-freqs[i]) * freqs[j]
            same = round(same, 3)
            diff = round(diff, 3)
            mat[2*i][2*j] = same
            mat[2*i+1][2*j+1] = same
            mat[2*i+1][2*j] = diff
            mat[2*i][2*j+1] = diff

            mat[2*j][2*i] = same
            mat[2*j+1][2*i+1] = same
            mat[2*j+1][2*i] = diff
            mat[2*j][2*i+1] = diff
    return mat

def get_hap(fiedler_vec):
    hap = []
    for i in range(0, len(fiedler_vec), 2):
        if fiedler_vec[i] > fiedler_vec[i+1]:
            hap.append(1)
        else:
            hap.append(0)
    return hap

def constr_graph(score_file, block_phase_file):
    frag_index = {}
    frag_list = []
    index = 0
    f = open(score_file)
    for line in f:
        line = line.strip()
        if line == '':
            continue
        if line[0] == '#':
            continue
        array = line.split()
        frag1 = array[0]
        frag2 = array[1]
        if frag1 not in frag_index:
            frag_index[frag1] = index
            frag_list.append(frag1)
            index += 1
        if frag2 not in frag_index:
            frag_index[frag2] = index
            frag_list.append(frag2)
            index += 1
    f.close()
    locus_num = len(frag_index)
    if locus_num == 0:
        os.system(":> %s"%(block_phase_file))
        return 0
    mat = np.zeros((2*locus_num, 2*locus_num))

    f = open(score_file)
    for line in f:
        line = line.strip()
        if line == '':
            continue
        if line[0] == '#':
            continue
        array = line.split()
        frag1 = array[0]
        frag2 = array[1]
        if len(array) < 3:
            score1 = 0
            score2 = 0
        else:
            new_array = array[2].split(";")
            score1 = float(new_array[0]) # 00 and 11
            score2 = float(new_array[1]) # 01 and 10
        frag1_index = frag_index[frag1]
        frag2_index = frag_index[frag2]

        mat[2*frag1_index][2*frag2_index] = score1
        mat[2*frag1_index+1][2*frag2_index+1] = score1
        mat[2*frag1_index+1][2*frag2_index] = score2
        mat[2*frag1_index][2*frag2_index+1] = score2

        mat[2*frag2_index][2*frag1_index] = score1
        mat[2*frag2_index+1][2*frag1_index+1] = score1
        mat[2*frag2_index+1][2*frag1_index] = score2
        mat[2*frag2_index][2*frag1_index+1] = score2
    # print ("block linkage graph:\n", mat)
    f.close()
    
    fiedler_vec = get_fiedler_vec(mat)
    hapotype = get_hap(fiedler_vec)
    # print ("fiedler_vec:\n", fiedler_vec, hapotype)
    output(frag_list, hapotype, block_phase_file)
    # print (frag_list)

def output(frag_list, hapotype, block_phase_file):
    f = open(block_phase_file, 'w')
    print ("#gene\t start\t end\t hap1\t hap2\t", file = f)
    for i in range(len(frag_list)):
        array = frag_list[i].split(":")
        gene = array[0]
        new_array = array[1].split("-")
        start = int(new_array[0])
        end = int(new_array[1])
        print (gene, start, end , hapotype[i] , 1-hapotype[i], sep = "\t", file = f)
    f.close()


if __name__ == "__main__":
    # score_file = "/mnt/d/HLAPro_backup/HLAPro/example/whole/output/HG00118/HLA_DQB1_break_points_score.txt"
    # block_phase_file = "/mnt/d/HLAPro_backup/HLAPro/example/whole/output/HG00118/HLA_DQB1_break_points_phased.txt"
    score_file = sys.argv[1]
    block_phase_file = sys.argv[2]
    constr_graph(score_file, block_phase_file)
    # output(frag_list, hapotype, block_phase_file)