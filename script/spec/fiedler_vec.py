import numpy as np
from scipy.sparse.linalg import eigsh 

def get_fiedler_vec(mat):
    D = np.diag(np.sum(mat, axis= 0))
    L = np.matrix(D - mat)
    vals, vecs = eigsh(L, k=2, which='SM')
    fiedler_vec = vecs[:,[1]]
    return fiedler_vec

def constr_graph(freqs): # construct the linkage graph with the allele frequencies at each locus.
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
        if fiedler_vec[i] >0:
            hap.append(0)
        else:
            hap.append(1)
    return hap


if __name__ == "__main__":
    # mat = [ [0,0,0.18,0.82,0,0],
    #         [0,0,0.82,0.18,0,0],
    #         [0.18,0.82,0,0,0.18,0.82],
    #         [0.82,0.18,0,0,0.82,0.18],
    #         [0,0,0.18,0.82,0,0],
    #         [0,0,0.82,0.18,0,0]
    #         ]
    # mat = np.array(mat)
    freqs = [0.9, 0.3, 0.8, 0.4, 0.9]
    mat = constr_graph(freqs)
    print (mat)
    fiedler_vec = get_fiedler_vec(mat)
    hap = get_hap(fiedler_vec)
    print (fiedler_vec)
    print (hap)
