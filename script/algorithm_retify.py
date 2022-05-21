#!/usr/bin/env python3
from my_imports import *

def table(k,allele_num):
    mytable=[]
    for j in range(allele_num):
        mytable.append([j])
    if k>1:
        for i in range(k-1):
            double_table=[]
            for j in range(allele_num):
                add_allele=[]
                for list in mytable:
                    newlist=list[:]
                    newlist.append(j)
                    add_allele.append(newlist)
                double_table+=add_allele
            mytable=double_table
        mytable.remove([0]*k)
        mytable.remove([1]*k)
    return mytable

def table_allele(k,max_allele):
    table_set=[]
    for allele_num in range(2,max_allele+1):
        table_set.append(table(k,allele_num))
    return table_set

def init_delta(points_num):
    delta_set=[]
    for i in range(points_num-1):
        delta_set.append([[0,0],[0,0]])
    return delta_set

def random_alpha(k):
    alpha = np.array([np.random.randint(1, 10) for i in range(k)])
    return alpha/sum(alpha)

def fixed(k):
    num_list=[]
    for i in range(k):
        num_list.append(i+1)
    list_sum=sum(num_list)
    fixed_alpha=[]
    for num in num_list:
        fixed_alpha.append(num/list_sum)
    return fixed_alpha

def alpha_step(geno_set,beta_set):
    # compute allele frequency with least square
    locus_num = 0
    alpha = np.array([0.0, 0.0])
    for i in range(len(beta_set)):
        beta = beta_set[i][1]
        if geno_set[i][0] == 0:
            alpha[0] += (1-beta)
            alpha[1] += beta
            locus_num += 1
        elif geno_set[i][0] == 1:
            alpha[0] += beta
            alpha[1] += (1-beta)
            locus_num += 1
    if locus_num > 0:
        return alpha/locus_num
    else:
        return [1, 0]

def index2seq(strain_number,locus_index,allele_set):
    geno_set=[]
    table_set=table_allele(strain_number,max(allele_set))
    # my_table=table(strain_number,allele_num)
    # for index in locus_index:
    for locus in range(len(locus_index)):
        index=locus_index[locus]
        geno_set.append(table_set[allele_set[locus]-2][index])
    geno_set=np.array(geno_set)
    seq_list=np.transpose(geno_set)
    return seq_list

class Phase_step(): #too many same step in the iteration
    def __init__(self,given_alpha,delta_set,beta_set,share_set,weight,allele_set):
        self.k=len(given_alpha)
        self.beta_set=beta_set
        # self.allele_num=len(beta_set[0])
        self.allele_set=allele_set
        self.given_alpha=given_alpha
        self.delta_set=delta_set
        self.points_num=len(beta_set)
        self.table_set=table_allele(self.k,max(allele_set))
        self.table_delta_set=self.all_delta_table()
        self.w=weight
        self.share_set=share_set
    def estimate(self,locus):
        estimated_beta=[]
        locus_table=self.table_set[self.allele_set[locus]-2]
        # for colum in self.table:
        for colum in locus_table:
            # result=[0]*self.allele_num
            result=[0]*self.allele_set[locus]
            for i in range(self.k):
                result[colum[i]]+=self.given_alpha[i]
            estimated_beta.append(result)
        return estimated_beta
    def delta_diff(self,delta_a,delta_b):
        diff=0
        #print (delta_a,delta_b)
        for i in range(len(delta_a)):
            diff+=abs(delta_a[i]-delta_b[i])
        return diff
    def delta_table(self,locus):
        # table_delta=[]
        # for i in range(len(self.table_set[self.allele_set[locus]-2])):
        #     middle_table=[]
        #     for j in range(len(self.table_set[self.allele_set[locus+1]-2])):
        #         # ij_delta=[0]*(self.allele_num**2)
        #         ij_delta=[0]*(self.allele_set[locus]*self.allele_set[locus+1])
        #         for l in range(self.k):
        #             index=self.table_set[self.allele_set[locus]-2][i][l]*self.allele_set[locus+1]+self.table_set[self.allele_set[locus+1]-2][j][l]
        #             # index=self.table[i][l]+self.table[j][l]*self.allele_num
        #             ij_delta[index] += self.given_alpha[l]
        #         middle_table.append(ij_delta)
        #     table_delta.append(middle_table)
        # return table_delta
        # print (self.allele_set[locus]-2,self.allele_set[locus+1]-2)
        return self.table_delta_set[self.allele_set[locus]-2][self.allele_set[locus+1]-2]
    def generate_deltas(self,pre_allele,fol_allele):
        table_delta=[]
        for i in range(len(self.table_set[pre_allele-2])):
            middle_table=[]
            for j in range(len(self.table_set[fol_allele-2])):
                # ij_delta=[0]*(self.allele_num**2)
                ij_delta=[0]*(pre_allele*fol_allele)
                for l in range(self.k):
                    index=self.table_set[pre_allele-2][i][l]*fol_allele+self.table_set[fol_allele-2][j][l]
                    # index=self.table[i][l]+self.table[j][l]*self.allele_num
                    ij_delta[index] += self.given_alpha[l]
                middle_table.append(ij_delta)
            table_delta.append(middle_table)
        return table_delta
    def all_delta_table(self):
        table_delta_set=[]
        allele_num=max(self.allele_set)
        for pre_allele in range(2,allele_num+1):
            one_locus=[]
            for fol_allele in range(2,allele_num+1):
                one_locus.append(self.generate_deltas(pre_allele,fol_allele))
            table_delta_set.append(one_locus)
        return table_delta_set
    def delta_phase(self):
        # self.delta_set[start:end+1],self.beta_set[start:end+2]=graph.lp_correct(self.delta_set[start:end+1],self.beta_set[start:end+2])
        save_table=[]
        # geno_num=len(self.table)
        for r in range(len(self.delta_set)+1): #delta_len = beta_len -1            
            point_table=[]
            if r == 0:
                geno_num=len(self.table_set[self.allele_set[r]-2])
                for m in range(geno_num):
                    point_table.append([0,0])
                save_table.append(point_table)
            else:
                # geno_num=len(self.table_set[self.allele_set[r-1]-2])
                #n is previous locus.
                for m in range(len(self.table_set[self.allele_set[r]-2])):
                    this_geno=[float('inf'),0]
                    for n in range(len(self.table_set[self.allele_set[r-1]-2])):
                        delta_loss=self.delta_diff(self.delta_set[r-1],self.delta_table(r-1)[n][m])
                        ratio=self.beta_set[r]
                        beta_loss=sum(abs(np.array(ratio)-np.array(self.estimate(r)[m])))
                        weight_loss=beta_loss*(self.w)+delta_loss*(1-self.w)
                        add_loss=save_table[r-1][n][0] + weight_loss
                        add_loss=round(add_loss,6)
                        if add_loss<this_geno[0]:
                            this_geno=[add_loss,n]
                    point_table.append(this_geno)
                save_table.append(point_table)
        frag_index,part_loss=self.backtrack(save_table)
        return frag_index,part_loss
    def backtrack(self,save_table):
        geno_num=len(self.table_set[self.allele_set[-1]-2])
        # geno_num=len(self.table)
        reverse_index=[]
        final_geno=[float('inf'),0]
        final_index=0
        for m in range(geno_num):
            this_geno=save_table[-1][m]
            if float(this_geno[0]) < float(final_geno[0]):
                final_geno=this_geno
                final_index=m
        # print ("delta loss",final_geno[0])
        part_loss=final_geno[0]
        reverse_index.append(final_index)
        for r in reversed(range(len(save_table)-1)):
            reverse_index.append(final_geno[1])
            final_geno=save_table[r][reverse_index[-1]]
        reverse_index.reverse()
        answer_index=reverse_index
        return answer_index,part_loss
    def breaks_phase(self):
        answer_index,phase_loss=self.delta_phase()
        geno_set=self.genotype(answer_index)
        return answer_index,geno_set,phase_loss
    def loss(self,answer_index):
        my_loss=0
        for i in range(len(answer_index)):
            my_loss+=abs(self.estimated_beta[answer_index[i]] - self.beta_set[i])
        return my_loss
    def genotype(self,answer_index):
        geno_set=[]
        for i in range(self.points_num):
            # print ('xx',i,len(answer_index),len(self.points_num))
            geno_set.append(self.table_set[self.allele_set[i]-2][answer_index[i]])
        return geno_set

class Workflow():
    def __init__(self,beta_set,delta_set,share_set,weight,elbow,allele_set):
        self.beta_set=beta_set
        self.delta_set=delta_set
        self.allele_set=allele_set
        self.share_set=share_set
        self.w=weight
        self.elbow=elbow
    def given_k(self):
        geno_index,corr_loss,final_alpha=self.multi_init()
        seq_list=index2seq(len(final_alpha),geno_index,self.allele_set)
        return final_alpha,seq_list,corr_loss
    def multi_init(self):
        mini_set=['',float('inf'),'']
        for i in range(10):
            geno_index,corr_loss,final_alpha=self.iteration()
            if corr_loss < mini_set[1]:
                mini_set=[geno_index,corr_loss,final_alpha]
        # print (mini_set[0],mini_set[2])
        # print (table_allele(3,2))
        return mini_set[0],mini_set[1],mini_set[2]
    def iteration(self):
        #T is supposed strain number
        times=0
        past_loss=float('inf')
        save_list=[]
        loss_list=[]
        # current_alpha=fixed(T)
        current_alpha=[0.5, 0.5]
        while True:
            ph=Phase_step(current_alpha,self.delta_set,self.beta_set,self.share_set,self.w,self.allele_set)
            answer_index,geno_set,phase_loss=ph.breaks_phase()
            # print ('phase step done')
            # print (phase_loss,current_alpha)
            save_list.append([answer_index,geno_set,phase_loss,current_alpha])
            loss_list.append(phase_loss)
            if abs(past_loss-phase_loss) < 0.000001 or times > 1:
                final_index,final_loss,final_alpha = answer_index,phase_loss,current_alpha
                break
            past_loss=phase_loss
            times+=1
            current_alpha=alpha_step(geno_set,self.beta_set)
            # print ('alpha step done')
            current_alpha=sorted(current_alpha)
        return final_index,final_loss,final_alpha


if __name__ == "__main__":
    beta_set=[[0.7,0.3],[0.8,0.2],[0.3,0.7],[0.5,0.2,0.3],[0.5,0.5],[0.2,0.3,0.5],[0.5,0.5],[0.3,0.7]]
    # delta_set=[[0.5,0.3,0,0.2,0.0,0,0,0,0],[0.3,0,0,0.5,0.2,0,0,0,0],[0,0.5,0,0,0.2,0,0.3,0,0],[0,0.2,0.3,0.5,0,0,0,0,0],\
    #     [0.2,0,0,0.3,0,0,0,0.5,0],[0,0,0.5,0.2,0.3,0,0,0,0],[0,0.3,0,0.5,0.2,0,0,0,0]]
    delta_set=[[0.5,0.2,0.3,0],[0.3,0.5,0,0.2],[0,0,0.3,0.5,0.2,0],[0,0.5,0.2,0,0.3,0],[0.2,0.3,0,0,0,0.5],[0,0.2,0,0.3,0.5,0],[0,0.5,0.3,0.2]]
    geno_set=[[0,1,0],[1,0,0],[1,0,1],[1,2,0],[0,0,1],[0,1,2],[1,1,0],[1,0,1]]
    allele_set=[2,2,2,3,2,3,2,2]
    wo=Workflow(beta_set,delta_set,allele_set)
    # wo.choose_k()
    final_alpha,seq_list=wo.given_k(3)
    print(np.array(geno_set))
    print (final_alpha)
    print(np.array(seq_list))
