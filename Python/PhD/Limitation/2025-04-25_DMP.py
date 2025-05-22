import numpy as np
import scipy as sp
from scipy.linalg import block_diag
from qutip import *
np.set_printoptions(legacy='1.13',linewidth=np.inf, threshold=np.inf, precision=3)

import matplotlib.pyplot as plt

#Gibbs energy dist from beta = 1 case, from the RSJJ model calculation value, single mode coup 1 gamma 1.4
Gibbs = np.array([[0.58412,0,0],[0,0.278235,0],[0,0,0.137646]])
degeneracy = np.array([4,4,4])

#Gibbs state with the degeneracy (Ergodic restriction concerned)
def degen_Gibbs_create(Gib, deg, TF):
    dim = Gib.shape[0]
    
    if dim == len(deg):
        sys_gibbs = []
        for i in range(dim):
        ###### debug ######
            #print(Gibbs[i][i])
            #print(degeneracy[i])
            #print("----------")
        ###### debug ######
            temp = (Gib[i][i], deg[i])
            div = temp[0]/temp[1]

            make_matrix_dim = np.identity(temp[1])
            make_matrix_dim *= div
            
            sys_gibbs.append(make_matrix_dim)
            ###### debug ######
            #print(sys_gibbs)
            ###### debug ######
        ret_matrix = block_diag(*sys_gibbs)
        if TF == False:
            return ret_matrix
        if TF == True:
            return sys_gibbs

    else:
        print("Not available input source")
        return None
    
'''
def slicing_Func(Mat):
    if Mat.shape[0] == Mat.shape[1]:
        arr = []
        for i in range(Mat.shape[0]):
            index = 0
            for j in range(Mat.shape[0]):
                if i==j and index < 2:
                    if Mat[index][index] != Mat[index+1][index+1]:
                        arr.append(index)
                        index += 1
                        break
    return(arr)
'''

Cal_Gibbs = degen_Gibbs_create(Gibbs,degeneracy,False)
Cal_Gibbs

Cal_Gibbs_T = degen_Gibbs_create(Gibbs,degeneracy,True)
Cal_Gibbs_T
#print(Cal_Gibbs)

#Create Qubit system with coherence
asdf = coherent_dm(2,1+0.2j)
Sys_block = asdf.full()

c = Bloch()

pointx = expect(sigmax(), asdf)
pointy = expect(sigmay(), asdf)
pointz = expect(sigmaz(), asdf)

c.add_points([pointx, pointy, pointz])
c.add_vectors([pointx, pointy, pointz])

#System construction, False case
Nup_asdf_qobj = Qobj(Sys_block)
Cal_Gibbs_qobj = Qobj(Cal_Gibbs)
rhoRS_qutip = tensor(Nup_asdf_qobj, Cal_Gibbs_qobj)
rhoRS=rhoRS_qutip.full()

def System_block_arrange(sys,bath,tot):
    sys_index = (sys.shape[0])
    bath_index = len(bath)

    #print("system_index",sys_index)

    index_comb = []

    for i in range(sys_index):
        for j in range(bath_index):
            comb = ((i,i),j)
            index_comb.append(comb)

    sorted_comb = sorted(index_comb, key=lambda x: x[0][0] + x[1])
    avail_values = [sorted_comb[i][0][0] + sorted_comb[i][1] for i in range(len(sorted_comb))] #possible combinations for indices corresponding energy conditions

    value_bound = (min(avail_values),max(avail_values))
    print("limited values : ", value_bound)

    avail_values = avail_values[1:-1]
    sorted_comb = sorted_comb[1:-1]
    print("sorted result : ",sorted_comb)
    print("avail_values : ",avail_values)

    sub_block_arr = []
    for k in range(len(sorted_comb)):
        sub_block = np.zeros(bath[sorted_comb[k][1]].shape[0])
        # This condition only acts under qubit case
        #print(sys[sorted_comb[k][0][0],sorted_comb[k][0][1]])
        sub_block = (sys[sorted_comb[k][0][0],sorted_comb[k][0][1]]) * bath[sorted_comb[k][1]]
        sub_block_arr.append(sub_block)
                
    print(sub_block_arr)
    #print(len(sub_block_arr))

    arrsize = 0 
    for l in range(len(sub_block_arr)):
        arrsize += sub_block_arr[l].shape[0]

    total_block = np.zeros((arrsize,arrsize))
    
    # 기존 오류 있는 루프 구조 대신 순차적 블록 배치
    row_offset = 0
    col_offset = 0

    for sub_block in sub_block_arr:
        size = sub_block.shape[0]
        total_block[row_offset:row_offset+size, col_offset:col_offset+size] = sub_block
        row_offset += size
        col_offset += size

    return total_block

testbed = System_block_arrange(Sys_block,Cal_Gibbs_T,rhoRS)

plt.imshow(testbed.real)
plt.colorbar()