'''
    REG: Random Effect General Method 
    Copyright(C) 2018 Cue Hyunkyu Lee

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''
#!/bin/python3
# dependency: numpy,scipy.stats, scipy.optimize,

import numpy as np
import pandas as pd
import multiprocessing as mp
from itertools import product
from decimal import *
from scipy.stats import multivariate_normal
from scipy.optimize import minimize
from meta_code.regeneral import REG_optim
from meta_code.LS import LS_chi, LS_apply

#os.path.basename(__file__);

### These definitions are necessary to estimate sample probability density functions for deterministic IS 
def generate_Pj (means, stders, cov, nstudy):
    Pj_list = [[],[]];
    for i in range(len(means)):
        mean_array = np.array([means[i]] * nstudy); stder_matrix = np.diag([stders[i]] * nstudy);
        cov_matrix = stder_matrix.dot(cov).dot(stder_matrix);
        Pj_list[0].append(mean_array); Pj_list[1].append(cov_matrix);
    return(Pj_list);

### sampling from mixture densities
def mixture_sampling (nsample, probs, Pj):
    nP = len(probs); choices = [i for i in range(nP)]; N = nsample; count = [0] * nP; df = pd.DataFrame();
    comp = np.random.choice(choices, N, replace = True, p = probs);
    for i in range(len(comp)): count[comp[i]] += 1;
    for i in range(nP):
        cmean = Pj[0][i]; ccov = Pj[1][i]; cN = count[i];
        adf = pd.DataFrame(np.random.multivariate_normal(mean = cmean, cov = ccov, size = cN));
        df = df.append(adf,ignore_index=True);
    return(df);

def get_H (tau, U, R): 
    return(tau * U + R );

def get_Cov(stders, H):
    cov=np.diag(stders).dot(H).dot(np.diag(stders));
    return(cov);

def h_t (ts,thres):
    return_vec = [0] * len(ts);
    for i in range(len(ts)):
        if(ts[i] >= thres): 
            return_vec[i] = 1;
    return(return_vec);

def estim_prob_Pj (Pj, X):
    nPj = len(Pj[0]); pdf_Pj = [];
    for i in range(nPj):
        cmean = Pj[0][i];
        ccov = Pj[1][i];
        cpdf = multivariate_normal.pdf(x = X, mean = cmean, cov = ccov);
        pdf_Pj.append(np.array(cpdf));
    return(pdf_Pj);

### These definitions are necessary for estimating I 
def const_mul(array, pdf_Pj):
    alist = [];
    for i in range(len(array)):
        amul = [array[i] * pdf_Pj[i][j] for j in range(len(pdf_Pj[i]))];
        alist.append(np.array(amul));
    return(alist);

def vector_sum(alist):
    sumvec = [0]*len(alist[0]);
    for i in range(len(alist)):
        sumvec = [sumvec[j] + alist[i][j] for j in range(len(alist[i]))];
    return(sumvec);

def estim_cov_tm(pdf_Pj, m):
    l = len(pdf_Pj);
    tm_vec = [0]*l;
    for i in range(l):
        tm_vec[i] = np.cov(m, pdf_Pj[i])[0][1];
    return(np.array(tm_vec));

def estim_cov_t(pdf_Pj, Palpha):
    l = len(pdf_Pj);
    t_mat = [];
    for i in range(l):
        array = np.array([pdf_Pj[i][j]/Palpha[j] for j in range(len(pdf_Pj[i]))])
        t_mat.append(array)
    return(np.cov(np.array(t_mat)));

def svd_inv(cov_t):
    u,s,v = np.linalg.svd(cov_t);

    ds = np.diag([1/s[j] for j in range(len(s)-1)]);
    us = np.matrix(u)[:,:-1];
    vs = np.matrix(np.transpose(v))[:,:-1];
    inv_cov_t = vs.dot(ds).dot(np.transpose(us));
    return(inv_cov_t)

def ims_estimate_statistics(df_data, se, Sg, Re, n):
    df_out = pd.DataFrame(index = df_data.index)
    df_out['stats'] =df_data.apply(lambda x: REG_optim(x.tolist(), se, Sg, Re, n), axis=1)
    return(df_out)

def ims_parallelize(df_input, func, cores, partitions, se, Sg, Re, n):
    data_split = np.array_split(df_input, partitions)
    iterable = product(data_split, [se], [Sg], [Re], [n])
    pool = mp.Pool(int(cores))
    df_output = pd.concat(pool.starmap(func, iterable))
    pool.close()
    pool.join()
    return(df_output)


def thres_estimate_pvalue(thres, REG_X, Palpha, probs, d_Q, d_Pj, nPj, N):
    h_reg = h_t(ts = REG_X, thres = thres);
    m_reg = [h_reg[i] * d_Q[i] / Palpha[i] for i in range(len(d_Q))];
    cov_tm_reg = estim_cov_tm(d_Pj, m_reg); 
    cov_t = estim_cov_t(d_Pj, Palpha);
    inv_cov_t = svd_inv(cov_t);
    denominator = vector_sum(const_mul(probs, d_Pj));
    betas_reg = [inv_cov_t.dot(cov_tm_reg)[0,i] for i in range(nPj)];
    control_variate_reg = vector_sum(const_mul(betas_reg, d_Pj));
    nominator_reg = [ h_reg[i] * d_Q[i] - control_variate_reg[i] for i in range(len(d_Q)) ];
    res = sum( nominator_reg[i] / denominator[i] for i in range(len(d_Q))) /N + np.sum(betas_reg);
    return(pd.DataFrame([res], columns = ['pvalue'], index = [thres]))


def thres_parallelize(thres_vec, func, cores, REG_X, Palpha, probs, d_Q, d_Pj, nPj, N):
    iterable = product(thres_vec, [REG_X], [Palpha], [probs], [d_Q], [d_Pj], [nPj], [N])
    pool = mp.Pool(int(cores))
    res_list = pd.concat(pool.starmap(func, iterable))
    pool.close()
    pool.join()
    return(res_list)

## GenCor and RECor: np.matrix, N: int, outfn: str
def importance_sampling(N, Sg, Re, outfn, mp_cores):
    
    ### set multi processing options 
    #mp.set_start_method('spawn')
    if(mp_cores == 0):
        cores = mp.cpu_count();
        partitions = cores;
    else:
        cores = mp_cores;
        partitions = cores;

    output_filename = outfn;
    
    n = nstudy = Sg.shape[0]; se=[1]*n;
   
    ## Deterministic importance sampling needs probability density functions for sampling purposes. They have means of zeros and standard errors(s_stders)
    std_P = [1,1.1,1.2,1.3,1.4,1.7,2,2.5,3,4,5];

    nPj = len(std_P); mean_P = [0]*nPj; probs= [1/nPj]*nPj; H = get_H(tau=0, U = Sg, R = Re);  
    Pj = generate_Pj(means = mean_P, stders = std_P, cov = H, nstudy = n);
    
    ## generate sample X
    df_input = mixture_sampling(nsample = N, probs=probs, Pj=Pj)
    print( "Generating {len_X} stats.".format( len_X=N ) ); 
    data = ims_parallelize(df_input, ims_estimate_statistics, cores, partitions, se, Sg, Re, n)
    REG_X = data['stats'].tolist()
   
    null_Re = Re; null_means = [0]*n; null_std = np.diag([1]*n);
    null_cov = null_std.dot(null_Re).dot(null_std);

    d_Q = multivariate_normal.pdf(x = df_input, mean = null_means, cov = null_cov);
    d_Pj = estim_prob_Pj (Pj, X = df_input);
   
    ### It is recommended to get tabulated pdf at 0.1, 0.2, 0.3 ... 1.0, 2.0, 3.0,.... 31.0.  
    thres_vec = [ float(i)/10 for i in range(10) ] + [ float(i+1) for i in range(40) ]; ntest = len(thres_vec);
    pvalue_estim = pd.DataFrame(index = thres_vec);

    Palpha = vector_sum(const_mul(probs,d_Pj));
   
    df_pvalue = thres_parallelize(thres_vec, thres_estimate_pvalue, cores, REG_X, Palpha, probs, d_Q, d_Pj, nPj, N)
    print('Completed estimating the inverse CDFs at different threshold values for DELPY\'s variance component model.'); 
    
    sorted_pvalue = df_pvalue.sort_index()
    print(sorted_pvalue)
    sorted_pvalue.to_csv(output_filename, header = False, index = True, sep = " ")
    print('Wrote tabulated inverse cdf');
