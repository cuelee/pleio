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

import os, sys, math
import numpy as np
from decimal import *
#import pandas as pd
from scipy.stats import multivariate_normal
from scipy.optimize import minimize
from meta_code.regeneral import REG_optim, REG_apply
from meta_code.LS import LS_chi, LS_apply
from importance_sam.etc.time import *

#os.path.basename(__file__);

### These definitions are necessary for estimating pdfs
def setPj (s_means, s_stders, H, nstudy):
    Pj_list = [[],[]];
    for i in range(len(s_means)):
        mean_array = np.array([s_means[i]] * nstudy);
        stder_mat = np.diag([s_stders[i]] * nstudy);
        cov_mat = stder_mat.dot(H).dot(stder_mat);
        Pj_list[0].append(mean_array);
        Pj_list[1].append(cov_mat);
    return(Pj_list);

## sampling from mixture densities
def mixt_sampling (N, probs, Pj):
    nPj = len(probs); choices = [i for i in range(nPj)];
    comp = np.random.choice(choices, N, replace = True, p = probs);
    array = np.empty((0,len(Pj[0][0])), int); count = [0] * nPj;
    for i in range(len(comp)): count[comp[i]] += 1;
    for i in range(nPj):
        cmean = Pj[0][i]; ccov = Pj[1][i]; cN = count[i];
        appended = gen_sample(N = cN, means = cmean, cov = ccov);
        array=np.append(array,appended,axis=0);
    return(array);

def get_H (tau, U, R): 
    return(tau * U + R );

def get_Cov(stders, H):
    cov=np.diag(stders).dot(H).dot(np.diag(stders));
    return(cov);

def gen_sample (N, means, cov):
    ## mean : ints, cov : 2D matrix, size : int
    array = np.random.multivariate_normal(mean = means, cov = cov, size = N);
    return(array);

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

## GenCor and RECor: np.matrix, N: int
def run_importance_sampling(N, GenCov, RECor, outfn):
    print('\n Cue Hyunkyu Lee: {}'.format(os.path.basename(__file__)));
    start_time = get_time();
    cur_time(' Job started at');
    outfile = outfn;

    n = nstudy = GenCov.shape[0];
    
    s_stders = [1,1.1,1.2,1.3,1.4,1.7,2,2.5,3,4,5];
    nPj = len(s_stders);
    s_means = [0]*nPj;
    probs= [1/nPj]*nPj;
    
    Sg = GenCov; Re = RECor;
    
    tau = 0;
    H = get_H(tau, Sg, Re);
    
    Pj = setPj(s_means = s_means, s_stders = s_stders, H = H, nstudy = n);
    
    null_cor = Re;
    null_means = [0]*n;
    null_stders_mat = np.diag([1]*n);
    null_cov = null_stders_mat.dot(null_cor).dot(null_stders_mat);
    
    X = mixt_sampling (N = N, probs = probs, Pj=Pj);
    
    print(" Finished sampling of {} variants".format(len(X)));
    
    REG_X = REG_apply(Sg=Sg, Re=Re, n=n, X=X, row_wise = True);
    
    pdf_Q = multivariate_normal.pdf(x = X, mean = null_means, cov = null_cov);
    pdf_Pj = estim_prob_Pj (Pj, X=X);
    
    thres_vec = [ float(i)/10 for i in range(10) ] + [ float(i+1) for i in range(30) ];
    no_test = len(thres_vec)
    reg_estim = [float(0)] * no_test;    

    k = 0;
    for thres in thres_vec:

        Palpha = vector_sum(const_mul(probs,pdf_Pj));
        h_reg = h_t(ts = REG_X, thres = thres);
        m_reg = [h_reg[i] * pdf_Q[i] / Palpha[i] for i in range(len(pdf_Q))];
        cov_tm_reg = estim_cov_tm(pdf_Pj, m_reg); 

        cov_t = estim_cov_t(pdf_Pj, Palpha);
        inv_cov_t = svd_inv(cov_t);
        denominator = vector_sum(const_mul(probs, pdf_Pj));

        betas_reg = [inv_cov_t.dot(cov_tm_reg)[0,i] for i in range(nPj)];
        control_variate_reg = vector_sum(const_mul(betas_reg, pdf_Pj));
        nominator_reg = [ h_reg[i] * pdf_Q[i] - control_variate_reg[i] for i in range(len(pdf_Q)) ];
        reg_estim[k] = sum( nominator_reg[i] / denominator[i] for i in range(len(pdf_Q))) /N + np.sum(betas_reg);
        print("\t{}: {}".format(thres,reg_estim[k]));

        k=k+1;

    ## print(output)
    print('\n Importance samping result');
    with open(outfile,'w') as fin:
        [print('{} {}'.format(thres_vec[i],reg_estim[i]),file=fin) for i in range(no_test)];

    run_time = get_time() - start_time;
    print(' Analysis time: {} sec'.format(round(run_time,3)));
