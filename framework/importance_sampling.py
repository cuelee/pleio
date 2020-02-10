'''
    Copyright(C) 2018 Cue Hyunkyu Lee
'''
#!/bin/python3
# dependency: numpy,scipy.stats, scipy.optimize,

import numpy as np
import pandas as pd
import multiprocessing as mp
import random
from itertools import product
from decimal import *
from scipy.stats import multivariate_normal
from meta_code.variance_component import vcm_optimization
from meta_code.LS import LS_chi, LS_apply

#os.path.basename(__file__);


def generate_P(means, stders, Re, n, nPj):
    '''
    P is the list contains 
    '''

    class Pj(object):
        '''
        Pj is a sampling density function of the deterministic importance sampling procedure 
        The class define the covariance matrix and means of Pj
        '''
        def __init__(self, mean, cov):
            self.means = [mean] * cov.shape[0];
            self.cov = cov;
            self.pdf = None;

    P = [Pj(means[i], np.diag( [stders[i]] * n ).dot( Re ).dot( np.diag([stders[i]] * n) )) for i in range(nPj) ];
    return(P) 

### These definitions are necessary to estimate sample probability density functions for deterministic IS 
def generate_Pj (means, stders, cov, nstudy):
    Pj_list = [[],[]];
    for i in range(len(means)):
        mean_array = np.array([means[i]] * nstudy); stder_matrix = np.diag([stders[i]] * nstudy);
        cov_matrix = stder_matrix.dot(cov).dot(stder_matrix);
        Pj_list[0].append(mean_array); Pj_list[1].append(cov_matrix);
    return(Pj_list);

### sampling from mixture densities
def mixture_sampling (nsample, alpha, P):
    J = len(alpha); choices = [i for i in range(J)]; N = nsample; count = [0] * J; input_df = pd.DataFrame();
    comp = np.random.choice(choices, N, replace = True, p = alpha);
    for i in range(len(comp)): count[comp[i]] += 1;
    for j in range(J):
        Pj_mean = P[j].means; Pj_cov = P[j].cov; Pj_N = count[j];
        Pj_df = pd.DataFrame(np.random.multivariate_normal(mean = Pj_mean, cov = Pj_cov, size = Pj_N));
        P[j].pdf = Pj_df
        input_df = input_df.append(Pj_df,ignore_index=True);
    return(input_df);

def h_t (ts,thres):
    return_vec = [0] * len(ts);
    for i in range(len(ts)):
        if(ts[i] >= thres): 
            return_vec[i] = 1;
    return(return_vec);

def P_density_estimation (P, input_df):
    nP = len(P); pdf_P = [];
    for i in range(nP):
        Pj_mean = P[i].means;
        Pj_cov = P[i].cov;
        Pj_pdf = multivariate_normal.pdf( x = input_df, mean = Pj_mean, cov = Pj_cov );
        pdf_P.append( np.array( Pj_pdf ) );
    return( pdf_P );

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
    df_out['LL_RTS'] = df_data.apply(lambda x: vcm_optimization(x.tolist(), se, Sg, Re, n), axis=1)
    return(df_out)

def ims_parallelize(df_input, func, cores, partitions, se, Sg, Re, n):
    data_split = np.array_split(df_input, partitions)
    iterable = product(data_split, [se], [Sg], [Re], [n])
    pool = mp.Pool(int(cores))
    df_output = pd.concat(pool.starmap(func, iterable))
    pool.close()
    pool.join()
    return(df_output)

def thres_estimate_pvalue(thres, Sdelpy, Palpha, alpha, d_Q, d_P, nPj, N):
    h = h_t(ts = Sdelpy, thres = thres);
    m = [h[i] * d_Q[i] / Palpha[i] for i in range(len(d_Q))];
    cov_tm = estim_cov_tm(d_P, m); 
    cov_t = estim_cov_t(d_P, Palpha);
    inv_cov_t = svd_inv(cov_t);
    denominator = vector_sum(const_mul(alpha, d_P));
    betas = [inv_cov_t.dot(cov_tm)[0,i] for i in range( nPj )];
    control_variate = vector_sum(const_mul(betas, d_P));
    nominator = [ h[i] * d_Q[i] - control_variate[i] for i in range(len( d_Q )) ];
    IS_estim = sum( nominator[i] / denominator[i] for i in range(len( d_Q ))) /N + np.sum( betas );
    return(pd.DataFrame([IS_estim], columns = ['pvalue'], index = [thres]))


def thres_parallelize(thres_vec, func, cores, Sdelpy, Palpha, alpha, d_Q, d_P, nPj, N):
    iterable = product(thres_vec, [Sdelpy], [Palpha], [alpha], [d_Q], [d_P], [nPj], [N])
    pool = mp.Pool(int(cores))
    res_list = pd.concat(pool.starmap(func, iterable))
    pool.close()
    pool.join()
    return(res_list)

## GenCor and RECor: np.matrix, N: int, outfn: str
def importance_sampling(N, Sg, Re, outfn, mp_cores):
    np.random.seed(1)
    
    ### set multi processing options 
    #mp.set_start_method('spawn')
    if(mp_cores == 0):
        cores = mp.cpu_count();
        partitions = cores;
    else:
        cores = mp_cores;
        partitions = cores;

    output_filename = outfn;
   
    nstudy = Sg.shape[0]; n = nstudy; Se = np.diag([1]*n).dot(Re).dot(np.diag([1]*n)); se = [1] * n;
   
    ## Deterministic importance sampling needs probability density functions for sampling purposes. They have means of zeros and standard errors(s_stders)
    c_Pj = [1,1.1,1.2,1.3,1.4,1.7,2,2.5,3,4,5];

    nPj = len(c_Pj); mean_P = [0] * nPj; alpha= [1 / nPj] * nPj; H = Se;  
    P = generate_P(means = [0] * nPj, stders = c_Pj, Re = Re, n = n, nPj = nPj)
    #Pj = generate_Pj(means = mean_P, stders = c_Pj, cov = H, nstudy = n);
    
    ## generate sample X
    input_df = mixture_sampling ( nsample = N, alpha = alpha, P = P )
    print( "Generating {len_X} stats.".format( len_X=N ) ); 
    data = ims_parallelize( input_df, ims_estimate_statistics, cores, partitions, se, Sg, Re, n )
    Sdelpy = data['LL_RTS'].tolist()

    null_means = [0] * n;   null_cov = Se;

    d_Q = multivariate_normal.pdf( x = input_df, mean = null_means, cov = null_cov );
    d_P = P_density_estimation( P = P, input_df = input_df );
   
    ### It is recommended to get tabulated pdf at 0.1, 0.2, 0.3 ... 1.0, 2.0, 3.0,.... 31.0.  
    thres_vec = [ float(i) / 10 for i in range(10) ] + [ float(i+1) for i in range(40) ];

    Palpha = vector_sum(const_mul(alpha,d_P));
   
    pvalue_df = thres_parallelize(thres_vec, thres_estimate_pvalue, cores, Sdelpy, Palpha, alpha, d_Q, d_P, nPj, N)
    print('Completed estimating the inverse CDFs at different threshold values for PLEIO\'s variance component model.'); 
    
    sorted_pvalue = pvalue_df.sort_index()
    print(sorted_pvalue)
    sorted_pvalue.to_csv(output_filename, header = False, index = True, sep = " ")
    print('Wrote tabulated inverse cdf');
