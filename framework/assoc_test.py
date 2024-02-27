from framework.significance_estimation import cof_estimation, pvalue_estimation, flattening_p_value 
from meta_code.variance_component import vcm_optimization
from framework.mtag_estimation import mtag
from framework.utilities import sqrt_ginv
from meta_code.LS import LS, LS_p 
from itertools import product
import multiprocessing as mp
import pandas as pd
import numpy as np

def parallel_computing(iterable, function, ncpu): 
    '''
    Parallel computation using starmap function 
    '''
    pool = mp.Pool(ncpu)
    results = pd.concat(pool.starmap(function, iterable))
    pool.close()
    pool.join()
    return(results)

def blup_optimization_old(y, G, R, c):
    pinv_R = np.linalg.pinv(R)
    pinv_G = np.linalg.pinv(G)
    V = np.linalg.pinv(pinv_R + pinv_G)
    u = V.dot(pinv_R).dot(y)
    return(pd.Series(u, index = c))

def blup_optimization(y, G, R, c):
    pinv_R = np.linalg.pinv(R)
    pinv_G = np.linalg.pinv(G)
    V = np.linalg.pinv(pinv_R + pinv_G)  # This is the key step for variance estimation
    u = V.dot(pinv_R).dot(y)
    
    # Calculating standard errors: sqrt of diagonal elements of V
    se_u = np.sqrt(np.diag(V))

    # Combine u and se_u into a single vector
    combined_values = np.empty(u.size + se_u.size, dtype=u.dtype)
    combined_values[0::2] = u
    combined_values[1::2] = se_u

    # Generate combined index strings
    combined_index = sum([[x, x + '_se'] for x in c], []) 

    return pd.Series(combined_values, index=combined_index)

def stat_estimation(sumstat, sg, ce, isf):
    def run_PLEIO(arr, ind):
        b=arr[ind]
        s=np.diag(arr[ind+1])
        return(vcm_optimization(np.transpose(u).dot(b), np.transpose(u).dot(s).dot(ce).dot(s).dot(u)))
        
    def run_LS(arr, ind, ce):
        b=arr[ind]
        s=arr[ind+1]
        return(LS(b, s, ce))
        
    n = len(np.diag(sg))
    ind = np.array([i*2 for i in range(n)])
    u = sqrt_ginv(sg); c = ce
    res = pd.DataFrame(index = sumstat.index, columns = ['tausq','pleio_stat','LS_stat','pleio_p','LS_p'])
    pleio_res = sumstat.apply(lambda x: run_PLEIO(x.to_numpy(), ind), result_type = 'expand', axis=1)
    res.loc[sumstat.index,['tausq','pleio_stat']] = pleio_res.loc[sumstat.index,['tausq','pleio_stat']]
    del pleio_res
    res.loc[:,'LS_stat'] = sumstat.apply(lambda x: run_LS(x.to_numpy(), ind, ce), axis=1)
    pfunction = cof_estimation(isf);
    res.loc[:,'pleio_p'] = res.loc[:,'pleio_stat'].apply(lambda x: pvalue_estimation(x, pfunction));
    res.loc[:,'LS_p'] = res.loc[:,'LS_stat'].apply(lambda x: LS_p(x));
    return(res)

def blup_estimation(d, sg, ce, trait_name):
    def run_BLUP(x, sg, ce, ind, c):
        y = x[ind]
        se = x[ind+1]
        tausq = np.max([10**-12,x[-1]])
        G = np.multiply(tausq, sg) 
        R = np.diag(se).dot(ce).dot(np.diag(se))
        return(blup_optimization(y,G,R,c))

    n = np.size(sg,1)
    t = trait_name
    ind = np.array([i*2 for i in range(n)], dtype = int)
    d = d.apply(lambda x: run_BLUP(x.to_numpy(dtype = float), sg, ce, ind, t), result_type = 'expand', axis=1)
    return(d)

