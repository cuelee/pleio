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
    def run_BLUP(x, Sg, Ce, ind, c):
        y = x[ind]
        se = x[ind+1]
        tausq = np.max([10**-12,x[-1]])
        G = np.multiply(tausq, Sg) 
        R = np.diag(se).dot(Ce).dot(np.diag(se))
        return(mtag(y,G,R,c))

    n = np.size(Sg,1)
    t = trait_name
    ind = np.array([i*2 for i in range(n)])
    d = df_data.apply(lambda x: run_BLUP(x.to_numpy(), Sg, Ce, ind, t), result_type = 'expand', axis=1)
    return(d)

