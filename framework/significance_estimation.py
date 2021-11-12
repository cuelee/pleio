from scipy.interpolate import splev, splrep
from scipy.stats import chi2
from decimal import *
import pandas as pd
import numpy as np
import multiprocessing as mp
from itertools import product
from numbers import Number
from framework.utilities import sqrt_ginv

def readf(f):
    d = pd.read_csv(f,names=['x','s'],sep=' ');
    d.x = d.x.round(8);
    d.s.iloc[0] = 1;
    return(d);

def manual_estimation(x1,y1,x0=0,y0=1):
    c = (y1-y0)/(x1-x0);
    mestim=lambda x: 1+x*c;
    return(mestim);

def interpolationf(d):
    tck = splrep(d.x.values,np.log(d.s.values), s = 0);
    return(tck);

def extrapolationf(d):
    x = d.x.values;
    y = np.log(d.s.values);
    b = np.cov(x,y)[1,0]/np.var(x);
    a = y[-1]-b*x[-1];
    mestim = lambda x: a + b * x;
    return(mestim);

def pvalue_estimation(s,iso):
    if(not isinstance(s,Number)):
        raise ValueError('The value of the input must be numeric')
    if (s <= iso.min):
        return(Decimal(iso.low(s)));
    elif (s <= iso.max):
        return(np.exp(Decimal(float(splev(s,iso.itck,der=0,ext=2)))));
    else:
        return(np.exp(Decimal(iso.tail(s))));

class cof_estimation(): 
    def __init__(self,isf):
        d = readf(isf);
        self.max = d.x.iloc[-1];
        self.min = d.x.iloc[1];
        self.low = manual_estimation(d.x.iloc[1], d.s.iloc[1]);
        self.itck = interpolationf(d.iloc[1:,:]);
        self.tail = extrapolationf(d.loc[d.x >= 20,:])


### vcm_optimization
def LL_fun(x,n,P_sq,w):
    return(-0.5*(n*np.log(2*np.pi)+sum(np.log(w+x))+sum(P_sq/(w+x))));

def LLp_fun(x,P_sq,w):
    return(0.5*(sum(1/(w+x))-sum(P_sq/(w+x)**2)));

def LLdp_fun(x,P_sq, w):
    return(-0.5*(sum(1/(w+x)**2)-2*sum(P_sq/(w+x)**3)));

def NR_root(f, df, x, P_sq, w, i = 0, iter_max = 10000, tol = 2.22044604925e-16**0.5):
    while ( abs(f(x,P_sq,w)) > tol ):
        x = x - f(x,P_sq,w) / df(x,P_sq,w);
        i = i + 1;
        if (i == iter_max):
            break;
    return(x)

def vcm_optimization (b, n, w, t_v):
    t = [10**(i/4) for i in range(-36,24,1)];
    crossP = t_v.dot(b);
    P_sq = crossP**2;
    init = t[np.argmax([LL_fun(i, n, P_sq, w) for i in t])];
    mle_tausq = NR_root(LLp_fun, LLdp_fun, init, P_sq, w);

    if (mle_tausq <0):
        mle_tausq = 0;
    null_ll = LL_fun(0, n, P_sq, w);
    alt_ll = LL_fun(mle_tausq, n, P_sq, w) ;
    if(alt_ll < null_ll):
        mle_tausq = 0;
        alt_ll = null_ll;
    return (- 2 * (null_ll - alt_ll))

def estimate_statistics(df_data, n, w, t_v):
    df_out = pd.DataFrame(index = df_data.index)
    df_out['null_stat'] = df_data.apply(lambda x: vcm_optimization(x.tolist(), n, w, t_v), axis=1)
    return(df_out)

def parallelize(df_input, func, cores, partitions, n, w, t_v):
    data_split = np.array_split(df_input, partitions)
    iterable = product(data_split, [n], [w], [t_v])
    pool = mp.Pool(int(cores))
    df_output = pd.concat(pool.starmap(func, iterable))
    pool.close()
    pool.join()
    return(df_output)

def flattening_p_value(summary, gwas_N, gencov, envcor, cores, isf, tol = 2.22044604925e-16**0.5):
    ### set multi processing options 
    if(cores == 0):
        cores = mp.cpu_count() - 1; partitions = cores;
    else:
        partitions = cores;

    U = gencov
    Ce = envcor
    se = 1/(np.array(gwas_N)**0.5)
    np.random.seed(1)
    n = len(se); nsim = 100000; 
    D = np.diag(se).dot(Ce).dot(np.diag(se));null_D = np.diag([1]*n).dot(Ce).dot(np.diag([1]*n));
    sqrt_U_inv = sqrt_ginv(U);
    K = sqrt_U_inv.dot(D).dot(sqrt_U_inv)
    w, v = np.linalg.eigh(K); t_v = np.transpose(v)
    pos = w > max(tol * w[0], 0)
    w_pos = w[pos]
    t_v_pos = t_v[pos]

    null_df = pd.DataFrame(np.random.multivariate_normal(mean = [0]*n, cov = null_D, size = nsim));
    eta_df = null_df.multiply(se, axis = 1)
    transformed_df = eta_df.apply(func = lambda x: sqrt_U_inv.dot(x), axis = 1, raw = True)

    res_out = parallelize(transformed_df, estimate_statistics, cores, partitions, n, w_pos, t_v_pos )
    p_functions = cof_estimation(isf);
    res_out['null_p'] = res_out.loc[:,'null_stat'].apply(lambda x: pvalue_estimation(x, p_functions));
    
    Nbin = 1000
    bin_average = np.array(nsim/Nbin, dtype = np.float)
    
    def find_num(p,res):
        inds = np.floor(p * 1000)
        for ind in inds:
            res[int(ind)-1] += 1
        return(res)
    
    res = np.array([0] * Nbin)
    p = np.array(res_out.null_p.values, dtype = np.float)
    res = find_num(p, res)
    
    bins = pd.DataFrame(index = [i for i in range(Nbin)], columns = ['start','end'])
    bins.start = bins.index / Nbin
    bins.end = (bins.index + 1) / Nbin
    bins.loc[:,'num'] = res
    bins.loc[:,'above_thres'] = bins.num > (bin_average * 0.1)
    
    for i in range(Nbin-2, 0, -1):
        if bins.above_thres[i]:
            target_i = i + 1
            ind = (p > bins.start[target_i]) & (p <= bins.end[target_i])
            target_val = max(p[ind])
            break
    
    random_unif_min = target_val
    random_unif_max = 1
    ind = summary.pleio_p > target_val
    summary.loc[ind, 'pleio_p'] = np.array(np.random.uniform(low = random_unif_min, high = random_unif_max, size = sum(ind)),dtype=np.dtype(Decimal))
        
    return(summary)
