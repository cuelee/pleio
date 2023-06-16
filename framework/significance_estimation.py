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
    df_out['null_stat'] = df_data.apply(lambda x: vcm_optimization(x.tolist(), n, w, t_v))
    return(df_out)

def parallelize(df_input, func, cores, partitions, n, w, t_v):
    data_split = np.array_split(df_input, partitions)
    iterable = product(data_split, [n], [w], [t_v])
    pool = mp.Pool(int(cores))
    df_output = pd.concat(pool.starmap(func, iterable))
    pool.close()
    pool.join()
    return(df_output)

def quantile_mapping(s, N_bin = 100, lower_bound = 0.3):
    INDEX = s.index
    s = s.to_numpy()
    # Define the bin edges
    bin_edges = np.linspace(lower_bound, 1, N_bin+1)

    # Get the counts in each bin for values in s that are between 0.6 and 1 (exclusive)
    hist, _ = np.histogram(s[(s > lower_bound) & (s <= 1)], bins=bin_edges)

    # The target number of values in each bin is the total count divided by the number of bins
    target_bin_count = len(s[s > lower_bound]) // N_bin
    count_ones = np.sum(s == 1)
    
    # Now we will replace the 1's in s with values that are distributed evenly across the bins
    for i in range(N_bin):
        bin_count = hist[i]
        values_to_add = target_bin_count - bin_count
        values_to_add = np.min([values_to_add, count_ones])
        if values_to_add > 0:
            # Generate 'values_to_add' random values within this bin
            bin_start = bin_edges[i]
            bin_end = bin_edges[i+1]
            new_values = np.random.uniform(bin_start, bin_end, size=values_to_add)

            # Find the first 'values_to_add' 1's in s and replace them with the new values
            ones_indices = np.flatnonzero(s == 1)
            indices_to_replace = ones_indices[:values_to_add]
            s[indices_to_replace] = new_values
            count_ones = len(ones_indices)
    return(pd.Series(s, index = INDEX, name = 'pleio_p'))

def flattening_p_value(summary):
    summary = summary.copy()  # create a copy of the DataFrame to avoid modifying the original data
    s = summary['pleio_p']
    summary.loc[:, 'pleio_p'] = quantile_mapping(s)

    return summary
