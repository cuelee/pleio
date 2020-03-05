## Dependancies : 
import numpy as np
from scipy.stats import multivariate_normal 
from scipy.optimize import fminbound

## T is the number of traits
## Sg is the Genetic covariance matrix with the size of [T x T];
## Se is the covariance matrix of error with the size of [T x T];

def LL_fun(x,n,sq,w):
    return(-0.5*(n*np.log(2*np.pi)+sum(np.log(w+x))+sum(sq/(w+x))));

def LLp_fun(x,sq,w):
    return(0.5*(sum(1/(w+x))-sum(sq/(w+x)**2)));

def LLdp_fun(x,sq, w):
    return(-0.5*(sum(1/(w+x)**2)-2*sum(sq/(w+x)**3)));

def NR_root(f, df, x, sq, w, i = 0, iter_max = 10000, tol = 2.22044604925e-16**0.5):
    '''
        f: LLp_fun
        df: LLdp_fun
    '''
    while ( abs(f(x,sq,w)) > tol ):
        x = x - f(x,sq,w) / df(x,sq,w);
        i = i + 1;
        if (i == iter_max): 
            break;
    return(x)

def vcm_optimization (b, K, tol = 2.22044604925e-16**0.5):
    '''
    
    inputs:
        b : betas, type: numpy array object, dtype: float 
        K : covariance matrix, type: numpy matrix object, dtype: float
    
    outputs:
        stats : PLEIO LRT statistic
        mle.tausq : maximum likelihood estimator (variance estimate)
    
    parameters description
        crossP : cross product of betas and transposed eigen vectors of K
        w : Eigen values of K
        v : right eigen vectors of K
        n : number of studies
    
    '''
    n = len(b);
    l, c = np.linalg.eigh(K);
    p = l > tol 
    if(all(p)):
        crossP = np.transpose(c).dot(b);
        sq = crossP**2;
        w = l;
    else:
        cr = np.transpose(c[:,p]).dot(b);
        sq = cr**2;
        w = l[p];
    t = [10**(i/4) for i in range(-36,24,1)];
    init = t[np.argmax([LL_fun(i, n, sq, w) for i in t])];
    tausq = NR_root(LLp_fun, LLdp_fun, init, sq, w);
    if (tausq <0): tausq = 0;
    null_ll = LL_fun(0, n, sq, w);
    alt_ll = LL_fun(tausq, n, sq, w) ;
    if(alt_ll < null_ll): tausq = 0; alt_ll = null_ll;
    return (- 2 * (null_ll - alt_ll))

def sqrt_ginv (X, tol = 2.22044604925e-16**0.5):
    u,s,vh = np.linalg.svd(X)
    Positive = s > max(tol * s[0], 0) 
    if (all(Positive)):
        res=np.transpose(vh).dot(np.diag(1/s)**0.5).dot(np.transpose(u))
    elif (not any(Positive)): 
        res=np.array([0]*np.prod(X.shape)).reshape(X.shape)
    else:
        res=vh[:, Positive].dot(np.diag(1/s[Positive])**0.5).dot(np.transpose(u[:, Positive]))
    return(res)
