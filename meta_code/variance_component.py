## Dependancies : 
import numpy as np
from scipy.stats import multivariate_normal 
from scipy.optimize import minimize

## Sg = GenCor[nonNA_array];
## Se = RECor[nonNA_array];

def likelihood_mvn_tp (pars, x, Sg, Sn, n):
    tau = pars[0];
    vcg = tau * Sg;
    vc = vcg + Sn;
    k = multivariate_normal.logpdf(x = x, mean = [0] * n, cov = vc);
    return(-1 * k);

def vc_optimization (b, se, Sg, Rn, n, bnds = [(0,200)]):
    Sn = np.diag(se).dot(Rn).dot(np.diag(se));
    res = minimize(likelihood_mvn_tp, x0 = np.random.uniform(0.001,.2,1), method ='SLSQP', bounds = bnds, args=(b, Sg, Sn, n), options = {'ftol' : 1e-9, 'disp':False});
    if (res.success != True):
        print(res.success)
    tau = res.x[0];
    alt_var = -1 * res.fun;
    nul_var = multivariate_normal.logpdf( x = b, mean = [0]*n, cov = Sn);
    return( max(2 * (alt_var - nul_var), 0.0) );

def delpy_apply(Sg, Rn, n, X, row_wise):
    stder = [1]*n;
    if(row_wise == True):
        niter = X.shape[0];
    else:
        niter = X.shape[1];
    delpy_stats = [0.0]*niter;
    if(row_wise == True):
        for i in range(niter):
            b = X[i,:]; 
            delpy_stats[i] = vc_optimization(b = b, se = se, Sg = Sg, Rn = Rn, n = n);
    else:
        for i in range(niter):
            b = X[:,i]; 
            delpy_stats[i] = vc_optimization(b = b, se = se, Sg = Sg, Rn = Rn, n = n);
    return(float(delpy_stats))
