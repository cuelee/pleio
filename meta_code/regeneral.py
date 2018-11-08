## Dependancies : 
import math
import numpy as np
from scipy.stats import multivariate_normal 
from scipy.optimize import minimize

## U = GenCor[nonNA_array];
## R = RECor[nonNA_array];

def likelihood_mvn_tp (pars, x, Sg, sigma_e, n):
    mu = [0]*n; ## optim pars 
    tau = pars[0];
    sigma_g = tau * Sg;
    sigma_ge = sigma_g + sigma_e;
    k = multivariate_normal.logpdf(x = x, mean = mu, cov = sigma_ge);
    return(-1 * k);

def REG_optim (beta, stder, Sg, Re, n, bnds = [(0,200)]):
    sigma_e = np.diag(stder).dot(Re).dot(np.diag(stder));
    res = minimize(likelihood_mvn_tp, x0 = np.random.uniform(0.001,.2,1), method ='SLSQP', bounds = bnds, args=(beta, Sg, sigma_e, n), options = {'ftol' : 1e-9, 'disp':False});
    if (res.success != True):
        print(res.success)
    ta = res.x[0];
    alt_var = -1 * res.fun;
    nul_var = multivariate_normal.logpdf( x = beta, mean = [0]*n, cov = sigma_e);
    return( max(2 * (alt_var - nul_var),0.0) );

def REG_apply(Sg, Re, n, X, row_wise):
    stder = [1]*n;
    if(row_wise == True):
        niter = X.shape[0];
    else:
        niter = X.shape[1];
    reg_vec = [0.0]*niter;
    if(row_wise == True):
        for i in range(niter):
            beta = X[i,:]; reg_vec[i] = REG_optim(beta = beta, stder = stder, Sg=Sg, Re=Re, n=n);
    else:
        for i in range(niter):
            beta = X[:,i]; reg_vec[i] = REG_optim(beta = beta, stder = stder, Sg=Sg, Re=Re, n=n);
    return(reg_vec)
