## Dependancies : 
import math
import numpy as np
from scipy.stats import multivariate_normal 
from scipy.optimize import minimize

## U = GenCor[nonNA_array];
## R = RECor[nonNA_array];

def likelihood_mvn_tp (pars, x, Rg, sigma_e, n):
    mu = pars[0]; ## optim pars 
    tau = pars[1];
    sigma_g = tau * Rg;
    sigma_ge = sigma_g + sigma_e;
    k = multivariate_normal.logpdf(x = x, mean = [mu]*n, cov = sigma_ge);
    return(-1 * k);

def REG_optim (beta, stder, Rg, Re, n, bnds = ((-200,200),(0,200))):
    sigma_e = np.diag(stder).dot(Re).dot(np.diag(stder));
    res = minimize(likelihood_mvn_tp, x0 = np.random.uniform(0.001,.2,2), method ='SLSQP', bounds = bnds, args=(beta, Rg, sigma_e, n), options = {'ftol' : 1e-9, 'disp':False});
    if (res.success != True):
        print(res.success)
    mu, ta = res.x;
    alt_var = -1 * res.fun;
    tau_var = multivariate_normal.logpdf( x = beta, mean = [mu]*n, cov = sigma_e);
    nul_var = multivariate_normal.logpdf( x = beta, mean = [0]*n, cov = sigma_e);
    return( [ max(2 * (alt_var - nul_var),0.0), max(2 * (alt_var - tau_var),0.0) ] );

