from scipy.stats import norm
from scipy.integrate import quad
from scipy.optimize import brentq as root
from scipy.optimize import newton as root2
import numpy as np
import time 

#Example
#OR= 1.237049 ;logOR=np.log(OR);V= 1.244427e-05 ;se=np.sqrt(V); Ncase =  2e+05 ; Ncont =  1500000; Fp = 0.03
#print(bineff2liabeff(V, Ncase, Ncont, 0.03, OR))

def posterior_nom(x):
    return(x*norm.pdf(x, loc = 0, scale = 1));


def posterior_denom(x):
    return(norm.pdf(x, loc = 0, scale = 1));


def generate_liability_posterior_mean(F):
    m = norm.ppf(q = (1-F), loc = 0, scale = 1);
    
    # estimate residual posterior probability of cases and controls 
    case_rpp = quad(func = posterior_nom, a = m, b = np.inf)[0] / quad(func = posterior_denom, a = m, b = np.inf)[0];
    cont_rpp = quad(func = posterior_nom, a = - np.inf, b = m)[0] / quad(func = posterior_denom, a = - np.inf, b = m)[0];
    
    # estimate mean liability of cases and controls 
    ml_case = - m + case_rpp;
    ml_cont = - m + cont_rpp;
    return(ml_case, ml_cont)


def search_cont_af(x, Ncase, Ncont, OR, V):
    cont_af = x;
    case_af = ( OR * cont_af ) / ( (OR - 1) * cont_af +1 );
    
    a = Ncase * (case_af * (1.0 - case_af) + case_af**2);
    b = Ncont * (cont_af * (1.0 - cont_af) + cont_af**2);
    c = Ncase * (1.0 - case_af) * (case_af + (1.0 - case_af));
    d = Ncont * (1.0 - cont_af) * (cont_af + (1.0 - cont_af));
    
    estim_V = 0.5 * (1/a + 1/b + 1/c + 1/d);
    return(V - estim_V)

def root_finding(Ncase, Ncont, OR, V, i):
    try:
        cont_af = root(f = search_cont_af, a = 0.001, b = 0.5, args = (Ncase, Ncont, OR, V*(1+i*0.01)));
        return(cont_af)
    except:
        return(-1)

## function binary effect estimates to liability scaled effect estimates 
## OR : odd ratio, var_OR : variance of OR, Ncase: number of cases, Ncont: number of controls, Fp: population prevalence 
def bineff2liabeff(signed_val, se, Ncase, Ncont, Fp, signed_name = None, SNP = None):
    if signed_name == None: 
        raise ValueError('no signed stats are assigned.')
    
    if signed_name == 'BETA':
        OR = np.exp(signed_val);
    elif signed_name =='OR':
        logOR = np.log(OR);

    V = se**2
    N = Ncase + Ncont;
    # Fcc: case control sample prevalence
    Fcc = Ncase/N;
    
    lib_case, lib_cont = generate_liability_posterior_mean(Fp);
    i=0; cont_af = -1;
    while(cont_af < 0 and i < 100):
        cont_af = root_finding(Ncase,Ncont,OR,V,i)
        i += 1;
    if i > 100 and root_finding(Ncase,Ncont,OR,V,i) < 0 or V > 1:
        print('SNP: {id} has non realistic OR: {Or} and SE: {se} values /nReturn NA'.format(id = SNP, Or = OR, se = se))
        return('NA', 'NA')
    case_af = OR * cont_af / (( OR - 1 ) * cont_af + 1);
    # smaf : sample maf
    sam_af = case_af * Fcc + cont_af * (1-Fcc);
    
    mean_std_case_X = (case_af - sam_af) / (sam_af * (1 - sam_af));
    mean_std_cont_X = (cont_af - sam_af) / (sam_af * (1 - sam_af));
    var_std_case_X = (2 * case_af * (1 - case_af)) / (2 * sam_af * (1 - sam_af))**2;
    var_std_cont_X = (2 * cont_af * (1 - cont_af)) / (2 * sam_af * (1 - sam_af))**2;
    
    cov_XY = ( mean_std_case_X * lib_case * Ncase + mean_std_cont_X * lib_cont * Ncont) / N - (mean_std_case_X * Ncase + mean_std_cont_X * Ncont) / N * (lib_case * Ncase + lib_cont * Ncont) / N;
    var_X = ((Ncase - 1) * var_std_case_X  + (Ncont - 1) * var_std_cont_X +  Ncase * Ncont / N * (mean_std_case_X - mean_std_cont_X)**2 ) / (N - 1);
    beta_estim = cov_XY / var_X;
    
    mean_Y = (lib_case * Ncase + lib_cont * Ncont) / N;
    mean_std_X = 0;
    
    alpha_estim = mean_Y - beta_estim * mean_std_X;
    mean_squared_simul_err = (      (lib_case - alpha_estim - beta_estim * (0 - 2 * sam_af) / (2 * sam_af * (1 - sam_af)))**2 * (1 - case_af)**2 * Fcc 
                                +   (lib_case - alpha_estim - beta_estim * (1 - 2 * sam_af) / (2 * sam_af * (1 - sam_af)))**2 * 2 * case_af * (1 - case_af) * Fcc
                                +   (lib_case - alpha_estim - beta_estim * (2 - 2 * sam_af) / (2 * sam_af * (1 - sam_af)))**2 * case_af**2 * Fcc
                                +   (lib_cont - alpha_estim - beta_estim * (0 - 2 * sam_af) / (2 * sam_af * (1 - sam_af)))**2 * (1 - cont_af)**2 * (1 - Fcc)
                                +   (lib_cont - alpha_estim - beta_estim * (1 - 2 * sam_af) / (2 * sam_af * (1 - sam_af)))**2 * cont_af * (1 - cont_af) * (1 - Fcc)
                                +   (lib_cont - alpha_estim - beta_estim * (2 - 2 * sam_af) / (2 * sam_af * (1 - sam_af)))**2 * cont_af**2 * (1 - Fcc)    )
    
    var_estim = np.sqrt(mean_squared_simul_err * N / (N - 2) / ((Ncase - 1) * var_std_case_X  + (Ncont - 1) * var_std_cont_X +  Ncase * Ncont / N * (mean_std_case_X - mean_std_cont_X)**2 ))
    return(beta_estim, var_estim)

