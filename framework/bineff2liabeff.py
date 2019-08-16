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


def search_cont_af(x, OR, V, Ncase, Ncont):
    cont_af = x;
    case_af = ( OR * cont_af ) / ( (OR - 1) * cont_af +1 );

    a = Ncase * 2 * (case_af * (1.0 - case_af) + case_af**2);
    b = Ncont * 2 * (cont_af * (1.0 - cont_af) + cont_af**2);
    c = Ncase * 2 * (1.0 - case_af) * (case_af + (1.0 - case_af));
    d = Ncont * 2 * (1.0 - cont_af) * (cont_af + (1.0 - cont_af));
    
    a,b,c,d = map(max, [ [i,1e-20] for i in [a,b,c,d] ])

    estim_V = 1/a + 1/b + 1/c + 1/d;
    return(V - estim_V)


def check_variance(af1, OR, N2, N1):
    af2 = ( OR * af1 ) / ( (OR - 1) * af1 +1 );

    a = N2 * 2 * (af2 * (1.0 - af2) + af2**2);
    b = N1 * 2 * (af1 * (1.0 - af1) + af1**2);
    c = N2 * 2 * (1.0 - af2) * (af2 + (1.0 - af2));
    d = N1 * 2 * (1.0 - af1) * (af1 + (1.0 - af1));
    
    a,b,c,d = map(max, [ [i,1e-20] for i in [a,b,c,d] ])

    V = 1/a + 1/b + 1/c + 1/d;
    return(V)

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
    if(se == 0 or Ncase == 0 or Ncont == 0):
        print(signed_val, se, Ncase, Ncont, Fp, signed_name, SNP)
        return(np.nan, np.nan)
    N = Ncase + Ncont;
    # Fcc: case control sample prevalence
    Fcc = Ncase/N;
    
    lib_case, lib_cont = generate_liability_posterior_mean(Fp);
    
    V_min = check_variance(0.5, OR, Ncase, Ncont)
    V_max = check_variance(0.00000001, OR, Ncase, Ncont)
    if V_min > V:
        cont_af = 0.5
        case_af = OR * cont_af / (( OR - 1 ) * cont_af + 1);
    elif V_max < V: 
        print(signed_val, se, Ncase, Ncont, Fp, signed_name,SNP);
        return(np.nan, np.nan)
    else:
        cont_af = root(f = search_cont_af, a = 0.00000001, b = 0.5, args = (OR, V, Ncase, Ncont));
        case_af = OR * cont_af / (( OR - 1 ) * cont_af + 1);

    if(cont_af > 1 or case_af > 1 or cont_af < 0 or case_af < 0):
        print(signed_val, se, Ncase, Ncont, Fp, signed_name,SNP, cont_af, case_af);
        raise ValueError('Convergence Error occurred: bineff2liabeff')
    print(OR, se, V, Ncase, Ncont, Fp, cont_af, case_af)

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
  
    if(mean_squared_simul_err < 0):
            print(signed_val, se, Ncase, Ncont, Fp, signed_name,SNP, i);
            raise ValueError('Found negative mean_squared_simul_err')
    se_estim = np.sqrt(mean_squared_simul_err * N / (N - 2) / ((Ncase - 1) * var_std_case_X  + (Ncont - 1) * var_std_cont_X +  Ncase * Ncont / N * (mean_std_case_X - mean_std_cont_X)**2 ))
    return(beta_estim, se_estim)

#SNP         CHR BP          A1  A2  Z           P           N       BETA           SE          N_CASE          N_CONTROL       INFO
#rs16945492  18  3839896     A   G   -0.35561    0.722133    361194  -8.05873e-05   0.000226617 1405.00371458   359788.996285   0.996052
#rs111726894 18  71796714    G   T   -0.168595   0.866115    361194  -0.333912      1.98056     1405.0002171    359788.999783   0.902917
#print(bineff2liabeff(signed_val=0.2129459, se=0.003527057, Ncase=200000, Ncont=1500000, Fp=0.03, signed_name = 'BETA', SNP = 'A'))
