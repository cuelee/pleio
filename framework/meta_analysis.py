from pval_estim.estim import regPestim
from meta_code.regeneral import REG_optim
from framework.parse import *
import pandas as pd

def meta_analysis(df, isres, args, log):
    '''
    The function integrates multiple summary statistics and provide data frames of pvalues and pooled log OR 
    df is a python class object which contains 1. summary statistics(dataframe), 
    2. Genetic Covariance matrix(numpy matrix) 3. Random errors(numpy matrix) 
    isres is a python class object which contains tabulated importance sampling results to estimate pvalues  
    '''
    Sg=df.Sg; Re=df.Re; dat = df.metain; 
    mtck, itck, etck, tck_lim = isres.mtck, isres.itck, isres.etck, isres.tck_lim

    res = pd.DataFrame(dat[['SNP']])
    
    def rowwise_optim(x, se, Sg, Re, n):
        s = REG_optim(beta=x, stder=se, Sg=Sg, Re=Re, n=n)
        return(s)
    
    def rowwise_pestim(s,mtck,itck,etck,tck_lim):
        p = regPestim(cstat=s, mtck=mtck, inter_tck=itck, extra_tck=etck, tck_lim=tck_lim)
        return(p)

    n = Sg.shape[1];
    se=[1]*n
    res['stat'] = dat.apply(lambda x: rowwise_optim(x.tolist()[1:], se, Sg, Re, n), axis=1)
    res['Pvalue'] = res.apply(lambda x: rowwise_pestim(x['stat'], mtck, itck, etck, tck_lim), axis=1)

    return(res)

