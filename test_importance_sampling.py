from framework.importance_sampling import importance_sampling as run_is
from framework.significance_estimation import pfun_estim, pvalue_estimation
import time
import numpy as np
import random

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    intervals = (('d', 86400), ('h', 3600), ('m', 60), ('s', 1), )
    f = ''
    for n, c in intervals:
        v = t // c
        if v:
            t -= v * c
            if t !=1:
                f += '{}{} '.format(round(v), n)
    return (f)

print('Beginning analysis at {T}'.format(T=time.ctime()));
start_time = time.time()

n = 7
rg = np.array([[0.9]*n]*n) 
for i in range(n):
    rg[i,i] = 1.0 
re = np.diag([1.0] * n)
#h2 = np.array([ np.random.uniform(0,0.4,1)[0] for i in range(n)])
h2 = [1]*n
#print(h2)
sg = np.diag(np.sqrt(h2)).dot(rg).dot(np.diag(np.sqrt(h2)))
#print(rg, re, sg)
ccov = sg+re


run_is(N=1000, gwas_N =np.array(random.sample([1000,1000000,100,100000]*100,n)), U = sg, Ce = re, outf = './isf.isf', mp_cores = 1)

#iso=pfun_estim('./isf.isf')

#for i in range(100):
#    print(pvalue_estimation(i,iso))

time_elapsed = round(time.time()-start_time,2)
print('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
