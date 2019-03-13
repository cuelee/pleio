from framework.importance_sampling import importance_sampling as run_is
from pval_estim.estim import pfun_estim, pestim
import time
import numpy as np

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

cov = np.diag([1]*7)
mean = np.array([1]*7)

run_is(N=1000, Sg = cov, Re = cov, outfn = './isf.isf')
iso=pfun_estim('./isf.isf')

for i in range(100):
    print(pestim(i*0.01,iso))

time_elapsed = round(time.time()-start_time,2)
print('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
