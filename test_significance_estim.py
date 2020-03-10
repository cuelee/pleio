import pandas as pd
from numpy.random import uniform
from framework.significance_estimation import cof_estimation, pvalue_estimation

isof = '/Users/cuelee/Dropbox/github/delpy/isf.isf'
cof = cof_estimation(isof)
for i in uniform(0, 10, size = 1000000):
    print(pvalue_estimation(float(i), cof))
