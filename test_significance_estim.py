import pandas as pd
from numpy.random import uniform
from framework.significance_estimation import cof_estimation, pvalue_estimation

isof = '/Users/cuelee/Dropbox/github/delpy/isf.isf'
cof = cof_estimation(isof)
print(pvalue_estimation(7.5, cof))
