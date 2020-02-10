import numpy as np
from numpy.linalg import svd
import pandas as pd
from scipy.stats import multivariate_normal as mvn

x=[1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1]
mean = [0]*18
rg_fn = '/Users/cuelee/Dropbox/PLEIO/Meeting/2019_07_24_RealDA/rg.txt.gz'
rg = pd.read_csv(rg_fn, sep='\t', compression = 'gzip')
rg.index = rg.columns
cov = rg.values

#a = mvn.logpdf(x=x, mean = mean, cov = cov)


def _makearray(a):
    new = np.core.asarray(a)
    wrap = getattr(a, "__array_prepare__", new.__array_wrap__)
    return new, wrap

def _isEmpty2d(arr):
    # check size first for efficiency
    return arr.size == 0 and product(arr.shape[-2:]) == 0



def pinv(a, rcond=1e-15, hermitian=False):
    """
    Compute the (Moore-Penrose) pseudo-inverse of a matrix.
    Calculate the generalized inverse of a matrix using its
    singular-value decomposition (SVD) and including all
    *large* singular values.

    Date at:2020.02.06
    Cue fixed numpy.linalg.pinv to generate a product of components in s

    """
    a, wrap = _makearray(a)
    rcond = np.core.asarray(rcond)
    if _isEmpty2d(a):
        m, n = a.shape[-2:]
        res = empty(a.shape[:-2] + (n, m), dtype=a.dtype)
        return wrap(res)
    a = a.conjugate()
    u, s, vt = svd(a, full_matrices=False, hermitian=hermitian)

    # discard small singular values
    cutoff = rcond[..., np.core.newaxis] * np.core.amax(s, axis=-1, keepdims=True)
    large = s > cutoff
    s = np.core.divide(1, s, where=large, out=s)
    s[~large] = 0

    res = np.matmul(np.core.transpose(vt), np.core.multiply(s[..., np.core.newaxis], np.core.transpose(u)))
    return wrap(res), 


def mvn_logpdf(x, mean, cov, tol = np.finfo(float).eps):
    if( len(cov.shape) > 2 or not np.issubdtype(cov.dtype, np.number) ):
        raise("`cov` must be a numeric or complex matrix")
    cov_inv = pinv(cov)
    return(cov_inv)

print(mvn_logpdf(x=x, mean=mean, cov=cov))
