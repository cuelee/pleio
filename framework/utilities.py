import numpy as np

### Code that can compute the root square pseudo inverse. Created by modifying numpy.linalg.pinv
def _makearray(a):
    new = np.asarray(a)
    wrap = getattr(a, "__array_prepare__", new.__array_wrap__)
    return new, wrap

def _is_empty_2d(arr):
    # check size first for efficiency
    return arr.size == 0 and np.product(arr.shape[-2:]) == 0

def sqrt_ginv(a, rcond=1e-15):
    '''
    This can compute the root square pseudo inverse.
    Dependencies: _is_empty_2d, _makearray
    '''
    a, wrap = _makearray(a)
    rcond = np.asarray(rcond)
    if _is_empty_2d(a):
        m, n = a.shape[-2:]
        res = np.empty(a.shape[:-2] + (n, m), dtype=a.dtype)
        return wrap(res)
    a = a.conjugate()
    u, s, vt = np.linalg.svd(a, full_matrices=False, hermitian=True)

    # discard small singular values
    cutoff = rcond[..., np.newaxis] * np.amax(s, axis=-1, keepdims=True)
    large = s > cutoff
    s = np.divide(1, s**0.5, where=large, out=s)
    s[~large] = 0

    res = np.matmul(np.transpose(vt), np.multiply(s[..., np.newaxis], np.transpose(u)))
    return wrap(res)


def is_pos_def(x):
    '''
    This function checks whether the input matrix is positive semi-definite (PSD).
    '''
    eigs = np.linalg.eigvals(x)
    return(np.all(eigs > 0))

