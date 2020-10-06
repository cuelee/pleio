import numpy as np
import pandas as pd

def blup_optimization(y, G, R, c):
    pinv_R = np.linalg.pinv(R)
    pinv_G = np.linalg.pinv(G)
    V = np.linalg.pinv(pinv_R + pinv_G)
    u = V.dot(pinv_R).dot(y)
    if(no.count_nonzero(pinv_G) > 0):
        u_se = np.sqrt(np.diag(V))
        return(pd.Series([item for sublist in zip(u, u_se) for item in sublist], index = [item for sublist in zip([a + '_beta' for a in c],[a + '_se' for a in c]) for item in sublist]) )
    else
        return(pd.Series([item for sublist in zip(u,[-9]*len(u)) for item in sublist], index = [item for sublist in zip([a + '_beta' for a in c],[a + '_se' for a in c]) for item in sublist]))

