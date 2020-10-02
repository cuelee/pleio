import numpy as np
import pandas as pd

def blup_optimization(y, G, R, c):
    pinv_R = np.linalg.pinv(R)
    pinv_G = np.linalg.pinv(G)
    return(pd.Series(np.linalg.pinv(pinv_R + pinv_G).dot(pinv_R).dot(y), index = c))
