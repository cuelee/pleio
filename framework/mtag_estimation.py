import pandas as pd
import numpy as np

def mtag(y, G, R, c):
    P = len(y)
    ind = [val+a for val in c for a in ('_beta', '_se')]
    res = pd.Series(np.zeros(P*2), index = ind)
    for p in range(P):
        omega_t = G[:,p]
        omega_tt = G[p,p]
        om_min_gam = G - np.outer(omega_t,omega_t)/omega_tt

        inv_mid = np.linalg.inv(om_min_gam + R)
        tdtt = omega_t/omega_tt
        top = np.transpose(tdtt).dot(inv_mid).dot(y)
        bot = np.transpose(tdtt).dot(inv_mid).dot(tdtt)
        res[c[p]+'_beta'] = top/bot 
        res[c[p]+'_se'] = 1/bot
    return(res)

