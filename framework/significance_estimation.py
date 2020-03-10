from scipy.interpolate import splev, splrep
from decimal import *
import pandas as pd
import numpy as np
from numbers import Number

def readf(f):
    d = pd.read_csv(f,names=['x','s'],sep=' ');
    d.x = d.x.round(8);
    d.s.iloc[0] = 1;
    return(d);

def manual_estimation(x1,y1,x0=0,y0=1):
    c = (y1-y0)/(x1-x0);
    mestim=lambda x: 1+x*c;
    return(mestim);

def interpolationf(d):
    tck = splrep(d.x.values,np.log(d.s.values), s = 0);
    return(tck);

def extrapolationf(d):
    x = d.x.values;
    y = np.log(d.s.values);
    b = np.cov(x,y)[1,0]/np.var(x);
    a = y[-1]-b*x[-1];
    mestim = lambda x: a + b * x;
    return(mestim);

def pvalue_estimation(s,iso):
    if(not isinstance(s,Number)):
        raise ValueError('The value of the input must be numeric')
    if (s <= iso.min):
        return(Decimal(iso.low(s)));
    elif (s <= iso.max):
        return(np.exp(Decimal(float(splev(s,iso.itck,der=0,ext=2)))));
    else:
        return(np.exp(Decimal(iso.tail(s))));

class cof_estimation(): 
    def __init__(self,isf):
        d = readf(isf);
        self.max = d.x.iloc[-1];
        self.min = d.x.iloc[1];
        self.low = manual_estimation(d.x.iloc[1], d.s.iloc[1]);
        self.itck = interpolationf(d.iloc[1:,:]);
        self.tail = extrapolationf(d.loc[d.x >= 20,:])
