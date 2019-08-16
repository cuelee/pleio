from scipy import interpolate
from decimal import *

def readf(afile):
    with open(afile,'r') as fin:
        n = 0;
        for line in fin:
            n = n + 1;

    with open(afile,'r') as fin:
        x = [0] * n;
        reg = [0] * n;
        i = 0;
        for line in fin:
            std = line.strip().split();
            x[i] = float(std[0]);
            reg[i] = float(std[1]);
            i = i + 1;
    return(x, reg)

def estim_interPfun(x,y):
    return(interpolate.splrep(x, y, s=0));

def estim_extraPfun(x,y,interval=10):
    dtab = dict()
    for i in range(len(x)):
        dtab[x[i]]=Decimal(y[i]);
    top=max(list(dtab.keys()))
    ran = [ top-interval+i for i in range(interval) ]
    diff = [ dtab[i+1] / dtab[i] for i in ran]
    ratio = Decimal(1 / len(ran)) * sum(diff);
    irange = [ i + 1 + int(top) for i in range(10001-int(top)) ];
    for i in irange:
        dtab[i]=Decimal(dtab[i-1])*Decimal(ratio)
    return(dtab);

def extrapolate_cue(x, dtab):
    res = estim_stat_extralinear(x,dtab);
    return(res);

def estim_stat_extralinear(aval, dtab):
    if(aval > 10000 ):
        aval = Decimal(10000);
    ivalue = int(round(aval));
    fvalue = Decimal(aval - ivalue);
    if(fvalue > 0):
        diff = Decimal(dtab[ivalue+1] - dtab[ivalue]);
        res = dtab[ivalue] + Decimal(diff * fvalue);
    else:
        res = dtab[ivalue];
    return(res);

#def regPestim(cstat, inter_tck, extra_tck, tck_lim):
#    if (cstat <= tck_lim):
#        res = Decimal( str( interpolate.splev(cstat, inter_tck, der=0)));
#    elif(cstat > tck_lim):
#        res = extrapolate_cue(cstat, extra_tck);
#    return(res);

def manPestim(cstat, mtck):
    x=mtck[0];y=mtck[1];
    c=float(cstat);
    ci=round(c); cf=c-ci;
    diff = y[ci+1] - y[ci];
    return( y[ci] + diff * cf );

def pvalue_estimation(cstat, iso):
    if (cstat <= 0.1):
        res = Decimal( manPestim(cstat, iso.mtck) );
    elif (cstat > 0.1 and cstat <= iso.tck_lim):
        res = Decimal( str( interpolate.splev(cstat, iso.itck, der=0)));
    elif (cstat > iso.tck_lim):
        res = extrapolate_cue( cstat, iso.dtab );
    return(res);

class pfun_estim(): 
    def __init__(self,isf):
        x,reg = readf(isf);
        self.tck_lim = x[-1];
        reg[0]=1;
        x1 = [ x[i] for i in range(len(x)) if x[i] <= 0.1 ];
        reg1 = [ reg[i] for i in range(len(x)) if x[i] <= 0.1 ];
        x2 = [ x[i] for i in range(len(x)) if x[i] >= 0.1 and x[i] <= self.tck_lim ];
        reg2 = [ reg[i] for i in range(len(x)) if x[i] >= 0.1 and x[i] <= self.tck_lim ];
        self.mtck = [x1,reg1];
        self.itck = estim_interPfun(x2,reg2);
        self.dtab = estim_extraPfun(x2,reg2);
