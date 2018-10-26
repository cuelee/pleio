from scipy import interpolate
from decimal import *
from pval_estim.readf import readf

def estim_interPfun(x,y):
    y[0]=1;
    return(interpolate.splrep(x, y, s=0));

def estim_extraPfun(x,y):
    diff = [ y[i+1] /  y[i] for i in range(len(y)-1) ];
    ratio = round( 1 / len(diff[-10:]) * sum(diff[-10:]),100 );
    newx = [ x[i] for i in range(len(x)) ];
    newy = [ Decimal(y[i]) for i in range(len(y)) ];
    irange = [ i + len(x) for i in range(10001-len(x)) ];
    for i in irange:
        newx.append(i);
        newy.append(Decimal(newy[i-1]*Decimal(ratio)));
    return([newx,newy]);    

def extrapolate_cue(x, extra_tck):
    res = estim_stat_extralinear(x,extra_tck);
    return(res);

def estim_stat_extralinear(aval, extra_tck):
    tck_x = extra_tck[0];
    tck_y = extra_tck[1];
    if(aval > 10000 ):
        aval = Decimal(10000);
    ivalue = int(round(aval));
    fvalue = Decimal(aval - ivalue);
    if(fvalue > 0):
        diff = Decimal(tck_y[ivalue+1] - tck_y[ivalue]);
        res = tck_y[ivalue] + Decimal(diff * fvalue);
    else:
        res = tck_y[ivalue];
    return(res);

def regPestim(cstat, inter_tck, extra_tck, tck_lim):
    if (cstat <= tck_lim):
        res = Decimal( str( interpolate.splev(cstat, inter_tck, der=0)));
    elif(cstat > tck_lim):
        res = extrapolate_cue(cstat, extra_tck);
    return(res);

def manPestim(cstat, mtck):
    x=mtck[0];y=mtck[1];
    c=float(cstat);
    ci=round(c); cf=c-ci;
    diff = y[ci+1] - y[ci];
    return( y[ci] + diff * cf );

def hetPestim(cstat, mtck, inter_tck, extra_tck, tck_lim):
    if (cstat <= 0.1):
        res = Decimal( manPestim(cstat, mtck) );
    elif (cstat > 0.1 and cstat <= tck_lim):
        res = Decimal( str( interpolate.splev(cstat, inter_tck, der=0)));
    elif (cstat > tck_lim):
        res = extrapolate_cue( cstat, extra_tck );
    return(res);

def pfun_estim(isf):
    x,reg,het = readf(isf);

    x1 = [ x[i] for i in range(len(x)) if x[i] <= 0.1 ];
    het1 = [ het[i] for i in range(len(x)) if x[i] <= 0.1 ];

    x2 = [ x[i] for i in range(len(x)) if x[i] > 0.1 and x[i] <= 30 ];
    het2 = [ het[i] for i in range(len(x)) if x[i] > 0.1 and x[i] <= 30 ];
    reg_itck = estim_interPfun(x,reg);
    reg_etck = estim_extraPfun(x,reg);

    het_mtck = [x1,het1];
    het_itck = estim_interPfun(x2,het2);
    het_etck = estim_extraPfun(x2,het2);
    return(reg_itck, reg_etck, het_mtck, het_itck, het_etck, x[-1]);
