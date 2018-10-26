#!/usr/bin/env python3
'''
REG: Random Effect General Method
Copyright(C) 2018 Cue Hyunkyu Lee 

REG is a command line tool for performing cross-disease meta-analysis
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

from pval_estim.estim import regPestim, hetPestim, pfun_estim
from importance_sam.importance_sampling_REG import run_importance_sampling as runis
from meta_code.regeneral import REG_optim
from meta_code.LS import LS_chi
from basic_tool.etc import *
import numpy as np
import os, sys, traceback, argparse, time
from decimal import *

__version__ = 'beta 1.0.0'
MASTHEAD = "************************************************************\n"
MASTHEAD += "* General Random Effect Model (REG)\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C)2018 Cue H. Lee\n"
MASTHEAD += "* Seoul National University\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "************************************************************\n"

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    intervals = (('d', 86400), ('h', 3600), ('m', 60), ('s', 1), )
    f = ''

    for n, c in intervals:
        v = t // c
        if v:
            t -= v * c
            if c !=1:
                f += '{}{} '.format(round(v), n)
    return (f)


class Logger(object):
    '''
    Lightweight logging.
    TODO: replace with logging module

    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'w');

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.

        '''
        print (msg, file = self.log_fh);
        print (msg);


def onis(N,Rg,Re,isf):
    runis(N= int(N), GenCor=Rg, RECor=Re, outfn = os.path.dirname(os.path.abspath(__file__))+"/"+isf);

parser = argparse.ArgumentParser()
parser.add_argument('--out', default='regen',type=str,
    help='Output filename prefix. If --out is not set, REG will use regen as the '
    'default output filename prefix.')
# Besic Random Effect Meta-analysis Flags
parser.add_argument('--input', default=None, type=str,
    help='input file: summary statistics data')
parser.add_argument('--cor', default=None, type=str,
    help='Prefix for genetic and environment correlation matrix files '
    '[Prefix].rg and [Prefix].re. ')
parser.add_argument('--isf', default=None, type=str,
    help='Filename Prefix for Estimated null-distribution cumulative densities. ')
parser.add_argument('--create', default = False, action='store_true',
    help='If this flag is set, REG will run importance sampling and create new isf. ')
parser.add_argument('--nis', default = 1e6, type = float,
    help='Number of samples for importance sampling. ')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.out is None:
        raise ValueError('--out is required');
    
    log = Logger(args.out+'.log');
    try:
        defaults = vars(parser.parse_args(''));
        opts = vars(args);
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]];
        header = MASTHEAD;
        header += 'Call: \n';
        header += './reg.py \\\n';
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults];
        header += '\n'.join(options).replace('True','').replace('Flase','');
        header = header[0:-1]+'\n';
        log.log(header);
        log.log('Beginning analysis at {T}'.format(T=time.ctime()));
        start_time = time.time()

        if args.input is None:
            raise ValueError('Must specify --input with data.');
        if args.cor is None:
            log.log('No inputs for --cor: Rg and Re are diagonal of ones.');
        if args.isf is None:
            log.log('No inputs for --isf: Use {}.isf instead.'.format(args.out));
        
        isf = os.path.basename(args.isf)+".is";
        Rg = load_mat(args.cor+'.rg');Re = load_mat(args.cor+'.rg');
        
        if args.create:
            log.log('Importance sampling will generate new isf with {} samples. '.format(args.nis));
            onis(N = args.nis, Rg = Rg, Re = Re, isf = isf);
            
        with open(args.input,"r") as fin:
            i = 0;
            for line in fin:
                i = i + 1;
        ninput = i;
        del i

        reg_itck, reg_etck, het_mtck, het_itck, het_etck, tck_lim = pfun_estim(os.path.dirname(os.path.abspath(__file__))+"/"+isf)
        
        with open(args.input,"r") as fin, open(args.out,"w") as fout:
            Sreg = [0.0]*ninput;
            Preg = [Decimal(0)]*ninput;
            Shet = [0.0]*ninput;
            nShet = [0]*ninput;
            Phet = [Decimal(0)]*ninput;
            i = 0;
            for line in fin:
                cin = line.strip().split();
                n = int(len(cin[1:])/2);
                cvar = cin[0];
                css = cin[1:];
                cbeta = np.array([float(css[i*2]) for i in range(n)]);
                cstder = np.array([float(css[i*2+1]) for i in range(n)]);
                Sreg[i], Shet[i] = REG_optim(beta=cbeta, stder=cstder, Rg=Rg, Re=Re, n=n)
                Preg[i] = regPestim(cstat=Sreg[i], inter_tck=reg_itck, extra_tck=reg_etck, tck_lim=tck_lim)
                Phet[i] = hetPestim(cstat=Shet[i], mtck = het_mtck, inter_tck=het_itck, extra_tck=het_etck, tck_lim=tck_lim)
                print(" ".join(map(str,[cvar, Sreg[i], Preg[i], Shet[i], Phet[i]])),file = fout);
                i = i + 1;

    except Exception:
        ex_type, ex, tb = sys.exc_info();
        log.log(traceback.format_exc(ex));
        raise
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))

