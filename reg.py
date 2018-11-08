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

from framework.matrixIO import *
from framework.meta_analysis import *
from framework.parse import *
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

def check_xtractmat_output(cur_path):
    if os.path.isdir(os.path.join(cur_path, '_out_')):
        if (os.path.exists(os.path.join(cur_path,'_out_','sumstats.list')) and os.path.exists(os.path.join(cur_path,'_out_','sumstats.re')) and os.path.exists(os.path.join(cur_path,'_out_','sumstats.sg'))): return(1);
        else: return(0);
    else: return (0);

def verify_test(ft, fl):
    for fi in ft:
        if not fi in fl: raise ValueError('{} is not included in the summary statistics data.'.format(fi) );
    new_ft = [x for x in fl if x in ft]
    if len(new_ft) < 1: raise ValueError ('There is only one sumstats to test');
    return(new_ft)


parser = argparse.ArgumentParser()
parser.add_argument('--out', default='regen',type=str,
    help='Output filename prefix. If --out is not set, REG will use regen as the '
    'default output filename prefix.')
# Besic Random Effect Meta-analysis Flags
parser.add_argument('--sumstats', default=None, type=str,
    help='input folder path: path to summary statistics data.')
parser.add_argument('--test', default=None, type=str,
    help='comma separated sumstats filename: Use the same file name from ./{sumstats}/_out_/sumstats.list ')
parser.add_argument('--isf', default=None, type=str,
    help='Filename Prefix for Estimated null-distribution cumulative densities. ')
parser.add_argument('--create', default = False, action='store_true',
    help='If this flag is set, REG will run importance sampling and create new isf. ')
parser.add_argument('--nis', default = 1e6, type = int,
    help='Number of samples for importance sampling. ')

### generate inputf parser
parser.add_argument('--snp', default='SNP', type=str,
    help='Name of SNP column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--z', default='Z', type=str,
    help='Name of z-score column (if not a name that ldsc understands). NB: case insensitive.')

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

        if args.sumstats is None:
            raise ValueError('Must specify --sumstats with data folder');
 
        if args.sumstats is not None: 
            if args.test is not None:
                fl = list_files_in_directory_slience(args.sumstats)
                ft = verify_test(args.test.split(','), fl)
            if not check_xtractmat_output(args.sumstats):
                raise ValueError('Cannot find xtract outputs')
            if (True):
                inputf = args.out+'.ss.gz'; Sgf = args.out+'.sg'; Ref = args.out+'.re'
                generate_meta_analysis_input_files(inputf, Sgf, Ref, ft, fl, args, log)
            if args.isf is None:
                log.log('No inputs for --isf: Use _isf_/{}.is instead.'.format(args.out));
                ensure_dir(os.path.join(args.sumstats,'_isf_')); isf = os.path.join(args.sumstats,'_isf_',args.out+'.is');
            else:
                isf = os.path.join(args.sumstats,'_isf_',args.isf+'.is');
                if not os.path.exists(isf):
                    raise ValueError('Can\'t find {}'.format(isf));
            if not float(args.nis) > float(1e5):
                log.log('We advise you to have --nis greater than 100000')

            meta_analysis(inputf, isf, Sgf, Ref, args, log)

    except Exception:
        ex_type, ex, tb = sys.exc_info();
        log.log(traceback.format_exc(ex));
        raise
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))

