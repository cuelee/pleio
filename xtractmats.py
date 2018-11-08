#!/usr/bin/env python
'''
REG : Pair Wise LDSC(getMartices)
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
from __future__ import print_function 
import framework.pairwise_ldsc as pairwise_ldsc
from framework.parse import *
import os, sys, traceback, argparse, time

MASTHEAD = "************************************************************\n"
MASTHEAD += "* Pair Wise LDSC (GetMatrices)\n"
MASTHEAD += "* (C)2018 Cue H. Lee\n"
MASTHEAD += "* Seoul National University\n"
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


parser = argparse.ArgumentParser()
parser.add_argument('--log', default='regmat',type=str,
    help='log filename prefix.')
parser.add_argument('--sumstats', default=None, type=str,
    help='input folder path: path to summary statistics data.')
parser.add_argument('--binary', default=False, action='store_true',
    help='If this flag is set, xtractmat will run binary traits pair wise LDSC. ')
parser.add_argument('--binary-prev', default=None, type=str,
    help='prevalence filename for binary traits pair wise LDSC. ')
parser.add_argument('--ldsc', default=os.path.join(os.path.dirname(os.path.realpath(__file__)),'ldsc'), type=str,
    help='LDSC path. If --ldsc is not set, xtractmat will use ./reg/ldsc.'
    'as the default ldsc package path.')
parser.add_argument('--ldsc-ld-chr', default=None, type=str,
    help='Filename prefix for LDSC --ref-ld-chr and --w-ld-chr flag.')

### munge_sumstat_options 
parser.add_argument('--N', default=None, type=float,
                    help="Sample size If this option is not set, will try to infer the sample "
                    "size from the input file. If the input file contains a sample size "
                    "column, and this flag is set, the argument to this flag has priority.")
parser.add_argument('--N-cas', default=None, type=float,
                    help="Number of cases. If this option is not set, will try to infer the number "
                    "of cases from the input file. If the input file contains a number of cases "
                    "column, and this flag is set, the argument to this flag has priority.")
parser.add_argument('--N-con', default=None, type=float,
                    help="Number of controls. If this option is not set, will try to infer the number "
                    "of controls from the input file. If the input file contains a number of controls "
                    "column, and this flag is set, the argument to this flag has priority.")
parser.add_argument('--info-min', default=0.9, type=float,
                    help="Minimum INFO score.")
parser.add_argument('--maf-min', default=0.01, type=float,
                    help="Minimum MAF.")
parser.add_argument('--daner', default=False, action='store_true',
                    help="Use this flag to parse Stephan Ripke's daner* file format.")
parser.add_argument('--daner-n', default=False, action='store_true',
                    help="Use this flag to parse more recent daner* formatted files, which "
                    "include sample size column 'Nca' and 'Nco'.")
parser.add_argument('--no-alleles', default=False, action="store_true",
                    help="Don't require alleles. Useful if only unsigned summary statistics are available "
                    "and the goal is h2 / partitioned h2 estimation rather than rg estimation.")
parser.add_argument('--merge-alleles', default=None, type=str,
                    help="Same as --merge, except the file should have three columns: SNP, A1, A2, "
                    "and all alleles will be matched to the --merge-alleles file alleles.")
parser.add_argument('--n-min', default=None, type=float,
                    help='Minimum N (sample size). Default is (90th percentile N) / 2.')
parser.add_argument('--chunksize', default=5e6, type=int,
                    help='Chunksize.')

# optional args to specify column names
parser.add_argument('--snp', default=None, type=str,
                    help='Name of SNP column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--N-col', default=None, type=str,
                    help='Name of N column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--N-cas-col', default=None, type=str,
                    help='Name of N column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--N-con-col', default=None, type=str,
                    help='Name of N column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--a1', default=None, type=str,
                    help='Name of A1 column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--a2', default=None, type=str,
                    help='Name of A2 column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--p', default=None, type=str,
                    help='Name of p-value column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--frq', default=None, type=str,
                    help='Name of FRQ or MAF column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--signed-sumstats', default=None, type=str,
                    help='Name of signed sumstat column, comma null value (e.g., Z,0 or OR,1). NB: case insensitive.')
parser.add_argument('--info', default=None, type=str,
                    help='Name of INFO column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--info-list', default=None, type=str,
                    help='Comma-separated list of INFO columns. Will filter on the mean. NB: case insensitive.')
parser.add_argument('--nstudy', default=None, type=str,
                    help='Name of NSTUDY column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--nstudy-min', default=None, type=float,
                    help='Minimum # of studies. Default is to remove everything below the max, unless there is an N column,'
                    ' in which case do nothing.')
parser.add_argument('--ignore', default=None, type=str,
                    help='Comma-separated list of column names to ignore.')
parser.add_argument('--a1-inc', default=False, action='store_true',
                    help='A1 is the increasing allele.')
parser.add_argument('--keep-maf', default=False, action='store_true',
                    help='Keep the MAF column (if one exists).')



if __name__ == '__main__':
    args = parser.parse_args()
    if args.log is None:
        raise ValueError('--log is required');
    
    log = Logger(args.log+'.log');
    start_time = time.time()
    try:
        if args.sumstats is None:
            raise ValueError('The --sumstats flag is required.')
        defaults = vars(parser.parse_args(''));
        opts = vars(args);
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]];
        munge_args = generate_munge_arguments(opts,non_defaults)
        header = MASTHEAD;
        header += 'Call: \n';
        header += './xtractmats.py \\\n';
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults];
        header += '\n'.join(options).replace('True','').replace('Flase','');
        header = header[0:-1]+'\n';

        if args.sumstats is not None:
            if args.binary is True:
                if args.binary_prev is None:
                    raise ValueError('Must specify --binary-prev with prevalence information.');
         
            log.log(header);
            log.log('Beginning analysis at {T}'.format(T=time.ctime()));

            pairwise_ldsc.pairwise_ldsc(args, munge_args, log);
        # bad flags
        else:
            print (header);
            print ('ERROR: no analysis selected.');
            print ('xtractmats.py -h describes options.');
    except Exception:
        ex_type, ex, tb = sys.exc_info();
        log.log(traceback.format_exc(ex));
        raise
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) );
        time_elapsed = round(time.time()-start_time,2);
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)));

