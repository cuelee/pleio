#!/usr/bin/env python
'''
Copyright(C) 2018 Cue Hyunkyu Lee 
'''
from __future__ import print_function
from framework.parse import *
from framework.bineff2liabeff import bineff2liabeff
import pandas as pd
import os, sys, traceback, argparse, time
import multiprocessing as mp 
from itertools import product

MASTHEAD = "************************************************************\n"
MASTHEAD += "* Summary statistics to Liability (sum2lib)\n"
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
        Print to log file and stdout with a single comman.
        '''
        print (msg, file = self.log_fh);
        print (msg);


parser = argparse.ArgumentParser()
parser.add_argument('--out', default='sum2lib', type=str, 
        help='out filename prefix.')
parser.add_argument('--gzip', default=False, action='store_true',
        help='gzip outfile.')
parser.add_argument('--log', default='sum2lib', type=str, 
        help='log filename prefix.')
parser.add_argument('--sumstats', default=None, type=str,
        help='input file paths, list all summary statistics with commas')
parser.add_argument('--prevs', default=None, type=str,
        help='input binary trait prevalences, list all prevalences of multiple traits with commas')
parser.add_argument('--ncase', default='NCASE', type=str,
        help='Name of Ncase column (if not a name that delpy understands). NB: case insensitive.')
parser.add_argument('--ncont', default='NCONT', type=str,
        help='Name of Ncont column (if not a name that delpy understands). NB: case insensitive.')
parser.add_argument('--signed-sumstats', default=None, type=str,
        help='Name of signed sumstat column, comma null value (e.g., Z,0 or OR,1)')
parser.add_argument('--Z', default='Z', type=str,
        help='Name of Z column (if not a name that delpy understands). NB: case insensitive.')
parser.add_argument('--OR', default='OR', type=str,
        help='Name of OR column (if not a name that delpy understands). NB: case insensitive.')
parser.add_argument('--se', default='SE', type=str, 
        help='Name of SE column (if not a name that delpy understands). NB: case insensitive.')
parser.add_argument('--beta', default='BETA', type=str,
        help='Name of log(OR) column (if not a name that delpy understands). NB: case insensitive.')
parser.add_argument('--snp', default='SNP', type=str,
        help='Name of SNP column (if not a name that delpy understands). NB: case insensitive.')
parser.add_argument('--a1', default='A1', type=str,
        help='Name of risk allele column (if not a name that delpy understands). NB: case insensitive.')
parser.add_argument('--a2', default='A2', type=str,
        help='Name of reference allele column (if not a name that delpy understands). NB: case insensitive.')
parser.add_argument('--delim', default=' ', type=str,
        help='The delimiter of the summary statistics text file. The default value is \' \'')


def generate_input(args, log, input_dict = dict()):
    fn_sumstats = args.sumstats.split(',');
    val_prevs = args.prevs.split(',');
    n = len(fn_sumstats);
    if( n == len(val_prevs) ):
        for i in range(n):
            input_dict[fn_sumstats[i]] = val_prevs[i];  
    if(len(input_dict) is not n ):
            raise ValueError('Detected an input with the same value ')
    return(input_dict)

def parse_delim(di):
    if di == 'w':
        d = ' '
    elif di == 'c':
        d = ','
    elif di == 't':
        d = '\t'
    else:
        raise ValueError('The value of the flag \'--delim\' cannot be recognized. Check sum2lib.py -h for details.');
    return(d)

def read_header(fn, args):
    data = pd.read_csv(fn, sep=args.delim, nrows = 1);
    in_str = [args.snp, args.a1, args.a2, args.signed_sumstats, args.se, args.ncase, args.ncont]
    col = list(data.keys())
    COL = [a.upper() for a in col];
    if all([a.upper() in COL for a in in_str]):
       return([COL.index(a.upper()) for a in in_str],col);
    else:
        print([a.upper() for a in in_str])
        print(COL)
        raise ValueError('Cannot parse the GWAS header. See sum2lib.py -h for details.')

def parse_signed_sumstats(s):
    a, b = s.split(',');
    if (a.upper() == 'Z' and int(b) == 0):
        r = args.Z;
        return(r,'Z')
    elif (a.upper() == 'OR' and int(b) == 1):
        r = args.OR
        return(r,'OR')
    elif (a.upper() == 'BETA' and int(b) == 0):
        r = args.beta;
        return(r,'BETA')
    else:
        raise ValueError('The value of --signed-sumstats cannot be recognized. Check sum2lib.py -h for details.');

def sum2lib_sub(df, prev, signed):
    applied_df = df.apply(lambda x: bineff2liabeff(float(x.iloc[2]),float(x.iloc[3]),float(x.iloc[4]),float(x.iloc[5]), float(prev), signed, x.name), axis=1, result_type = 'expand') 
    applied_df.columns=['LIB_BETA','LIB_SE']
    df = pd.concat([df,applied_df], axis = 'columns')
    return(df)

def parallelize (df_input, func, prev, signed):
    cores = mp.cpu_count()
    data_split = np.array_split(df_input, cores)
    iterable = product(data_split, [prev], [signed]);
    pool = mp.Pool(int(cores))
    df_output = pd.concat(pool.starmap(func, iterable))
    pool.close()
    pool.join()
    return(df_output)

def sum2lib(input_dic, args, log):
    args.delim = parse_delim(args.delim);
    if len(args.signed_sumstats.split(',')) > 2:
        raise ValueError('The value of --signed-sumstats cannot be recognized. Check sum2lib.py -h for details.');
    args.signed_sumstats, signed = parse_signed_sumstats(args.signed_sumstats);
    
    fs = list(input_dic.keys());
    for i in range(len(fs)):
        if os.path.exists(fs[i]):
            prev = input_dic[fs[i]];
            fn = fs[i];
        log.log('Read summary statistics from: {GWAS}'.format(GWAS=fn));
        fn_col, col= read_header(fn, args)
        df = pd.read_csv(fn, sep=args.delim, usecols = fn_col)
        df = df[[col[a] for a in fn_col]]
        df = df.set_index(df.columns[0])
        df = parallelize(df, sum2lib_sub, prev, signed)
        df.columns = ['A1','A2','OLD_'+signed,'OLD_SE','Ncase','Ncont',signed,'SE']
        col=list(df.columns)
        df = df[ [col[a] for a in [0,1,6,7,4,5] ]]
    return(df)

if __name__ == '__main__':
    args = parser.parse_args();
    if args.log is None:
        raise ValueError('--log is required');
    log = Logger(args.log+'.log');
    start_time = time.time();
    try:
        default = vars(parser.parse_args(''));
        opts = vars(args);
        non_defaults = [x for x in opts.keys() if opts[x] != default[x]];
        munge_args = generate_munge_arguments(opts, non_defaults);
        header = MASTHEAD;
        header += 'Call: \n';
        header += './sum2lib.pu \\\n';
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults];
        header += '\n'.join(options).replace('True','').replace('False','');
        header = header[0:-1]+'\n';

        if args.sumstats is None or args.prevs is None or args.signed_sumstats is None:
            raise ValueError('The --sumstats, --prevs, --signed-sumstats flags are required');

        log.log(header);
        log.log('Beginning analysis at {T}'.format(T=time.ctime()));

        input_dic = generate_input(args, log);
        if (len(input_dic) > 0):
            res_df = sum2lib(input_dic, args, log);
            if args.gzip:
                res_df.to_csv(args.out+'.txt.gz', index = True, sep = '\t', compression = 'gzip')
            else:
                res_df.to_csv(args.out+'.txt', index = True, sep = '\t')
        else:
            print(header);
            print('ERROR: no analysis selected.');
            print('sum2lib.py -h for more details.');
    except Exception:
        ex_type, ex, tb = sys.exc_info();
        formatted_lines = traceback.format_exc().splitlines()
        for i in range(len(formatted_lines)): log.log(formatted_lines[i]); 
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()));
        time_elapsed = round(time.time()-start_time,2);
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)));

