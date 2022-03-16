#!/usr/bin/env python3

'''
PLEIO: Pleiotropic Locus Exploration and Interpretation using Optimal test
Copyright(C) 2018 Cue Hyunkyu Lee 

PLEIO is a summary-statistic-based framework to map and interpret pleiotropic loci in a joint analysis of multiple diseases and complex traits. Our method maximizes power by systematically accounting for genetic correlations and heritabilities of the traits in the association test.
'''

from framework.importance_sampling import importance_sampling
from framework.assoc_test import parallel_computing, stat_estimation, blup_estimation
from framework.utilities import is_pos_def
import numpy as np
import pandas as pd
import os, sys, traceback, argparse, time
import multiprocessing as mp
import logging
from itertools import product

#logger = mp.log_to_stderr()
#logger.setLevel(mp.SUBDEBUG)

#print(mp.cpu_count())

codename = 'PLEIO'
__version__ = '2.0'
MASTHEAD = "************************************************************\n"
MASTHEAD += "* PLEIO({c})\n".format(c=codename)
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C)2018 Cue H. Lee\n"
MASTHEAD += "* Seoul National University\n"
MASTHEAD += "* Unlicensed software\n"
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
        
    def mlog(self, msg):
        '''
        [ mutelog ] Print to log file and without stdout with a single command.

        '''
        print (msg, file = self.log_fh);

class generate_data(object):
    def __init__(self, args):
        self.output = args.out;
        self.ncpu = int(args.ncores)
        
        self.pleio_output = os.path.join(self.output + '.txt.gz');
        self.sumstat_fn = args.metain;
        self.sumstat = self.read_metain(args.snp);
        self.trait_name, self.column_beta, self.column_se, self.N_trait = self.read_column()
        self.N_sample = (self.sumstat.loc[:,self.column_se]).apply(lambda x: np.mean(1/x**2))

        self.sg_fn = args.sg;
        self.ce_fn = args.ce;
        self.gencov, self.envcor = self.read_mat();
        
        self.isf_output = os.path.join(self.output + '.isf'); 
        self.run_isf = args.create
        if self.run_isf is False:
            self.isf_output = args.isf
        self.isf_Nsim = args.nis
        self.isf_zero_prob = float(0.195);
        
    def read_column(self):
        c = self.sumstat.columns.to_numpy()
        ind = np.array([i for i in range(len(c)) if i % 2 == 0])
        b = c[ind] 
        s = c[ind+1]
        t = np.vectorize(lambda x: x.split('_beta')[0])(b)
        n = int(len(t))
        return(t,b,s,n)
        
    def read_metain(self,s):
        try:
            col = pd.read_csv(self.sumstat_fn, header = 'infer',compression = 'infer', index_col = s, sep='\t', nrows = 0)
            dtype_dict = dict(); dtype_dict[s] = 'str'
            for c in col: dtype_dict[c] = 'float'
            d = pd.read_csv(self.sumstat_fn, header = 'infer',compression = 'infer', index_col = s, sep='\t', dtype = dtype_dict, na_filter = False)
        except Exception:
            raise('All data in the summary statistics file should be numeric. NA values are not allowed')
        return(d)
    
    def read_mat(self):
        df = pd.read_csv(self.sg_fn, sep='\t', compression = 'infer', header = 'infer'); 
        df.index = df.columns
        sg = df.loc[self.trait_name,self.trait_name]
        
        df = pd.read_csv(self.ce_fn, sep='\t', compression = 'infer', header = 'infer'); 
        df.index = df.columns 
        ce = df.loc[self.trait_name,self.trait_name]
        return(sg, ce)
    
def pleio(args,log):
    ### Read inputs as a class_array_object 
    data = generate_data(args)
    
    ### PLEIO itself modifies the PSD matrix to PD. The code below will tell you whether this fix has been made or not.
    if not is_pos_def(data.gencov):
        log.log('Genetic covariance matrix is not a positive definite matrix')
    if not is_pos_def(data.envcor):
        log.log('Environmental correlation matrix is not a positive definite matrix')

    ### importance sampling procedure
    if data.run_isf:
        data.isf_zero_prob = importance_sampling(data.isf_Nsim, data.N_sample, data.gencov, data.envcor, data.isf_output, data.ncpu);
    
    ### PLEIO analysis
    log.log('Running PLEIO...');

    parallel_input = product(np.array_split(data.sumstat, data.ncpu), [data.gencov.values], [data.envcor.values], [data.isf_output])
    res = parallel_computing(parallel_input, stat_estimation, data.ncpu);

    if args.flattening_p_values:
        res = flattening_p_value(res, data.N_sample, data.gencov, data.envcor, data.ncpu , data.isf_output )
    
    log.log('Writing PLEIO output file');
    res.to_csv(data.pleio_output, index = True, sep='\t', compression='gzip');

    if args.blup:
        log.log('Estimating BLUP...');
        dat = data.sumstat.merge(res.loc[:,'tausq'], how = 'inner', left_index = True, right_index = True)
        iterable = product(np.array_split(dat, data.ncpu), [data.gencov.values], [data.envcor.values],[data.trait_name])
    
        blup_estimates = parallel_computing(iterable, blup_estimation, data.ncpu);
        blup_estimates.to_csv(data.output + '.blup.gz', index = True, sep ='\t', compression = 'gzip')

    log.log('Work done');
    

parser = argparse.ArgumentParser()
parser.add_argument('--out', default='pleio',type=str,
    help='Path to the output. If --out is not set, PLEIO will use pleio as the default output directory.')
parser.add_argument('--metain', default=None, type=str,
    help='input file: file prefix of the meta input data.')
parser.add_argument('--sg', default=None, type=str,
    help='input file: file prefix of the genetic covariance matrix.')
parser.add_argument('--ce', default=None, type=str,
    help='input file: file prefix of the non-genetic correlation matrix.')
parser.add_argument('--isf', default=None, type=str,
    help='Filename Prefix for Estimated null-distribution cumulative densities. ')
parser.add_argument('--create', default = False, action='store_true',
    help='If this flag is set, PLEIO will run importance sampling and create new isf. ')
parser.add_argument('--nis', default = int(100000), type = int,
    help='Number of samples for importance sampling. ')
parser.add_argument('--blup', default = False, action='store_true',
    help='If this flag is set, PLEIO will estimate Best Linear Unbiased Prediction (BLUP)'
    'and write [output].blup.gz')
parser.add_argument('--flattening_p_values', default = False, action= 'store_true',
    help='Flattening p-value distribution. This is necessary if you want to check lambda GC')
parser.add_argument('--parallel', default = False, action='store_true',
    help='If this flag is set, PLEIO will run parallel computing ')
parser.add_argument('--ncores', default = mp.cpu_count() - 2, type = int, 
    help='Number of cpu cores for parallel computing. If --ncores is not set, PLEIO will use n_cores - 1 as the default.')
parser.add_argument('--snp', default='SNP', type=str,
    help='Name of SNP column (if not a name that ldsc understands). NB: case insensitive.')

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
        header += './pleio.py \\\n';
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults];
        header += '\n'.join(options).replace('True','').replace('False','');
        header = header[0:-1]+'\n';
        log.log(header);
        log.log('Beginning analysis at {T}'.format(T=time.ctime()));
        start_time = time.time()
        
        if args.metain is not None:
            if args.metain is None:
                raise ValueError('--metain was not specified')
            if args.metain is not None and args.sg is None or args.ce is None:
                raise ValueError('--sg and --ce flags are required.')
            if args.create is not True and args.isf is None:
                raise ValueError('--create or --isf is required') 
            if (args.parallel == False):
                args.ncores = 1;
            
            pleio(args,log)

        else:
            print ('Error: no analysis selected.')
            print ('Use pleio -h for details.')

    except Exception:
        log.mlog(traceback.format_exc());
        raise

    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
