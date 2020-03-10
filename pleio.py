#!/usr/bin/env python3
'''
PLEIO: PLEIOtropy
Copyright(C) 2018 Cue Hyunkyu Lee 

PLEIO is command line framework to perform cross-disease meta-analyses
'''

#from framework.meta_analysis import *
from framework.parse import *
from framework.importance_sampling import importance_sampling
from framework.significance_estimation import pfun_estim, pvalue_estimation
from meta_code.variance_component import vcm_optimization
from meta_code.LS import LS, LS_p 
import numpy as np
import pandas as pd
import os, sys, traceback, argparse, time
import multiprocessing as mp
from itertools import product

codename = 'PLEIO'
__version__ = 'beta 1.0.0'
MASTHEAD = "************************************************************\n"
MASTHEAD += "* PLEIOtropy({c})\n".format(c=codename)
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

def setup_input_files(args ,log, outd = '_out_', outname = 'sumstats' ):

    class _out_(object):
        def read_list(self, fl = []):
            with open(self.list_path,'r') as fin:
                for line in fin:
                    fl.append(line.strip());
            return(fl)
                    
        def read_matrices(self):
            Sg=pd.DataFrame(load_mat(self.sg_path), columns = self.LIST, index = self.LIST); 
            Rn=pd.DataFrame(load_mat(self.rn_path), columns = self.LIST, index = self.LIST); 
            return(Sg, Rn)

        def __init__(self, args, outd, outname):
            try:
                ssd = args.sumstats;
                self.sumstat_dir = ssd;
                self.default_nis = args.nis;
                self.output_dir = os.path.join(ssd,outd);
                self.sg_path = os.path.join(self.output_dir,outname+'.sg');
                self.rn_path = os.path.join(self.output_dir,outname+'.rn');
                self.list_path = os.path.join(self.output_dir,outname+'.list');
                self.LIST = self.read_list()
                self.Sg, self.Rn = self.read_matrices();
            except:
                log.log('FATAL ERROR: Exception occurred while parsing sumstats folder [_out_]');
                raise

    return(_out_(args, outd, outname))

def generate_input_class_array_object(ssin, log, issd = '_ind_ldsc_'):

    class input_class_array_object_generator(object):

        def __init__ (self, ssn, issf):
            self.ssn = ssn;
            self.issf = issf;
            self.annot = self.read_df(issf);

        def read_df(self,issf):
            df = pd.read_csv(issf,sep='\t',nrows =100000)
            df.columns = ['SNP', 'A1', 'A2', 'N', self.ssn]
#            print('{} has {} lines'.format(self.ssn,df.shape[0]))
#            with pd.option_context('display.max_rows', 1, 'display.max_columns', None):
#                print(df)
            return (df)

    try:
        issl = [ os.path.join(ssin.sumstat_dir, issd, ssin.LIST[i]+'.sumstats.gz') for i in range(len(ssin.LIST))]; nss = len(issl);
        input_carray = [input_class_array_object_generator(ssn = ssin.LIST[i], issf=issl[i]) for i in range(nss)]
#        input_class_array_objects = [input_class_array_object_generator(ssf=ssl[i],issf=issl[i],index=i) for i in range(2)]

    except:
        log.log('FATAL ERROR: Exception occurred while parsing sumstats folder [_ind_ldsc_]')
        raise 

    return(input_carray)

def generate_mode_class_array_object(ssin, cain, log):

    def generate_df_from_cain(cain,ncain,log):
        merged_df = pd.DataFrame(columns=['SNP'])
        for i in range(ncain):
            merged_df = pd.merge(merged_df, cain[i].annot[['SNP',cain[i].ssn]], how='outer', on = ['SNP'])
        merged_df = merged_df.dropna(axis = 0, how = 'any', inplace = False)
        merged_df = merged_df.set_index('SNP')
        return(merged_df)

    class meta_class_array_object_generator(object):
        def __init__ (self, ssin, merged_df, ncain):
            self.metain = merged_df;
            self.Sg = ssin.Sg;
            self.Rn = ssin.Rn;
            self.n = ncain;
             
    try:
        ncain = len(cain);
        merged_df = generate_df_from_cain(cain, ncain, log);
        meta_cain = meta_class_array_object_generator(ssin, merged_df, ncain);
    except:
        log.log('FATAL ERROR: {} cannot parse the input_class_array_object'.format(codename)) 
        raise
    return(meta_cain)

def sqrt_ginv (X, tol = 2.22044604925e-16**0.5):
    u,s,vh = np.linalg.svd(X)
    Positive = s > max(tol * s[0], 0)  
    if (all(Positive)):
        res=np.transpose(vh).dot(np.diag(1/s)**0.5).dot(np.transpose(u))
    elif (not any(Positive)): 
        res=np.array([0]*np.prod(X.shape)).reshape(X.shape)
    else:
        res=vh[:, Positive].dot(np.diag(1/s[Positive])**0.5).dot(np.transpose(u[:, Positive]))
    return(res)

def run_vc_optimizer(x, Ut, Ct):
    b = np.array([x[i*2] for i in range(n)]);
    se = np.array([x[i*2+1] for i in range(n)]);
    h = Ut.dot(b)
    K = np.diag(se).dot(Ct).dot(np.diag(se))
    return(vcm_optimization(h, K))

def LS_input_parser(x, Ce, n):
    b= [x[i*2] for i in range(n)];
    se = [x[i*2+1] for i in range(n)];
    return(LS(b, se, Ce))

def _estimate_statistics(df_data, Sg, Ce):
    n = np.size(Sg,1)
    sqrt_Sg_ginv = sqrt_ginv(Sg) 
    trans_ce = sqrt_Sg_ginv.dot(Ce).dot(sqrt_Sg_ginv)
    df_out = pd.DataFrame(index = df_data.index)
    df_out['DELPY_stat'] = df_data.apply(lambda x: run_vc_optimizer(x.tolist(), sqrt_Sg_ginv, trans_ce), axis=1)
    df_out['LS_stat'] = df_data.apply(lambda x: LS_input_parser(x.tolist(), Ce, n), axis=1)
    return(df_out)

def _parallelize(meta_cain, func, args): 
    data_split = np.array_split(meta_cain.metain, args.ncores)
    iterable = product(data_split, [meta_cain.Sg.values], [meta_cain.Ce.values])
    pool = mp.Pool(int(args.ncores))
    df_output = pd.concat(pool.starmap(func, iterable))
    pool.close()
    pool.join()
    return(df_output)

def find_closed_PSD(x):
    w, v = np.linalg.eig(x.values)
    if (w>0).all():
        return(x)
    else: 
        w_correct = [max(i,1e-6) for i in w]
        rho = sum(w)/sum(w_correct)
        s=v.dot(np.diag(w_correct)).dot(np.transpose(v))
        x = pd.DataFrame((s+np.transpose(s))/2,index = x.index, columns = x.columns, copy = True  )
        return(x); 

def is_pos_def(x):
    eigs = np.linalg.eigvals(x) 
    return np.all(eigs > 0)


def read_filetype(filename, magic_dict = { "\x1f\x8b\x08": "gzip", "\x42\x5a\x68": "bz2", "\x50\x4b\x03\x04": "zip" }):

    max_len = max(len(x) for x in magic_dict)

    with open(filename, encoding="latin-1") as f:
        file_start = f.read(max_len)
        print(file_start)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return None


def read_metain(metainf, sgf, cef):
    class meta_class_array_object_generator(object):
        def __init__(self, metain_df, sg_df, ce_df, N_gwas, tlist):
            self.metain = metain_df;
            self.Sg = sg_df;
            self.Ce = ce_df;
            self.N_gwas = N_gwas;
            self.n = len(tlist);

    def read_mat(fn, LIST, name):
        df=pd.read_csv(fn,sep='\t', compression = 'infer'); 
        df.index = df.columns 
        df = df.loc[LIST,LIST]
        if not is_pos_def(df):
            print(name,'is not a positive definite matrix')
        return(df)

    metain_df = pd.read_csv(metainf, sep='\t', compression = 'infer', index_col = ['SNP','A1','A2'])
    col = metain_df.columns.tolist()
    include = [col[i].split('_beta')[0] for i in range(len(col)) if i % 2 == 0]
    N_gwas = (metain_df.loc[:,[col[i] for i in range(len(col)) if i % 2 != 0]]).apply(lambda x: np.mean(1/x**2))
    error_correlation_df = read_mat(cef,include,'error correlation matrix');
    genetic_covariance_df = read_mat(sgf,include,'genetic covariance matrix');
    meta_cain = meta_class_array_object_generator(metain_df, genetic_covariance_df, error_correlation_df, N_gwas, include)
    return(meta_cain)

def delpy(args,log):

   ### Read inputs as a class_array_object 
    if args.sumstats is not None:
        ssin = setup_input_files(args ,log);
        cain = generate_input_class_array_object(ssin,log);
        log.log('Read {} Summary statistics from the summary statistics directory '.format(len(cain)));
        #log.log('Read {} Summary statistics'.format(len(cain)));
        meta_cain = generate_mode_class_array_object(ssin, cain, log);
        
        del ssin, cain
    
    if args.metain is not None:
        meta_cain = read_metain(args.metain, args.sg, args.ce)

    if args.create:
        importance_sampling_fn = 'isf.txt';
        outdir = args.out; ensure_dir(outdir);
        importance_sampling_f = os.path.join(outdir,importance_sampling_fn); 
        importance_sampling(args.nis, meta_cain.N_gwas, meta_cain.Sg, meta_cain.Ce, importance_sampling_f, args.ncores);
        iso = pfun_estim(importance_sampling_f);
    else:
        quit('The feature has not been supported yet');
    
    summary = _parallelize(meta_cain, _estimate_statistics, args);
    summary['DELPY_p'] = summary.apply(lambda x: pvalue_estimation(x['DELPY_stat'], iso), axis=1);
    summary['LS_p'] = summary.apply(lambda x: (x['LS_stat'], ), axis=1);
    
    out_path = os.path.join(outdir+'.delpy.sum.gz');
    summary.to_csv(out_path, index = True, sep='\t', compression='gzip');

    quit('work done')


parser = argparse.ArgumentParser()
parser.add_argument('--out', default='out',type=str,
    help='Output file directory name. If --out is not set, PLEIO will use out as the '
    'default output directory.')
# Besic Random Effect Meta-analysis Flags
parser.add_argument('--sumstats', default=None, type=str,
    help='input folder path: path to summary statistics data.')
parser.add_argument('--metain', default=None, type=str,
    help='input file: file prefix of the meta input data.')
parser.add_argument('--sg', default=None, type=str,
    help='input file: file prefix of the genetic covariance matrix.')
parser.add_argument('--ce', default=None, type=str,
    help='input file: file prefix of the non-genetic correlation matrix.')
parser.add_argument('--test', default=None, type=str,
    help='comma separated sumstats filename: Use the same file name from ./{sumstats}/_out_/sumstats.list ')
parser.add_argument('--all', default=False, action='store_true',
    help='Test all possible combinations in /{sumstats}/_out_/sumstats.list ')
parser.add_argument('--isf', default=None, type=str,
    help='Filename Prefix for Estimated null-distribution cumulative densities. ')
parser.add_argument('--create', default = False, action='store_true',
    help='If this flag is set, DELPY will run importance sampling and create new isf. ')
parser.add_argument('--nis', default = int(100000), type = int,
    help='Number of samples for importance sampling. ')
parser.add_argument('--parallel', default = False, action='store_true',
    help='If this flag is set, DELPY will run parallel computing ')
parser.add_argument('--ncores', default = mp.cpu_count(), type = int, 
    help='Number of cpu cores for parallel computing. DELPY uses max(CPUs) as default.')
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
        header += './delpy.py \\\n';
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults];
        header += '\n'.join(options).replace('True','').replace('False','');
        header = header[0:-1]+'\n';
        log.log(header);
        log.log('Beginning analysis at {T}'.format(T=time.ctime()));
        start_time = time.time()

        if args.all is not False or args.test is not None:
            if args.sumstats is None and args.metain is None:
                raise ValueError('--sumstats should be defined')
            if args.metain is not None and args.sg is None or args.ce is None:
                raise ValueError('--metain need --sg and --ce flags to be defined')
            if args.create is not True and args.isf is None:
                raise ValueError('--create or --isf flag should be defined') 
            if (args.parallel == False):
                args.ncores = 1;
            
            delpy(args,log)

        else:
            print (header)
            print ('Error: no analysis selected.')
            print ('Use delpy.py -h for details.')

    except Exception:
        log.mlog(traceback.format_exc());
        raise

    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))

## Sg and Rn should have same dimension
