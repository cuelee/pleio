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

#from framework.meta_analysis import *
from framework.parse import *
from framework.importance_sampling import importance_sampling as imsa
from pval_estim.estim import pfun_estim, pestim
from meta_code.regeneral import REG_optim
import numpy as np
import pandas as pd
import os, sys, traceback, argparse, time




codename = 'REG'
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
            Re=pd.DataFrame(load_mat(self.re_path), columns = self.LIST, index = self.LIST); 
            return(Sg, Re)

        def __init__(self, args, outd, outname):
            try:
                ssd = args.sumstats;
                self.sumstat_dir = ssd;
                self.default_nis = args.nis;
                self.output_dir = os.path.join(ssd,outd);
                self.sg_path = os.path.join(self.output_dir,outname+'.sg');
                self.re_path = os.path.join(self.output_dir,outname+'.re');
                self.list_path = os.path.join(self.output_dir,outname+'.list');
                self.LIST = self.read_list()
                self.Sg, self.Re = self.read_matrices();
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
            df = pd.read_csv(issf,sep='\t')
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

    def generate_mod_from_cain(cain,ncain,log):
        merged_df = pd.DataFrame(columns=['SNP'])
        for i in range(ncain):
            merged_df = pd.merge(merged_df, cain[i].annot[['SNP',cain[i].ssn]], how='outer', on = ['SNP'])
        cmod_df = pd.DataFrame(merged_df[['SNP']])
        cmod_df['mod'] = 0;
        for i in range(ncain):
            idf = pd.DataFrame(cain[i].annot[['SNP']])
            idf['mod'] = int(2**(i));
            cmod_df = pd.concat([cmod_df,idf]).groupby('SNP')['mod'].sum().reset_index()
        merged_df = pd.merge(merged_df, cmod_df, how='outer', on=['SNP'])
        mod_count = cmod_df.groupby('mod').count().to_dict()['SNP']
        return(merged_df,mod_count)

    class mode_class_array_object_generator(object):
        def __init__ (self, ssin, merged_df, index, nmod):
            self.index = index;  self.nmod = nmod;
            self.metain, self.Sg, self.Re = self.generate_input(ssin, merged_df);
             
        def generate_input(self, ssin, merged_df):
            on = self.read_binstr(ssin.LIST, bin(self.index));
            row_idx = col_idx = [i for i in range(len(ssin.LIST)) if ssin.LIST[i] in on];
            mod_Sg = ssin.Sg.loc[on, on];
            mod_Re = ssin.Re.loc[on, on];
            selected_df = merged_df.loc[merged_df['mod'] == self.index][['SNP'] + on];
            return(selected_df, mod_Sg, mod_Re)

        def read_binstr(self, ssl, binstr):
            b = True; vec = [];
            for c in binstr:
                if c == 'b':
                    b = False;
                    continue
                if b:
                    continue
                vec.append(c);
            reversed_vec = [ i for i in reversed(vec) ]
            j = len(reversed_vec);
            on = [ ssl[i] for i in range(j) if reversed_vec[i] == '1' ]
            return( on )

    def generate_mode_class_array(ssin, merged_df, mod_count):
        sorted_mod = sorted(mod_count, key = mod_count.get, reverse=True)
#        mode_carray = [None]*len(sorted_mod);
#        for i in range(len(sorted_mod)):
#            index = sorted_mod[i];
        mode_carray = [None];
        for i in range(1):
            index = 2**len(ssin.LIST)-1;
            mode_carray[i] = mode_class_array_object_generator(ssin, merged_df, int(index), int(mod_count[index])) 
            merged_df = merged_df.loc[merged_df['mod'] != index]
        return(mode_carray)

    try:
        ncain = len(cain);
        merged_df, mod_count = generate_mod_from_cain(cain, ncain, log);
        mode_carray = generate_mode_class_array(ssin, merged_df, mod_count);
    except:
        log.log('FATAL ERROR: {} cannot parse the input_class_array_object'.format(codename)) 
        raise
    return(mode_carray)


def meta_analysis(df, isres, args, log):
    '''
    The function integrates multiple summary statistics and provide data frames of pvalues and pooled log OR 
    df is a python class object which contains 1. summary statistics(dataframe), 
    2. Genetic Covariance matrix(numpy matrix) 3. Random errors(numpy matrix) 
    isres is a python class object which contains tabulated importance sampling results to estimate pvalues  
    '''
    Sg=df.Sg; Re=df.Re; dat = df.metain;

    res = pd.DataFrame(dat[['SNP']])

    def rowwise_optim(x, se, Sg, Re, n):
        s = REG_optim(beta=x, stder=se, Sg=Sg, Re=Re, n=n)
        return(s)

    def rowwise_pestim(s,mtck,itck,dtab,tck_lim):
        p = pestim(cstat=s, iso = isres)
        return(p)

    n = Sg.shape[1];
    se=[1]*n
    res['stat'] = dat.apply(lambda x: rowwise_optim(x.tolist()[1:], se, Sg, Re, n), axis=1)
    res['Pvalue'] = res.apply(lambda x: rowwise_pestim(x['stat'], mtck, itck, etck, tck_lim), axis=1)

    return(res)


def mvp(args,log):

    def print_camod(camod):
        for i in range(len(camod)):
            print(camod[i].metain)
            print(camod[i].Sg)
            print(camod[i].Re)

   ### Read inputs as a class_array_object 
    ssin = setup_input_files(args ,log);
    cain = generate_input_class_array_object(ssin,log);
    log.log('Read {} Summary statistics from the summary statistics directory '.format(len(cain)));
    #log.log('Read {} Summary statistics'.format(len(cain)));
    camod = generate_mode_class_array_object(ssin, cain, log);

    outdir = args.out; ensure_dir(outdir);
    for i in range(len(camod)):
        is_dir = os.path.join(outdir,'_isfs_'); ensure_dir(is_dir)
        importance_sampling_fn = 'is_mod_{mode_number}_result'.format(mode_number=str(i));
        is_path = os.path.join(is_dir,importance_sampling_fn); 
        
        if args.create:
            imsa(N = ssin.default_nis ,Sg = camod[i].Sg, Re = camod[i].Re, outfn = is_path);
            isres = pfun_estim(isf);
        else:
            quit('The feature has not been supported yet');
 
        res = meta_analysis(camod[i], isres, args, log)
    
        result_fn = 'res_mod_{mode_number}_result.gz'.format(mode_number=str(i));
        res_path = os.path.join(outdir,result_fn);
        res.to_csv(res_path, index = False, sep='\t', compression='gzip');

    quit('work done')



parser = argparse.ArgumentParser()
parser.add_argument('--out', default='out',type=str,
    help='Output file directory name. If --out is not set, REG will use out as the '
    'default output directory.')
# Besic Random Effect Meta-analysis Flags
parser.add_argument('--sumstats', default=None, type=str,
    help='input folder path: path to summary statistics data.')
parser.add_argument('--test', default=None, type=str,
    help='comma separated sumstats filename: Use the same file name from ./{sumstats}/_out_/sumstats.list ')
parser.add_argument('--all', default=False, action='store_true',
    help='Test all possible combinations in /{sumstats}/_out_/sumstats.list ')
parser.add_argument('--isf', default=None, type=str,
    help='Filename Prefix for Estimated null-distribution cumulative densities. ')
parser.add_argument('--create', default = False, action='store_true',
    help='If this flag is set, REG will run importance sampling and create new isf. ')
parser.add_argument('--nis', default = int(100000), type = int,
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
        header += '\n'.join(options).replace('True','').replace('False','');
        header = header[0:-1]+'\n';
        log.log(header);
        log.log('Beginning analysis at {T}'.format(T=time.ctime()));
        start_time = time.time()

        if args.all is not False or args.test is not None:
            if args.sumstats is None:
                raise ValueError('--sumstats should be defined')
            if args.create is not True and args.isf is None:
                raise ValueError('--create or --isf flag should be defined') 

            mvp(args,log)

        else:
            print (header)
            print ('Error: no analysis selected.')
            print ('reg.py -h describes options.')

    except Exception:
        log.mlog(traceback.format_exc());
        raise

    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))

## Rg and Re should have same dimension
