#!/usr/bin/env python2

'''
PLEIO: Pleiotropic Locus Exploration and Interpretation using Optimal test
Copyright(C) 2018 Cue Hyunkyu Lee 

PLEIO is a summary-statistic-based framework to map and interpret pleiotropic loci in a joint analysis of multiple diseases and complex traits. Our method maximizes power by systematically accounting for genetic correlations and heritabilities of the traits in the association test.
'''

from __future__ import division
import os, sys, traceback, argparse, time
import pandas as pd
import numpy as np

from itertools import product
from scipy import stats

codename = 'ldsc_preprocess'
__version__ = 1.0
MASTHEAD = "************************************************************\n"
MASTHEAD += "* ldsc_preprocess({c})\n".format(c=codename)
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C)2018 Cue H. Lee\n"
MASTHEAD += "* Seoul National University\n"
MASTHEAD += "* Unlicensed software\n"
MASTHEAD += "************************************************************\n"

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

def _remove_dtype(x):
    '''Removes dtype: float64 and dtype: int64 from pandas printouts'''
    x = str(x)
    x = x.replace('\ndtype: int64', '')
    x = x.replace('\ndtype: float64', '')
    return x

class Logger(object):
    '''
    Lightweight logging.

    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'wb')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.

        '''
        print >>self.log_fh, msg
        print msg

class read_input_file(object):
    def __init__(self,args):
        self.args = args
        self.out = args.out
        self.input = args.input
        self.temp_dir = self.set_out('temp')
        self.snplist_out = self.set_out('snps.txt.gz')
        self.sg_out = self.set_out('sg.txt.gz')
        self.ce_out = self.set_out('ce.txt.gz')
        self.res_out = self.set_out('metain.txt.gz')
        self.read_input() 
        self.sumstat = dict.fromkeys(self.trait_id)
        self.create_folder(self.out)
        self.create_folder(self.temp_dir)
        
    def set_out(self, p):
        return(os.path.join(self.out,p))
    
    def read_input(self):
        try:
            self.info = pd.read_csv(self.input, header = 'infer', compression = 'infer', sep = '\t', na_filter = True, index_col = None, dtype = {'FILE': str,'TYPE': str,'SPREV': float,'PPREV': float, 'NAME': str}, usecols=['FILE', 'TYPE', 'SPREV', 'PPREV', 'NAME']) #info = ldsc_preprocess_input
            if 'NAME' in self.info.columns:
                self.trait_id = self.info.NAME.to_numpy()
            else:
                self.trait_id = np.array(['trait_' + str(i) for i in self.info.shape[0]])
                self.info.insert(self.info.shape[1], 'NAME', self.trait_id, True)
            if not np.all( (self.info.shape[1] == 5) and np.array_equal(set(self.info.columns.to_numpy()), set(np.array(['FILE', 'TYPE', 'SPREV', 'PPREV','NAME',])) )): 
                raise ValueError('Header should contains following columns: FILE, TYPE, SPREV, PPREV')
            self.info.rename(columns={"FILE": "f", "TYPE": "phe", "SPREV": "sprev", "PPREV": "pprev", "NAME": "trait_id"}, inplace = True) 
            self.info.set_index('trait_id',inplace=True)
            n_col = self.info.shape[1]
            self.info.insert(n_col, 'fn',  np.vectorize(lambda x: x.split('/')[-1])(self.info.f), True)
            sumstat = self.vectorize_join_fn(self.info.fn.to_numpy(),'.sumstats.gz', fn = lambda x, y : "".join([x,y]))
            self.info.insert(n_col, 'sumstat', sumstat, True)
            temp = self.vectorize_join_fn(self.temp_dir, sumstat, fn = lambda x, y : os.path.join(x,y))
            self.info.insert(n_col, 'temp', temp, True)
            out = self.vectorize_join_fn(self.out, sumstat, fn = lambda x, y : os.path.join(x,y))
            self.info.insert(n_col, 'out', out, True)
            for f in self.info.f:
                if not os.path.isfile(f):
                    raise OSError('The file does not exist in the following path. Please consider providing the full path to the file in the input.')

        except ValueError:
            log.log ('Failed reading the input file. The input should be a T x 4(5) matrix with a header where T is the number of traits.')
            raise
        else:
            log.log ('Read %d traits from input' % self.info.shape[0])
    
    def vectorize_join_fn(self,a,b,fn):
        return(np.vectorize(fn)(a,b))
    
    def create_folder(self, p):
        try:
            os.mkdir(p)
        except OSError:
            log.log ("Failed to create a directory at : %s" % p)
        else:
            log.log ("Successfully created a directory at : %s " % p)

    def remove_folder(self):
        p = self.temp_dir
        try:
            for root, dirs, files in os.walk(p, topdown=False):
                for name in files:
                    os.remove(os.path.join(root, name))
                for name in dirs:
                    os.rmdir(os.path.join(root, name))
            os.rmdir(p)
        except OSError:
            log.log ("Failed to remove the directory at : %s" % p)
        else:
            log.log ("Successfully removed the directory at : %s " % p)
            
            
class generate_sumstat_data(object):
    def __init__(self, data):
        def sanity_check(fs):
            for f in fs:
                column_names = pd.read_csv(f, sep ='\t', index_col='SNP', header = 'infer', usecols = ['SNP','A1','A2','Z','N'], dtype = {'SNP': str,'A1': str,'A2': str,'Z': float,'N': float}, nrows = 0).columns
                if (len(column_names) != 4) or (not np.array_equal(set(column_names), set(np.array(['A1','A2','Z','N'])))):
                    raise ValueError('The summary statistics data file must contain the following six columns: SNP,A1,A2,Z,N')
            
        self.data = data
        self.sumstats = dict()
        self.snps = []        
        sanity_check(data.info.f.to_numpy())
        for trait in data.info.index:
            self.sumstats[trait] = self.h2_analysis(trait)
        for trait in data.info.index:
            self.snps += [self.sumstats[trait].index.to_list()]
        self.common_snps = np.array(list(set.intersection(*map(set,self.snps))))
        log.log('Generate input files (sumstats.txt.gz) for LDSC --rg analysis')
        for trait in data.info.index:
            self.sumstats[trait] = self.sumstats[trait].loc[self.common_snps,:]
            self.sumstats[trait].to_csv(self.data.info.temp[trait], compression = 'gzip', index = True, header = True, sep = '\t')
            log.log('Generated %s' %self.data.info.temp[trait])
        self.rg_analysis()
        self.check_sumstats()
        log.log('Final step: combine several summary statistical data')
        self.res = pd.DataFrame(index = self.common_snps)
	self.res.index.rename(name='SNP', inplace=True)
        for trait in data.info.index:
            sprev, pprev, ptype = self.data.info.loc[trait,['sprev','pprev','phe']]
            eta, se = generate_pleio_input(self.sumstats[trait], trait, sprev, pprev, ptype)
            eta.name = str(trait) + '_beta'
            se.name = str(trait) + '_se'
            self.res = self.res.merge(eta, left_index = True, right_index = True)
            self.res = self.res.merge(se, left_index = True, right_index = True)
            self.sumstats[trait] = None
        self.res.to_csv(self.data.res_out, index = True, header = True, sep = '\t', compression = 'gzip')

    def h2_analysis(self, t):
        def read_intercept(f):
            with open(f, 'rb') as fin:
                for line in fin:
                    if 'Intercept: ' in line:
                        intercept = float(line.split()[1])
            return(intercept)            

        def correction_of_sumstat_data_using_icoef (i, c): 
            d = pd.read_csv(i, sep ='\t', index_col='SNP', usecols = ['SNP','A1','A2','Z','N'], dtype = {'SNP': str,'A1': str,'A2': str,'Z': float,'N': float})
            d.Z = d.Z / np.sqrt(c)
            return(d)

        try:
            h2_args = self.data.args
            h2_args.h2_log = h2_args.h2 = h2_args.samp_prev = h2_args.pop_prev = None;
            h2_args.h2 = self.data.info.f[t] 
            h2_args.phe = self.data.info.phe[t]
            h2_args.h2_log = os.path.join(self.data.temp_dir, t + '_h2.log') 
            
            if h2_args.phe == 'binary' or h2_args.phe == 'bin' or h2_args.phe == 'b':
                h2_args.samp_prev, h2_args.pop_prev = np.vectorize(lambda x: float(x) if not np.isnan(x) else None)(self.data.info.loc[t,['sprev','pprev']].to_numpy())
                if((h2_args.samp_prev is None) != (h2_args.pop_prev is None)):
                    raise ValueError('sprev and pprev must both be floating-point values, or NA values.')
            else: 
                h2_args.samp_prev = h2_args.pop_prev = None

            ldsc_h2_log = Logger(h2_args.h2_log)
            sumstats.estimate_h2(h2_args, ldsc_h2_log)
            ldsc_h2_log.log_fh.close()
            
            if os.path.isfile(h2_args.h2_log):
                intercept = read_intercept(args.h2_log)
            else:
                raise RuntimeError('{} not exist'.format(h2_args.h2_log))
            
            self.data.info.loc[t,'temp'] 
            log.log('Dividing z-scores with the correction factor (the squared root of the LDSC h2 analysis intercept value).')
            sumstat = correction_of_sumstat_data_using_icoef(args.h2, intercept)
            return(sumstat)
        except RuntimeError:
            raise('Interrupted during LDSC h2 analysis. %s' % t)    

    def rg_analysis(self):
        def read_gencov(f):
            if not os.path.isfile(f):
                raise OSError('log file is missing: %s' % f)
            with open(f, 'rb') as fin:
                while 1:
                    line = fin.readline()
                    if not line:
                        break
                    elif 'Total Liability scale gencov:' in line:
                        genetic_covariance = float(line.split()[4])
                    elif 'Summary of Genetic Correlation Results' in line:
                        fin.readline()
                        genetic_covariance_intercept = fin.readline().split()[10]
            return(genetic_covariance_intercept, genetic_covariance)

        try:
            traits = self.data.info.index.to_list()  
            traits_to_be_tested = list(traits) 
            T = len(traits)
            sg = pd.DataFrame([[0.0]*T]*T, index = traits, columns = traits)
            ce = pd.DataFrame([[0.0]*T]*T, index = traits, columns = traits)
            for A in traits:
                for B in traits_to_be_tested:
                    rg_args = self.data.args
                    rg_args.no_check_alleles = True
                    rg_args.rg = rg_args.r2_log = rg_args.samp_prev = rg_args.pop_prev = None
                    rg_args.rg = ','.join([self.data.info.temp[A],self.data.info.temp[B]])
                    rg_args.r2_log = os.path.join(self.data.temp_dir, ','.join([A,B]) + '.log')
                    rg_args.samp_prev = ','.join(map(str,self.data.info.loc[[A,B],'sprev'].to_list()))
                    rg_args.pop_prev = ','.join(map(str,self.data.info.loc[[A,B],'pprev'].to_list()))
                    ldsc_rg_log = Logger(args.r2_log)
                    sumstats.estimate_rg(rg_args, ldsc_rg_log)
                    ldsc_rg_log.log_fh.close()
                    gencov_int, gencov = read_gencov(args.r2_log)
                    sg.loc[A,B]=sg.loc[B,A] = gencov
                    ce.loc[A,B]=ce.loc[B,A] = gencov_int
                traits_to_be_tested.remove(A)
            sg.to_csv(self.data.sg_out, index = False, header = True, sep ='\t', compression = 'gzip')        
            ce.to_csv(self.data.ce_out, index = False, header = True, sep ='\t', compression = 'gzip')
            return(None)
        except RuntimeError:
            raise('Interrupted during LDSC rg analysis.')         

    def check_sumstats(self):
        def find_duplicated_snps(i):
            return(i[i.duplicated()].to_list())

        def check_allele_mismatch(df):
            res = ~df.eq(df.iloc[:, 0],axis = 0).all(axis = 1)
            if np.sum(res) > 0:
                s = df.index[res].to_numpy(str)
                raise ValueError('Found Allele mismatch: {}'.format(s))
        
        log.log('Number of variants in common: {}'.format(len(self.common_snps)))

        duplicated_snps = []
        for trait in self.data.info.index:
            duplicated_snps += (find_duplicated_snps(self.sumstats[trait].index))
        log.log('Found %d duplicated variants' %len(duplicated_snps))

        for trait in self.data.info.index:
            self.sumstats[trait].drop(duplicated_snps,inplace=True)
            
        for trait in self.data.info.index:
            Nvariants = self.sumstats[trait].shape[0]
            self.sumstats[trait].dropna(axis = 0, how = 'any', inplace = True)    
            diff = Nvariants - self.sumstats[trait].shape[0]
            if diff > 0:
                raise ValueError('Found {} variants with NA values: {}'.format(diff, trait))

        self.common_snps = self.sumstats[trait].index.to_numpy()
        
        A1 = pd.DataFrame(index = self.common_snps)
        A2 = pd.DataFrame(index = self.common_snps)
        for trait in self.data.info.index:
            A1[trait] = self.sumstats[trait].A1
            A2[trait] = self.sumstats[trait].A2
        check_allele_mismatch(A1)
        check_allele_mismatch(A2)
        log.log('Found no Allele mismatch')

def generate_pleio_input(sumstat, trait, sprev, pprev, ptype, eta = 'eta', se = 'se', Z = 'Z', N = 'N'):
    def estimate_coef(sprev, pprev):
        '''
        This function estimate the alpha and theta coefficients.
        These coefficients will be used in the standardization of summary statistics.
        The procedure converts observed beta and SE(beta) estimates of binary traits to standardized scales; eta and SE(eta) 
        '''
        K = float(pprev)
        P = float(sprev)

        t = stats.norm.isf(K)
        z = stats.norm.pdf(t,0,1)
        i = z/K # i2 = -i*K / (1-K)
        lamb = (P-K)/(1-K)
        theta = (i*lamb)*(t-i*lamb)
        alpha = (K*(1-K))**2/(P*(1-P)*z**2)
        return(alpha, theta)    
    
    def gen_beta(Z, N, alpha, theta):
        if (Z == 0):
            return(float(0.0))
        se = np.sqrt(1/N)
        beta = Z * se
        eta = np.sqrt(alpha) * beta / np.sqrt(1-alpha*theta*beta**2)
        return(eta)
    
    def gen_se(eta, Z, N):
        if (Z == 0):
            return(np.sqrt(1/N))
        se = eta / Z
        return(se)
    
    log.log('Summary of {}: '.format(trait))
    if ptype == 'binary' or ptype == 'bin' or ptype == 'b': 
        alpha, theta = estimate_coef(sprev, pprev);
        log.log('Estimated alpha and theta for {}: {} and {}'.format(trait,alpha,theta))
        sumstat[eta] = sumstat.loc[:,[Z,N]].apply(lambda x: gen_beta(x.Z,x.N,alpha,theta),axis=1) 
        sumstat[se] = sumstat.loc[:,[eta,Z,N]].apply(lambda x: gen_se(x.eta,x.Z,x.N),axis=1)
    else: 
        sumstat[se] = np.sqrt(1/sumstat.N)
        sumstat[eta] = sumstat.Z*sumstat.se
    log.log('z-score (mean) {} \nsample number (mean) {}'.format(np.mean(sumstat.Z), np.mean(sumstat.N)))
    return(sumstat[eta], sumstat[se])
    
    
def preprocess(args, log):

    data = read_input_file(args)
    sumstat_data = generate_sumstat_data(data)
    data.remove_folder()
    
    return(None)


parser = argparse.ArgumentParser()
parser.add_argument('--out', default='input_folder',type=str,
    help='Output file directory name. If --out is not set, PLEIO will use input_folder as the '
    'default output directory.')
# Basic Flags for PLEIO-preprocess
parser.add_argument('--input', default=None, type=str,
    help='The input file for the preprocessing. '
    'This file is a T x 4 matrix with header row (tab delimited) where T is number of traits. Header should contains following: FILE: filename, TYPE: either binary or quantitative, SPREV: sample prevalence = cases_samples/ control_samples, PPREV: disease prevalence')
# Basic Flags for Working with Variance Components
parser.add_argument('--ref-ld', default=None, type=str,
    help='Use --ref-ld to tell LDSC which LD Scores to use as the predictors in the LD '
    'Score regression. '
    'LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix.')
parser.add_argument('--ref-ld-chr', default=None, type=str,
    help='Same as --ref-ld, but will automatically concatenate .l2.ldscore files split '
    'across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz '
    'to the filename prefix. If the filename prefix contains the symbol @, LDSC will '
    'replace the @ symbol with chromosome numbers. Otherwise, LDSC will append chromosome '
    'numbers to the end of the filename prefix.'
    'Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gz'
    'Example 2: --ref-ld-chr ld/@_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz')
parser.add_argument('--w-ld', default=None, type=str,
    help='Filename prefix for file with LD Scores with sum r^2 taken over SNPs included '
    'in the regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz.')
parser.add_argument('--w-ld-chr', default=None, type=str,
    help='Same as --w-ld, but will read files split into 22 chromosomes in the same '
    'manner as --ref-ld-chr.')
parser.add_argument('--overlap-annot', default=False, action='store_true',
    help='This flag informs LDSC that the partitioned LD Scores were generates using an '
    'annot matrix with overlapping categories (i.e., not all row sums equal 1), '
    'and prevents LDSC from displaying output that is meaningless with overlapping categories.')
parser.add_argument('--print-coefficients',default=False,action='store_true',
    help='when categories are overlapping, print coefficients as well as heritabilities.')
parser.add_argument('--frqfile', type=str,
    help='For use with --overlap-annot. Provides allele frequencies to prune to common '
    'snps if --not-M-5-50 is not set.')
parser.add_argument('--frqfile-chr', type=str,
    help='Prefix for --frqfile files split over chromosome.')
parser.add_argument('--no-intercept', action='store_true',
    help = 'If used with --h2, this constrains the LD Score regression intercept to equal '
    '1. If used with --rg, this constrains the LD Score regression intercepts for the h2 '
    'estimates to be one and the intercept for the genetic covariance estimate to be zero.')
parser.add_argument('--intercept-h2', action='store', default=None,
    help = 'Intercepts for constrained-intercept single-trait LD Score regression.')
parser.add_argument('--intercept-gencov', action='store', default=None,
    help = 'Intercepts for constrained-intercept cross-trait LD Score regression.'
    ' Must have same length as --rg. The first entry is ignored.')
parser.add_argument('--M', default=None, type=str,
    help='# of SNPs (if you don\'t want to use the .l2.M files that came with your .l2.ldscore.gz files)')
parser.add_argument('--two-step', default=None, type=float,
    help='Test statistic bound for use with the two-step estimator. Not compatible with --no-intercept and --constrain-intercept.')
parser.add_argument('--chisq-max', default=None, type=float,
    help='Max chi^2.')

# Flags for both LD Score estimation and h2/gencor estimation
parser.add_argument('--print-cov', default=False, action='store_true',
    help='For use with --h2/--rg. This flag tells LDSC to print the '
    'covaraince matrix of the estimates.')
parser.add_argument('--print-delete-vals', default=False, action='store_true',
    help='If this flag is set, LDSC will print the block jackknife delete-values ('
    'i.e., the regression coefficeints estimated from the data with a block removed). '
    'The delete-values are formatted as a matrix with (# of jackknife blocks) rows and '
    '(# of LD Scores) columns.')

# Flags you should almost never use
parser.add_argument('--invert-anyway', default=False, action='store_true',
    help="Force LDSC to attempt to invert ill-conditioned matrices.")
parser.add_argument('--n-blocks', default=200, type=int,
    help='Number of block jackknife blocks.')
parser.add_argument('--not-M-5-50', default=False, action='store_true',
    help='This flag tells LDSC to use the .l2.M file instead of the .l2.M_5_50 file.')
parser.add_argument('--return-silly-things', default=False, action='store_true',
    help='Force ldsc to return silly genetic correlation estimates.')
parser.add_argument('--no-check-alleles', default=False, action='store_true',
    help='For rg estimation, skip checking whether the alleles match. This check is '
    'redundant for pairs of chisq files generated using munge_sumstats.py and the '
    'same argument to the --merge-alleles flag.')

if __name__ == '__main__':
    
    
    ldsc_dir = os.path.dirname(os.path.abspath(__file__))+'/ldsc'    
    if os.path.isdir(ldsc_dir):
        open(ldsc_dir+'/__init__.py', 'a').close()
    else: raise(ValueError('directory name {} cannot be found in directory {}. Install ldsc in the same directory as the \'preprocess.py\''.format('ldsc', os.path.dirname(os.path.abspath(__file__)))))

    import ldsc.ldscore.ldscore as ld
    import ldsc.ldscore.parse as ps
    import ldsc.ldscore.sumstats as sumstats
    import ldsc.ldscore.regressions as reg

    args = parser.parse_args()
    log = Logger(args.out+'.log');

    if args.out is None:
        raise ValueError('--out is required');

    try:
        defaults = vars(parser.parse_args(''));
        opts = vars(args);
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]];
        header = MASTHEAD;
        header += 'Call: \n';
        header += './ldsc_preprocess.py \\\n';
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults];
        header += '\n'.join(options).replace('True','').replace('False','');
        header = header[0:-1]+'\n';
        log.log(header);
        log.log('Beginning analysis at {T}'.format(T=time.ctime()));
        start_time = time.time()

        if args.n_blocks <= 1:
            raise ValueError('--n-blocks must be an integer > 1.')
        elif args.input and (args.ref_ld or args.ref_ld_chr) and (args.w_ld or args.w_ld_chr):
            if args.ref_ld and args.ref_ld_chr:
                raise ValueError('Cannot set both --ref-ld and --ref-ld-chr.')
            if args.w_ld and args.w_ld_chr:
                raise ValueError('Cannot set both --w-ld and --w-ld-chr.')
        
            if not args.overlap_annot or args.not_M_5_50:
                if args.frqfile is not None or args.frqfile_chr is not None:
                    log.log('The frequency file is unnecessary and is being ignored.')
                    args.frqfile = None
                    args.frqfile_chr = None
            if args.overlap_annot and not args.not_M_5_50:
                if not ((args.frqfile and args.ref_ld) or (args.frqfile_chr and args.ref_ld_chr)):
                    raise ValueError('Must set either --frqfile and --ref-ld or --frqfile-chr and --ref-ld-chr')
            preprocess(args,log)
        else:
            print header
            print 'Error: no analysis selected.'
            print 'preprocess.py -h describes options.'
    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.log( traceback.format_exc(ex) )
        raise
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))

## Sg and Rn should have same dimension
