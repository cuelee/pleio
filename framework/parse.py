''' 
(C) 2018 Cue Hyunkyu Lee

This module contains functions for parsing various reg-defined file formats.

'''

from __future__ import division, print_function
import numpy as np
#import pandas as pd
import os
import gzip
import bz2
from meta_code.LS import LS

_munge_args_list_ = ['N_con_col', 'signed_sumstats', 'daner', 'info_list', 'keep_maf', 'merge_alleles', 'a1_inc', 'N_cas_col', 'frq', 'n_min', 'N_col', 'info_min', 'chunksize', 'maf_min', 'N_cas', 'N', 'a1', 'a2', 'snp', 'N_con', 'info', 'nstudy', 'ignore', 'p', 'no_alleles', 'nstudy_min', 'daner_n']

def ensure_dir(file_path):
    try:
        os.mkdir(file_path);
    except OSError:
        if not os.path.isdir(file_path):
            raise

def separate_extension(fh):
    '''Which sort of compression should we use with read_file'''
    spd = fh.split(".");
    filename = spd[0]; ext = "."+".".join(map(str,spd[1:]));
    return(filename, ext)

def get_compression(fh):
    '''Which sort of compression should we use with read_file'''
    if fh.endswith('gz'):
        compression = 'gzip';
        openfunc = gzip.open;
    elif fh.endswith('bz2'):
        compression = 'bz2';
        openfunc = bz2.BZ2File;
    else:
        compression = None;
        openfunc = open;

    return (openfunc,compression)

def list_files_in_directory(path, log):
    '''Return all files in a directory'''
    files = [file for file in os.listdir(path)
        if os.path.isfile(os.path.join(path, file))]
    f_txt = '\nRead total {} summary statistics\n'.format(len(files))
    for cfile in files:
        f_txt += '{} '.format(cfile)
    log.log(f_txt)
    return(files)

def list_files_in_directory_slience(path):
    '''Return all files in a directory'''
    files = [file for file in os.listdir(path)
        if os.path.isfile(os.path.join(path, file))]
    return(files)

def generate_munge_arguments(opts,non_defaults):
    '''Generate arguments for mungesumstats.py'''
    munge_non_defaults = list(set(_munge_args_list_)&set(non_defaults))
    munge_args=dict((x,opts[x]) for x in munge_non_defaults)
    return(munge_args)

def checkeq_list_dic(l,d):
    for fn in l:
        if fn in d:
            continue;
        else:
            raise ValueError('Can\'t find {} in --binary-prev'.format(fn));
    return( len(l) == len(d) )
 
def read_prev_file(fh, inputs, log):
    fin = open(fh, 'r');
    pprev_dic = dict(); sprev_dic = dict();
    prev_txt = '\nSample and populaiton prevalences: \n'
    for line in fin:
        std = line.strip().split();
        prev_txt += '{} {} {}\n'. format(std[0],std[1],std[2])
        sprev_dic[std[0]] = float(std[1]); pprev_dic[std[0]] = float(std[2]);
    if (checkeq_list_dic(inputs,sprev_dic) and checkeq_list_dic(inputs,pprev_dic)):
        log.log(prev_txt); return(sprev_dic, pprev_dic);
    else:
        raise ValueError('Inequality found in --sumstats and --binary-prev ')

def make_blank_dics(inputs):
    pprev_dic = dict(); sprev_dic = dict();
    for cinput in inputs:
        sprev_dic[cinput] = ""; pprev_dic[cinput] = "";
    return(sprev_dic, pprev_dic);

def get_prev_for_sumstats(fn, sprev, pprev):
    sprev_fn = sprev[fn]; pprev_fn = pprev[fn];
    if sprev_fn != '' and pprev_fn != '':
        return(' --samp-prev {} --pop-prev {}'.format(sprev_fn, pprev_fn));
    else:
        return('');

def generate_ldsc_rginput(f1, f2, sprev, pprev, stat_path, ldsc_ld_chr):
    sp_f1 = sprev[f1]; pp_f1 = sprev[f1]; sp_f2 = sprev[f2]; pp_f2 = sprev[f2];
    sg = '.sumstats.gz'; rg_fp = os.path.join(stat_path, f1+sg)+','+os.path.join(stat_path, f2+sg);
    if sp_f1 != '' and pp_f1 != '' and sp_f2 != '' and pp_f2 != '':
        o = (' --samp-prev {},{} --pop-prev {},{}'.format(sp_f1, sp_f2, pp_f1, pp_f2))
    else:
        o = '';
    c = '/ldsc.py --w-ld-chr {} --ref-ld-chr {} --rg {} {}'.format(ldsc_ld_chr, ldsc_ld_chr, rg_fp, o)
    return(c)

def read_ldsc_intercept(fh):
    sub_str = 'Intercept:'
    with open(fh,'r') as fin:
        for line in fin:
            if sub_str in line:
                return(float(line.strip().split()[-2]))
    raise 

def read_ldsc_covariance(fh):
    sub_str = 'scale gencov:'
    with open(fh,'r') as fin:
        for line in fin:
            if sub_str in line:
                return(float(line.strip().split()[-2]))
    raise 

def write_mat(X,fname):
    np.savetxt(fname=fname,X=X,fmt='%1.5f',delimiter=' ',newline='\n',header='',footer='',comments='#')
    return();

def load_mat(fname):
    X = np.loadtxt(fname=fname);
    return(X);

def write_vec(l,fname):
    with open(fname,'w') as fout:
        for ind in l:
            print(ind, file = fout);
    return();

def parse_indiv_file(fn):
    snp = [];
    data = dict();
    necessary = ['SNP', 'A1', 'A2', 'N', 'Z'];
    index = ['']*5;
    with gzip.open(fn, 'r') as fin:
        flr = True
        for line in fin:
            if flr:
                flr = False; col = line.strip().split();
                for i in range(len(necessary)):
                    if necessary[i] in col:
                       index[i] = col.index(necessary[i]);
                    else:
                        raise
                continue
            std = line.strip().split();
            if len(std) > 1:
                data[std[index[0]]] = [std[x] for x in index];
                snp.append(std[index[0]]);
    return(data,snp)

def correct_inflation_in_sumstats(intercept, ofh, nfh):
    '''Will perform Z/alpha of all rows in ofh and save them in nfh'''
    flr = True;
    with gzip.open(ofh) as fin, gzip.open(nfh,'w') as fout:
        for line in fin:
            if flr:
                flr = False; col = line.strip().split(); z = col.index('Z'); lcol = len(col);
                print(line,file=fout);
                continue
            std = line.strip().split();
            if len(std) == lcol:
                zs = float(std[z]) / np.sqrt(intercept); std[z] = str(zs);
                print(' '.join(std), file=fout);
            else:
                continue

def fixed_effect_meta_analysis(d1, d2, common, fp12):
    col = 'SNP A1 A2 N Z';
    cor = np.diag([1.0]*2); stders = [1.0]*2;
    with gzip.open(fp12, 'w') as fout:
        print(col,file = fout);
        for s in common:
            snp1, snp2 = d1[s], d2[s]; betas = [snp1[4], snp2[4]]; nz = LS(betas, stders, cor);
            if snp1[:3] == snp2[:3]:
                nn = str(float(snp1[3])+float(snp2[3]));
                line = '{} {} {} {} {}'.format(snp1[0],snp1[1],snp1[2],nn,nz);
                print(line,file = fout)
            else:
                raise ValueError('filed snp parsing {} - {}'.format(snp1,snp2));

def write_output(sumfiles, gen_cov, env_cor, prefix):
    write_vec(sumfiles, prefix + '.list');
    write_mat(gen_cov, prefix + '.sg')
    write_mat(env_cor, prefix + '.re')

def write_mat(X,fname):
    np.savetxt(fname=fname,X=X,fmt='%1.5f',delimiter=' ',newline='\n',header='',footer='',comments='#')
    return();

def load_mat(fname):
    X = np.loadtxt(fname=fname);
    return(X);
