'''
(C) 2018 Cue Hyunkyu Lee

This module deals with creating genetic covariance matrix and environmental correlation matrix
by performing pairwise LDSC implemented in LDSC software.
'''

from __future__ import print_function
from framework.parse import *
import math
import os,sys
import subprocess 

def run_munge_sumstats(sumstats_path, munge_path, files, munge_command_pre, log):
    '''we run LDSC munge sumstats to ensure the quality of summary statistics data'''
    txt = 'Generate munge.sumstats.gz....'
    ensure_dir(munge_path);
    for ff in files:
        fn, ext = separate_extension(ff);
        comm = '{} --out {} --sumstats {}'.format(munge_command_pre, os.path.join(munge_path, fn),os.path.join(sumstats_path, ff));
        subprocess.call(comm, shell = True, stdout=subprocess.PIPE);
    txt += '\tDone.'
    log.log(txt);

def intercept_corrected_sumstats(munge_path, indiv_ldsc_path, inputs, sprev, pprev, ldsc_code_path, ldsc_ld_chr, log):
    ''' We run LDSC on each summary statistics and correct the inflation of chi^2 statistics using inflaction factor alpha ''' 
    ensure_dir(indiv_ldsc_path);
    log.log('\nEstimating LDSC intercept...');
    txt = '';
    for fn in inputs:
        fp = fn + '.sumstats.gz'; munge_fp = os.path.join(munge_path, fp); ldsc_fp = os.path.join(indiv_ldsc_path, fp);
        prevs = get_prev_for_sumstats(fn, sprev, pprev);
        comm = '{}/ldsc.py --w-ld-chr {} --ref-ld-chr {} --h2 {} --out {} {}'.format(ldsc_code_path, ldsc_ld_chr, ldsc_ld_chr, munge_fp, ldsc_fp, prevs);
        subprocess.call(comm, shell = True, stdout=subprocess.PIPE);
        fp_intercept = read_ldsc_intercept(ldsc_fp+'.log');
        txt += '{}: {}, '.format(fn, fp_intercept);
        correct_inflation_in_sumstats(max(1.0,fp_intercept), munge_fp, ldsc_fp);
    log.log(txt[:-2]);

def generate_genetic_covariance_matrix(indiv_ldsc_path, gen_pair_path, inputs, sprev, pprev, ldsc_code_path, ldsc_ld_chr, log):
    ''' We run LDSC on each pair of summary statistics and get genetic covariance matrix '''
    ensure_dir(gen_pair_path);
    log.log('\nEstimating genetic covariance using Pair Wise LDSC...')
    n = len(inputs);
    m = np.diag([0.0]*n);
    for i in range(n):
        f1 = inputs[i];
        txt = '{} - '.format(f1);
        for j in range(n-i):
            f2 = inputs[j+i];
            comm_pre = generate_ldsc_rginput(f1,f2,sprev,pprev,indiv_ldsc_path,ldsc_ld_chr); outf = os.path.join(gen_pair_path,'{},{}'.format(f1,f2))
            comm = '{}{} --out {}'.format(ldsc_code_path, comm_pre, os.path.join(gen_pair_path,'{},{}'.format(f1,f2)),);
            subprocess.call(comm, shell = True, stdout=subprocess.PIPE);
            c = read_ldsc_covariance(outf + '.log');
            txt += '{}: {}, '.format(f2, c);
            m[i,j+i] = c;
            m[j+i,i] = c;
        log.log(txt[:-2]);
    return(m)

def fixed_effect_framework(input_path, fn1, fn2, new_path, log):
    sg = '.sumstats.gz';
    fp1 = os.path.join(input_path, fn1+sg); fp2 = os.path.join(input_path, fn2+sg);
    d1, s1 = parse_indiv_file(fp1); d2, s2 = parse_indiv_file(fp2);
    common = set(s1) & set(s2);
    fp12 = os.path.join(new_path, fn1+'_'+fn2+'.sumstats.gz')
    fixed_effect_meta_analysis(d1, d2, common, fp12);
    return(fp12)

def generate_environment_correlation_matrix(indiv_ldsc_path, env_pair_path, inputs, ldsc_code_path, ldsc_ld_chr, log):
    ''' We apply fixed effect meta-analysis on each pair of summary statistics and run LDSC to get environmental correlation matrix from intercept of LDSC '''
    ensure_dir(env_pair_path);
    log.log('\nEstimating environment correlation matrix...');
    n = len(inputs);
    m = np.diag([1.0]*n);
    for i in range(n-1):
        f1 = inputs[i];
        txt = '{} - '.format(f1);
        for j in range(n-i-1):
            f2 = inputs[j+i+1];
            fp12 = fixed_effect_framework(indiv_ldsc_path, f1, f2, env_pair_path, log); 
            comm = '{}/ldsc.py --w-ld-chr {} --ref-ld-chr {} --h2 {} --out {}'.format(ldsc_code_path, ldsc_ld_chr, ldsc_ld_chr, fp12, fp12);
            subprocess.call(comm, shell = True, stdout=subprocess.PIPE);
            intercept = read_ldsc_intercept(fp12+'.log');
            c = intercept - 1
            txt += '{}: {}, '.format(f2, c);
            m[i,j+i+1] = c;
            m[j+i+1,i] = c;
        log.log(txt[:-2]);
    return(m)


def pairwise_ldsc(args, munge_args, log):
    ''' Will perform pairwise LDSC '''
    _cur_path_ = args.sumstats; _file_list_ = list_files_in_directory(_cur_path_, log); _input_list_ = [separate_extension(x)[0] for x in _file_list_]; 
    _ldsc_code_path_ = args.ldsc; _ldsc_ld_chr_ = args.ldsc_ld_chr;
    _munge_dir_='_munged_'; _alpha_dir_ = '_ind_ldsc_'; _gencov_dir_ = '_gen_cov_'; _envcor_dir_ = '_env_cor_'; _meta_dir_='_meta_';
    _out_dir_ = os.path.join(_cur_path_,'_out_'); ensure_dir(_out_dir_);

    _munge_command_='{}/munge_sumstats.py {}'.format(_ldsc_code_path_,' '.join(map(str,['--'+x.replace('_','-')+' '+str(munge_args[x]) for x in munge_args])));
    if args.binary:
        _sprev_dic_, _pprev_dic_ = read_prev_file(args.binary_prev, _input_list_ ,log);
    else:
        _sprev_dic_, _pprev_dic_ = make_blank_dics(_input_list_);
    run_munge_sumstats(_cur_path_,os.path.join(_cur_path_,_munge_dir_), _file_list_ , _munge_command_, log);
    intercept_corrected_sumstats(os.path.join(_cur_path_,_munge_dir_), os.path.join(_cur_path_,_alpha_dir_), _input_list_, _sprev_dic_, _pprev_dic_, _ldsc_code_path_, _ldsc_ld_chr_, log);
    genetic_covariance_matrix = generate_genetic_covariance_matrix(os.path.join(_cur_path_,_alpha_dir_), os.path.join(_cur_path_,_gencov_dir_), _input_list_, _sprev_dic_, _pprev_dic_, _ldsc_code_path_, _ldsc_ld_chr_, log);
    environment_correlation_matrix = generate_environment_correlation_matrix(os.path.join(_cur_path_,_alpha_dir_), os.path.join(_cur_path_, _envcor_dir_), _input_list_, _ldsc_code_path_, _ldsc_ld_chr_, log);

    write_output(_input_list_, genetic_covariance_matrix, environment_correlation_matrix, os.path.join(_out_dir_,'sumstats'));
