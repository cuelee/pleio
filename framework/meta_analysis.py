from pval_estim.estim import regPestim, pfun_estim
from importance_sam.importance_sampling_REG import run_importance_sampling as runis
from meta_code.regeneral import REG_optim
from meta_code.LS import LS_chi
from framework.parse import *
from decimal import *
import os, sys
import numpy as np
import gzip

def onis(N,Sg,Re,isf):
    runis(N= int(N), GenCov=Sg, RECor=Re, outfn = isf);

def read_snp_and_z(fh,s,z):
    open_, comp_ = get_compression(fh); data = dict(); flr = True;
    with open_(fh, 'rt') as fin:
        for line in fin:
            if flr:
                col = [x.upper() for x in  line.strip().split()];
                si = col.index(s); zi = col.index(z); flr = False;
                continue
            std = line.strip().split(); data[std[si]] = std[zi];
    return(data)

def find_common_snps(data,il):
    for i in range(len(il)):
        if i == 0:
            common = set(list(data[il[i]].keys()));
        common = set(common) & set(list(data[il[i]].keys()));
    return(common)

def _generate_matrix_(matrix, fl, ft, new_file):
    m = load_mat(matrix);
    sub_ind  = [ i for i in range(len(fl)) if fl[i] in ft]
    write_mat(m[np.ix_(sub_ind,sub_ind)], new_file)

    
def generate_meta_analysis_input_files(inputf, Sgf, Ref, ft, fl,  args, log):
    log.log('Generate meta-anlaysis input file format....')
    if args.snp is not None: snp_colname = 'SNP'; 
    if args.z is not None: z_colname = 'Z'; 
    cur_path = args.sumstats; test_list = ft; input_list = [separate_extension(x)[0] for x in test_list]; data = dict();
    sum_sg = os.path.join(args.sumstats,'_out_/sumstats.sg'); sum_re = os.path.join(args.sumstats,'_out_/sumstats.re')
    Sg=_generate_matrix_(sum_sg, fl, ft, Sgf); Re=_generate_matrix_(sum_re, fl, ft, Ref);
    for i in range(len(input_list)): data[input_list[i]] = read_snp_and_z(os.path.join(cur_path,test_list[i]),snp_colname, z_colname)
    common_snps = find_common_snps(data,input_list)
    log.log('Total {} snps are in common between {} summary statistics'.format(len(common_snps),len(input_list)));
    with gzip.open(inputf, 'wt') as fout:
        for snp in common_snps:
            txt = snp;
            for l in input_list: txt += ' {} 1'.format(str(data[l][snp]))
            print(txt, file = fout); 

def meta_analysis(inputf, isf, Sgf, Ref, args, log):
    Sg=load_mat(Sgf); Re=load_mat(Ref);

    if args.create:
            log.log('Importance sampling will generate new isf with {} samples. '.format(args.nis));
            onis(N = args.nis, Sg = Sg, Re = Re, isf = isf);

    with gzip.open(inputf,"r") as fin:
        i = 0;
        for line in fin:
            i = i + 1;
        ninput = i;
        del i
        
    reg_mtck, reg_itck, reg_etck, tck_lim = pfun_estim(isf)
    
    with gzip.open(inputf,"rt") as fin, gzip.open(args.out+'.txt.gz',"wt") as fout:
        Sreg = [0.0]*ninput;
        Preg = [Decimal(0)]*ninput;
        i = 0;
        for line in fin:
            cin = line.strip().split();
            n = int(len(cin[1:])/2);
            cvar = cin[0];
            css = cin[1:];
            cbeta = np.array([float(css[i*2]) for i in range(n)]);
            cstder = np.array([float(css[i*2+1]) for i in range(n)]);
            Sreg[i] = REG_optim(beta=cbeta, stder=cstder, Sg=Sg, Re=Re, n=n)
            Preg[i] = regPestim(cstat=Sreg[i], mtck=reg_mtck, inter_tck=reg_itck, extra_tck=reg_etck, tck_lim=tck_lim)
            print(" ".join(map(str,[cvar, Sreg[i], Preg[i]])),file = fout);
            i = i + 1;

