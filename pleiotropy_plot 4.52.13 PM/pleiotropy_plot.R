`%notin%` <- Negate(`%in%`)
load('/Users/cuelee/Dropbox/github/pleio/pleiotropy_plot/R/circoplot.Rdata')

## Dependencies
require(circlize)
require(ComplexHeatmap)

#### read inputs
if (!exists('datf')){
  datf = gzfile("/Users/cuelee/Dropbox/PLEIO/Meeting/2019_10_30_CIRCLIZE/metain.txt.gz")
  dat = read.table(datf,header = T, check.names = F)
  ddf = dat[,-c(1:3)]
  rownames(ddf) = dat[,1]
  
  gwas_datf = gzfile("/Users/cuelee/Dropbox/PLEIO/Meeting/2019_10_30_CIRCLIZE/gwasin.txt.gz")
  gwas_dat = read.table(gwas_datf,header = T, check.names = F)
  gwasin = gwas_dat
  
  gwasResults <- read.table('/Users/cuelee/Dropbox/PLEIO/Meeting/2019_10_15_RealAR/results.txt',header=T, sep='\t')
  
  novel = read.table(gzfile('/Users/cuelee/Dropbox/PLEIO/Meeting/2019_10_15_RealAR/results_NOVEL.txt.gz'),header = T, sep='\t')
  thres = 30
  maxp = function(x){
    return(max(10^(-thres),x))
  }
  gwasResults$DELPY_P <- sapply(gwasResults$DELPY_P,maxp)
  
  h2 = (as.vector(read.table('/Users/cuelee/Dropbox/PLEIO/Meeting/2019_10_30_CIRCLIZE/heritability.txt',header= T,check.names=F)))
  
  in_mat=as.matrix(read.table('/Users/cuelee/Dropbox/PLEIO/Meeting/2019_10_30_CIRCLIZE/rg.txt.gz',header= T,check.names=F))
  
  traits = as.matrix(read.table('/Users/cuelee/Dropbox/PLEIO/Meeting/2019_10_30_CIRCLIZE/factors.txt', sep='\t',header= F,check.names=F))
  
  rownames(in_mat) = colnames(in_mat) = sapply(colnames(in_mat), function(x) traits[which(traits[,1]==x),2])
  colnames(h2) = sapply(colnames(h2), function(x) traits[which(traits[,1]==x),2])
  
  #ind=corrMatOrder(in_mat)
  #in_mat = in_mat[ind,ind]
  manual_order = c('High density lipoprotein','Essential hypertension','Hypertension','Obesity','Fasting glucose','Type 2 diabetes','Ischemic heart disease','Coronary atherosclerosis','Coronary artery disease','Major coronary heart disease','Heart attack','Myocardial infarction','A. myocardial infarction','Heart failure','Type 1 diabetes','Triglycerides','Low density lipoprotein','Total cholesterol' )
  manual_ind = rep(0,nrow(in_mat))
  for ( i in 1:nrow(in_mat)){manual_ind[i] = which(colnames(in_mat) == manual_order[i])} 
  in_mat = in_mat[manual_ind,manual_ind]
  factors = colnames(in_mat)
}


#### Analysis 
k=0
for (i in 1:(length(novel$SNP))){
  snp = as.character(novel$SNP[i])
  snp_delpy = gwasResults[gwasResults$SNP == snp,]
  bp = snp_delpy$BP
  bp_s = bp-1000000
  bp_e = bp+1000000
  chr = snp_delpy$CHR
  p = snp_delpy$DELPY_P
  
  subset = gwasResults[gwasResults$CHR == chr &gwasResults$BP > bp_s &  gwasResults$BP < bp_e, c('SNP','BP')]
  subset_ordered = subset[order(subset$BP),]
  
  points = snp_info = NULL;
  gwasin_subset =  gwasin[rownames(gwasin)%in%subset_ordered$SNP,]
  
  for (fac in factors){
    gwasin_tn = traits[which(traits[,2]==fac),1]
    x = (subset_ordered$BP - bp)/(1000000)*0.8
    
    ps = NULL;
    p_fac = paste(gwasin_tn,'_P',sep='')
    for (s in subset_ordered$SNP){
      ps = c(ps, gwasin_subset[s,p_fac])          
    }
    y = -log10(ps)
    
    snp_info = c(snp_info, list(list(snp = snp,x = 0, y = -log10(gwasin_subset[snp,p_fac]) ) ) ) 
    points = c(points,list(cbind(x,y)) )
  }
  names(points) = factors
  names(snp_info) = factors
  
  gwas_beta = gwasin_subset[snp,paste(sapply(factors, function(x) traits[which(traits[,2]==x),1]),'_BETA',sep='')] 
  names(gwas_beta)=sapply(colnames(gwas_beta),function(x) traits[which(traits[,1]==strsplit(x,'_BETA')[[1]][1]),2])
  eta = ddf[snp,paste(sapply(factors, function(x) traits[which(traits[,2]==x),1]),'_beta',sep='')]
  names(eta)=sapply(colnames(eta),function(x) traits[which(traits[,1]==strsplit(x,'_beta')[[1]][1]),2])
  eta_se = ddf[snp,paste(sapply(factors, function(x) traits[which(traits[,2]==x),1]),'_se',sep='')]
  names(eta_se)=sapply(colnames(eta_se),function(x) traits[which(traits[,1]==substr(x,1,nchar(x)-3)),2])
  
  lci = eta - eta_se * 1.96 
  uci = eta + eta_se * 1.96 
  eta_ci = list(lci, uci)
  
  snp_h2 = eta^2
  prop_h2 = matrix(rep(0,length(factors)), nrow = 1)
  colnames(prop_h2) = factors
  for (fac in factors){
    prop_h2[1,fac] = as.numeric(snp_h2[fac])
  }
  
  ldf = gen_link(in_mat, prop_h2[1,])
  k=k+1
  jpeg(paste('/Users/cuelee/Dropbox/PLEIO/Meeting/2019_10_30_CIRCLIZE/out/',k,'_',snp,'.tiff',sep=''), width = 7, height = 7, units = 'in', res = 300)
  par(mar =c(0,0,0,0))
  ## lmp: local Manhattan plot input, ci: confidence intervals of joint_analysis input, beta: effect sizes of joint_analysis input
  
  circlize_plot(factors = factors , ldf = ldf, gwas_beta = gwas_beta, eta = eta, eta_ci = eta_ci, lmp_input = list(snp_info = snp_info, points = points))
  dev.off()
}



