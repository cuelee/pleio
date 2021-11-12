ARGS = commandArgs(trailingOnly = T)
if (length(ARGS) == 0 ){
  example_folder = '/Users/cuelee/Dropbox/github/example'
} else {example_folder = ARGS[1]}

library(pleiotropyPlot)

pleioin = gzfile(file.path(example_folder, 'input.txt.gz'))
pleiores = gzfile(file.path(example_folder, 'output.txt.gz'))
sumstats = gzfile(file.path(example_folder, 'sumstats.txt.gz'))
snps = gzfile(file.path(example_folder, 'snp.txt.gz'))
h2 = gzfile(file.path(example_folder, 'heritability.txt.gz'))
rg = gzfile(file.path(example_folder, 'rg.txt.gz'))
traits = file.path(example_folder, 'traits.txt.gz')
ref = gzfile(file.path(example_folder, '1000G_SNP_CHR_BP_A1_A2_CHR1.txt.gz'))
outf = example_folder

input_parser = function(pleioin, pleiores, sumstats, snps, h2, rg, traits, ref){
  dat = list()

  dat$pleioin = read.table(pleioin, header = T, check.names = F, row.names = 'SNP', sep='\t')
  dat$pleiores = read.table(pleiores,header = T, check.names = F, row.names = 'SNP', sep='\t')
  dat$sumstats = read.table(sumstats,header = T, check.names = F, row.names = 'SNP', sep='\t')
  dat$h2 = as.vector(read.table(h2, header= T, check.names=F))
  dat$rg = as.matrix(read.table(rg, header= T, check.names=F))
  dat$traits = as.matrix(read.table(traits, sep='\t',header= T,check.names=F))
  dat$snps = as.vector(as.matrix(read.table(snps, sep='\t',header= F,check.names=F)))
  dat$ref = read.table(ref,header = T, check.names = F, row.names = 'RSID', sep='\t')
  dat$outf = outf

  change_colname= function(col_name){
    new_colname = NA
    for (i in 1:nrow(dat$traits)){
      asubstr = dat$traits[i,'filename']
      if (grepl(asubstr, col_name)){
        new_colname = gsub(asubstr,dat$traits[i,'traitname'],col_name)
        break
      }
    }
    return(new_colname)
  }

  rownames(dat$rg) = colnames(dat$rg) = sapply(colnames(dat$rg), function(x) dat$traits[which(dat$traits[,'filename']==x),2])
  colnames(dat$h2) = sapply(colnames(dat$h2), function(x) dat$traits[which(dat$traits[,'filename']==x),2])
  colnames(dat$sumstats) = sapply(colnames(dat$sumstats), function(x) change_colname(x) )
  colnames(dat$pleioin) = sapply(colnames(dat$pleioin), function(x) change_colname(x) )

  return(dat)
}

if (!exists('dat')){ dat = input_parser(pleioin, pleiores, sumstats, snps, h2, rg, traits, ref) }

main = function(dat){

  for ( i in 1:length(dat$snps)){

    pleioin = dat$pleioin
    pleiores = dat$pleiores
    snp_reference = dat$ref
    sumstats = dat$sumstats
    rg_matrix = dat$rg
    h2 = dat$h2
    outf = dat$out
    traits = colnames(dat$rg)

    snp = as.character(dat$snps[i])

    jpeg(file.path(outf,paste(snp,'.jpeg',sep='')), width = 7, height = 7, units = 'in', res = 300)
    par(mar =c(0,0,0,0))
    pleioplot(snp, traits, rg_matrix, sumstats, pleioin, pleiores, h2, snp_reference)
    dev.off()

  }

}


main(dat)
