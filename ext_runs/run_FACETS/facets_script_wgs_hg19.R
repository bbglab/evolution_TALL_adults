#!/home/isentis/.conda/envs/copynumber/bin/Rscript

##################################### LIBRARIES USED ###############################################

library("optparse")
library("facets")
library(R.utils)

########################## PARSING COMMAND-LINE ARGUMENTS ##########################################

option_list <- list(
  make_option(c("-i", "--input"), 
              action="store_true", 
              dest="input", 
              type="character",
              default=FALSE,
              help="Path with the SNP file created by snp-pileup "),
  make_option(c("-o", "--output"), 
              action="store_true",
              dest="output", 
              type="character",
              default=FALSE,
              help="Path of the output directory, default: working directory"),
  make_option(c("-s", "--seed"),
              action="store_true",
              dest="seed",
              type="integer",
              default=FALSE,
              help="Seed random number for reproducible results"))

opt <- parse_args(OptionParser(option_list=option_list), args = commandArgs(trailingOnly = TRUE))

###################################### RUNNING FACETS ###############################################

in_file = opt$input
out_dir = opt$output
seed_num = as.integer(opt$seed)

set.seed(seed_num)

if (endsWith(in_file, 'gz')){
    print("unzipped")
    data = gunzip(in_file)
} else {
    data = in_file
}

#data = gunzip(paste(in_file, ".gz", sep=""))
#data = in_file
rcmat = readSnpMatrix(data)
xx = preProcSample(rcmat, snp.nbhd = 5000, ndepth=5, cval=75, gbuild="hg19")
oo = procSample(xx, cval=800, min.nhet = 25)
fit = emcncf(oo)

file.create(paste(out_dir, 'metrics.txt' , sep="/"))
write(paste("purity", fit$purity, sep=':'), file=paste(out_dir, 'metrics.txt', sep="/"), append=TRUE)
write(paste("ploidy",fit$ploidy, sep=':'), file=paste(out_dir, 'metrics.txt', sep="/"), append=TRUE)

setwd(out_dir)
jpeg("plot_sample.jpg")
plotSample(x=oo, emfit=fit)
dev.off()
write.table(fit$cncf, file= "results_x_seg",sep="\t", col.names=TRUE)


fit2 = emcncf2(oo)
write.table(fit2$cncf, file= "results_global",sep="\t", col.names=TRUE)

sessionInfo()
