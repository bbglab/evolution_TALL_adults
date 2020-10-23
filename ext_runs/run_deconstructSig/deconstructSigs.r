#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## load libraries
library(deconstructSigs)
library(gplots)

# read our input file
ttype=args[1]
x<-read.table(paste(ttype, sep=""), header=T, sep="\t", colClasses=c(NA,NA,NA,NA,NA))
colnames(x) <- c("chr", "pos", "ref", "alt", "Sample")
x$chr <- as.factor(x$chr)
head(x)

# read signatures file 
path_signatures = args[2]
signatures_used <- read.table(path_signatures, row.names=1,sep="\t", header=T)
cols.cosmic <- colnames(signatures.cosmic)
colnames(signatures_used) <- cols.cosmic
print(head(signatures_used))

# subset to use
subset = args[3]

if (subset == 'BALL_subset') { 
  list_blood_sig <- c("SBS1","SBS5","SBS18", "SBS2", "SBS13","SBS36", "SBS37", "SBS9") # run_subset_28012020 and run_subset_30012020
} else if (subset == 'TALL_subset') {
  list_blood_sig <- c("SBS1","SBS5","SBS18") # run_21012019 
} else if (subset == 'ALL_primary_subset') {
  list_blood_sig <- c("SBS1", "SBS2","SBS5", "SBS6","SBS9", "SBS13", "SBS17a", "SBS17b","SBS18","SBS34","SBS36", "SBS37")
} else if  (subset == 'TALL_HSCP_comparative') {
  list_blood_sig <- c("SBS1","SBS5","SBS18", "SBS_hscp") # run run_subsets_hemato
} else if (subset == 'TALL_relapse_subset') {
  list_blood_sig <- c("SBS1","SBS5", "SBS18","SBS32", "SBSA_new", "SBSB_new") #relapse samples #"SBSA_new", "SBSB_new # blood paper
}

# list of signatures from blood tumors of pcawg without the ones that Maura et al. 2019 do not recomend 
#list_blood_sig <- c("SBS1", "SBS2","SBS5", "SBS6","SBS9", "SBS13", "SBS17a", "SBS17b","SBS18","SBS34",
#                     "SBS36", "SBS37") #BALL

# list of signatures that had mutations asigned after first run
#list_blood_sig <- c("SBS1","SBS5","SBS18", "SBS2", "SBS13", "SBS37", "SBS9") #BALL
#list_blood_sig <- c("SBS1","SBS5","SBS18", "SBS37") #TALL

signatures_used_blood <- signatures_used[list_blood_sig,]
# signatures.pcawg <- read.table("/workspace/projects/all_aecc/scripts_general/SigProfiler_SBS_signatures_2018_03_28.deconstructsigs.tsv",
#                          row.names=1,sep="\t", header=T)
# 


# Convert to deconstructSigs input
# this step generates the matrix suitable for the program
sigs.input <- mut.to.sigs.input(mut.ref = x, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")
write.table(sigs.input, file="mut_count_96_ch.tsv", sep="\t", col.names = TRUE, row.names=TRUE)
# now we run deconstructSigs for each sample in our input list
flag = 0
all_tests <- list()
for (sample in unique(x$Sample))
{

    if (nrow(x[which(x$Sample==sample),]) > 30) 
    {
        test = whichSignatures(tumor.ref = sigs.input, 
                               signatures.ref = signatures_used_blood, 
                               sample.id = sample,
                               contexts.needed = TRUE,
                               tri.counts.method = 'default',
                               signature.cutoff = 0.1)
        all_tests[[sample]] <- test
        a = test$weights # save the weights for each signature. 
        a['SSE']  = round(sqrt(sum(test$diff * test$diff)), digits = 3) # compute the error rate
        a['unknown'] = test$unknown 
        a['mutation_count'] = nrow(x[which(x$Sample==sample),]) # number of mutations
        # append the results of each sample in to dataframe
        if (flag == 0){total = a; flag=1}
        else{total <- rbind(total, a)}
    }
}

save(all_tests, file = "all_tests.RData")
# prepare heatmap

if (endsWith(path_signatures, "cosmic")) {
  if (length(unique(x$Sample)) != 1) {
      new = as.matrix(total[ , grepl( "Signature" , names( total ) ) ])
      pdf(paste("heatmap_deconstruct", ".pdf", sep=""))
      heatmap.2(new,dendrogram='none', Rowv=TRUE, Colv=TRUE,trace='none', cexRow=0.5, cexCol=0.5)
      dev.off()
  }
}
if (endsWith(path_signatures, "pcawg")) {
  if (length(unique(x$Sample)) != 1) {
    new = as.matrix(total[ , grepl( "SBS" , names( total ) ) ])
    pdf(paste("heatmap_deconstruct", ".pdf", sep=""))
    heatmap.2(new,dendrogram='none', Rowv=TRUE, Colv=TRUE,trace='none', cexRow=0.5, cexCol=0.5)
    dev.off()
  }
}

# prepare CSV file
myDF <- cbind(sample_id = rownames(total), total) # assign row names
rownames(myDF) <- NULL
# write the output to a file
write.table(myDF, file=paste("signatures_weight", ".csv", sep=""), sep="\t", col.names = TRUE, row.names=FALSE)
