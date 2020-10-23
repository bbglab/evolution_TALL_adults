# This piece of code relies on a workspace directory structure such as
#
# cohort/
#     patientID/
#         DxTumorID_vs_normalID/
#         ReTumorID_vs_normalID/ (sometimes)
# patientID, DxTumorID etc can be found in ../ext_files/all_cohort_clinical_groups.tsv
#
# Be aware that the filtered mafs with clonal classification and joined mutations after running the scripts in filter/
# have the following file name: TumorID_vs_normalID + _strelka_uniq_all_anno_vep92_categories_filt_snps_cluster.maf.
# This file name is used in the following code.

c_path = "" # path to the adult T-ALL cohort folder

obs_ccf_batch1 <- list()
obs_ccf_batch1 <- lapply(c("PAT1","PAT2", "PAT3", "PAT4", "PAT5", "PAT6", "PAT7", "PAT8", "PAT9", "PAT10",
                    "PAT11", "PAT12", "PAT13", "PAT14", "PAT15"), function (c_pat){

  print(c_pat)
  c_files <- list.files(path = paste(c_path, c_pat, sep = "/") , pattern="_strelka_uniq_all_anno_vep92_categories_filt_snps_cluster.maf", recursive=TRUE, full.names=TRUE)
  if(length(c_files) != 2) return(NA)
  c_pri <- read.table(c_files[1], header = TRUE, sep = "\t", comment.char = "" )
  c_rel <- read.table(c_files[2], header = TRUE, sep = "\t", comment.char = "" )

  c_rel <- c_rel[c_rel$ccf > 0.9,]
  c_rel <- c_rel[,c("X.CHROM", "REF", "POS", "ALT")]
  c_pri <- c_pri[,c("X.CHROM", "REF", "POS", "ALT", "ccf")]
  c_rel <- merge(c_rel, c_pri, all.x = TRUE)
  c_rel
})
names(obs_ccf_batch1) <- c("PAT1","PAT2", "PAT3", "PAT4", "PAT5", "PAT6", "PAT7", "PAT8", "PAT9", "PAT10",
                    "PAT11", "PAT12", "PAT13", "PAT14", "PAT15")


c_path = "" # path to the adult T-ALL cohort folder (this was a different batch that's whay there are two paths to the mafs)

obs_ccf_batch2 <- list()
obs_ccf_batch2<- lapply(c("PAT16", "PAT17", "PAT18", "PAT19"), function (c_pat){

  print(c_pat)
  c_files <- list.files(path = paste(c_path, c_pat, sep = "/"), pattern="_strelka_uniq_all_anno_vep92_categories_filt_snps_cluster.maf", recursive=TRUE, full.names=TRUE)
  if(length(c_files) != 2) return(NA)
  c_pri <- read.table(c_files[1], header = TRUE, sep = "\t", comment.char = "" )
  c_rel <- read.table(c_files[2], header = TRUE, sep = "\t", comment.char = "" )

  c_rel <- c_rel[c_rel$ccf > 0.9,]
  c_rel <- c_rel[,c("X.CHROM", "REF", "POS", "ALT")]
  c_pri <- c_pri[,c("X.CHROM", "REF", "POS", "ALT", "ccf")]
  c_rel <- merge(c_rel, c_pri, all.x = TRUE)
  c_rel
})
names(obs_ccf_batch2) <- c("PAT16", "PAT17", "PAT18", "PAT19")

obs_ccf <- c(obs_ccf_batch1, obs_ccf_batch2)


pdf(file = "relapse_fixation.pdf", width = 9, height = 5)

c_lty <- c("solid", "longdash", "dotted", "dotdash")
c_color <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
par(mfrow = c(1,1), mar = c(4,4,4,4), bty ="l", las = 1)
c_type <- 0


plot(NA, xlim = c(0,1.1), ylim = c(0,1.1), xlab = "Primary CCF bins",ylab = "Frequency",
     main = "Fixed mutations (CCF > 90%) at relapse", axes = FALSE)
for (i in names(obs_ccf)){
  obs_ccf[[i]]$ccf[is.na(obs_ccf[[i]]$ccf)] <- runif(sum(is.na(obs_ccf[[i]]$ccf)), 0, 0.1)
  obs_ccf[[i]]$ccf[obs_ccf[[i]]$ccf > 0.9] <- 0.9
  obs_ccf[[i]]$ccf <- floor(obs_ccf[[i]]$ccf * 10)/10
  obs_ccf[[i]]$ccf <- factor(obs_ccf[[i]]$ccf, levels = seq(0, max(obs_ccf[[i]]$ccf), 0.1))
  c_table <- table(obs_ccf[[i]]$ccf)
  lines(x = as.numeric(names(c_table)), y = c_table/ sum(c_table),
        col = c_color[(c_type %% 6) + 1], lty = c_lty[(c_type %/% 6)+1])
  points(x = as.numeric(names(c_table)), y = c_table/ sum(c_table),
         pch = 16, col = c_color[(c_type %% 6) + 1], cex = 0.8)
  lines(x= c(0.2 + (c_type %/% 7)/7, 0.25 + (c_type %/% 7)/7), y= c(0.9 - (c_type %% 7) /20, 0.9 - (c_type %% 7)/20),
        col = c_color[(c_type %% 6) + 1], lty = c_lty[(c_type %/% 6)+1])
  points(x= 0.225 + (c_type %/% 7)/7, y = 0.9 - (c_type %% 7) /20,
         pch = 16, col = c_color[(c_type %% 6) + 1], cex = 0.8)
  text(x= 0.25 + (c_type %/% 7)/7, y = 0.9 - (c_type %% 7) /20, labels = i, cex = 0.8, pos = 4)
  c_type <- c_type + 1
}
axis(1, at = seq(0,.9, .1), labels = paste(seq(0,.9,.1), seq(.1,1,.1), sep = "-"))
axis(2, at = seq(0,1, 0.2), labels = seq(0,1, 0.2), las = 1)

dev.off()

