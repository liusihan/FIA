library(ggplot2)
library(reshape2)
library(grid)
library(RColorBrewer)
library(dplyr)

##################

args <- commandArgs(TRUE)

workdir = args[1] #result/stat/
inputdir = args[2] # dir of dat.RData
p_info = args[3] #pub.txt

load(paste(inputdir,"/dat.RData",sep =""))
dat <- get("dat")
setwd(workdir)

## 1. max af
af <- dat[,c(2:5, 9:16)]
af <- distinct(af, chrom,pos,ref,alt, .keep_all = T)
ans <- t(apply(af, 1, function(x){
    gnomad_af <- x[6:12]
    gnomad_af[which(gnomad_af == ".")] <- 0
    gnomad_af <- as.numeric(gnomad_af,na.rm = T)
    c(x[1:5], max(gnomad_af))
}))
af_max <- as.data.frame(ans, stringsAsFactors = F)
af_max$pos <- as.numeric(af_max$pos)
af_max$AF_control <- as.numeric(af_max$AF_control)
af_max$V6 <- as.numeric(af_max$V6)
colnames(af_max) <- c(colnames(af_max)[1:5],"max_gnomad_af")
af_max <- mutate(af_max, varid = paste(chrom, pos, ref, alt, sep = "_"))

## 2. impact
#vep <- read.delim(featureinfo, head = T, stringsAsFactors = F)
#vep <- vep %>% mutate(varid = paste(CHROM, POS, REF, ALT, sep = "_"))
#vep <- vep %>% distinct(varid, transcript, .keep_all = TRUE)
#dat <- left_join(dat, vep[,-c(1:4,6:8)])
dat <- dat %>% filter(impact == "HIGH" | impact == "MODERATE" | grepl("splice_region",consequence) | grepl("extended_intronic_splice_region_variant",splice_region))

## 3. diagnosed P/B/V
b_sites = dat %>% filter(AR_classify == "Likely benign" | AR_classify == "Benign")
b_sites = unique(b_sites$varid)

diagnosed_dat <- read.delim(p_info, head = T, stringsAsFactors = F)
colnames(diagnosed_dat) <- tolower(colnames(diagnosed_dat))
#colnames(diagnosed_dat) <- c("chrom", colnames(diagnosed_dat)[-1])
diagnosed_dat <- mutate(diagnosed_dat, varid = paste(chrom, pos, ref, alt, sep = "_"))
diagnosed_dat <- distinct(diagnosed_dat, varid, .keep_all = T)

## 4. melt
p_var <- dat[c(6,30:32,34:59)]
#p_var <- p_var[,-5]
p_var_melt <- melt(p_var, id.vars = c("varid","transcript","consequence","impact"),variable.name = "evidence", na.rm = T)
p_var_melt$evidence <- as.character(p_var_melt$evidence)
#p_var_melt <- p_var_melt %>% distinct(varid, transcript, evidence, .keep_all = TRUE)
## 5. merge af,pathogenicity
union_evidence <- left_join(p_var_melt, af_max[,6:7])
union_evidence <- union_evidence %>% mutate(maxaf_break = cut(max_gnomad_af,breaks = c(-Inf,0, 0.0001,0.005,0.01,1)))
union_evidence <- mutate(union_evidence, diagnosis = if_else(varid %in% diagnosed_dat$varid, "P", "NonP"))
union_evidence$pathogenicity <- "V"
union_evidence <- union_evidence %>% mutate(pathogenicity = ifelse(varid %in% b_sites, "B", pathogenicity))
union_evidence <- union_evidence %>% mutate(pathogenicity = ifelse(varid %in% diagnosed_dat$varid, "P", pathogenicity))

# pvb_info <- data.frame(varid = as.character(unique(union_evidence$varid)))
# pvb_info <- data.frame(pvb_info, pathogenicity = rep())
save(union_evidence,b_sites, diagnosed_dat, file = "melted_evidence.RData")

#######  stat  #######
source("/public/home/swgenetics_3/zengyuanyuan/Projects/ACMG/intervar/pipeline_test20181107/scripts/APpipeline/stat/evidence_stat.r")


dat1 = filter(union_evidence, grepl("missense_variant", consequence))
dat2 = filter(union_evidence, impact == "HIGH")
source("/public/home/swgenetics_3/zengyuanyuan/Projects/ACMG/intervar/pipeline_test20181107/scripts/APpipeline/stat/function_forstat.r")
ans <- generate_table(dat1,impact = "missense")
ans <- generate_table(dat2,impact = "impactHIGH")

###### for random forest and logistic regression
#==============================================================
## evidence association analysis  
## all variants
evidence_dcast <- union_evidence
evidence_dcast$value <- 1
evidence_dcast <- dcast(evidence_dcast, varid+diagnosis~evidence)
evidence_dcast[is.na(evidence_dcast)] <- 0
write.table(evidence_dcast, file = "ACMG_evidence_matrix_all_variants.txt", sep = "\t", col.names = T, row.names = F, quote = F)

## missense variants
evidence_dast2 <- filter(union_evidence, grepl("missense_variant", consequence))
evidence_dast2 <- dcast(evidence_dast2, varid+diagnosis~evidence)
temp <- evidence_dast2[,3:ncol(evidence_dast2)]
temp[!is.na(temp)] <- 1
temp[is.na(temp)] <- 0
evidence_dast2 <- data.frame(evidence_dast2[,1:2], temp)
write.table(evidence_dast2, file = "ACMG_evidence_matrix_missense_variants.txt", sep = "\t", col.names = T, row.names = F, quote = F)

#========with evidence strength====================
## evidence association analysis  
## all variants
evidence_dcast <- union_evidence
#evidence_dcast$value <- 1
evidence_dcast <- dcast(evidence_dcast, varid+diagnosis~evidence)
evidence_dcast[is.na(evidence_dcast)] <- 999
write.table(evidence_dcast, file = "ACMG_evidence_matrix_all_variants.strength.txt", sep = "\t", col.names = T, row.names = F, quote = F)

## missense variants
evidence_dast2 <- filter(union_evidence, grepl("missense_variant", consequence))
evidence_dast2 <- dcast(evidence_dast2, varid+diagnosis~evidence)
temp <- evidence_dast2[,3:ncol(evidence_dast2)]
#temp[!is.na(temp)] <- 1
temp[is.na(temp)] <- 999
evidence_dast2 <- data.frame(evidence_dast2[,1:2], temp)
write.table(evidence_dast2, file = "ACMG_evidence_matrix_missense_variants.strength.txt", sep = "\t", col.names = T, row.names = F, quote = F)




