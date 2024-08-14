#=====================  main  =================================================
library(parallel)

### main


args <- commandArgs(TRUE)

evidence_db_both_file = args[1]
evidence_db_ar_file = args[2]
evidence_db_ad_file = args[3]
evidence_db_other_dir = args[4]
#evidence_sample_list = args[5]
# evidence_sample_file = args[5]
# outdir <- dirname(evidence_sample_file)
# out_prefix <- gsub(".sample_evidence.tsv","", basename(evidence_sample_file))
outdir <- gsub("/w+$","", dirname(evidence_db_both_file))
outdir <- paste(outdir,"/data/",sep = "")
out_prefix <- "nonsample_evidence"

# workdir = args[1]
# inputdir = args[2]
# out_prefix = args[3]
# PP5 = args[4]

setwd(outdir)

### 1 all evidence
absolute_path <- normalizePath(evidence_db_ar_file, mustWork = FALSE)
file_dir <- dirname(absolute_path)

all_pre_temp <- read.delim(paste(file_dir,"/all_prediction_20181105.txt", sep=""),head = F, stringsAsFactors = F)      ####
rownames(all_pre_temp) <- all_pre_temp[,2]
all_pre <- unlist(all_pre_temp[,1])
names(all_pre) <- all_pre_temp[,2]

### 2 

    evidence_db_both <- read.delim(evidence_db_both_file, head = T, stringsAsFactors = F)      ####
    evidence_db_both$id <- paste(evidence_db_both$chrom, evidence_db_both$pos, evidence_db_both$ref, evidence_db_both$alt, evidence_db_both$transcript,sep = "_")
    evidence_db_both <- evidence_db_both[!duplicated(evidence_db_both$id),]

    evidence_db_ad <- read.delim(evidence_db_ad_file, head = T, stringsAsFactors = F)      ####
    evidence_db_ad$id <- paste(evidence_db_ad$chrom, evidence_db_ad$pos, evidence_db_ad$ref, evidence_db_ad$alt, evidence_db_ad$transcript,sep = "_")
    evidence_db_ad <- evidence_db_ad[!duplicated(evidence_db_ad$id),]

    evidence_db_ar <- read.delim(evidence_db_ar_file, head = T, stringsAsFactors = F)      ####
    evidence_db_ar$id <- paste(evidence_db_ar$chrom, evidence_db_ar$pos, evidence_db_ar$ref, evidence_db_ar$alt, evidence_db_ar$transcript,sep = "_")
    evidence_db_ar <- evidence_db_ar[!duplicated(evidence_db_ar$id),]


    n <- which(colnames(evidence_db_both) %in% c(all_pre_temp[,2] , "id"))
    evidence_db_both <- evidence_db_both[,c(1:5,n)]
    ##
    file_list <- list.files(path = evidence_db_other_dir, pattern = ".t[sx][tv]$",full.names = T)   
    bs2_ar_file <- file_list[grep("BS2", file_list)]
    file_list <- file_list[-grep("BS2", file_list)]
    for(i in 1:length(file_list)){
        dat <- read.delim(file_list[i], head = T, stringsAsFactors = F)
        dat$id <- paste(dat[,1], dat[,2], dat[,3], dat[,4], dat[,5], sep = "_")
        dat <- dat[!duplicated(dat$id),]
        evidence_db_both <- merge(evidence_db_both, dat[,-c(1:5)], by = 'id', all.x = T, sort = F)
    }

#    dat_merge_ar <- merge(evidence_db_ar, evidence_db_ar_pm3, by = c("chrom", "pos", "ref", "alt", "transcript"), all.x = T)    

    dat_merge_ar <- merge(evidence_db_both, evidence_db_ar, by = 'id', all.x = T, sort = F)
    bs2_ar <- read.delim(bs2_ar_file, head = T, stringsAsFactors = F)
    bs2_ar$id <- paste(bs2_ar[,1], bs2_ar[,2], bs2_ar[,3], bs2_ar[,4], bs2_ar[,5], sep = "_")
    bs2_ar <- bs2_ar[!duplicated(bs2_ar$id),]
    dat_merge_ar <- merge(dat_merge_ar, bs2_ar, by = 'id', all.x = T, sort = F)

    dat_merge_ad <- merge(evidence_db_both, evidence_db_ad, by = 'id', all.x = T, sort = F) 
    if(!all(dat_merge_ad$id == dat_merge_ar$id)){
        n <- match(dat_merge_ar$id, dat_merge_ad$id)
        dat_merge_ad <- dat_merge_ad[n,]
    }

    n1 <- which(colnames(dat_merge_ar) %in% all_pre_temp[grep("AR", all_pre_temp[,3]),2])
    n2 <- which(colnames(dat_merge_ad) %in% all_pre_temp[grep("AD", all_pre_temp[,3]),2])
    dat_merge_ar <- dat_merge_ar[,c(2:6,n1)]    
    dat_merge_ad <- dat_merge_ad[,c(2:6,n2)]
    colnames(dat_merge_ar) <- gsub(".x", "",colnames(dat_merge_ar))
    colnames(dat_merge_ad) <- gsub(".x", "",colnames(dat_merge_ad))

    ## PVS1 and PP3, only PVS1 keep
    n1 <- which(!is.na(dat_merge_ar[,"PP3"]) & !is.na(dat_merge_ar[,"PVS1"]))
    dat_merge_ar[n1,"PP3"] <- NA
    n2 <- which(!is.na(dat_merge_ad[,"PP3"]) & !is.na(dat_merge_ad[,"PVS1"]))
    dat_merge_ad[n2,"PP3"] <- NA
    save(dat_merge_ad, dat_merge_ar, file = paste(outdir,"/db_ad_ar.RData",sep = ""))
    write.table(dat_merge_ad, file = paste(outdir,"/test_dat_ad.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
    write.table(dat_merge_ar, file = paste(outdir,"/test_dat_ar.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)









