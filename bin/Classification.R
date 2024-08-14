library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 2) stop("Incorrect number of arguments, usage> Rscript Classification.R input output"))
input = args[1]
Output = args[2]

Classification<-function(evidence){
    p_evidence<-gsub(pattern = "[B][A-Z]*\\d[=][0-9]", replacement = "", x = evidence)
    #remove evidence for pathogenic
    b_evidence<-gsub(pattern = "[P][A-Z]*\\d[=][0-9]", replacement = "", x = evidence)
    p_VeryStrong=str_count(string = p_evidence, pattern = "=0")
    p_Strong=str_count(string = p_evidence, pattern = "=1")
    p_Moderate=str_count(string = p_evidence, pattern = "=2")
    p_Supporting=str_count(string = p_evidence, pattern = "=3")
    Stand_alone=str_count(string = b_evidence, pattern = "=0")
    b_Strong=str_count(string = b_evidence, pattern = "=1")
    b_Moderate=str_count(string = b_evidence, pattern = "=2")
    b_Supporting=str_count(string = b_evidence, pattern = "=3")
    PVS1=str_count(string = p_evidence, pattern = "PVS1=VS")
    PM2=str_count(string = p_evidence, pattern = "PM2=P")
    p_class=""
    b_class=""
    if(p_VeryStrong>=1 && p_Moderate==1){p_class="LP"}
    if(p_Strong==1 && (p_Moderate==1 || p_Moderate==2)){p_class="LP"}
    if(p_Strong==1 && p_Supporting>=2){p_class="LP"}
    if(p_Moderate>=3){p_class="LP"}
    if(p_Moderate==2 && p_Supporting>=2){p_class="LP"}
    if(p_Moderate==1 && p_Supporting>=4){p_class="LP"}
    if(p_VeryStrong>=1){
      if(p_Strong>=1 || p_Moderate>=2 || p_Supporting>=2 || (p_Moderate==1 && p_Supporting==1)){
        p_class="P"
      }
    }
    if(p_VeryStrong>=2){p_class="P"}
    if(p_Strong>=2){p_class="P"}
    if(p_Strong==1){
      if(p_Moderate>=3 || (p_Moderate==2 && p_Supporting>=2) || (p_Moderate==1 && p_Supporting>=4)){
        p_class="P"
      }
    }
    if((b_Strong==1 && b_Supporting==1) || (b_Supporting>=2) || (b_Strong==1 && b_Moderate==1) || ((b_Supporting>=1 && b_Moderate>=1))){b_class="LB"}
    if(Stand_alone>=1 || b_Strong>=2){b_class="B"}
    if(b_class=="" && p_class==""){return("VUS")}
    else if(b_class!="" && p_class!=""){return("VUS")}
    else if(b_class=="" && p_class!=""){return(p_class)}
    else{return(b_class)}
}

input_result <- read.table(input, head = T, stringsAsFactors = F, sep="\t")

for(i in 1:nrow(input_result)){
    input_result$Classify_AD[i]="VUS"
    input_result$Classify_AR[i]="VUS"
    if(is.na(input_result$AR_evidence[i]) && is.na(input_result$AD_evidence[i])){next}
    if(!is.na(input_result$AR_evidence[i])){input_result$Classify_AR[i]=Classification(input_result$AR_evidence[i])}
    if(!is.na(input_result$AD_evidence[i])){input_result$Classify_AD[i]=Classification(input_result$AD_evidence[i])}
}

write.table(input_result, Output,sep="\t", row.name=FALSE, col.names=TRUE, quote=FALSE)
