# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 16:20:10 2018

@author: Yuanyuan Zeng
"""

import os
import re
import vcf as pyvcf
import sys

#os.chdir('D:\\工作\\ACMG\\scripts')

#### data
## P 位点信息存入字典
## 按行读入VCF
## 拆分转录本、HGVSp蛋白信息、Consequence信息 存为list

# P_file = open("Public_db_20181026.txt")
# #transcript_length = open("ensemble_transcript_length.txt")
# transcript_length = open("ensemble_transcript_length_20181113.txt")
# link_variants = open("linkage_in_samples_20181015.txt")
# dvpred_file = open("test.dvpred.tsv")
#vcf = open("test_pm5_ps1.vcf")
#vcf = open("test_samples.vcf")
#vcf = open("")

### main
if __name__ == '__main__':  
    inputvcf=sys.argv[1]
    output=sys.argv[2]


vcf_reader = pyvcf.Reader(filename = inputvcf, strict_whitespace = True, encoding = 'utf-8')
#outpu = 
with open(output,'w')as f:
    f.write('\t'.join(["CHROM", "POS", "REF", "ALT", "transcript", "genesymbol","hgvsc","hgvsp", "impact", "consequence","splice_region"])+'\n')



########
#flag = {}

for record in vcf_reader:
    if str(record.ALT[0]) == "*":
        continue
    selectinfo = {}
    selectinfo = {}
    for i in ['VEP_SYMBOL','VEP_Feature','VEP_HGVSc','VEP_HGVSp','VEP_SpliceRegion','VEP_Consequence','VEP_IMPACT']:
        if isinstance(record.INFO[i], list):
            selectinfo[i] = record.INFO[i][0].split("|")
        else:
            selectinfo[i] = ['']
    # selectinfo["genesymbol"] = record.INFO['VEP_SYMBOL'][0].split("|")
    # selectinfo["transcript"] = record.INFO['VEP_Feature'][0].split("|")
    # selectinfo["hgvsc"] = record.INFO['VEP_HGVSc'][0].split("|")
    # selectinfo["hgvsp"] = record.INFO['VEP_HGVSp'][0].split("|")
    # selectinfo["consequence"] = record.INFO['VEP_Consequence'][0].split("|")
    # selectinfo["domain"] = record.INFO['VEP_DOMAINS'][0].split("|")
    # selectinfo["impact"] = record.INFO['VEP_IMPACT'][0].split("|")
    
    for i in range(len(selectinfo["VEP_Feature"])):
        with open(output, 'a')as f:
            f.write('\t'.join([record.CHROM, str(record.POS), record.REF, str(record.ALT[0]), selectinfo["VEP_Feature"][i], selectinfo["VEP_SYMBOL"][i], selectinfo['VEP_HGVSc'][i], selectinfo['VEP_HGVSp'][i], selectinfo["VEP_IMPACT"][i],selectinfo["VEP_Consequence"][i], selectinfo["VEP_SpliceRegion"][i]])+'\n')



