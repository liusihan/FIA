# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 16:20:10 2018

@author: Yuanyuan Zeng
"""

import os
import re
import vcf as pyvcf
import sys

### main
if __name__ == '__main__':  
    inputvcf=sys.argv[1]
    output=sys.argv[2]
#    genelist_file=sys.argv[3]


vcf_reader = pyvcf.Reader(filename = inputvcf, strict_whitespace = True, encoding = 'utf-8')
#outpu = 
with open(output,'w')as f:
    f.write('\t'.join(["CHROM", "POS", "REF", "ALT","Max_AF", "gnomAD_AF_SAS_WES", "gnomAD_grpmax_WES", "gnomAD_grpmax_WGS", "gnomAD_AF_EAS_WES", "gnomAD_AF_ASJ_WES", "gnomAD_AF_MID_WES", "gnomAD_AF_FIN_WES","gnomAD_AF_NFE_WES","gnomAD_AF_AMI_WES"])+'\n')
#    f.write('\t'.join(["CHROM", "POS", "REF", "ALT", "transcript", "genesymbol", "AF", "impact", "consequence"])+'\n')

for record in vcf_reader:
    if str(record.ALT[0]) == "*":
        continue
#     selectinfo = {}
#     selectinfo["genesymbol"] = record.INFO['VEP_SYMBOL'][0].split("|")
#     selectinfo["transcript"] = record.INFO['VEP_Feature'][0].split("|")
# #    selectinfo["hgvsc"] = record.INFO['VEP_HGVSc'][0].split("|")
# #    selectinfo["hgvsp"] = record.INFO['VEP_HGVSp'][0].split("|")
#     selectinfo["consequence"] = record.INFO['VEP_Consequence'][0].split("|")
# #    selectinfo["domain"] = record.INFO['VEP_DOMAINS'][0].split("|")
#     selectinfo["impact"] = record.INFO['VEP_IMPACT'][0].split("|")

    af = 0
    aflist = []
#    aflist = []
    # if record.CHROM == "1" and record.POS == 40773149:
    #     for i in ["gnomAD_AF", "gnomAD_AF_AFR", "gnomAD_AF_AMR", "gnomAD_AF_ASJ", "gnomAD_AF_EAS", "gnomAD_AF_FIN", "gnomAD_AF_NFE","gnomAD_AF_SAS"]:
    #         if i in record.INFO and (not record.INFO[i][0]):
    #             print(i, record.INFO[i][0])
    #             break

    #for i in ["AF_control","gnomAD_AF", "gnomAD_AF_AFR", "gnomAD_AF_AMR", "gnomAD_AF_ASJ", "gnomAD_AF_EAS", "gnomAD_AF_FIN", "gnomAD_AF_NFE","gnomAD_AF_SAS"]:
    for i in ["gnomAD_AF_SAS_WES", "gnomAD_grpmax_WES", "gnomAD_grpmax_WGS", "gnomAD_AF_EAS_WES", "gnomAD_AF_ASJ_WES", "gnomAD_AF_MID_WES", "gnomAD_AF_FIN_WES","gnomAD_AF_NFE_WES","gnomAD_AF_AMI_WES"]:

        if i in record.INFO and record.INFO[i][0]:
            aflist.append(str(record.INFO[i][0]))
            af_temp = float(record.INFO[i][0])
            #aflist.append(i+"="+str(af_temp))
            if(af_temp > af):
                af = af_temp
        else:
            aflist.append("0")
#    print(aflist)
    with open(output,'a')as f:
        f.write('\t'.join([record.CHROM, str(record.POS), record.REF, str(record.ALT[0]), str(af)]) + '\t' + '\t'.join(aflist)+'\n')
    
    # if af <= 0.01:
    #     #print(aflist)
    #     for i in range(len(selectinfo["transcript"])):
    #         if selectinfo["genesymbol"][i] in gene_list:
    #             if re.match("splice",selectinfo["consequence"][i]) or selectinfo["impact"][i] == "HIGH" or selectinfo["impact"][i] == "MODERATE":
    #                 with open(output,'a')as f:
    #                     f.write('\t'.join([record.CHROM, str(record.POS), record.REF, str(record.ALT[0]), selectinfo["transcript"][i], selectinfo["genesymbol"][i], str(af), selectinfo["impact"][i], selectinfo["consequence"][i]])+'\n')


