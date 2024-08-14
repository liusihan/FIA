import os, re, sys, argparse, configparser, logging, subprocess
from AutoPVS1.autoPVS1 import AutoPVS1
import vcf as pyvcf
description = 'Python script for PVS1 annotation with AutoPVS1'
parser = argparse.ArgumentParser(description = description)
parser.add_argument("--input","-i",required = True,help = "Input of VCF file.")
parser.add_argument("--genome","-g",required = True,help = "Genome version. {hg19,hg38}")
parser.add_argument("--output","-o",required = True,help = "Output directory")
args = parser.parse_args()

vcf_reader = pyvcf.Reader(filename = args.input, strict_whitespace = True, encoding = 'utf-8')

output_dir = os.path.abspath(args.output)

with open(output_dir+"/PVS1.tsv",'w') as f:
        f.write('\t'.join(["chrom","pos","ref","alt","transcript","PVS1"])+'\n')


for record in vcf_reader:
        if str(record.ALT[0]) == "*":
            continue
        transcripts = str(record.INFO["VEP_Feature"][0]).split("|")
        for transcript in transcripts:
                var = "-".join([record.CHROM, str(record.POS), record.REF, str(record.ALT[0])])
                #print(f"'{var}'",f"'{args.genome}'",f"'{transcript}'")
                pvs1_output=AutoPVS1(var,args.genome,user_trans=str(transcript))
                pvs1_strength = "NA"
                if pvs1_output.islof:
                        #print(pvs1_output.pvs1.strength.value)
                        if pvs1_output.pvs1.strength.name == "VeryStrong":
                                pvs1_strength = 0
                        elif pvs1_output.pvs1.strength.name == "Strong":
                                pvs1_strength.name = 1
                        elif pvs1_output.pvs1.strength.name == "Moderate":
                                pvs1_strength = 2
                        elif pvs1_output.pvs1.strength.name == "Supporting":
                                pvs1_strength = 3
                        else:
                                pvs1_strength = "NA"
                with open(output_dir+"/PVS1.tsv", 'a') as f:
                        f.write("\t".join([record.CHROM, str(record.POS), record.REF, str(record.ALT[0]),transcript,str(pvs1_strength)])+"\n")

