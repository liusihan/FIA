#coding:utf-8
"""按照VEP_Feature那一列拆分转录本，并保留'MODERATE' or HIGH 的位点生成新的文件
"""

"""[1]CHROM        [2]POS  [3]REF  [4]ALT  [5]VEP_Feature  
[6]VEP_IMPACT   [7]VEP_PICK     [8]VEP_Gene     [9]VEP_SYMBOL   [10]AH-001:GT  

"""
import sys
import re

input=sys.argv[1]
output=sys.argv[2]


if __name__ == '__main__':
	with open(input,'r')as f:
		linenum=0
		header=[]
		out={}
		for line in f:
			linenum+=1
			if(linenum==1):
				temp=line.strip().split('\t')
				header = [i.split("]")[1].split(":")[0] for i in temp]
				with open(output,'w')as f:
					f.write('\t'.join(temp[0:8]) + '\t'+'\t'.join(header[17:])+'\n')			
			else:
				#sampleIDs = header[17:]
				tmp = line.strip().split('\t')
				genotypes=tmp[17:]
				#snp = "_".join(tmp[0:4])
				##判断intergenic_variant位点
				if tmp[4] == 1:
					continue
				## 拆 Feature
				VEP_Feature=tmp[4].split('|')
				VEP_IMPACT=tmp[5].split('|')
				VEP_SYMBOL=tmp[6].split('|')
				VEP_Gene=tmp[7].split('|')
				#VEP_IMPACT=tmp[5].split('|')
				for i in range(len(VEP_Feature)):
					if re.findall(r"ENS", VEP_Feature[i]):
						continue
					key='\t'.join(tmp[0:4]) + '\t' + VEP_Feature[i] + '\t'+ VEP_IMPACT[i] + '\t'+ VEP_SYMBOL[i] + '\t'+VEP_Gene[i] + '\t'+'\t'.join(genotypes)
					with open(output, 'a')as f:
						f.write(key+'\n')
