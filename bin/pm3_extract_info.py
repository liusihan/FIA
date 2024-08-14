#coding:utf-8
"""按照VEP_Feature那一列拆分转录本，并保留'MODERATE' or HIGH 的位点生成新的文件
"""

"""[1]CHROM        [2]POS  [3]REF  [4]ALT  [5]VEP_Feature  
[6]VEP_IMPACT   [7]VEP_PICK     [8]VEP_Gene     [9]VEP_SYMBOL   [10]AH-001:GT  

"""
import sys
import re
input=sys.argv[1]#testdata_split.tsv
input_p=sys.argv[2]#20181206_z_AR_case_p.sites
linkage=sys.argv[3]#remove.txt
output=sys.argv[4]
inputgene=sys.argv[5]#gene_symbol


#outputfile=open(output, "w")
#read p sites
"""chrom   pos     ref     alt
10      102789780       C       A
10      102789810       C       CG

"""

if __name__ == '__main__':
	#read gene_symbol 对应到VEP_Symbol
	with open(inputgene,'r')as f:
		linenum=0
		gene_symbol=[]
		for line in f:
			head=line.strip().split('\t')[0]
			gene_symbol.append(head)
			
	with open(input_p,'r')as f:
		linenum=0
		snp_p=[]
		for line in f:
			linenum+=1
			if(linenum==1):
				continue
			else:
				head=line.strip().split('\t')[0:4]
				snp_p.append('_'.join(head))
	#read likage sites
	"""sampleid        chr     pos     ref     alt
	AH-002  17      18055194        G       A
	AH-048  17      18045546        C       T
	"""
	with open(linkage,'r')as f:
		linenum=0
		link={}
		for line in f:
			linenum+=1
			if(linenum==1):
				continue
			else:
				head=line.strip().split('\t')
				snp = '_'.join(head[1:5])
				if snp not  in link:
					link[snp] = [head[0]]
				else:
					link[snp].append(head[0])
		

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
					f.write('\t'.join(temp[0:8]) + '\t'+'\t'.join(header[19:])+'\n')
				
			else:
				sampleIDs = header[19:]
				tmp = line.strip().split('\t')
				#genotypes=[i.split(";")[0] for i in tmp[17:]]
				genotypes=tmp[19:]
				snp = "_".join(tmp[0:4])
				## 标记连锁位点
				if snp in link.keys():
					#print(snp)
					for i in range(len(sampleIDs)):
						sid = sampleIDs[i]
						if sid in link[snp]:
							#print(genotypes[i])
							#print(sid)
							genotypes[i] = 'rm'
							#print(genotypes[i])
				## 拆 Feature,判断基因
				VEP_Feature=tmp[4].split('|')
				VEP_IMPACT=tmp[5].split('|')
				VEP_SYMBOL=tmp[6].split('|')
				#VEP_Gene=tmp[7].split('|')
				VEP_splice=tmp[17].split('|')
				VEP_consequence=tmp[18].split('|')
				#VEP_IMPACT=tmp[5].split('|')
				if snp in snp_p:
					for i in range(len(VEP_Feature)):
						if VEP_SYMBOL[i] in gene_symbol:
							#print(snp)
						#print('\t'.join(tmp[0:4]) + '\t' + VEP_Feature[i] + '\t' +'\t'.join(genotypes))
							#with open(output, "a") as f:
							key='\t'.join(tmp[0:4]) + '\t' + VEP_Feature[i] + '\t'+ VEP_IMPACT[i] + '\t'+ VEP_SYMBOL[i] + '\t'+VEP_SYMBOL[i] + '\t'+'\t'.join(genotypes)
							with open(output, 'a')as f:
								f.write(key+'\n')

							#out[key]=VEP_IMPACT[i] + '\t'+ VEP_SYMBOL[i] + '\t'+VEP_Gene[i] + '\t'+'\t'.join(genotypes)

								#f.write('\t'.join(tmp[0:4]) + '\t' + VEP_Feature[i] + '\t' +VEP_IMPACT[i] + '\t' +VEP_SYMBOL[i] + '\t'+VEP_Gene[i] + '\t'+'\t'.join(genotypes)+'\n')
				else:
					af = tmp[7:17]    ### 判断频率
					filter_flag = 1
					for i in af:
						if i == ".":
							continue
						elif float(i) > 0.005:#任意一个频率＞0.005都过滤掉
							filter_flag = 0
							break
					if filter_flag==1:
						for i in range(len(VEP_Feature)):
							if (VEP_SYMBOL[i] in gene_symbol) and (VEP_IMPACT[i]=='MODERATE' or VEP_IMPACT[i]=='HIGH' or re.findall(r"extended_intronic_splice_region_variant",VEP_splice[i]) or re.findall(r"splice_region",VEP_consequence[i]) ):
							# and VEP_IMPACT[i]=='MODERATE' or VEP_IMPACT[i]=='HIGH':
								key='\t'.join(tmp[0:4]) + '\t' + VEP_Feature[i] + '\t' + VEP_IMPACT[i] + '\t' +VEP_SYMBOL[i] + '\t'+VEP_SYMBOL[i] + '\t'+'\t'.join(genotypes)
								with open(output, 'a')as f:
									f.write(key+'\n')
							#	out[key]=VEP_IMPACT[i] + '\t' +VEP_SYMBOL[i] + '\t'+VEP_Gene[i] + '\t'+'\t'.join(genotypes)

							# with open(output, "a") as f:
								# #print('\t'.join(tmp[0:4]) + '\t' + VEP_Feature[i] + '\t' + '\t'.join(tmp[15:])+'\n')
								# f.write('\t'.join(tmp[0:4]) + '\t' + VEP_Feature[i] + '\t' +VEP_IMPACT[i] + '\t' +VEP_SYMBOL[i]  + '\t'+VEP_Gene[i] + '\t'+ '\t'.join(genotypes)+'\n')
	# for i in out.keys():
	# 	tmp1=str(i)
	# 	#print(tmp1)
	# 	tmp2=str(out[i])
	# 	#print(tmp1,tmp2)
	# 	with open(output,'a')as f:
	# 		f.write(tmp1+'\t'+tmp2+'\n')
