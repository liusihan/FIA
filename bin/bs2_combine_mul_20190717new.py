#coding:utf-8
import numpy as np
import random
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import os
import multiprocessing
import sys
sns.set_style('whitegrid')


#inputdir='/public/home/swgenetics_3/zengyuanyuan/Projects/ACMG/intervar/pipeline_test20181107/scripts/ACMG_hmm/scripts/data/network.tsv'
#inputtrans='/public/home/swgenetics_5/hmm/Rawdata/scripts/ACMG_v3/1108uniqe_trans.tsv'
#inputz='/public/home/swgenetics_3/zengyuanyuan/Projects/ACMG/intervar/pipeline_test20181107/scripts/ACMG_hmm/scripts/data/input/20181206_z_AR_case_p.sites'
#output='/public/home/swgenetics_3/zengyuanyuan/Projects/ACMG/intervar/pipeline_test20181107/scripts/ACMG_hmm/scripts/data/pm3_result.tsv'
#os.chdir("D:\\工作\\ACMG\\scripts\\ba2")

def featureTypes():
#split trans
	featureTypes=[]
	news_ids = []
	with open(inputtrans, 'r') as f:
		lineNum = 0
		for line in f:
			tmp = line.replace('\n','').split('\t')
			featureTypes.append(tmp[0])
		return featureTypes		
def readp(inputz):
	with open(inputz,'r') as f:
		sites=[]
		num=0
		for line in f:
			num+=1
			if(num==1):
				continue
			tmp=line.replace('\n', '').split('\t')
			sites.append(tmp[0] + '_' + tmp[1] + '_' + tmp[2] + '_' + tmp[3])
	return sites
		
def readTable(feature,sites,inputdir):
	"""compute nodeList, sampleIDs, totalList one gene
		input sourcefile:gene_source(include head)
		return snp sampleIDs genoTypes 
	"""
	#sites=readp()
	#print(len(sites))
	sampleIDs=[]
	nodeList = []
	tmp0=[]
	totalList = []
	tmp1=[]
	tmpsample=[]
	sampleTmp=[]
	label='NA'
	with open( inputdir, 'r') as f:#open sourcefile
		lineNum = 0
		for line in f:
			lineNum += 1
			dat = line.strip().split('\t')
			if (lineNum == 1):
				sampleIDs = list(dat[9:])#sampleIDs of GT
				for i in range(len(sampleIDs)):
					totalList.append([]) #define list to save GT
			if dat[4] != feature:
				continue
			snp=dat[0] + '_' + dat[1] + '_' + dat[2] + '_' + dat[3]
			label0='V_'
			if snp in sites:
				label0='P_'
			vus1=label0+snp
				
#			homoflag=1
			homoN=0
			genotype = list(dat[9:])#get GT per line
			for i in range(len(genotype)):
				if genotype[i] == '1/1' or genotype[i] == './1':
					homoN = homoN + 1
#					homoflag = 0
			# if homoN > 0:
			# 	with open(output,'a')as fout:
			# 		fout.write('\t'.join(dat[:4])+'\t'+str(feature)+'\t'+str('1')+'\n')
					#break
			vus1 = vus1 +'_'+ str(homoN)
			#print(sampleTmp)
			#nodeList.append(vus1)
			#if homoflag:
			nodeList.append(vus1)
			#print(len(genotype), totalList)
			for i in range(len(genotype)):# save GT to tatalList
				totalList[i].append(genotype[i]) 
	
	#nodeList：
	#P_chr_pos_ref_alt_1
	#v_chr_pos_ref_alt_0
	#v_chr_pos_ref_alt_0.5
	return nodeList, sampleIDs,totalList

def buildG(feature):
	nodeList, sampleIDs, totalList= readTable(feature,sites,inputdir)
	#print('nodeList',nodeList)
	nodeNum = len(nodeList)
	#print('nodeNum',nodeNum)
	matrix = np.zeros((nodeNum, nodeNum))
	sampleDict = {}

	for i in range(len(totalList)):
		sampleID = sampleIDs[i]
		colList = totalList[i]
		
		for j in range(nodeNum):
			nodeStr = nodeList[j]
			if colList[j] == '0/1' or colList[j] == '1/0' or colList[j]=='1/1' or colList[j] == './1':
				#print(colList[j])
				if sampleID not in sampleDict.keys():
					nodeNumDict = {}
					nodeNumDict[nodeStr] = 1
					sampleDict[sampleID] = nodeNumDict
				else:
					nodeNumDict = sampleDict[sampleID]
					if nodeStr not in nodeNumDict.keys():
						nodeNumDict[nodeStr] = 1
					else:
						nodeNumDict[nodeStr] += 1
					sampleDict[sampleID] = nodeNumDict
			else:
				continue
	for key, nodeNumDict in sampleDict.items():
		#print(sampleDict)
		nowNodes = list(nodeNumDict.keys())
		#print(len(nowNodes))
		for i in range(len(nowNodes)):
			node1 = nowNodes[i]
			index1 = nodeList.index(node1)#查找node1的索引位置
			#print(len(nowNodes))
			n1 = nodeNumDict[node1]
			#print(n1)
			for j in range(i + 1, len(nowNodes)):
				node2 = nowNodes[j]
				index2 = nodeList.index(node2)
				n2 = nodeNumDict[node2]
				#print(n2)
				matrix[index1, index2] += min(n1, n2)
				matrix[index2, index1] += min(n1, n2)
				#print('matrix',matrix[index1, index2] )
	G = nx.Graph()
	for i in range(nodeNum):
		for j in range(i,nodeNum):
			if matrix[i][j] != 0:
				#print(i,j,matrix[i][j])
				G.add_edge(nodeList[i], nodeList[j], weight=matrix[i][j])
	
	wDegs = []
	Degs = []
	for node in G.nodes():
		wDegs.append(G.degree(node,weight='weight'))
		Degs.append(G.degree())

	return G,feature,nodeList
#cal weighted nb num ,nowSyb is P/B/V
#output:savefig in (graphFile+geneType+'WNbNum-'+nowSyb+'.png')
def findWgtNbNum(feature):
	G,feature,nodeList = buildG(feature)
	#homosites = []
	for i in nodeList:
		if i not in G.nodes():
			node=i.split('_')[1:5]#get node
			label1 = int(i.split('_')[5])
			if label1 > 0:
				node='\t'.join(node)
				#homosites.append(node)
				with open(output,'a')as f:
					f.write(str(node)+'\t'+str(feature)+'\t'+str('1')+'\n')
	
	for node in G.nodes():
		BA2='NA'
		label2 = node.split('_')[0]#get node label(p,v,b)
		label1 = int(node.split('_')[5])
		snplist=node.split('_')[1:5]#get node
		snp='\t'.join(snplist)
		nbs = G.adj[node]#访问node的邻居节点和边
		pNum = 0
#		Degs=(G.degree(node))#get degree
		for nb, wgt in nbs.items():#统计node的邻居节点nb和边wgt
			nbLabel = nb.split('_')[0]
			if nbLabel == 'P':
				pNum = 1
				break
		if(pNum>0 and label2 != "P"):
			BA2="1"
		if label1 > 0 and label2 != "P":
			BA2="1"
		if BA2 == "1":
			with open(output,'a')as f:
				f.write(str(snp)+'\t'+str(feature)+'\t'+str(BA2)+'\n')

def bk(m):
	m=m+1
	
if __name__ == '__main__':	
	
	inputdir=sys.argv[1]#network.tsv
	inputz=sys.argv[2]#p.tsv
	output=sys.argv[3]#result
	feature=sys.argv[4]#feature
	#if(feature=="NM"):
	# cofflist=[]
	# cofflist.append('chrom\tpos\tref\talt\ttranscript\tPM3\n')
	sites = readp(inputz)
	findWgtNbNum(feature)
	#print ("Sub-process(es) done.")

		
