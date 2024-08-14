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
def readp():
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
		
def readTable(feature):
	"""compute nodeList, sampleIDs, totalList one gene
		input sourcefile:gene_source(include head)
		return snp sampleIDs genoTypes 
	"""
	sites=readp()
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
			tmp = line.replace('\n', '').split('\t')
			#print(tmp)
			if (lineNum == 1):
				sampleIDs = list(tmp[9:])#sampleIDs of GT
				# print(sampleIDs)
				for i in range(len(sampleIDs)):
					totalList.append([]) #define list to save GT
			if tmp[4] != feature:
				continue
			
			snp=tmp[0] + '_' + tmp[1] + '_' + tmp[2] + '_' + tmp[3]
			#print(snp)
			label0='V_'
			for i in sites:
				if (i==snp):
					label0='P_'
				else:
					continue
			vus1=label0+ snp 
		#print(tmp[0])
			sampleTmp = list(tmp[9:])#get GT per line
			mod=0
			for i in sampleTmp:
				if i == '1/1':
					#print(i)
					mod+=1
			#print(mod)
			
			if(mod==1):
				label1='_0.5'
			elif(mod>1):
				label1='_1'
			else:
				label1='_0'
			#print(vus1+label1)
			nodeList.append(vus1+label1)
		#print(sampleTmp)
			for i in range(len(sampleTmp)):# save GT to tatalList
				totalList[i].append(sampleTmp[i]) 
	
	#nodeList：
	#P_chr_pos_ref_alt_1
	#v_chr_pos_ref_alt_0
	#v_chr_pos_ref_alt_0.5
	return nodeList, sampleIDs,totalList

def buildG(feature):
	nodeList, sampleIDs, totalList= readTable(feature)
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
			if colList[j] == '0/1' or colList[j] == '1/0' or colList[j]=='1/1':
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
	#G,feature = buildG(feature)
	G,feature,nodeList = buildG(feature)
	homozygosity=[]
	for i in nodeList:
		source0= float(i.split('_')[-1])
		if (source0>0):
			homozygosity.append(i)
	for i in homozygosity:
		if i not in G.nodes():
			#print(i)
			node=i.split('_')[1:5]#get node
			label1 = i.split('_')[0]
			sourcelabel=i.split('_')[-1]
			node='\t'.join(node)
			with open(output,'a')as f:
				if (sourcelabel=='0.5'):
					f.write(str(node)+'\t'+str(feature)+'\t'+str('3')+'\n')
				if (sourcelabel=='1'):
					f.write(str(node)+'\t'+str(feature)+'\t'+str('2')+'\n')
	print(feature)
	pNums = []
	bNums = []
	vNums = []
	#print(nowSyb)
	for node in G.nodes():
		label2 = node.split('_')[0]#get node label(p,v,b)
		source= float(node.split('_')[-1])
		#print(label2)
		snplist=node.split('_')[1:5]#get node
		snp='\t'.join(snplist)
		nbs = G.adj[node]#访问node的邻居节点和边
		pNum = 0
		#wDegs=(G.degree(node,weight='weight'))#get weight
		Degs=(G.degree(node))#get degree
		#print('wDegs:',wDegs)
	
		for nb, wgt in nbs.items():#统计node的邻居节点nb和边wgt
			nbLabel = nb.split('_')[0]
			#snplink=nb.split('_')[1:5]
			if nbLabel == 'P':
				#pNum += wgt['weight']
				pNum += 1
		pNum+=source
			#if nbLabel == 'B' :
				#bNum += wgt['weight']
			#elif nbLabel == 'V':
				#vNum += wgt['weight']
		if(pNum==4 or pNum>4):
			PM3="0"
		elif(pNum==2 or pNum==3 or pNum==2.5 or pNum==3.5):
			PM3="1"
		elif(pNum==1 or pNum==1.5):
			PM3="2"
		elif(pNum==0.5):
			PM3="3"
		else:
			PM3="NA"
		#print('PM3:',PM3)
		#print('PM3:',str(snp)+'\t'+str(feature)+'\t'+str(PM3)+'\t'+str(snplink)+'\n')
		with open(output,'a')as f:
			f.write(str(snp)+'\t'+str(feature)+'\t'+str(PM3)+'\n')
	#print('%.2f' % basStats(pNums),'%.2f' % basStats(bNums),'%.2f' % basStats(vNums))

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
	findWgtNbNum(feature)
	#print ("Sub-process(es) done.")

		
