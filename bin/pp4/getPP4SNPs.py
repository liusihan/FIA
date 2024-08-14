import os 
import shutil

def detectExistDirectory(tempDir):
	if os.path.exists(tempDir):
		shutil.rmtree(tempDir)
	os.makedirs(tempDir)


def buildPP4Files(filepath, fileNames):
    for i in range(len(fileNames)):
        fileName = filepath + fileNames[i] + '.tsv'
        with open(fileName, 'w') as fobjw:
            fobjw.write("chrom\tpos\tref\talt\ttranscript\tsymbol\tgene\tgnomAD_AF_SAS_WES\tgnomAD_grpmax_WES\tgnomAD_AF_grpmax_WGS\tgnomAD_AF_EAS_WES\tgnomAD_AF_ASJ_WES\tgnomAD_AF_MID_WES\tgnomAD_AF_FIN_WES\tgnomAD_AF_NFE_WES\tGT\tAD\tDP\tGQ\tPL\tPP4\n")
        
def getPP4Files(SamplePheno_dict, GenePheno_dict, inputFile, HLsampleNames, outputFilePath):
    inputFile = inputFile
    HLsampleNames = HLsampleNames
    buildPP4Files(outputFilePath, HLsampleNames)
    with open(inputFile,'r',encoding='UTF-8') as fobjr:
        for eachline in fobjr:
            eachline = eachline.strip().split("\t")
#             if eachline[0] == 'Sample':
# #                fobjw.write("\t".join(eachline)+"\tPP4\n")
#                 eachline= fobjr.readline()
#                 eachline = eachline.strip().split("\t")
#             if eachline[4] == "*":
#                 pass
            if eachline[0] == 'Sample' or eachline[4] == "*":
                continue
            else:
                sample = eachline[0]
                sampleFile = outputFilePath + sample + '.tsv'
                with open(sampleFile, 'a') as fobja:
                    transList = eachline[5].split('|')
                    symbolList = eachline[6].split('|')
                    geneList = eachline[7].split('|')

                    transNum = len(transList)
                    for i in range(transNum):
                        symbol = symbolList[i]
                        if (sample in SamplePheno_dict) and (SamplePheno_dict[sample] in GenePheno_dict) and (GenePheno_dict[SamplePheno_dict[sample]]==symbol):
                            fobja.write("\t".join(eachline[1:5])+'\t'+transList[i]+'\t'+symbol+'\t'+geneList[i]+'\t'+"\t".join(eachline[8:len(eachline)])+'\t'+"3\n")
                        else:
                            fobja.write("\t".join(eachline[1:5])+'\t'+transList[i]+'\t'+symbol+'\t'+geneList[i]+'\t'+"\t".join(eachline[8:len(eachline)])+'\t'+"NA\n")

                        # if (symbol in sgeneSamDict) and (sample in sgeneSamDict[symbol].strip(",").split(",")):
                        #     fobja.write("\t".join(eachline[1:5])+'\t'+transList[i]+'\t'+symbol+'\t'+geneList[i]+'\t'+eachline[8]+'\t'+eachline[9]+'\t'+"1"+"\t".join(eachline[9:len(eachline)])+"\n")
                        # elif (symbol in mgeneSamDict) and (sample in mgeneSamDict[symbol].strip(",").split(",")):
                        #     fobja.write("\t".join(eachline[1:5])+'\t'+transList[i]+'\t'+symbol+'\t'+geneList[i]+'\t'+eachline[8]+'\t'+eachline[9]+'\t'+"2"+"\t".join(eachline[9:len(eachline)])+"\n")
                        # else:
                        #     fobja.write("\t".join(eachline[1:5])+'\t'+transList[i]+'\t'+symbol+'\t'+geneList[i]+'\t'+eachline[8]+'\t'+eachline[9]+'\t'+"NA"+"\t".join(eachline[9:len(eachline)])+"\n")
 
