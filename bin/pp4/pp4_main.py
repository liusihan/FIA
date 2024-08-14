import getPP4SNPs as gps
import sys
import time



 
if __name__ == '__main__':   
    SamplePheno_dict = {}
    with open(sys.argv[1], 'r') as file:
        for line in file:
            key, value = line.strip().split("\t")
            if key != "ID":
                SamplePheno_dict[key] = value
    GenePheno_dict = {}
    with open(sys.argv[2], 'r') as file:
        for line in file:
            value, key = line.strip().split("\t")
            GenePheno_dict[key] = value
    pp4ResultPath = sys.argv[3]
    tsvFile = sys.argv[4]
    VCFfile = sys.argv[5]
    HLsampleNames = list(SamplePheno_dict.keys())
    begin_time=int(round(time.time()))
    gps.getPP4Files(SamplePheno_dict,GenePheno_dict, tsvFile, HLsampleNames, pp4ResultPath)
    end_time=int(round(time.time()))
    print("--------------------------pp4 tag time： %f min-----------------------" % ((end_time-begin_time)/60))
    
