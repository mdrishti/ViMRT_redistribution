import os
import sys, getopt
from package import Batch
os.chdir("./")
def main(argv):
    inputfile = ''
    outputpath= '.'
    source= 'PubMed_PMC'
    try:
        opts, args = getopt.getopt(argv,"h:i:o:s:",["help","ifile=","opath=","source="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            sys.exit()
        if opt in ("-o", "--ofile"):
            outputpath = arg
        if opt in ("-i", "--ipath"):
            inputfile = arg      
        if opt in ("-s", "--source"):
            source = arg
    return outputpath,inputfile,source
          
if __name__ == "__main__":
    outputpath,inputfile,source=main(sys.argv[1:])
    if source=="PubMed":  
        Batch.BiocPubmed(outputpath,inputfile).getBatch()
    if source=="PMC":  
        Batch.BiocPmcoa(outputpath,inputfile).getBatch()
    if source=="PubTator":  
        Batch.PubTatorTiab(outputpath,inputfile).getBatch()
    if source=="PubMed_PMC":  
        Batch.BiocPubmed(outputpath,inputfile).getBatch()
        Batch.BiocPmcoa(outputpath,inputfile).getBatch()
