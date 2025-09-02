import os
import sys, getopt
from package import Batch
from package import found
os.chdir("./")
def main(argv):
    inputpath = '.'
    outputpath= '.'
    form= 'BioC'
    virus= 'Unknown'
    method= "regex+rules"
    concat= "no"
    try:
        opts, args = getopt.getopt(argv,"h:i:o:v:f:m:c:",["help","ipath=","opath=","virus=","form=","method=","concat="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            sys.exit()
        if opt in ("-o", "--opath"):
            outputpath = arg
        if opt in ("-i", "--ipath"):
            inputpath = arg      
        if opt in ("-v", "--virus"):
            virus = arg
        if opt in ("-f", "--form"):
            form= arg
        if opt in ("-m", "--method"):
            method = arg
        if opt in ("-c", "--concat"):
            concat = arg
    return inputpath,outputpath,virus,form,method,concat
          
if __name__ == "__main__":
    inputpath,outputpath,virus,form,method,concat=main(sys.argv[1:])        
    if concat=="concat":
        result = found.Result(inputpath,outputpath,form,concat)
        result._rulesAppendRegex()
    else:
        if method == "rules":
            f = found.Found(inputpath,outputpath,virus,form,"rules")
            f.save()
        elif method == "regex": 
            f = found.Found(inputpath,outputpath,virus,form,"regex")
            f.save()
        else: 
            f = found.Found(inputpath,outputpath,virus,form,"regex+rules")
            inputpath=outputpath
            result = found.Result(inputpath,outputpath,form)
            f.save()
            result._rulesAppendRegex()
        
    

    


