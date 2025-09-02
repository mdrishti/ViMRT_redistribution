import os
from superior.package.sub1 import initiation
def initiate(root):
    run = initiation.Initiate2(root)
    run.getPMIDlist()
    print ('tiab')
    run.tiabFile()
    print ('PMCID')
    run.getPMCID()
    print ('xml')
    run.xmlFile()
    print ('DOI')
    run.getDOI()
    print ('down')
    run.downPDF()
    run.mkTXT()
    
if __name__ == "__main__":
    root = "VIRUS/version_1"
    for r,d,f in os.walk(root):
        rootlist = d
        break
    for ro in rootlist:
        print (ro)
        r = root+'/'+ro
        initiate(r)
    