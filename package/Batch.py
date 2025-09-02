import json
import requests
import subprocess
import os,re
import xml.etree.ElementTree as ET
class Batch:

    steps = [5000, 2000, 1000, 200, 100, 10, 1]

    def __init__(self, root:str, inputfile:str):
        self.outdir = root
        self.root = root
        if not os.path.exists(self.outdir):os.makedirs(self.outdir)
    def _getFunc(self):
        pass

    def _getlines(self):
        pass

    def getBatch(self):
        for batch_size in self.steps:
            print('batch_size: {}'.format(batch_size))
            self._batchRequest(self._getlines(), batch_size)

    def _batchRequest(self, lines:list, batch_size:int):
        lines = self._checkExist(lines)
        if lines == [] or lines == ['']:
            return
        batch_num = len(lines)//batch_size
        func = self._getFunc()
        for i in range(batch_num):
            print(i, end=' \n')
            try:
                func(lines[i*batch_size:(i+1)*batch_size])
            except Exception as e:
                print ('error:',e, ' failed', end='!')
        in_batch_num = batch_size*batch_num
        if len(lines) > in_batch_num:
            func(lines[in_batch_num:])

    def _checkExist(self, lines:list):
        if self.outdir == '':
            return lines
        print (lines)
        new_lines = [line for line in lines if not os.path.exists(self.outdir.format(pmid=str(line)))]
        # Rule out the files that are already downloaded 

        return new_lines

    def _writeEachFile(self, nameList:list, contentList:list):
        for i in range(len(nameList)):
            outfile = self.outdir.format(pmid=nameList[i])
            try:
                with open(outfile, 'w') as f:
                    f.write(contentList[i])
            except:
                return
class BioC(Batch):

    """https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pubmed.cgi/BioC_xml/20869540/unicode"""
    """https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/31601912/ascii"""
    base_url = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/{source}.cgi/BioC_xml/{pmid}/{encode}"

    def __init__(self, root:str, inputfile:str):
        super().__init__(root,inputfile)
        self.steps = [1]
        self.pmidlist = inputfile
        self.bioc = self.root+'/bioc'
        if not os.path.exists(self.bioc):
            os.mkdir(self.bioc)

    def _getFunc(self):

        return self._biocFile

    def _getlines(self):
        with open(self.pmidlist) as f:
            pmidlist = f.read().strip().split('\n')

        return pmidlist

    def _request(self, id:str):

        return requests.get(self.url.format(pmid=id), timeout = 5).text

    def _biocFile(self, pmidlist:list):
        assert len(pmidlist) == 1, "{pmidlist} must have one pmid, and only one".format(pmidlist=pmidlist)
        pmid = pmidlist[0]
        xml = self._request(pmid)
        if xml[:7] == "[Error]":
            with open(self.bioc+"/nopmcxml.txt", 'a') as Errorf:
                Errorf.write(pmid+"\n")
        else:
            xml = xml.replace('\n','').replace('  ','')
            with open(self.outdir.format(pmid=pmid), 'w') as f:
                f.write(xml)

class BiocPubmed(BioC):
   
    url = BioC.base_url.format(source="pubmed", pmid="{pmid}", encode="ascii")

    def __init__(self, root:str, inputfile:str):
        super().__init__(root, inputfile)
        print('BioC Pubmed Open Access')
        self.pubmed = self.bioc + '/pubmed'
        if not os.path.exists(self.pubmed):
            os.mkdir(self.pubmed)
        self.outdir = self.pubmed+"/{pmid}.pubmed.xml"
        
class BiocPmcoa(BioC):
    
    url = BioC.base_url.format(source="pmcoa", pmid="{pmid}", encode="ascii")
    print (url)
    def __init__(self,  root:str, inputfile:str):
        super().__init__(root, inputfile)
        print('BioC PMC Open Access')
        self.pmcoa = self.bioc + '/pmcoa'
        if not os.path.exists(self.pmcoa):
            os.mkdir(self.pmcoa)
        self.outdir = self.pmcoa+"/{pmid}.pmcoa.xml"



if __name__ == "__main__":
    #BiocPmcoa('IV','./IV/pmidlist.txt',).getBatch()
    BiocPubmed('Virus','pmidlist.txt',).getBatch()
