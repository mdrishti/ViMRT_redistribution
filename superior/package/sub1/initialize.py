import pandas as pd
import os
import requests
import json
import time
import re
import copy
from bs4 import BeautifulSoup
import subprocess
VirusNames = ['EBV','HBV','HIV','HPV','HTLV1','XMRV']
BaseUrl = "https://www.ncbi.nlm.nih.gov"
class Begin:
    @staticmethod
    def down(info):
        crawl = CrawlerMonitor(info)
        download_ist = crawl.identify_download()
        download_ist.download()
    def __init__(self) -> None:
        self.pmidDict = {}
    def _getPMID(self):
        for virusName in VirusNames:
            filename = 'GoldStandard/%slast1.xlsx' % virusName
            df = pd.read_excel(filename)
            pmidlist = df.pmid.unique().astype('str').tolist()
            self.pmidDict[virusName] = pmidlist
    def getPubTatorAbstract(self):
        for virusName, pmidlist in self.pmidDict.items():
            info = {
                'directory':{'outdir':'download/'+virusName},
                'IDlist':{'pmid':pmidlist},
                'source':{'api':'PubTator','content':'Abstract'}
                }
            Begin.down(info)
    def getBionlpAbstract(self):
        for virusName, pmidlist in self.pmidDict.items():
            info = {
                'directory':{'outdir':'download/'+virusName},
                'IDlist':{'pmid':pmidlist},
                'source':{'api':'Bionlp','content':'Abstract'}
                }
            Begin.down(info)
    def _getPMCID(self):
        self.pmcidDict = {}
        for virusName, pmidlist in self.pmidDict.items():
            pmcidlist = ID.getBatchPMCID(pmidlist)
            self.pmcidDict[virusName] = pmcidlist
    def getPubTatorPMCtext(self):
        for virusName, pmcidlist in self.pmcidDict.items():
            info = {
                'directory':{'outdir':'download/'+virusName},
                'IDlist':{'pmid':pmcidlist},
                'source':{'api':'PubTator','content':'Full-text'}
            }
            Begin.down(info)
    def getBionlpPMCtext(self):
        for virusName, pmcidlist in self.pmcidDict.items():
            info = {
                'directory':{'outdir':'download/'+virusName},
                'IDlist':{'pmid':pmcidlist},
                'source':{'api':'Bionlp','content':'Full-text'}
            }
            Begin.down(info)
    def _getDOI(self):
        self.doiDict = {}
        for virusName, pmidlist in self.pmidDict.items():
            doilist = ID.getBatchDOI(pmidlist)
            self.doiDict[virusName] = doilist
    def getEFetchMedline(self):
        for virusName, pmidlist in self.pmidDict.items():
            info = {
                'directory':{'outdir':'download/'+virusName},
                'IDlist':{'pmid':pmidlist},
                'source':{'api':'EFetch','content':'Abstract'}
            }
            Begin.down(info)
    
    
            
class ID:
    bash_line_DOI = """
    esearch -db pubmed -query "{}" |
    efetch -format xml |
    xtract -head '<html><body>' -tail '</body></html>' \
        -pattern PubmedArticle -PMID MedlineCitation/PMID \
        -block ArticleId -if @IdType -equals doi \
        -tab '\n' -pfx '<p><a href="http://dx.doi.org/' \
        -sep '">' -sfx '</a></p>' -encode ArticleId,"&PMID"
            """
    def __init__(self) -> None:
        pass
    @staticmethod
    def getDOI(pmidlist):
        pmid = '\n'.join(pmidlist).replace("\n", "[pmid] or ") + "[pmid]"
        a, b = subprocess.getstatusoutput(ID.bash_line_DOI.format(pmid))
        doi = b.split("\n")
        newdoi = []
        for i in range(1, len(doi)-1):
            d = doi[i]
            newd = d.replace("<p><a href=\"","").replace("\">","\t").replace("</a></p>","").replace("http://dx.doi.org","https://sci-hub.se")
            newdoi.append(newd)
            print (newd)
        return newdoi
    @staticmethod
    def getBatchDOI(pmidlist):
        doiList = []
        N = int(len(pmidlist) / 1000)
        if N > 0:
            for i in range(N):
                doiList += ID.getDOI(pmidlist[1000*i : 1000*(i+1)])
        doiList += ID.getDOI(pmidlist[1000*N:])
        return doiList
    @staticmethod
    def getPMCID(pmid):
        """transfer PMID to PMCID and DOI"""
        url = BaseUrl + "/pmc/utils/idconv/v1.0/?tool=my_tool&email=tanfanglin@tongji.edu.cn&ids=%s&format=json&versions=no" % pmid
        return requests.get(url).text
    @staticmethod
    def getBatchPMCID(pmidlist):
        pmcidList = []
        for pmid in pmidlist:
            if pmid == "": continue
            dctr = json.loads(ID.getPMCID(pmid))
            record = dctr["records"][0]
            if "pmcid" in record.keys():
                pmcid = record["pmcid"]
                pmcidList.append(pmid+","+pmcid)
        return pmcidList

class Download:
    def __init__(self, Download_data) -> None:
        self.Download_data = Download_data
        self.outdir = Download_data['directory']['outdir']
        if not os.path.exists(self.outdir): os.mkdir(self.outdir)
        self.pmidlist = Download_data['IDlist']['pmid']
        self.dirname = Download_data['source']['api']
        self.suffix = Download_data['source']['content']
        self.outfile = "{0}/{1}/%s.{2}".format(self.outdir, self.dirname, self.suffix)
    @staticmethod
    def meets_condition(Download_data: dict):
        return False
    @staticmethod
    def meets_condition_pre(Download_data: dict):
        """Precondition of the contract of this interface.
        Validate that the ``Download_data`` parameter is properly formed.
        """
        assert isinstance(Download_data, dict), f"{Download_data!r} is not a dict"
        for info in ("directory", "IDlist", "source"):
            assert info in Download_data, f"{info} not in {Download_data}"
            assert isinstance(Download_data[info], dict)
    def download(self): pass

class UnknownDownload(Download):
    """A type of Download that cannot be identified from its data"""

class PubTatorAbstractDownload(Download):
    @staticmethod
    def meets_condition(Download_data: dict):
        return (
            Download_data["source"]["api"] == "PubTator"
            and Download_data["source"]["content"] == "Abstract"
        )
    def _pubtator(self, pmid):
        url = BaseUrl + "/research/pubtator-api/publications/export/pubtator?pmids=%s&concepts=mutation" % pmid
        return requests.get(url, timeout = 5).text
    def download(self):
        for pmid in self.pmidlist:
            outfile = self.outfile % pmid
            if os.path.exists(outfile): continue
            try:
                text = self._pubtator(pmid)
                with open(outfile,'w') as f:
                    f.write(text)
            except Exception as e:
                print (pmid, e, end=" ")

class PubTatorFullTextDownload(Download):
    @staticmethod
    def meets_condition(Download_data: dict):
        return (
            Download_data["source"]["api"] == "PubTator"
            and Download_data["source"]["content"] == "Full-text"
        )
    def _pubtator(self, pmcid):
        url = BaseUrl + "/research/pubtator-api/publications/export/biocxml?pmcids=%s&concepts=mutation" % pmcid
        return requests.get(url, timeout = 5).text
    def download(self):
        for pm_pmc in self.pmidlist:
            pmid, pmcid = pm_pmc.split(',')
            outfile = self.outfile % pmid
            if os.path.exists(outfile): continue
            try:
                text = self._pubtator(pmcid)
                with open(outfile,'w') as f:
                    f.write(text)
            except Exception as e:
                print (pmid, e, end=" ")
    

class EFetchAbstractDownload(Download):
    bash_line = """
    esearch -db pubmed -query "{}" | 
    efetch -format "medline"
    """
    @staticmethod
    def meets_condition(Download_data: dict):
        return (
            Download_data["source"]["api"] == "EFetch"
            and Download_data["source"]["content"] == "Abstract"
        )
    def _efetch(self):
        pmids_line = "[pmid] or ".join(self.pmidlist) + "[pmid]"
        a, b = subprocess.getstatusoutput(self.bash_line.format(pmids_line))
        return b
    def download(self):
        all_medlines = self._efetch()
        with open(self.outfile % 'test','w') as f:
            f.write(all_medlines)

class BionlpFullTextDownload(Download):
    @staticmethod
    def meets_condition(Download_data: dict):
        return (
            Download_data["source"]["api"] == "Bionlp"
            and Download_data["source"]["content"] == "Full-text"
        )

class BionlpFullTextDownload(Download):
    @staticmethod
    def meets_condition(Download_data: dict):
        return (
            Download_data["source"]["api"] == "Bionlp"
            and Download_data["source"]["content"] == "Full-text"
        )

class PDFFullTextDownload(Download):
    @staticmethod
    def meets_condition(Download_data: dict):
        return Download_data["source"]["api"] == "scihub"
    
class CrawlerMonitor:
    def __init__(self, download_data):
        self.download_data = download_data

    def identify_download(self):
        Download.meets_condition_pre(self.download_data)
        download_cls = next(
            (
                download_cls
                for download_cls in Download.__subclasses__()
                if download_cls.meets_condition(self.download_data)
            ),
            UnknownDownload,
        )
        return download_cls(self.download_data)
