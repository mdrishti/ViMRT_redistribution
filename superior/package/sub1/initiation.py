import pandas as pd
import os
import requests
import json
import time
import re
import copy
from bs4 import BeautifulSoup
import subprocess

class Initiate:
    def __init__(self, root):
        self.pmcidlist = []
        self.doi = []
        self.pmid_noxml = []
        self.noInfo = []
        self.root = root
        self.outPM = root + "/pmid.txt"
        self.outPMC = root + "/pmcid.csv"
        self.pmid_all = root+ "/pmid.txt"
        self.xmlFail = root + "/xmlfail.csv"
        self.xmlOutDir = root + "/xml"
        self.pmid_filter = root + "/pmid-filter.txt"
        self.outDOI = root + "/doi.txt"
        self.tiabOutDir = root + "/tiab"
        self.pdfOutDir = root + "/pdf"
        self.nofile = root + "/nofile.txt"
        self.txtDir = root + "/txt"
        self.pubDir = root + "/pubtator"
    def getPMIDlist(self):
        infile = self.root + "/"+ self.root+"mutation.xlsx"
        if os.path.exists(infile):
            df = pd.read_excel(infile)
            self.pmidlist = df.pmid.unique().astype('str').tolist()
        else:
            with open(self.pmid_all, 'r') as f:
                self.pmidlist = f.read().split('\n')
    def getPMCID(self):
        base_url = "https://sci-hub.se/"
        if os.path.exists(self.nofile): return
        for i in range(len(self.pmidlist)):
            pmid = self.pmidlist[i]
            if pmid == "": continue
            dctr = json.loads(self.pm_pmc(pmid))
            record = dctr["records"][0]
            if "pmcid" in record.keys():
                pmcid = record["pmcid"]
                self.pmcidlist.append(pmid+","+pmcid)
            if "doi" in record.keys():
                doi = record["doi"]
                self.doi.append(base_url+doi+"\t"+pmid)
            if "pmcid" not in record.keys() and "doi" not in record.keys():
                self.noInfo.append(pmid)
        self.writeIDlist()
    def pm_pmc(self, pmid):
        """transfer PMID to PMCID and DOI"""
        url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=tanfanglin@tongji.edu.cn&ids=%s&format=json&versions=no" % pmid
        return requests.get(url).text
    def writeIDlist(self):
        with open(self.outPM,'w') as f: f.write("\n".join(self.pmidlist))
        with open(self.outPMC,'w') as f: f.write("PMID,PMCID\n"+"\n".join(self.pmcidlist))
        with open(self.nofile,'w') as f: f.write("\n".join(self.noInfo))
    def xmlFile(self):
        if not os.path.exists(self.xmlOutDir): os.mkdir(self.xmlOutDir)
        if os.path.exists(self.xmlFail): return
        dfpmc = pd.read_csv(self.outPMC, header=0, sep = ',')
        pmidlist = copy.deepcopy(self.pmidlist)
        with open(self.xmlFail,'w') as F:
            F.write("PMID,PMCID\n")
            for i in range(len(dfpmc)):
                pmcid = dfpmc['PMCID'][i]
                pmid = str(dfpmc['PMID'][i])
                outfile = self.xmlOutDir + '/' + pmid + ".xml"
                if not os.path.exists(outfile):
                    try:
                        xml = self.getPMC(pmcid)
                        if len(xml) < 200:
                            F.write(pmid+","+pmcid+"\n")
                            continue
                        else: 
                            with open(outfile,'w',encoding='utf-8') as f: f.write(xml)
                            pmidlist.remove(pmid)
                    except Exception as e:
                        print (pmcid, e, end=" ")
                else: 
                    pass
                    # try: pmidlist.remove(pmid)     
                    # except: continue   
        self.pmid_noxml = pmidlist
        with open(self.pmid_filter,'w') as f: f.write("\n".join(pmidlist))
    def getPMC(self, pmcid):
        url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmcids=%s&concepts=mutation" % pmcid
        return requests.get(url, timeout = 5).text
    def getBioC(self, pmid):
        url = 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pubmed.cgi/BioC_xml/%s/ascii' % pmid
        URL = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/%s/ascii" % pmid
        return requests.get(URL).text
    def pubtator(self, pmid):
        format = "pubtator" 
        concepts = "&concepts=mutation"
        url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=%s&concepts=mutation" % pmid
        return requests.get(url, timeout = 5).text
    def tiabFile(self):
        if not os.path.exists(self.tiabOutDir): os.mkdir(self.tiabOutDir)
        for i in range(len(self.pmidlist)):
            pmid = self.pmidlist[i]
            outfile = self.tiabOutDir + "/" +pmid+".tiab.txt"
            if os.path.exists(outfile): continue
            try:
                tiab = self.pubtator(pmid)
                with open(outfile,'w') as f:
                    f.write(tiab)
            except Exception as e:
                print (pmid, e, end=" ")
    def downPDF(self):
        with open(self.outDOI,'w') as f: f.write("\n".join(self.doi))
        outdir = self.pdfOutDir
        if not os.path.exists(outdir): os.mkdir(outdir)
        if len(self.pmid_noxml) == 0: 
            with open(self.pmid_filter) as f:
                self.pmid_noxml = f.read().split("\n")
        if len(self.doi) == 0:
            with open(self.outDOI) as f:
                self.doi = f.read().split("\n")
        for i in range(len(self.doi)):
            line1 = self.doi[i].split("\t")
            if len(line1) < 2: continue
            line =line1[0]
            pmid = line1[1]
            out_fname =outdir+"/"+pmid+'.pdf'
            # if pmid not in self.pmid_noxml: 
            #     if os.path.exists(out_fname): os.remove(out_fname)
            #     continue
            if os.path.exists(out_fname): continue
            print (pmid+" downloading")
            try:
                res = requests.get(line)
                res.encoding = 'utf-8'
                soup = BeautifulSoup(res.text, 'html.parser')
                news = soup.select('iframe')
            except Exception as e:
                print (e)
                time.sleep(5)
                continue
            if len(news) == 0: 
                print (pmid +" no news")
                continue
            src=news[0]['src']
            if re.match(r"//sci-hub|//dacemir", src[0:9]):
                pdf = "https:"+src
                self.download(pdf, out_fname)
                continue
            elif re.match(r"https://mosc|https://zero|https://sci-|https://twin|https://cybe|https://dace|https://zero", src[0:12]) or re.match(r"http://www.hxkq", src[0:15]):
                pdf = src
                self.download(pdf, out_fname)
                continue
            print (pmid+" failed")
    def mkTXT(self):
        if not os.path.exists(self.txtDir): os.mkdir(self.txtDir)
        if not os.path.exists(self.pubDir): os.mkdir(self.pubDir)
    def download(self, pdf, out_fname):
        try:
            r = requests.get(pdf)
            with open(out_fname,'wb') as f2:
                f2.write(r.content)
        except Exception as e:
            print (e)
            time.sleep(5)
            r = requests.get(pdf)
            with open(out_fname,'wb') as f2:
                f2.write(r.content)

                
class Initiate1(Initiate):
    def getPMIDlist(self):
        for r,d,f in os.walk(self.root):
            for fi in f:
                if "add" in fi:
                    infile = r+'/'+fi
                    break
        with open(infile) as f:
            pmid = f.read().strip('\n').split('\n')
        self.pmidlist = pmid


class Initiate2(Initiate1):
    bash_line = """
    esearch -db pubmed -query "{}" |
    efetch -format xml |
    xtract -head '<html><body>' -tail '</body></html>' \
        -pattern PubmedArticle -PMID MedlineCitation/PMID \
        -block ArticleId -if @IdType -equals doi \
        -tab '\n' -pfx '<p><a href="http://dx.doi.org/' \
        -sep '">' -sfx '</a></p>' -encode ArticleId,"&PMID"
            """
    def getDOI(self):
        def new(pmidlist):
            pmid = '\n'.join(pmidlist).replace("\n", "[pmid] or ") + "[pmid]"
            a, b = subprocess.getstatusoutput(Initiate2.bash_line.format(pmid))
            doi = b.split("\n")
            newdoi = []
            for i in range(1,len(doi)-1):
                d = doi[i]
                newd = d.replace("<p><a href=\"","").replace("\">","\t").replace("</a></p>","").replace("http://dx.doi.org","https://sci-hub.se")
                newdoi.append(newd)
                print (newd)
            return newdoi
        N = int(len(self.pmidlist) / 1000)
        if N > 0:
            for i in range(N):
                self.doi += new(self.pmidlist[1000*i : 1000*(i+1)])
        self.doi += new(self.pmidlist[1000*N:])
            



