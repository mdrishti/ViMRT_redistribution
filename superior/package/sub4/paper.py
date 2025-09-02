import re
from superior.package.sub3.virus import VirusName
from superior.package.sub3 import tmVar
class PaperFile:
    def __init__(self, filename:str, virus:str=None, file_format='pubtator'):
        
        if 'xml' in filename and 'tor' not in filename:
            file_format = 'xml'
        if file_format in ['pubtator', 'tiab', 'PubTator']:
            self.paper = tmVar.TXT(filename)
        elif file_format in ['xml', 'bioc','pubmed', 'XML']:
            self.paper = tmVar.XML(filename)
        else:
            print(file_format)
        self.paper.all()
        self.pmid = self.paper.pmid
        if self.pmid == "": print(filename)
        self.virusinput = virus
        self.virus ='Unknown'
    def virususer(self):
        if (self.virusinput in ['EBV', "SARS-CoV-2",'HBV',"IV", 'HIV', 'HPV', 'HTLV1', 'XMRV', 'MCV']):
            virus_cls = VirusName(self.virusinput).identify_virus()
            print (virus_cls,2)
            return virus_cls({
                'pmid': self.paper.pmid,
                'passage': self.paper.passage,
                'annotation': self.paper.annotation,
                'useinput': self.virusinput
            })

    def processPaper(self):
        #if (self.virus not in ['EBV', "Unknown","SARS-CoV-2",'HBV',"IV", 'HIV', 'HPV', 'HTLV1', 'XMRV', 'MCV']):
        self._getVirus()
        
        virus_cls = VirusName(self.virus).identify_virus()
        return virus_cls({
            'pmid': self.paper.pmid,
            'passage': self.paper.passage,
            'annotation': self.paper.annotation,
            'useinput': self.virusinput
        })
    def _getVirus(self):
        IdentifyVirus._getVirusNamePatterns()
        idVirus = IdentifyVirus(self.paper.passage)
        idVirus._virusDict()
        virusDict = idVirus.Virus
        self.virusDict = virusDict
        self._getVirusMatch()
    def _getVirusMatch(self):
        titleVirus = []
        abstractVirus = []
        fulltextVirus = []
        for key,value in self.virusDict.items():
            if key[0].lower() == "t":
                titleVirus = value
            elif key[0].lower() == "a":
                abstractVirus = value
            else:
                fulltextVirus.extend(value)
        for found in (titleVirus, abstractVirus, fulltextVirus):
            if len(found) > 0:
                self.virus = found[0][0]
                return
        print(self.paper.pmid, self.virusDict, " no virus")
        self.virus = 'UKnown'   
    def _processVirusDict(self):
        flag = False
        vd = self.virusDict
        if not('t' in vd.keys() and 'a' in vd.keys()): return
        aV, tV = vd['a'], vd['t']
        def process(vt, virus):
            if len(vt) == 2 and vt[0][1] == vt[1][1]:
                return 2
            elif len(vt) == 2 and vt[0][0] == virus:
                return False
            elif len(vt) == 2 and vt[0][0] != virus:
                return 2
            elif len(vt) == 1 and vt[0][0] != virus:
                return True
            elif len(vt) == 1 and vt[0][0] == virus:
                return False
            else:
                return 2
        flagT = process(tV, self.virus)
        flagA = process(aV, self.virus)
        if flagT == 2:
            if flagA == 2:
                flag = True
            elif flagA:
                flag = True
            else:
                flag = False
        else:
            flag = flagT
        if flag: 
            pass# (self.pmid, 'flagT=%s'%flagT, 'flagA=%s'%flagA, vd)
        IdentifyVirus.count(flag)
class IdentifyVirus:
    VirusNameDict = {
        'EBV':([r'Epstein-Barr virus',r'human herpesvirus 4'],['EBV','HHV-4']),
        'HBV':([r'Hepatitis B Virus',r'hepatitis B'],['HBV']),
        'HIV':([r'Human Immunodeficiency Virus',r'\bhiv\b'],['HIV']),
        'HPV':([r'human papillomavirus'],['HPV']),
        'HTLV1':([r'Human T-lymphotropic virus type-1',r'Human T-cell leukemia virus type 1',
            r'Human T-lymphotropic virus type 1', r'Human T-cell leukemia virus type-1'],
            ['HTLV1','HTLV-1']),
        'SARS-CoV-2':([r'Severe Acute Respiratory Syndrome Coronavirus 2',r'COVID-19'],['SARS-CoV-2']),
        'IV':([r'Influenza Virus'],['IV']),
        'MCV':([r'Merkel cell virus'],['MCV']),
        'XMRV':([r'xenotropic murine leukemia virus-related virus'],['XMRV']),
        'ALL':([r'Merkel cell virus'],['ALL']),
        }
    VirusNamePatternDict = {}
    right, wrong = 0, 0
    @classmethod
    def _getVirusNamePatterns(cls):
        for virus, nameTuple in cls.VirusNameDict.items():
            fullnameList, acronymList = nameTuple
            fullPattern = re.compile(r'|'.join(fullnameList), re.IGNORECASE)
            acronymPattern = re.compile(r'|'.join(acronymList))
            cls.VirusNamePatternDict[virus] = (fullPattern, acronymPattern)
    @classmethod
    def count(cls, flag:bool):
        if flag:
            cls.wrong += 1
        else:
            cls.right += 1
    def __init__(self, passage:dict):
        self.passage = passage
    @property
    def Virus(self):
        return self._Virus
    def _virusDict(self):
        self._Virus = {}
        for k, v in self.passage.items():
            region, offset = k
            found = self._findVirus(v)
            if len(found) > 0:
                self._Virus[region] = found
                break
    def _findVirus(self, text):
        foundVirusDict = {}
        for virus, patternTuple in self.VirusNamePatternDict.items():
            fullPattern, acronymPattern = patternTuple
            fullFound = fullPattern.findall(text)
            acronymFound = acronymPattern.findall(text)
            found = fullFound + acronymFound
            if len(found) == 0: continue
            foundVirusDict[virus] = len(found)
            break
        foundVirus = sorted(foundVirusDict.items(),  key=lambda d: d[1], reverse=True)
        return foundVirus