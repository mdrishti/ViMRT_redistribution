from this import d
import requests
import os
import xml.etree.ElementTree as ET
class VirusPaperAll:
    ROOT = 'Virus'
    def __init__(self, virus):
        self.base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        self.virus_url = "{0}?db=pubmed&&term={1}&&retmax=15000000".format(self.base_url, '{}')
        self.virus = virus
        self.IDfile = self.ROOT+'/'+self.virus+'/listFile.txt'
    def _tree(self):
        xml = self._getXml()
        return ET.fromstring(xml)
    def ID(self):
        if os.path.exists(self.IDfile):
            return
        tree = self._tree()
        IDlist = tree.find('IdList')
        searchList = []
        for Id in IDlist:
            searchList.append(Id.text)
        with open(self.IDfile,'w') as f:
            f.write('\n'.join(searchList))

    def _getXml(self):
        url = self._getUrl()
        print(url)
        return requests.get(url, timeout=20).text

    def _getUrl(self):
        UrlDict = {
            'SARS_CoV_2': self._SARS_CoV_2_Url(),
            'SARS_CoV': self._SARS_CoV_Url(),
            'MERS_CoV': self._MERS_CoV_Url(),
            'IV': self._IV_Url(),
            'HSV_1': self._HSV_1_Url(),
            'HSV_2': self._HSV_2_Url(),
            'VZV': self._VZV_Url(),
            'HCMV': self._HCMV_Url(),
            'KSHV': self._KSHV_Url(),
            'EBV': self._EBV_Url(),
            'HBV': self._HBV_Url(),
            'HIV': self._HIV_Url(),
            'HPV': self._HPV_Url(),
            'HTLV_1': self._HTLV_1_Url(),
            'MCV': self._MCV_Url(),

        }
        return UrlDict[self.virus]  
    def _EBV_Url(self):
        self.term_query = "EBV[tiab]+or+Epstein-Barr+Virus[tiab]+or+HHV-4[tiab]+or+human+herpesvirus+4[tiab]"
        return self.virus_url.format(self.term_query)
    
    def _HBV_Url(self):
        self.term_query = "HBV[tiab]+or+Hepatitis+B+Virus[tiab]"
        return self.virus_url.format(self.term_query)

    def _HIV_Url(self):
        self.term_query = "HIV[tiab]+or+Human+Immunodeficiency+Virus[tiab]"
        return self.virus_url.format(self.term_query)
    
    def _HPV_Url(self):
        self.term_query = "HPV[tiab]+or+Human+Papilloma+Virus[tiab]"
        return self.virus_url.format(self.term_query)
        
    def _HTLV_1_Url(self):
        self.term_query = """
        HTLV1[tiab]+or+HTLV-1[tiab]+or+
        Human+T-lymphotropic+Virus+Type-1[tiab]+or+
        Human+T-lymphotropic+virus+type+1[tiab]+or+
        Human+T-cell+leukemia+virus+type+1[tiab]+or+
        Human+T-cell+leukemia+virus+type-1[tiab]
        """
        return self.virus_url.format(self.term_query)

    def _MCV_Url(self):
        self.term_query = "MCV[tiab]+or+Merkel+cell+virus[tiab]"
    
    def _SARS_CoV_2_Url(self):
        SARS_CoV_2_query_term = "(SARS-CoV-2[tiab]+or+SARS-CoV2[tiab]+or+COVID-19[tiab]+or+covid-19[tiab]+or+severe+acute+respiratory+syndrome+coronavirus+2[tiab]+or+coronavirus+disease[tiab])"
        date_query_term = "(2020/01/01:2021/10/01[dp])"
        self.term_query = "{}+AND+{}".format(
            SARS_CoV_2_query_term, date_query_term)
        return self.virus_url.format(self.term_query)

    def _SARS_CoV_Url(self):
        SARS_CoV_query_term = "(SARS[tiab]+or+SARS-CoV[tiab]+or+severe+acute+respiratory+syndrome+coronavirus[tiab])"
        date_query_term = "2003/01/01:2019/12/31[dp]"
        self.term_query = "{}+AND+{}".format(SARS_CoV_query_term, date_query_term)
        return self.virus_url.format(self.term_query)


    def _MERS_CoV_Url(self):
        self.term_query = "MERS-CoV[tiab]+or+middle+east+respiratory+syndrome+coronavirus[tiab]"
        return self.virus_url.format(self.term_query)

    def _IV_Url(self):
        self.term_query = "(IAV[tiab]+or+influenza+virus[tiab]+or+influenza+A+virus[tiab]+or+influenza+B+virus[tiab]+or+influenza+C+virus[tiab])"
        return self.virus_url.format(self.term_query)
    
    def _HSV_1_Url(self):
        self.term_query = "(HSV-1[tiab]+or+herpes+simplex+virus+type+1[tiab])"
        return self.virus_url.format(self.term_query)

    def _HSV_2_Url(self):
        self.term_query = "(HSV-2[tiab]+or+herpes+simplex+virus+type+2[tiab])"
        return self.virus_url.format(self.term_query)

    def _VZV_Url(self):
        self.term_query = "VZV[tiab]+varicella+zoster+virus[tiab]"
        return self.virus_url.format(self.term_query)
        
    def _HCMV_Url(self):
        self.term_query = "HCMV[tiab]+or+human+cytomegalovirus[tiab]"
        return self.virus_url.format(self.term_query)

    def _KSHV_Url(self):
        self.term_query = "KSHV[tiab]+or+Kaposi's+sarcoma-associated+herpesvirus[tiab]"
        return self.virus_url.format(self.term_query)



class MakeDirs:
    ROOT = 'Virus'
    dirnames = ['', 'all_abstract', 'output', 'pubtator', 'pubtator/tiab', 'pubtator/xml', 'with_mutation']
    def __init__(self, root):
        self.root = self.ROOT+'/'+root
    def makeDirs(self):
        for dirname in self.dirnames:
            fulldir = '{0}/{1}'.format(self.root, dirname)
            if os.path.exists(fulldir):
                continue
            else:
                os.mkdir(fulldir)

                
if __name__ == "__main__":
    
    VirusPaperAll('IV').ID()
    MakeDirs('IV').makeDirs()