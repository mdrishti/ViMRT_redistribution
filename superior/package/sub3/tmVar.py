import re
from functools import reduce
import xml.etree.ElementTree as ET
class SentenceOfWordInText:
    breakPattern = re.compile(r'')
    def __init__(self):
        pass
    @classmethod
    def getSentence(cls, text:str, loc:int, leg:int, offset:int=0):
        """get whole sentence from a passage;
        based on the location of mutation;
        offset is the location of text in the whole paper
        """
        if ". " not in text:
            return text
        assert ". " in text, f"{text} is one sentence."
        s, e = loc + leg - offset, loc - offset
        start = text.rfind('. ', 0, s) + 2
        if start == 1:
            start = 0
        end = text.find('. ', e, len(text)) + 1
        if end == 0:
            end = len(text)
        return text[start:end]
class LocationInPaper:
    @staticmethod
    def getRegion(location:int, offset:list):
        """
        >>> length = [10, 20, 30, 40, 70]
        >>> offset = 55
        >>> loc = getRegion(offset, length)
        >>> 4
        """
        return (0 
            if location < offset[0] 
            else LocationInPaper.getRegion(location, offset[1:]) + 1)

class NLP:
    '''defines regex func used by all other classes
    Parameters
    ----------
    : a file with annotated mutations
    '''
    def __init__(self, filename:str):
        assert filename != '', f'{filename} is empty'
        self.filename = filename
        self.annotation = {}
        self.passage = {}
    def _getPmid(self): pass
    def all(self):pass

class TXT(NLP):
    def __init__(self, filename:str):
        super().__init__(filename)
        if len(filename) > 100:
            self._lines = filename.split('\n')
        else:
            with open(filename,'r') as f:
                self._lines = f.read().strip().split("\n")

    def all(self):
        self._getPmid()
        if self.pmid == "": return "error, no pmid"
        self._getText()

    def _getPmid(self):
        self.pmid = self._lines[0].split("|")[0]

    def _getText(self):
        self._textList = []
        for line in self._splitPassageAnnotation():
            try:
                pmid, region, text = line.split("|",2)
                self._textList.append(text)
                if region == 't' or region == 'title': # {'t': <title>, 'a': <abstract>}
                    offset = 0
                
                self.passage[(region, offset)] = text
                offset = len(text) + 1
            except:
                print("error: ", self.filename, self.pmid, self._lines)
                return

    def _splitPassageAnnotation(self):
        for i in range(len(self._lines)):
            if self._lines[i].startswith(self.pmid+'\t'):
                self._getOffsetList()
                self._getAnnotation(self._lines[i:])
                break
            yield self._lines[i]

    def _getOffsetList(self):
        '''get the offset of each region'''
        lines = self._textList
        ls = [len(line) for line in lines]# list of ints
        t = [lambda x: reduce(lambda x,y: x+y+1, ls[:x]) for i in range(1, len(ls)+2)]# generate functions to add lengths
        self._offset = [t[i](i) for i in range(1, len(ls)+1)]# a list of length accumulative sums

    def _getAnnotation(self, lines:list):
        self._text = " ".join(self._textList)
        for line in lines:
            try:
                start, end, epr, typ, std = line.split('\t')[1:]
            except Exception as e:
                print (e, line, self.pmid)
                return
            start, end = int(start), int(end)
            region = list(self.passage.keys())[LocationInPaper.getRegion(start, self._offset)][0]
            sentence = SentenceOfWordInText.getSentence(self._text, start, end-start)
            offset = self._text.find(sentence)
            sentence = sentence.strip(' ')
            key = (sentence, offset, region)
            information = (epr, std, typ, start-offset, end-offset)
            if key in self.annotation:
                self.annotation[key].append(information)
            else:
                self.annotation[key] = [information]# The element is a tuple

class XML(NLP):
    _abandon = re.compile(r"ABBR|COMP_INT|AUTH_CONT|SUPPL|ACK_FUND")
    def __init__(self, filename):
        super().__init__(filename)
        if type(filename) == str:
            self.pmid = filename.rsplit('/',1)[1].split('.')[0]
            xml = ET.parse(filename).getroot()
            self._document = xml.find('document')
            if self._document == None:
                self._document = xml
        else:
            self._document = filename
            info = self._document[1]
            for i in info:
                print (i)
                if i.attrib["key"] == "article-id_pmid": 
                    self.pmid = i.text
                    break
        try:
            self._length = len(self._document)
        except TypeError as e:
            print(e, filename)

    def all(self):
        self._getAnnotated()
        self._getPassage()

    def _getAnnotated(self):
        '''get passages that contain annotations'''
        self._annotatedDict = {}
        for i in range(1, self._length):
            for j in range(len(self._document[i])):
                node = self._document[i][j]
                if node.tag == "annotation":
                    self._annotatedDict[i] = j
                    break

    def _getPassage(self):
        self._offset = []
        textDict = {}
        for i in range(1, self._length):
            node = self._document[i]
            region = ""
            for child in node:
                tag, attrib = child.tag, child.attrib
                if attrib.get("key") == "section_type":
                    region = child.text
                elif attrib.get("key") == "type" and region == "":
                    region = child.text
                elif tag == "offset":
                    offset = int(child.text)
                elif tag == "text":
                    passage = child.text
            if 'REF' in region: break
            elif self._abandon.findall(region): continue
            if region != 'TABLE':
                region = region[0].lower()
                textDict[(region, offset)] = passage
            else:
                textDict[(region, offset)] = 'TABLE'
            if i not in list(self._annotatedDict): continue
            annotationLoc = self._annotatedDict[i]
            node, annotation = node[:annotationLoc], node[annotationLoc:]
            for result in self._processNode(annotation, passage, offset, region):
                sentence, region, *rest, start, end = result
                sen_offset = passage.find(sentence) + offset
                key = (sentence, sen_offset, region)
                value = (*rest, start-sen_offset, end-sen_offset)
                self.annotation[key] = self.annotation[key] + [value] if key in self.annotation else [value]
        self.passage = textDict
        
    def _processNode(self, root, passage, offset:int, region):
        for node in root:
            standard, typ, expression = [child.text for child in (node[0], node[1], node[3])]
            varOffset, length = [int(x) for x in node[2].attrib.values()]
            if region == 'TABLE': 
                sentence = region
            else:
                sentence = SentenceOfWordInText.getSentence(passage, varOffset, length, offset)
                sentence = sentence.strip()
            yield (sentence, region, expression, standard, typ, varOffset, varOffset+length)

