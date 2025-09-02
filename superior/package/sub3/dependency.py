import re
import pandas as pd
aa_triple = '''
(Ala|Arg|As[np]|Cys|Gl[nuy]|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|T[hy]r|Trp|Val)
'''
aa_exp = '''
[ACGHILMPSTV][aehilrsy][aeglnoprstuy]'''
def tmVarColumns():
    return pd.DataFrame(columns=["processed","original","standard",
            "sentence","sen_offset","region","type","start_in_sen", "end_in_sen"])
def rulesColumns():
    return pd.DataFrame(columns=["processed","original","standard",
            "start","end",
            "sentence","sen_offset","region","type",#"start":[], "end":[],
            "fullword","start_in_sen","end_in_sen"])
def regexColumns():
    return pd.DataFrame(columns=["processed","original",
            "start","end",
            "sentence","sen_offset","region","gene",#"start":[], "end":[],
            "fullword","start_in_sen","end_in_sen"])
def AminoAcids():
    '''
    Not included: J, O, U
    Possible degeneration: Z, B
    '''
    return {
        'Alanine':'A','alanine':'A','ALANINE':'A','Ala':'A','ala':'A','ALA':'A',
        'Arginine':'R','arginine':'R','ARGININE':'R','Arg':'R','arg':'R','ARG':'R',
        'Asparagine':'N','asparagine':'N','ASPARAGINE':'N','Asn':'N','asn':'N','ASN':'N',
        'Asparticacid':'D','asparticacid':'D','ASPARTICACID':'D','Asp':'D','asp':'D','ASP':'D',
        'Cysteine':'C','cysteine':'C','CYSTEINE':'C','Cys':'C','cys':'C','CYS':'C',
        'Glutamine':'Q','glutamine':'Q','GLUTAMINE':'Q','Gln':'Q','gln':'Q','GLN':'Q',
        'Glutamicacid':'E','glutamicacid':'E','GLUTAMICACID':'E','Glu':'E','glu':'E','GLU':'E',
        'Glycine':'G','glycine':'G','GLYCINE':'G','Gly':'G','gly':'G','GLY':'G',
        'Histidine':'H','histidine':'H','HISTIDINE':'H','His':'H','his':'H','HIS':'H',
        'Isoleucine':'I','isoleucine':'I','ISOLEUCINE':'I','Ile':'I','ile':'I','ILE':'I',
        'Leucine':'L','leucine':'L','LEUCINE':'L','Leu':'L','leu':'L','LEU':'L',
        'Lysine':'K','lysine':'K','LYSINE':'K','Lys':'K','lys':'K','LYS':'K',
        'Methionine':'M','methionine':'M','METHIONINE':'M','Met':'M','met':'M','MET':'M',
        'Phenylalanine':'F','phenylalanine':'F','PHENYLALANINE':'F','Phe':'F','phe':'F','PHE':'F',
        'Proline':'P','proline':'P','PROLINE':'P','Pro':'P','pro':'P','PRO':'P',
        'Serine':'S','serine':'S','SERINE':'S','Ser':'S','ser':'S','SER':'S',
        'Threonine':'T','threonine':'T','THREONINE':'T','Thr':'T','thr':'T','THR':'T',
        'Tryptophan':'W','tryptophan':'W','TRYPTOPHAN':'W','Trp':'W','trp':'W','TRP':'W',
        'Tyrosine':'Y','tyrosine':'Y','TYROSINE':'Y','Tyr':'Y','tyr':'Y','TYR':'Y',
        'Valine':'V','valine':'V','VALINE':'V','Val':'V','val':'V','VAL':'V',
        'Stop':'X', 'STOP':'X','stop':'X','St.':'X','*':'X','Ter':'X',
   }
def Bases():
    return {
        'W': ['A','T'],
        'R': ['A','G'],
        'V': ['A','C','G'],
        'S': ['C','G'],
        'M': ['A','C'],
        'K': ['G','T'],
        'Y': ['C','T'],
        'B': ['C','G','T'],
        'D': ['A','G','T'],
        'H': ['A','C','T'],
        'N': ['A','C','G','T']
    }
def order(string:str):
    current_s = string[0]
    for s in string:
        if ord(s) > ord(current_s):
            return False
        else:
            current_s = s
    return True

def getFullWord(sentence:str, span:tuple):
    s, e = span
    try:
        start = s - re.search(r',(?=\))|\((?=[ \+])|[\[\{](?= )|[\]\[: ]', sentence[s::-1]).span()[0] + 1
    except Exception as E:
        start = sentence.rfind(r' ',0,s)
        if start == -1:
            start = 0
    try:
        end = re.search(r'(?<![\]])[\)\}](?=[ ,\.\]])|[,;\}](?= )|[\]\[\(: ]', sentence[e:]).span()[0] + e
    except Exception as E:
        end = sentence.find(r' ',e,len(sentence))
        if end == -1:
            end = len(sentence)-1
    if start != s or end != e:
        pass
        # print(sentence[s:e], sentence[start:end], sentence)
    full_span = (start, end)
    return full_span, sentence[start:end]


class MutationWord:
    '''
    a class for CLASSIC mutation forms <capital letter><digit><capital letter>
    '''
    DIGITPATTERN = re.compile(r'\d{1,}')# the digit can be very long
    CAPITALPATTERN = re.compile(r'[AC-IK-NP-TV-Y\*]')# the letter only exist individually
    @staticmethod
    def possible(string:str):
        uppers = re.findall(r'[A-Z]', string)
        return (
            (len(uppers) > 1 and
            (
                bool(MutationWord.DIGITPATTERN.findall(string)) 
                or order(''.join(uppers)) 
            )) 
            or re.findall(r'fs|Fs|X|p\.|c\.', string) 
            or re.findall(r'[ACGT]\d{1,4}[acgt]', string)
            )
    @staticmethod
    def fitPattern(string:str):
        return bool(
            re.findall(r'[AC-IK-NP-TV-Y\*]\d{1,}[AC-IK-NP-TV-Y\*]', string)
            )
    def __init__(self, string:str):
        self.varWord = string
        assert self.possible(string), f"{string} is not a mutation"
    def _getSiteDigit(self):
        if not self.DIGITPATTERN.findall(self.varWord):
            self.digit = ""
        else:
            self.digit = self.DIGITPATTERN.search(self.varWord).group()
        return self.digit # the first digit string found in varWord
    def _getCapitalLettersList(self):
        if self._getSiteDigit() == "":
            return list(self.varWord[:2])
        if self.varWord.isupper():
            return [self.CAPITALPATTERN.findall(string) for string in self.varWord.split(self.digit, 1)]
        if ReplaceAminoAcids().hasAminoAcids(self.varWord):
            string = ReplaceAminoAcids().replaceAminoAcids(self.varWord)
            return [self.CAPITALPATTERN.findall(s) for s in string.split(self.digit, 1)]
        if 'Delta' in self.varWord:
            self.varWord = self.varWord.replace('Delta','-')
        return [s.upper() for s in self.varWord.split(self.digit, 1)]
    def getSplitMutationWordList(self):
        # all Capital letters BEFORE digit in the FIRST element, AFTER the SECOND
        beforeLetters , afterLetters = self._getCapitalLettersList()
        for bL in beforeLetters:
            for aL in afterLetters:
                if aL.isalpha():
                    pass
                elif aL == '-':
                    aL = 'Delta'
                elif aL == '/':
                    continue
                yield bL+self.digit+aL


class ReplaceAminoAcids:
    def __init__(self):
        pass
    @staticmethod
    def replaceAminoAcids(varWord:str):
        for key in AminoAcids():
            if key in varWord:
                if key in ['ala', 'alanine']:
                    varWord = re.sub(r'\b'+key, AminoAcids()[key], varWord)
                else:
                    varWord = varWord.replace(key, AminoAcids()[key])
        return varWord
    @staticmethod
    def hasAminoAcids(varWord:str):
        for key in AminoAcids():
            if key in varWord:
                return True
    @staticmethod
    def matchAminoAcids(varWord:str):
        aa_long = r'Ala|Arg|As[np]|Cys|Gl[nuy]|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val'
        aa_after_long = r'{}|Stop'.format(aa_long)
        aa_site = r'[1-9]\d{1,2}'
        one_Pattern = r'(\b{0})(-?{1}-?)({2})(?![-\d])'.format(aa_long, aa_site, aa_after_long)
        return bool(re.findall(one_Pattern, varWord))
    
