'''
The difference between rules and regex, is that rules are used to judge whether a mutation
found by tmVar is a correct mutation, and try to normalize it, while regex are used to 
look for unknown mutations from the text.'''
import re
from superior.package.sub3.dependency import MutationWord, ReplaceAminoAcids, getFullWord
class RegExSentence:
    aa_single_wt = r'[AC-IK-NP-TV-Y]'
    aa_single_mt = r'([*]|STOP|stop|Stop|St\.|[AC-IK-NP-TV-Y])'
    aa_site = r'[1-9]\d{0,2}'
    aa_Pattern = r'\b{0}(/{0}){{0,}}({1})({2})(/{2}){{0,}}\b'.format(aa_single_wt, aa_site, aa_single_mt)
    nt_letter = r'[ACGT]'
    nt_site = r'[1-9]\d{2,6}'
    nt_Pattern = r'\b{0}(/{0}){{0,}}({1}){0}(/{0}){{0,}}\b'.format(nt_letter, nt_site)
    FigureSerial = r'(\bS\d|\bS1\d)[A-H]\b[^/]'
    patternList = [
        re.compile(aa_Pattern),
        re.compile(nt_Pattern)
        ]
    def __init__(self, sentence:str, patternList:list=[]):
        self.sentence = sentence
        self.patternList = [re.compile(pattern) for pattern in list(set(patternList))]

    def _processMatch(self, match:str):
        return ('' 
            if (len(match) == 3 or len(match) == 4) and re.search(self.FigureSerial, match)
            else ReplaceAminoAcids().replaceAminoAcids(match)
            )

    def _getFullword(self, span:tuple):
        return getFullWord(self.sentence, span)

    def _getVarMatch(self, found_length:int, pattern:re.Pattern):
        '''different in each class'''
        varDict = {}
        span,loc  = (0, 0), 0
        for i in range(found_length):
            loc = span[1]
            span = pattern.search(self.sentence[loc:]).span()
            span = tuple([loc+s for s in span])
            match = self.sentence[span[0]:span[1]]
            new = self._processMatch(match)
            if new == '': continue
            full_span, fullword = self._getFullword(span)
            key = (match, span, fullword, full_span)
            if 'SSS' not in new:
                varDict[key] = MutationWord(new).getSplitMutationWordList()
                # print(match, self.sentence, pattern)
            else:
                varDict[key] = [new]
        return varDict

    def _getVarList(self):
        '''same in every class'''
        self.varDict = {}
        for pattern in self.patternList:
            found = pattern.findall(self.sentence)
            if found == []: continue
            else: self.varDict.update(self._getVarMatch(len(found), pattern))

    def getVarDict(self):
        self._getVarList()
        return self.varDict

class AminoAcidsRegEx(RegExSentence):
    """zero width negative assertion: (?![-\d])"""
    aa_letter = RegExSentence.aa_single_wt
    aa_long = r'A(la|rg)|As[np]|Cys|Gl[nuy]|His|Ile|L(eu|ys)|Met|P(he|ro)|Ser|T[hy]r|Trp|Val'
    aa_after_long = r'{}|Stop'.format(aa_long)
    aa_site = r'[1-9]\d{1,2}'
    aa_short = r'[ACGHILMPSTV][aehilrsy][aeglnoprstuy]'
    # one_Pattern = r'(\b{3}|{0})(-?{1}-?)({2})(?![-\d])'.format(aa_long, aa_site, aa_after_long, aa_letter)
    aa_after_short = r'{}|Stop'.format(aa_short)
    multi_Pattern = r'{0}(/{0}){{1,}}(-?{1}-?)({2})(/{2}){{1,}}'.format(aa_short, aa_site, aa_after_short)
    one_Pattern = r'(\b{0})(-?{1}-?)({2})(?![-\d])'.format(aa_short, aa_site, aa_after_short, aa_letter)
    one_Pattern_Upper = r'(\b{3})(-?{1}-?)({2})(?![-\d])'.format(aa_short, aa_site, aa_after_short, aa_letter)

    aa_letter = RegExSentence.aa_single_wt
    stop_codon = r'STOP|Stop|stop|St\.|\*'
    stop_Pattern = r'{0}{1} ?({2})'.format(aa_letter, aa_site, stop_codon)
    
    arrow_Pattern = r'({1}-?{0}|{0}{1})--(&gt;|>){0}'.format(aa_letter, aa_site)
    arrow_Pattern_2 = r'({0})-{1} --(&gt;|>) ({0})'.format(aa_long, aa_site)
    connector = r'[- ]to[- ]|--(&gt;|>)'
    oneSite_Pattern = r'({0})(\(?)({1})(\))({2})({0})'.format(aa_short, aa_site, connector)
    doubleSite_Pattern = r'({0})(?P<site>{1})({2})({0})(?P=site)'.format(aa_short, aa_site, connector)
    '''555-V--(&gt;|>)I
    Leu(526)-to-Met
    F514--(&gt;|>)L, L528--(&gt;|>)M
    C144Arg and C144Lys'''
    patternList = [
            re.compile(multi_Pattern, re.IGNORECASE), 
            re.compile(one_Pattern, re.IGNORECASE),
            re.compile(one_Pattern_Upper),
            re.compile(stop_Pattern, re.IGNORECASE),
            re.compile(arrow_Pattern),
            re.compile(arrow_Pattern_2, re.IGNORECASE),
            re.compile(oneSite_Pattern, re.IGNORECASE),
            re.compile(doubleSite_Pattern, re.IGNORECASE),
            re.compile(RegExSentence.aa_Pattern)
            ]

    def __init__(self, sentence:str):
        self.sentence = sentence

    def _processMatch(self, match:str):
        site = re.findall(self.aa_site, match)
        new = ReplaceAminoAcids.replaceAminoAcids(match)
        new = re.sub(r'to|for|with|intron', '', new)
        return new
        

class NucleotideRegEx(RegExSentence):

    nt_letter, nt_site = RegExSentence.nt_letter, RegExSentence.nt_site
    nt_pattern = r'nt{0}{1}{0}(/{0}){{0,}}'.format(nt_letter, nt_site)
    dash_pattern = r'{1} \({0}\-{0}\)'.format(nt_letter, nt_site)
    nala_pattern = r'{0}[- ]to[- ]{0} .*?at .*?nt.*?{1}'.format(nt_letter, nt_site)
    arrow_pattern = r'{1}[a-z\( ]*?{0}--(&gt;|>){0}[\)]?'.format(nt_letter, nt_site)
    '''G-to-A mutation at nucleotide (nt) 1896
    1764(G--(&gt;|>)T)/1766(C--(&gt;|>)G)'''
    patternList = [
        re.compile(nt_pattern),
        re.compile(nala_pattern),
        re.compile(arrow_pattern),
        re.compile(dash_pattern),
        re.compile(RegExSentence.nt_Pattern)
        ]


    def __init__(self, sentence:str):
        self.sentence = sentence

    def _processMatch(self, match:str):
        if re.findall('[ -]', match):
            letters = re.findall(self.nt_letter, match)
            site = re.findall(self.nt_site, match)
            return site[0].join(letters)
        else:
            return match.replace('nt', '')

class EBVlower(RegExSentence):

    nt_letter = '[acgt]'
    nt_site = '[1-9]\d{3,5}'
    nt_Pattern = '{0}{1}{0}'.format(nt_letter, nt_site)

    def __init__(self, sentence, patternList=[]):
        super().__init__(sentence, patternList=patternList)
        self.patternList = [re.compile(self.nt_Pattern),]

    def _processMatch(self, match:str):
        return match.upper()

class HIVjoint(RegExSentence):

    aa_letter = RegExSentence.aa_single_wt
    aa_site = RegExSentence.aa_site
    aa_multi = r'{0}{{1,7}}'.format(aa_letter)
    aa_joint = r'\b({0}{1}{2})((/?{2}){{1,}})\b'.format(aa_letter, aa_site, aa_multi)

    def __init__(self, sentence, patternList=[]):
        super().__init__(sentence, patternList=patternList)
        if "motif" in sentence:
            self.sentence = ""
        self.patternList = [re.compile(self.aa_joint),]
    def _processMatch(self, match: str):
        if 'RT' in match:
            return match.replace('RT','')
        else:
            return match

class HPVdash(RegExSentence):
    letter = RegExSentence.aa_single_mt
    dash = r"[1-9]\d{1,3}[A-Z]-[A-Z]"

    def __init__(self, sentence, patternList=[]):
        super().__init__(sentence, patternList=patternList)
        self.patternList = [re.compile(self.dash),]

    def _processMatch(self, match:str):
        site = re.findall(r'\d{1,4}', match)[0]
        letters = re.findall(r'[A-Z]', match)
        return site.join(letters)

class Affix(RegExSentence):
    AFFIX = []
    aa_pattern = r'{0}(/{0}){{0,}}{1} ?{2}(/{2}){{0,}}\b'.format(
        RegExSentence.aa_single_wt, RegExSentence.aa_site, RegExSentence.aa_single_mt)
    
    def __init__(self, sentence:str):
        super().__init__(sentence)
        self.patternList = {}
        for A in self.AFFIX:
            self.patternList[A] = re.compile(r'\b{0}{1}'.format(A, self.aa_pattern))

    def _getVarMatch(self, le:int, A:str, pattern):
        varDict = {}
        span, loc  = (0, 0), 0
        for i in range(le):
            loc = span[1]
            span = pattern.search(self.sentence[loc:]).span()
            span = [loc+s for s in span]
            match = self.sentence[span[0]:span[1]]
            if match.startswith(A):
                span = tuple([span[0]+len(A), span[1]])
            elif match.endswith(A):
                span = tuple([span[0], span[1]-len(A)])
            full_span, fullword = self._getFullword(span)
            new = self.sentence[span[0]:span[1]]
            new = re.sub(r'STOP|stop|Stop|St\.|\*', 'X', new)
            mutationWord = MutationWord(new)
            key = (A, match, span, fullword, full_span)
            varDict[key] = mutationWord.getSplitMutationWordList()
        return varDict

    def _getVarList(self):
        self.varDict = {}
        for A, pattern in self.patternList.items():
            found = pattern.findall(self.sentence)
            if found == []: continue
            else: self.varDict.update(self._getVarMatch(len(found), A, pattern))

    def getVarDict(self):
        self._getVarList()
        return self.varDict


class HBVaffix(Affix):
    AFFIX = [
        r"s",r"x",r"c",r"rt",
        r"S",
        r"preS1",r"preS2",r"ps",r'ps2',r'ps1',
        r"preS1 ",r"preS2 ",r"ps ",r'ps2 ',r'ps1 ',
        r'preC-',r'C-',r'S-',r'X-',r'PC-',
        r"core",r"pol",r"sp",r"pc",r"RT",
        r"core ",r"pol ",r"sp ",r"pc ",r"RT ",]

class HIVaffix(Affix):
    AFFIX = [
        r"PR",r"RT",r"PI",
        r"PR ",r"RT ",r"PI "]
    SUFFIX = [
        r'rev',
    ]
    aa_pattern = HIVjoint.aa_joint
    
    def __init__(self, sentence: str):
        super().__init__(sentence)
        for S in self.SUFFIX:
            self.patternList[S] = re.compile(r'\b{1}{0}'.format(S, self.aa_pattern))

class HTLV1affix(Affix):
    AFFIX = [r"Tax", r"Tax "]

class HPVaffix(Affix):
    AFFIX = [
        r'E6',r'E7',r'URR',
        r'E6 ',r'E7 ',r'URR ']

class VirusRegex:

    def __init__(self, sentence:str):
        self.sentence = sentence
        self.varDict = {}

    def getRegex(self):
        pass

class HBVregex(VirusRegex):

    def getRegex(self):
        return HBVaffix(self.sentence).getVarDict()
        

class EBVregex(VirusRegex):

    def getRegex(self):
        return EBVlower(self.sentence).getVarDict()


class HIVregex(VirusRegex):

    def getRegex(self):
        self.varDict.update(HIVaffix(self.sentence).getVarDict())
        self.varDict.update(HIVjoint(self.sentence).getVarDict())
        return self.varDict


class HPVregex(VirusRegex):

    def getRegex(self):
        self.varDict.update(HPVaffix(self.sentence).getVarDict())
        self.varDict.update(HPVdash(self.sentence).getVarDict())
        return self.varDict

class HTLV1regex(VirusRegex):

    MutationDict = {"Tax M47": ("Tax", ["L319R", "L320S"])}

    def _getKnownMutation(self):
        for expression, information in self.MutationDict.items():
            match = re.findall(expression, self.sentence, re.IGNORECASE)
            if match:
                span = re.search(expression, self.sentence, re.IGNORECASE).span()
                fullspan = span
                fullword = match[0]
                gene, mutation = information
                self.varDict[(gene, match[0], span, fullword, fullspan)] = mutation

    def getRegex(self):
        self._getKnownMutation()
        self.varDict.update(HTLV1affix(self.sentence).getVarDict())
        return self.varDict


class PatternFromMutation:

    sub_dict = {
    r'(?P<s>[\(\)\-\+\.\>])': r'\\\g<s>',
    r'[1-9][0-9]{0,3}': '[1-9][0-9]{0,3}',
    r'[A-Z]': '[A-Z]',
    # r'[acgt]': '[acgt]',
    }

    def __init__(self, varWord:str):
        self.varWord = varWord

    def getPattern(self):
        varWord = self.varWord
        if re.findall(r'ins|del|Delta|fs|Non|intron|c\.|p\.|transition|at|bp|pc|rt|pr', varWord):
            return
        elif ReplaceAminoAcids.hasAminoAcids(varWord):
            return
        for k, v in self.sub_dict.items():
            if varWord.isupper() or 'acgt' not in v:
                varWord = re.sub(k, v, varWord)
            elif 'acgt' in v:
                varWord = re.sub(k, v, varWord, count=1)
        return r'\b{}\b(?!\/)'.format(varWord)
if __name__ == "__main__":
    RegExSentence('')