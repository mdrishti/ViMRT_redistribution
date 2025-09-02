import re
import nltk
from . import dependency, regularexpression, rules
class VarRules:
    def __init__(self, description):
        self.description = description
    def identify_rule(self):
        rules.Rules.meets_condition_pre(self.description)
        try:
            rule_cls = next(
                (
                    rule_cls
                    for rule_cls in rules.Rules.__subclasses__()
                    if rule_cls.meets_condition(self.description)
                ),
                rules.UnknownRule,
            )
        except RuntimeError as e:
            return rules.UnknownRule
        return rule_cls(self.description)

class VirusPaper:
    def __init__(self, description:dict):
        self._pmid = description['pmid']
        self._passage = description['passage']
        assert type(self._passage) == dict
        self._annotation = description['annotation']
        self._useinput =description["useinput"]
        self.regex = dependency.regexColumns()
        self.raw = dependency.tmVarColumns()
        self.rules = dependency.rulesColumns()
        self.total = {'rules': self.rules, 'raw': self.raw, 'regex': self.regex}
        self._loc = 0
    @staticmethod
    def meets_condition(virusName:str):
        return False
    def all(self):
        self._processAnnotation()
        self._getRegex()
        for k, v in self.total.items():
            if set((self.virus,self._useinput))=={"Unknown"}:
                v['virus'] ="Unknown"
            else:
                v['virus'] = ";".join({self.virus,self._useinput}-{'Unknown'})
            v['pmid'] = self._pmid
    def _processAnnotation(self):
        i, j = 0, 0
        for location, infoList in self._annotation.items():
            self.patternList = []
            for info in infoList:
                expression, standard, *rest = info
                if re.findall(r'rs\d+', standard, re.IGNORECASE): continue
                information = (expression, standard, *location, *rest)
                i = self._appendRaw(i, information)
                j = self._appendRules(j, information)
            self._processSentenceWithMutation(location)
    def _processSentenceWithMutation(self, location):
        sentence, *_ = location # offset, region
        self.patternList = list(filter(None, self.patternList))
        if self.patternList == []: return
        result = regularexpression.RegExSentence(sentence, self.patternList).getVarDict()
        if result == {}: return
        self._appendRegex(result, location)    
    def _appendRaw(self, i:int, information):
        expression, standard, *rest = information
        var = rules.PubTatorWord(standard)
        try:
            for p in var.getIdentifier():
                self.raw.loc[i] = [p, expression, standard, *rest]
                i += 1
        except Exception as e:
            print(e, expression, rest)
        return i
    def _appendRules(self, i:int, information):
        expression, standard, sentence, offset, region, *rest = information
        if not re.findall(r'\d', expression): return i
        typ, start, end = rest
        span = (start, end)
        description = {
                'virus': self.virus,
                'varWord': expression,
                'sentence': sentence,
                'span': span,
            }
        rule_cls = VarRules(description).identify_rule()
        processed = rule_cls.getRight()
        try:
            fullword = rule_cls.fullWord
            fullspan = rule_cls.fullspan
            # fullspan = [int(offset)+s for s in fullspan]
        except Exception as e:
            print(e, rule_cls, processed, description)
        # print(rule_cls, processed, description)
        try:
            for p in processed:
                if p == '': continue
                self.rules.loc[i] = [p, expression, standard, *span, sentence, offset, region, typ, fullword, *fullspan]
                self.patternList.append(regularexpression.PatternFromMutation(expression).getPattern())
                i += 1
        except Exception as e:
            print(e, information)
        return i
    def _getRegex(self):
        for info, text in self._passage.items():
            region, offset = info
            sentenceList = [sen for sen in text.split('. ') if sen!= '']
            offsetList = [text.find(sentence) + offset for sentence in sentenceList]
            list_len = len(sentenceList)
            for i in range(list_len):
                sentence = sentenceList[i]
                sentence = sentence if sentence.endswith('.') else sentence + '.'
                location = (sentence, offsetList[i], region)
                self._processSentence(location)
    def _processSentence(self, location):
        sentence, *_ = location
        Aa = regularexpression.AminoAcidsRegEx(sentence).getVarDict()
        Nt = regularexpression.NucleotideRegEx(sentence).getVarDict()
        result = Aa
        result.update(Nt)
        self._appendRegex(result, location)# location: region, typ
    def _appendRegex(self, result:dict, location:tuple):
        sentence, offset, region = location #rest: sentence, region
        for key, varList in result.items():
            assert type(key) == tuple, "{key} is not a tuple"
            if len(key) == 5:
                gene, match, span, fullword, fullspan = key
            elif len(key) == 4:
                match, span, fullword, fullspan = key
                gene = ""
            else:
                print(key, location)
            # paper_span = [int(offset)+s for s in span] # span: start, end
            # fullspan = [int(offset)+s for s in fullspan]
            for var in varList:
                row = [var, match, *span, sentence, offset, region, gene, fullword, *fullspan]
                self.regex.loc[self._loc] = row# *paper_span
                self._loc += 1
    

class HBVPaper(VirusPaper):
    virus = 'HBV'
    @staticmethod
    def meets_condition(virusName:str):
        return virusName == "HBV"   
    def _processSentence(self, location):
        super()._processSentence(location)
        sentence, *_ = location
        result = regularexpression.HBVregex(sentence).getRegex()
        self._appendRegex(result, (location))         
    

class HIVPaper(VirusPaper):
    virus = 'HIV'
    @staticmethod
    def meets_condition(virusName:str):
        return virusName == "HIV"
    def _processSentence(self, location):
        super()._processSentence(location)
        sentence, *_ = location
        result = regularexpression.HIVregex(sentence).getRegex()
        self._appendRegex(result, (location))

class HPVPaper(VirusPaper):
    virus = 'HPV'
    @staticmethod
    def meets_condition(virusName:str):
        return virusName == "HPV"
    def _processSentence(self, location):
        super()._processSentence(location)
        sentence, *_ = location
        result = regularexpression.HPVregex(sentence).getRegex()
        self._appendRegex(result, (location))

class EBVPaper(VirusPaper):
    virus = 'EBV'
    @staticmethod
    def meets_condition(virusName:str):
        return virusName == "EBV"
    def _processSentence(self, location):
        super()._processSentence(location)
        sentence, *_ = location
        result = regularexpression.EBVregex(sentence).getRegex()
        self._appendRegex(result, (location))

class HTLV1Paper(VirusPaper):
    virus = 'HTLV1'
    @staticmethod
    def meets_condition(virusName:str):
        return virusName == "HTLV1"
    def _processSentence(self, location):
        super()._processSentence(location)
        sentence, *_ = location
        result = regularexpression.HTLV1regex(sentence).getRegex()
        self._appendRegex(result, (location))

class XMRVPaper(VirusPaper):
    virus = 'XMRV'
    @staticmethod
    def meets_condition(virusName:str):
        return virusName == "XMRV"
class IVPaper(VirusPaper):
    virus = 'SARS-CoV-2'
    @staticmethod
    def meets_condition(virusName:str):
        return virusName == "SARS-CoV-2"
class IVPaper(VirusPaper):
    virus = 'IV'
    @staticmethod
    def meets_condition(virusName:str):
        return virusName == "IV"
class UnknownPaper(VirusPaper):
    """unknown virus of paper"""
    virus = 'ALL'

class NewPaper(VirusPaper):
    @staticmethod
    def meets_condition(virusName:str):
        NewPaper.virus = virusName
        return 'V' in virusName

class VirusName:
    def __init__(self, name):
        self.name = name
        print (name)
    def identify_virus(self):
        try:
            virus_cls = next(
                (
                    rule_cls
                    for rule_cls in VirusPaper.__subclasses__()
                    if rule_cls.meets_condition(self.name)
                ),
                UnknownPaper,
            )
        except RuntimeError as e:
            return UnknownPaper
        return virus_cls