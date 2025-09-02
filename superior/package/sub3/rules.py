import re
from superior.package.sub3.dependency import MutationWord, ReplaceAminoAcids, getFullWord


class DNAWord(MutationWord):
    '''Nucleotides are separated by >, not the digit'''
    def _getCapitalLettersList(self):
        self.varWord = re.sub(r'[cg]\.','',self.varWord)
        self.digit = re.findall(r'\d{1,}', self.varWord)[0]
        if '>' in self.varWord:
            sep = '>'
        else:
            sep = re.findall(r'[^\w\.]', self.varWord)[0]
        return [self.CAPITALPATTERN.findall(string) for string in self.varWord.split(sep)]


class IndelWord:
    def __init__(self, string):
        self.varWord = string

    def getVarWord(self):
        return [self.varWord]


class PubTatorWord:

    def __init__(self, string):
        self.varWord = string.split(';')[0]
        self.varWord = ReplaceAminoAcids.replaceAminoAcids(self.varWord)

    def _preAnnotated(self):
        if re.findall(r'ins|del|dup|Delta', self.varWord):
            return IndelWord(self.varWord).getVarWord()
        elif re.findall(r'c\.|g\.', self.varWord):
            try:
                dnaWord = DNAWord(self.varWord)
            except AssertionError as e:
                return [self.varWord]
            return dnaWord.getSplitMutationWordList()
        elif re.findall(r'p.', self.varWord):
            return MutationWord(self.varWord.replace('p.','')).getSplitMutationWordList()
        else:
            return [self.varWord]

    def _newAnnotation(self):
        cgp, varTyp, *rest = self.varWord.split("|")
        varWord = "".join(rest)
        if varTyp == "SUB":
            try:
                return MutationWord(varWord).getSplitMutationWordList()
            except AssertionError as e:
                return [varWord]
        else:
            return [varTyp+varWord]

    def getIdentifier(self):
        if "|" in self.varWord:
            return self._newAnnotation()
        else:
            return self._preAnnotated()

        

class Rules:
    '''Includes a sentence and a variation in the sentence.
    Use some rules to judge whether it's a true mutation.
    '''
    virus_dict = {
        'HIV': ['PR', 'RT', 'Tat', 'Y712HRL', 'A3C'],
        'HPV': ['E1', 'E2', 'E7M', '350G', 'R10'],
        'EBV': ['Z'],
        'HBV': ['PC','PS','BCP','core','rt', 'RT', 'HIV', 'FPolE', 'pc'],
        'HTLV1': [],
        'XMRV': [],
        'HCMV': ['UL97', 'UL54', 'UL32']
    }
    for k, v in virus_dict.items():
        v.append('EGFP')
        virus_dict[k] = v

    @staticmethod
    def meets_condition(description:dict):
        return False
    @staticmethod
    def meets_condition_pre(description:dict):
        assert isinstance(description, dict), f"{description!r} is not a dict"
        for aspect in ('virus', 'varWord', 'sentence'):
            assert aspect in description, f"{aspect} not in {description}"
            assert isinstance(description[aspect], str)

    def __init__(self, description:dict):
        self.description = description
        self.sentence = description['sentence']
        self.varWord = description['varWord']
        self.virus = description['virus']
        self.span = description['span']
        self.fullspan = None
    def __str__(self):
        return self.__class__.__name__

    def _replaceAbbr(self, fullWord=None):
        if fullWord == None:
            fullWord = self.fullWord
        varWord = ReplaceAminoAcids.replaceAminoAcids(self._getFullWord()) 
        return varWord

    def _getStart(self):
        return self.sentence.find(self.varWord)

    def _getEnd(self):
        return self._getStart() + len(self.varWord)

    def _getVarSite(self):
        return re.findall(r'\d{1,}', self.varWord)

    def _getVarLetter(self):
        return re.findall(r'[AC-IK-NP-TV-Y]', self.varWord)

    def _VirusWrongWordsPattern(self):
        wrongWordList = self.virus_dict.get(self.virus, [])
        wrongWordList.append(self.virus)
        pattern = r'\b{}\b'.format('|'.join(wrongWordList))

        return pattern

    def _getFullWord(self):
        if self.sentence == "TABLE":
            self.fullspan = self.span
            self.fullWord = self.varWord
            return self.varWord
        if self.fullspan == None:
            self.fullspan, self.fullWord = getFullWord(self.sentence, self.span)
        s, e = self.fullspan

        return self.sentence[s:e]

    def _fullWordRight(self):
        pattern = self._VirusWrongWordsPattern()
        fullWord = re.sub(pattern, '', self._getFullWord(), re.IGNORECASE)
        self.fullWord = fullWord
        if fullWord == '' or (
            len(self._getFullLetter()) < 2 
            and not ReplaceAminoAcids().hasAminoAcids(self.fullWord)
            and not re.findall(r'ins|del', self.fullWord)):

            return True  
        return False

    def _getFullSite(self, fullWord=None):
        if fullWord == None:
            fullWord = self.fullWord

        return re.findall(r'\d{1,}', fullWord)

    def _getFullLetter(self, fullWord=None):
        if fullWord == None:
            fullWord = self.fullWord

        return re.findall(r'[AC-IK-NP-TV-Y]', fullWord)

    def isWrong(self):
        return self._fullWordRight()

    def getRight(self): 
        return []

    def _split(self):
        return MutationWord(self.varWord).getSplitMutationWordList()
        
    def _hasSlash(self, word):
        slash_split_wordlist = [w for w in word.split('/') if w!= ""]
        allMutation = [MutationWord.fitPattern(w) for w in slash_split_wordlist]

        if all(allMutation):
            new_list = []
            for string in slash_split_wordlist:
                if not re.findall(r'[A-Z]{2,}', string): 
                    new_list.append(string)
                else:
                    for var in MutationWord(string).getSplitMutationWordList():
                        new_list.append(var)
            return new_list

        newWordList = []
        temp = []
        for w in slash_split_wordlist:
            if MutationWord.fitPattern(w):
                if temp == []:
                    newWordList.append(w)
                else:
                    if len(newWordList) > 0:
                        last = newWordList.pop()
                        newWordList.append(last+'/'+'/'.join(temp))
                        newWordList.append(w)
                    else:
                        newWordList.append('/'.join(temp)+'/'+w)
                    temp = []
                    
            elif not MutationWord.possible(w):
                temp.append(w)
        if temp != []:
            last = newWordList.pop()
            newWordList.append(last+'/'+'/'.join(temp))
        result_list = []
        for w in newWordList:
            for var in MutationWord(w).getSplitMutationWordList():
                result_list.append(var)

        return result_list

    def _wordToVarList(self, wordList):
        varList = []
        wrongWordList = self.virus_dict.get(self.virus, [])
        wrongWordList.append(self.virus)
        pattern = re.compile(r'{}'.format('|'.join(wrongWordList)))
        for word in wordList:
            word = pattern.sub('', word)
            if not MutationWord.fitPattern(word): continue
            if '/' in word:
                for var in self._hasSlash(word):
                    varList.append(var)
            else:
                varList.append(word)

        return varList

class UnknownRule(Rules):
    '''unknown rules'''

class NaturalLanguage(Rules):
    wrongWords = [
        'C-A1',
        'C8166 T'
    ]
    @staticmethod
    def meets_condition(description:dict):
        varWord = description['varWord']
        return (
            bool(re.findall(r'to|with|by|[- ]|for|ins|del|Non|Delta|of|at', varWord) 
            and '>' not in varWord)
            )

    def isWrong(self):
        if self.virus == 'HIV' and self.varWord in self.wrongWords:
            return True
        if not re.findall(r'\d', self.varWord):
            return True
        if len(self._getFullSite()) > 1:
            return True
        self.fullWord = ReplaceAminoAcids.replaceAminoAcids(self.fullWord)
        if len(re.findall(r'[A-Z]', self.fullWord)) < 2 and not re.findall(r'ins|del', self.fullWord):
            return True
        return False

    def getRight(self):
        if super().isWrong() or self.isWrong(): 
            return []
        if re.findall(r'ins|del|Non|Delta', self.varWord):
            return [ReplaceAminoAcids().replaceAminoAcids(self.fullWord)]
        elif '&gt' in self.varWord:
            return []
        # elif '>' in self.varWord:
        #     return [self._getVarSite()[0].join(self._getVarLetter())]
        elif '-' in self.varWord:
            pass
        elif ' ' in self.varWord and self.varWord.isupper():
            return self._hasSpace(self.varWord)
        varWord = self._replaceAbbr().replace(',','')
        self.varWord = ReplaceAminoAcids.replaceAminoAcids(self.varWord) if re.findall(r'[\(\)]', varWord) else varWord
        try:
            site = self._getVarSite()[0]
        except Exception as e:
            return []
        letters = self._getVarLetter()
        if len(letters) != 2: return []
        if (re.findall(r'of|rather|for', self.varWord) 
        and not re.findall(r'substituted', self.varWord)):
            return [site.join(letters[::-1])]
        return [site.join(letters)]

    def _hasSpace(self, word): 
        allMutation = [MutationWord.fitPattern(word) for word in word.split(' ')]
        var = MutationWord(word.replace(' ', ''))
        return word.split(' ') if all(allMutation) else var.getSplitMutationWordList()

class MultiVar(Rules):
    '''if re.findall(r'\d{1,}.{1,}\d{1,}|[-\(\)/]', varWord): '''
    
    @staticmethod
    def meets_condition(description:dict):
        pattern = re.compile(r'\d{1,}\D{1,}\d{1,}|[/,]')
        varWord = description.get('varWord')
        return (varWord[-1] != 'S' 
            and bool(pattern.findall(varWord)) 
            and varWord.isupper()
            )

    def _oneSiteVar(self):
        return ('/'.join(self._getVarLetter()) + self._getVarSite()[0] == self.varWord)

    def isWrong(self):
        if not self._fullWordRight():
            return False
        if len(self._getVarSite()) == len(self._getVarLetter()):
            return True
        elif len(self._getFullSite()) == len(self._getFullLetter()):
            return True
        elif len(self._getVarSite()) == 1:
            return self._oneSiteVar()
        else:
            return False

    def getRight(self):
        if self.isWrong(): return []
        fullWord = self._getFullWord()
        if 'Delta' not in fullWord:
            fullWord = re.sub(r'[\(\)a-z]', '', fullWord)
        wordList = re.compile(r'[_:\(\)\[\]+\-,\.]').split(fullWord)
        return self._wordToVarList(wordList)

class DnaNalaRule(Rules):

    @staticmethod
    def meets_condition(description:dict):
        varWord = description['varWord']
        return ('>' in varWord
        and bool(re.findall(r'to|with|by|for|ins|del|Non|Delta|of|at', varWord))
        )
    def _getFullWord(self):
        self.fullspan = self.span
        self.fullWord = self.varWord
    def getRight(self):
        self._getFullWord()
        nt_list = re.findall(r'[ACGT]', self.varWord)
        if len(nt_list) == 2:
            site = self._getVarSite()[0]
            return [site.join(nt_list)]
        else:
            print(nt_list)

class DnaRule(Rules):

    @staticmethod
    def meets_condition(description:dict):
        varWord = description['varWord']
        return ('>' in varWord
        and not re.findall(r'[lpus]', varWord))

    def getRight(self):
        if self.isWrong(): return []
        new_word_list = []
        if '/' in self.fullWord:
            wordList = self.fullWord.split('/')
            for word in wordList:
                new_word_list.extend([w for w in DNAWord(word).getSplitMutationWordList()])
        elif ',' in self.fullWord:
            fullLetter = self._getFullLetter(self.fullWord)
            fullSite = self._getFullSite(self.fullWord)
            if len(fullLetter) == 2 and len(fullSite) == 1:
                new_word_list.append(fullSite[0].join(fullLetter))
            else:
                wordList = self.fullWord.split(',')
                for word in wordList:
                    if MutationWord.fitPattern(word) and '>' not in word:
                        continue
                    elif '>' in word:
                        new_word_list.extend([w for w in DNAWord(word).getSplitMutationWordList()])
                    else:
                        print(word, self.fullWord)
        else:
            try:
                new_word_list = [w for w in DNAWord(self.fullWord).getSplitMutationWordList()]
            except:
                print(self.fullWord, self.varWord)
        return new_word_list


class AaRule(Rules):

    @staticmethod
    def meets_condition(description:dict):
        varWord = description['varWord']
        return ReplaceAminoAcids().hasAminoAcids(varWord)

    def getRight(self):
        fullWord = ReplaceAminoAcids().replaceAminoAcids(self._getFullWord())
        full_list = re.split(r'\.|:', fullWord)
        new_list = []
        for full in full_list:
            if MutationWord.fitPattern(full):
                for var in MutationWord(full).getSplitMutationWordList():
                    new_list.append(var)
            else:
                letters = self._getFullLetter(full)
                site = self._getFullSite(full)
                if len(letters) == 2 and len(site) == 1:
                    new_list.append(site[0].join(letters))
        return new_list


class PrototypeNucleotide(Rules):
    '''if virus.isHPV() and varWord[0] in list('ACEGT'): '''
    @staticmethod
    def meets_condition(description:dict):
        varWord = description['varWord']
        return (
            not re.findall(r'[\(\)]', varWord)
            and description.get('virus') == 'HPV'
            and varWord[0] in list('ACEGT')
            and varWord[-1] in list('ACGT')
            )
    Epattern = re.compile(r'E\d{2,5}[ACGT]')# E350G
    NTpattern = re.compile(r'[ACGT]\d{2,5}[ACGT]')# E-C109G
    keyWords = ['European', 'prototype']

    def _hasKeyWords(self):
        ls = [True if word in self.sentence else False for word in self.keyWords]
        return any(ls)

    def _hasPattern(self):
        if bool(self.Epattern.findall(self.varWord)):
            return True
        elif bool(self.NTpattern.findall(self.varWord)):
            return self._hasE_ahead()
        else: 
            return False

    def _hasE_ahead(self):
        lo = self._getStart()
        return True if self.sentence[lo-2:lo] == 'E-' else False

    def isWrong(self):
        return (self._hasKeyWords() and self._hasPattern())

    def getRight(self):
        self._getFullWord()
        return [] if self.isWrong() else [self.varWord]

class StandardRule(Rules):

    @staticmethod
    def meets_condition(description:dict):
        varWord = description['varWord']
        return (MutationWord.fitPattern(varWord) 
            and not bool(re.findall(r'STOP|stop|Stop|[\(\)/ ]|[pcg]\.', varWord))
            )

    def _isFull(self):
        return (self.varWord == self._getFullWord())

    def isWrong(self):
        if self._isFull():
            return False
        elif re.findall(r'[\[\]/\(\)Z\-+:;,\._\>]|&gt;|pc|rt', self._getFullWord()):
            return False
        else:
            return True

    def _fullWordToList(self, Word):
        fullWord = self._replaceAbbr(Word)
        wordList = re.compile(r'[_:\(\)\[\]+\-,\.?]|&gt;|[a-z]').split(fullWord)

        return self._wordToVarList(wordList)
    def getRight(self):
        if Rules.isWrong(self) or self.isWrong(): return []
        if self.virus == 'HIV' and 'SSS' in self.varWord:
            return [self.varWord]
        if self.varWord.isupper():
            if self._isFull():
                return self._split()
            elif self.varWord == self._getFullWord().strip('\*'):
                return self._split()
            else:
                return [self.varWord]
        else:
            return self._fullWordToList(self.fullWord)

class ParenthesisRule(Rules):
    pattern = re.compile(r'[ACGT]\(\d{1,5}\)[ACGT]')

    @staticmethod
    def meets_condition(description:dict):
        varWord = description['varWord']
        return bool(ParenthesisRule.pattern.findall(varWord))

    def isWrong(self):
        if len(self._getFullSite()) == 2 and len(self._getVarSite()) == 1:
            return True
        else: return False

    def getRight(self):
        self._getFullWord()
        if super().isWrong() or self.isWrong(): return []
        site = self._getVarSite()[0]
        letter = self._getVarLetter()
        return [site.join(letter)]

class StopEnd(Rules):
    '''if varWord[-1] == 'S':'''
    
    @staticmethod
    def meets_condition(description:dict):
        varWord = description.get('varWord')
        return (
            varWord[-1] == 'S' or
            bool(re.findall(r'STOP|stop|Stop',varWord))
        )

    def isWrong(self):
        end = self._getEnd()
        return bool(re.findall(r'STOP|stop|Stop', self._getFullWord()))

    def getRight(self):
        fullWord = self._getFullWord()
        if self.varWord == fullWord:
            return self._split()
        elif self.isWrong():
            self.varWord = re.sub(r'STOP|stop|Stop|S', 'X', self.varWord)
            return self._split()
        elif self.virus == 'EBV' and 'SPM' in fullWord:
            return [self.varWord.split('/')[0]]
        elif '/' in self.varWord:
            return MultiVar(description=self.description).getRight()
        else:
            self.varWord = re.sub(r'[a-z]', '', self.varWord)
            return self._split()

class Lower(Rules):

    @staticmethod
    def meets_condition(description:dict):
        varWord = description['varWord']
        return bool(re.findall(r'[ACGT]\d{1,5}[acgt]', varWord))

    def getRight(self):
        self._getFullWord()
        return [self.varWord.upper()]

class PubTatorFormat(Rules):

    @staticmethod
    def meets_condition(description:dict):
        return re.findall(r'[pc]\.', description['varWord'][:2])

    def isWrong(self):
        return False

    def _getFullWord(self):
        self.fullspan = self.span
        self.fullWord = self.varWord

    def getRight(self):
        self._getFullWord()
        var = PubTatorWord(self.varWord)
        return var.getIdentifier()

