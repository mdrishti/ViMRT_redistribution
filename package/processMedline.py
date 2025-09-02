import os
import re
class Medline:
    pattern = re.compile(r'\b(mutation|mutant|variant|variat|deletion|insertion|transvertion|transition|heterozygous|homozygote|homozygous|substitut|SNP)', re.IGNORECASE)
    ROOT = 'Virus'
    def __init__(self, paperdir):
        self.root = self.ROOT+'/'+ paperdir.split('/')[1]
        self.passage = {}
        self.pmid = paperdir.rsplit('/')[-1].strip('.txt')
        with open(paperdir) as f:
            self.lines = f.readlines()
    def process(self):
        self._getPassage()
        self._getTerm()
    def _getPassage(self):
        key, OT = '', []
        for line in self.lines:
            line = line.strip('\n')
            if line[:2] == 'TI' or line[:2] == 'AB':
                key = line[0].lower()
                self.passage[key] = line[6:]
            elif line[:6] == ' '*6 and key != '':
                self.passage[key] += line[5:]
            elif line[:3] == 'OT ':
                OT.append(line[6:])
            else:
                key = ''
        if OT != []:
            self.passage['o'] = ';'.join(OT)
    def _getTerm(self):
        assert 't' in self.passage.keys() and 'a' in self.passage.keys(), "{}{}".format(self.lines, self.pmid)
        outdir = "{0}/{1}/{2}.txt".format(self.root, 'with_mutation', self.pmid)
        if self.pattern.findall(self.passage['a']):
            with open(outdir, 'w') as f:
                for k,v in self.passage.items():
                    f.write(k+'|'+self.pmid+'|'+v+'\n')
        else:
            if os.path.exists(outdir):
                os.remove(outdir)


class WithMutation:
    ROOT = 'Virus'
    def __init__(self, root):
        self.root = self.ROOT+'/'+ root
        self.medline = self.root+'/all_abstract/'
        self.pmidlist = self.root+'/pmidlist.txt'
    def _processMedline(self):
        for paperFile in os.listdir(self.medline):
            try:
                Medline(self.medline + paperFile).process()
            except AssertionError as e:
                continue
    def getList(self):
        self._processMedline()
        filelist = [name.strip('.txt') for name in os.listdir(self.root+'/with_mutation')]
        with open(self.pmidlist, 'w') as f:
            f.write('\n'.join(filelist))

if __name__ == "__main__":
    WithMutation('IV').getList()