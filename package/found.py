import os
import pandas as pd
from superior.package.sub4.paper import PaperFile

class Found:
    
    def __init__(self, root:str,outpath:str,Virus="Unknown",form="BioC", method="regex_rules"):
        self.virus = Virus
        self.form = form
        self.outpath = outpath+'/'+form
        self.root = root
        self.outpathroot=outpath
        if method == "regex":
            self.total = {'regex': pd.DataFrame()}
        elif method == "rules":
            self.total = {'rules': pd.DataFrame()}                
        else:
            self.total = {'regex': pd.DataFrame(),'rules': pd.DataFrame()}
            
    def _sum(self):

        if os.path.isfile(self.root):
            self._find(self.root)
            return
        i = 0
        for fi in os.listdir(self.root):
            filename = self.root + '/'+ fi
            if fi.endswith(".xml") or fi.endswith(".XML") or fi.endswith(".PubTator"):
                filename = self.root + '/'+ fi
                print(i,"filename: ", filename)
                try:
                    self._find(filename)
                except Exception as e:
                    self.save(outfile='temp')
            else:
                print (self.root+ '/'+ fi +" - It's not either BioC or PubTor format.")
                break
            i += 1
    def _find(self, filename):
        if os.path.getsize(filename) == 0:
            os.remove(filename)
            return
        self.paper = PaperFile(filename, self.virus, self.form).processPaper()
        self.paper.all()
        for key, value in self.total.items():
            df = self.paper.total[key].drop_duplicates() #  ['pmid', 'sentence', 'processed']
            self.total[key] = pd.concat([value, df])

    def save(self, outfile:str=''):
        def region(data):
            region_dict = {
            't': "Title",
            'title': "Title",
            "a": "Abstract",
            "abstract": "Abstract",
            'i': "Introduction",
            'm': "Method",
            "r": "Result",
            "c": "Conclusion",
            'd': "Discussion",
            'f': "Figure",
            'TABLE': "Table",
            }  
            for key in region_dict.keys():
                data.loc[data["region"]==key,"region"]=region_dict[key]
            return data
        if outfile != '':
            for key, value in self.total.items():
                outdir = f'{outfile}_{key}.csv'
                return
        self._sum()
        for key, value in self.total.items():
            outdir = f'{self.outpath}/{self.form}_{key}.csv'
            currdir = self.outpath
            if not os.path.exists(self.outpathroot):
                os.mkdir(self.outpathroot)
            if not os.path.exists(currdir):
                os.mkdir(currdir)
            if value.shape[0]==0:
                print ("No virus mutations are identified")
            else:   
                value = value.sort_values(['pmid', 'sen_offset', 'start_in_sen'])
                value=region(value)
                value.rename(columns={"processed":'Standardized_mutation',"original":'Original_mutation',"start_in_sen":"start_in_sentence","end_in_sen":"end_in_sentence"},inplace=1)
                value=value[["Standardized_mutation","Original_mutation","pmid","region","sentence","start_in_sentence","end_in_sentence","virus"]]
                value.to_csv(outdir, index=0)
    
    def output(self, mut=""):
        self._sum()
        df_rules = self.total['rules']
        df_regex = self.total['regex']
        df_rules['source'] = 'rules'
        df_regex['source'] = 'regex'
        df_rr = pd.concat([df_rules, df_regex]).drop_duplicates(['pmid', 'processed', 'sentence', 'start_in_sen', 'end_in_sen'], keep= 'first')
        df_rr = df_rr.sort_values(['pmid', 'sen_offset', 'start_in_sen', 'source'], ascending=[True, True, True, False])



class Result:
    def __init__(self, root:str,outpath:str,form="BioC",concat="no"):
        self.form = form
        if concat=="no":
            self.root = root+'/'+form
            self.output = outpath+ '/'+form
        else:
            self.root = root
            self.output = outpath
        if not os.path.exists(self.output):
            os.mkdir(self.output)

    def _rulesAppendRegex(self):
        def share(series):
            if series['source'] == "rules" and series['shared'] == "right":
                return "rules_regex"
            elif series['source'] == "rules" and series['shared'] == "left":
                return "rules"
            elif series['source'] == "regex" and series['shared'] == "left":
                return "rules_regex"
            elif series['source'] == "regex" and series['shared'] == "right":
                return "regex"
            return 
        def region(data):
            region_dict = {
            't': "Title",
            'title': "Title",
            "a": "Abstract",
            "abstract": "Abstract",
            'i': "Introduction",
            'm': "Method",
            "r": "Result",
            "c": "Conclusion",
            'd': "Discussion",
            'f': "Figure",
            'TABLE': "Table",
            }  
            for key in region_dict.keys():
                data.loc[data["region"]==key,"region"]=region_dict[key]
            return data
        main_key = ['pmid', 'Standardized_mutation', 'sentence', 'start_in_sentence', 'end_in_sentence']
        if os.path.exists(self.root+f'/{self.form}_rules.csv') or os.path.exists(self.root+f'/{self.form}_regex.csv'):  
            df_rules = pd.read_csv(self.root+f'/{self.form}_rules.csv')    
            df_regex = pd.read_csv(self.root+f'/{self.form}_regex.csv')
            df_rules['source'] = 'rules'
            df_regex['source'] = 'regex'
            df_rules = df_rules.sort_values(['pmid', 'start_in_sentence']).drop_duplicates(main_key)
            df_regex = df_regex.sort_values(['pmid', 'start_in_sentence']).drop_duplicates(main_key)
            df_rules=region(df_rules)
            df_regex=region(df_regex)
            outer = pd.concat([df_rules, df_regex]).drop_duplicates(main_key, keep=False)
            right = pd.concat([outer, df_rules]).drop_duplicates(main_key, keep=False)
            left = pd.concat([outer, df_regex]).drop_duplicates(main_key, keep=False)
            right['shared'] = 'right'
            left['shared'] = 'left'
            df_total = pd.concat([right, left])
            df_total['found_by'] = df_total.apply(share, axis=1)
            self.df_total = df_total.drop_duplicates(main_key)
            self.df_total=self.df_total[["Standardized_mutation","Original_mutation","pmid","region","sentence","start_in_sentence","end_in_sentence","virus","found_by"]]
            self.df_total= self.df_total.sort_values(['pmid', 'start_in_sentence'])
            self.df_total.to_excel(self.output+f'/{self.form}_rules_regex.xlsx', index=0)


if __name__ == "__main__":
    found.save()