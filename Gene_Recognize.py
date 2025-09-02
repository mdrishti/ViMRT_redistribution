# *_* coding: utf-8 *_*
# @Time     : 10/21/2021 2:45 PM
# @Author   : zong hui
# @object   : some words here

import sys,os,getopt
import re
import json
import pandas as pd
from tqdm import tqdm
import string
os.chdir("./")

def main(argv):
    inputfile = './'
    outputfile= 'generesult.txt'
    virus= 'fullvirus'
    try:
        opts, args = getopt.getopt(argv,"i:o:v:",["ifile=","ofile=","virus="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            sys.exit()
        if opt in ("-o", "--opath"):
            outputfile = arg
        if opt in ("-i", "--ipath"):
            inputfile = arg      
        if opt in ("-v", "--virus"):
            virus = arg  
    return inputfile,outputfile,virus

class Gene_Recognize():

    def __init__(self,vocabulary:'pandas.core.frame.DataFrame',virus="fullvirus"):
        self.vocabulary = vocabulary.loc[vocabulary.gene == vocabulary.gene].copy()
        self.vocabulary.loc[:,"length"]=self.vocabulary.loc[:,"gene"].str.len()
        self.vocabulary=self.vocabulary.sort_values(by = ['length'], ascending =False)
        if virus!="fullvirus":
            virus=virus.split(";") 
            self.vocabulary = self.vocabulary[self.vocabulary["virus"].isin(virus)]
        else:
            self.vocabulary=self.vocabulary
    def str_insert(self,str_origin, pos, str_add):
        str_list = list(str_origin)    
        str_list.insert(pos, str_add)  
        str_out = ''.join(str_list)    
        return  str_out
    def quoteadd(self,quotecheck):
        ###处理引号,括号
        for quote in ["'",'"',"\(",'\)',"\[",'\]']:
            quotecheck=quotecheck.strip(quote)
            quotesiteadd=0
            for quotesite in re.finditer(quote,quotecheck, re.IGNORECASE):
                quotesite=int(list(quotesite.span())[0])+quotesiteadd
                quotecheck=self.str_insert(quotecheck,quotesite,"\\")
                quotesiteadd+=1   
        return  quotecheck
    def Match(self,description:dict,gene_vocabulary_error:"pandas.core.frame.DataFrame"):
        generesult={}
        single_wt = r'[AC-IK-NP-TV-YZ]'
        single_mt = r'(|[*]|STOP|stop|Stop|non-stop|St\.|[AC-IK-NP-TV-Y]{0,4}|\s)'
        site = r'[1-9]\d{0,5}'
        nublist=list(map(str,range(10)))
        nubletter=list(string.ascii_lowercase)
        for txt in description["text"]:
            lineinfo=txt.split("|pmid|")[0]
            print (txt)
            linetext=txt.split("|pmid|")[1]
            resultlast=[]
            sitefirst=set()
            errorsitesure=[]
            siteall=set()
            for key in set(self.vocabulary["virus"]):
                case=list(gene_vocabulary_error["gene"][(gene_vocabulary_error["error_type"]=="upper_lowercase") & (gene_vocabulary_error["virus"]==key)])
                gene_repeat=gene_vocabulary_error[(gene_vocabulary_error["error_type"]=="gene_repeat") & (gene_vocabulary_error["virus"]==key)]
                gene_repeat = gene_repeat.loc[gene_repeat.gene == gene_repeat.gene].copy()
                gene_repeat.loc[:,"gene"]=gene_repeat.loc[:,"gene"].str.upper()
                Special=gene_vocabulary_error[(gene_vocabulary_error["error_type"]=="same_word") & (gene_vocabulary_error["virus"]==key)]
                Special = Special.loc[Special.gene == Special.gene].copy()
                Special.loc[:,"gene"]=Special.loc[:,"gene"].str.upper()
                genetxt=self.vocabulary["gene"][self.vocabulary["virus"]==key]
                symbol=self.vocabulary["genesymbol"][self.vocabulary["virus"]==key]
                Patternraw=r'(vt){{0,1}}(wt){{0,1}}{3}(_){{0,1}}{0}(/{0}){{0,}}({1})({2})(/{2}){{0,}}(mutation){{0,1}}|{3}({1})({2})(/{2}){{0,}}|{3}({1})|(vt){{0,1}}(wt){{0,1}}{3}(s){{0,1}}'
                for regex in zip(genetxt,symbol):
                    regex=list(regex)
                    ###处理引号,括号
                    regex[0]=self.quoteadd(regex[0])
                    ####in的大小写
                    if regex[0].upper() in case:
                        Pattern = Patternraw.format(single_wt, site, single_mt,regex[0])
                        regexlast='\\b(' + Pattern + ')\\b'
                        for match in re.finditer(regexlast,linetext):
                            matchsite=list(match.span())
                            matchtext=match.group()
                            resultlast.append((key,matchtext,
                                str(matchsite[0]),str(matchsite[1]),regex[1]))
                    else:
                        ###有和其他单词重叠的
                        if regex[0].upper() in list(Special["gene"]):                            
                            Specialgeneall=list(Special["error_word"][Special["gene"]==regex[0].upper()])
                            for Specialgene in Specialgeneall:
                              Specialgene=self.quoteadd(Specialgene)
                           
                              #Patternerror = Patternraw.format(single_wt, site, single_mt,Specialgene)
                              regexerror='\\b' + Specialgene+'(s{0,})\\b'
                              
                              for error in re.finditer(regexerror,linetext, re.IGNORECASE):
                                    value=range(error.span()[0],error.span()[1]+1)
                                    errorsitesure.extend(list(value))
                  
                        ####处理"A -to- C"
                        if len(regex[0].upper())==1:
                            Patternto=r'({1}){{1}}(-){{0,1}}(\s){{0,1}}to{{1}}(-){{0,1}}(\s){{0,1}}({1}){{1}}'
                            Patternto = Patternto.format(regex[0],single_wt)
                            Patternto='\\b(' + Patternto + ')\\b'            
                            for error in re.finditer(Patternto,linetext, re.IGNORECASE): 
                                value=range(error.span()[0],error.span()[1])
                                errorsitesure.extend(list(value))                                          
                        Pattern = Patternraw.format(single_wt, site, single_mt,regex[0])
                        regexlast='\\b(' + Pattern + ')\\b'

                        for match in re.finditer(regexlast,linetext, re.IGNORECASE):
                            matchsite=list(match.span())
                            matchtext=match.group()
                            Patternmut=r'(vt){{0,1}}(wt){{0,1}}{0}(/{0}){{0,}}({1})({2})(/{2}){{0,}}(mutation){{0,1}}|{3}({1})({2})(/{2}){{1,}}|{3}{1}'
                            Patternmut = Patternmut.format(single_wt, site, single_mt,regex[0])
                            Patternmut='\\b(' + Patternmut + ')\\b'
                            number=matchsite[0]+len(regex[0])
                            
                            if matchtext.startswith("vt") and not regex[0].endswith("vt"):
                                matchtext=matchtext[2:(len(regex[0])+2)]
                                matchsite[1]=matchsite[0]+len(regex[0])+2
                                matchsite[0]=matchsite[0]+2
                            elif matchtext.startswith("wt") and not regex[0].endswith("wt"):
                                matchtext=matchtext[2:(len(regex[0])+2)]
                                matchsite[1]=matchsite[0]+len(regex[0])+2
                                matchsite[0]=matchsite[0]+2
                            ####处理gene_repeat
                            elif len(matchtext)>len(regex[0]) and regex[0].upper() in list(gene_repeat["gene"]):
                                matchtext=matchtext[:len(regex[0])]
                                matchsite[1]=matchsite[0]+len(regex[0])
                                rightword=list(gene_repeat["error_word"][gene_repeat["gene"]==matchtext.upper()])[0]                         
                                if matchtext==rightword:
                                  continue
                                else:
                                    matchtext=matchtext[:len(regex[0])]
                                    matchsite[1]=matchsite[0]+len(regex[0])                       
                            ###判断s是基因还是单复数
                            elif matchtext=="s":                                
                                if linetext[(matchsite[0]-1):(matchsite[0]+2)]=="(s)" and linetext[(matchsite[0]-2)] in nubletter:      
                                    continue
                                elif linetext[(matchsite[0]-1)]=="'" and not re.findall(r'\d+',matchtext, re.IGNORECASE):  
                                    continue
                                else:
                                    matchtext=matchtext[:len(regex[0])]
                                    matchsite[1]=matchsite[0]+len(regex[0])
                            ###处理s匹配上的
                            elif matchtext.endswith("s") and not regex[0].endswith("s") and matchtext!="s":
                                if len(matchtext)<len(regex[0]):
                                  matchsite[1]=matchsite[0]+len(matchtext)
                                else: 
                                  matchsite[1]=matchsite[0]+len(regex[0])+1
                            ###处理突变中的错误基因 rtA181S 
                            elif len(regex[0])==1 and len(matchtext)!=1 and re.findall(Patternmut,matchtext, re.IGNORECASE) and linetext[matchsite[0]] in [regex[0].upper()]: 
                                continue
 
                            ###处理突变中的错误基因 rtA181V/T/S 
                            elif len(regex[0])==1 and len(matchtext)==1:
                                Patternmut1=r'{0}{{0,}}(/{0}){{0,}}({1}){{0,}}({2}){{0,}}(/{2}){{0,}}'
                                Patternmut1 = Patternmut1.format(single_wt, site, single_wt)
                                Patternmut1='\\b(' + Patternmut1 + ')\\b'
                                Patternmut2=r'{0}{{1,}}(/{0}){{0,}}({1})'
                                Patternmut2 = Patternmut2.format(single_wt, site)
                                Patternmut2='\\b(' + Patternmut2 + ')\\b'
                                nogene=[]
                                ###判断基因P<0.009
                                if regex[0].upper()=="P" or regex[0].upper()=="N":
                                    if matchsite[1]<=len(linetext) and linetext[matchsite[1]] in [">","<","="]:           
                                      nogene.append("N")
                                    elif matchsite[1]+1<=len(linetext) and linetext[matchsite[1]:(matchsite[1]+2)] in [" >"," <"," =",".1"," 0"]:           
                                      nogene.append("N")
                                
                                ###判断x表示乘号的基因  
                                if regex[0].upper()=="X":
                                   if linetext[matchsite[0]-2] in nublist and linetext[matchsite[1]+1] in nublist: 
                                    continue
                                for errormut in re.finditer(Patternmut1, linetext, re.IGNORECASE):
                                    if len(errormut.group())!=1 and matchsite[0] in list(range(errormut.span()[0],errormut.span()[1]+1)):
                                        errortxt=errormut.group() 
                                        if re.findall(Patternmut2, errortxt, re.IGNORECASE): 
                                            matchtext=matchtext[:len(regex[0])]
                                            matchsite[1]=matchsite[0]+len(regex[0])   
                                        else: 
                                            nogene.append("N")
                                if "N" in nogene:    
                                    continue
                            
                            elif len(matchtext)<len(regex[0]):
                                matchtext=matchtext
                                matchsite[1]=matchsite[0]+len(matchtext)
                            ###判断基因最后是不是数字
                            elif regex[0][-1] in nublist and linetext[number] in nublist: 
                                continue
                            ###判断X
                            elif regex[0].upper()=="X":
                              if linetext[matchsite[0]-2] in nublist and linetext[matchsite[1]+1] in nublist: 
                                continue
                              else:
                                matchtext=matchtext[:len(regex[0])]
                                matchsite[1]=matchsite[0]+len(regex[0])
                            ###判断单个字母匹配到突变如T204V/S
                            elif len(matchtext)>len(regex[0]) and len(regex[0])==1 and matchtext[0] not in nubletter:
                                continue 
                                                                ######判断p53 nogene
                            elif regex[0].upper()=="P" and re.findall(r'\bp[1-9]\d{0,2}\b',matchtext, re.IGNORECASE):   
                                continue     
                            else:
                                matchtext=matchtext[:len(regex[0])]
                                matchsite[1]=matchsite[0]+len(regex[0])                            
                            if len (set(matchsite) & siteall)==0 and len (set(errorsitesure) & set(matchsite))==0:                               
                                siteall=siteall | set(range(matchsite[0],matchsite[1]))                                
                                resultlast.append((key,matchtext,
                                str(matchsite[0]),str(matchsite[1]),regex[1]))
            generesult[txt]=resultlast
        return(generesult)

if __name__ == "__main__":
    file_name,out_file,virus=main(sys.argv[1:])
    error_file="./genecorpus/gene_vocabulary_error.txt"
    gene_vocabularyname="./genecorpus/gene_vocabulary.txt"
    infile= pd.read_csv(file_name,sep="\t",header=None,index_col=False)
    vocabulary= pd.read_csv(gene_vocabularyname,sep="\t",header=0,keep_default_na=False,index_col=False,encoding='gbk')
    gene_vocabulary_error= pd.read_csv(error_file,sep="\t",header=0,index_col=False,encoding='gbk')
    description={}
    description['text']=infile[0]
    infile= pd.read_csv(file_name,sep="\t",header=None,index_col=False)
    t = Gene_Recognize(vocabulary,virus)
    resultgene=t.Match(description,gene_vocabulary_error)
    with open(out_file,"a") as outf:
        outf.write("\t".join(("pmid","sentence","virus","abstract_text","start_site","end_site","standardized_name"))+"\n") 
        for sentence in resultgene:
            for gene in resultgene[sentence]: 
                outf.write(sentence.split("|pmid|")[0]+"\t"+sentence.split("|pmid|")[1]+"\t")           
                outf.write("\t".join(gene)+"\n")    


