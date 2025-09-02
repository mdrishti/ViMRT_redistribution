# *_* coding: utf-8 *_*
import stanza
import pandas as pd
import string
import math
import difflib
import sys,os,getopt
import re
import json
os.chdir("./")
def main(argv):
    inputfile = './'
    outputfile= 'disease_stanza.txt'
    try:
        opts, args = getopt.getopt(argv,"i:o:v:",["ifile=","ofile="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            sys.exit()
        if opt in ("-o", "--opath"):
            outputfile = arg
        if opt in ("-i", "--ipath"):
            inputfile = arg  
    return inputfile,outputfile

class Disease_Recognize():
    def __init__(self,entity="ncbi_disease;bc5cdr"):
        self.disentity=entity.split(";") 
    def Match(self,description:dict):
        resultdisease={}
        disease_bc5cdr_model = stanza.Pipeline('en', package='mimic', processors={'ner': "bc5cdr"}, verbose=False)
        disease_ncbi_model = stanza.Pipeline('en', package='mimic', processors={'ner': "ncbi_disease"}, verbose=False)
        for txt in description["text"]:
            lineinfo=txt.split("|pmid|")[0]
            linetext=txt.split("|pmid|")[1]
            resultlast=[]
            if  "bc5cdr" in self.disentity:
                disease_docs_bc5cdr = disease_bc5cdr_model(linetext)
                for ent in disease_docs_bc5cdr.entities:
                    if  ent.type!="CHEMICAL":
                        resultlast.append(("bc5cdr",ent.text,ent.type,str(ent.start_char),str(ent.end_char)))                
            if  "ncbi_disease" in self.disentity:
                disease_ncbi_bc5cdr = disease_ncbi_model(linetext)
                for ent in disease_ncbi_bc5cdr.entities:
                    resultlast.append(("ncbi_disease",ent.text,ent.type,str(ent.start_char),str(ent.end_char)))
            resultdisease[txt]=resultlast 
        return(resultdisease)

class Deal_Disease_Recognize():
    def __init__(self):
       pass
    def str_insert(self,str_origin, pos, str_add):
        str_list = list(str_origin)    
        str_list.insert(pos, str_add)  
        str_out = ''.join(str_list)    
        return  str_out
    def quoteadd(self,quotecheck):
        ###处理引号,括号
        for quote in ["'",'"',"\(",'\)',"\[",'\]','\/']:
            quotecheck=quotecheck.strip(quote)
            quotesiteadd=0
            for quotesite in re.finditer(quote,quotecheck, re.IGNORECASE):
                quotesite=int(list(quotesite.span())[0])+quotesiteadd
                quotecheck=self.str_insert(quotecheck,quotesite,"\\")
                quotesiteadd+=1   
        return  quotecheck

    def Deal(self,output:"pandas.core.frame.DataFrame",vocabulary_error:"pandas.core.frame.DataFrame",CTD_disease:"pandas.core.frame.DataFrame",vocabulary_disease:"pandas.core.frame.DataFrame"):
        deletion=vocabulary_error[vocabulary_error["error_type"]=="deletion"]
        deletion = deletion.loc[deletion.entity == deletion.entity].copy()
        deletion.loc[:,"entity"]=deletion.loc[:,"entity"].str.upper()
        finalword=vocabulary_error[vocabulary_error["error_type"]=="finalword"]
        CTD_disease=CTD_disease.loc[:,["DiseaseName","DiseaseID","Synonyms"]]
        CTD_disease['Synonyms']=CTD_disease['Synonyms'].map(lambda x:str(x).split('|'))
        CTD_Synonyms=CTD_disease.explode('Synonyms')
        CTD_disease["Synonyms"]=CTD_disease["DiseaseName"]
        CTD_Synonyms=CTD_Synonyms[CTD_Synonyms['Synonyms']!="nan"].loc[:,["Synonyms","DiseaseID","DiseaseName"]]
        CTD_disease=pd.concat([CTD_disease,CTD_Synonyms],ignore_index=True)
        vocabulary_disease1=vocabulary_disease.loc[:,["Standardized_Disease","Standardized_Disease"]]
        vocabulary_disease1.columns=["disease","Standardized_Disease"]
        vocabulary_disease=pd.concat([vocabulary_disease,vocabulary_disease1],ignore_index=True)
        vocabulary_disease.loc[:,"length"]=vocabulary_disease.loc[:,"disease"].str.len()
        vocabulary_disease=vocabulary_disease.sort_values(by = ['length'], ascending =False)
        vocabulary_error = vocabulary_error.loc[vocabulary_error.entity == vocabulary_error.entity].copy()
        vocabulary_error.loc[:,"entity"]=vocabulary_error.loc[:,"entity"].str.upper()
        delloc=output.index
        nrow=list(range(output.shape[0]))
        dellist=[]
        output["difflib_value"]=""
        output["Standardized_Disease"]=""
        for i in nrow:
            key=list(output.loc[i])
            lineword=str(key[3])
            sentence=key[1]
            ###Delete some obvious errors
            if lineword.upper() in list(deletion["entity"]) or lineword.upper() in ['{}S'.format(a) for a in list(deletion["entity"])] or "R"+lineword.upper() in list(deletion["entity"]) or lineword.upper()+"R" in list(deletion["entity"]):
                dellist.append(i)
            for word in finalword["entity"]:
                if lineword.endswith(word) or lineword.endswith(word+"s") or lineword.endswith(word+"es"):
                    dellist.append(i)
                    
            maxv=[]
            ###Delete some low diff
            if str(lineword)!="nan" and i not in dellist:
                if lineword.upper() not in list(vocabulary_error["entity"]) and lineword.upper()[:-1] not in list(vocabulary_error["entity"]):
                    diff=difflib.get_close_matches(lineword,list(CTD_disease["Synonyms"]),1, cutoff=0.1)
                    if diff!=[]:                    
                        diffword=CTD_disease["Synonyms"][CTD_disease["Synonyms"]==diff[0]]
                        diffSeq=difflib.SequenceMatcher(None,lineword, diff[0]).quick_ratio()
                        output.loc[[i],["difflib_value"]]=diffSeq
                        output.loc[[i],["Standardized_Disease"]]=list(CTD_disease["DiseaseName"][CTD_disease["Synonyms"]==diff[0]])[0]
                        if diffSeq<0.7:
                            output.loc[[i],["difflib_value"]]=diffSeq
                            dellist.append(i)
            ###Corpus recognition
            copy=[]
            for word_id in range(vocabulary_disease.shape[0]):
                word=vocabulary_disease["disease"][word_id]
                word=self.quoteadd(word)
                regexerror='\\b' + word+'(s{0,})\\b'
                for error in re.finditer(regexerror,sentence, re.IGNORECASE):
                  if error.span()[0] not in copy and error.span()[0]+len(error.group()) not in copy:
                    keynew=[key[0],key[1],"match",error.group(),'DISEASE',error.span()[0],error.span()[0]+len(error.group()),"",vocabulary_disease["Standardized_Disease"][word_id]]
                    output.loc[output.shape[0]+1]=keynew
                    copy.extend(list(range(error.span()[0],error.span()[0]+len(error.group())+1)))
                    
        stanzadel=output.loc[dellist]
        output=output.drop(dellist)
        return output,stanzadel
if __name__ == "__main__":
    file_name,out_file=main(sys.argv[1:])
    CTD_name="./diseasecorpus/CTD_diseases.txt"
    vocabulary_error_name="./diseasecorpus/disease_vocabulary_error.txt"
    vocabulary_disease_name="./diseasecorpus/disease_vocabulary.txt"
    description={}
    infile= pd.read_csv(file_name,sep="\t",header=None,index_col=False)
    description['text']=infile[0]
    t = Disease_Recognize("ncbi_disease;bc5cdr")
    resultdisease=t.Match(description)
    with open(out_file,"a") as outf:
        outf.write("\t".join(["PMID","Sentence","Models","Disease","Type","Start_site","End_site","\n"]))
        for sentence in resultdisease:
            if  resultdisease[sentence]==[]:
                outf.write(sentence.split("|pmid|")[0]+"\t"+sentence.split("|pmid|")[1]+"\n")
            else:
                for disease in resultdisease[sentence]: 
                    outf.write(sentence.split("|pmid|")[0]+"\t"+sentence.split("|pmid|")[1]+"\t")      
                    outf.write("\t".join(disease)+"\n") 
    vocabulary_error= pd.read_csv(vocabulary_error_name,sep="\t",header=0,index_col=False)
    vocabulary_disease=pd.read_csv(vocabulary_disease_name,sep="\t",header=0,index_col=False)
    CTD_disease= pd.read_csv(CTD_name,sep="\t",header=0,index_col=False)
    t = Deal_Disease_Recognize()
    input_file=pd.read_csv(out_file,sep="\t",header=0,index_col=False)
    input_file=input_file.loc[:, ~input_file.columns.str.contains('^Unnamed')]
    diseaseresult,diseaseresultdel=t.Deal(input_file,vocabulary_error,CTD_disease,vocabulary_disease)
    diseaseresult=diseaseresult.loc[:,["PMID","Sentence","Disease","Start_site","End_site","Standardized_Disease","difflib_value"]]
    diseaseresultdel=diseaseresultdel.loc[:,["PMID","Sentence","Disease","Start_site","End_site","Standardized_Disease","difflib_value"]]
    diseaseresult=diseaseresult.drop_duplicates()
    out_filelast=out_file.replace(".txt","_")+"optimize.txt"
    stanzadel=out_file.replace(".txt","_")+"del.txt"
    diseaseresult.to_csv(out_filelast,sep="\t",index=False)
    diseaseresultdel.to_csv(stanzadel,sep="\t",index=False)