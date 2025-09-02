## Purpose of this GitHub repo

This is a redistribution of ViMRT software. The original implementation provided by the authors is here: http://bmtongji.cn:1225/mutation/download/?file_name=ViMRT.zip . This git repo is for playing around with code and repurpose for personal projects only. 

The readme contents below is from the original readme file provided by the authors.

---------------------------------------------------------------------------------------------------------

######   Using ViMRT
ViMRT is a text-mining tool and search engine for automated virus mutation recognition by rule patterns and regular expression patterns for different written forms of virus mutation based on natural language processing. It can also quickly and accurately search virus mutation-related information including virus genes and related diseases.

#### Getting Started with ViMRT
System Python
1.Installing Python


```bash
conda create --name py3.8 python=3.8
conda activate py3.8
```


2.Installing python package

```
pip3 install -r requirements.txt
```

####  Virus mutation recognition
1.Downloading BioC_xml format

```
python Bio_download.py -i [input] -o [output] -s [sources]
```
input: user can provide the input file, such as PMIDlist.txt (for PMIDlist type). 
output: user can provide the output folder path. 
sources: user can choose the output file source: PubMed | PMC | PubMed_PMC, which will obtain the PubMed abstracts, PMC full text literatures or both, respectively)

!Example: python Bio_download.py -i PMIDlist.txt -o ./tmVar/tmvar_input -s PubMed

2.Optimizing results of tmVar by rule patterns
2.1 Identifying mutation by tmVar
In this step, user firstly needs to download tmVar(https://www.ncbi.nlm.nih.gov/research/bionlp/Tools/tmVar/), User can also downlaod the concise version min-tmVar.zip(http://bmtongji.cn:1225/mutation/download/?file_name=min-tmVar.zip). tmVar can only run in window environment or linux environment.
```
java -Xmx5G -Xms5G -jar tmVar.jar [input] [output]
```
input: user can provide the input folder path.
output: user can provide the output folder path. 

2.2.Optimizing results of tmVar by rules patterns

```
python ViMRT.py -i [input] -o [output] -v [virus] -f [formart] -m rules
```
input: input folder path includes output result files of tmVar.
output: user can provide the output folder path. The default path is the current path.
virus: user can provide one virus name. The default parameter is 'Unknown'.
formart: user can choose one input file formart: PubTator or BioC. The default parameter is 'BioC'.


Example: python ViMRT.py -i ./tmVar_result/ -o ./ViMRT_result/ -v HBV -f BioC -m rules



3.Identifying mutation by regular expression patterns
```
python ViMRT.py -i [input] -o [output] -v [virus] -f [formart] -m regex
```
input: input folder path includes output result files of tmVar or BioC format files by downloading.
output: user can provide the output folder path. The default path is the current path.
virus: user can provide one virus name. The default parameter is 'Unknown'.
formart: user can choose one input file format: PubTator or BioC. The default parameter is 'BioC'.

Example: python ViMRT.py -i ./tmVar_result/ -o ./ViMRT_result/ -v HBV -f BioC -m regex


4.Integrating the result of rules and regex

```
python ViMRT.py -i [input] -o [output] -f [formart] -c concat
```


input: input path includes identification result files by both rule and regex.
output: user can provide the output folder path. The default path is the current path.
formart: user can choose one input file format: PubTator or BioC. The default parameter is 'BioC'.


or

Identifying directly the mutation by rules and regex
```
python ViMRT.py -i ./tmVar_result/ -o ./ViMRT_mutation/ -v HBV -f BioC
python ViMRT.py -i ./tmVar_result/ -o ./ViMRT_mutation/ -v HBV -f PubTator
```


####   Virus gene recognition
```
python Gene_Recognize.py  -i [input] -o [output] -v [virus]
```

input: user can provide input file (gene_example.txt). The input file format: PMID+"|pmid|"+sentence
output: user can provide output folder path. The default path is the current path.
virus: user can choose virus name, such as HBV, HBV;HPV, etc. The default parameter is "fullvirus"


Example: python Gene_Recognize.py -i ./gene_example.txt -o ./ViMRT_gene/ -v HBV


Note: we recommend selecting one virus name. Becuase the default parameter will match genes of all viruses in turn, which will run for a long time.

####   Disease recognition

1.Installing stanza
```
pip install stanza
```
2.Download stanza disease models
```
import stanza 
stanza.download('en', package='mimic', processors={'ner': 'bc5cdr'}, verbose=False)
stanza.download('en', package='mimic', processors={'ner': 'ncbi_disease'}, verbose=False)
```
3.Identifying and optimizing disease
```
python Disease_Recognize.py -i [input] -o [output]
```

input: user can provide input file (disease_example.txt). The input file format: PMID+"|pmid|"+sentence.
output: user can provide output folder path. The default path is the current path.

Example: python Disease_Recognize.py -i ./disease_example.txt -o ./ViMRT_disease/ -v HBV



