# VO-SemRep Installation

Installation details to be provided here ....

- Prepared by: Nor Danish Imran
- Helped by Oliver He, Jie Zheng. 

## Outline:
- General SemRep installation. 
   - https://github.com/lhncbc/SemRep 
   - Give some links and summarize the important points. 
- VO customization for SemRep. 
   - Introduce Background with EFO customization (https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/additional-tools/DataFileBuilder.html), and why/how it does not work for VO.  
   - How to prepare VO-related files for SemRep to use. 
   - How to use these files.
   - Possible bugs and how to solve them. 

## General SemRep installation:
The general steps of installing SemRep can be found in the SemRep repo [here](https://github.com/lhncbc/SemRep). The SemRep linked provide installation steps for 2006AA and 2018AA of UMLS release. If you wish to use newer release of UMLS, navigate below.

### Using SemRep with updated UMLS Metathesaurus:
Because SemRep hasn't been updated in a while, we need to do some workaround for SemRep to work with newer version of UMLS Metathesaurus. 
   1. Make sure you have installed the UMLS release you wanted to use it with.  
   2. Navigate to <SemRep dir>/DATA/UMLS_HIERARCHY folder, and make a copy of each hierarchy file from 2006AA and 2018AA. Then for each copy, rename it to the {YEAR}{VERSION} of your chosen UMLS release, as shown below.   
```
hrel-UMLS_btree.2006AA_1
hrel-UMLS_btree.2006AA_2
...
hrel-UMLS_btree.2006AA_5
hrel-UMLS_btree.2018AA_1

// example: 2022AB
// rename the files to

hrel-UMLS_btree.2022AB_1
hrel-UMLS_btree.2022AB_2
...
hrel-UMLS_btree.2022AB_5
hrel-UMLS_btree.2022AB_6

// You should have 6 new hierarchy files.
```
   3. Add {YEAR} to "AVAILABLE DATA" variable in <SemRep dir>/bin/install.sh file.
   4. Follow general steps of installation from SemRep repo. 

## VO customization for SemRep. 

### Introduce Background with EFO customization
Introduction: 
Link: https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/additional-tools/DataFileBuilder.html
However, it does not work for VO.  
Possible reasons: 
- The EFO instruction was prepared a long time ago, and ...
- other reasons ... 


### Preparation and Usage of VO-related files for SemRep to use. 

More clear details 


### Possible bugs and how to solve them. 

