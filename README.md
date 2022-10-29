# Code Repo for Cell Village Project for Xu  
Wilson wrote it on 10/29/2022  
  
### List of scripts used in the following orders:  
1. scRNAseq.sh  
Purpose: To run cellranger analysis on the FASTQ files. Only for one run.
How to run:  
```
sh scRNAseq.sh <samplename>  
  
Details:
<samplename> : The name of the sample used to identify the specific pairs of fastq files. See FASTQ File naming.

Within the shell script:
REF: reference genome built for cellrange (hg38).
FQ: the full path to the fastq file folder. Therefore, $S1 refers to the folder containing the fastq file.
```  
  
Output:
A folder with complete cellranger run.  

Note:
Possible to rerun in case the analysis was cutoff. cellranger is able to continue the run from where it stops.  


### List of R object and its meaning:  

