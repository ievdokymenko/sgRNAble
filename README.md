# CRISPR_Guide_RNA

# About
This program is a combination of the Azimuth machine learning software developed by Microsoft https://www.microsoft.com/en-us/research/project/crispr/ and a thermodynamic model constructed by the Salis Lab (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004724). 

Uses the Azimuth model to determine the strongest RNA binding sequences, and passes the sequences to the model developed by Farasat/Salis to run a script that finds the best guide RNAs for a given gene (Target Sequence), Genome and any additional DNA sequences. The model is based on the 
information based on CAS9 and would be different for other CRISPR systems. 

The script is written using python v3 and cannot run on python v2. 

## Setup

To run the script, ensure you have Anaconda installed
  
Place the Azimuth Folder in the Anaconda Folder in /lib/python3.6/site-packages/ 

If Azimuth is previously installed, replace the folder (Converting the program from Python v2 to v3)

## How to use
Ensure that you have a Fasta or Genebank file containing the gene of interest (Target Sequence), the genome of the organism (Genome), and any additional DNA present(Plasmids/Fosmids). The gene of interest must be present in the genome or the other additional DNA added to the script.

navigate to the folder in terminal and type in

```
python Optimal_Guide.py
```

Follow the prompts on screen. 

  
