# sgRNAble

Brief Outline

## Installation

### Prerequisites

What things you need to install the software and how to install them

* Python3
* Environment Manager (Anaconda is used)


### Installation Guide

Prior to installation,it is best practise to create a new enviroment to store the program and dependencies locally. This setup will create an conda environment with the name sgRNAble and install all required dependencies. 

```
conda create --name sgRNAble python=3.7
conda activate sgRNAble
pip install sgRNAble
conda deactivate
```

In the future, the program can be run by activating the python env and running the program.

```
conda activate sgRNAble
sgRNAble -t TARGET_FILE -g GENOME_FILE
conda deactivate
```

## Quick Run Guide

Ensure that you have a file containing the gene of interest (Target Sequence), the genome of the organism (Genome), and
any additional DNA present. The gene of interest must be present in the genome or the other additional DNA added to the script.

navigate to the folder in terminal and type in

```
python optimal_guide_finder/guide_finder.py -t data/Fasta_Files/GFP.fasta -g data/Fasta_Files/E_coli_MG1655_genome.fasta data/Fasta_Files/GFP.fasta
```

## Authors
Siddarth Raghuvanshi, Ahmed Abdelmoneim, and Avery Noonan

## Contact

Need something? Send me an email at Raghuvanshi.Siddarth@gmail.com
