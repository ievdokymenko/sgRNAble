#!/usr/bin/env python3

import numpy as np
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
from Cas9_Calculator import *
from File_GUI import Get_Sequence

#Length of the Guide RNA desired
Guide_RNA_length = 20

#Whether the CRISPR program is being used for CRIi (1) or CRIa (n) or normal/wt (0)
Usage = 1 # length of Amplification range user would like

#Combine the Coding and Template Strands into a single strand
def CombinetoStr (Template_Guides, Coding_Guides):
  Guides = []

  for i in range (len(Template_Guides)):
    if (i < len(Template_Guides)):
      Guides.append(Template_Guides[i])

  for i in range (len(Coding_Guides)):
    if (i < len(Coding_Guides)):
      Guides.append(Coding_Guides[i])

  return Guides

def ReverseComplement(nucleotide_sequence):
  comp = []
  for c in nucleotide_sequence:
    if c == 'A' or c == 'a':
      comp.append('T')
    if c == 'G' or c == 'g':
      comp.append('C')
    if c == 'U' or c == 'u' or c == 'T' or c == 't':
      comp.append('A')
    if c == 'C' or c == 'c':
      comp.append('G')
  rev_comp = ''.join(reversed(comp))
  return rev_comp

Target_Seqs, Genome = Get_Sequence()

Genome = Genome + Genome.reverse_complement()
SeqIO.write(Genome, "Total_Genome_Plus_RC", "fasta")

#Obtain the Guide RNAs from the Target Sequence
from Azimuth_Finder import * #placed here to avoid multithreading issues with Tkinter (I think)
PAM = "GG"
Guide_Num_Cutoff = 20 # Cutoff per strand

Target_Guides, Position_List, Direction_List = Azimuth_Guides(Target_Seqs, PAM,Guide_Num_Cutoff, Usage )

#Combine the information from the Guide RNAs into one single array.
Guide_Info = np.vstack((Position_List, Direction_List)).T

#Send Data to the model
Cas9Calculator=clCas9Calculator(['Total_Genome_Plus_RC'])
sgRNA1 = sgRNA(Target_Guides,Guide_Info, Cas9Calculator)
sgRNA1.run()
