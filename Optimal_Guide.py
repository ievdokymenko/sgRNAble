#!/usr/bin/env python3

import numpy as np
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
from Cas9_Calculator import *
from File_GUI import Get_Sequence


#Length of the Guide RNA desired
Guide_RNA_length = 20

#Find the Guide RNAs in a Sequence
def PAM_Finder(Sequence, PAM, Direction):
 Guide_RNAs = []
 Location = []
 Strand = []

 Position = 0
 Temp_Sequence = Sequence
 if(Direction < 0):
   PAM = str(Seq(PAM).reverse_complement())
   print(PAM)
 while True:
   i = Temp_Sequence.find(PAM)
   if(i == -1):
       break
   if((len(Temp_Sequence)-i) < Guide_RNA_length):
       break
   Position = Position + i + 1
   if(Position > Guide_RNA_length):
       if(Direction > 0):
           Location.append(Position - 2)
           Strand.append(Direction)
           Guide_RNAs.append(Sequence[Position-26:Position+4])
       if(Direction < 0):
           Location.append(Position + 1)
           Strand.append(Direction)
           Guide_RNAs.append(Sequence[Position-4:Position+26])
   Temp_Sequence = Temp_Sequence[i+1:]

 if(Direction < 0):
     for i in range(len(Guide_RNAs)):
         Guide_RNAs[i] = str(Guide_RNAs[i].reverse_complement())
 #predictions = model_comparison.predict(np.array(Guide_RNAs))

 #Sorting the Guides based on the Azimuth scores and then only selecting the
 #top number of guides as defined by the Azimuth Cutoff
 #Guide_RNAs = [x for y,x in sorted(zip(predictions,Guide_RNAs))][:Azimuth_Cutoff]
 #Location = [x for y,x in sorted(zip(predictions,Location))][:Azimuth_Cutoff]
 #Strand = [x for y,x in sorted(zip(predictions,Strand))][:Azimuth_Cutoff]

 #Trimming the extra basepairs need for the Azimuth Model from the Guide RNAs
 #so that they fit into the Salis Lab's modelself.

 for i in range(len(Guide_RNAs)):
     Guide_RNAs[i] = Guide_RNAs[i][5:25]

 return Guide_RNAs,Location,Strand

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

Target_Seq, Genome = Get_Sequence()

Genome = Genome + Genome.reverse_complement()
SeqIO.write(Genome, "Total_Genome_Plus_RC", "fasta")

#Obtain the Guide RNAs from the Target Sequence
from Azimuth_Finder import * #placed here to avoid multithreading issues with Tkinter *I think
PAM = "GG"
Guide_Num_Cutoff = 20

T_Guides_Pos, Position_Pos, Direction_Pos,T_Guides_Neg, Position_Neg, Direction_Neg  = Azimuth_Guides(Target_Seq, PAM,Guide_Num_Cutoff )

#Combine the information from the Guide RNAs into one single array.
Position_List = Position_Pos + Position_Neg
Direction_List = Direction_Pos + Direction_Neg
Guide_Info = np.vstack((Position_List, Direction_List)).T

#Combine the two guides
Target_Guides = CombinetoStr(T_Guides_Pos, T_Guides_Neg)

#Send Data to the model
Cas9Calculator=clCas9Calculator(['Total_Genome_Plus_RC'])
sgRNA1 = sgRNA(Target_Guides,Guide_Info, Cas9Calculator)
sgRNA1.run()
