#!/usr/bin/env python3

import numpy as np
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
from Cas9_Calculator import *
from File_GUI import Get_Sequence

#Length of the Guide RNA desired
Guide_RNA_length = 20

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

Target_Seqs, Genome, Usage = Get_Sequence()

#Obtain the Guide RNAs from the Target Sequence
from Azimuth_Finder import * #placed here to avoid multithreading issues with Tkinter (I think)
PAM = "GG"
Guide_Num_Cutoff = 20 # Cutoff per strand

Target_Guides, Position_List, Direction_List, Con_Ids, Met_Distance, ORF = Azimuth_Guides(Target_Seqs, PAM,Guide_Num_Cutoff, Usage )

if Usage == "CRISPRi Screening":
    workbook = xlsxwriter.Workbook('Guide RNAs.xlsx')
    worksheet = workbook.add_worksheet()
    row = 0
    col = 0
    for Name in Con_Ids:
        worksheet.write(row,col, "Contig ID:")
        worksheet.write(row,col + 1 , Name)
        for i in range(len(ORF)):
            worksheet.write(row + 1,col + 1 , "ORF Sequence: ")
            worksheet.write(row + 1,col + 9, str(ORF[i]))
            for j in range(len(Target_Guides[i])):
                worksheet.write(row + 2, col +2, "Guide:")
                worksheet.write(row + 2, col +3 , Target_Guides[i][j])
                worksheet.write(row + 2, col +4, "Position: ")
                worksheet.write(row + 2, col +5, Position_List[i][j])
                worksheet.write(row + 2, col +6, "Strand: ")
                worksheet.write(row + 2, col +7, Direction_List[i][j])
                worksheet.write(row + 2, col +8, "Distance to Previous Methanine")
                worksheet.write(row + 2, col +9, Met_Distance[i][j][0])
                worksheet.write(row + 2, col +10, "Distance to Next Methanine")
                worksheet.write(row + 2, col +11, Met_Distance[i][j][1])
                row += 2
        row += 2

    workbook.close()

else:
    Genome = Genome.upper() + Genome.upper().reverse_complement()
    SeqIO.write(Genome, "Total_Genome_Plus_RC", "fasta")

    #Combine the information from the Guide RNAs into one single array.
    Guide_Info = np.vstack((Position_List, Direction_List)).T

    #Send Data to the model
    Cas9Calculator=clCas9Calculator(['Total_Genome_Plus_RC'])
    sgRNA1 = sgRNA(Target_Guides,Guide_Info, Cas9Calculator)
    sgRNA1.run()
