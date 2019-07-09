#!/usr/bin/env python3

import numpy as np
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
from Cas9_Calculator import *
from File_GUI import Get_Sequence

if __name__ == "__main__":

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

    if Usage == "Guide Screening":
        #Convert all Sequences into Strings
        Target_Seqs = list(map(lambda Target_Seqs:str(Target_Seqs.seq),Target_Seqs))
        Guide_Info = [["NA","NA"]] *len(Target_Seqs)

        #Combine the Genome
        Genome = Genome.upper() + Genome.upper().reverse_complement()
        SeqIO.write(Genome, "Total_Genome_Plus_RC", "fasta")

        #Send Data to the model
        Cas9Calculator=clCas9Calculator(['Total_Genome_Plus_RC'])
        sgRNA1 = sgRNA(Target_Seqs,Guide_Info, Cas9Calculator)
        sgRNA1.run()
        sys.exit()

    #Obtain the Guide RNAs from the Target Sequence
    from Azimuth_Finder import * #placed here to avoid multithreading issues with Tkinter (I think)
    PAM = "GG"
    Guide_Num_Cutoff = 20 # Cutoff per strand

    Target_Guides, Position_List, Direction_List, Con_Ids, Met_Distance, Con_ORFs, Azimuth_Score = Azimuth_Guides(Target_Seqs, PAM,Guide_Num_Cutoff, Usage )

    if Usage == "CRISPRi Screening":

        Guide_Cutoff = 5
        workbook = xlsxwriter.Workbook('Guide RNAs.xlsx')
        worksheet = workbook.add_worksheet()
        row = 0
        col = 0
        ORF_Num = 0
        Maximum_Met_Distance = 500
        Minimum_Gene_Distance = 600
        for k in range(len(Con_Ids)):
            worksheet.write(row,col, "Contig ID:")
            worksheet.write(row,col + 1 , Con_Ids[k])
            if (Con_ORFs[k] == 0):
                continue
            for i, ORF in enumerate(Con_ORFs[k]):
                worksheet.write(row + 1,col + 1 , "ORF Length: ")
                worksheet.write(row + 1,col + 2, len(ORF))
                worksheet.write(row + 1,col + 3 , "ORF Sequence: ")
                worksheet.write(row + 1,col + 10, str(ORF))
                Guides_Listed = 0
                for j in range(len(Target_Guides[ORF_Num + i])):
                    if (- Met_Distance[i][j][0] > Position_List[i][j]):
                        continue
                    if (- Met_Distance[i][j][0] > Maximum_Met_Distance):
                        continue
                    if (len(ORF) + Met_Distance[i][j][0] < Minimum_Gene_Distance):
                        continue
                    if(Guides_Listed > Guide_Cutoff):
                        break
                    else:
                        Guides_Listed += 1
                    worksheet.write(row + 2, col +2, "Guide:")
                    worksheet.write(row + 2, col +3 , Target_Guides[ORF_Num + i][j])
                    worksheet.write(row + 2, col +4, "Position: ")
                    worksheet.write(row + 2, col +5, Position_List[ORF_Num + i][j])
                    worksheet.write(row + 2, col +6, "Strand: ")
                    worksheet.write(row + 2, col +7, Direction_List[ORF_Num + i][j])
                    worksheet.write(row + 2, col +8, "Distance to Previous Methionine (bp)")
                    worksheet.write(row + 2, col +9, Met_Distance[ORF_Num + i][j][0])
                    worksheet.write(row + 2, col +10, "Distance to Next Methionine (bp)")
                    worksheet.write(row + 2, col +11, Met_Distance[ORF_Num + i][j][1])
                    worksheet.write(row + 2, col +12, "Azimuth Model Score")
                    worksheet.write(row + 2, col +13, Azimuth_Score[ORF_Num + i][j])
                    row += 2
            ORF_Num += i
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
