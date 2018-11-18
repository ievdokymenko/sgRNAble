from Azimuth_Model import model_comparison
import numpy as np
from Bio.Seq import Seq



def PAM_Finder(Sequence, PAM, Direction, Cutoff):
  Guide_RNAs = []
  Location = []
  Strand = []

  Azimuth_Cutoff = Cutoff
  Guide_RNA_length = 30

  Position = 0
  Temp_Sequence = Sequence
  if(Direction < 0):
    PAM = str(Seq(PAM).reverse_complement())
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
  if(Direction > 0):
      for i in range(len(Guide_RNAs)):
           Guide_RNAs[i] = str(Guide_RNAs[i])

  predictions = model_comparison.predict(np.array(Guide_RNAs))
  #Sorting the Guides based on the Azimuth scores and then only selecting the
  #top number of guides as defined by the Azimuth Cutoff
  Guide_RNAs = [x for y,x in sorted(zip(predictions,Guide_RNAs))][:Azimuth_Cutoff]
  Location = [x for y,x in sorted(zip(predictions,Location))][:Azimuth_Cutoff]
  Strand = [x for y,x in sorted(zip(predictions,Strand))][:Azimuth_Cutoff]

  #Trimming the extra basepairs need for the Azimuth Model from the Guide RNAs
  #so that they fit into the Salis Lab's modelself.
  for i in range(len(Guide_RNAs)):
      Guide_RNAs[i] = Guide_RNAs[i][4:24]

  return Guide_RNAs,Location,Strand

def Azimuth_Guides(Gene_Seq, PAM_Seq, Guides_Cutoff):

    Gene1,Location1,Direction1 = PAM_Finder(Gene_Seq, PAM_Seq, 1,Guides_Cutoff)
    Gene2,Location2,Direction2 = PAM_Finder(Gene_Seq, PAM_Seq, -1,Guides_Cutoff)


    return Gene1,Location1,Direction1,Gene2,Location2,Direction2
