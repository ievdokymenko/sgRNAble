from Azimuth_Model import model_comparison
import numpy as np
from Bio.Seq import Seq
from bisect import bisect_left

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
  #so that they fit into the Salis Lab's model.
  for i in range(len(Guide_RNAs)):
      Guide_RNAs[i] = Guide_RNAs[i][4:24]

  return Guide_RNAs,Location,Strand

def Azimuth_Guides(Gene_Seq, PAM_Seq, Guides_Cutoff, Usage):

    if (Usage == "CRISPRi Screening"):
        Position = []
        Direction = []
        Contig_Name = []
        Guides = []
        Frames = [1,2,3,-1,-2,-3]
        Stop_Codons = ["TAG","TGA","TAA"]
        Met = "ATG"
        Reverse_Met = "TAC"
        Start_Codons = [0]
        Reverse_Start_Condons = [0]
        ORFs = []
        ORF_Cutoff = 1000
        Location_in_Contig = 0
        Distance_to_Nearest_Met = []

        Guides_Cutoff = 10

        for Contig in Gene_Seq:
            for Frame in Frames:
                Contig_Frame = Contig[abs(Frame)-1:]
                if Frame < 0:
                    Contig_Frame = Contig[abs(Frame)-1:].reverse_complement()
                for i in range(len(Contig_Frame)):
                    if str(Contig_Frame[i*3:i*3+3].seq) in Met:
                        Start_Codons.append(i*3)

                    if str(Contig_Frame[i*3:i*3+3].seq) in Reverse_Met:
                        Reverse_Start_Condons.append(i*3)

                    if Contig_Frame[i*3:(i*3)+3].seq in Stop_Codons:

                        if(i*3-Location_in_Contig) > ORF_Cutoff and Location_in_Contig != 0:
                            Gene1, Location1, Direction1 = PAM_Finder(Contig_Frame[Location_in_Contig:(i*3)+3].seq, PAM_Seq, 1,Guides_Cutoff)
                            Position.append(Location1)
                            Direction.append(Direction1)
                            Meth_Distances = []
                            for Loc in Location1:
                                if Frame > 0:
                                    Meth_Distances.append(ClosestDiff(Start_Codons,Loc+Location_in_Contig))
                                else:
                                    Meth_Distances.append(ClosestDiff(Reverse_Start_Condons,Loc+Location_in_Contig))
                            Distance_to_Nearest_Met.append(Meth_Distances)
                            Guides.append(Gene1)
                            ORFs.append(Contig_Frame[Location_in_Contig:(i*3)+3].seq)
                        Location_in_Contig = i*3
                Location_in_Contig = 0
            Contig_Name.append(Contig.id)

        return Guides, Position, Direction, Contig_Name, Distance_to_Nearest_Met, ORFs

    elif (Usage == "CRISPRi on a Gene"):
        for Gene in Gene_Seq:
            Gene,Location,Direction = PAM_Finder(Gene.seq.upper(), PAM_Seq, -1,Guides_Cutoff)

        return Gene, Location, Direction, 0 ,0, 0

    else:

        Gene1,Location1,Direction1 = PAM_Finder(Gene_Seq, PAM_Seq, 1,Guides_Cutoff)
        Gene2,Location2,Direction2 = PAM_Finder(Gene_Seq, PAM_Seq, -1,Guides_Cutoff)

        for i in range(len(Location1)):
            if (Location1[i] > Usage):
                Location1 = Location1[:i]
                Gene1 = Gene1[:1]
                Direction1 = Direction1[:i]

        for i in range(len(Location2)):
            if (Location2[i] > Usage):
                Location2 = Location2[:i]
                Gene2 = Gene2[:1]
                Direction2 = Direction2[:i]

        Position = Location1 + Location2
        Direction = Direction1 + Direction2

        #Combine the two guides
        Guides = CombinetoStr(Gene1, Gene2)

        return Guides, Position, Direction, 0, 0, 0

def ClosestDiff(myList, myNumber):
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return ["NA",(myNumber - myList[pos + 1])]
    if pos == len(myList):
        return [(myNumber - myList[pos - 1 ]),"NA"]
    before = myList[pos - 1]
    after = myList[pos]
    return [(myNumber - before),(myNumber - after)]
