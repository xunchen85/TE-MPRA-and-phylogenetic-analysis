import sys
import getopt
import re

##### Author: xunchen85@gmail.com
##### This script is used to extract sequences in a given position and length of one TE subfamily for MPRA experiment
##### Date: 2021/1/24

try:
    opts,args = getopt.getopt(sys.argv[1:], '-h:-i:-S:-E:-b:-n:-c:-o:', ['help', 'iseqFile=', 'Start=', 'End=', 'bedFile=', 'nameSubfamily=','consensusSubfamily=','outFile='])
except getopt.GetoptError:
    print ('python command variables')
    sys.exit()

for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print("[*] Help info")
        sys.exit()
    if opt_name in ('-i','--iseqFile'):
        iseqFile = opt_value
    if opt_name in ('-S','--Start'):
        Start = int(opt_value)
    if opt_name in ('-E','--End'):
        End = int(opt_value)
    if opt_name in ('-b','--bedFile'):
        bedFile = opt_value
    if opt_name in ('-n','--nameSubfamily'):
        nameSubfamily = opt_value
    if opt_name in ('-c','--consensusSubfamily'):
        consensusSubfamily = opt_value
    if opt_name in ('-o','--outFile'):
        outFile = opt_value

#iseqFile = "/home/xchen/Desktop/Overall_Data/TE_rmsk_profilling.Oct.2/hg19.fa.align.seq"
#bedFile = "/home/xchen/OneDriveU/Projects_Bourque/EMC_project_2020_12_1/Refs/hg19_rmsk_TE_1bp.bed"
#bedFile = "MER41B.bed"
#Start = 1
#End = 100
#nameSubfamily = "MER41B"
#consensusSubfamily = "MER41B"

#####
TElist = {}

##### read bedFile
with open(bedFile, 'rt') as b:
    for line in b:
        #print (line)
        line2 = re.split(r'\s+',line)
        TEinfo = re.split(r':',line2[3])
        if TEinfo[0] == nameSubfamily:
            TEcoordinate_1bp = line2[0] + ":" + str(line2[1]) + "-" + str(line2[2])
            TElist[TEcoordinate_1bp] = TEinfo[1]
            #print (TEcoordinate_1bp)
print("TEinstance_coordinate\tDirectionAgainstConsensus\tTEConsensusName\tTEgroup\tPosiConsensus_1\tPosiConsensus_2\tPosiConsensus_3\tFullSeqHuman\tFullSeqConsensus\tTEsubfamily\tTEinstanceName\tFrameRefCoordinate\tFrameRefSeq\tFrameRefLen\tFrameConsensusStart\tFrameConsensusEnd\tFrameConsensusSeq\tFrameConsensusLen")
with open(iseqFile, 'rt') as i:
    for line in i:
        #print (line)
        line2 = re.split(r'\s+',line)
		#print (line2)
        coordinate_tmp = line2[0].replace(":","-")
        coordinates = re.split(r'[-:]',coordinate_tmp)
        coordinates[1] = int(coordinates[1])
        coordinates[2] = int(coordinates[1])
        #print (line2[0],line2[2],consensusSubfamily)
        if line2[0] in TElist and line2[2] == consensusSubfamily:
            seq_r = line2[7]
            seq_c = line2[8]
            seq_r2 = ""
            seq_c2 = ""
            posi = 0
            if line2[1] == "+":
                Start1 = int(line2[4])
                End1 = int(line2[4])
            else:
                Start1 = int(line2[5])
                End1 = int(line2[5])
            for i in range(0,len(line2[8]),1):
                if End1 >= Start and End1 <= End:
                    seq_r2 = seq_r2 + line2[7][i]
                    seq_c2 = seq_c2 + line2[8][i]
                    if line2[1] == "+" and line2[8][i] != "-":
                        End1 = End1 + 1
                    elif line2[1] == "C" and line2[8][i] != "-":
                        End1 = End1 - 1
                    if line2[7][i] != "-":
                        coordinates[2] = coordinates[2] + 1
                elif Start1 < Start or Start1 > End:
                    if line2[1] == "+" and line2[8][i] != "-":
                        Start1 = Start1 + 1
                        End1 = End1 + 1
                    elif line2[1] == "C" and line2[8][i] != "-":
                        Start1 = Start1 - 1
                        End1 = End1 - 1
                    if line2[7][i] != "-":
                        coordinates[1] = coordinates[1] + 1
                        coordinates[2] = coordinates[1]
            if len(seq_r2) == 0:
                seq_r2 = "N"
                seq_c2 = "N"
                Len_r2 = 0
                Len_c2 = 0
                if line2[1] == "+":
                    Start1 = int(line2[4])
                    End1 = int(line2[5])
                else:
                    Start1 = int(line2[5])
                    End1 = int(line2[6])
            else:
                if line2[1] == "+":
                    End1 = End1 - 1
                else:
                    End1 = End1 + 1
                coordinates[2] = coordinates[2] - 1
                Len_r2 = len(seq_r2.replace("-",""))
                Len_c2 = len(seq_c2.replace("-",""))
            print (line.rstrip() + "\t" + nameSubfamily + "\t" + TElist[line2[0]] + "\t" + coordinates[0] + ":" + str(coordinates[1]) + "-" + str(coordinates[2]) + "\t" + seq_r2 + "\t" + str(Len_r2) + "\t" + str(Start1) + "\t" + str(End1) + "\t" + seq_c2 + "\t" + str(Len_c2))
