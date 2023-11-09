import sys
import getopt
import re

#####
try:
    opts,args = getopt.getopt(sys.argv[1:], '-h:-i:-l:', ['help', 'inputFile=', 'listFile='])
except getopt.GetoptError:
    print ('python Association_Genes_TEs.py -i <inputFile> -l <listFile>')
    sys.exit()

    
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print("[*] Help info")
        sys.exit()
    if opt_name in ('-i','--inputFile'):
        inputFile = opt_value
        #print('Input file: ' + inputFile)
    if opt_name in ('-l','--listFile'):
        listFile = opt_value
        #print('Input file: ' + inputFile)

list1 = {}
with open(listFile, 'rt') as f:
    for line in f:
        line2 = re.split(r'\s+',line)
        list1[line2[0]] = [0,0,0]
        list1[line2[0]][0] = line2[1]
        list1[line2[0]][1] = 0
        list1[line2[0]][2] = 0

best_matched = {}
#####
with open(inputFile, 'rt') as f:
    for line in f:
        line2 = re.split(r'\s+',line)
        if int(line2[10]) > 0:
            line_original = re.split(r':',line2[3])
            line_liftover = re.split(r':',line2[9])
            line_liftover_name = ":".join(line_liftover[0:4])
            family_name = line_original[0] + ":" + line_original[2] + ":" + line_original[3]
            TE_len = int(line2[2]) - int(line2[1])
            if line2[3] == line_liftover_name:
                if line2[3] not in best_matched:
                    best_matched[line2[3]] = TE_len
                    list1[family_name][1] += 1
                    if TE_len >= 200:
                        list1[family_name][2] += 1
                else:
                    if best_matched[line2[3]] < 200 and TE_len >=200:
                        best_matched[line2[3]] = TE_len
                        list1[family_name][2] += 1
fileName = inputFile.replace("hg19_rmsk_TE_0bp.","").replace(".back.intersect.bed","")
#print ("fileName\tTEfamily\ttotalCount\tliftoverintersect_all\tliftoverintersect_200bp")
for TE in list1:
    if TE != "TE_family":
        print(fileName + "\t" + TE + "\t" + "\t".join(map(str,list1[TE])))
