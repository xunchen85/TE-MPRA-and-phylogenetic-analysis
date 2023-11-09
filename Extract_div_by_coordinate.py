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
        list1[line2[0]] = line2

#####
with open(inputFile, 'rt') as f:
    for line in f:
        line2 = re.split(r'\s+',line)
        line3 = re.split(r':',line2[3])
        coordinate = line2[0] + ":" + line2[1] + "-" + line2[2]
        if coordinate in list1:
            print ("\t".join(map(str,line2)) + "\t" + list1[coordinate][4])
        else:
            print ("\t".join(map(str,line2)) + "\tNA")
