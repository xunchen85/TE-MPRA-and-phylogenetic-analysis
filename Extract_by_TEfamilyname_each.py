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

#####
Order = 1
OrderS = 1
with open(inputFile, 'rt') as f:
    for line in f:
        line2 = re.split(r'\s+',line)
        line3 = re.split(r':',line2[3])
        if line3[0] == listFile and int(line2[2]) - int(line2[1]) >= 200:
            line2[4] = line3[3]
            line2[3] = "T_" + str(Order)
            print("\t".join(line2))
            Order = Order +1
        elif line3[0] == listFile and int(line2[2]) - int(line2[1]) < 200:
            line2[4] = line3[3]
            line2[3] = "S_" + str(OrderS)
            print("\t".join(line2))
            OrderS = OrderS +1
