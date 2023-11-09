import sys
import getopt
import re

#####
try:
    opts,args = getopt.getopt(sys.argv[1:], '-h:-i:-l:-g:-n:', ['help', 'inputFile=', 'listFile=', 'groupFile=', 'nameFile='])
except getopt.GetoptError:
    sys.exit()

    
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print("[*] Help info")
        sys.exit()
    if opt_name in ('-i','--inputFile'):
        inputFile = opt_value
    if opt_name in ('-l','--listFile'):
        listFile = opt_value
    if opt_name in ('-g','--groupFile'):
        groupFile = opt_value  
    if opt_name in ('-n','--nameFile'):
        nameFile = opt_value

##### list of family
list1 = {}
with open(listFile, 'rt') as f:
    for line in f:
        line2 = re.split(r'\s+',line)
        list1[line2[0]] = [0,0,0]
        list1[line2[0]][0] = line2[1]
        list1[line2[0]][1] = 0
        list1[line2[0]][2] = 0

##### list of instance name
name1 = {}
with open(nameFile, 'rt') as n:
    for line in n:
        line2 = re.split(r'\s+',line)
        coordinate = line2[0] + ":" + line2[1] + "-" + line2[2]
        name1[coordinate] = line2[3]
        #print (coordinate,name1[coordinate])

##### lift of group
group1 = {}
list_group1 = {}
with open(groupFile, 'rt') as g:
    for line in g:
        line2 = re.split(r'\s+',line)
        if line2[20] == line2[0]:
            line2[24] = line2[24].replace("macFas5_rmsk_","")
            list_group1[line2[20]] = [line2[24],line2[21],line2[23],line2[1],line2[22],0,0,line2[15]]
            rename = line2[24] + ":" + line2[7]
            group1[rename] = line2[20]
        elif line2[19] == line2[0]:
            line2[23] = line2[23].replace("macFas5_rmsk_","")
            list_group1[line2[19]] = [line2[23],line2[20],line2[22],line2[1],line2[21],0,0,line2[15]]
            rename = line2[23] + ":" + line2[7]
            group1[rename] = line2[19]
        else:
            line2[24] = line2[24].replace("macFas5_rmsk_","")
            list_group1[line2[20]] = [line2[24],line2[21],line2[23],line2[1],line2[22],"error","error",line2[15]]
        #print (rename,group1[rename])
        #print(line[20],list_group1[line2[20]])

##### liftover file
best_matched = {}
with open(inputFile, 'rt') as f:
    for line in f:
        line2 = re.split(r'\s+',line)
        if line2[10] == "":
            next
        elif int(line2[10]) > 0:
            line_original = re.split(r':',line2[3])
            line_liftover = re.split(r':',line2[9])
            line_liftover_name = ":".join(line_liftover[0:4])
            family_name = line_original[0] + ":" + line_original[2] + ":" + line_original[3]
            TE_len = int(line2[2]) - int(line2[1])
            if line2[3] == line_liftover_name:
                # coordinate
                coordinate = line2[0] + ":" + line2[1] + "-" + line2[2]
                # rename
                if coordinate in name1:
                    rename = line_original[0] + ":" +name1[coordinate]
                    #print (rename)
                else:
                    rename = "not_tested"
                if line2[3] not in best_matched:
                    best_matched[line2[3]] = TE_len
                    list1[family_name][1] += 1
                    if rename in group1:
                        list_group1[group1[rename]][5] += 1
                    if TE_len >= 200:
                        list1[family_name][2] += 1
                        if rename in group1:
                            list_group1[group1[rename]][6] += 1
                else:
                    if best_matched[line2[3]] < 200 and TE_len >=200:
                        best_matched[line2[3]] = TE_len
                        list1[family_name][2] += 1
                        if rename in group1:
                            list_group1[group1[rename]][6] += 1

fileName = inputFile.replace("hg19_rmsk_TE_0bp.","").replace(".back.intersect.bed","")
fileName = fileName.replace("macFas5_rmsk_TE_0bp.","")
#print ("fileName\tTEfamily\ttotalCount\tliftoverintersect_all\tliftoverintersect_200bp")
#for TE in list1:
#    if TE != "TE_family":
#        print(fileName + "\t" + TE + "\t" + "\t".join(map(str,list1[TE])))
for TE in list_group1:
    if TE != "TE_family":
        print(fileName + "\t" + "\t".join(map(str,list_group1[TE])))
