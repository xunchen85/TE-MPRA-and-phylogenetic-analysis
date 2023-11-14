#######################
### Author: Xun Chen, Ph.D.
### Email: xunchen85@gmail.com or xchen@outlook.com
### ORCID: https://orcid.org/0000-0003-0327-1888
### Date: 2021/3/16
###
#######################

#!/usr/bin/python
# -*- coding: UTF-8 -*-

import gzip
import re
import sys
import getopt
import numpy as np

### define variables

try:
	opts,args = getopt.getopt(sys.argv[1:], '-h:-i:-b:-c:-t:-m:', ['help', 'insert=','bcFile=','maxbc=','cell=','typeCol='])
except getopt.GetoptError:
	sys.exit()

for opt_name,opt_value in opts:
	if opt_name in ('-h','--help'):
		sys.exit()
	if opt_name in ('-i','--insert'):
		insertFile = opt_value
	if opt_name in ('-b','--bcFile'):
		bcFile = opt_value
	if opt_name in ('-m','--maxbc'):
		maxbc = int(opt_value)
	if opt_name in ('-c','--cell'):
		cellFile = opt_value
	if opt_name in ('-t','--typeCol'):
		typeCol = int(opt_value)

#print("Insert\tTotalBCperInsertOriginal\tBC\tDNAcount\tRNAcount")

#### insert list 
insertList = {}
with open (insertFile, 'rt') as f:
	for line in f:
		line2 = re.split(r'\s+',line.rstrip())
		ID = line2[0]
		insertList[ID] = [0] * (maxbc+1)
		insertList[ID][0] = ID
		#print (ID)

#### bc list
bcList = {}
with open (bcFile, 'rt') as f:
	for line in f:
		line2 = re.split(r'\s+',line.rstrip())
		if line2[0] == "Insert":
			next
		else:
			ID = line2[0] + "|" + line2[1]
			bcList[ID] = line2[2]
		#print (ID)

#### cell list
cellList = {}
with open (cellFile, 'rt') as f:
	for line in f:
		line2 = re.split(r'\s+',line.rstrip())
		if line2[0] == "Insert":
			next
		else:
			ID = line2[0] + "|" + line2[2]
			insertList[line2[0]][int(bcList[ID])] = line2[typeCol]
		#	print (insertList[line2[0]])
print("insert" + "\t" + "\t".join(map(str,range(1,(maxbc+1)))))
for ID in insertList:
	print ("\t".join(map(str,insertList[ID])))
