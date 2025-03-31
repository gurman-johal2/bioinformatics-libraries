#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 20:02:38 2021

@author: gurmanjohal
"""

#! /usr/bin/python

import sys # needed for command line arguement access 
from Bio import SeqIO #file processing method taken from lecture 1 - biopython 

#---- Parameter Check --------------

fileName = sys.argv[1] 

#check if the file can be opened, if not exit system
try:
    fastaFile = SeqIO.parse(fileName, "fasta") #check if file can open
except:
    print("ERROR: Unable to Open File :(")
    sys.exit()

#Check if threshold is a number
try:
    threshold = float(sys.argv[2])
except:
    print("Threshold is not a number :(")
    sys.exit()

#check threshold range
if threshold > 1 or threshold <0 :
    print("Invalid Threshold (must be value between 0 and 1 ")
    sys.exit()


print("No errors with the parameters")
print(f"File Name: {fileName} Threshold: {threshold}")


#process file into 2D list ------------------------------

data = [] #data file will be 2D list (1D: record, 2D: characters)

for record in fastaFile:
    seq =list(record.seq) #convert the sequence into a list 
    data.append(seq) #add sequence to data list 
    
# misc info ------------------------------

nRecords = len(data)
nBases=len(data[0])

print("Number of Reads: ",nRecords, "Read Lenght: ", nBases)

# print header ------------------------
print ("\n{:<10} {:^15} {:^8} {:^15} {:^8} {:^15}".format("Postion", "Major Allele", "count", "Minor Allele", "Count", "Frequency"))


countG = 0
countC = 0
countA = 0
countT = 0

postion = 0 #nucletoide postion 
record = 0 # record number 


#get base frequency data based on postion ---------------------------------------------------
# sort through postion then by each record collecting the base frequecy values 
for postion in range (nBases):
    for record in range (nRecords):
        base = data[record][postion] 
        if base == "G":
           countG += 1
        if base == "C":
           countC += 1
        if base == "A":
           countA += 1
        if base == "T":
           countT += 1
           
    total = countG + countC + countA + countT #total number of bases counted 
    #total = nRecords # total value used for the frequnecy calculation is equal to the number of reads (read with _ or x inculded)
    
    basefrequency = {"G":countG,"C":countC,"A":countA,"T":countT} # put information into a dictonary 
    
    basefrequency = sorted(basefrequency.items(), key=lambda x: x[1], reverse=True) #sort bases based on count
    
    basefrequency_list = []
    
    for key, value in basefrequency: # convert sotrted base and base count values into list from dictionary 
        basefrequency_list.append(key)
        basefrequency_list.append(int(value))
    
    
    if ((basefrequency_list[3]/total) > (threshold)):
        print ("{:<10} {:^15} {:^8} {:^15} {:^8} {:^15.2f}".format(postion,basefrequency_list[0],basefrequency_list[1],basefrequency_list[2],basefrequency_list[3],(basefrequency_list[3]/total)))
    
    #reset counters 
    countG = 0
    countC = 0
    countA = 0
    countT = 0
    total = 0

