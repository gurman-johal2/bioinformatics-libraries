#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 10:37:39 2021

@author: gurmanjohal
"""
#imports 
import sys # needed for command line arguement access 
import csv # nneeded for CSV output 

#comand Line arguement 
fileName = sys.argv[1] #unzip file first 
windowSize = int(sys.argv[2])
increment = int(sys.argv[3])
count = 0

print("FileName: ", fileName, "windowSize: ", windowSize, "Increment: ", increment)
#lists 
line = []
data = []

#set default values for density values based on patient ID
I=0
C=0
B=0
A=0
P=0
D=0
F=0
E=0


#check if the file can be opened, if not exit system
try:
    file = open(fileName,mode= "r") 
except:
    print("ERROR: Unable to Open File :(")
    sys.exit()


#process file and create a 2D array 1D: various records, 2D: polymorphsim data for the patients 
for record in file:
        if record.find("#")==-1: #if the data is actually data and not a header (# or ##)
            dataSpilt = record.split("\t")
            
            line.append(int(dataSpilt[1])) # add the postion
         
            for i in range(9,17): #loop through patient data 
                    if dataSpilt[i].find("0/0")!=-1:
                        line.append(0) #zero if no polymorishim detected
                    else:
                        line.append(1) #one if detected
            data.append(line)
            line = [] #blank the line 
            

# inital setup for writing loop      
start = data[0][0] #lowest postion in file 
end = start + windowSize #window range 
absEnd = data[len(data)-1][0] #highest postion in file 


#CSV OUTPUT FILE LOOP --------------------------------------
with open("output-assignment2-GJ.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    header = ["Postion","Patient I","Patient C","Patient B","Patient A","Patient P","Patient D","Patient F","Patient E"]
    writer.writerow(header)
    while end < absEnd: # only run if the end value of iteration is below to last postion mapped 
        line.append(start) #add the start value for the iteration
        for z in range (0,len(data)-1):#loop through all lines  
            if data[z][0] >= start and data[z][0] <= end: #only process those whose postion is within the iteration range
                I+=data[z][1] #count number polymohsims detected for the patient witin the window 
                C+=data[z][2]
                B+=data[z][3]
                A+=data[z][4]
                P+=data[z][5]
                D+=data[z][6]
                F+=data[z][7]
                E+=data[z][8]
                
        #add N polymorphsims to the list 
        line.append(I)
        line.append(C)
        line.append(B)
        line.append(A)
        line.append(P)
        line.append(D)
        line.append(F)
        line.append(E)
        
        #reset values for next iteration 
        I=0
        C=0
        B=0
        A=0
        P=0
        D=0
        F=0
        E=0
        
        writer.writerow(line) #write the line containg frequency and postion data 
        
        line = [] #clear exisitng values 
        start+=increment #update new start value 
        end = start + windowSize #update new end value 
        
print("Code Complete: Output file Generated")     
csvfile.close() #close file 


  
  
    

