#! /usr/bin/python3

from Bio import SeqIO
from Bio.Seq import MutableSeq
import sys

#Error Handling
if len(sys.argv) != 5:
    sys.exit(sys.argv[0]+" <kmer count file> <file to correct> <kmer size> <threshold>")

if int(sys.argv[4]) <= 0:
    sys.exit ("Invalid: threshold size")
if int(sys.argv[3]) <= 0:
    sys.exit ("Invalid Kmer size")
    
# Global dictionary for our kmers
kmers = {}


def isErroneous(pos, size, threshold, s):
   
 
    maxKmerCount = 0 #max kmer frequency 
    
    start = pos-(size-1) # start postion
    if start < 0: # if start postion is less than zero default to zero 
        start = 0
    
    end = pos # end is equal of postion
    if end > len(s)-(size): #if postion is greater than read size - kmer lenght default to size-kmer lenght 
        end = len(s)-(size)
        
        
    for x in range (start,end+1): # excute this loop for all possible kmers 
        kmer_seq = s[x:x+size] #window
        kmer_seq = str(kmer_seq) #conver seq to string 
        
        try: # try as kmer might not be in dictionary 
            if kmers[kmer_seq] > maxKmerCount:
                maxKmerCount = kmers[kmer_seq]
            
        except:
            pass 
            
    if maxKmerCount > threshold: 
        return(False,maxKmerCount)
    else:
        return(True,maxKmerCount) #modified output to display max frequency of kmer with base 


def getAlternativeBase(pos,size,threshold,s):
    bases = ['A','G','C','T']
    
    newseq = list(s)
    count = 0
    maxCount = 0 #count of the base
    maxCountBase = "" #base that has the max count 

   
    for base in bases:
        newseq[pos] = base #replace possible error base with another base 
        strseq = ''.join(str(i) for i in newseq) #turn the list into a string 
        
        bol,count = isErroneous(pos, size, threshold, strseq)
        
        if count > maxCount:
            maxCountBase = base
            maxCount = count
   
    return (maxCountBase) # return the base with the highest count 


def setUpKmers():
   
    fileName = sys.argv[1]
    
    try:
        fastaFile = SeqIO.parse(fileName, "fasta") # open kmer count file 
    except:
       print("ERROR: Unable to Open File :(") #if the file cannont be opened exit 
       sys.exit()

    for record in fastaFile: # cycle through each kmer record and put the data into the kmers dictionary 
        kmers[str(record.seq)] = int(record.id) #key = kmer seq, value = kmer count

def main():
    setUpKmers()
    kmer_size = int(sys.argv[3])
    threshold = int(sys.argv[4])
    outfile = open("outfile.fastq","w")
    with open(sys.argv[2],'r') as f:
        seqIO = SeqIO.parse(f,'fastq')
        for sr in seqIO:
            s = MutableSeq(sr.seq)
            for pos in range(0,len(s)):
                bol,count = isErroneous(pos, kmer_size, threshold, s)
                if bol:
                    newbase = getAlternativeBase(pos,kmer_size,threshold,s)
                    oldbase = s[pos]
                    s[pos] = ("-")
                    print("Sequence:",s, "Position:",pos,"Erroneous Base:", oldbase, "Corrected Base:", newbase)
                    s[pos] = newbase
            sr.seq = s
            SeqIO.write(sr,outfile,'fastq')
    outfile.close()
	
if __name__ == "__main__":
    main()
