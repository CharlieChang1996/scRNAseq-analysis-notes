#In this example we load two example files
#"Oct4.pos.fasta" contains sequence enriched from ChIP-seq experiment (peak centres) using Oct4, a transcription factor
#"bground.fasta" contains the background (non-enriched) sequences from the same experiment
#This code reads these files and searches for enrichment using Python dictionaries
#Complete the analysis and output the results of sorted 8bp enriched sequences
#Note we use 8bp because Oct4 is an "Octamer binding protein", it is known to bind 8mers
#The most common known binding site for Oct4 is ATCGCAAAT

# Complete the missing code in the compareDictionaries(kmer) method- each missing statement is indicated by ???
#Open the files and interpret the results-do you find the expected enrichment patterns
#How might you combine the results of the foreground and background searches?
#What happens if you change the kmer size?

#!/usr/python3
from operator import itemgetter
import re

def normaliseDict(mydict):

    #count total reads
    #express as percent total

    count=0

    for myseq in mydict.keys():
        count =count +mydict[myseq]
    count =count/100

    for myseq in mydict.keys():
        mydict[myseq]=mydict[myseq]/count



def processSeq(mydict,myline,kmer):
    #go though the line chopping sequence
    p=re.compile('[^ATGC]',re.IGNORECASE)

    for i in range(0, len( myline)-kmer):
        #do something
        myword=myline[i:i+kmer]

        #we ignore issues of case here
        if myword.find("N")!=-1:
            continue
        if len(myword)!=kmer:
            continue
        if p.match(myword):
            continue;


        #print(myword)

        if myword not in mydict:
            mydict[myword] = 1
        else:
            mydict[myword] += 1

def buildDictionary(fastafile,kmer):
    f1 = open(fastafile, 'r')
    f1lines=f1.readlines()
    readBody=False
    myline=""
    mydict ={}

    for line in f1lines:
        if line.startswith('>'):
            if readBody==True:
                readBody=False
                #we are ending the processing of a seq
                processSeq(mydict,myline,kmer)
            else:
                #we are starting the processing of a seq
                readBody=True
            continue

        if readBody==False:
            myline=""
        else:
            if len(line)>0:
                myline=myline+line;
    f1.close()
    return mydict

def compareDictionaries(kmer):
    #just load the two files and count each kmer
    dictforeground= buildDictionary("Oct4.pos.fasta",8)
    dictbackground = buildDictionary("bground.fasta",8)

    #normalise by size
    normaliseDict(dictforeground)
    normaliseDict(dictbackground)

    F = open("testfile_fg.txt", "w")
    for myseq in sorted(dictforeground.items(), key=itemgetter(1), reverse=True):
      #  print(myseq[0])
        if(len(myseq[0])!=kmer):
            continue
        F.write(myseq[0])
        F.write("\t")
        F.write(str(myseq[1]))
        F.write("\n")
    F.close()

    F = open("testfile_bg.txt", "w")
    for myseq in sorted(dictbackground.items(), key=itemgetter(1), reverse=True):
        #  print(myseq[0])
        if (len(myseq[0]) != kmer):
            continue
        F.write(myseq[0])
        F.write("\t")
        F.write(str(myseq[1]))
        F.write("\n")
    F.close()

#actually run the programme
compareDictionaries(8)
