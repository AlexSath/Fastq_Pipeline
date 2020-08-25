#this main script needs to import everything that was imported within other scripts for it to work. Thus the import is a collection of all the imports needed for the script to run. 
import argparse
import createfastq as cf
import os
import part23 as p23 
import part45 as gm
#Example use is 
# python parseFastq.py --fastq /home/rbif/week6/hawkins_pooled_sequences.fastq

#defining variables filePath and filedir while using .os command for absolutle path while interfacing operating system
filePath = os.path.abspath(__file__) 
filedir = os.path.dirname(filePath)


#defining the functino of main: which will be used to execute all other scripts in this master script.
def main():
    fastqPath = os.path.join(filedir,"DATA", "hawkins_pooled_sequences.fastq") #variable while using OS to join the original filedir variable with the hawkins pooled seq data
    clinicalpath = os.path.join(filedir,"DATA", "harrington_clinical_data.txt") #variable while using OS to join the original filedir variable with the hawkins pooled seq data
    cf.createfastqs(fastqPath, clinicalpath) #above we imported this script as cf. we are executing the script upon the variables of fastqPath and clinicalpath
    
    p23.createSams() 
    p23.processBams() 
    bamdir = os.path.join(filedir,"sortedBAM")
    outfile = os.path.join(filedir, "tempout.tsv")
    outfile2 = os.path.join(filedir, "report.txt")
    gm.popstats(bamdir,outfile)
    gm.compilereport(outfile,outfile2)

if __name__ == '__main__':
    main()



