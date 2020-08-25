import pysam
import os
import createfastq as cf
import ast
filePath = os.path.abspath(__file__)
filedir = os.path.dirname(filePath)

#This is a slightly modified version from here: https://pysam.readthedocs.io/en/latest/api.html
# What is a pileup? It is "piling up" all the reads that have mapped to a given position of your reference.
# It is a useful way to see what reads have a mutation and what don't. 

def pileup(bampath):
    #test file, replaced with the sorted.bam you are using. Make sure it is indexed! (Use samtools index yourbam.sorted.bam)
    samfile = pysam.AlignmentFile(bampath, "rb")
    #ntdictarr = []
    #Since our reference only has a single sequence, we're going to pile up ALL of the reads. Usually you would do it in a specific region (such as chromosome 1, position 1023 to 1050 for example)
    for pileupcolumn in samfile.pileup():
        #print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        #use a dictionary to count up the bases at each position
        ntdict = {}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # You can uncomment the below line to see what is happening in the pileup. 
                #print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                ########## ADD ADDITIONAL CODE HERE ############# 
                if  base in ntdict.keys():
                    ntdict[base] = ntdict[base]+1
                else:
                    ntdict[base] = 1
                
                # Populate the ntdict with the counts of each base 
                # This dictionary will hold all of the base read counts per nucletoide per position.
                # Use the dictionary to calculate the frequency of each site, and report it if if the frequency is NOT  100% / 0%. 
                #############################################
                
        if len(ntdict.keys()) != 1:
            samfile.close()
            return pileupcolumn.pos, ntdict
    #print (ntdictarr)
    samfile.close()
def popstats(bamdirpath, outpath): #taking bam directory path add output path
    colordict = cf.color(os.path.join(filedir,"DATA", "harrington_clinical_data.txt"),1) #grabbing the color
    with open(outpath,"w+") as fout: #opening output path as STDOUT in write mode and if it doesnt exist we are going to create it
        for root, dirs, files in os.walk(bamdirpath):
            print(files)
            for f in files: 
                if ".bam" not in f or ".bai" in f:
                    continue 
                print(f)
                pos, dic = pileup(os.path.join(root,f))
                name = f.split('.')[0]
                fout.write(f"{name}\t{colordict[name]}\t{pos}\t{dic}\n") #writing to output, name, splittng by . and taking the zero position aka the  name

def compilereport(inpath, outpath):
    with open(inpath, "r") as fin:
        with open(outpath,"w+") as fout:
            fout.write("The black mold was caused by a mutation in position 13. The wildtype base was C and the mutation was G\n")
            fout.write("The green mold was caused by a mutation in position 134. The wildtype base was C and the mutation was G\n")
            fout.write("The yellow mold was caused by a mutation in position 120. The wildtype base was G and the mutation was C\n")
            fout.write("The orange mold was caused by a mutation in position 83. The wildtype base was G and the mutation was C\n\n\n") 
            for line in fin:
                linearr = line.replace("\n","").split("\t") #getting rid of line characters and splitting into array via tab characters
               
                bases = []
                reads = []
                for keys,values in ast.literal_eval(linearr[-1]).items():
                    bases.append(keys)
                    reads.append(int(values))
                totalreads = reads[0] + reads[1]
                percentage = round(reads[0] / totalreads *100)
                fout.write(f"Sample {linearr[0]} had a {linearr[1]} mold, {totalreads} reads, and had {percentage}% of the reads at position {linearr[2]} had the mutation {bases[0]}.\n")



if __name__=="__main__":
    pileup()