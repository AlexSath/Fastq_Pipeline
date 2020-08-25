import subprocess
import os 
filePath = os.path.abspath(__file__)
filedir = os.path.dirname(filePath)

def createSams():
	referencepath = os.path.join(filedir, "DATA", "dgorgon_reference.fa")
	subprocess.call(["bwa", "index", referencepath])


	samdirpath = os.path.join(filedir, "samfiles")
	fastqdir = os.path.join(filedir, "fastqs") 
	for root,dirs,files in os.walk(fastqdir): 
		for f in files:
			if ".fastq" not in f: #checks to make sure the file we are using is a .fastq file
				continue
			filepath = os.path.join(root,f)
			samname = f.split("_")[0]+".sam" 
			sampath = os.path.join(samdirpath, samname)
			fout = open(sampath, "w+")
			subprocess.call(["bwa", "mem", referencepath, filepath],stdout = fout)

def processBams():
	samdirpath = os.path.join(filedir, "samfiles")
	unsortedBAM = os.path.join(filedir, "unsortedBAM")
	sortedBAM = os.path.join(filedir, "sortedBAM")
	for root,dirs,files in os.walk(samdirpath):
		for f in files:
			if ".sam" not in f:
				continue
			name = f.split(".")[0]
			sampath = os.path.join(root,f)
			unsortedBAMpath = os.path.join(unsortedBAM,f"{name}.bam")
			sortedBAMpath = os.path.join(sortedBAM,f"{name}.sorted.bam")

			#samtools view -bS {name}.sam > {name}.bam
			fout = open(unsortedBAMpath,"w+")
			subprocess.call(["samtools", "view", "-bS", sampath],stdout = fout)
			#samtools sort -m 100M {name}.bam > {name}.sorted.bam
			fout = open(sortedBAMpath, "w+")
			subprocess.call(["samtools", "sort", "-m", "100M", unsortedBAMpath],stdout = fout)
			#samtools index {name}.sorted.bam 
			subprocess.call(["samtools", "index", sortedBAMpath])

