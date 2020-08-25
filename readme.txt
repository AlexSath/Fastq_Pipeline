1) For the pipeline.py script to work please download Assignment3 directory within home/pated/week6
2) Within Assignment3 absolute paths were created within "DATA" for harrington_clinical_data.txt, hawkins_pooled_sequences.fasta and dgorgon_reference.fa
3) Additional requirements: installation of samtools, pysam, python3, Burrows-Wheeler Aligner tool
4)Pipeline.py is the master script and broken up into createfastq.py which gives the trimmed fasta files, part23.py with will incorporate the BWA tool while also converting the sam files to bam files which are then sorted and part45.py which will produce our variant discovery and out final output report. 


P.S. for the first section of part 5, I was not able to figure out how to code to output what nucleotide position and mutation is responsible for each color of the mold. I think I am overthinking into statistics but cant wrap my head around how to get this part done. Nonetheless when looking out the output created in part 4 (in this case tempout.csv), it seemed like there were only 4 different iterations, thus I used fout and just wrote them into the report. 
