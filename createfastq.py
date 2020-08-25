import argparse
import os

#defining variables filePath and filedir while using .os command for absolutle path while interfacing operating system
filePath = os.path.abspath(__file__) 
filedir = os.path.dirname(filePath)
##############################################
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)


def getbarcode(filepath, index): #defining a function for grabbing the barcode from the clinical data
    barcodekey = {} #creates and empty dictionary for the barcodes
    with open(filepath,"r") as fin: #opens the filepath in read mode as fin and reads every line
        for line in fin:
            linearr = line.split("\t") #line array takes the lines and splits by tab
            barcodekey[linearr[index].replace("\n","")] = linearr[0]
    return barcodekey

#this next set of lines basically does the same thing as the previous lines of code, but instead of the barcodes we are grapping the color of the mold. 
def color(filepath, index):
    barcodekey = {}
    with open(filepath,"r") as fin:
        for line in fin:
            linearr = line.split("\t")
            barcodekey [linearr[0]]= linearr[index].replace("\n","") #this line reads the opposite from the corresponding code above because of its position in the clinical data files. 
    return barcodekey



def createfastqs(fastqpath, barcodepath):
    hawkins_seq = ParseFastQ(fastqpath)
    barcodedict = getbarcode(barcodepath,2) 
    for root, dirs, files in os.walk(os.path.join(filedir,"fastqs")):
        for f in files:
            os.remove(os.path.join(root, f))
    for x in range(100000):
        try:
            nextdata = hawkins_seq.next()
        except:
            break
        

        tracker = 0
        prevchar = ""
        for char in nextdata[3]:
            if (prevchar == "D" or prevchar == "F")\
                and (char =="D" or char == "F"):
                break
            prevchar = char
            tracker +=1 

        
        filename = f"{barcodedict[nextdata[1][0:5]]}_trimmed.fastq"
        filepath = os.path.join(filedir, "fastqs", filename)

        with open (filepath, "a+") as fout:
            fout.write(nextdata[0]+"\n")
            fout.write(nextdata[1][5:tracker-1]+"\n")
            fout.write(nextdata[2]+"\n")
            fout.write(nextdata[3][5:tracker-1]+"\n")


    print("done")





