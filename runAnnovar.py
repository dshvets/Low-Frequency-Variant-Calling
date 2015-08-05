import sys
import subprocess
import re

sample = sys.argv[1]

#Function for regular expression search for the ID name from the input .vcf file from the command line
def captureName(x):
	sampleSearch = re.search("([0-9]+)",x)
	sampleCapture=sampleSearch.group(1)
	return sampleCapture

#Function for regular expression search of gene name
def captureGene(x):
        search=re.search("(^[A-Za-z0-9\-\.]*):",x)
        capture=search.group(1)
        return capture

#Function for regular expression search for transcript names
def captureTranscript(x):
	search=re.findall(":([A-Za-z0-9\.]*):exon",x)
	capture=','.join(search)
	return capture

#Function for regular expression search of cDNA
def captureCDNA(x):
	search=re.search("(c\.[0-9]*[AGCT]>[AGCT]):p",x)
	capture=search.group(1)
	return capture

#Function for regular expression search for amino acid change
def captureAminoAcid(x):
	search=re.search("p\.([A-Z][0-9]*[A-Z])",x)
	capture=search.group(1)
	return capture

#Function for extracting read depth
def getDepth(x):
        sampleSearch=re.search("DP=([0-9]+);",x)
        capture=sampleSearch.group(1)
        return capture

#Function for extracting allele frequency
def getFreq(x):
        sampleSearch=re.search("AF=([0-9]\.[0-9]*);",x)
        capture=sampleSearch.group(1)
        return capture

#Function for extracting bias
def getBias(x):
        sampleSearch=re.search("SB=([0-9]*);",x)
        capture=sampleSearch.group(1)
        return capture

#Function for extracting DP4 values
def getFour(x):
        sampleSearch=re.search("DP4=([0-9]*,[0-9]*,[0-9]*,[0-9]*)",x)
        capture=sampleSearch.group(1)
        return capture


def indelTranscript(x):
        search=re.findall(":([A-Za-z0-9\.]*):exon",x)
        capture=','.join(search)
        return capture



def indelcDNA(x):
        search=re.search("(c\.[0-9]+_[0-9]+[a-zA-Z]+):",x)
        capture=search.group(1)
        return capture


def indelProt(x):
        search=re.search(":(p\.[0-9A-Za-z]+)",x)
        capture=search.group(1)
        return capture

def getHrun(x):
	search=re.search("HRUN=([0-9]+)",x)
	capture=search.group(1)
	return capture


#Setting to variables all of the needed paths for the files and scripts used when calling Annovar
database="/cluster/directoy/location/annovar_db"
annovar_script="/usr/directory/location/annotate_variation.pl"
conversion_script="/usr/directory/location/convert2annovar.pl"

#Subprocess call for converting the filtered vcf files that contain lofreq variants into the correct format for Annovar
id = captureName(sample)
fileName=id+".avinput"
outputFile=open(fileName,"w")
run=subprocess.Popen([conversion_script,"-format","vcf4", sample,"-includeinfo"],stdout=subprocess.PIPE)
result = run.stdout.read()
run.wait()
outputFile.write(result)
outputFile.close()


#Subprocess call for running the Annovar annotation using the hg19 database, this creates two resulting files in the directory with the extensions avinput.exonic_variant_function and avinput.variant_function. We only need to use the avinput.exonic_variant_function file because we only care about the variants in the exon regions rather than intron regions
subprocess.call([annovar_script, "-geneanno","-dbtype","knownGene","-hgvs","-buildver", "hg19",fileName, database])


#finalFile represents the final output file to which our desired data is going to be written to in the proper tab delimited format
finalFile = id+"final.annovar"
finalOutput=open(finalFile,"w")
#Writing the header to the final output file
finalOutput.write("Sample\tChromosome\tFrameShift\tPosition\tRef.Base\tAlt.Base\tGeneName\tTranscript\tcDNA\tAminoAcid\tLofreqQuality\tReadDepth\tAlleleFrequency\tStrandBias\tForward.Ref\tReverse.Ref\tForward.NonRef\tReverse.NonRef\tHRUN(if indel)\n")
#exonFile is the file created by running Annovar, it contains the desired amino acid change and annotation information(avinput.exonic_variant_function). It is being opened here to extract the needed information from it and write it to the final file in the correct format for later analysis in R
exonFile = fileName+".exonic_variant_function"
readFile=open(exonFile,"r")
#This loop goes through every line in the file containing exonic variant annotation and extracts the needed data from each line and writes it in the correct format to the final output file. It also skips lines that return an UNKNOWN gene name from annovar
for line in readFile:
	try:
		frameShift=line.split()[1]
		frameRest=line.split()[2]
		geneData = line.split()[3]
		chr=line.split()[4]
		position=line.split()[10]
		refBase=line.split()[12]
		altBase=line.split()[13]
		lofreqQual=line.split()[14]
		vcfInfo=line.split()[16]
		depth=getDepth(vcfInfo)
		freq=getFreq(vcfInfo)
		bias=getBias(vcfInfo)
		dp4=getFour(vcfInfo)
		dp4=dp4.split(',')
		gene=captureGene(geneData)
		geneArray=geneData.split(',')
		geneArray.pop()
		totFrameShift=frameShift+' '+frameRest
		for x in geneArray:
			if totFrameShift == "frameshift insertion" or totFrameShift=="nonframeshift insertion" or totFrameShift=="frameshift deletion" or totFrameShift=="nonframeshift deletion":
				trans=indelTranscript(x)
				indeldna=indelcDNA(x)	
				indelaa=indelProt(x) 
				hrun=getHrun(vcfInfo)
				newLine=id,chr,totFrameShift,position,refBase,altBase,gene,trans,indeldna,indelaa,lofreqQual,depth,freq,bias,dp4[0],dp4[1],dp4[2],dp4[3],hrun
				finalOutput.write('\t'.join(map(str,newLine))+'\n')
			else:
				transcript=captureTranscript(x)
				basechange =captureCDNA(x)
				aminoacid=captureAminoAcid(x)
				newLine=id,chr,totFrameShift,position,refBase,altBase,gene,transcript,basechange,aminoacid,lofreqQual,depth,freq,bias,dp4[0],dp4[1],dp4[2],dp4[3],0
				finalOutput.write('\t'.join(map(str,newLine))+'\n')
	except:
		if "UNKNOWN" in geneData:
			continue

finalOutput.close()
readFile.close()






