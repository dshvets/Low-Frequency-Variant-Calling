import subprocess
import sys
import re

#Receive command line input of the sample ID
annovarName=sys.argv[1]

def captureName(x):
        sampleSearch = re.search("([0-9]+)",x)
        sampleCapture=sampleSearch.group(1)
        return sampleCapture

id=captureName(annovarName)

#Provide names for the final file to be created, the file containing the binomial and BH correction pvalues(somatic.txt) and the R script
finalName=id+"finalBinomialTESTING.txt"
somatic=id+"pvals.txt"
BinomScript="binomialTesting.R"

#Run the R script using the annotated file, the output is saved to a file called idsomatic.txt. Using subprocess.call is considered "blocking" so this script should wait for the Rscript to finish executing before continuing further down the code. 
runBinom =subprocess.call(["Rscript",BinomScript,annovarName])

#Create empty directory for pushing the data from the binomial results
dict={}

#Open the binomial results file and push data into dictionary
somaticFile=open(somatic,'r')
for line2 in somaticFile:
                chr2=line2.split()[0]
                pos2=line2.split()[1]
                base2=line2.split()[2]
                gene2=line2.split()[3]
		pVal=line2.split()[7]
		bhVal=line2.split()[8]
		key = chr2+","+pos2+","+base2+","+gene2
		dict[key]=pVal +","+ bhVal		

#Close results file and open original annotated file as well as new final file for writing that will include all annotation information as well as the p- and q-values
somaticFile.close()
annovarFile=open(annovarName,'r')
finalBinomial=open(finalName,'w')
finalBinomial.write("Sample\tChromosome\tFrameShift\tPosition\tRef.Base\tAlt.Base\tGeneName\tTranscript\tcDNA\tAminoAcid\tLofreqQuality\tReadDepth\tAlleleFrequency\tStrandBias\tForward.Ref\tReverse.Ref\tForward.NonRef\tReverse.NonRef\tHRUN(if indel)\tBinomialPvalue\tBHCorrection\n")

#For each line in annotated file we use four identifiers(chromosome,position,base,and gene) as the key to retrieve the proper p- and q-values and write these to the final file
for line1 in annovarFile:
	if not line1.startswith("Sample"):
		chr1=line1.split()[1]
		pos1=line1.split()[4]
		base1=line1.split()[5]
		gene1=line1.split()[7]
		newKey = chr1+","+pos1+","+base1+","+gene1
		data = dict[newKey]
		data=data.split(",")
		pValue=data[0]
		qValue=data[1]
		line1=line1.rstrip()
		finalLine=line1+'\t'+data[0]+'\t'+data[1]
		finalBinomial.write(finalLine+'\n')

annovarFile.close()
finalBinomial.close()

