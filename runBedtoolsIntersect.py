#!/opt/sw/software/Python/2.7.3-GCC-4.1.2/bin/python
import subprocess
import glob
import re
import sys

sample = sys.argv[1]

#Function for regular expression match of the sample ID 
def captureName(x):
	sampleSearch = re.search("([0-9]+)",x)
	sampleCapture = sampleSearch.group(1)
	return sampleCapture


#Obtaining all BED files from directory where they are located. All six of these files are necessary for filtering out all of the less desirable regions before conducting further analysis. 
filterFiles =[]
for x in glob.glob('/directory/location/Genovese_BedFiles/*.bed'):
       filterFiles.append(x)


#The below code recreates the command line "pipe" using subprocess, each previous subprocess run is saved as "stdin" and gets fed into the following subprocess run. The final bedtools intersect filter does not have a "-v" because it is the Strict Mask filter which requires the opposite of the others. 
x = captureName(sample)	
fileName = x + "BedFilter.vcf"
outputFile = open(fileName,"w")
run1 = subprocess.Popen(["bedtools","intersect","-a",sample,"-b", filterFiles[0], "-v"], stdout=subprocess.PIPE)
run2 = subprocess.Popen(["bedtools","intersect","-a","stdin","-b", filterFiles[1], "-v"],stdin=run1.stdout,stdout=subprocess.PIPE)
run3 = subprocess.Popen(["bedtools","intersect","-a","stdin","-b", filterFiles[2],"-v"],stdin=run2.stdout,stdout = subprocess.PIPE)
run4 = subprocess.Popen(["bedtools","intersect","-a","stdin","-b", filterFiles[3],"-v"],stdin=run3.stdout, stdout= subprocess.PIPE)
run5 = subprocess.Popen(["bedtools","intersect","-a","stdin","-b", filterFiles[5],"-v"],stdin=run4.stdout, stdout= subprocess.PIPE)
run6 = subprocess.Popen(["bedtools","intersect","-a","stdin","-b", filterFiles[4]],stdin=run5.stdout, stdout= outputFile)
run6.wait()
outputFile.close()



 
