import re
import glob

#Glob all desired files from specified location. 
fileNames=[]
for fileName in glob.glob('../Allele_Filter/*filterAF.txt'):
	fileNames.append(fileName)

#Regex to save sample ID
def captureName(x):
	search=re.search("([0-9]+)",x)
	capture=search.group(1)
	return capture

#Initialize the dictionary that is going to be storing all unique variants in the entire cohort. 
#The count for any variant is increased by one any time that same variant is seen again in a different sample. 
#The dictionary is stored in memory until the last for loop of this script where each variant present only once or twice in the cohort is saved and all others are deleted from each sample.  
dict={}


#This first for loop ensures that variants don't get double counted in each individual sample due to multiple transcript. 
#This is done using the keyDict dictionary which is re-set for each new sample. 
for x in fileNames:
	openFile=open(x,'r')
	keyDict={}
	for line in openFile:
		if not line.startswith('Sample'):
			chr=line.split('\t')[1]
			pos=line.split('\t')[3]
			ref=line.split('\t')[4]
			alt=line.split('\t')[5]
			gene=line.split('\t')[6]
			key=chr+','+pos+','+ref+','+alt+','+gene
			if key not in keyDict:
				keyDict[key]=1
				if key in dict:
					dict[key]+=1
				else:
					dict[key]=1
	openFile.close()


for key in dict:
	val=dict[key]
	print key,val

for x in fileNames:
	openAgain=open(x,'r')
	id=captureName(x)
	newName=id+'filter.txt'
	newFile=open(newName,'w')
	for line in openAgain:
		if not line.startswith('Sample'):
			chr=line.split('\t')[1]
	               	pos=line.split('\t')[3]
                	ref=line.split('\t')[4]
                	alt=line.split('\t')[5]
                	gene=line.split('\t')[6]
			key=chr+','+pos+','+ref+','+alt+','+gene
			value=dict[key]
			if (value <= 2):
				newFile.write(line)	

	openAgain.close()
	newFile.close()
