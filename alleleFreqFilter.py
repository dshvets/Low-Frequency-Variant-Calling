import sys
import re

inputName=sys.argv[1]

def captureName(x):
        sampleSearch = re.search("([0-9]+)",x)
        sampleCapture=sampleSearch.group(1)
        return sampleCapture

id=captureName(inputName)

newName=id+'filterAF.txt'

openFile=open(inputName,'r')
newFile=open(newName,'w')

for line in openFile:
	if not line.startswith("Sample"):
		af=line.split('\t')[12]
		af=float(af)
		if (af > 0.10):
			newFile.write(line)

openFile.close()
newFile.close()
