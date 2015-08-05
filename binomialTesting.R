fullName <- commandArgs(trailingOnly=TRUE)

id <- gsub("[A-Za-z/_.]+","",fullName)

data<-read.table(fullName,header=TRUE,sep="\t")

#Unique column names from Original sample table
uniqueCols <- c("Chromosome","Position","Ref.Base","GeneName","ReadDepth","Forward.NonRef","Reverse.NonRef")

#Get the unique data necessary for binomial test and BH correction
uniqueData <- unique(data[,uniqueCols])

#Create the function called test that performs the binomial test and extracts the p-value
test <- function(k,n){binom.test(k,n,p=0.5,alternative="less")$p.value}

#Apply the binomial test to each unique row of data and save results to new column BinomialPvalue
uniqueData$BinomialPvalue <- mapply(test,uniqueData$Forward.NonRef+uniqueData$Reverse.NonRef,uniqueData$ReadDepth)

#Save the p-values to an array in order to perform BH correction
pValsArray <- uniqueData$BinomialPvalue

#Perform the BH correction and save results to new column called BHcorrection
uniqueData$BHcorrection <- p.adjust(pValsArray,method="BH")

newName <- paste(id,"pvals.txt",sep="")

#Remove quotes and row name count that R automatically populates
write.table(uniqueData,newName,sep="\t",row.names=FALSE,quote=FALSE)




