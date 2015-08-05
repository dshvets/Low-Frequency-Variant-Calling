#This code was written as we aimed to solve our quality problem related to high strand bias in certain capture kits. 
#Our cohort under analysis used three different capture kits and it was noted that in variants where the base was changing from A to C (or T to G), the strand bias was exceptionally high.
#To solve which capture kit had caused this error, this script creates a plot of every possible base change, with strand bias on the x-axis and the number of variants seen on the y-axis. 
#Change the file name in the first line of this script for read.table in order to change the input file. 
#The format of the input file is a tab delimited file with the first column being the base change (ex: A,C), the second column is strand bias, and the third column is count. 

indexTable <- read.table("~/Desktop/NIH-CH/truseq.txt",stringsAsFactors=FALSE,header=TRUE)
library(plyr)
library(reshape)
library(ggplot2)

acIndexTable <-data.frame()
agIndexTable <-data.frame()
atIndexTable <-data.frame()
caIndexTable <-data.frame()
cgIndexTable <-data.frame()
ctIndexTable <-data.frame()
gaIndexTable <-data.frame()
gcIndexTable <-data.frame()
gtIndexTable <-data.frame()
taIndexTable <-data.frame()
tcIndexTable <-data.frame()
tgIndexTable <-data.frame()

for(i in 1:length(indexTable$Base.Change)){
  if(indexTable$Base.Change[i] == 'A,C'){
    acIndexTable <- rbind(acIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='A,G' ){
    agIndexTable <- rbind(agIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='A,T'){
    atIndexTable <- rbind(atIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='C,A'){
    caIndexTable <- rbind(caIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='C,G'){
    cgIndexTable <- rbind(cgIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='C,T'){
    ctIndexTable <- rbind(ctIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='G,A'){
    gaIndexTable <- rbind(gaIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='G,C'){
    gcIndexTable <- rbind(gcIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='G,T'){
    gtIndexTable <- rbind(gtIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='T,A'){
    taIndexTable <- rbind(taIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='T,C'){
    tcIndexTable <- rbind(tcIndexTable,indexTable[i,])
  }else if(indexTable$Base.Change[i] =='T,G'){
    tgIndexTable <- rbind(tgIndexTable,indexTable[i,])
  }
}

row.names(acIndexTable) <- NULL
row.names(agIndexTable) <- NULL
row.names(atIndexTable) <- NULL
row.names(caIndexTable) <- NULL
row.names(cgIndexTable) <- NULL
row.names(ctIndexTable) <- NULL
row.names(gaIndexTable) <- NULL
row.names(gcIndexTable) <- NULL
row.names(gtIndexTable) <- NULL
row.names(taIndexTable) <- NULL
row.names(tcIndexTable) <- NULL
row.names(tgIndexTable) <- NULL

acIndexTable$Base.Change <- NULL
agIndexTable$Base.Change <- NULL
atIndexTable$Base.Change <- NULL
caIndexTable$Base.Change <- NULL
cgIndexTable$Base.Change <- NULL
ctIndexTable$Base.Change <- NULL
gaIndexTable$Base.Change <- NULL
gcIndexTable$Base.Change <- NULL
gtIndexTable$Base.Change <- NULL
taIndexTable$Base.Change <- NULL
tcIndexTable$Base.Change <- NULL
tgIndexTable$Base.Change <- NULL

acIndexTable <-acIndexTable[!(acIndexTable$Strand.Bias >= 75),]
tgIndexTable <-tgIndexTable[!(tgIndexTable$Strand.Bias >= 75),]
agIndexTable <-agIndexTable[!(agIndexTable$Strand.Bias >= 75),]
atIndexTable <-atIndexTable[!(atIndexTable$Strand.Bias >= 75),]
tcIndexTable <-tcIndexTable[!(tcIndexTable$Strand.Bias >= 75),]
taIndexTable <-taIndexTable[!(taIndexTable$Strand.Bias >= 75),]
caIndexTable <-caIndexTable[!(caIndexTable$Strand.Bias >= 75),]
cgIndexTable <-cgIndexTable[!(cgIndexTable$Strand.Bias >= 75),]
ctIndexTable <-ctIndexTable[!(ctIndexTable$Strand.Bias >= 75),]
gtIndexTable <-gtIndexTable[!(gtIndexTable$Strand.Bias >= 75),]
gcIndexTable <-gcIndexTable[!(gcIndexTable$Strand.Bias >= 75),]
gaIndexTable <-gaIndexTable[!(gaIndexTable$Strand.Bias >= 75),]

#AC and TG pair

for (i in 1:length(acIndexTable$Strand.Bias)){
  if (is.element(acIndexTable$Strand.Bias[i], tgIndexTable$Strand.Bias)==FALSE){
    newrow <-c(acIndexTable$Strand.Bias[i],0)
    tgIndexTable <- rbind(tgIndexTable,newrow)
  }
}

for (i in 1:length(tgIndexTable$Strand.Bias)){
  if (is.element(tgIndexTable$Strand.Bias[i], acIndexTable$Strand.Bias)==FALSE){
    newrow <-c(tgIndexTable$Strand.Bias[i],0)
    acIndexTable <- rbind(acIndexTable,newrow)
  }
}


#AG and TC pair

for (i in 1:length(agIndexTable$Strand.Bias)){
  if (is.element(agIndexTable$Strand.Bias[i], tcIndexTable$Strand.Bias)==FALSE){
    newrow <-c(agIndexTable$Strand.Bias[i],0)
    tcIndexTable <- rbind(tcIndexTable,newrow)
  }
}

for (i in 1:length(tcIndexTable$Strand.Bias)){
  if (is.element(tcIndexTable$Strand.Bias[i], agIndexTable$Strand.Bias)==FALSE){
    newrow <-c(tcIndexTable$Strand.Bias[i],0)
    agIndexTable <- rbind(agIndexTable,newrow)
  }
}

#AT and TA pair


for (i in 1:length(atIndexTable$Strand.Bias)){
  if (is.element(atIndexTable$Strand.Bias[i], taIndexTable$Strand.Bias)==FALSE){
    newrow <-c(atIndexTable$Strand.Bias[i],0)
    taIndexTable <- rbind(taIndexTable,newrow)
  }
}

for (i in 1:length(taIndexTable$Strand.Bias)){
  if (is.element(taIndexTable$Strand.Bias[i], atIndexTable$Strand.Bias)==FALSE){
    newrow <-c(taIndexTable$Strand.Bias[i],0)
    atIndexTable <- rbind(atIndexTable,newrow)
  }
}

#CA and GT pair

for (i in 1:length(caIndexTable$Strand.Bias)){
  if (is.element(caIndexTable$Strand.Bias[i], gtIndexTable$Strand.Bias)==FALSE){
    newrow <-c(caIndexTable$Strand.Bias[i],0)
    gtIndexTable <- rbind(gtIndexTable,newrow)
  }
}

for (i in 1:length(gtIndexTable$Strand.Bias)){
  if (is.element(gtIndexTable$Strand.Bias[i], caIndexTable$Strand.Bias)==FALSE){
    newrow <-c(gtIndexTable$Strand.Bias[i],0)
    caIndexTable <- rbind(caIndexTable,newrow)
  }
}

#CG and GC pair

for (i in 1:length(cgIndexTable$Strand.Bias)){
  if (is.element(cgIndexTable$Strand.Bias[i], gcIndexTable$Strand.Bias)==FALSE){
    newrow <-c(cgIndexTable$Strand.Bias[i],0)
    gcIndexTable <- rbind(gcIndexTable,newrow)
  }
}

for (i in 1:length(gcIndexTable$Strand.Bias)){
  if (is.element(gcIndexTable$Strand.Bias[i], cgIndexTable$Strand.Bias)==FALSE){
    newrow <-c(gcIndexTable$Strand.Bias[i],0)
    cgIndexTable <- rbind(cgIndexTable,newrow)
  }
}

#CT and GA pair

for (i in 1:length(gaIndexTable$Strand.Bias)){
  if (is.element(gaIndexTable$Strand.Bias[i], ctIndexTable$Strand.Bias)==FALSE){
    newrow <-c(gaIndexTable$Strand.Bias[i],0)
    ctIndexTable <- rbind(ctIndexTable,newrow)
  }
}

for (i in 1:length(ctIndexTable$Strand.Bias)){
  if (is.element(ctIndexTable$Strand.Bias[i], gaIndexTable$Strand.Bias)==FALSE){
    newrow <-c(ctIndexTable$Strand.Bias[i],0)
    gaIndexTable <- rbind(gaIndexTable,newrow)
  }
}

#Order all newly created tables by Strand Bias column in numerically increasing order
acIndexTable <-arrange(acIndexTable, Strand.Bias)
agIndexTable <-arrange(agIndexTable,Strand.Bias)
atIndexTable <-arrange(atIndexTable,Strand.Bias)
caIndexTable <-arrange(caIndexTable,Strand.Bias)
cgIndexTable <-arrange(cgIndexTable,Strand.Bias)
ctIndexTable <-arrange(ctIndexTable,Strand.Bias)
gaIndexTable <-arrange(gaIndexTable,Strand.Bias)
gcIndexTable <-arrange(gcIndexTable,Strand.Bias)
gtIndexTable <-arrange(gtIndexTable,Strand.Bias)
taIndexTable <-arrange(taIndexTable,Strand.Bias)
tcIndexTable <-arrange(tcIndexTable,Strand.Bias)
tgIndexTable <-arrange(tgIndexTable,Strand.Bias)


plotAC.TG <- data.frame(x=acIndexTable$Strand.Bias, y1=acIndexTable$Count, y2=tgIndexTable$Count)
melt_AC.TG <- melt(plotAC.TG,id='x')
graph_AC.TG <- ggplot(melt_AC.TG,aes(x=x,y=value,fill=variable))+geom_bar(stat='identity',position='identity',alpha=.7)+xlab("Strand Bias")+ylab("Count")+ggtitle("AC and TG")+scale_fill_discrete(name="Base Change",labels=c("AC","TG"))

plotAG.TC <- data.frame(x=agIndexTable$Strand.Bias, y1=agIndexTable$Count, y2=tcIndexTable$Count)
melt_AG.TC <- melt(plotAG.TC,id='x')
graph_AG.TC <- ggplot(melt_AG.TC,aes(x=x,y=value,fill=variable))+geom_bar(stat='identity',position='identity',alpha=.7)+xlab("Strand Bias")+ylab("Count")+ggtitle("AG and TC")+scale_fill_discrete(name="Base Change",labels=c("AG","TC"))

plotAT.TA <- data.frame(x=atIndexTable$Strand.Bias, y1=atIndexTable$Count, y2=taIndexTable$Count)
melt_AT.TA <- melt(plotAT.TA,id='x')
graph_AT.TA <- ggplot(melt_AT.TA,aes(x=x,y=value,fill=variable))+geom_bar(stat='identity',position='identity',alpha=.7)+xlab("Strand Bias")+ylab("Count")+ggtitle("AT and TA")+scale_fill_discrete(name="Base Change",labels=c("AT","TA"))


plotCA.GT <- data.frame(x=caIndexTable$Strand.Bias, y1=caIndexTable$Count, y2=gtIndexTable$Count)
melt_CA.GT <- melt(plotCA.GT,id='x')
graph_CA.GT <- ggplot(melt_CA.GT,aes(x=x,y=value,fill=variable))+geom_bar(stat='identity',position='identity',alpha=.7)+xlab("Strand Bias")+ylab("Count")+ggtitle("C->A and G->T")+scale_fill_discrete(name="Base Change",labels=c("C->A","G->T"))


plotCG.GC <- data.frame(x=cgIndexTable$Strand.Bias, y1=cgIndexTable$Count, y2=gcIndexTable$Count)
melt_CG.GC <- melt(plotCG.GC,id='x')
graph_CG.GC <- ggplot(melt_CG.GC,aes(x=x,y=value,fill=variable))+geom_bar(stat='identity',position='identity',alpha=.7)+xlab("Strand Bias")+ylab("Count")+ggtitle("C->G and G->C")+scale_fill_discrete(name="Base Change",labels=c("C->G","G->C"))

plotCT.GA <- data.frame(x=ctIndexTable$Strand.Bias, y1=ctIndexTable$Count, y2=gaIndexTable$Count)
melt_CT.GA <- melt(plotCT.GA,id='x')
graph_CT.GA <- ggplot(melt_CT.GA,aes(x=x,y=value,fill=variable))+geom_bar(stat='identity',position='identity',alpha=.7)+xlab("Strand Bias")+ylab("Count")+ggtitle("C->T and G->A")+scale_fill_discrete(name="Base Change",labels=c("C->T","G->A"))

graph_AC.TG
graph_AG.TC
graph_AT.TA
graph_CA.GT
graph_CG.GC
graph_CT.GA
multiplot(graph_AC.TG,graph_AG.TC,graph_AT.TA,graph_CA.GT,graph_CG.GC,graph_CT.GA,cols=2)



