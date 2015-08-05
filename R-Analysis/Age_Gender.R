#This script creates a histogram of a set of data containing age and gender information.
#The histogram separates the two genders and plots age on the x-axis, with the number of people in the cohort having that age on the y-axis. 

librar(plyr)

age <- read.table("~/Desktop/NIH-CH/AgeAtEnrollment.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)
gender <- read.table("~/Desktop/NIH-CH/gender.txt",stringsAsFactors=FALSE,sep="\t",header=TRUE)

genderSamples <- gender[,1]
newAge <- data.frame()

for(i in 1:length(age$Sample)){
  if(is.element(age$Sample[i],genderSamples)){
    newAge <- rbind(newAge,age[i,])
  }
}

row.names(newAge) <- NULL
newAge <- arrange(newAge,Sample)
gender <- arrange(gender,Sample)

maleAge <-vector()
femAge <- vector()
for(i in 1:length(newAge$Sample)){
  if(gender[i,2]=='M'){
    maleAge <- c(maleAge,newAge[i,2])
  }else{
    femAge <- c(femAge,newAge[i,2])
  }
}



df <- rbind(data.frame(fill="orange",obs=maleAge),data.frame(fill="green",obs=femAge))
plotHist <-ggplot(df,aes(x=obs,fill=fill))+geom_histogram(binwidth=1,position="dodge",alpha=0.8)+xlab("Age")+ylab("Number of Patients")+ggtitle("Clinseq Cohort Age & Gender")+theme(axis.title=element_text(size=15),plot.title=element_text(size=20,face='bold'))+scale_fill_discrete(name="Gender",labels=c("Male","Female"))

