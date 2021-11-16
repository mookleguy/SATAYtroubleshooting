###GRAVEYARD***
#pergenefiles<-list.files(pattern="pergene.txt$")
#contset<-(1:length(voldatatn))[!(1:length(voldatatn)) %in% testset]
#testName<-gsub(".NlaIII_S.{1,2}_L001_R1_001.fastq_pergene.txt","",testName)
#setname<-paste0(testName,"_vs_all")
#datasets <-names(voldatatn[testset])
#names(voldataread[testset])
###GRAVEYARD###



workpath<-"E:/satay/allPerGenes"
setwd(workpath)

# The fold change is determined by the ratio between the reference and the experimental dataset.
# When one of the datasets is 0, this is false results for the fold change.
# To prevent this, the genes with 0 insertions are set to have 5 insertions, and the genes with 0 reads are set to have 25 reads.
# These values are determined in dicussion with the Kornmann lab.
noisetn=5 
noiseread=25

#define list variable for transposons and reads
voldatatn<-list()
voldataread<-list()



#iterate through directory and read pergene.txt files
#normalize each genes fraction of transposons by the total number of transposons for that data set. Add 5 to each gene in denominator
#normalize each genes fraction of reads by the total number of reads for that data set. Add 25 to each gene in denominator
for (i in list.files(pattern="pergene.txt$")) {
  voldata<-read.table(i, header = TRUE)
  voldata$number_of_transposon_per_gene<-(voldata$number_of_transposon_per_gene+noisetn)/sum(voldata$number_of_transposon_per_gene+noisetn)
  voldata$number_of_read_per_gene<-(voldata$number_of_read_per_gene+noiseread)/sum(voldata$number_of_read_per_gene+noiseread)
  voldatatn[[i]]<-voldata$number_of_transposon_per_gene
  voldataread[[i]]<-voldata$number_of_read_per_gene
}



#make tns and reads data frame type
voldatatn<-data.frame(voldatatn)
voldataread<-data.frame(voldataread)

#definre your testset
testset<-c(2)
#define control set
contset<-c(3)

#variable to name file. Original script did this poorly. 
testName<-paste(names(voldatatn[testset]), "_vs_ ",collapse = '')
contName <- paste(names(voldatatn[contset]), collapse = '')
setname<-paste0(testName, contName)

#want another placeholder for this variable so we can mess with the two voldata variables
dtn<-voldatatn
drd<-voldataread

#if more than 1 control set chosen, average the reads and transposons
if(length(contset)>1){
  dtnsum<-rowSums(dtn[,contset])/length(contset)
  drdsum<-rowSums(drd[,contset])/length(contset)
} else {
  dtnsum<-dtn[,contset]/length(contset)
  drdsum<-drd[,contset]/length(contset)
}

#if more than 1 test set chosen, average the reads and transposons
#compute the fold change of tns and reads compared to cont
if(length(testset)>1){
  volcanotn<-rowSums(dtn[,testset])/(dtnsum*length(testset))
  volcanoread<-rowSums(drd[,testset])/(drdsum*length(testset))
} else {
  volcanotn<-(dtn[,testset])/(dtnsum*length(testset))
  volcanoread<-(drd[,testset])/(drdsum*length(testset))
}

#more variable duplication
studenttn<-volcanotn
studentread<-volcanotn


pb<-txtProgressBar(min=0, max = length(dtn[[1]]))
if(length(testset)>1){
  for (j in 1:length(dtn[[1]])) {
    tt<-t.test(as.numeric(dtn[j,contset]), as.numeric(dtn[j,testset] ), var.equal = TRUE)
    studenttn[j]<-tt$p.value
    tt<-t.test(as.numeric(drd[j,contset]), as.numeric(drd[j,testset] ), var.equal = TRUE)
    studentread[j]<-tt$p.value
    setTxtProgressBar(pb,j)
  }
}else{
  # for (j in 1:length(dtn[[1]])) {
  #   tt<-t.test(dtn[j,contset], mu=dtn[j,testset] )
  #   #print (dtn[j,testset])
  #   studenttn[j]<-tt$p.value
  #   tt<-t.test(drd[j,contset], mu=drd[j,testset] )
  #   studentread[j]<-tt$p.value
  #   setTxtProgressBar(pb,j)
  # }
  for (j in 1:length(dtn[[1]])) {
    tt<-t.test(dtn[j,contset], dtn[j,testset] )
    #print (dtn[j,testset])
    studenttn[j]<-tt$p.value
    tt<-t.test(drd[j,contset], drd[j,testset] )
    studentread[j]<-tt$p.value
    setTxtProgressBar(pb,j)
  }
}

library(plotly)
df<-list()
df$xx<-log2(volcanotn)
df$yy<--log(studenttn) #+0.00000000000001)
df$name<-voldata$gene_name
#df$pathway<-factor(voldata$pathway)
df<-data.frame(df)

p<-plot_ly(df,
           x=~xx,
           y=~yy,
           text=~name,
           #color=~pathway,
           colors = "Set1"
)
htmlwidgets::saveWidget(p,paste0(workpath,setname,"TN.html"))

print(p)

write.csv(df, file = paste0(workpath,setname,"TN.csv"))

library(plotly)
df<-list()
df$xx<-log2(volcanoread)
df$yy<--log(studentread)
df$name<-voldata$gene_name
df<-data.frame(df)

p<-plot_ly(df,
           x=~xx,
           y=~yy,
           text=~name
)
htmlwidgets::saveWidget(p,paste0(workpath,setname,"read.html"))

print(p)

write.csv(df, file = paste0(workpath,setname,"read.csv"))