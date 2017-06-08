### This code includes functions for managing de datasets from the hierarchical models

library(Hmisc)

##spshort creates an abbreviation of six letters for each species
spshort=function()
{
   fulldata=read.delim(file=nourage)
   spnames=sort(unique(fulldata$sp))
   spshort=character()
   counter=1
   for (i in 1:length(spnames))  
   {
     spc<-gsub(" ","-",spnames[i])
     spc1<-regexpr("-",spc)
     genuscode<-substr(spc,1,3)
     spcode<-substr(spc,as.numeric(spc1[1])+1,as.numeric(spc1[1])+3)
     spshort[counter]=paste(genuscode,spcode,sep="")
    counter=counter+1
   } 
 return=spshort
}   
   

datasetnrg<-function(dataset=nourage)
{
nrg<-read.delim(dataset)
spnames=sort(unique(nrg$sp))
for (i in 1:length(spnames))
 {
  spdata=extract.seedfall_onesp(latin = spnames[i], file = nourage)
  spc<-gsub(" ","-",spnames[i])
  spc1<-regexpr("-",spc)
  genuscode<-substr(spc,1,3)
  spcode<-substr(spc,as.numeric(spc1[1])+1,as.numeric(spc1[1])+3)
  spshort=paste(genuscode,spcode,sep="")
  spshort=spdata
  write.table(spdata,file=paste(spshort,"txt",sep="."),row.names=F,sep="\t") 
 }

}

initial.peak=function(dataset)
{
 nyear=unique(dataset$year)
 peaksyr<-vector()
 for (i in 1:length(nyear))
  {
    oneyr<-subset(dataset, year==nyear[i])
    peaksyr[i]= round(sum(oneyr$quantity, na.rm=T),0)
  }
 return (peaksyr)
}



## do a loop for all species in Nouragues

 hierarchical.NRG<-function(file=nourage)
 
{
   fulldata=read.delim(file)
   spnames=sort(unique(fulldata$sp))
   results=list()
   spshort=spshort()
   startv=vector()
   for (i in 1:length(spnames))

  {  
     spdata<-extract.seedfall_onesp(latin =spnames[i], file = file)
     maxdate=which.max(spdata$quantity)
     startv[1]=spdata$yday[maxdate]
     startv[2]= 10
     startv[3]=round(log(spdata$quantity[maxdate]),1)
     startv[4]= 1.2
     startv[5]= 25
     fit=model.seedfall.Gibbs(data=spdata,start=startv,startyear=2001,endyear=2010,steps=4000,show=10,burn=1000)
     results[[spshort[i]]]= fit
   }
 save(results, file="hierarchical models Nouragues.RData")
 return(results)
} 

## do a loop for all species in BCI


 hierarchical.BCI<-function(file=bci)
 
{
   bci1=read.delim(file)   
   inc<-getseeds(bci1$part,bci1$sp,bci1$trap)
   bcifruits<-bci1[inc,]  
   maspp=mostabundspp.bci()
   ma<-data.frame(maspp)
   bcishort<-merge(bcifruits,ma,by.x="sp",by.y="maspp", all.y=T)
   spnames=as.vector(sort(unique(bcishort$sp)))
   #spnameshort=spnames[c(157,149,81)]
   
   results=list()          
   startv=vector()
   for (i in 1:length(spnames))

  {     
     onesp<-subset(bcishort,sp==spnames[i])
     cat("sp:",as.character(spnames[i]), "\n" )
     censsum<-bcicens(infile=onesp)
     bcianal=bciwithzeros(infile=censsum)    
     spdata<-extract.seedfall_bci(dataset=bcianal)
     #spdata=subset(spdata,spdata$year==1998|spdata$year==1999|spdata$year==2000|spdata$year==2001|spdata$year==2002)
     maxdate=which.max(spdata$quantity)
     startv[1]=spdata$yday[maxdate]
     startv[2]= 10
     startv[3]=3.7
     startv[4]= 1.2
     startv[5]= 25
     #fit=model.seedfall.Gibbs(data=spdata,start=startv,startyear=1998,endyear=2002,steps=2000,show=10,burn=1000)
     fit=model.seedfall.Gibbs(data=spdata,start=startv,startyear=1987,endyear=2011,steps=200,show=10,burn=10)
     #save(fit,file=paste(spnames[i],".RData"))
     results[[spnames[i]]]= fit
   }
 save(results, file="hierarchical models BCI.RData")
 return(results)
}


graphierspp<-function(sp="Arrchi", file=nourage)
{
spsh=spshort()
inc<-which(spsh==sp)
fulldata=read.delim(file)
spnames=sort(unique(fulldata$sp))
dataname=extract.seedfall_onesp(latin = spnames[inc], file = nourage)
fitset=y[[inc]]
graph.fitseedrain(spdata=dataname,fit=fitset)

}

#function for constructing a table with the hyperparameters of peak and peakday for all species
resultstablenrghyp<-function(fit=results,file=nourageanalyses,beginyearfile=beginyearfile)
{

beginyr=read.delim(file=beginyearfile)
datafile=read.delim(file)
spnames=names(fit)
species=character()
mu=CImu2=CImu97=SD=CISD2=CISD97=logmu=CIlogmu2=CIlogmu97=logSD=CIlogSD2=CIlogSD97=bestSD=maxlik=vector()
counter=1

for (i in 1:length(spnames))
{ 
f=fit[[spnames[i]]]
speciesused=as.character(spnames[i])  
splong= as.character(splong(infile=file,shortname=speciesused))
bgy= unique(beginyr$beginyr[beginyr$sp==splong])
fitset=retranslate.seedfalldate(fit=f,beginyear=bgy)

bgy= unique(beginyr$beginyr[beginyr$sp==splong])
fitset=retranslate.seedfalldate(fit=f,beginyear=bgy)

species[counter]=splong
mu[counter]=fitset$besthyper[1]
CImu2[counter]= fitset$CIhyper[1]
CImu97[counter]= fitset$CIhyper[2]
SD[counter]=fitset$besthyper[2]
CISD2[counter]= fitset$CIhyper[3]
CISD97[counter]= fitset$CIhyper[4]
logmu[counter]=fitset$besthyper[3]
CIlogmu2[counter]= fitset$CIhyper[5]
CIlogmu97[counter]= fitset$CIhyper[6]
logSD[counter]=fitset$besthyper[4]
CIlogSD2[counter]= fitset$CIhyper[7]
CIlogSD97[counter]= fitset$CIhyper[8]
bestSD[counter]=fitset$bestSD
maxlik[counter]=max(fitset$llike)
counter=counter+1
}
results=data.frame(species,mu,CImu2,CImu97,SD,CISD2,CISD97,logmu,CIlogmu2,CIlogmu97,logSD,CIlogSD2,CIlogSD97,bestSD,maxlik)
return(results)
}

#table for extracting parameters of each year in Nouragues
resultstablenrgyears<-function(fit=results,file=nourageanalyses,beginyearfile=beginyearfile)
{
 
   nrgdata=read.delim(file)
   beginyr=read.delim(file=beginyearfile)
   fulldata=merge(nrgdata,beginyr,by="sp", all.x=TRUE)
   spnames=names(fit)



counter=1

for (i in 1:length(spnames))
{ 
f=fit[[spnames[i]]]
speciesused=as.character(spnames[i])  
splong= as.character(splong(infile=file,shortname=speciesused))
bgy= unique(beginyr$beginyr[beginyr$sp==splong])
fitset=retranslate.seedfalldate(fit=f,beginyear=bgy)

spdata<-extract.seedfall_onesp2(latin =splong, file = file) 
spdata2<-create.rep.year(data=spdata,beginyear=bgy)
nyear= unique(spdata2$year)    
maxyday=numfiles=dif=numeric()
     for (j in 1:length(nyear))
     { 
        oneyr=subset(spdata2,spdata2$year==nyear[j])
        maxyday[j]<-max(oneyr$yday)
        dif[j]<-max(oneyr$yday)-min(oneyr$yday)
        numfiles[j]=dim(oneyr)[1]  
     }
     
     yearselection=data.frame(maxyday,dif,numfiles)
     numyear=dim(yearselection)[1]
     startv=endv=numeric()
     startv= ifelse(yearselection$maxyday[1]>=320&dif[1]>300,nyear[1],nyear[2])
     endv=ifelse(yearselection$maxyday[numyear]>=320&dif[numyear]>300,nyear[numyear],nyear[numyear-1])     
     inc2=which(spdata2$year>=startv&spdata2$year<=endv)
 
  yearsused=unique(spdata$year[inc2])  [1: length(fitset$bestpeakday)]

 peak=peakday=CI2peak=CI97peak=CI2peakday=CI97peakday=year=vector()      
 species=splong
 
 for (k in 1:length(yearsused))
 {
  year[k]=yearsused[k]
  peak[k]=fitset$bestpeak[k]
  first=seq(1,2*length(yearsused),2)
  last=seq(2,2*length(yearsused),2)
  CIpeak2= fitset$CIpeak[first]
  CIpeak97= fitset$CIpeak[last]
  CIpeakday2= fitset$CIpeakday[first]
  CIpeakday97= fitset$CIpeakday[last]
  CI2peak[k]= CIpeak2[k]
  CI97peak[k]= CIpeak97[k]
  peakday[k]=fitset$bestpeakday[k] 
  CI2peakday[k]= CIpeakday2[k]
  CI97peakday[k]= CIpeakday97[k]
  }
  results=data.frame(species, year, peak,CIpeak2,CIpeak97,peakday,CIpeakday2, CIpeakday97)
  if (i==1) allresults=results else allresults=rbind(allresults,results)
  } 

 return(allresults)
}


#results10b=list(results10[["FITR"]],results10[["FIYO"]],results10[["GUA1"]],results10[["GUA2"]],results10[["GUAD"]],results10[["GUAS"]]) 

#this function serves for extracting the hyperparameters of peak and peakday for each species
resultstablebcihyper<-function(fit=results,beginyearfile=bcibeginyear)
{

beginyr=read.delim(file=beginyearfile)
spnames=names(fit)
species=character()
mu=CImu2=CImu97=SD=CISD2=CISD97=logmu=CIlogmu2=CIlogmu97=logSD=CIlogSD2=CIlogSD97=bestSD=maxlik=vector()
counter=1

for (i in 1:length(spnames))
{ 
f=fit[[spnames[i]]]

bgy= unique(beginyr$beginyr[beginyr$sp==spnames[i]])
fitset=retranslate.seedfalldate(fit=f,beginyear=bgy)

species[counter]=spnames[i]
mu[counter]=fitset$besthyper[1]
CImu2[counter]= fitset$CIhyper[1]
CImu97[counter]= fitset$CIhyper[2]
SD[counter]=fitset$besthyper[2]
CISD2[counter]= fitset$CIhyper[3]
CISD97[counter]= fitset$CIhyper[4]
logmu[counter]=fitset$besthyper[3]
CIlogmu2[counter]= fitset$CIhyper[5]
CIlogmu97[counter]= fitset$CIhyper[6]
logSD[counter]=fitset$besthyper[4]
CIlogSD2[counter]= fitset$CIhyper[7]
CIlogSD97[counter]= fitset$CIhyper[8]
bestSD[counter]=fitset$bestSD
maxlik[counter]=max(fitset$llike)
counter=counter+1
}
results=data.frame(species,mu,CImu2,CImu97,SD,CISD2,CISD97,logmu,CIlogmu2,CIlogmu97,logSD,CIlogSD2,CIlogSD97,bestSD, maxlik)
return(results)
}


 
#this function extracts parameters of peak and peakday for each species and year
resultstablebciyears<-function(fit=results,beginyearfile=bcibeginyear,shortfile=bcishort)
{
 beginyr=read.delim(file=beginyearfile)
 bcishort=read.delim(file=shortfile)
 spnames=names(fit)

 for (i in 1:length(spnames))
{ 
  f=fit[[spnames[i]]]

  bgy= unique(beginyr$beginyr[beginyr$sp==spnames[i]])
  fitset=retranslate.seedfalldate(fit=f,beginyear=bgy)
#this part is for the extraction of the years that were used for the model running
 onesp<-subset(bcishort,sp==spnames[i]) 
 censsum<-bcicens(infile=onesp) #this function sums the quantity values of the traps of the same censuses 
 bcianal=bciwithzeros(infile=censsum) #bciwithzeros adds all the censuses to the dataset for each single species   
 spdata<-extract.seedfall_bci(dataset=bcianal)
 bgy= unique(bcishort$beginyr[bcishort$sp==spnames[i]])
 spdata2<-create.rep.year(data=spdata,beginyear=bgy)
 nyear= unique(spdata2$year)    
 maxyday=numfiles=dif=numeric()
     for (j in 1:length(nyear))
     { 
        oneyr=subset(spdata2,spdata2$year==nyear[j])
        maxyday[j]<-max(oneyr$yday)
        dif[j]<-max(oneyr$yday)-min(oneyr$yday)
        numfiles[j]=dim(oneyr)[1]  
     }
     
  yearselection=data.frame(maxyday,dif,numfiles)
  numyear=dim(yearselection)[1]
  startv=endv=numeric()
  startv= ifelse(yearselection$maxyday[1]>=320&dif[1]>300,nyear[1],nyear[2])
  endv=ifelse(yearselection$maxyday[numyear]>=320&dif[numyear]>300,nyear[numyear],nyear[numyear-1])
    
  inc2=which(spdata2$year>=startv&spdata2$year<=endv)
  yearsused=unique(spdata$year[inc2])  [1: length(fitset$bestpeakday)]

 peak=peakday=CI2peak=CI97peak=CI2peakday=CI97peakday=year=vector()         
 species=spnames[i]
 
 for (k in 1:length(yearsused))
 {
  year[k]=yearsused[k]
  peak[k]=fitset$bestpeak[k]
  first=seq(1,2*length(yearsused),2)
  last=seq(2,2*length(yearsused),2)
  CIpeak2= fitset$CIpeak[first]
  CIpeak97= fitset$CIpeak[last]
  CIpeakday2= fitset$CIpeakday[first]
  CIpeakday97= fitset$CIpeakday[last]
  CI2peak[k]= CIpeak2[k]
  CI97peak[k]= CIpeak97[k]
  peakday[k]=fitset$bestpeakday[k] 
  CI2peakday[k]= CIpeakday2[k]
  CI97peakday[k]= CIpeakday97[k]
  }
  results=data.frame(species, year, peak,CIpeak2,CIpeak97,peakday,CIpeakday2, CIpeakday97)
  if (i==1) allresults=results else allresults=rbind(allresults,results)
  } 

 return(allresults)
}

#this function serves for extracting the hyperparameters of peak and peakday for each species
resultsbcihypertrans<-function(fit=results,beginyearfile=bcibeginyear)
{

beginyr=read.delim(file=beginyearfile)
spnames=names(fit)
species=character()
mu=CImu2=CImu97=SD=CISD2=CISD97=logmu=CIlogmu2=CIlogmu97=logSD=CIlogSD2=CIlogSD97=bestSD=maxlik=vector()
counter=1

for (i in 1:length(spnames))
{ 
fitset=fit[[spnames[i]]]

bgy= unique(beginyr$beginyr[beginyr$sp==spnames[i]])


species[counter]=spnames[i]
mu[counter]=fitset$besthyper[1]
CImu2[counter]= fitset$CIhyper[1]
CImu97[counter]= fitset$CIhyper[2]
SD[counter]=fitset$besthyper[2]
CISD2[counter]= fitset$CIhyper[3]
CISD97[counter]= fitset$CIhyper[4]
logmu[counter]=fitset$besthyper[3]
CIlogmu2[counter]= fitset$CIhyper[5]
CIlogmu97[counter]= fitset$CIhyper[6]
logSD[counter]=fitset$besthyper[4]
CIlogSD2[counter]= fitset$CIhyper[7]
CIlogSD97[counter]= fitset$CIhyper[8]
bestSD[counter]=fitset$bestSD
maxlik[counter]=max(fitset$llike)
counter=counter+1
}
results=data.frame(species,mu,CImu2,CImu97,SD,CISD2,CISD97,logmu,CIlogmu2,CIlogmu97,logSD,CIlogSD2,CIlogSD97,bestSD, maxlik)
return(results)
}



#this function extracts parameters of peak and peakday for each species and year
resultsbciyearstrans<-function(fit=results,beginyearfile=bcibeginyear,shortfile=bcishort)
{
 beginyr=read.delim(file=beginyearfile)
 bcishort=read.delim(file=shortfile)
 spnames=names(fit)

 for (i in 1:length(spnames))
{ 
  fitset=fit[[spnames[i]]]

  bgy= unique(beginyr$beginyr[beginyr$sp==spnames[i]])
#this part is for the extraction of the years that were used for the model running
 onesp<-subset(bcishort,sp==spnames[i]) 
 censsum<-bcicens(infile=onesp) #this function sums the quantity values of the traps of the same censuses 
 bcianal=bciwithzeros(infile=censsum) #bciwithzeros adds all the censuses to the dataset for each single species   
 spdata<-extract.seedfall_bci(dataset=bcianal)
 bgy= unique(bcishort$beginyr[bcishort$sp==spnames[i]])
 spdata2<-create.rep.year(data=spdata,beginyear=bgy)
 nyear= unique(spdata2$year)    
 maxyday=numfiles=dif=numeric()
     for (j in 1:length(nyear))
     { 
        oneyr=subset(spdata2,spdata2$year==nyear[j])
        maxyday[j]<-max(oneyr$yday)
        dif[j]<-max(oneyr$yday)-min(oneyr$yday)
        numfiles[j]=dim(oneyr)[1]  
     }
     
  yearselection=data.frame(maxyday,dif,numfiles)
  numyear=dim(yearselection)[1]
  startv=endv=numeric()
  startv= ifelse(yearselection$maxyday[1]>=320&dif[1]>300,nyear[1],nyear[2])
  endv=ifelse(yearselection$maxyday[numyear]>=320&dif[numyear]>300,nyear[numyear],nyear[numyear-1])
    
  inc2=which(spdata2$year>=startv&spdata2$year<=endv)
  yearsused=unique(spdata$year[inc2])  [1: length(fitset$bestpeakday)]

 peak=peakday=CI2peak=CI97peak=CI2peakday=CI97peakday=year=vector()         
 species=spnames[i]
 
 for (k in 1:length(yearsused))
 {
  year[k]=yearsused[k]
  peak[k]=fitset$bestpeak[k]
  first=seq(1,2*length(yearsused),2)
  last=seq(2,2*length(yearsused),2)
  CIpeak2= fitset$CIpeak[first]
  CIpeak97= fitset$CIpeak[last]
  CIpeakday2= fitset$CIpeakday[first]
  CIpeakday97= fitset$CIpeakday[last]
  CI2peak[k]= CIpeak2[k]
  CI97peak[k]= CIpeak97[k]
  peakday[k]=fitset$bestpeakday[k] 
  CI2peakday[k]= CIpeakday2[k]
  CI97peakday[k]= CIpeakday97[k]
  }
  results=data.frame(species, year, peak,CIpeak2,CIpeak97,peakday,CIpeakday2, CIpeakday97)
  if (i==1) allresults=results else allresults=rbind(allresults,results)
  } 

 return(allresults)
}

#this function pastes all the results of the hierarchical models for BCI in a single table
allresultsbci=function(file=bci,beginyearfile=bcibeginyear)
{
  exc=c(130,54,71,21,27,44,85,77)    #this is fhe number of species included in results17b (it is needed to replace these species in the other files by the new ones)
  
  thyper7=resultstablebcihyper(fit=results7, beginyearfile=bcibeginyear)
  thyper13=resultstablebcihyper(fit=results13, beginyearfile=bcibeginyear)
  thyper15=resultstablebcihyper(fit=results15, beginyearfile=bcibeginyear)
  thyper9=resultstablebcihyper(fit=results9, beginyearfile=bcibeginyear)
  thyper14=resultstablebcihyper(fit=results14, beginyearfile=bcibeginyear)
  thyper11=resultstablebcihyper(fit=results11, beginyearfile=bcibeginyear)
  thyper12=resultstablebcihyper(fit=results12, beginyearfile=bcibeginyear)
  thyper16=resultstablebcihyper(fit=results16, beginyearfile=bcibeginyear)
  thyper17b=resultstablebcihyper(fit=results17b, beginyearfile=bcibeginyear)
  tothyp=rbind(thyper7,thyper13,thyper15,thyper9,thyper14,thyper11,thyper12,thyper16)
  tothyp2= tothyp[-exc,]
  tothypfinal=rbind(tothyp2,thyper17b)
  o=order(tothypfinal$species)
  write.table(tothypfinal[o,], file="results hyperparameters BCI.txt", sep="\t", row.names=FALSE)
  
  
  tyears7=resultstablebciyears(fit=results7,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears13=resultstablebciyears(fit=results13,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears15=resultstablebciyears(fit=results15,beginyearfile=bcibeginyear,shortfile=bcishort) 
  tyears9=resultstablebciyears(fit=results9,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears14=resultstablebciyears(fit=results14,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears11=resultstablebciyears(fit=results11,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears12=resultstablebciyears(fit=results12,beginyearfile=bcibeginyear,shortfile=bcishort) 
  tyears16=resultstablebciyears(fit=results16,beginyearfile=bcibeginyear,shortfile=bcishort)  
  tyears17b=resultstablebciyears(fit=results17b,beginyearfile=bcibeginyear,shortfile=bcishort)
  totyears=rbind(tyears7,tyears13,tyears15,tyears9,tyears14,tyears11,tyears12,tyears16)
  sora=which(totyears$species=="SORA")
  fic2=which(totyears$species=="FIC2")
  hipv=which(totyears$species=="HIPV")
  astg=which(totyears$species=="ASTG")
  ceip=which(totyears$species=="CEIP")
  desp=which(totyears$species=="DESP")
  masn=which(totyears$species=="MASN")
  jacc=which(totyears$species=="JACC")
  exc2=c(sora, fic2,hipv,astg,ceip,desp,masn,jacc)
  totyears2=totyears[-exc2,]
  totyearfinal=rbind(totyears2,tyears17b)
  o=order(totyearfinal$species)
 
 write.table(totyearfinal[o,], file="results parameters per year BCI.txt", sep="\t", row.names=FALSE)
}

#this function pastes all the results of the hierarchical models for BCI in a single table
allresultsbcitrans=function(file=bci,beginyearfile=bcibeginyear)
{
  exc=c(130,54,71,21,27,44,85,77)    #this is fhe number of species included in results17b (it is needed to replace these species in the other files by the new ones)
  
  thyper7=resultsbcihypertrans(fit=results7, beginyearfile=bcibeginyear)
  thyper13=resultsbcihypertrans(fit=results13, beginyearfile=bcibeginyear)
  thyper15=resultsbcihypertrans(fit=results15, beginyearfile=bcibeginyear)
  thyper9=resultsbcihypertrans(fit=results9, beginyearfile=bcibeginyear)
  thyper14=resultsbcihypertrans(fit=results14, beginyearfile=bcibeginyear)
  thyper11=resultsbcihypertrans(fit=results11, beginyearfile=bcibeginyear)
  thyper12=resultsbcihypertrans(fit=results12, beginyearfile=bcibeginyear)
  thyper16=resultsbcihypertrans(fit=results16, beginyearfile=bcibeginyear)
  thyper17b=resultsbcihypertrans(fit=results17b, beginyearfile=bcibeginyear)
  tothyp=rbind(thyper7,thyper13,thyper15,thyper9,thyper14,thyper11,thyper12,thyper16)
  tothyp2= tothyp[-exc,]
  tothypfinal=rbind(tothyp2,thyper17b)
  o=order(tothypfinal$species)
  write.table(tothypfinal[o,], file="results hyperparameters BCI translated.txt", sep="\t", row.names=FALSE)
  
  
  tyears7=resultsbciyearstrans(fit=results7,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears13=resultsbciyearstrans(fit=results13,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears15=resultsbciyearstrans(fit=results15,beginyearfile=bcibeginyear,shortfile=bcishort) 
  tyears9=resultsbciyearstrans(fit=results9,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears14=resultsbciyearstrans(fit=results14,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears11=resultsbciyearstrans(fit=results11,beginyearfile=bcibeginyear,shortfile=bcishort)
  tyears12=resultsbciyearstrans(fit=results12,beginyearfile=bcibeginyear,shortfile=bcishort) 
  tyears16=resultsbciyearstrans(fit=results16,beginyearfile=bcibeginyear,shortfile=bcishort)  
  tyears17b=resultsbciyearstrans(fit=results17b,beginyearfile=bcibeginyear,shortfile=bcishort)
  totyears=rbind(tyears7,tyears13,tyears15,tyears9,tyears14,tyears11,tyears12,tyears16)
  sora=which(totyears$species=="SORA")
  fic2=which(totyears$species=="FIC2")
  hipv=which(totyears$species=="HIPV")
  astg=which(totyears$species=="ASTG")
  ceip=which(totyears$species=="CEIP")
  desp=which(totyears$species=="DESP")
  masn=which(totyears$species=="MASN")
  jacc=which(totyears$species=="JACC")
  exc2=c(sora, fic2,hipv,astg,ceip,desp,masn,jacc)
  totyears2=totyears[-exc2,]
  totyearfinal=rbind(totyears2,tyears17b)
  o=order(totyearfinal$species)
 
 write.table(totyearfinal[o,], file="results parameters per year BCI translated.txt", sep="\t", row.names=FALSE)
}



graph.allpeakdays=function()
{
t=resultstable(fit=results,nrgfile=nourageanalyses, beginyearfile=beginyearfile)

}

modeltable<-function(fit=results,nrgfile=nourageanalyses,beginyearfile=beginyearfile)
{

beginyr=read.delim(file=beginyearfile)
datafile=read.delim(file=nrgfile)
spnames=unique(datafile$sp)
species=character()
model=vector()
counter=1

for (i in 1:length(spnames))
{ 
species=as.character(spnames[i])
sphort=spshort2(spname=species)
f=fit[[spshort]]

bgy= unique(beginyr$beginyr[beginyr$sp==species])
fitset=retranslate.seedfalldate(fit=f,beginyear=bgy)

species[counter]=splong
mu[counter]=fitset$besthyper[1]
CImu2[counter]= fitset$CIhyper[1]
CImu97[counter]= fitset$CIhyper[2]
SD[counter]=fitset$besthyper[2]
CISD2[counter]= fitset$CIhyper[3]
CISD97[counter]= fitset$CIhyper[4]
logmu[counter]=fitset$besthyper[3]
CIlogmu2[counter]= fitset$CIhyper[5]
CIlogmu97[counter]= fitset$CIhyper[6]
logSD[counter]=fitset$besthyper[4]
CIlogSD2[counter]= fitset$CIhyper[7]
CIlogSD97[counter]= fitset$CIhyper[8]
counter=counter+1
 

}                               
results=data.frame(species,mu,CImu2,CImu97,SD,CISD2,CISD97,logmu,CIlogmu2,CIlogmu97,logSD,CIlogSD2,CIlogSD97)
return(results)
}

graph.allpeakdays=function()
{
t=resultstable(fit=results,nrgfile=nourageanalyses, beginyearfile=beginyearfile)


}


#some graphs integrating hyperparameters of species

histhyper=function()
{
  jpeg(filename="histhyperparameters.jpg",width = 600, height = 600,quality=100)
  par(mfrow=c(2,2))
  t=resultstable(fit=results,nrgfile=nourageanalyses, beginyearfile=beginyearfile)
  t=t[-18,]
  hist(t$mu,xlab="day of the year",ylab="number of species",main="")
  hist(t$SD,xlab="SD of day of the year",ylab="number of species",main="")
  hist(exp(t$logmu)+1,xlab="peak of seed production (number seeds/year)",ylab="number of species",main="")
  hist(exp(t$logSD)+1,xlab="SD of peak of seed production (number seeds/year)",ylab="number of species",main="")
   dev.off()
}

histhyp=function(file="BCI results hyperparameters.txt", site="BCI")
{
 t=read.delim(file)
 # nrghyp=read.delim(file="Nouragues results hyperparameters.txt")
  pdf(file="hist paryear.pdf")
  par(mfrow=c(2,2))
  hist(t$mu,xlab="day of the year",ylab="number of species",main=site)
  hist(t$SD,xlab="SD of day of the year",ylab="number of species",main=site)
  hist(exp(t$logmu)+1,xlab="peak of seed production (number seeds/year)",ylab="number of species",main=site)
  hist(exp(t$logSD)+1,xlab="SD of peak of seed production (number seeds/year)",ylab="number of species",main=site)
   dev.off()
}
#curvesallspp(file="Nouragues results hyperparameters.txt",site="Nouragues")
curvesallspp=function()

{
setwd("C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
tr=read.delim(file="BCI results hyperparameters.txt")
#tr=read.delim(file="table some results BCI.txt")
negval= which(tr$logmu<=0)
negdate=which(tr$mu<=0)
if(length(negval)>=1) {tr=tr[-negval,]} 
if(length(negdate)>0) {tr=tr[-negval,]}

x11(height=9,width=9)
par(mfrow=c(2,1))

 maxyr=which.max(tr$logmu)

   oneyrpred=14*tr$logmu[maxyr]*dnorm(1:400,mean=tr$mu[maxyr],sd=tr$bestSD[maxyr])
   plot(1:400,oneyrpred,type='l',ylab='Weekly seedfall',xlab='day of year',ylim=c(0,max(oneyrpred+1)),lwd=1.25)
   polygon(c(0,120,120,0),c(0,0,3,3),col = "#FF000020", border = NA)
   polygon(c(120,400,400,120),c(0,0,3,3),col = "#0000FF20", border = NA)
   title(main="BCI (1987-2011)")
    for (i in 1:dim(tr)[1])
{
oneyrpred=14*tr$logmu[i]*dnorm(1:400,mean=tr$mu[i],sd=tr$bestSD[i])
lines(1:400,oneyrpred,col="blue",lwd=1.25) 
}



tr=read.delim(file="Nouragues results hyperparameters.txt")
negval= which(tr$logmu<=0)
negdate=which(tr$mu<=0)
if(length(negval)>=1) {tr=tr[-negval,]} 
if(length(negdate)>0) {tr=tr[-negval,]}

maxyr=which.max(tr$logmu)


 
   oneyrpred=14*tr$logmu[maxyr]*dnorm(1:400,mean=tr$mu[maxyr],sd=tr$bestSD[maxyr])
   plot(1:400,oneyrpred,type='l',ylab='Biweekly seedfall',xlab='day of year',ylim=c(0,max(oneyrpred+1)),lwd=1.25)
   polygon(c(248,341,341,248),c(0,0,3,3),col = "#FF000020", border = NA)
   polygon(c(0,248,248,0),c(0,0,3,3),col = "#0000FF20", border = NA)
   polygon(c(341,400,400,341),c(0,0,3,3),col = "#0000FF20", border = NA)
   title(main="Nouragues (2001-2011)")
    for (i in 1:dim(tr)[1])
{
oneyrpred=14*tr$logmu[i]*dnorm(1:400,mean=tr$mu[i],sd=tr$bestSD[i])
lines(1:400,oneyrpred,col="blue", lwd=1.25) 
}

#lines(1:400, 14*mean(tr$logmu)*dnorm(1:400,mean=mean(tr$mu),sd=mean(tr$bestSD)),col="black",lheight=2)
}

nrg.graphs.pdf<-function(file=nourage, fit=results, beginyearfile=beginyearfile)
{
   nrgdata=read.delim(file)
   beginyr=read.delim(file=beginyearfile)
   fulldata=merge(nrgdata,beginyr,by="sp", all.x=TRUE)
   spnames=sort(unique(fulldata$sp))
   
    
   for (i in 1:length(spnames))
  {     
     species=as.character(spnames[i])
     bgy=beginyr$beginyr[beginyr$sp==species]
     spdata<-extract.seedfall_onesp2(latin =species, file = file) 
     spdata2<-create.rep.year(data=spdata,beginyear=bgy)
     nyear= unique(spdata2$year)    
     maxyday=numfiles=dif=numeric()
     for (j in 1:length(nyear))
     { 
        oneyr=subset(spdata2,spdata2$year==nyear[j])
        maxyday[j]<-max(oneyr$yday)
        dif[j]<-max(oneyr$yday)-min(oneyr$yday)
        numfiles[j]=dim(oneyr)[1]  
     }
     
     yearselection=data.frame(maxyday,dif,numfiles)
     numyear=dim(yearselection)[1]
     startv=endv=numeric()
     startv= ifelse(yearselection$maxyday[1]>=320&dif[1]>300,nyear[1],nyear[2])
     endv=ifelse(yearselection$maxyday[numyear]>=320&dif[numyear]>300,nyear[numyear],nyear[numyear-1])     
     inc2=which(spdata2$year>=startv&spdata2$year<=endv)
     
     spsh=spshort2(spname=species)
     fitset=fit[[spsh]] 
     #fittrans= retranslate.seedfalldate2(data=spdata,fit=fitset,beginyear=bgy,startyear=startyear,endyear=endyear)
     fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
      
    filename=paste(species,".pdf")   
    pdf(file=filename)
     
    par(mfrow=c(2,1))

    #plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,max(fittrans$bestpeak)))
    plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="years",ylab="number of seeds")
    jan=tojulian(pst('01/01/',min(unique(spdata$year)):(max(unique(spdata$year))+1)))
    #axis(side=1, at=jan, labels= c(min(unique(spdata$year)):max(unique(spdata$year))))
    axis(side=1, at=jan, labels= jan)
    axis(2)
    Dj=round(ifelse(bgy<182.5,bgy,bgy-366),0)
    lines(spdata2$julian[inc2]+Dj,fitset$model)
    title(main=species)  
    newdates=fromjulian(spdata$julian[inc2]+Dj,dateform="%Y-%m-%d")
    fulldate=create.fulldate(newdates,format='%Y-%m-%d')
    startyear=min(unique(fulldate$year))
    endyear=max(unique(fulldate$year))
    jan1=tojulian(pst('01/01/',startv:endv))
    fitpeakjulian=jan1+fittrans$bestpeakday 
    points(fitpeakjulian,fittrans$bestpeak,col='red',pch=16)
    
   
   maxyr=which.max(fittrans$bestpeak)
   oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(1:400,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
   plot(1:400,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year')
   for(k in 1:length(fittrans$bestpeak)) 
    { 
     oneyrpred=14*fittrans$bestpeak[k]*dnorm(1:400,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
     lines(1:400,oneyrpred,col=k) 
    }
 dev.off()
 }
}

 
bci.graph.pdf<-function(file=bcishort, fit=bciresults2, beginyearfile=bcibeginyear,spstart=101,spend=166)
{
  
   bcishort<-read.delim(file)
   #sprange=spstart:spend
   #spnames=as.vector(sort(unique(bcishort$sp)))
   spnames=names(fit)
   #spnameshort=spnames[c(157,149,81)]
   
   results=list()          
   startv=vector()
   #for (i in 1:length(sprange))
   for (i in 1:length(spnames))
                       
  {     
     #species=as.character(spnames[sprange][i])
     species=as.character(spnames[i])
     onesp<-subset(bcishort,sp==species)
     cat("sp:",as.character(species), "\n" )
     censsum<-bcicens(infile=onesp)
     bcianal=bciwithzeros(infile=censsum)    
     spdata<-extract.seedfall_bci(dataset=bcianal)
     bgy= unique(bcishort$beginyr[bcishort$sp==species])
     spdata2<-create.rep.year(data=spdata,beginyear=bgy)
     #spdata<-create.rep.year(data=spdata,beginyear=bgy)
    
     Dj=round(ifelse(bgy<182.5,bgy,bgy-366),0)
     nyear= unique(spdata2$year)    
     maxyday=numfiles=dif=numeric()
     for (j in 1:length(nyear))
     { 
        oneyr=subset(spdata2,spdata2$year==nyear[j])
        maxyday[j]<-max(oneyr$yday)
        dif[j]<-max(oneyr$yday)-min(oneyr$yday)
        numfiles[j]=dim(oneyr)[1]  
     }
     
     yearselection=data.frame(maxyday,dif,numfiles)
     numyear=dim(yearselection)[1]
     startv=endv=numeric()
     startv= ifelse(yearselection$maxyday[1]>=320&dif[1]>300,nyear[1],nyear[2])
     endv=ifelse(yearselection$maxyday[numyear]>=320&dif[numyear]>300,nyear[numyear],nyear[numyear-1]) 
     inc2=which(spdata2$year>=startv&spdata2$year<=endv)
     fitset=fit[[species]] 
     #fittrans= retranslate.seedfalldate2(data=spdata,fit=fitset,beginyear=bgy,startyear=startyear,endyear=endyear)
     fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
     
     filename=paste(species,".pdf")   
     pdf(file=filename)

     
    par(mfrow=c(2,1))
    
    plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="years",ylab="number of seeds") 
    jan=tojulian(pst('01/01/',min(unique(spdata$year)):(max(unique(spdata$year))+1)))
    janlab=pst('1Jan',min(unique(spdata$year)):(max(unique(spdata$year))+1))
    axis(side=1, at=jan, labels= janlab)
    axis(2)  
    ejex=as.vector(spdata2$julian[inc2]+Dj)
    lines(ejex[1:length(fitset$model)],fitset$model)  
    title(main=species)
    newdates=fromjulian(spdata$julian[inc2]+Dj,dateform="%Y-%m-%d")
    fulldate=create.fulldate(newdates,format='%Y-%m-%d')
    startyear=min(unique(fulldate$year))
    endyear=max(unique(fulldate$year))
    jan1=tojulian(pst('01/01/',startv:endv))
    fitpeakjulian=jan1+fittrans$bestpeakday 
    points(fitpeakjulian,fittrans$bestpeak,col='red',pch=16)
   
   maxyr=which.max(fittrans$bestpeak)
   oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(1:400,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
   plot(1:400,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year')
   for(k in 1:length(fittrans$bestpeak)) 
    { 
     oneyrpred=14*fittrans$bestpeak[k]*dnorm(1:400,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
     lines(1:400,oneyrpred,col=k) 
    }
 dev.off()
}
}  
#this function repeat models just for a subset of species that were not working with other beginyear priors 
somerepetBCI<-function(file=bci,beginyearfile=bcibeginyear, steps=10000, burnin=1000,outfile)
 
{
   bci1=read.delim(file)  
   beginyr=read.delim(file=beginyearfile) 
   inc<-getseedsbci200traps(bci1$part,bci1$sp,bci1$trap)
   bcifruits<-bci1[inc,]  
  #mostabundsppbci contains the most abundant species according to the criterion in Wright & Calderon: seeds or fruits captured in 10 or more traps in any single year 
   ma<-read.delim(file=mostabundsppbci)
   ma2<-as.data.frame(ma)
   bcimaspp<-merge(bcifruits,ma2,by.x="sp",by.y="spma", all.y=T)
   bcishort<-merge(bcimaspp,beginyr,by.x="sp",by.y="sp", all.x=T)
   inc2=which(is.na(bcishort$beginyr)==T)
   bcishort=bcishort[-c(inc2),]
   spnames=c("SORA","FIC2","HIPV","ASTG","CEIP")
   #spnameshort=spnames[c(157,149,81)]
   
   results17=list()          
   startv=vector()
   for (i in 1:length(spnames))
                       
  {     
     species=as.character(spnames[i])
     onesp<-subset(bcishort,sp==species)
     cat("sp:",as.character(species), "\n")
     #cat("number:",i, "\n" )
     censsum<-bcicens(infile=onesp) #this function sums the quantity values of the traps of the same censuses 
     bcianal=bciwithzeros(infile=censsum) #bciwithzeros adds all the censuses to the dataset for each single species   
     spdata<-extract.seedfall_bci(dataset=bcianal)
     bgy= unique(bcishort$beginyr[bcishort$sp==species])
     spdata2<-create.rep.year(data=spdata,beginyear=bgy)
     #spdata=subset(spdata,spdata$year==1998|spdata$year==1999|spdata$year==2000|spdata$year==2001|spdata$year==2002)
     nyear= unique(spdata2$year)    
     maxyday=numfiles=dif=numeric()
     for (j in 1:length(nyear))
     { 
        oneyr=subset(spdata2,spdata2$year==nyear[j])
        maxyday[j]<-max(oneyr$yday)
        dif[j]<-max(oneyr$yday)-min(oneyr$yday)
        numfiles[j]=dim(oneyr)[1]  
     }
     
     yearselection=data.frame(maxyday,dif,numfiles)
     numyear=dim(yearselection)[1]
     startv=endv=numeric()
     startv= ifelse(yearselection$maxyday[1]>=320&dif[1]>300,nyear[1],nyear[2])
     endv=ifelse(yearselection$maxyday[numyear]>=320&dif[numyear]>300,nyear[numyear],nyear[numyear-1])
     peakdays= rep(182,length(startv:endv))
     peaks= log(initial.peak(data =spdata2,startyear=startv,endyear=endv)+1)+1
     fit=model.seedfall.Gibbs2(data=spdata2,peakdays=peakdays,peaks=peaks,SDvector=c(10,1,35),startyear=startv,endyear=endv,steps=steps,showstep=200,burnin=burnin)
     #fit=model.seedfall.Gibbs(data=spdata,start=startv,startyear=1998,endyear=2002,steps=2000,show=10,burn=1000)
     results17[[spshort2(spname=species)]]= fit

     #save(fit,file=paste(spnames[i],".RData"))
     results17[[species]]= fit
   }
 save(results17, file=outfile)
 return(results17)
}

 
bci.graphsmoved.arr1=function(file=bcishort, fit=results14,species="ARR1")         #this graph put 0 days of the year in the center: it is only valid for LOW numbers of peakday
{
   
   bcishort<-read.delim (file) 
   spnames=as.vector(sort(unique(bcishort$sp)))
   
     results=list()          
     startv=vector()

     onesp<-subset(bcishort,sp==species)
     cat("sp:",as.character(species), "\n" )
     censsum<-bcicens(infile=onesp)
     bcianal=bciwithzeros(infile=censsum)    
     spdata<-extract.seedfall_bci(dataset=bcianal)
     bgy= unique(bcishort$beginyr[bcishort$sp==species])
     spdata2<-create.rep.year(data=spdata,beginyear=bgy)
     #spdata<-create.rep.year(data=spdata,beginyear=bgy)
    
     Dj=round(ifelse(bgy<182.5,bgy,bgy-366),0)
     nyear= unique(spdata2$year)    
     maxyday=numfiles=dif=numeric()
     for (j in 1:length(nyear))
     { 
        oneyr=subset(spdata2,spdata2$year==nyear[j])
        maxyday[j]<-max(oneyr$yday)
        dif[j]<-max(oneyr$yday)-min(oneyr$yday)
        numfiles[j]=dim(oneyr)[1]  
     }
     
     yearselection=data.frame(maxyday,dif,numfiles)
     numyear=dim(yearselection)[1]
     startv=endv=numeric()
     startv= ifelse(yearselection$maxyday[1]>=320&dif[1]>300,nyear[1],nyear[2])
     endv=ifelse(yearselection$maxyday[numyear]>=320&dif[numyear]>300,nyear[numyear],nyear[numyear-1]) 
     inc2=which(spdata2$year>=startv&spdata2$year<=endv)
     fitset=fit[[species]] 
     #fittrans= retranslate.seedfalldate2(data=spdata,fit=fitset,beginyear=bgy,startyear=startyear,endyear=endyear)
     fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
     
     filename=paste(species,"moved2.pdf")   
     pdf(file=filename)

     
    par(mfrow=c(2,1))
    inc=which(spdata$year>=startv&spdata$year<=endv)
    plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="years",ylab="number of seeds")
    jan=tojulian(pst('01/01/',min(unique(spdata$year)):(max(unique(spdata$year))+1)))
    janlab=pst('1Jan',min(unique(spdata$year)):(max(unique(spdata$year))+1))
    axis(side=1, at=jan, labels= janlab)
    axis(2)
    ejex=as.vector(spdata2$julian[inc2]+Dj)
    lines(ejex[1:length(fitset$model)],fitset$model)  
    title(main=species)
    newdates=fromjulian(spdata$julian[inc2]+Dj,dateform="%Y-%m-%d")
    fulldate=create.fulldate(newdates,format='%Y-%m-%d')
    startyear=min(unique(fulldate$year))
    endyear=max(unique(fulldate$year))
    jan1=tojulian(pst('01/01/',startv:endv))
    fitpeakjulian=jan1+fittrans$bestpeakday 
    points(fitpeakjulian,fittrans$bestpeak,col='red',pch=16)
   
   maxyr=which.max(fittrans$bestpeak)
   oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(-120:245,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
   plot(-120:245,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year',axes=FALSE)
   axis(side=1, at=c(-120,-60,1,60,120,180,240),labels=c("245","305","0","60","120","180","240") )
   axis(side=2)
   for(k in 1:length(fittrans$bestpeak)) 
    { 
     oneyrpred=14*fittrans$bestpeak[k]*dnorm(-120:245,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
     lines(-120:245,oneyrpred,col=k) 
    }
 dev.off()
}


boxplot.BCI=function(parameter="peakday")                           #this function draws the plots for the estimated parameters of peak or peakday for BCI. The desired parameter must be indicated
{
  setwd("C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  bcy=read.delim(file="BCI results parameters per year.txt")
  spnames=sort(unique(bcy$species))

    for (i in 1:length(spnames))
                       
  {     
  sp=as.character(spnames[i])
   cat("sp:",as.character(sp), "\n" ) 
  onesp=subset(bcy,bcy$species==sp)
  filename=paste(sp,"year.pdf")   
  pdf(file=filename)
  color=c("red","blue","grey","grey","red","red","red","grey","grey","grey","red","red",
  "grey","grey","grey","red","grey","grey","grey","grey","blue","blue","purple","blue","blue")
  barplot(onesp$peak, ylim=c(0,max(onesp$CIpeak97)), space=1,axes=FALSE, col=color[1:dim(onesp)[1]],xlab="years",ylab="estimated number of seeds")
  title(main=sp)
  xpos=seq(1.5,(dim(onesp)[1])*2,2)
  axis(side =1, at=xpos,lab=as.character(c(min(bcy$year):max(bcy$year))[1:length(xpos)]))
  axis(side=2)
  errbar(x=xpos,y=onesp$peak,yplus=onesp$CIpeak2,yminus=onesp$CIpeak97,add=TRUE)
  dev.off()
  }
}

boxplot.NRG=function(parameter="peak")           #this function draws the plots for the estimated parameters of peak or peakday for Nouragues. The desired parameter must be indicated
{
  setwd("C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  bcy=read.delim(file="Nouragues results parameters per year.txt")
  spnames=sort(unique(bcy$species))

    for (i in 1:length(spnames))
                       
  {     
  sp=as.character(spnames[i])
   cat("sp:",as.character(sp), "\n" ) 
  onesp=subset(bcy,bcy$species==sp)
  
  color=c("grey","red","grey","grey","grey","grey","blue","blue","purple","blue","blue")
  if(parameter=="peak") 
  {
    filename=paste(sp,"peak year.pdf")   
  pdf(file=filename)
  barplot(onesp$peak, ylim=c(0,max(onesp$CIpeak97)), space=1,axes=FALSE, col=color[1:dim(onesp)[1]],xlab="years",ylab="estimated number of seeds")
  title(main=sp)
  xpos=seq(1.5,(dim(onesp)[1])*2,2)
  axis(side =1, at=xpos,lab=as.character(c(min(bcy$year):max(bcy$year))[1:length(xpos)]))
  axis(side=2)
  errbar(x=xpos,y=onesp$peak,yplus=onesp$CIpeak2,yminus=onesp$CIpeak97,add=TRUE)
  }
   
  if(parameter=="peakday") 
  {
  filename=paste(sp,"peakday year.pdf")   
  pdf(file=filename)
  barplot(onesp$peakday, ylim=c(min(onesp$peakday),max(onesp$CIpeakday97)), space=1,axes=FALSE, col=color[1:dim(onesp)[1]],xlab="years",ylab="estimated day of the year")
  title(main=sp)
  xpos=seq(1.5,(dim(onesp)[1])*2,2)
  axis(side =1, at=xpos,lab=as.character(c(min(bcy$year):max(bcy$year))[1:length(xpos)]))
  axis(side=2)
  errbar(x=xpos,y=onesp$peakday,yplus=onesp$CIpeakday2,yminus=onesp$CIpeakday97,add=TRUE)
   
  }
  dev.off()
  }
}