
setwd("C:/Irene/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models")
load("C:\\Irene\\Brunoy\\Base Datos Nouragues\\Metadata Joe Wright\\modeling max likelihood\\CTFSRPackage.rdata")
load("C:\\Irene\\Brunoy\\hierarchical models\\working files hierarchical models\\NRG-results.Rdata")
#source("C:\\Brunoy\\Base Datos Nouragues\\Metadata Joe Wright\\modeling max likelihood\\BCI-max likelihood.r")
#attach('CTFSRPackage.rdata')

library(date)

#### NOURAGUES DATASETS####
nourage<-"nouragues.txt"
beginyearfile<-"beginyearseeds 2011 newfecha_NRG.txt"

#### BCI DATASETS####
bcibeginyear<-"bci_beginyearseeds.txt"
bcishort<-"bcishort.txt"
mostabundsppbci<-"BCI-most abundant spp.txt"
bci<-"bci seed traps long.txt"  ##this is the original dataset of BCI


#### key for variable "part" (within bci file) follows:########
## 0 = Reproductive buds (only recorded for selected species) (first recorded on X)
## 1 = Mature fruit
## 2 = Single diaspores
## 3 = "Capsules" (This is a part that vertebrates never eat. Botanically this might be a capsule, pedicel, bract, etc.)
## 4 = Fragments of fruit dropped by vertebrates (we try to record the number of fruit represented by counting pedicels or the points of attachment of fruits to the mother plant)
## 5 = Immature fruit (endosperm of seeds is not filled; falls later during development than #8 below)
## 6 = Perfect and female flowers
## 7 = Fruit with insect emergence hole (only recorded for selected species)
## 8 = Aborted fruits (fall soon after flowering, have a swollen ovule, and often have some flower parts attached)
## 9 = Male flowers



lengthunique = function(x) return(length(unique(x)))


####FUNCTIONS FOR HIERARCHICAL MODELS####

# Hierarchical model of seed fall across many years. Params are a hypermean and hyperSD for date (as yday)
# of peak seed fall, hyperlogmean and hyperlogSD of annual production, each varying across years, then a single
# within-year SD. Error function is dpois, so no other parameter needed. The data table consists of a julian date, yday,
# quantity (seed count), already extracted for a single species. The startyear and endyear are submitted; one peak
# and one peak date will be fitted per year.


#### BASIC MODEL OF SEED PRODUCTION####
# A function to produce a predicted seed count on every day across many years.
# Pass a vector of julian dates x spanning N years, and a vector of N peakdays (vector of julian dates, for instance, the output of function 
# initial.peakdays), N peaks, but a single
# standard deviation, and a single startyear to give the initial year.

model.seedrain=function(data,peakdays,peaks,SD,startyear,endyear,debugmode=FALSE)  
{ 
  nyear=startyear:endyear
  if(length(peakdays)==1) peakdays=rep(peakdays,nyear)
  if(length(peakdays)!=length(peaks)) return(rep(NA, length(x)))
  #if(length(peakdays)!=numberpeaks) return(rep(NA, length(x)))  ## for bimodal species, I should include another parameter with the number of peak
  ## (= twice number of years)
  
  
  inc=which(data$year>=startyear&data$year<=endyear)
  x=data$julian
  pred.quantity=pred.trapcount=rep(NA, length(x)) 
  cldyear=unique(data$calendaryear[inc])
  jan1= tojulian(pst(as.character(nyear),"-01-01"),dateform="%Y-%m-%d")
  peakjulian=jan1+peakdays
  midjulian=diff(peakjulian)/2+peakjulian[-length(nyear)]
  
  savetonextyr=0
  
  for (i in 1:length(nyear))
  {                
    include=which(data$year==nyear[i]) 
    
    if(length(include)>0)
    {
      pred.quantity[include]=peaks[i]*pnorm(q=x[include],mean=peakjulian[i],sd=SD)
      pred.trapcount[include]=c(0,diff(pred.quantity[include]))
      pred.trapcount[include[1]]=savetonextyr + pred.trapcount[include[1]]
      
      # browser()
      savetonextyr=peaks[i]-(pred.quantity[max(include)])
    }
    else savetonextyr=savetonextyr+peaks[i]
  }
  
  
  if(length(which(pred.trapcount<0))>0) browser()
  return(pred.trapcount[inc])
}

#### LIKELIHOOD FUNCTIONS####

# Take one of two hyperparams, an integer to identify which it is, then the vector of both hyperparams (mu and SD), 
# plus the vector of peaks (or peakdays). Then there is a flag to indicate whether it is peak or peakdays, since they have
# different likelihood functions. Return the llikelihood of observing the vector of peaks given the hyper mean and SD.
hyper.llike=function(testparam,whichtest,hyper,P,type='dnorm')
{
  if(whichtest==2 & testparam<=0) return(-Inf)
  hyperparam=arrangeParam.llike(testparam,hyper,whichtest)
  
  if(type=='dnorm') hyper.llike=dnorm(P,mean=hyperparam[1],sd=hyperparam[2],log=TRUE)
  else if(type=='dlnorm') hyper.llike=dlnorm(P,meanlog=hyperparam[1],sdlog=hyperparam[2],log=TRUE)
  
  if(length(which(is.na(hyper.llike)))>0) browser()
  
  return(sum(hyper.llike))
}


# Given a table of data, with julian and seed fall, a start year, and full parameters. There must be a vector of peakdays
# and a vector of peaks, both the same length as the number of years.
full.llike_seedrain=function(data,peakdays,peaks,SD,startyear,endyear,hyperparam=NULL)
{
  inc=which(data$year>=startyear&data$year<=endyear)
  pred=model.seedrain(data=data,peakdays=peakdays,peaks=peaks,SD=SD,startyear=startyear,endyear=endyear)
  model.llike=dpois(round(x=data$quantity[inc],0),lambda=pred,log=TRUE)
  
  if(!is.null(hyperparam))
  {                                 
    hyper.llike1=hyper.llike(testparam=hyperparam[1],whichtest=1,hyper=hyperparam[1:2],P=peakdays,type='dnorm')
    hyper.llike2=hyper.llike(testparam=hyperparam[3],whichtest=1,hyper=hyperparam[3:4],P=peaks,type='dlnorm')
  }
  else hyper.llike1=hyper.llike2=0
  
  if(length(which(is.na(model.llike)))>0) browser()
  return(sum(model.llike)+hyper.llike1+hyper.llike2)
}


# Likelihood function for a single one of the peakdays for seed production. Pass the one peakday to be updated, an integer to indicate
# which is being updated, then full vectors of peaks, peakdays, and the SD, plus data and startyear.
llike.peakday=function(onepeakday,whichpeak,data,peakdays,peaks,SD,hyper,startyear,endyear)
{
  if(onepeakday<=0 | onepeakday>366) return(-Inf)
  
  peakdays[whichpeak]=onepeakday
  model.llike=full.llike_seedrain(data=data,peakdays=peakdays,peaks=peaks,SD=SD,startyear=startyear,endyear=endyear)
  hyper.llike=dnorm(onepeakday,mean=hyper[1],sd=hyper[2],log=TRUE)
  
  return(model.llike+hyper.llike)
}


# Likelihood function for a single one of the peaks for seed production. Pass the one peak to be updated, an integer to indicate
# which is being updated, then full vectors of peaks, peakdays, and the SD, plus data, startyear and endyear.
llike.peak=function(onepeak,whichpeak,data,peakdays,peaks,SD,hyper,startyear,endyear)
{
  if(onepeak<=0) return(-Inf)
  
  peaks[whichpeak]=onepeak
  
  model.llike=full.llike_seedrain(data=data,peakdays=peakdays,peaks=peaks,SD=SD,startyear=startyear,endyear=endyear)
  hyper.llike=dlnorm(onepeak,meanlog=hyper[1],sdlog=hyper[2],log=TRUE)
  
  return(model.llike+hyper.llike)
}


# Likelihood function for the SD. Pass full vectors of peaks, peakdays, and the SD, plus data and startyear.
llike.SD=function(SD,data,peakdays,peaks,startyear,endyear)
{
  if(SD<=0) return(-Inf)
  return(full.llike_seedrain(data=data,peakdays=peakdays,peaks=peaks,SD=SD,startyear=startyear,endyear=endyear))
}

#initial.peak function calculates inital parameters for the peaks of seed production, making the sum of the seeds produced for calender years
initial.peak=function(dataset,startyear,endyear)
{
  nyear=startyear:endyear
  peaksyr<-vector()
  year=vector()
  for (i in 1:length(nyear))
  {
    oneyr<-subset(dataset, year==nyear[i])
    peaksyr[i]= round(sum(oneyr$quantity, na.rm=T),0)
    #year[i]=nyear[i]
  }
  #return (data.frame(year,peaksyr))
  return (peaksyr)
}

#initial.peak.bci function calculates inital parameters for the peaks of seed production, making the sum of the seeds produced for calender years
initial.peak.bci=function(dataset,startyear,endyear)
{
  nyear=startyear:endyear
  peaksyr<-vector()
  year=vector()
  for (i in 1:length(nyear))
  {
    oneyr<-subset(dataset, year==nyear[i])
    peaksyr[i]= round(mean(oneyr$quantity, na.rm=T),0)
    #year[i]=nyear[i]
  }
  #return (data.frame(year,peaksyr))
  return (peaksyr)
}

####Bayesian model####
#model.seedfall.Gibbs2 function calculates by turn the posterior distribution of the species-specific parameters of the seed production model
# SDvector has initial values for standard deviations 1) peakday across years, 2) production, 3) day within year
#peakdays and peaks must be vectors with length=total number of years
#noyear=length(peaks)
#allyear=startyear:(startyear+noyear-1)

#numpeaks = tallies of the total number of peaks of seed production (double number of years for biannual species)

model.seedfall.Gibbs2=function(data,peakdays,peaks,SDvector=c(10,1,35),startyear=1,endyear=10, steps=4000,showstep=200,burnin=1000) 
{
  
  noyear=endyear-startyear+1
  peak=peakday=matrix(nrow=steps,ncol=noyear)
  peakstep=peakdaystep=numeric(noyear)
  peakday[1,]=peakdaystep=peakdays
  peak[1,]=peakstep=peaks 
  
  hyper=data.frame(matrix(nrow=steps,ncol=4))
  colnames(hyper)=c('mu','SD','logmu','logSD')
  hyperstep=numeric(4)         
  
  hyper$mu[1]=hyperstep[1]=mean(peakdays)
  hyper$SD[1]=hyperstep[2]=SDvector[1]
  hyper$logmu[1]=hyperstep[3]=mean(log(peaks))
  hyper$logSD[1]=hyperstep[4]=SDvector[2] 
  
  SD=llike=numeric()
  SD[1]=SDstep=SDvector[3]
  
  i=1
  llike[i]=full.llike_seedrain(data=data,SD=SD[i],peaks=peak[i,],peakdays=peakday[i,],startyear=startyear,endyear=endyear,hyperparam=drp(hyper[i,]))
  # browser()
  
  for(i in 2:steps)
  {
    # Update each of the peakdays in turn
    for(j in 1:noyear)
    {
      useparam=arrangeParam.Gibbs(i,j,peakday)
      
      metropResult=metrop1step(func=llike.peakday,start.param=peakday[i-1,j],scale.param=peakdaystep[j],adjust=1.02,target=0.25,
                               whichpeak=j,data=data,peakdays=useparam,peaks=peak[i-1,],SD=SD[i-1],hyper=drp(hyper[i-1,]),startyear=startyear,endyear=endyear)
      
      peakday[i,j]=metropResult[1]
      peakdaystep[j]=metropResult[2]
    }
    
    # Update each of the peaks in turn
    for(j in 1:noyear)
    {
      useparam=arrangeParam.Gibbs(i,j,peak)
      
      metropResult=metrop1step(func=llike.peak,start.param=peak[i-1,j],scale.param=peakstep[j],adjust=1.02,target=0.25,
                               whichpeak=j,data=data,peakdays=peakday[i,],peaks=useparam,SD=SD[i-1],hyper=drp(hyper[i-1,]),startyear=startyear,endyear=endyear)
      
      peak[i,j]=metropResult[1]
      peakstep[j]=metropResult[2]
    }
    
    # Update the SD
    metropResult=metrop1step(func=llike.SD,start.param=SD[i-1],scale.param=SDstep,adjust=1.02,target=0.25,data=data,
                             peakdays=peakday[i,],peaks=peak[i,],startyear=startyear,endyear=endyear)
    SD[i]=metropResult[1]
    SDstep=metropResult[2]
    
    # Update the two hyper parameters for peakday
    for(j in 1:2)
    {
      useparam=drp(arrangeParam.Gibbs(i,j,hyper))
      
      metropResult=metrop1step(func=hyper.llike,start.param=hyper[i-1,j],scale.param=hyperstep[j],adjust=1.02,target=0.25,
                               whichtest=j,P=peakday[i,],hyper=useparam[1:2],type='dnorm')
      
      hyper[i,j]=metropResult[1]
      hyperstep[j]=metropResult[2]
    }
    
    # Update the two hyper parameters for peak
    for(j in 3:4)
    {
      useparam=drp(arrangeParam.Gibbs(i,j,hyper))
      
      metropResult=metrop1step(func=hyper.llike,start.param=hyper[i-1,j],scale.param=hyperstep[j],adjust=1.02,target=0.25,
                               whichtest=j-2,P=peak[i,],hyper=useparam[3:4],type='dlnorm')
      
      hyper[i,j]=metropResult[1]
      hyperstep[j]=metropResult[2]
    }
    
    
    llike[i]=full.llike_seedrain(data=data,SD=SD[i],peaks=peak[i,],peakdays=peakday[i,],startyear=startyear,endyear=endyear)
    
    if(i%%showstep==0 | i==2) # browser()
      cat("Step ", i, "Mu : ", round(hyper$mu[i],1), "SD: ", round(hyper$SD[i],1), 
          "LogMu : ", round(hyper$logmu[i],1), "LogSD: ", round(hyper$logSD[i],1), 
          "llike=", round(llike[i],1), "\n")
    
  }
  keep=(burnin+1):steps 
  besthyper=colMeans(hyper[keep,])
  CIhyper=apply(hyper[keep,],2,CI)
  bestpeak=colMeans(peak[keep,]) 
  bestpeakday=colMeans(peakday[keep,])
  CIpeak=apply(peak[keep,],2,CI)
  CIpeakday=apply(peakday[keep,],2,CI)
  bestSD=mean(SD[keep])
  model=model.seedrain(data,peaks=bestpeak,peakdays=bestpeakday,SD=bestSD,startyear=startyear,endyear=endyear)
  
  
  return(list(hyper=hyper,besthyper=besthyper,bestpeak=bestpeak,bestpeakday=bestpeakday,CIpeak=CIpeak,CIpeakday=CIpeakday,CIhyper=CIhyper,
              fullpeak=peak,fullpeakday=peakday,SD=SD,bestSD=bestSD,model=model,llike=llike))
}


#### GRAPHIC REPRESENTATION OF MODELS####
# graph.fitseedrain2=function(spdata=spdata,fit=fit,beginyear,startyear=2002,endyear=2010)

graph.fitseedrain2=function(spdata,fit,startyear=2002,endyear=2010)
  {
   nyear=startyear:endyear
   inc=which(spdata$year>=startyear&spdata$year<=endyear)

   x11(xpos=400,ypos=0,height=5.5,width=9)
   
   #### (Rick) I changed to spdata$julian. Use different dataframe for converted or non-converted dates; and year not calendaryear
   plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,max(fit$bestpeak)))
   lines(spdata$julian[inc],fit$model)
   ncldyear=unique(spdata$year[inc])
   ####
   
   jan1=tojulian(pst('01/01/',startyear:endyear))
   fitpeakjulian=jan1+fit$bestpeakday 
   points(fitpeakjulian,fit$bestpeak,col='red',pch=16)
   
   x11(xpos=400,ypos=700,height=5,width=9)
   maxyr=which.max(fit$bestpeak)
   oneyrpred=14*fit$bestpeak[maxyr]*dnorm(1:400,mean=fit$bestpeakday[maxyr],sd=fit$bestSD)
   plot(1:400,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year')
   for(i in 1:length(fit$bestpeak)) 
    { 
     oneyrpred=14*fit$bestpeak[i]*dnorm(1:400,mean=fit$bestpeakday[i],sd=fit$bestSD)
     lines(1:400,oneyrpred,col=i) 
    }
  }
   
### ADING PRIORS TO CALENDAR YEARS 
## Because reproductive years are not concordant with calendar years, 
## we used priors for estimating the starting day of the reproductive year of each species  (= "beginyear").
## Using this prior, all Julian dates were transformed adding a constant value Dj. This way, all reproductive years were centered on the 2nd July. )
## # retranslate.seedfalldate function adds a constant, Dj, to every date in the output of the previous function.  The function needs a prior "beginyear" for working

#original retranslate function
retranslate.seedfalldate=function(fit,beginyear)
{
 Dj=round(ifelse(beginyear<182.5,beginyear,beginyear-366),0)
 fit$hyper$mu=fit$hyper$mu+Dj
 fit$besthyper[1]=fit$besthyper[1]+Dj
 fit$CIhyper[,1]=fit$CIhyper[,1]+Dj

 fit$bestpeakday=fit$bestpeakday+Dj
 fit$CIpeakday=fit$CIpeakday+Dj
 
 fit$fullpeakday=fit$fullpeakday+Dj
 
 return(fit)
}


#create.rep.year function takes the output of "extract.seedfall_onesp2" for translating dates to species-specific reproductive years using a prior coming from
#the "beginyear" priors. Dates are moved a fixed number of days towards the middle of the year (2 July)

## data for Virola micheli: beginyr=218
## vm=extract.seedfall_onesp2(latin = "Virola michelii", file = nourage)
## vm2= create.rep.year(vm,218)
     
create.rep.year=function(data,beginyear)
{
 #### (Rick) I rewrote this, removing many unnecessary steps, returning dataframe exactly like input but with shifted date
  newdata=data
  
  Dj=round(ifelse(beginyear<182.5,beginyear,beginyear-366),0)
  
  newdata$julian=data$julian-Dj
  newdata$Date=fromjulian(newdata$julian,dateform="%Y-%m-%d")
  fulldate=create.fulldate(newdata$Date,format='%Y-%m-%d')
  newdata[,c('year','yday')]=fulldate[,c('year','yday')]
  
  return(newdata)

}   
  

extract.seedfall_onesp2=function(latin='Guettarda acreana',file=nourage)
{
 fulldata=read.delim(file)
 spdata=subset(fulldata,sp==latin & !is.na(quantity))

#### (Rick) The very first census has to be set to zero, otherwise it can't be used. 
 spdata$quantity[1]=0
 
 spdata$julian=as.integer(tojulian(spdata$Date,dateform='%Y-%m-%d') )
 fulldate=create.fulldate(spdata$Date,format='%Y-%m-%d')
 
### (Rick) no need for second dataframe, and I want to keep the date for debugging
 spdata[,c('year','yday')]=fulldate[,c('year','yday')]
 spdata$fecha=NULL
 spdata$calendaryear=NULL
 
 return(spdata)

}



##spshort creates an abbreviation of six letters for each species (useful for Nouragues dataset)
##### (Rick) functions strsplit and left make this much easier
spshort=function(infile=nourage)
{
   fulldata=read.delim(infile,as.is=TRUE)
   spnames=sort(unique(fulldata$sp))
   genuscode=spcode=character()
   
   latin=strsplit(spnames,split=' ')
   for(i in 1:length(latin))
    {
     genuscode[i]=left(latin[[i]][1],3)
     spcode[i]=left(latin[[i]][2],3)
    }
    
   return(pst(genuscode,spcode))
}

##spshort2 creates an abbreviation of six letters for a given name of a species
#spshort2("Virola michelii")

spshort2=function(spname)
{
  genuscode=spcode=character()
  latin=strsplit(spname,split=' ')
  genuscode=left(latin[[1]][1],3)
  spcode=left(latin[[1]][2],3)
  return(pst(genuscode,spcode))
}

##splong reconverts an abbreviation of six letters to the full name of a species
#splong(infile=nourage, shortname="Virmic")
splong=function(infile, shortname)
{
   fulldata=read.delim(infile,as.is=TRUE)
   spnames=sort(unique(fulldata$sp))
   splong=character()
   genuscode=spcode=character()
   
   latin=strsplit(spnames,split=' ')
   for(i in 1:length(latin))
    {
     genuscode[i]=left(latin[[i]][1],3)
     spcode[i]=left(latin[[i]][2],3)
     splong[i]=spnames[i]
    }
    
   spshort=pst(genuscode,spcode)
   spp=data.frame(splong,spshort)
   result=which(spp$spshort==shortname)
   return(spp[result,1])
}   
  

###### APLICATION OF MODELS####

## example of a loop for Nouragues data:
#  nrgresults= hierarchical.NRG.repyr(file=nourage, beginyearfile=beginyearfile,spstart=1,spend=45, steps=10000, burnin=1000,outfile="nrg-results.RData")
#  bciresults= hierarchical.BCI(file=bci,beginyearfile=bcibeginyear,spstart=1,spend=149, steps=10000, burnin=1000,outfile="BCI-results.RData")
##graphics: nrg.graphierspp(file=nourage, fit=results.nrg2, beginyearfile=beginyearfile,spstart=11,spend=20)


#### NOURAGUES DATASET####
##"hierarchical.NRG.repyr" do a loop for all species in Nouragues
## nrg.results2<-hierarchical.NRG.repyr(file=nourage, beginyearfile=beginyearfile,spstart=5,spend=10, steps=4000, burnin=1000,outfile="nrg.results2.RData")

 hierarchical.NRG.repyr<-function(file=nourage, beginyearfile=beginyearfile,spstart=1,spend=10, steps=4000, burnin=1000,outfile)
 
{
   sprange=spstart:spend
   nrgdata=read.delim(file)
   beginyr=read.delim(file=beginyearfile)
   fulldata=merge(nrgdata,beginyr,by="sp", all.x=TRUE)
   spnames=sort(unique(fulldata$sp))
   results=list()
   spshort=spshort()
                          
   for (i in 1:length(sprange))

  {  
     species=as.character(spnames[sprange][i])
     spdata<-extract.seedfall_onesp2(latin =species, file = file)
     bgy= unique(fulldata$beginyr[fulldata$sp==species])
     spdata2<-create.rep.year(data=spdata,beginyear=bgy)
     
     #this is for the selection of the years for model running
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
     results[[spshort2(spname=species)]]= fit
   }
 save(results, file=outfile)
 return(results)
} 



nrg.graphierspp2<-function(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=11,spend=20, filename="Nouragues seed production all spp2.pdf",longnames="total number of seeds per species.txt")
{
   sprange=spstart:spend
   nrgdata=read.delim(file)
   beginyr=read.delim(file=beginyearfile)
   fulldata=merge(nrgdata,beginyr,by="sp", all.x=TRUE)
   totseed=read.delim(file=longnames)
   names(totseed)=c("species","longname","totseed","form","disp","fruit")
   spnames=sort(unique(fulldata$sp))
   spnames2=totseed$longname
   #filename=paste(species,".jpg")   
   #jpeg(filename="Nouragues seed production all spp.jpg",width = 1000, height = 650,quality=100) 
   pdf(file=filename) 
    
   for (i in 1:length(sprange))
  {     
     species=as.character(spnames[sprange][i])
     splong=as.character(spnames2[sprange][i])
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
     fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
     
     par(mfrow=c(2,1), las=1, mar=c(4.5,4,4,1),cex=1)

    #plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,max(fittrans$bestpeak)))
    plot(spdata$julian,spdata$quantity,pch=16,cex=1,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="years",ylab="number of seeds")
    jan=tojulian(pst('01/01/',min(unique(spdata$year)):(max(unique(spdata$year))+1)))
    #axis(side=1, at=jan, labels= c(min(unique(spdata$year)):max(unique(spdata$year))))
    axis(side=1, at=jan, labels= jan)
    axis(2)
    Dj=round(ifelse(bgy<182.5,bgy,bgy-366),0)
    lines(spdata2$julian[inc2]+Dj,fitset$model)  
    newdates=fromjulian(spdata$julian[inc2]+Dj,dateform="%Y-%m-%d")
    fulldate=create.fulldate(newdates,format='%Y-%m-%d')
    startyear=min(unique(fulldate$year))
    endyear=max(unique(fulldate$year))
    jan1=tojulian(pst('01/01/',startv:endv))
    fitpeakjulian=jan1+fittrans$bestpeakday 
    points(fitpeakjulian,fittrans$bestpeak,col='red',pch=16,cex=1)
     mtext(side=3, text=splong, font=3,line=0.75, cex=2) 
   
   maxyr=which.max(fittrans$bestpeak)
   oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(1:400,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
   plot(1:400,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year', lwd=2)
   for(k in 1:length(fittrans$bestpeak)) 
    { 
     oneyrpred=14*fittrans$bestpeak[k]*dnorm(1:400,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
     lines(1:400,oneyrpred,col=k, lwd=2) 
    }
   }
 dev.off()   
}

nrg.graphsmoved.low<-function(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(45,11,40,23,25,26,28,32,11,18,19,4,34,43), filename="Nouragues seed production moved low.pdf",longnames="total number of seeds per species.txt")
{
  
  nrgdata=read.delim(file)
  beginyr=read.delim(file=beginyearfile)
  fulldata=merge(nrgdata,beginyr,by="sp", all.x=TRUE)
  totseed=read.delim(file=longnames)
  names(totseed)=c("species","longname","totseed","form","disp","fruit")
  spnames=sort(unique(fulldata$sp))
  spnames2=totseed$longname
  
  
  pdf(file=filename) 
  for (i in 1:length (spnumber))
    
  {
    species=as.character(spnames[spnumber[i]])
    splong=as.character(spnames2[spnumber[i]])
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
    fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
    
    filename=paste(species,"moved.jpg")   
    #jpeg(filename=filename,width = 1000, height = 650,quality=100)
    
    par(mfrow=c(2,1), las=1, mar= c(4.5,4,4,1),cex=1)
    
    #plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,max(fittrans$bestpeak)))
    plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="years",ylab="number of seeds")
    jan=tojulian(pst('01/01/',min(unique(spdata$year)):(max(unique(spdata$year))+1)))
    #axis(side=1, at=jan, labels= c(min(unique(spdata$year)):max(unique(spdata$year))))
    axis(side=1, at=jan, labels= jan)
    axis(2)
    Dj=round(ifelse(bgy<182.5,bgy,bgy-366),0)
    lines(spdata2$julian[inc2]+Dj,fitset$model)  
    newdates=fromjulian(spdata$julian[inc2]+Dj,dateform="%Y-%m-%d")
    fulldate=create.fulldate(newdates,format='%Y-%m-%d')
    startyear=min(unique(fulldate$year))
    endyear=max(unique(fulldate$year))
    jan1=tojulian(pst('01/01/',startv:endv))
    fitpeakjulian=jan1+fittrans$bestpeakday 
    points(fitpeakjulian,fittrans$bestpeak,col='red',pch=16)
    mtext(side=3, text=splong, font=3,line=0.75, cex=2)
    
    maxyr=which.max(fittrans$bestpeak)
    oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(-185:185,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
    plot(-185:185,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year',axes=FALSE,lwd=2)
    axis(side=1, at=c(-180,-120,-60,1,60,120,180),labels=c("185","245","305","0","60","120","180") )
    axis(side=2)
    for(k in 1:length(fittrans$bestpeak)) 
    { 
      oneyrpred=14*fittrans$bestpeak[k]*dnorm(-185:185,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
      lines(-185:185,oneyrpred,col=k,lwd=2) 
    }
  }
  dev.off()
}

#this graph puts 0 days of the year in the center: it is only valid for HIGH numbers of peakday
nrg.graphsmoved.high<-function(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(5,33,35,37), filename="Nouragues seed production moved high.pdf",longnames="total number of seeds per species.txt")
{
  
  nrgdata=read.delim(file)
  beginyr=read.delim(file=beginyearfile)
  fulldata=merge(nrgdata,beginyr,by="sp", all.x=TRUE)
  totseed=read.delim(file=longnames)
  names(totseed)=c("species","longname","totseed","form","disp","fruit")
  spnames=sort(unique(fulldata$sp))
  spnames2=totseed$longname
  
  pdf(file=filename)
  for (i in 1:length (spnumber))
    
  {
    species=as.character(spnames[spnumber[i]])
    splong=as.character(spnames2[spnumber[i]])
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
    fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
    
    #filename=paste(species,"moved.jpg")   
    #jpeg(filename=filename,width = 1000, height = 650,quality=100)
    
    par(mfrow=c(2,1), las=1,mar= c(4.5,4,4,1),cex=1)
    
    #plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,max(fittrans$bestpeak)))
    plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="years",ylab="number of seeds")
    jan=tojulian(pst('01/01/',min(unique(spdata$year)):(max(unique(spdata$year))+1)))
    #axis(side=1, at=jan, labels= c(min(unique(spdata$year)):max(unique(spdata$year))))
    axis(side=1, at=jan, labels= jan)
    axis(2)
    Dj=round(ifelse(bgy<182.5,bgy,bgy-366),0)
    lines(spdata2$julian[inc2]+Dj,fitset$model)  
    newdates=fromjulian(spdata$julian[inc2]+Dj,dateform="%Y-%m-%d")
    fulldate=create.fulldate(newdates,format='%Y-%m-%d')
    startyear=min(unique(fulldate$year))
    endyear=max(unique(fulldate$year))
    jan1=tojulian(pst('01/01/',startv:endv))
    fitpeakjulian=jan1+fittrans$bestpeakday 
    points(fitpeakjulian,fittrans$bestpeak,col='red',pch=16)
    mtext(side=3, text=splong, font=3,line=0.75, cex=2)
    
    maxyr=which.max(fittrans$bestpeak)
    oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(185:550,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
    plot(185:550,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year',axes=FALSE,lwd=2)
    axis(side=1, at=c(185,245,305,366,425,485,545),labels=c("185","245","305","0","60","120","180"))
    axis(side=2)
    for(k in 1:length(fittrans$bestpeak)) 
    { 
      oneyrpred=14*fittrans$bestpeak[k]*dnorm(185:550,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
      lines(185:550,oneyrpred,col=k,lwd=2) 
    }
  }
  dev.off()
}

nrg.graphsmoved.other<-function(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(1), lowvalue=100, highvalue=465,filename="Acacia tenuifolia moved.pdf")
{
  
  nrgdata=read.delim(file)
  beginyr=read.delim(file=beginyearfile)
  fulldata=merge(nrgdata,beginyr,by="sp", all.x=TRUE)
  spnames=sort(unique(fulldata$sp))
  pdf(file=filename)
  for (i in 1:length (spnumber))
    
  {
    species=as.character(spnames[spnumber[i]])
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
    fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
    
    #filename=paste(species,"moved.jpg")   
    #jpeg(filename=filename,width = 1000, height = 650,quality=100)
    
    par(mfrow=c(2,1), las=1,mar= c(4.5,4,4,1),cex=1)
    
    #plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,max(fittrans$bestpeak)))
    plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="years",ylab="number of seeds")
    jan=tojulian(pst('01/01/',min(unique(spdata$year)):(max(unique(spdata$year))+1)))
    #axis(side=1, at=jan, labels= c(min(unique(spdata$year)):max(unique(spdata$year))))
    axis(side=1, at=jan, labels= jan)
    axis(2)
    Dj=round(ifelse(bgy<182.5,bgy,bgy-366),0)
    lines(spdata2$julian[inc2]+Dj,fitset$model)  
    newdates=fromjulian(spdata$julian[inc2]+Dj,dateform="%Y-%m-%d")
    fulldate=create.fulldate(newdates,format='%Y-%m-%d')
    startyear=min(unique(fulldate$year))
    endyear=max(unique(fulldate$year))
    jan1=tojulian(pst('01/01/',startv:endv))
    fitpeakjulian=jan1+fittrans$bestpeakday 
    points(fitpeakjulian,fittrans$bestpeak,col='red',pch=16)
    mtext(side=3, text=species, font=3,line=0.75, cex=2)
    
    maxyr=which.max(fittrans$bestpeak)
    oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(lowvalue:highvalue,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
    plot(lowvalue:highvalue,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year',axes=FALSE,lwd=2)
    axis(side=1, at=c(seq(lowvalue, highvalue, 60)),labels=TRUE)
    axis(side=2)
    for(k in 1:length(fittrans$bestpeak)) 
    { 
      oneyrpred=14*fittrans$bestpeak[k]*dnorm(lowvalue:highvalue,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
      lines(lowvalue:highvalue,oneyrpred,col=k,lwd=2) 
    }
  }
  dev.off()
}



#### BCI dataset#####

#getseedsbci200traps function subsets the dataset of BCI for fruits and single seeds and the original 200 seed traps (excluding those of gaps)
getseedsbci200traps=function(part,species,trap) {
  inc = ((part==1 & species!="ALSB") | part==2) & trap<=200
  return(inc)
} 

#mostabundspp.bci function selects for the most abundant species according to the criterion in Wright & Calderon: seeds or fruits captured in 10 or more traps in any single year 


mostabundspp.bci<-function(file=bci)
{ 
  bci1=read.delim(file)
  inc<-getseedsbci200traps(bci1$part,bci1$sp,bci1$trap)
  bcifruits<-bci1[inc,] 
  bcifruits$fecha2<-strptime(as.character(bcifruits$fecha), format="%Y-%m-%d")                
  bcifruits$Year<-bcifruits$fecha2$year+1900
  
  trapfreq=aggregate(data.frame(ntraps=bcifruits$trap),by=list(sp=bcifruits$sp,year=bcifruits$Year),lengthunique)
  inc=which(trapfreq$ntraps>=10)
  tr<-trapfreq[inc,]
  unk=which(tr$sp=="UNK?")
  maso=which(tr$sp=="MASO")
  tet1=which(tr$sp=="TET1")
  
  tr<-tr[-c(unk,maso,tet1),]
  spma<-unique(tr$sp)
  return(spma)
} 


##bci.results10<-hierarchical.BCI(file=bci,beginyearfile=bcibeginyear, spstart=61,spend=66, steps=10000,burnin=1000,outfile="bci.results10.RData")
hierarchical.BCI<-function(file=bci,beginyearfile=bcibeginyear,spstart=1,spend=10, steps=2000, burnin=1000,outfile)
 
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
   sprange=spstart:spend
   spnames=as.vector(sort(unique(bcishort$sp)))
   #spnameshort=spnames[c(157,149,81)]
   
   results14=list()          
   startv=vector()
   for (i in 1:length(sprange))
                       
  {     
     species=as.character(spnames[sprange][i])
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
     results14[[spshort2(spname=species)]]= fit

     #save(fit,file=paste(spnames[i],".RData"))
     results14[[species]]= fit
   }
 save(results14, file=outfile)
 return(results14)
}

 
#bcicens function sums the quantity values of the traps of the same census 

bcicens<-function(infile)
 {                                       
  bci=infile
  #cens<-c(3:61,63:1287)
  cens<-unique(bci$census)
  fecha=census=quantity=numeric()
 
 for (i in 1:length(cens))
 {
    onecs<-bci[bci$census==cens[i],]
    quantity <-sum(onecs$quantity)
    fecha<-unique(onecs$fecha)[1]    
    census<-cens[i]
    results<-data.frame(fecha,census,quantity)
    if (i==1) allresults=results else allresults=rbind(allresults,results)
 }   
 allresults$sp=unique(bci$sp)  
 return(allresults) 
}  

#bciwithzeros adds all the censuses to the dataset for each single species
bciwithzeros<-function(infile)
 {                                        
 dataset=infile
 cens<-read.delim (file="bci-census list.txt",header=T,as.is=T)
 #cens=cens[1:1284,] #this is for the old dataset from BCI
 #cens$julian=as.integer(tojulian(cens$fecha, dateform="%Y%m%d"))
 cens$fecha=as.character(cens$fecha)
 dataset$fecha=as.character(dataset$fecha) 
 allcens<-merge(dataset, cens, by="census",all.y=T)
 allcens$quantity<-ifelse(is.na(allcens$quantity),0,allcens$quantity)
 allcens$fecha=allcens$fecha.y   
 results<-allcens[,c(6,1,3)]
 #results$Date<-strptime(as.character(results$fecha), format="%Y%m%d")   #this is for the old version of data from BCI (where date is a expressed as yyyymmdd)
 results$Date<-strptime(as.character(results$fecha), format="%Y-%m-%d") 
 results$year<-(results$Date$year)+1900
 return(results)
#write.table(results,outfile,row.names=F,sep="\t")

} 

#this function prepares data for hierarchical models of reproductive year for a given dataset
extract.seedfall_bci=function(dataset) #dataset=tri
{
 spdata=dataset
 spdata$quantity[1]=0
 spdata$julian=as.integer(tojulian(spdata$Date,dateform='%Y-%m-%d') )
 fulldate=create.fulldate(spdata$Date,format='%Y-%m-%d')
 spdata[,c('year','yday')]=fulldate[,c('year','yday')]
 
 return(spdata)
}

 
 
bci.graphierspp<-function(file=bci, fit=bciresults2, beginyearfile=bcibeginyear,spstart=101,spend=166)
{
   bci1=read.delim(file)  
   beginyr=read.delim(file=beginyearfile) 
   inc<-getseedsbci200traps(bci1$part,bci1$sp,bci1$trap)
   bcifruits<-bci1[inc,]  
   maspp=read.delim(file=mostabundsppbci) #this line selects the most abundant species according to the criterion in Wright & Calderon: seeds or fruits captured in 10 or more traps in any single year 
   ma<-data.frame(maspp)
   bcimaspp<-merge(bcifruits,ma,by.x="sp",by.y="spma", all.y=T)
   bcishort<-merge(bcimaspp,beginyr,by.x="sp",by.y="sp", all.x=T)
   inc2=which(is.na(bcishort$beginyr)==T)
   bcishort=bcishort[-c(inc2),]
   sprange=spstart:spend
   spnames=as.vector(sort(unique(bcishort$sp)))
   #spnameshort=spnames[c(157,149,81)]
   
   results=list()          
   startv=vector()
   for (i in 1:length(sprange))
                       
  {     
     species=as.character(spnames[sprange][i])
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
     fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
     
     filename=paste(species,"2.jpg")   
     jpeg(filename=filename,width = 1200, height = 850,quality=100)

     
     par(mfrow=c(2,1))
     inc=which(spdata$year>=startv&spdata$year<=endv)
     #plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,max(fittrans$bestpeak)))
     plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="years",ylab="number of seeds")
    #plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,30),axes=F,xlab="years",ylab="number of seeds")
    #jan=tojulian(pst('01/01/',startyear:endyear))
    jan=tojulian(pst('01/01/',min(unique(spdata$year)):(max(unique(spdata$year))+1)))
    #axis(side=1, at=jan, labels= c(min(unique(spdata$year)):max(unique(spdata$year))))
    axis(side=1, at=jan, labels= jan)
    axis(2)
    #model=model.seedrain(data=spdata,peaks=fittrans$bestpeak,peakdays=fittrans$bestpeakday,SD=fittrans$bestSD,startyear=startyear,endyear=endyeard)
    lines(spdata2$julian[inc2]+Dj,fitset$model)  
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

#this graph puts 0 days of the year in the center: it is only valid for LOW numbers of peakday
bci.graphsmoved.low=function(file=bcishort, fit=results14,species="HAM2")        
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
     fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
     
     filename=paste(species,"moved.pdf")   
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
   oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(-185:185,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
   plot(-185:185,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year',axes=FALSE)
   axis(side=1, at=c(-180,-120,-60,1,60,120,180),labels=c("185","245","305","0","60","120","180") )
   axis(side=2)
   for(k in 1:length(fittrans$bestpeak)) 
    { 
     oneyrpred=14*fittrans$bestpeak[k]*dnorm(-185:185,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
     lines(-185:185,oneyrpred,col=k) 
    }
 dev.off()
}


#this graph puts 0 days of the year in the center: it is only valid for HIGH numbers of peakday
bci.graphsmoved.high=function(file=bcishort, fit=results14,species="HAM2")         
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
     fittrans= retranslate.seedfalldate(fit=fitset,beginyear=bgy)
     
     filename=paste(species,"moved.pdf")   
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
   oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(185:550,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
   plot(185:550,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year',axes=FALSE)
   axis(side=1, at=c(185,245,305,366,425,485,545),labels=c("185","245","305","0","60","120","180"))
   axis(side=2)
   for(k in 1:length(fittrans$bestpeak)) 
    { 
     oneyrpred=14*fittrans$bestpeak[k]*dnorm(185:550,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
     lines(185:550,oneyrpred,col=k) 
    }
 dev.off()
}



