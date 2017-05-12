
nourage<-"seeddatesp all spp 2011 with zeros NA and data grouped by censuses.txt"
#attach_if_needed("C:\\Brunoy\\Base Datos Nouragues\\Metadata Joe Wright\\modeling max likelihood\\CTFSRPackage.rdata")
library(date)

# Hierarchical model of seed fall across many years. Params are a hypermean and hyperSD for date (as yday)
# of peak seed fall, hyperlogmean and hyperlogSD of annual production, each varying across years, then a single
# within-year SD. Error function is dpois, so no other parameter needed. The data table consists of a julian date, yday,
# quantity (seed count), already extracted for a single species. The startyear and endyear are submitted; one peak
# and one peak date will be fitted per year.

# Sample for 2 species
# spdataGa=extract.seedfall_onesp(latin = "Guettarda acreana", file = nourage)
# fitGa=model.seedfall.Gibbs(data=spdataGa,start=c(225,10,3.7,1.2,25),startyear=2001,endyear=2010,steps=2000,show=10,burn=1000)
# graph.fitseedrain(spdata=spdataGa,fit=fitGa)
# attach("C:\\Brunoy\\Base Datos Nouragues\\Metadata Joe Wright\\modeling max likelihood\\CTFSRPackage.rdata")
# fitTp=model.seedfall.Gibbs(data=spdata,start=c(225,10,3.7,1.2,25),startyear=2001,endyear=2010,steps=2000,show=10,burn=1000)
# graph.fitseedrain(spdata=spdataTp,fit=fitTp)

graph.fitseedrain=function(spdata=spdata,fit=fit,startyear=2001,endyear=2010)
  {
   graphics.off()
   x11(xpos=200,ypos=0,height=5.5,width=9)
   plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fit$bestpeak)))
   lines(spdata$julian,fit$model)
   jan1=tojulian(pst('01/01/',startyear:endyear))
   fitpeakjulian=jan1+fit$bestpeakday
   points(fitpeakjulian,fit$bestpeak,col='red',pch=16)
   
   x11(xpos=200,ypos=700,height=5,width=9)
   maxyr=which.max(fit$bestpeak)
   oneyrpred=14*fit$bestpeak[maxyr]*dnorm(1:400,mean=fit$bestpeakday[maxyr],sd=fit$bestSD)
   plot(1:400,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year')
   for(i in 1:length(fit$bestpeak)) 
    { 
     oneyrpred=14*fit$bestpeak[i]*dnorm(1:400,mean=fit$bestpeakday[i],sd=fit$bestSD)
     lines(1:400,oneyrpred,col=i) 
    }
  }
  
  
model.seedfall.Gibbs=function(data,start=c(180,10,5,3,40),startyear=2001,endyear=2011,steps=4000,showstep=200,burnin=1000)
{
 
 noyear=endyear-startyear+1
 peak=peakday=matrix(nrow=steps,ncol=noyear)
 peakstep=peakdaystep=numeric(noyear)
 peakday[1,]=peakdaystep=rep(start[1],noyear)
 peak[1,]=peakstep=rep(exp(start[3]),noyear)
 
 hyper=data.frame(matrix(nrow=steps,ncol=4))
 colnames(hyper)=c('mu','SD','logmu','logSD')
 hyperstep=numeric(4)
 
 hyper$mu[1]=hyperstep[1]=start[1]
 hyper$SD[1]=hyperstep[2]=start[2]
 hyper$logmu[1]=hyperstep[3]=start[3]
 hyper$logSD[1]=hyperstep[4]=start[4]

 SD=llike=numeric()
 SD[1]=SDstep=start[5]
 
 i=1
 llike[i]=full.llike_seedrain(data=data,SD=SD[i],peaks=peak[i,],peakdays=peakday[i,],startyear=startyear,hyperparam=drp(hyper[i,]))
 # browser()
 
 for(i in 2:steps)
  {
   # Update each of the peakdays in turn
   for(j in 1:noyear)
    {
     useparam=arrangeParam.Gibbs(i,j,peakday)

     metropResult=metrop1step(func=llike.peakday,start.param=peakday[i-1,j],scale.param=peakdaystep[j],adjust=1.02,target=0.25,
                              whichpeak=j,data=data,peakdays=useparam,peaks=peak[i-1,],SD=SD[i-1],hyper=drp(hyper[i-1,]),startyear=startyear)
                              
     peakday[i,j]=metropResult[1]
     peakdaystep[j]=metropResult[2]
    }
   
   # Update each of the peaks in turn
   for(j in 1:noyear)
    {
     useparam=arrangeParam.Gibbs(i,j,peak)

     metropResult=metrop1step(func=llike.peak,start.param=peak[i-1,j],scale.param=peakstep[j],adjust=1.02,target=0.25,
                              whichpeak=j,data=data,peakdays=peakday[i,],peaks=useparam,SD=SD[i-1],hyper=drp(hyper[i-1,]),startyear=startyear)
                              
     peak[i,j]=metropResult[1]
     peakstep[j]=metropResult[2]
    }
  
   # Update the SD
   metropResult=metrop1step(func=llike.SD,start.param=SD[i-1],scale.param=SDstep,adjust=1.02,target=0.25,data=data,
                            peakdays=peakday[i,],peaks=peak[i,],startyear=startyear)
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
   
   
   llike[i]=full.llike_seedrain(data=data,SD=SD[i],peaks=peak[i,],peakdays=peakday[i,],startyear=startyear)
                       
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
 model=model.seedrain(data$julian,peaks=bestpeak,peakdays=bestpeakday,SD=bestSD,startyear=startyear)
 
 return(list(hyper=hyper,besthyper=besthyper,bestpeak=bestpeak,bestpeakday=bestpeakday,CIpeak=CIpeak,CIpeakday=CIpeakday,CIhyper=CIhyper,
             fullpeak=peak,fullpeakday=peakday,SD=SD,bestSD=bestSD,model=model,llike=llike))
}


# Likelihood function for a single one of the peaks for seed production. Pass the one peak to be updated, an integer to indicate
# which is being updated, then full vectors of peaks, peakdays, and the SD, plus data and startyear.
llike.peak=function(onepeak,whichpeak,data,peakdays,peaks,SD,hyper,startyear)
{
 if(onepeak<=0) return(-Inf)
 
 peaks[whichpeak]=onepeak
 
 model.llike=full.llike_seedrain(data=data,peakdays=peakdays,peaks=peaks,SD=SD,startyear=startyear)
 hyper.llike=dlnorm(onepeak,meanlog=hyper[1],sdlog=hyper[2],log=TRUE)
 
 return(model.llike+hyper.llike)
}

# Likelihood function for a single one of the peakdays for seed production. Pass the one peakday to be updated, an integer to indicate
# which is being updated, then full vectors of peaks, peakdays, and the SD, plus data and startyear.
llike.peakday=function(onepeakday,whichpeak,data,peakdays,peaks,SD,hyper,startyear)
{
 if(onepeakday<=0 | onepeakday>366) return(-Inf)
 
 peakdays[whichpeak]=onepeakday
 model.llike=full.llike_seedrain(data=data,peakdays=peakdays,peaks=peaks,SD=SD,startyear=startyear)
 hyper.llike=dnorm(onepeakday,mean=hyper[1],sd=hyper[2],log=TRUE)

 return(model.llike+hyper.llike)
}


# Likelihood function for the SD. Pass full vectors of peaks, peakdays, and the SD, plus data and startyear.
llike.SD=function(SD,data,peakdays,peaks,startyear)
{
 if(SD<=0) return(-Inf)
 return(full.llike_seedrain(data=data,peakdays=peakdays,peaks=peaks,SD=SD,startyear=startyear))
}

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
full.llike_seedrain=function(data,peakdays,peaks,SD,startyear,hyperparam=NULL)
{
 pred=model.seedrain(x=data$julian,peakdays=peakdays,peaks=peaks,SD=SD,startyear=startyear)
 model.llike=dpois(x=data$quantity,lambda=pred,log=TRUE)
 
 if(!is.null(hyperparam))
  {                                 
   hyper.llike1=hyper.llike(testparam=hyperparam[1],whichtest=1,hyper=hyperparam[1:2],P=peakdays,type='dnorm')
   hyper.llike2=hyper.llike(testparam=hyperparam[3],whichtest=1,hyper=hyperparam[3:4],P=peaks,type='dlnorm')
  }
 else hyper.llike1=hyper.llike2=0
 
 if(length(which(is.na(model.llike)))>0) browser()
 return(sum(model.llike)+hyper.llike1+hyper.llike2)
}

# A function to produce a predicted seed count on every day across many years.
# Pass a vector of julian dates x spanning N years, and a vector of N peakdays, N peaks, but a single
# standard deviation, and a single startyear to give the initial year.

model.seedrain=function(x,peakdays,peaks,SD,startyear,debugmode=FALSE)

{
 noyear=length(peaks)
 if(length(peakdays)==1) peakdays=rep(peakdays,noyear)
 if(length(peakdays)!=length(peaks)) return(rep(NA, length(x)))

 allyear=startyear:(startyear+noyear-1)
 jan1= tojulian(pst(as.character(allyear),"-01-01"),dateform="%Y-%m-%d")
 peakjulian=jan1+peakdays
 # if(min(x)<peakjulian[1]-250) return(rep(NA, length(x)))
 # if(max(x)>peakjulian[noyear]+250) return(rep(NA, length(x)))
 midjulian=diff(peakjulian)/2+peakjulian[-noyear]
 savetonextyr=0
 pred.quantity=pred.trapcount=rep(NA, length(x))

 for (i in 1: noyear)
  {
    if(i==1) startjulian=min(x)
    else startjulian=midjulian[i-1]
    if (i<noyear) endjulian=midjulian[i]
    else endjulian=max(x)+1
  
    include=which(x>=startjulian&x<endjulian)
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
 return(pred.trapcount)
}



extract.seedfall_onesp=function(latin='Guettarda acreana',file=nourage)
{
 fulldata=read.delim(file)
 spdata=subset(fulldata,sp==latin & !is.na(quantity))
 
 spdata$julian=as.integer(tojulian(spdata$Date,dateform='%Y-%m-%d') )
 fulldate=create.fulldate(spdata$Date,format='%Y-%m-%d')
 result=subset(spdata,select=c('quantity','julian'))
 result$year=fulldate$year
 result$yday=fulldate$yday
 result$quantity=round(result$quantity,0)
 
 return(result)
}
