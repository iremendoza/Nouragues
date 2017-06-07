lengthunique <- function(x) return(length(unique(x)))
lengthisna <- function(x) return(length(which(is.na(x))))
cvLognormal <- function(logSD,...) sqrt(exp((logSD^2))-1)


runningmean = function(dataset, k = 3) #this function calculates the running mean for a k number of observations; 
  #dataset must be a vector including the variable of interest  
{
  N<-length(dataset)
  total<-N-k+1
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))
  
  counter=1
  val=vector()
  
  for (i in 1:total)
  {
    val[i]<-mean(dataset[i:(k+i-1)])
    counter=counter+1
  }  
  return(val)
  
}




###Getting the correct number of estimated number of seeds per month ###############

####1- First, extract the original parameters of peak, peakday, and SD for each species.#### 
# Which years were included?
#the function includes the possibility of correcting peakday dates with a factor Dj
#file=nourage; beginyearfile=beginyearfile; spstart=1; spend=45; fit=results; dj=FALSE
#parameters=parametersyr(file=nourage, beginyearfile=beginyearfile, spstart=1, spend=45, fit=results, dj=TRUE)
#pp=aggregate(data.frame(peakday=parameters$peakday,  CI2peakday=parameters$CI2peakday,  CI97peakday=parameters$CI97peakday,peak=parameters$peak,  CI2peak=parameters$CI2peak,  CI97peak=parameters$CI97peak), by=list(year=parameters$cycle,sp=parameters$sp),mean)
#hyper=aggregate(data.frame(peakdaymu=parameters$peakdaymu,  CI2hypermu=parameters$CI2hypermu,  CI97hypermu=parameters$CI97hypermu,peakdaysd=parameters$peakdaysd,  CI2hyperSD=parameters$CI2hyperSD,  CI97hyperSD=parameters$CI97hyperSD,peaklogmu=parameters$peaklogmu,CI2hyperlogmu=parameters$CI2hyperlogmu,  CI97hyperlogmu=parameters$CI97hyperlogmu,peaklogsd=parameters$peaklogsd,  CI2hyperlogSD=parameters$CI2hyperlogSD,CI97hyperlogSD=parameters$CI97hyperlogSD), by=list(sp=parameters$sp),unique)  
#hyper2=data.frame(sp=hyper$sp,peakdaymu=round(hyper$peakdaymu,2),  CIhypermu=paste(round(hyper$CI2hypermu,2),"-", round(hyper$CI97hypermu,2)),peakdaysd=round(hyper$peakdaysd,2),CIhyperSD=paste(round(hyper$CI2hyperSD,2),"-",round(hyper$CI97hyperSD,2)), peaklogmu=round(hyper$peaklogmu,2),CIhyperlogmu=paste(round(hyper$CI2hyperlogmu,2),"-",round(hyper$CI97hyperlogmu)),peaklogsd=round(hyper$peaklogsd,2),CIhyperlogSD=paste(round(hyper$CI2hyperlogSD,2),"-",round(hyper$CI97hyperlogSD,2)))  
#yrmin=aggregate(data.frame(min=parameters$year), by=list(year=parameters$cycle,sp=parameters$sp),min)
#yrmax=aggregate(data.frame(max=parameters$year), by=list(year=parameters$cycle,sp=parameters$sp),max)
#pp2=data.frame(sp=pp$sp, cycle=pp$year,years=paste(yrmin$min,"-",yrmax$max), peakdays=round(pp$peakday,2), CIpeakday= paste(round(pp$CI2peakday,2),"-" ,round(pp$CI97peakday,2)), peak=round(pp$peak,2), CIpeak=paste(round(pp$CI2peak,2),"-", round(pp$CI97peak,2)))

months = list(jan = 1:31, feb=32:59, mar=60:90, apr=91:120,
            may=121:151, jun=152:181, jul=182:212, aug=213:243, sep=244:273, oct=274:304, nov=305:334, dec=335:365)
rmonths=unlist(months)

#peak=ifelse(pp2$peakdays<0,pp2$peakdays+365,pp2$peakdays)
#peak2=ifelse(peak>365,peak-365,peak)
#pp2$month=names(rmonths[trunc(peak2)])  

parametersyr = function(file = nourage, beginyearfile = beginyearfile, spstart = 1, spend = 45, fit = results, dj = FALSE)
  
{
  
  sprange = spstart:spend
  nrgdata = read.delim(file)
  beginyr = read.delim(file = beginyearfile)
  fulldata = merge(nrgdata,beginyr,by="sp", all.x=TRUE)
  spnames=sort(unique(fulldata$sp))
  results=list()
  spshort = spshort()
  
  for (i in 1:length(sprange))
    
  {  
    species = as.character(spnames[sprange][i])
    bgy = beginyr$beginyr[beginyr$sp==species]
    spdata <- extract.seedfall_onesp2(latin = species, file = file)
    spdata2 <- create.rep.year(data = spdata, beginyear = bgy)
    
    #this is for the selection of the years for model running
    nyear = unique(spdata2$year)    
    maxyday = numfiles = dif = numeric()
    for (j in 1:length(nyear))
    { 
      oneyr = subset(spdata2,spdata2$year == nyear[j])
      maxyday[j] <- max(oneyr$yday)
      dif[j] <- max(oneyr$yday)-min(oneyr$yday)
      numfiles[j] = dim(oneyr)[1]  
    }
    
    yearselection = data.frame(maxyday,dif,numfiles)
    numyear = dim(yearselection)[1]
    startv = endv = numeric()
    startv = ifelse(yearselection$maxyday[1] >= 320 & dif[1]>300, nyear[1], nyear[2])
    endv = ifelse(yearselection$maxyday[numyear] >= 320 & dif[numyear] > 300, nyear[numyear], nyear[numyear-1])
    inc2 = which(spdata2$year>= startv & spdata2$year <= endv)
    
    spsh = spshort2(spname=species)
    fitset = fit[[spsh]] 
    #fittrans= retranslate.seedfalldate2(data=spdata,fit=fitset,beginyear=bgy,startyear=startyear,endyear=endyear)
    spdata3 = data.frame(spdata2[inc2,], model = fitset$model)
    Dj = round(ifelse(bgy<182.5, bgy, bgy-366), 0)
    cycle=numeric()
    
    for (k in 1:nrow(spdata3))
    {
      cycle[k]=as.numeric(which(spdata3$year[k]==unique(spdata3$year)))
      
    }
    spdata3$cycle = cycle
    if (dj==TRUE) spdata3$julianew = spdata3$julian+Dj 
    if(dj==FALSE) spdata3$julianew=spdata3$julian    #### I have included this for corrections
    fulldate = create.fulldate(fromjulian(spdata3$julianew,dateform='%Y-%m-%d'),format='%Y-%m-%d')
    spdata3$yearn=fulldate$year
    spdata3$monthn=fulldate$month
    spdata3$ydayn=fulldate$yday
    
    newmod=data.frame(spdata3[,1:5], Date=fromjulian(spdata3$julianew,dateform='%Y-%m-%d'),julian=spdata3$julianew, year=spdata3$yearn, month=spdata3$monthn, yday=spdata3$ydayn, model=spdata3$model, cycle=spdata3$cycle)
    first=seq(1,2*length(unique(spdata3$cycle)),2)
    last=seq(2,2*length(unique(spdata3$cycle)),2)
    peakday=peak=CI2peak=CI97peak=CI2peakday=CI97peakday=numeric()
    for (l in 1:nrow(newmod))
    {
      if (dj== TRUE) peakday[l]=fitset$bestpeakday[newmod$cycle[l]]+Dj 
      if (dj ==FALSE) peakday[l]=fitset$bestpeakday[newmod$cycle[l]]
      peak[l]=fitset$bestpeak[newmod$cycle[l]]
      CI2peak[l]= fitset$CIpeak[first[newmod$cycle[l]]]
      CI97peak[l]= fitset$CIpeak[last[newmod$cycle[l]]]
      if (dj== TRUE) CI2peakday[l]= fitset$CIpeakday[first[newmod$cycle[l]]]+Dj
      if (dj== FALSE) CI2peakday[l]= fitset$CIpeakday[first[newmod$cycle[l]]]
      if (dj== TRUE) CI97peakday[l]= fitset$CIpeakday[last[newmod$cycle[l]]]+Dj
      if (dj== FALSE) CI97peakday[l]= fitset$CIpeakday[last[newmod$cycle[l]]]
      
    }
    if (dj==TRUE) {results=data.frame(newmod, peakday, CI2peakday,CI97peakday, peak, CI2peak, CI97peak, sd=fitset$bestSD, peakdaymu=fitset$besthyper[[1]]+Dj,CI2hypermu=fitset$CIhyper[[1]]+Dj,CI97hypermu=fitset$CIhyper[[2]]+Dj,peakdaysd=fitset$besthyper[[2]], CI2hyperSD=fitset$CIhyper[[3]],CI97hyperSD=fitset$CIhyper[[4]],
                                      peaklogmu=exp(fitset$besthyper[[3]]),CI2hyperlogmu=exp(fitset$CIhyper[[5]]),CI97hyperlogmu=exp(fitset$CIhyper[[6]]),peaklogsd=exp(fitset$besthyper[[4]]),CI2hyperlogSD=exp(fitset$CIhyper[[7]]),CI97hyperlogSD=exp(fitset$CIhyper[[8]]) )
    }                
    
    if (dj==FALSE) {results=data.frame(newmod, peakday, CI2peakday,CI97peakday, peak, CI2peak, CI97peak, sd=fitset$bestSD, peakdaymu=fitset$besthyper[[1]],CI2hypermu=fitset$CIhyper[[1]],CI97hypermu=fitset$CIhyper[[2]],peakdaysd=fitset$besthyper[[2]], CI2hyperSD=fitset$CIhyper[[3]],CI97hyperSD=fitset$CIhyper[[4]],
                                       peaklogmu= exp(fitset$besthyper[[3]]),CI2hyperlogmu=exp(fitset$CIhyper[[5]]),CI97hyperlogmu=exp(fitset$CIhyper[[6]]),peaklogsd=exp(fitset$besthyper[[4]]),CI2hyperlogSD=exp(fitset$CIhyper[[7]]),CI97hyperlogSD=exp(fitset$CIhyper[[8]]) )
    }                
    
    
    if (i==1) allresults=results else allresults=rbind(allresults,results)
  }
  
  return(allresults)
}


####2- Second, calculate the accumulated amounts for each month of the dataset####


#dd=monthlyvalues(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", beginyearfile=beginyearfile,spstart=1,spend=45)

## write.table(dd, file="monthly values of seed model.txt", sep="\t", row.names=F)
monthlyvalues=function(file="Nouragues model all spp.txt", fileno="NRG model all spp without Dj.txt",beginyearfile=beginyearfile,spstart=1,spend=45){
  
  sprange=spstart:spend
  nrgdata=read.delim(file)
  nrgdatano=read.delim(fileno)
  spnames=sort(unique(nrgdata$sp))
  beginyr=read.delim(file=beginyearfile)
  
  for (i in 1:length(sprange))
    
  {  
    species=as.character(spnames[sprange][i])
    onesp=nrgdata[nrgdata$sp==species,]
    onespno=nrgdatano[nrgdatano$sp==species,]
    startyear=min(unique(onesp$year))
    endyear=max(unique(onesp$year))
    startyearno=min(unique(onespno$year))
    endyearno=max(unique(onespno$year))
    bgy=beginyr$beginyr[beginyr$sp==species]
    Dj=round(ifelse(bgy<182.5,bgy,bgy-366),0)
    
    oneyr=numeric()  
    #for(k in 1:length((startyear+1):(endyear-1)))#this is for calculating the julian of the first day of each month
    for(k in 1:length((startyear):(endyear)))
    {
      oneyr=rep(startyear+(k-1),12)  
      if (k ==1) yrvector=oneyr else yrvector=c(yrvector, oneyr)
    }
    
    firstmonth=as.numeric(unique(onesp$month)[1])
    lastmonth=onesp$month[onesp$julian==sort(unique(onesp$julian))[length(unique(onesp$julian))]]
    #day1=c(tojulian(pst(firstmonth:12,'/01/',startyear)), c(tojulian(pst(rep(1:12, length((startyear+1):(endyear-1))),"/01/",yrvector))),c(tojulian(pst(1:lastmonth,'/01/',endyear))))
    day1=c(tojulian(pst(rep(1:12, length(startyear:endyear)),"/15/",yrvector)))
    day1min=which(day1==tojulian(paste(firstmonth,"/15/",startyear),dateform = "%m /%d/ %Y"))
    day1max=which(day1==as.numeric(tojulian(paste(lastmonth,"/15/",endyear),dateform = "%m /%d/ %Y")))
    day1=day1[day1min:day1max]
    day1=day1-Dj
    day1=day1[which(day1>=onespno$julian[1])]
    #which(day1>onespno$julian[length(onesp$julian)])
    #datanew=data.frame(julian=day1)
    #datanew$year=create.fulldate(fromjulian(day1, dateform="%Y-%m-%d"),format="%Y-%m-%d")$year
    #mm=model.seedrain(data=datanew, peaks=unique(onesp$peak), peakdays=unique(onesp$peakday), SD=unique(onesp$sd), startyear=min(unique(onesp$year)),endyear=max(unique(onesp$year)))
    #data1=create.fulldate(fromjulian(day1, dateform="%Y-%m-%d"),format="%Y-%m-%d")$year
    datamonth=seq(1,length(unique(onespno$cycle))*12,12)
    
    jan1=tojulian(pst('01/01/',startyearno:endyearno))   
    meanvector=numeric()
    for (p in 1:length(unique(onespno$year)))
    {
      meanvector[p]=jan1[p]+unique(onespno$peakday)[p] 
    }
    
    savetonextyr=0
    pred.quantity= pred.trapcount=numeric()
    for(l in 1: length(unique(onespno$cycle)))
    {
      include=which(onesp$cycle==l)
      dayinclude=datamonth[l]:(datamonth[l]+11)
      pred.quantity[dayinclude]=unique(onespno$peak[onesp$cycle==l])*pnorm(q=day1[dayinclude],mean=meanvector[l],sd=unique(onespno$sd))
      pred.trapcount[dayinclude]=c(0,diff(pred.quantity[dayinclude]))
      pred.trapcount[dayinclude[1]]=savetonextyr + pred.trapcount[dayinclude[1]]
      savetonextyr=unique(onespno$peak[onespno$cycle==l])-(pred.quantity[max(dayinclude)])
      
    }
    pred.trapcount=pred.trapcount[which(is.na(pred.trapcount)==FALSE)]
    results=data.frame(species,create.fulldate(fromjulian(day1+Dj, dateform="%Y-%m-%d")),model=pred.trapcount)
    if (i ==1) allresults = results else allresults = rbind (allresults, results)
  }
  return(allresults)
}

####3- Third, calculate the accumulated amounts for each year of the dataset####

#dd=yearvalues(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", beginyearfile=beginyearfile,spstart=1,spend=45)

## write.table(dd, file="monthly values of seed model.txt", sep="\t", row.names=F)
yearvalues=function(file="Nouragues model all spp.txt", fileno="NRG model all spp without Dj.txt",beginyearfile=beginyearfile,spstart=1,spend=45){
  
  sprange=spstart:spend
  nrgdata=read.delim(file)
  nrgdatano=read.delim(fileno)
  spnames=sort(unique(nrgdata$sp))
  beginyr=read.delim(file=beginyearfile)
  
  for (i in 1:length(sprange))
    
  {  
    species=as.character(spnames[sprange][i])
    onesp=nrgdata[nrgdata$sp==species,]
    onespno=nrgdatano[nrgdatano$sp==species,]
    startyear=min(unique(onesp$year))
    endyear=max(unique(onesp$year))
    startyearno=min(unique(onespno$year))
    endyearno=max(unique(onespno$year))
    bgy=beginyr$beginyr[beginyr$sp==species]
    Dj=round(ifelse(bgy<182.5,bgy,bgy-366),0)
    
    oneyr=numeric()  
    #for(k in 1:length((startyear+1):(endyear-1)))#this is for calculating the julian of the first day of each month
    for(k in 1:length((startyear):(endyear)))
    {
      oneyr=rep(startyear+(k-1),12)  
      if (k ==1) yrvector=oneyr else yrvector=c(yrvector, oneyr)
    }
    
    firstmonth=as.numeric(unique(onesp$month)[1])
    lastmonth=onesp$month[onesp$julian==sort(unique(onesp$julian))[length(unique(onesp$julian))]]
    #day1=c(tojulian(pst(firstmonth:12,'/01/',startyear)), c(tojulian(pst(rep(1:12, length((startyear+1):(endyear-1))),"/01/",yrvector))),c(tojulian(pst(1:lastmonth,'/01/',endyear))))
    day1=c(tojulian(pst(rep(1:12, length(startyear:endyear)),"/15/",yrvector)))
    day1min=which(day1==tojulian(paste(firstmonth,"/15/",startyear),dateform = "%m /%d/ %Y"))
    day1max=which(day1==as.numeric(tojulian(paste(lastmonth,"/15/",endyear),dateform = "%m /%d/ %Y")))
    day1=day1[day1min:day1max]
    day1=day1-Dj
    day1=day1[which(day1>=onespno$julian[1])]
    #which(day1>onespno$julian[length(onesp$julian)])
    #datanew=data.frame(julian=day1)
    #datanew$year=create.fulldate(fromjulian(day1, dateform="%Y-%m-%d"),format="%Y-%m-%d")$year
    #mm=model.seedrain(data=datanew, peaks=unique(onesp$peak), peakdays=unique(onesp$peakday), SD=unique(onesp$sd), startyear=min(unique(onesp$year)),endyear=max(unique(onesp$year)))
    #data1=create.fulldate(fromjulian(day1, dateform="%Y-%m-%d"),format="%Y-%m-%d")$year
    datamonth=seq(1,length(unique(onespno$cycle))*12,12)
    
    jan1=tojulian(pst('01/01/',startyearno:endyearno))   
    meanvector=numeric()
    for (p in 1:length(unique(onespno$year)))
    {
      meanvector[p]=jan1[p]+unique(onespno$peakday)[p] 
    }
    
    savetonextyr=0
    pred.quantity= pred.trapcount=numeric()
    for(l in 1: length(unique(onespno$cycle)))
    {
      pred.quantity[l]=unique(onespno$peak[onesp$cycle==l])*pnorm(q=meanvector[l],mean=meanvector[l],sd=unique(onespno$sd))
      pred.trapcount[include]=c(0,diff(pred.quantity[include]))
      pred.trapcount[dayinclude[1]]=savetonextyr + pred.trapcount[dayinclude[1]]
      savetonextyr=unique(onespno$peak[onespno$cycle==l])-(pred.quantity[max(dayinclude)])
      
    }
    pred.trapcount=pred.trapcount[which(is.na(pred.trapcount)==FALSE)]
    results=data.frame(species,create.fulldate(fromjulian(day1+Dj, dateform="%Y-%m-%d")),model=pred.trapcount)
    if (i ==1) allresults = results else allresults = rbind (allresults, results)
  }
  return(allresults)
}

