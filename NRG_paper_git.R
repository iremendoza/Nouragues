
#load(".\\NRG.results2.Rdata")
#load(".\\EstimatedCV.Rdata")
#load(".\\tetpan2.Rdata")

#setwd("C:\\Projects\\Nouragues")
#setwd("C:\\Brunoy\\hierarchical models\\working files hierarchical models\\")

nrgresults = read.delim(file = "nouragues results parameters per year.txt")

#source("C:\\Irene\\Brunoy\\hierarchical models\\working files hierarchical models\\hierarchical-enso-local climate.r")
source(".\\hierarchical models_git.r")
source(".\\Nouragues_functions_git.r")
source(".\\hierarchicalModelsRepyear_git.r")
source(".\\dataset management_git.r")
#source("C:\\Irene\\Brunoy\\Base Datos Nouragues\\Metadata Joe Wright\\NRG\\Rcode_CV_PhenYr_20110509.r")
#source("C:\\Irene\\Brunoy\\Base Datos Nouragues\\Metadata Joe Wright\\modeling max likelihood\\BCI-max likelihood.r")


#attach('CTFSRPackage.rdata')

library(ggplot2)
library(dplyr)

####DATASETS####

nourage <- "nouragues.txt" ## raw data with seed counts and density of seeds/m2 per census and species
beginyearfile <- "beginyearseeds 2011 newfecha_NRG.txt"

#biomassraw = read.delim(file="biomass raw.txt")
#biomassraw$fecha = create.fulldate(biomassraw$date, format="%d/%m/%Y")
#summonth = aggregate(data.frame(fruits=biomassraw$fruits, flowers=biomassraw$flowers), by=list(month=biomassraw$month, year=biomassraw$year),sum)
#nbcens = aggregate(data.frame(nbcens=biomassraw$Census), by=list(month=biomassraw$month, year=biomassraw$year),lengthunique)
#summonth$nbcens=nbcens$nbcens
#summonth$fr=summonth$fruits/summonth$nbcens
#summonth$fl=summonth$flowers/summonth$nbcens
#biomass = read.delim(file="biomass all months.txt") ##this dataset includes the biomass values per month and standarized by the number of censuses
#ensoanomalies <- read.delim(file="ENSO anomalies - 2012.txt")
#nrgraw = read.delim("nrg 2011.txt",header=T)
#NRGallspp = read.delim(file="Nouragues model all spp.txt")
#NRGallsppnoDj = read.delim("NRG model all spp without Dj.txt") #this dataset includes all the parameters models without the correction for date 
#NRGmonthly = read.delim(file="monthly values of seed model.txt")
#nrghyper = read.delim(file="Nouragues results hyperparameters.txt")
#est<-read.delim("number total spp per month estimated.txt",header=T)
#totseed<-read.delim(file="total number of seeds per species.txt")

clim <- read.table("local climate data Nouragues.txt", header = T)
#nrgprior<-read.delim(file="Nouragues priors per year.txt")
clim$dates = strptime(paste("1-",clim$month,"-",clim$year), format="%d - %m - %Y")
clim$julian = tojulian(clim$dates,dateform = "%Y-%m-%d")
clim$yday = clim$dates$yday+1
clim$month = clim$dates$mon+1
clim$dates = strptime(paste("1-",clim$month,"-",clim$year), format="%d - %m - %Y")
clim$julian = tojulian(clim$dates,dateform = "%Y-%m-%d")
clim$yday = clim$dates$yday+1

newman = read.delim(file = "Nouragues manual estimated new.txt")


#census<-read.delim(file="census list.txt")

#climraw<-read.delim(file="Nouragues climate 2003-2012_manual.txt") ##this includes raw data of local climate at Nouragues (manual station)
#climraw$dates=strptime(paste(climraw$Day,"-",climraw$Month,"-",climraw$Year), format="%d - %m - %Y")
#climraw$julian=tojulian(climraw$dates,dateform = "%Y-%m-%d")
#is.na(climraw$TempMin)<-which(climraw$TempMin<19)

#autoraw<-read.delim(file="Nouragues_automatique climat journalieres 2006-2009.txt")
#autoraw$dates=strptime(paste(autoraw$Day,"-",autoraw$Month,"-",autoraw$Year), format="%d - %m - %Y")
#autoraw$julian=tojulian(autoraw$dates,dateform = "%Y-%m-%d")

#newmansummary=aggregate(data.frame(tmin=newman$tmine,tmax=newman$tmaxe,rain=newman$raine),by=list( month=newman$Month,year=newman$Year),mean, na.rm=T)


####FIGURES OF THE PAPER####

####  FIGURE 1 OF THE PAPER   ######################################
### beautiful graph of mean number of species  and climatogram ###

#figure1(datfile="new",graphname="figure1.tif",est=est,biomass=biomass)
figure1 = function(datfile = c("new", "old"), graphname="figure1.tif", est = est, biomass = biomass, newman = "Nouragues manual estimated new.txt")
  #newman is the estimated dataset of local climate from the climatic station of Nouragues
  #old datafile referes to "local climate data Nouragues.txt"
{
  #dat <- read.table("local climate data Nouragues.txt", header = T) #this was for the previous version of the paper
  if (datfile == "old") {
    dat=read.delim("local climate data Nouragues.txt")
    means <- aggregate(dat[, 3:5], by = list(month = dat$month), mean)
    sds  <- aggregate(dat[, 3:5], by = list(month = dat$month), sd)
  }
  
  if (datfile == "new") {
    dat = read.delim(newman)
    dry = aggregate(data.frame(rain=dat$raine), by = list(Month=dat$Month,Year=dat$Year), sum,na.rm=T)
    for (k in 1:dim(dry)[1]){   
      dry$lack[k]=length(which(is.na(dat$raine[dat$Year==dry$Year[k]&dat$Month==dry$Month[k]])==T))
    }
    # we establish a threshold of 4 days for keeping rainfall data of that month
    dry2=dry[-which(dry$lack>4),]
    monthrain=aggregate(data.frame(rain=dry2$rain),by=list(month=dry2$Month),mean)
    sdrain=aggregate(data.frame(rain=dry2$rain),by=list(month=dry2$Month),sd)
    
    means <- data.frame(aggregate(data.frame(tmin=dat[,16],tmax=dat[,17]), by = list(month = dat$Month), mean,na.rm=T),rain=monthrain$rain)
    sds  <- data.frame(aggregate(data.frame(tmin=dat[,16],tmax=dat[,17]), by = list(month = dat$Month), sd,na.rm=T),rain=sdrain$rain)
  }
  
  labels <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
  
  #est<-read.delim("number total spp per month estimated.txt",header=T)
  sppmeans=aggregate( est$estnumbspp, by = list(month = est$month), mean)
  names(sppmeans)=c("month","nbspp")
  sppsds=aggregate( est$estnumbspp, by = list(month = est$month), sd)
  names(sppsds)=c("month","nbspp")
  
  #biomass=read.delim(file="biomass all months.txt",as.is=T)
  biomass=biomass[-1,]
  meanbio=aggregate(data.frame(fruit=biomass$fruits,flower=biomass$flowers), by=list(month=biomass$month), mean)
  sdbio=aggregate(data.frame(fruit=biomass$fruits,flower=biomass$flowers), by=list(month=biomass$month), sd)
  
  #x11(height=10,width=6) 
  tiff(filename=graphname,width = 1900, height = 3000,pointsize=12, res=300)
  par(las = 1, bty = "o", tcl = 0.2, mar = c(2, 5, 1, 5), mgp = c(0.25, 0.25, 0),cex.axis=1.5,lwd=1.5)
  par(mfrow=c(3,1))
  
  plot(meanbio$fruit,type="o",axes=FALSE, xlab="",ylab="", col="black", ylim=c(0,15) )
  segments(1:12, meanbio$fruit- sdbio$fruit, 1:12, meanbio$fruit + sdbio$fruit, lwd = 2, lend = 3, col = "black")
  axis(1, at = 1:12, labels = NA)
  axis(2) 
  mtext(text=expression(paste("Biomass (g ",m^-2, ")", sep = " ")), 2, las = 3, line = 3) 
  text(2,15,labels="A",pos=2, offset=2, cex=2)
  par(new=TRUE)
  ji12=jitter(c(1:12))
  plot(ji12,flowermeans$flower, type="o",axes=FALSE, xlab="", ylab="", col="grey70", ylim=c(0,15) ,lwd=2)
  segments(ji12, flowermeans$flower- flowersds$flower, ji12, flowermeans$flower + flowersds$flower, lwd = 2, lend = 3, col = "grey70")
  legend(7,15,legend=c("fruits","flowers"), col=c("black","grey70"),lty=1,lwd=3,bty="n",cex=2)
  
  print(cor.test(means$rain, meanbio$flower))
  print(cor.test(means$rain, meanbio$fruit))
  
  plot(sppmeans$nbspp,type="o",axes=FALSE, ylim=c(10,40), xlab="",ylab="",col="black")
  text(2,40,labels="B",pos=2, offset=2, cex=2)
  axis(1, at = 1:12, labels = NA)
  axis(2)
  segments(1:12, sppmeans$nbspp - sppsds$nbspp, 1:12, sppmeans$nbspp + sppsds$nbspp, lwd = 2, lend = 3, col = "black")
  mtext("Number of fruiting species", 2, las = 3, line = 3)
  
  print(cor.test(means$rain,sppmeans$nbspp))
  
  #plot(1:12, means$rain, axes = F, type = "n", xlab = "", ylab = "", ylim = range(means$rain - sds$rain, means$rain+ sds$rain))
  plot(1:12, means$rain, axes = F, type = "n", xlab = "", ylab = "", ylim = c(0,800))
  text(2,800,labels="C",pos=2, offset=2, cex=2)
  polygon(c(1:12, 12:1), c(means$tmin, rev(means$tmax))*20, col = "#FF000050", border = NA)
  segments( 1:12, 0, 1:12, means$rain, lwd = 40, lend = 3, col = "blue",lend=3)
  segments(1:12, means$rain - sds$rain, 1:12, means$rain+ sds$rain, lwd = 2, lend = 3, col = "black")
  
  
  axis(1, at = 1:12, labels = labels)
  axis(2, at = 100*1:5, col.axis = "blue")
  axis(4, at = 200*1:4, labels = 10*1:4, col.axis = "red")
  mtext("Temperature (ºC)", 4, las = 3, line = 3, col = "red")
  mtext("Precipitation (mm)", 2, las = 3, line = 3, col = "blue")
  
  dev.off() 
}


####FIGURE 2 OF THE PAPER ###############

#figure2 plots the hyperparameters for each species in Nouragues (in an exponential scale)
#figure2(trfile="Nouragues results hyperparameters.txt",longnames="total number of seeds per species.txt",filename="figure2.tif")
figure2 = function(trfile = "Nouragues results hyperparameters.txt", longnames = "total number of seeds per species.txt", filename = "figure2.tif")
  
{
  
  #x11(height=9,width=9)
  tiff(filename=filename,height=1200,width=1900,res=300)
  par(mar=c(3,3,3,1),oma=c(1,1,1,1))
  tr = read.delim(file = trfile)
  tr$mu2=tr$mu
  tr$mu2[which(tr$mu<=0)]= tr$mu[which(tr$mu<=0)]+365.25
  tr$mu2[which(tr$mu>365)]= tr$mu[which(tr$mu>365)]-365.25
  negdate=which(tr$mu<=0)
  
  
  totseed = read.delim(file = longnames)
  names(totseed)=c("species","longname","totseed","form","disp","fruit", "length", "width", "Smythe")
  spnames2 = totseed$longname
  dat <- merge(tr,totseed,by="species")
  zoo = which(totseed$disp=="zoo")
  bal = which(totseed$disp=="bal")
  ane = which(totseed$disp=="ane")
  
  colores=1:45
  colores[zoo]="red"
  colores[ane]="blue"
  colores[bal]="blue"
  
  
  maxyr=which.max(tr$logmu)
  
  oneyrpred=14*exp(tr$logmu[maxyr])*dnorm(1:365,mean=tr$mu2[maxyr],sd=tr$bestSD[maxyr])
  plot(1:365,oneyrpred,type='l',ylab='',xlab='',ylim=c(-0.5,max(oneyrpred+1)),lwd=1.25, cex=1, axes=F)
  axis(side=2, las=2)
  axis(side=1, at=c(1,32,60,92,122,152,183,214,245,275,306,336,366), labels=c("J", "F", "M", "A","M", "J", "J", "A", "S","O","N","D","J"))
  #polygon(c(248,341,341,248),c(0,0,3,3),col = "#FF000020", border = NA)
  #polygon(c(0,248,248,0),c(0,0,3,3),col = "#0000FF20", border = NA)
  polygon(c(212,334,334,212),c(0,0,30,30),col = "grey80", border = NA)
  #polygon(c(341,400,400,341),c(0,0,3,3),col = "#0000FF20", border = NA)
  mtext(3, text="Nouragues (2001-2011)",line=2,cex=1.5)
  mtext(2,text="biweekly seedfall",line=3,cex=1)
  mtext(1,text="month of the year",line=3,cex=1)
  for (i in 1:dim(tr)[1])
  {
    oneyrpred=14*exp(tr$logmu[i])*dnorm(1:365,mean=tr$mu2[i],sd=tr$bestSD[i])
    lines(1:365,oneyrpred,col=colores[i], lwd=2) 
  }
  
  lines(1:365, 14*mean(exp(tr$logmu))*dnorm(1:365,mean=mean(tr$mu2),sd=mean(tr$bestSD)),col="black",lwd=4)
  legend(10,37,lty=c(1,1),col=c("black","red","blue"),legend=c("community","biotic","abiotic"),bty="n",cex=0.8,horiz=T,lwd=2)
  
  #### contingency test associated with figure 2: we test the Ho that the peaks of fruit production were uniform through the year####
  # March-May is 25% of the year. Does 21 species differ significantly from 25% of 45 species?
  
  ## Species having peaks during the rainy peak (Mar-May)
  dim(tr[which(tr$mu2 >= 60 & tr$mu2 <= 151),])[1]
  chisq.test(x=c(22,23),p=c(0.25,0.75))
  
  # contingency test associated with figure 2: we test the Ho that the peaks of fruit production were uniform through the year for all dispersal modes
  biotic = which(totseed$disp=="zoo")
  abiotic = which(totseed$disp=="bal"| totseed$disp=="ane")
  
  
  ## Species having peaks during the rainy peak (Mar-May) 
  rainypeak <- dat[which(dat$mu2 >= 60 & dat$mu2 <= 151),]
  dim(rainypeak)
  length(which(rainypeak$disp == "zoo"))
  rainy2<-dat[which(dat$mu2>=152&dat$mu2<=243),]
  dim(rainy2)
  length(which(rainy2$disp=="zoo"))
  dry <- dat[which(dat$mu2 >= 212 & dat$mu2 <= 334),]
  rest <- dat[which(dat$mu2 < 212 | dat$mu2 > 334),]
  dim(dry)
  length(which(dry$disp == "zoo"))
  length(which(dry$disp == "ane"))
  
  rainy1<-dat[which(dat$mu2>=335|dat$mu2<=59),]
  dim(rainy1)
  length(which(rainy1$disp=="zoo"))
  
  drydisp = matrix(ncol=2,nrow=2)
  colnames(drydisp)=c("dryseason","rest")
  rownames(drydisp)=c("biotic","abiotic")
  drydisp[,1]=c(length(which(dry$disp == "zoo")),length(which(dry$disp == "ane")))
  drydisp[,2]=c(length(which(rest$disp == "zoo")),length(which(rest$disp == "ane"|rest$disp == "bal")))
  
  chisq.test(disptest) 
  
  
  disptest=matrix(ncol=2,nrow=2)
  colnames(disptest)=c("rainypeak","rest")
  rownames(disptest)=c("biotic","abiotic")
  disptest[1,]=c(11,13)
  disptest[2,]=c(11,10)
  biotic<-c(11,1,5,7)
  abiotic<-c(11,3,3,4)
  chisq.test(disptest) ## no signficant differences between dispersal mode of species
  
  dev.off()
  
}

#figure2nodispersal plots the hyperparameters for each species in Nouragues (in an exponential scale) without dispersal modes (for ATBC presentation)
#figure2nodispersal(trfile="Nouragues results hyperparameters.txt",longnames="total number of seeds per species.txt",filename="NRG for ATBC2016.tif")

figure2nodispersal=function(trfile="Nouragues results hyperparameters.txt",longnames="total number of seeds per species.txt",filename="figure2.tif")
  
{
  
  #x11(height=9,width=9)
  tiff(filename=filename,height=1200,width=1900,res=300)
  par(mar=c(3,3,3,1),oma=c(1,1,1,1))
  tr=read.delim(file=trfile)
  tr$mu2=tr$mu
  tr$mu2[which(tr$mu<=0)]= tr$mu[which(tr$mu<=0)]+365.25
  tr$mu2[which(tr$mu>365)]= tr$mu[which(tr$mu>365)]-365.25
  negdate=which(tr$mu<=0)
  
  
  totseed=read.delim(file=longnames)
  names(totseed)=c("species","longname","totseed","form","disp","fruit")
  spnames2=totseed$longname
  dat<-merge(tr,totseed,by="species")
  zoo=which(totseed$disp=="zoo")
  bal=which(totseed$disp=="bal")
  ane=which(totseed$disp=="ane")
  
  colores=1:45
  colores[zoo]="blue"
  colores[ane]="blue"
  colores[bal]="blue"
  
  
  maxyr=which.max(tr$logmu)
  
  oneyrpred=14*exp(tr$logmu[maxyr])*dnorm(1:365,mean=tr$mu2[maxyr],sd=tr$bestSD[maxyr])
  plot(1:365,oneyrpred,type='l',ylab='',xlab='',ylim=c(-0.5,max(oneyrpred+1)),lwd=1.25, cex=1, axes=F)
  axis(side=2, las=2)
  axis(side=1, at=c(1,32,60,92,122,152,183,214,245,275,306,336,366), labels=c("J", "F", "M", "A","M", "J", "J", "A", "S","O","N","D","J"))
  #polygon(c(248,341,341,248),c(0,0,3,3),col = "#FF000020", border = NA)
  #polygon(c(0,248,248,0),c(0,0,3,3),col = "#0000FF20", border = NA)
  polygon(c(212,334,334,212),c(0,0,30,30),col = "#CCCCCC42", border = NA)
  #polygon(c(341,400,400,341),c(0,0,3,3),col = "#0000FF20", border = NA)
  mtext(3, text="Nouragues (2001-2011)",line=2,cex=1.5)
  mtext(2,text="biweekly seedfall",line=3,cex=1)
  mtext(1,text="month of the year",line=3,cex=1)
  for (i in 1:dim(tr)[1])
  {
    oneyrpred=14*exp(tr$logmu[i])*dnorm(1:365,mean=tr$mu2[i],sd=tr$bestSD[i])
    lines(1:365,oneyrpred,col=colores[i], lwd=2) 
  }
  
  lines(1:365, 14*mean(exp(tr$logmu))*dnorm(1:365,mean=mean(tr$mu2),sd=mean(tr$bestSD)),col="black",lwd=4)
  dev.off()
  
}


####FIGURE 3 OF THE PAPER ###############
##how many species had their peaks during the wet season?

figure3 = function(file = "Nouragues results hyperparameters.txt", longnames = "total number of seeds per species.txt", filename = "figure3.tif"){
  
  tr = read.delim(file)
  tr$mu2 = tr$mu
  tr$mu2[which(tr$mu<=0)] = tr$mu[which(tr$mu<=0)]+365.25
  tr$mu2[which(tr$mu>365)] = tr$mu[which(tr$mu>365)]-365.25
  
  tr$CImu2.2 = tr$CImu2
  #tr$CImu2.2[which(tr$CImu2 <= 0)] = tr$CImu2[which(tr$CImu2 <= 0)]+365.25
  tr$CImu2.2[which(tr$CImu2 > 365)] = tr$CImu2[which(tr$CImu2 > 365)]-365.25
 
  tr$CImu97.2 = tr$CImu97
  #tr$CImu97.2[which(tr$CImu97 <= 0)] = tr$CImu97[which(tr$CImu97 <= 0)]+365.25
  tr$CImu97.2[which(tr$CImu97 > 365)] = tr$CImu97[which(tr$CImu97 > 365)]-365.25
  
   
  totseed = read.delim(file = longnames)
  names(totseed) = c("species", "longname", "totseed", "form", "disp", "fruit", "length", "width", "Smythe")
  long = merge(tr, totseed, by ="species")
  type1 = subset(long, long$Smythe == "Type 1")
  type2 = subset(long, long$Smythe == "Type 2")
  long[long$longname == "Serjania paucidentata", 2:4] = long[long$longname == "Serjania paucidentata", 2:4] - 365.25
  long[long$longname == "Tetragastris panamensis", 2:4] = long[long$longname == "Tetragastris panamensis", 2:4] + 365.25
  long[long$longname == "Licania membranacea", 2:4] = long[long$longname == "Licania membranacea", 2:4] + 365.25
  
  dat <- data.frame(sp = factor(long$longname[order(long$mu)], levels = long$longname[order(long$mu)]), mu = long$mu[order(long$mu)], CImu2 = long$CImu2[order(long$mu)], CImu97 = long$CImu97[order(long$mu)], 
                    SD = long$SD[order(long$mu)], Smythe = long$Smythe[order(long$mu)], disp = long$disp[order(long$mu)])
  
  dat$CV <- dat$SD/dat$mu
  zoo = which(dat$disp == "zoo")
  bal = which(dat$disp == "bal")
  ane = which(dat$disp == "ane")
  dat$disp2 <- character(length = dim(dat)[1])
  dat$disp2[zoo] = "biotic"
  dat$disp2[ane] = "abiotic"
  dat$disp2[bal] = "abiotic"
  
  colores = character(length = 45)
  colores[zoo] = "gray48"
  colores[ane] = "white"
  colores[bal] = "white"
  
  tiff(filename = filename, height = 2000, width = 3000, res = 300)
  par(las = 1, mar = c(3,13,1,1), oma = c(2,4,1,1), cex = 0.6)
  
  plot(dat$mu, 1:45, ylim = c(0, 50), xlim = c(-80, 450), pch = 21, bg = colores, cex = 1.5, ylab = "", xlab = "", axes = F)
  axis(side=1, at=c(1,32,60,92,122,152,183,214,245,275,306,336,366), labels = c("J", "F", "M", "A","M", "J", "J", "A", "S","O","N","D","J"), cex.axis = 1.5)
  #segments(x0 = dat$mu, x1 = dat$CImu2, y0 = 1:45, y1 = 1:45)
  #segments(x0 = dat$mu, x1 = dat$CImu97, y0 = 1:45, y1 = 1:45)
  segments(x0 = dat$mu-dat$SD, x1 = dat$mu+dat$SD, y0 = 1:45, y1 = 1:45)
  #segments(x0 = dat$mu, x1 = dat$CImu97, y0 = 1:45, y1 = 1:45)
  
  
  polygon(c(212, 334, 334, 212),c(0, 0,50,50),col = "#CCCCCC42", border = NA)
  axis(side = 2, at = 1:45, labels = dat$sp[seq(1,45,1)],las = 2, font = 3, cex.axis = 1.2)
  legend(10, 44, lty = c(1,1), pch = 21, legend = c("biotic","abiotic"), pt.bg = c("gray48", "white"), bty = "n", cex = 1.5, horiz = T, lwd = 2)
  mtext(side = 1, line = 3, text = "date of peaks (month of the year)")
  
  # The errorbars overlapped, so use position_dodge to move them horizontally
 # pd <- position_dodge(0.1) # move them .05 to the left and right
  #ggplot(dat, aes (x = mu, y = sp)) +
   # geom_point(position = pd) +
  #  geom_errorbar(aes(ymin = CImu2, ymax = CImu97), stat = "identity", width=.1, position=pd) 
  #+
    #geom_errorbar(aes(ymin = CImu2, ymax = CImu2), width=.1, position = )
  dev.off()
 
  ####classification of species according to hyperSD####
  length(which(dat$SD <= 45))
  summary(dat$SD)
  low <- dat[which(dat$SD <= 24.250),]
  medium <- dat[which(dat$SD > 24.250 & dat$SD <= 63.580),]
  high <- dat[which(dat$SD > 63.580),]
  
  variability <- matrix(ncol = 3, nrow = 2)
  colnames(variability) = c("low","medium", "high")
  rownames(variability) = c("biotic","abiotic")
  variability[,1] = c(length(which(low$disp == "zoo")), length(which(low$disp == "ane")))
  variability[,2] = c(length(which(medium$disp == "zoo")), length(which(medium$disp == "ane"|medium$disp == "bal")))
  variability[,3] = c(length(which(high$disp == "zoo")), length(which(high$disp == "ane")))
  chisq.test(variability)
  
  cvariability <- cbind(variability[,"low"] + variability[,"medium"] , variability[,"high"])
  chisq.test(cvariability) ## not significant differences in variability according to a chi-squared test
  
  glm1 <- glm(dat$SD ~ dat$disp2)
  summary(glm1)
  #boxplot(dat$SD ~ dat$disp2)
  }  

figure3.old = function(file = "Nouragues results hyperparameters.txt", filename = "figure3_old.tif"){
  
  tr = read.delim(file)
  str(tr)
  tr$mu2 = tr$mu
  tr$mu2[which(tr$mu<=0)] = tr$mu[which(tr$mu<=0)]+365.25
  tr$mu2[which(tr$mu>365)] = tr$mu[which(tr$mu>365)]-365.25
  
  totseed = read.delim(file = longnames)
  names(totseed) = c("species", "longname", "totseed", "form", "disp", "fruit", "length", "width", "Smythe")
  long = merge(tr, totseed, by ="species")
  type1 = subset(long, long$Smythe == "Type 1")
  type2 = subset(long, long$Smythe == "Type 2")
  
  sort(long$mu2)
  
  tiff(filename = filename, height = 1000, width = 2000, res = 300)
  par(mfrow = c(1,2), las = 1, mar = c(3,3,1,0), oma = c(2,2,1,1), cex = 0.6)
  #tt=hist(tr$mu2, breaks=c(0,30,60,90,120,150,180,210,240,270,300,330,360,390), xlab="",ylab="",right=FALSE,  axes=F, main="", ylim=c(0,9))
  tt = hist(tr$mu2, breaks=c(0,60,120,180,240,300,360,420), xlab="",ylab="",right=FALSE,  axes=F, main="",ylim=c(0,15),cex=0.5)
  mtext(1,text="peakday (hypermean)",line=2.5,cex=1)
  mtext(2,text="number of species",line=3,cex=1.5,las=0)
  axis(side=1, at=c(0,60,120,180,240,300,360,420), labels=c("0","60","120","180","240","300","360","365"))
  axis(side=2, at=c(seq(0,15,5)), labels=T,las=1)
  data.frame(breaks=tt$breaks[1:length(tt$breaks)-1], counts=tt$counts)
  tsd=hist(tr$SD,breaks=c(0,15,30,45,60,75,90,105,120),xlab="",ylab="",main="", axes=F, right=F)
  mtext(1,text = "variation of peakday (number of days)",line = 2.5, cex = 1)
   axis(side=2, las = 2)
  axis(side=1, at=c(0,15,30,45,60,75,90,105,120), labels=T, las = 1)
  #c("0","","","45","","","90","","","135")
  
  data.frame(breaks=tsd$breaks[1:length(tsd$breaks)-1], counts=tsd$counts)
  #hist(exp(t$logmu)+1,xlab="peak of seed production (number seeds/year)",ylab="number of species",main=site)
  #hist(exp(t$logSD)+1,xlab="SD of peak of seed production (number seeds/year)",ylab="number of species",main=site)
  # hist(tr$logmu,main="", ylab="number of species", xlab="hyperpeak")
  
  dev.off()
  
}  


####FIGURE 4 OF THE PAPER ###############
### coefficient of variation of seed production for the 45 species

figure4 = function(file = "nouragues results parameters per year.txt", hyper = "Nouragues results hyperparameters.txt", EstimatedCV = EstimatedCV, longnames = "seed size.txt", filename = "figure4.tif") {
  
  nrg = read.delim(file)
  #spnames=sort(unique(nrg$sp))
  totseed = read.delim(file = longnames)
  names(totseed) = c("Taxa", "species", "longname", "family", "Poncy", "adult", "observations", "form", "disp", "fruit", "totseed", "measured", "length", "width", "D3", "comments.measure")
  spnames2 = totseed$longname
  spmeans = aggregate(data.frame(peak = nrg$peak),by=list(species = nrg$species),mean)
  spsd = aggregate(data.frame(peak = nrg$peak), by = list(species = nrg$species),sd)
  condensed = merge(spmeans, spsd, by = "species")
  names(condensed)=c("species","mean","sd")
  condensed2 = merge(totseed, condensed, by="species")
  condensed2$CV = condensed$sd/condensed$mean
  orco = order(condensed2$CV)
  
  ##this is for the second way of calculating the CV from a lognormal
  tr <-read.delim(hyper)
  tr$CV <- sqrt(exp(tr$logSD^2)-1 )
  orCV <- order(log(tr$CV+1))
  condensedtr = merge(tr, totseed, by = "species")
  
  #this is the third and correct way, using the table calculated by Rick Condit
  orest = order(EstimatedCV$meancv)
  EstimatedCV$species <- rownames(EstimatedCV)
  condensedest = merge(EstimatedCV, totseed, by = "species")
  
  #split.screen(c(1,2))
  tiff(filename = filename, height = 1600, width = 2500, pointsize = 24) #
  par(mar = c(20,5,12,1), cex = 1)
  zoo = which(condensedest$disp[orest] == "zoo")
  bal = which(condensedest$disp[orest] == "bal")
  ane = which(condensedest$disp[orest] == "ane")
  colores=1:45
  colores[zoo] = "gray68"
  colores[ane] = "white"
  colores[bal] = "white"
  #barplot(condensed2$CV[orco], main="", space=0,las=2,axes=F, ylim=c(0,3),col=colores) #this is for the old way of calculating the CV of peak
  #barplot(log(condensedtr$CV + 1)[orCV], main="", space = 0, las = 2, axes = F, ylim = c(0,30), col = colores)
  #ggplot(condensedest, aes(x = longname[orest], y = as.numeric(meancv[orest]))) +
           #geom_bar(stat = "identity")
  barplot(as.numeric(condensedest$meancv[orest]), main="", space=0,las=2, ylim=c(0,40),col=colores, axes = F)
  segments(y0 = as.numeric(condensedest$meancv[orest]), y1 = as.numeric(EstimatedCV$upper[orest]), x0 = seq(0.5, 44.5, 1), x1 = seq(0.5, 44.5, 1), lwd = 2, col = "black")
  segments(y0 = as.numeric(condensedest$lower[orest]), y1 = as.numeric(condensedest$meancv[orest]), x0 = seq(0.5, 44.5, 1), x1 = seq(0.5, 44.5, 1), lwd = 2, col = "black")
  
  arrows(y0 = as.numeric(condensedest$meancv[orest]), y1 = as.numeric(EstimatedCV$upper[orest]), x0 = seq(0.5, 44.5, 1), x1 = seq(0.5, 44.5, 1), lwd = 2, col = "black", angle = 90,
         code = 3, length = 0.05)
  arrows(y0 = as.numeric(condensedest$lower[orest]), y1 = as.numeric(condensedest$meancv[orest]), x0 = seq(0.5, 44.5, 1), x1 = seq(0.5, 44.5, 1), lwd = 2, col = "black", angle = 90,
         code = 3, length = 0.05)
  
  
  axis(side=2,las=2, cex.axis=1.75)
  mtext(side=2, "Coefficient of variation of seed production", line=3,cex=2, las=3)
  #axis(side=1, at=seq(0.5,45,1), labels = condensed2$longname[orco][seq(1,45,1)],las=2,font=3, cex.axis=1.5) #this is for the old way of calculating the CV of peak
  axis(side=1, at=seq(0.5,45,1), labels = condensedest$longname[orest][seq(1,45,1)],las=2,font=3, cex.axis=1.5) #this is for the new way of calculating the CV of peak
  
  abline(h=1,lwd=2)
  text(0,28, "A", cex=2)
  legend(4,28,fill=c("gray68","white"),legend=c("biotic","abiotic"),bty="n",cex=2)
  par(new=T, mar = c(2,2,2,2))
  m=matrix(ncol=3,nrow=4)
  m[1,3]=2
  m[2:4,]=1
  m[1,c(1,2)]=1
  layout(m, widths=c(0.9,0.9))
  #CVspp(file=file)
  h <- hist(condensedest$meancv,  breaks = c(0,1,5,10,15,20,25), plot = F)
  barplot(h$counts, space = 0, axes =F)
  axis(side = 1, las = 1, cex.axis = 1.75, at = 1:6, labels = c( "1", "5", "10", "15", "20", "25"))
  axis(side = 2, las = 1, cex.axis = 1.75)
  abline(v=1,lwd=3)
  text(3,15, "B", cex=2)
  mtext(side=1,"coefficient of variation of seed production",line=3,cex=1.5)
  mtext(side=2, "number of species",line=3,cex=1.5)
  dev.off()
}


###FIGURE 5 OF THE PAPER ###############

figure5=function(file=nourage, fit=results, beginyearfile=beginyearfile)
{
  nrgdata=read.delim(file)
  beginyr=read.delim(file=beginyearfile)
  fulldata=merge(nrgdata,beginyr,by="sp", all.x=TRUE)
  spnames=sort(unique(fulldata$sp))
  jpeg(filename="Figure5.jpg",width = 900, height = 1050,quality=100)
  par(mfrow=c(6,2),mar=c(2.5,1,2,5), oma=c(4,5,1,1),las=1, cex=1)
  sprange=c(45,10,23,8,18,20)
  for (i in 1:length(sprange)){
    species=as.character(spnames[sprange][i])
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
    
    #plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,max(fittrans$bestpeak)))
    
    if (i ==1|i==5){
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
      mtext(side=3, text=species, font=3,line=0.5, cex=1.5)
      
      maxyr=which.max(fittrans$bestpeak)
      oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(-185:185,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
      plot(-185:185,oneyrpred,type='l',ylab='biweekly seedfall',xlab='day of year',axes=FALSE,lwd=2, bty="l")
      axis(side=1, at=c(-180,-120,-60,1,60,120,180),labels=c("185","245","305","0","60","120","180") )
      axis(side=2)
      for(k in 1:length(fittrans$bestpeak)) 
      { 
        oneyrpred=14*fittrans$bestpeak[k]*dnorm(-185:185,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
        lines(-185:185,oneyrpred,col=k,lwd=2) 
      }
      #mtext(side=3, text=species, font=3,line=0.75, cex=1.5)
    }
    if(i ==6){
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
      mtext(side=3, text=species, font=3,line=0.5, cex=1.5)
      mtext(side=1, text="years", font=1,line=3, cex=2) 
      
      maxyr=which.max(fittrans$bestpeak)
      oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(1:400,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
      plot(1:400,oneyrpred,type='l',ylab='',xlab='',lwd=2,ylim=c(0,max(oneyrpred+1)),bty="l")
      for(k in 1:length(fittrans$bestpeak)) 
      { 
        oneyrpred=14*fittrans$bestpeak[k]*dnorm(1:400,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
        lines(1:400,oneyrpred,col=k,lwd=2) 
      } 
      mtext(side=1, text="day of year", font=1,line=3, cex=2)   
    }  
    
    if(i == 2|i ==3){
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
      mtext(side=3, text=species, font=3,line=0.5, cex=1.5)
      maxyr=which.max(fittrans$bestpeak)
      oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(1:400,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
      plot(1:400,oneyrpred,type='l',ylab='',xlab='',lwd=2,ylim=c(0,max(oneyrpred+1)),bty="l")
      for(k in 1:length(fittrans$bestpeak)) 
      { 
        oneyrpred=14*fittrans$bestpeak[k]*dnorm(1:400,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
        lines(1:400,oneyrpred,col=k,lwd=2) 
      } 
      #mtext(side=3, text=species, font=3,line=0.75, cex=1.5)   
    }  
    if(i == 4){
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
      mtext(side=2, text="number of seeds", font=1,line=3,adj=0, las=0,cex=2)   
      mtext(side=3, text=species, font=3,line=0.5, cex=1.5)
      
      maxyr=which.max(fittrans$bestpeak)
      oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(1:400,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
      plot(1:400,oneyrpred,type='l',ylab='',xlab='',lwd=2,ylim=c(0,max(oneyrpred+1)),bty="l")
      mtext(side=2, text="biweekly seedfall", font=1,line=3,adj=0, las=0,cex=2)
      for(k in 1:length(fittrans$bestpeak)) 
      { 
        oneyrpred=14*fittrans$bestpeak[k]*dnorm(1:400,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
        lines(1:400,oneyrpred,col=k,lwd=2) 
      } 
      
    }  
    
  }  
  
  dev.off()  
}


####FIGURE 5b OF THE PAPER ###############
## figure 5b is a short version of figure5 for the FAPESP relatorio cientifico

figure5b=function(file=nourage, fit=results, beginyearfile=beginyearfile)
{
  nrgdata=read.delim(file)
  beginyr=read.delim(file=beginyearfile)
  fulldata=merge(nrgdata,beginyr,by="sp", all.x=TRUE)
  spnames=sort(unique(fulldata$sp))
  jpeg(filename="Figure5b.jpg",width = 900, height = 500,quality=100)
  par(mfrow=c(2,2),mar=c(2.5,4,4,5), oma=c(4,5,1,1),las=1, cex=1)
  sprange=c(45,20)
  for (i in 1:length(sprange)){
    species=as.character(spnames[sprange][i])
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
    
    #plot(spdata$julian[inc],spdata$quantity[inc],pch=16,ylim=c(0,max(fittrans$bestpeak)))
    
    if (i ==1){
      plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="",ylab="")
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
      mtext(side=3, text=species, font=3,line=0.5, cex=1.5)
      mtext(side=2, text="number of seeds", font=1,line=3,adj=0, las=0,cex=1.5)
      
      
      maxyr=which.max(fittrans$bestpeak)
      oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(-185:185,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
      plot(-185:185,oneyrpred,type='l',ylab='',xlab='',axes=FALSE,lwd=2, bty="l")
      axis(side=1, at=c(-180,-120,-60,1,60,120,180),labels=c("185","245","305","0","60","120","180") )
      axis(side=2)
      for(k in 1:length(fittrans$bestpeak)) 
      { 
        oneyrpred=14*fittrans$bestpeak[k]*dnorm(-185:185,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
        lines(-185:185,oneyrpred,col=k,lwd=2) 
      }
      mtext(side=2, text="biweekly seedfall", font=1,line=3,adj=0, las=0,cex=1.5)
    }
    if(i ==2){
      plot(spdata$julian,spdata$quantity,pch=16,ylim=c(0,max(fittrans$bestpeak)),axes=F,xlab="",ylab="")
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
      mtext(side=3, text=species, font=3,line=0.5, cex=1.5)
      mtext(side=1, text="years", font=1,line=3, cex=2) 
      mtext(side=2, text="number of seeds", font=1,line=3,adj=0, las=0,cex=1.5)
      
      maxyr=which.max(fittrans$bestpeak)
      oneyrpred=14*fittrans$bestpeak[maxyr]*dnorm(1:400,mean=fittrans$bestpeakday[maxyr],sd=fittrans$bestSD)
      plot(1:400,oneyrpred,type='l',ylab='',xlab='',lwd=2,ylim=c(0,max(oneyrpred+1)),bty="l")
      for(k in 1:length(fittrans$bestpeak)) 
      { 
        oneyrpred=14*fittrans$bestpeak[k]*dnorm(1:400,mean=fittrans$bestpeakday[k],sd=fittrans$bestSD)
        lines(1:400,oneyrpred,col=k,lwd=2) 
      } 
      mtext(side=1, text="day of year", font=1,line=3, cex=2)  
      mtext(side=2, text="biweekly seedfall", font=1,line=3,adj=0, las=0,cex=1.5)
    }  
    
    
  }  
  
  dev.off()  
}


#### FIGURE 6 ################################
##figure6 plots the contribution of each species to the total amount of seed production per year

#figure6(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", beginyearfile=beginyearfile, spstart=1, spend=45, graph=5,longnames="total number of seeds per species.txt",graphname="figure6.tif")
figure6=function(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", beginyearfile=beginyearfile, spstart=1, spend=45,graph=5,longnames="total number of seeds per species.txt",graphname="figure6.tif")
  
{
  data=monthlyvalues(file=file,fileno=fileno, beginyearfile=beginyearfile,spstart=spstart,spend=spend)
  totseed=read.delim(file=longnames)
  names(totseed)=c("species","longname","totseed","form","disp","fruit")
  dispersal=data.frame(species=totseed$longname,disp=totseed$disp)
  datayr<-merge(data,totseed, by="species",all.x=T)
  sumyr=aggregate(data.frame(sumseed=datayr$model), by=list(species=datayr$longname, year=datayr$year), sum)
  years=split(sumyr$sumseed, sumyr$year)
  tiff(filename=graphname,width = 4000, height = 3100,res=300)
  par(mfrow=c(5,2), mar=c(2,2,1,0),oma=c(3,5,3,3))  
  for (i in c(1:10))
  {
    oneyr=sumyr[sumyr$year==unique(sumyr$year)[i],]
    oneyr$perc=(oneyr$sumseed/sum(oneyr$sumseed))*100
    orspp=order(oneyr$sumseed, decreasing=T)
    spplist=as.vector(oneyr$species[orspp][1:15])
    perc=as.vector(oneyr$perc[orspp][1:15])
    spshort=strsplit(as.character(oneyr$species[orspp][1:graph]),' ')
    genuscode=spcode=dispersalmode=character()
    for(k in 1:length(spshort))
    {
      genuscode[k]=left(spshort[[k]][1],1)
      #spcode[k]=left(spshort[[k]][2],3)
      spcode[k]=spshort[[k]][2]
      dispersalmode[k]=as.character(dispersal$disp[which(dispersal$species==spplist[k])])
      
    }
    spnames=pst(genuscode,"."," ", spcode)
    splist=data.frame(spnames,dispersalmode)
    zoo=which(splist$dispersalmode=="zoo")
    bal=which(splist$dispersalmode=="bal")
    ane=which(splist$dispersalmode=="ane")
    colores=1:graph
    colores[zoo]="gray28"
    colores[ane]="white"
    colores[bal]="white"
    if (i==1|i==9) barplot(oneyr$perc[orspp][1:graph],width=0.7,space=0.2,axes=F, las=1, col=colores,ylim=c(0,60)) else barplot(oneyr$perc[orspp][1:graph],width=0.7,space=0.2,axes=F, las=1, col=colores,ylim=c(0,40))
    axis(side=1,at=c(0.5,1.35,2.2,3.05,3.90),labels=spnames, cex.axis=1.25,font=3)
    if (i==1|i==9) axis(side=2,las=1,at=c(0,20,40,60),cex.axis=1.5) else axis(side=2,las=1,at=c(0,10,20,30,40),cex.axis=1.5)
    if (i ==2) legend(3,43,fill=c("gray28","white"),legend=c("biotic","abiotic"),bty="n",cex=1.8,horiz=F)
    mtext(side=3, text=unique(sumyr$year)[i], line=-2.5,cex=1.5)
    results=data.frame( year=unique(sumyr$year)[i], species=spplist[1:graph], perc[1:graph],disp=splist$dispersalmode[1:graph])  
    if (i==1) allresults=results else allresults=rbind(allresults,results)
    if (i==5) mtext(side=2, text="number of seeds (%)", line=3.5, las=3,cex=2)
    
  }
  ## chi-squared analyses of biotic vs abiotic dispersal dominance
  abiotic<-c(4,4,3,5,5,4,3,4,4,4)
  biotic<-c(1,1,2,0,0,1,2,1,1,1)
  chisq.test(c(40,10),p=c(0.47,0.53))
  dev.off()
  return(allresults)
  
}

#barplot(oneyr$perc[orspp][1:5],names.arg=oneyr$species[orspp][1:5], las=2, main=unique(sumyr$year)[i])


#### FIGURE 7 ################################
## NOT INCLUDED IN THE LAST VERSION OF THE PAPER

#figure 7 includes three graphs: 1) values of fruit biomass; 2) number of fruiting species; 
#3)summed values of our model of seed production vs. MEI 
#figure7(biomass=biomass,estfile="number total spp per month estimated.txt",k=3,file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", beginyearfile=beginyearfile, spstart=1, spend=45)
figure7=function(biomass=biomass,estfile="number total spp per month estimated.txt",k=3,file="Nouragues model all spp.txt",fileno="NRG model all spp without Dj.txt", beginyearfile=beginyearfile, spstart=1, spend=45)
  
{
  par(mfrow=c(3,1),las=1,mar=c(2,3,2,3),oma=c(2,2,2,2),cex=1)
  
  biomass$fecha=strptime(paste(biomass$Year2,"-",biomass$month,"-",1), format="%Y - %m - %d")
  start=which(biomass$Year2==2001 & biomass$month==2)
  end=which(biomass$Year2==2011 & biomass$month==2)
  N=length(biomass$fruitsb[start:end])
  biomass=biomass[start:end,]
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))
  
  plot(biomass$fruitsb,type="n", xlab="time", ylab="", ylim=c(0,6),bty="l", axes=F, las=1)
  axis(side=2, las=1,col="black")
  jan1=which(biomass$month==1)
  #axis(side=1, ,at=c(jan1[1]-12,jan1), labels=as.character(c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")))
  axis(side=1, at=c(jan1[1]-12,jan1), labels=F)
  mtext(side=2, "fruit biomass (g/m?)", line=3, cex=1.1,las=0)
  lines(c(start:(end-2)),runningmean(biomass$fruitsb,k), col="blue3", lwd=3)
  legend(end-125,8, lty=c(1,1,1,2), c("fruit biomass","# fruiting species","# seeds","MEI"), bty="n",col=c("blue3", "purple","black" ,"red"), horiz=T,lwd=2,xpd=T)
  #lines(biomass$fecha[start:(end-2)],runningmean(biomass$flowersb[start:end],k)/(160*0.5), lty=2)
  
  startmei=which(ensoanomalies$YEAR==2001 & ensoanomalies$MONTH==2)
  endmei=which(ensoanomalies$YEAR==2011 & ensoanomalies$MONTH==2)
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2)
  polygon(c(16,25,25,16),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(42,48,48,42),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(68,71,71,68),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(102,111,111,102),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(59,62,62,59),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(79,89,89,79),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(94,98,98,94),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(114,120,120,114),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, "MEI", line=3, cex=1.1, las=0)
  
  
  est=read.delim(estfile)
  est$fecha=strptime(paste(est$year,"-",est$month,"-",1), format="%Y - %m - %d")
  start=which(est$year==2001 & est$month==2)
  end=which(est$year==2011 & est$month==2)
  N=length(biomass$fruitsb[start:end])
  est=est[start:end,]
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))  
  
  plot(est$estnumbspp,type="n", xlab="time", ylab="", ylim=c(0,45),bty="l", axes=F, las=1)
  axis(side=2, las=1, col="black")
  jan1=which(est$month==1)
  #axis(side=1, at=c(jan1[1]-12,jan1), labels=as.character(c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011")))
  axis(side=1, at=c(jan1[1]-12,jan1), labels=FALSE)
  mtext(side=2, "number of fruiting species", line=3, cex=1.1,las=0)
  lines(c(start:(end-2)),runningmean(est$estnumbspp,k), col="purple", lwd=3)
  #legend(end-103,5, lty=c(1,2), c("number of fruiting species","MEI"), col=c("blue","red"),bty="n", cex=0.75, horiz=T,lwd=2)
  #lines(biomass$fecha[start:(end-2)],runningmean(biomass$flowersb[start:end],k)/(160*0.5), lty=2)
  
  startmei=which(ensoanomalies$YEAR==2001 & ensoanomalies$MONTH==2)
  endmei=which(ensoanomalies$YEAR==2011 & ensoanomalies$MONTH==2)
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2)
  polygon(c(16,25,25,16),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(42,48,48,42),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(68,71,71,68),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(102,111,111,102),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(59,62,62,59),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(79,89,89,79),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(94,98,98,94),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(114,120,120,114),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, "MEI", line=3, cex=1.1,las=0)  
  
  nrg=monthlyvalues(file=file,fileno=fileno, beginyearfile=beginyearfile,spstart=spstart,spend=spend)
  species=unique(nrg$species)
  sphigh=unique(nrg$sp[which.max(nrg$model)])
  onesphigh=nrg[nrg$species==sphigh,]
  onesphigh=onesphigh[onesphigh$julian>=15355&onesphigh$julian<=18246,]
  y2001=as.numeric(tojulian(c("01/01/2001","02/01/2001","03/01/2001","04/01/2001","05/01/2001","06/01/2001","07/01/2001","08/01/2001","09/01/2001","10/01/2001","11/01/2001","12/01/2001")))
  y2010=as.numeric(tojulian(c("01/01/2010", "02/01/2010","03/01/2010","04/01/2010","05/01/2010","06/01/2010","07/01/2010","08/01/2010","09/01/2010","10/01/2010","11/01/2010","12/01/2010","01/01/2011","02/01/2011")))
  tt=create.fulldate(fromjulian(c(y2001,onesphigh$julian, y2010), dateform="%Y-%m-%d"))
  plot(c(y2001,onesphigh$julian, y2010), c(rep(0,12),onesphigh$model, rep(0,14)),type="n", las=1, ylab="",xlab="time",bty="l", ylim=c(0,2000), axes=FALSE)
  #jan1=which(onesphigh$month==1)
  jan1=which(tt$month==1)
  axis(side=1,  at=c(y2001,onesphigh$julian, y2010)[jan1], labels=as.character(c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011")))
  axis(side=2, las=1)
  mtext(side=2, "total number of seeds",line=3, cex=1.1,las=0)
  polygon(c(onesphigh$julian[5], onesphigh$julian[14], onesphigh$julian[14], onesphigh$julian[5]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  polygon(c(onesphigh$julian[31], onesphigh$julian[37], onesphigh$julian[37], onesphigh$julian[31]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  polygon(c(onesphigh$julian[57], onesphigh$julian[60], onesphigh$julian[60], onesphigh$julian[57]),c(0,0,2000,2000),col = "#FF000050", border = NA)  
  #polygon(c(onesphigh$julian[91], onesphigh$julian[96], onesphigh$julian[96], onesphigh$julian[91]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  polygon(c(onesphigh$julian[48], onesphigh$julian[51], onesphigh$julian[51], onesphigh$julian[48]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  polygon(c(onesphigh$julian[68], onesphigh$julian[78], onesphigh$julian[78], onesphigh$julian[68]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  polygon(c(onesphigh$julian[83], onesphigh$julian[87], onesphigh$julian[87], onesphigh$julian[83]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  polygon(c(c(y2001,onesphigh$julian, y2010)[115],c(y2001,onesphigh$julian, y2010)[121],c(y2001,onesphigh$julian, y2010)[121],c(y2001,onesphigh$julian, y2010)[115]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  polygon(c(c(y2001,onesphigh$julian, y2010)[103],c(y2001,onesphigh$julian, y2010)[112],c(y2001,onesphigh$julian, y2010)[112],c(y2001,onesphigh$julian, y2010)[103]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  
  for (i in 1:length(species)){
    onesp=nrg[nrg$species==species[i],]
    onesp=onesp[onesp$julian>=15355&onesp$julian<=18246,]
    lines(onesp$julian, onesp$model, col=i, lwd=2)
  }
  sumseed=aggregate(data.frame(peak=nrg$model), by=list(month=nrg$month,year=nrg$year ), sum)
  sumseed=sumseed[sumseed$year>=2002 & sumseed$year < 2010,]  
  sumseed$julian=tojulian(paste("15/",sumseed$month,"/",sumseed$year), dateform="%d/ %m / %Y")
  lines(sumseed$julian, sumseed$peak, col="black",lwd=3)
  startmei=which(ensoanomalies$YEAR==2001 & ensoanomalies$MONTH==2)
  endmei=which(ensoanomalies$YEAR==2011 & ensoanomalies$MONTH==2)
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F,las=2, col="red", lwd=2)
  abline(h=0, lty=2)
  axis(side=4, las=2, col="red")
  mtext(side=4, "MEI", line=3, las=0, cex=1.1)
  
}



####FIGURE 7b OF THE PAPER ###############
#figure7b function only includes the graph of biomass and estimated number of species against MEI
figure7b=function(biomass=biomass,estfile="number total spp per month estimated.txt", k=3)
  
{
  par(mfrow=c(2,1),las=1,mar=c(2,3,2,3),oma=c(2,2,2,2),cex=1)
  
  biomass$fecha=strptime(paste(biomass$Year2,"-",biomass$month,"-",1), format="%Y - %m - %d")
  start=which(biomass$Year2==2001 & biomass$month==2)
  end=which(biomass$Year2==2011 & biomass$month==2)
  N=length(biomass$fruitsb[start:end])
  biomass=biomass[start:end,]
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))
  
  plot(biomass$fruitsb,type="n", xlab="time", ylab="", ylim=c(0,6),bty="l", axes=F, las=1)
  axis(side=2, las=1,col="black")
  jan1=which(biomass$month==1)
  #axis(side=1, ,at=c(jan1[1]-12,jan1), labels=as.character(c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")))
  axis(side=1, at=c(jan1[1]-12,jan1), labels=F)
  mtext(side=2, "fruit biomass (g/m?)", line=3, cex=1.1,las=0)
  lines(c(start:(end-2)),runningmean(biomass$fruitsb,k), col="blue3", lwd=2)
  legend(end-125,7.5, lty=c(1,1,2), c("fruit biomass","number of fruiting species","MEI"), bty="n",col=c("blue3", "black", "red"), horiz=T,lwd=2,xpd=T)
  #lines(biomass$fecha[start:(end-2)],runningmean(biomass$flowersb[start:end],k)/(160*0.5), lty=2)
  
  startmei=which(ensoanomalies$YEAR==2001 & ensoanomalies$MONTH==2)
  endmei=which(ensoanomalies$YEAR==2011 & ensoanomalies$MONTH==2)
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2)
  polygon(c(16,25,25,16),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(42,48,48,42),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(68,71,71,68),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(102,111,111,102),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(59,62,62,59),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(79,89,89,79),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(94,98,98,94),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(114,120,120,114),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, "MEI", line=3, cex=1.1, las=0)
  
  
  est=read.delim(estfile)
  est$fecha=strptime(paste(est$year,"-",est$month,"-",1), format="%Y - %m - %d")
  start=which(est$year==2001 & est$month==2)
  end=which(est$year==2011 & est$month==2)
  N=length(biomass$fruitsb[start:end])
  est=est[start:end,]
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))  
  
  plot(est$estnumbspp,type="n", xlab="time", ylab="", ylim=c(0,45),bty="l", axes=F, las=1)
  axis(side=2, las=1, col="black")
  jan1=which(est$month==1)
  axis(side=1, at=c(jan1[1]-12,jan1), labels=as.character(c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011")))
  mtext(side=2, "number of fruiting species", line=3, cex=1.1,las=0)
  lines(c(start:(end-2)),runningmean(est$estnumbspp,k), col="black", lwd=2)
  #legend(end-103,5, lty=c(1,2), c("number of fruiting species","MEI"), col=c("blue","red"),bty="n", cex=0.75, horiz=T,lwd=2)
  #lines(biomass$fecha[start:(end-2)],runningmean(biomass$flowersb[start:end],k)/(160*0.5), lty=2)
  
  startmei=which(ensoanomalies$YEAR==2001 & ensoanomalies$MONTH==2)
  endmei=which(ensoanomalies$YEAR==2011 & ensoanomalies$MONTH==2)
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2)
  polygon(c(16,25,25,16),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(42,48,48,42),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(68,71,71,68),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(102,111,111,102),c(-2.5,-2.5,6,6),col = "#FF000050", border = NA)
  polygon(c(59,62,62,59),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(79,89,89,79),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(94,98,98,94),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  polygon(c(114,120,120,114),c(-2.5,-2.5,6,6),col = "#0000FF50", border = NA)
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, "MEI", line=3, cex=1.1,las=0)  
  
}

###### FIGURE 7c ####################################################

#figure7c(file="NRG model all spp without Dj.txt",fileno="NRG model all spp without Dj.txt", beginyearfile=beginyearfile, spstart=1, spend=45, k=3)
#this function plots the summed number of seeds of the models vs. MEI 
figure7c=function(file="NRG model all spp without Dj.txt",fileno="NRG model all spp without Dj.txt", beginyearfile=beginyearfile, spstart=1, spend=45, k=3)
{
  nrg=monthlyvalues(file=file,fileno=fileno, beginyearfile=beginyearfile,spstart=spstart,spend=spend)
  species=unique(nrg$species)
  
  #go=which(NRGallspp$sp=="Gouania blanchetiana")#this removes Gouania blanchetiana because it is the most abundant spp and marks all the fruiting pattern
  #NRGallspp=NRGallspp[-go,]
  sumseed=aggregate(data.frame(peak=nrg$model), by=list(month=nrg$month,year=nrg$year ), sum)
  sumseed=sumseed[sumseed$year>=2002 & sumseed$year < 2010,]  
  sumseed$julian=tojulian(paste("15/",sumseed$month,"/",sumseed$year), dateform="%d/ %m / %Y")
  N=length(sumseed$julian)
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))
  
  par(mar=c(5,5,3,4), cex=1.5,las=1)
  sphigh=unique(nrg$sp[which.max(nrg$model)])
  onesphigh=nrg[nrg$species==sphigh,]
  onesphigh=onesphigh[onesphigh$julian>=15355&onesphigh$julian<=18246,]
  plot(onesphigh$julian, onesphigh$model, type="n", las=1, ylab="",axes=F,xlab="time",bty="l", ylim=c(0,2000))
  jan1=which(sumseed$month==1)
  axis(side=1, at=onesphigh$julian[c(jan1, jan1[length(jan1)]+11)], labels=as.character(c("2002","2003","2004","2005","2006","2007","2008","2009","2010")))
  mtext(side=2, "monthly number of seeds", line=3, cex=1.5,las=0)
  axis(side=2)
  
  polygon(c(onesphigh$julian[5], onesphigh$julian[14], onesphigh$julian[14], onesphigh$julian[5]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  polygon(c(onesphigh$julian[31], onesphigh$julian[37], onesphigh$julian[37], onesphigh$julian[31]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  polygon(c(onesphigh$julian[57], onesphigh$julian[60], onesphigh$julian[60], onesphigh$julian[57]),c(0,0,2000,2000),col = "#FF000050", border = NA)  
  polygon(c(onesphigh$julian[91], onesphigh$julian[96], onesphigh$julian[96], onesphigh$julian[91]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  polygon(c(onesphigh$julian[48], onesphigh$julian[52], onesphigh$julian[52], onesphigh$julian[48]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  polygon(c(onesphigh$julian[68], onesphigh$julian[78], onesphigh$julian[78], onesphigh$julian[68]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  polygon(c(onesphigh$julian[83], onesphigh$julian[87], onesphigh$julian[87], onesphigh$julian[83]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  
  lines(sumseed$julian[comienzo:(final-2)],runningmean(sumseed$peak[comienzo:final],k),col="blue", lwd=2)
  legend(sumseed$julian[final-90],2180, lty=c(1,2), col=c("blue", "red"),c("number of seeds","MEI"), bty="n", cex=1, horiz=T, lwd=2)
  #lines(biomass$fecha[start:(end-2)],runningmean(biomass$flowersb[start:end],k)/(160*0.5), lty=2)
  
  startmei=which(ensoanomalies$YEAR==2002 & ensoanomalies$MONTH==1)
  endmei=which(ensoanomalies$YEAR==2009 & ensoanomalies$MONTH==12)
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F,las=2, col="red", lwd=2)
  abline(h=0, lty=2)
  axis(side=4, las=2)
  mtext(side=4, "MEI", line=3, las=0, cex=1.5)
  
}


#### FIGURE 7d ####################################################
##Figure 7d includes a line per species for the monthly model of seed production

#figure7d(file="NRG model all spp without Dj.txt",fileno="NRG model all spp without Dj.txt", beginyearfile=beginyearfile, spstart=1, spend=45)
figure7d=function(file="Nouragues model all spp.txt",fileno="NRG model all spp without Dj.txt", beginyearfile=beginyearfile, spstart=1, spend=45)
{
  par(las=1,mar=c(2,3,2,3),oma=c(2,2,2,2),cex=1)
  nrg=monthlyvalues(file=file,fileno=fileno, beginyearfile=beginyearfile,spstart=spstart,spend=spend)
  species=unique(nrg$species)
  sphigh=unique(nrg$sp[which.max(nrg$model)])
  onesphigh=nrg[nrg$species==sphigh,]
  onesphigh=onesphigh[onesphigh$julian>=15355&onesphigh$julian<=18246,]
  plot(onesphigh$julian, onesphigh$model, type="n", las=1, ylab="",xlab="time",bty="l", ylim=c(0,2000), axes=FALSE)
  jan1=which(onesphigh$month==1)
  axis(side=1,  at=onesphigh$julian[c(jan1, jan1[length(jan1)]+11)], labels=as.character(c("2002","2003","2004","2005","2006","2007","2008","2009","2010")))
  axis(side=2, las=1)
  mtext(side=2, "total number of seeds",line=3, cex=1.1,las=0)
  polygon(c(onesphigh$julian[5], onesphigh$julian[14], onesphigh$julian[14], onesphigh$julian[5]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  polygon(c(onesphigh$julian[31], onesphigh$julian[37], onesphigh$julian[37], onesphigh$julian[31]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  polygon(c(onesphigh$julian[57], onesphigh$julian[60], onesphigh$julian[60], onesphigh$julian[57]),c(0,0,2000,2000),col = "#FF000050", border = NA)  
  polygon(c(onesphigh$julian[91], onesphigh$julian[96], onesphigh$julian[96], onesphigh$julian[91]),c(0,0,2000,2000),col = "#FF000050", border = NA)
  polygon(c(onesphigh$julian[48], onesphigh$julian[52], onesphigh$julian[52], onesphigh$julian[48]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  polygon(c(onesphigh$julian[68], onesphigh$julian[78], onesphigh$julian[78], onesphigh$julian[68]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  polygon(c(onesphigh$julian[83], onesphigh$julian[87], onesphigh$julian[87], onesphigh$julian[83]),c(0,0,2000,2000),col = "#0000FF50", border = NA)
  
  
  for (i in 1:length(species)){
    onesp=nrg[nrg$species==species[i],]
    onesp=onesp[onesp$julian>=15355&onesp$julian<=18246,]
    lines(onesp$julian, onesp$model, col=i, lwd=2)
  }
  sumseed=aggregate(data.frame(peak=nrg$model), by=list(month=nrg$month,year=nrg$year ), sum)
  sumseed=sumseed[sumseed$year>=2002 & sumseed$year < 2010,]  
  sumseed$julian=tojulian(paste("15/",sumseed$month,"/",sumseed$year), dateform="%d/ %m / %Y")
  lines(sumseed$julian, sumseed$peak, col="black",lwd=3)
  startmei=which(ensoanomalies$YEAR==2002 & ensoanomalies$MONTH==1)
  endmei=which(ensoanomalies$YEAR==2009 & ensoanomalies$MONTH==12)
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F,las=2, col="red", lwd=2)
  abline(h=0, lty=2)
  axis(side=4, las=2, col="red")
  mtext(side=4, "MEI", line=3, las=0, cex=1.1)
  
}


###### FIGURE APPENDIX 1 ####################################################

appendix1 = function(file = nourage, fit= results, beginyearfile = beginyearfile)
{
  high=c(5,33,35,37)
  low=c(4,11,18,19,23,25,26,28,32,34,40,45)
  special=sort(c(high,low))
  resto=c(1:45)[-special]
  
  nrg.graphierspp2(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=1,spend=3, filename="appendix1_1.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.low(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(4), filename="appendix1-2b.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.high(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(5), filename="appendix1-3.pdf",longnames="total number of seeds per species.txt")
  nrg.graphierspp2(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=6,spend=17, filename="appendix1-4.pdf",longnames="total number of seeds per species.txt")  
  nrg.graphsmoved.low(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(18,19), filename="appendix1-5.pdf",longnames="total number of seeds per species.txt")
  nrg.graphierspp2(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=20,spend=22, filename="appendix1-6.pdf",longnames="total number of seeds per species.txt")  
  nrg.graphsmoved.low(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(23), filename="appendix1-7.pdf",longnames="total number of seeds per species.txt")
  nrg.graphierspp2(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=24,spend=24, filename="appendix1-8.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.low(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(25,26), filename="appendix1-9.pdf",longnames="total number of seeds per species.txt")
  nrg.graphierspp2(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=27,spend=27, filename="appendix1-10.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.low(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(28), filename="appendix1-11.pdf",longnames="total number of seeds per species.txt")
  nrg.graphierspp2(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=29,spend=31, filename="appendix1-12.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.low(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(32), filename="appendix1-13.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.high(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(33), filename="appendix1-14.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.low(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(34), filename="appendix1-15.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.high(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(35), filename="appendix1-16.pdf",longnames="total number of seeds per species.txt")
  nrg.graphierspp2(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=36,spend=36, filename="appendix1-17.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.high(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(37), filename="appendix1-18.pdf",longnames="total number of seeds per species.txt")
  nrg.graphierspp2(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=38,spend=39, filename="appendix1-19.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.low(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(40), filename="appendix1-20.pdf",longnames="total number of seeds per species.txt")
  nrg.graphierspp2(file=nourage, fit=results, beginyearfile=beginyearfile,spstart=41,spend=44, filename="appendix1-21.pdf",longnames="total number of seeds per species.txt")
  nrg.graphsmoved.low(file=nourage, fit=results, beginyearfile=beginyearfile,spnumber=c(45), filename="appendix1-22.pdf",longnames="total number of seeds per species.txt")
}

#### APPENDIX 2 (for the ATBC 2016 talk) ###############  
##What was the trend of peak parameter?

appendix2 = function(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", beginyearfile=beginyearfile, spstart=1, spend=45,graph=5,longnames="total number of seeds per species.txt",startyear=2001,endyear=2010,filename="appendix2.tiff"){
  
  data=monthlyvalues(file=file,fileno=fileno, beginyearfile=beginyearfile,spstart=spstart,spend=spend)
  totseed=read.delim(file=longnames)
  names(totseed)=c("species","longname","totseed","form","disp","fruit")
  dispersal=data.frame(species=totseed$longname,disp=totseed$disp)
  datayr<-merge(data,totseed, by="species",all.x=T)
  sumyr=aggregate(data.frame(sumseed=datayr$model), by=list(species=datayr$longname, year=datayr$year), sum)
  filno=read.delim(file=fileno)
  uniquefilno<-aggregate(data.frame(peak=filno$peak), by=list(species=filno$sp, year=filno$year), unique)
  sndmean=aggregate(data.frame(mean=uniquefilno$peak),by=list(species=uniquefilno$species),mean)
  sndsd=aggregate(data.frame(sd=uniquefilno$peak),by=list(species=uniquefilno$species),sd)
  
  snd=numeric() 
  for (i in 1:nrow(uniquefilno)){
    uniquefilno$snd[i]= (log(uniquefilno$peak[i]+1,10)-log(sndmean$mean[sndmean$species==uniquefilno$species[i]]+1,10))/log(sndsd$sd[sndsd$species==uniquefilno$species[i]]+1,10)
  }  
  ds=data.frame(years=startyear:endyear,mean=aggregate(uniquefilno$snd,by=list(year=uniquefilno$year),mean)[,2],se=as.numeric(aggregate(uniquefilno$snd,by=list(year=uniquefilno$year),sd)[,2]/sqrt(lengthunique(uniquefilno$species))))
  ggplot(ds, aes(x=years,y=mean))+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1,position = "identity")+
    geom_line() +
    geom_point(size=2, shape=21, fill="black")+
    xlab("years") +
    ylab("Seed production (SND)") +
    #expand_limits(y=0) +                        # Expand y range
    scale_y_continuous(breaks=seq(-0.5,0.5,0.1)) +         # Set tick every 4
    scale_x_continuous(breaks=seq(startyear,endyear,1)) +
    theme_classic()+
    theme(text = element_text(size=20),axis.line = element_line(lineend="square"),axis.text.x = element_text(angle=90, vjust=1),axis.ticks.length = unit(.25, "cm"))
  ggsave(filename,width=20, height=15, units="cm")
  
  
  #tiff(filename=filename,height=1000,width=2000,res=300)
  #par(mfrow=c(1,1),las=1,mar=c(3,3,1,0),oma=c(2,2,1,1),cex=0.6)
  #meanyr=aggregate(data.frame(mean=uniquefilno$peak),by=list(year=uniquefilno$year),mean)
  #plot(meanyr$year, meanyr$mean,type="o")
  #plot(uniquefilno$year,uniquefilno$peak)  
  #dev.off()
}  

#### SUPPLEMENTARY TABLE 1 ####

stable1 = function(seed.data = "seed size.txt"){
ss = read.delim(file = seed.data)
ss_df = group_by(tbl_df(ss), Taxa, Family)

Smythe <- character(length = dim(ss)[1])
Smythe[which(ss$ave.L > 15)] <- "TypeI"
Smythe[which(ss$ave.L <= 15 & ss$Dispersal.syndrome == "zoo")] <- "TypeII"

ss = data.frame(ss, Smythe)
write.table(ss, file = "supplementary table 1.txt", row.names = F, sep = "\t")
return(ss)
}

#### SUPPLEMENTARY TABLE 2 ####

stable2 = function(file = "Nouragues results hyperparameters.txt", longnames = "total number of seeds per species.txt", filename = "supplementary Table 2.txt"){
  
  tr = read.delim(file)
  tr$mu2 = tr$mu
  tr$mu2[which(tr$mu <= 0)] = tr$mu[which(tr$mu <= 0)] + 365.25
  tr$mu2[which(tr$mu > 365)] = tr$mu[which(tr$mu > 365)] - 365.25
  
  
  months = list(jan = 1:31, feb = 32:59, mar = 60:90, apr = 91:120,
              may = 121:151, jun = 152:181, jul = 182:212, aug = 213:243, sep = 244:273, oct = 274:304, nov = 305:334, dec = 335:365)
  rmonths=unlist(months)
  month = names(rmonths[trunc(tr$mu2)]) 
  
  tr$CImu2.2 = tr$CImu2
  tr$CImu2.2[which(tr$CImu2 <= 0)] = tr$CImu2[which(tr$CImu2 <= 0)]+365.25
  tr$CImu2.2[which(tr$CImu2 > 365)] = tr$CImu2[which(tr$CImu2 > 365)]-365.25
  month2 <- names(rmonths[trunc(tr$CImu2.2)])  
  
  tr$CImu97.2 = tr$CImu97
  tr$CImu97.2[which(tr$CImu97 <= 0)] = tr$CImu97[which(tr$CImu97 <= 0)]+365.25
  tr$CImu97.2[which(tr$CImu97 > 365)] = tr$CImu97[which(tr$CImu97 > 365)]-365.25
  month97 <- names(rmonths[trunc(tr$CImu97.2)])
  
  tr$pmu = exp(tr$logmu)
  tr$pCI2 = exp(tr$CIlogmu2)
  tr$pCI97 = exp(tr$CIlogmu97)
  tr$pSD = exp(tr$logSD)
  tr$pSDCI2 = exp(tr$CIlogSD2)
  tr$pSDCI97 = exp(tr$CIlogSD97)
  
  totseed = read.delim(file = longnames)
  names(totseed) = c("species", "longname", "totseed", "form", "disp", "fruit", "length", "width", "Smythe")
  long = merge(tr, totseed, by ="species")
  dataset = data.frame(species = long$longname, peakdays = paste(round(long$mu2,2)," (",month,")",sep=""), CIpeakday = paste(round(long$CImu2.2,2),"-" ,round(long$CImu97.2,2), " (", month2, "-", month97, ")", sep =""), SDpeakday =  round(long$SD,2), CISDpeakday = paste(round(long$CISD2,2),"-", round(long$CISD97,2), sep =""), peakmu = round(long$pmu,2), CIpeakmu = paste(round(long$pCI2,2),"-", round(long$pCI97,2), sep =""), peakSD = round(long$pSD, 2), CIpeakSD = paste(round(long$pSDCI2,2),"-", round(long$pSDCI97,2), sep =""), SD = round(long$bestSD, 2))
  write.table(dataset, file = filename, sep = "\t", row.names = F)
}

#### SUPPLEMENTARY TABLE 3 ##################
#file = nourage; beginyearfile = beginyearfile; spstart = 1; spend = 45; fit = results; dj=FALSE; filename = "Supplementary Table 3.txt"
stable3 = function(file = nourage, beginyearfile = beginyearfile, spstart = 1, spend = 45, fit = results, dj = TRUE, filename = "Supplementary Table 3.txt")
{
  parameters = parametersyr(file = file, beginyearfile = beginyearfile, spstart = spstart, spend = spend, fit = results, dj = dj)
  pp = aggregate(data.frame(peakday = parameters$peakday,  CI2peakday=parameters$CI2peakday,  CI97peakday=parameters$CI97peakday,peak=parameters$peak,  CI2peak=parameters$CI2peak,  CI97peak=parameters$CI97peak), by=list(year=parameters$cycle,sp=parameters$sp),mean)
  yrmin = aggregate(data.frame(min=parameters$year), by=list(year=parameters$cycle,sp=parameters$sp),min)
  yrmax = aggregate(data.frame(max=parameters$year), by=list(year=parameters$cycle,sp=parameters$sp),max)
  months = list(jan = 1:31, feb = 32:59, mar = 60:90, apr = 91:120,
                may = 121:151, jun = 152:181, jul = 182:212, aug = 213:243, sep = 244:273, oct = 274:304, nov = 305:334, dec = 335:365)
  rmonths = unlist(months)
 
  
  peak = ifelse(pp$peakday< 0, pp$peakday+365, pp$peakday)
  peak2 = ifelse(peak > 365, peak - 365, peak)
  month = names(rmonths[trunc(peak2)]) 
  
  CI2peak = ifelse(pp$CI2peakday < 0, pp$CI2peakday + 365, pp$CI2peakday)
  CI2peak2 = ifelse(CI2peak > 365, CI2peak - 365, CI2peak)
  CI2peak3 = ifelse(CI2peak2 < 1 & CI2peak2 > 0, 1, CI2peak2)
  month2 = names(rmonths[trunc(CI2peak3)])
  
  CI97peak = ifelse(pp$CI97peakday < 0, pp$CI97peakday + 365, pp$CI97peakday)
  CI97peak2 = ifelse(CI97peak > 365, CI97peak - 365, CI97peak)
  CI97peak3 = ifelse(CI97peak2 < 1 & CI97peak2 > 0, 1, CI97peak2)
  month97 = names(rmonths[trunc(CI97peak3)])
  
  pp2 = data.frame(sp = pp$sp, cycle = pp$year, years = paste(yrmin$min,"-",yrmax$max), peakdays = paste(round(peak2,2)," (",month,")",sep=""), CIpeakday = paste(round(CI2peak2,2),"-" , round(CI97peak2,2), " (", month2, "-", month97, ")", sep =""), peak = round(pp$peak,2), CIpeak = paste(round(pp$CI2peak,2),"-", round(pp$CI97peak,2)))
  pp3 = data.frame(sp=pp$sp, cycle=pp$year,year1=yrmin$min, year2=yrmax$max, peakdays=pp$peakday, CI2peakday= pp$CI2peakday, CI97peakday=pp$CI97peakday, peak=pp$peak, CI2peak=pp$CI2peak, CI97peak=pp$CI97peak)
  
  write.table(pp2, file = filename, sep = "\t", row.names = F)
  return(pp2)
 
}

#### SUPPLEMENTARY FIGURE 1 ################

SFig1 = function(file = "Nouragues results hyperparameters.txt", graphname = "SFig1.tif") {
  
  tr <- read.delim(file)
  sp<- read.delim("total number of seeds per species.txt")
  sp2 = data.frame(species = sp$sp, totseed = sp$totseed)
  sp3<- merge(sp2,tr, by = "species")
  plot(log(sp3$totseed), sp3$logSD)
  cor(log(sp3$totseed), sp3$logSD)
  plot(log(sp3$totseed), sp3$SD)
  cor(log(sp3$totseed), sp3$SD)
  tiff(filename = graphname,width = 1500, height = 1000,pointsize=12, res=300)
  par(las = 1, bty = "o", tcl = 0.2, mar = c(5, 5, 2,2), mgp = c(0.25, 0.25, 0),cex.axis=1.2,lwd=1.5)
  plot(tr$logmu, tr$SD,xlab="",ylab="",las=1,bty="l",pch=19)
  plot(tr$logSD, tr$SD,xlab="",ylab="",las=1,bty="l",pch=19)
  mtext(side=2,text = expression(paste (sigma, " of peakday", sep ="")),line=2.5,las=0,cex=1.2)
  mtext(side=1,text = expression(paste ("log", sigma, " of P", sep ="")),line=2.5,las=0,cex=1.2)
  abline(lm(tr$SD ~ tr$logSD))
  summary(lm(tr$SD ~ tr$logSD))
  cor.test(tr$SD,tr$logSD,method="pearson")
  
  dev.off()
}


#### SUPPLEMENTARY FIGURE 3 ################
#relationship between seed size and variation of peakday (hyper-SD of peakday)

supplementary3 = function(file = "Nouragues results hyperparameters.txt", seedfile = "seed size.txt", graphname = "SFig3.tif") {
  
  tr <- read.delim(file)
  ssize <- read.delim(seedfile)
  
  condensed <- CVpeakdays()
  combined <- merge(tr, ssize, by = "species")
  #combined <- combined[-4,] 
  type = character(length = dim(combined)[1])
  combined$type[which(combined$ave.L >= 15)] <- "Type 1"
  combined$type[which(combined$ave.L < 15)] <- "Type 2"
  boxplot(combined$logSD ~ combined$type)
  
  #combined <- merge(condensed, ssize, by = "species")
  
  tiff(filename = graphname, width = 1500, height = 1000, pointsize = 12, res = 300)
  par(las = 1, bty = "o", tcl = 0.2, mar = c(5, 5, 2,2), mgp = c(0.25, 0.25, 0),cex.axis=1.2,lwd=1.5)
  plot(combined$ave.L, combined$SD, xlab="", ylab="", las=1, bty="l", pch=19)
  mtext(side = 1,text="seed length (mm)",line=2.5,las=0,cex=1.2)
  mtext(side = 2,text= "SD of peakday",line=2.5,las=0,cex=1.2)
  abline(lm(combined$SD ~ combined$ave.L))
  summary(lm(combined$SD ~ combined$ave.L))
  cor.test(combined$ave.L, combined$SD, method="pearson")
  
  dev.off()
}

#### MISCELLANEOUS ####

### Smythe´s hypothesis###
smythe = function(file = "nouragues results parameters per year.txt", longnames="total number of seeds per species.txt") {
  
  nrg = read.delim(file)
  #spnames=sort(unique(nrg$sp))
  totseed = read.delim(file = longnames)
  names(totseed) = c("species", "longname", "totseed", "form", "disp", "fruit", "length", "width", "Smythe")
  spnames2 = totseed$longname
  sptot = aggregate(data.frame(ptot = nrg$peak),by = list(species = nrg$species),sum)
  
  condensed = merge(nrg, sptot, by = "species", all.x = T)
  condensed$p = condensed$peak/condensed$ptot
  condensed$Di <- condensed$p*log(condensed$p)
  
  spD <- aggregate(data.frame(spD = condensed$Di),by = list(species = condensed$species),sum)
  spD$maxdiv <- -spD$spD/log(10)
  
  condensed2 = merge(totseed, spD, by="species")
  type = character(length = dim(condensed2)[1])
  condensed2$type[which(condensed2$length >= 15)] <- "Type 1"
  condensed2$type[which(condensed2$length < 15)] <- "Type 2"
  condensed2 <- condensed2[- which(condensed2$species == "Bombacaceae sp1"),]
  orco = order(condensed2$maxdiv)
  
  ggplot(condensed2, aes(length, maxdiv)) + 
    geom_point() +
    geom_smooth(aes(x = length, y = maxdiv), span = 0.4, method = "loess", se = F,  col = "black", alpha = 0.4) 
  
  #abline(lm(condensed2$maxdiv ~ condensed2$length))
  summary(lm(condensed2$maxdiv ~ condensed2$length))
  summary(glm(condensed2$maxdiv ~ condensed2$length), family = "binomial")
  boxplot(condensed2$maxdiv ~ condensed2$type)
  lm1 <- lm(condensed2$maxdiv ~ condensed2$type)
  summary(lm1)
}

#stats with table 1: how many species were masting; do fruit variability change according to the dispersal mode?#

table1stats <- function (table1 = "table 1.txt"){
  tab1 <- read.delim(table1)
  table(tab1$masting)
  (18/45)
}

# using the bestSD as proxy of fit of the model

fitmodel = function(file = "Nouragues results hyperparameters.txt"){
  
  tr <-read.delim(file)
  bestSD <- tr$bestSD
  hist(bestSD)
  boxplot(bestSD)
  CI(bestSD)
  tr[which(bestSD >=94),]
  
}

## calculating CV for Peak
#load("hyper all species Nouragues.RData")

CVlog <- function (hyper = hyper){
  acacia = hyper[hyper$species == "Acacia tenuifolia",]
  CV = geommean = geomsd = vector(length = 9000)
  for (i in 1:9000) {
    lognormal <- rlnorm(1000, hyper$logmu[i], hyper$logSD[i])
    CV[i]= sd(lognormal)/mean(lognormal)
   geommean[i] = mean(lognormal)
   geomsd[i] = sd(lognormal)
  }
  logSD = mean(acacia$logmu)
  meancv = mean(CV)
  uppercv = as.numeric(CI(CV)[1])
  lowercv = as.numeric(CI(CV)[2])
  
}

#this function plots fruit biomass for Nouragues in relation to MEI (it doesn't work well)
biomasstime=function(biomass=biomass, k=3) {
  
  biomass$fecha=strptime(paste(biomass$Year2,"-",biomass$month,"-",1), format="%Y - %m - %d")
  start=which(biomass$Year2==2001 & biomass$month==2)
  end=which(biomass$Year2==2012 & biomass$month==1)
  N=length(biomass$fruitsb[start:end])
  biomass=biomass[start:end,]
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))
  
  par(mar=c(4,4,2,4),cex=1)
  plot(biomass$fruitsb,type="n", xlab="time", ylab="fruit biomass (g/m?)", ylim=c(0,6),bty="l", axes=F, las=1)
  axis(side=2, las=1,col="blue3")
  jan1=which(biomass$month==1)
  axis(side=1, at=c(jan1[1]-12,jan1), labels=as.character(c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")))
  lines(c(start:(end-2)),runningmean(biomass$fruitsb,k), col="blue3", lwd=2)
  legend(end-113,0.8, lty=c(1,2), c("fruit biomass","MEI"), bty="n",col=c("blue3", "red"), horiz=T,lwd=2)
  #lines(biomass$fecha[start:(end-2)],runningmean(biomass$flowersb[start:end],k)/(160*0.5), lty=2)

  startmei=which(ensoanomalies$YEAR==2001 & ensoanomalies$MONTH==2)
  endmei=which(ensoanomalies$YEAR==2012 & ensoanomalies$MONTH==1)
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2)
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, "MEI", line=3, cex=1.5)

  }

#esttime plots the estimated number of fruiting species vs MEI
esttime=function(estfile="number total spp per month estimated.txt", k=3) {
  
  est=read.delim(estfile)
  est$fecha=strptime(paste(est$year,"-",est$month,"-",1), format="%Y - %m - %d")
  start=which(est$year==2001 & est$month==2)
  end=which(est$year==2011 & est$month==2)
  N=length(biomass$fruitsb[start:end])
  est=est[start:end,]
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))
  
  par(mar=c(4,4,2,4), cex=1.5)
  plot(est$estnumbspp,type="n", xlab="time", ylab="number of fruiting species", ylim=c(0,45),bty="l", axes=F, las=1)
  axis(side=2, las=1, col="blue")
  jan1=which(est$month==1)
  axis(side=1, at=c(jan1[1]-12,jan1), labels=as.character(c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011")))
  lines(c(start:(end-2)),runningmean(est$estnumbspp,k), col="blue", lwd=2)
  legend(end-103,5, lty=c(1,2), c("number of fruiting species","MEI"), col=c("blue","red"),bty="n", cex=0.75, horiz=T,lwd=2)
  #lines(biomass$fecha[start:(end-2)],runningmean(biomass$flowersb[start:end],k)/(160*0.5), lty=2)
  
  startmei=which(ensoanomalies$YEAR==2001 & ensoanomalies$MONTH==2)
  endmei=which(ensoanomalies$YEAR==2011 & ensoanomalies$MONTH==2)
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2)
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, "MEI", line=3, cex=1.5)
  
}

#this graphs plots the values of a climatic variable (rainfall, tmin or tmax) versus MEI
#graphrainenso(clim=clim, ensoanomalies=ensoanomalies, climvariable="tmin", k=3)
graphrainenso=function(clim=clim,climvariable="rain", ensoanomalies=ensoanomalies, k=3, unit="(mm)"){
  
  ind=which(colnames(clim)==climvariable)
  rain=clim[,ind]
  
  par(mar=c(5,5,3,5))
  plot(clim$julian, rain,  ylab=paste(climvariable, unit), xlab="time", las=1, type="n", axes=T)
  lines(clim$julian[1:(length(clim$julian)-2)],runningmean(rain,k), col="black")
  startmei=which(ensoanomalies$YEAR==2001 & ensoanomalies$MONTH==1)
  endmei=which(ensoanomalies$YEAR==2011 & ensoanomalies$MONTH==7)
  legend(clim$julian[70],max(rain),lty=c(1,2), c(climvariable,"MEI"), bty="n")
  par(new=T)
  plot(ensoanomalies$MEI[startmei:endmei],type="l", lty=2,ylab="",xlab="",axes=F)
  abline(h=0, lty=2)
  axis(side=4, las=1)
  mtext(side=4, "MEI", line=3)
  
}


ccfbiomassmei=function(biomass, ensoanomalies){
  
  start=which(biomass$Year2==2001 & biomass$month==2)
  end=which(biomass$Year2==2011 & biomass$month==4)
  biomass=biomass[start:end,]
  biomassmean=aggregate(data.frame(fruit=biomass$fruitsb),by=list(month=biomass$month), mean)
  estmean=aggregate(data.frame(est=est$estnumbspp),by=list(month=est$month), mean)
  newbio=newest=vector()
  
  for ( i in 1:length(biomass$fruitsb))
  {
         #esto sirve para quitar a cada valor la media del mes
    newbio[i]<- biomass$fruitsb[i]-biomassmean$fruit[biomassmean$month==biomass$month[i]] 
  }
  
  for ( i in 1:length(biomass$fruitsb))
  {
    #esto sirve para quitar a cada valor la media del mes
    newest[i]<- est$estnumbspp[i]-estmean$est[estmean$month==est[i]] 
  }
  
  startmei=which(ensoanomalies$YEAR==2001 & ensoanomalies$MONTH==2)
  endmei=which(ensoanomalies$YEAR==2011 & ensoanomalies$MONTH==4)
  
  for ( i in 1:length(ensoanomalies$MEI[startmei:endmei]))###this is for making a correlation at 5 month lag
  {
    #esto sirve para quitar a cada valor la media del mes
    newbio[i]<- biomass$fruitsb[i]-biomassmean$fruit[biomassmean$month==biomass$month[i]] 
  }
  
  ccf(ensoanomalies$MEI[startmei:endmei],newbio)
  
}


crosscorrbioclimate<-function()
  
{
  est<-read.delim("number total spp per month estimated.txt",header=T)
  estt<-est[,c(1,2,4)]
  soi<-read.delim(file="ENSO anomalies - 2012.txt",as.is=T)
  mei<-soi[,c(1,2,8)]
  mei10<-subset(mei,mei$YEAR>=2001)
  mei<-mei10[2:122,]
  biotot<-read.delim(file="biomass all months.txt",as.is=T)
  bio<-biotot[,1:3]
  mediamei<-c(0.011,-0.082, -0.185, -0.110,  0.217,  0.180,  0.150,  0.166, -0.027, -0.039,  0.085,  0.018) 
  mediaest<-c(21.13000, 23.20909, 32.20000, 31.40000, 28.85000, 19.75000, 19.40000, 18.30000, 14.65000, 17.45000, 20.80000, 15.63000)
  mediabio<-c(342.6270, 361.7227, 494.0609, 531.3882, 450.0148, 395.8168, 244.5356, 159.5320, 139.1103, 185.0658, 176.2310, 237.1320)
  newdata2=newdata=newbio=numeric()
  
  for ( i in 1:121)
  {
    #mei$MEI[i]<- mei$MEI[i]-mediamei[(i+1)-(trunc(i/12)*12)] 
    newdata2[i]<- estt$estnumbspp[i]-mediaest[(i+1)-(trunc(i/12)*12)]       #esto sirve para quitar a cada valor la media del mes
    newbio[i]<- bio$fruitsb[i]-mediabio[(i+1)-(trunc(i/12)*12)] 
  }
  
  for ( i in 1:(121-12))
  {
    newdata[i]<- estt$estnumbspp[i]-0.3*estt$estnumbspp[i+12] 
  }
  
  
  
  xt<-ccf(mei$MEI, est$estnumbspp)
  xy<-ccf(mei$MEI, bio$fruitsb)
  
  t05<-qt(0.025,119)
  t05/sqrt(119+(t05*t05))
  xt<-ccf(mei$MEI, newdata2,lag.max=30,main="cross correlation number of species and MEI")
  len=length(newdata2)
  newdatalag<-numeric()
  
  for ( i in 1:len-7)
  {
    newdatalag[i]<- newdata2[i+7] 
  }
  
  for ( i in 1:len-8)
  {
    newdatalag[i]<- newdata2[i+8] 
  }
  c<-cor(mei$MEI[1:114],newdatalag[1:114],method="pearson")
  c<-cor.test(mei$MEI[1:113],newdatalag[1:113],method="pearson")
  ccbio<-ccf(mei$MEI, newbio,lag.max=30,main="cross correlation fruit biomass and MEI")
  
  xt12<-ccf(nino12$NINO12, newdata2,lag.max=30,main="cross correlation number of species and NINO 12")
  xt3<-ccf(nino3$NINO3, newdata2,lag.max=30,main="cross correlation number of species and NINO 3")
  xt4<-ccf(nino4$NINO4, newdata2,lag.max=30,main="cross correlation number of species and NINO 4")
  xt34<-ccf(nino34$NINO34, newdata2,lag.max=30,main="cross correlation number of species and NINO 34") 
  xtSOI<-ccf(soindex$SOI, newdata2,lag.max=30,main="cross correlation number of species and SOI") 
  
  rainnrg<-read.delim(file="rain estimation nouragues.txt")
  rain<-rainnrg[2:122,]
  
  newrain<-numeric()
  for ( i in 1:(121-12))
  {
    newrain[i]<- rain[i]-0.3*rain[i+12] 
  }
  
  ccrain<-ccf(newrain, newdata2,lag.max=30,main="cross correlation number of species and rainfall") 
  
}

#this graph plots the mean and sd values of the peak hyperparameters per year (not very interesting)
graphmeanprod=function(nrgresults="nouragues results parameters per year.txt"){
 
  nrgresults=read.delim(nrgresults)
  peakyear=aggregate(data.frame(peak=log(nrgresults$peak+1)), by=list(year=nrgresults$year),mean)
  sdyear=aggregate(data.frame(peak=log(nrgresults$peak+1)), by=list(year=nrgresults$year),sd)
  par(mar=c(2,2,2,2), oma=c(4,3,2,2))
  plot(peakyear$year[1:9],peakyear$peak[1:9], las=1, ylim=c(0,6), main="Nouragues", type="o", bty="l",axes="F")
  #segments(c(0.75,2,3.25,4.5,5.75,7,8.25,9.5,10.75), peakyear$peak - sdyear$peak, c(0.75,2,3.25,4.5,5.75,7,8.25,9.5,10.75), peakyear$peak +sdyear$peak, lwd = 2, lend = 3, col = "black")
  segments(peakyear$year[1:9], peakyear$peak - sdyear$peak, peakyear$year[1:9], peakyear$peak +sdyear$peak, lwd = 2, lend = 3, col = "black")
  axis(side=2)
  axis(side=1, at=peakyear$year[1:9], labels=as.character(peakyear$year[1:9]))
  mtext(side=2, "log of mean number of seeds", line=3)
  mtext(side=1, "years", line=3)
  
  par(mar=c(2,2,2,2), oma=c(4,3,2,2))
  peakdayyear=aggregate(data.frame(peak=nrgresults$peakday), by=list(year=nrgresults$year),mean)
  sdpeakday=aggregate(data.frame(peak=nrgresults$peakday), by=list(year=nrgresults$year),sd)
  #barplot(peakdayyear$peak[1:9], las=1,names.arg=peakdayyear$year[1:9],ylim=c(1,365),space=0.25, main="Nouragues")
  #segments(c(0.75,2,3.25,4.5,5.75,7,8.25,9.5,10.75), peakdayyear$peak - sdpeakday$peak, c(0.75,2,3.25,4.5,5.75,7,8.25,9.5,10.75), peakdayyear$peak +sdpeakday$peak, lwd = 2, lend = 3, col = "black")
  plot(peakdayyear$year[1:9],peakdayyear$peak[1:9], las=1, ylim=c(1,365), main="Nouragues", type="o", bty="l",axes="F")
  segments(peakdayyear$year[1:9], peakdayyear$peak - sdpeakday$peak, peakdayyear$year[1:9], peakdayyear$peak +sdpeakday$peak, lwd = 2, lend = 3, col = "black")
  axis(side=2)
  axis(side=1, at=peakyear$year[1:9], labels=as.character(peakyear$year[1:9]))
  mtext(side=2, "mean day of year", line=3)
  mtext(side=1, "years", line=3)

  #par(new=T)
  #plot(ensoanomalies$MEI[ensoanomalies$YEAR>=2001& ensoanomalies$YEAR<=2009], type="l", axes=F)
  }
  
###this function plots the estimated model of seed production for all the species across time and the mean values
spptime=function(nrgallspp="Nouragues model all spp.txt", k=3){
  
 nrg=read.delim(file=nrgallspp) 
 spnames=sort(unique(nrg$sp))
 par(mar=c(3,3,3,2))
 sphighest=nrg$sp[which.max(nrg$model)]
 dates=sort(unique(nrg$julian))
 N=length(nrg$julian[nrg$sp==sphighest])
 comienzo<-floor((k/2)+1)     
 final<-floor(N-(k/2))+1
 plot(nrg$julian[nrg$sp==sphighest][comienzo:final], runningmean(log(nrg$model[nrg$sp==sphighest]+1),k), bty="l", type="n", las=1,axes=F,las=1, main="Nouragues", ylab="log of biweekly seed fall")
 axis(side=2,las=2)
 xlabels=seq(nrg$julian[nrg$sp==sphighest][1],nrg$julian[nrg$sp==sphighest][length(nrg$julian[nrg$sp==sphighest])]+365,365)
 #axis(side=1, at=xlabels,labels=fromjulian(xlabels))
 axis(side=1, at=xlabels,labels=c("2002","2003","2004","2005","2006","2007","2008","2009","2010","2011"))
 for (i in 1:length(spnames))
 {
   N=length(nrg$julian[nrg$sp==spnames[i]])
   comienzo<-floor((k/2)+1)     
   final<-floor(N-(k/2))+1
   lines(nrg$julian[nrg$sp==spnames[i]][comienzo:final], runningmean(log(nrg$model[nrg$sp==spnames[i]]+1),k),col="grey")
   
 }
 meanjulian=aggregate(data.frame(model=log(nrg$model+1)),by=list(julian=nrg$julian), mean)
 N=length(meanjulian$julian)
 comienzo<-floor((k/2)+1)     
 final<-floor(N-(k/2))+1
 lines( meanjulian$julian[comienzo:final],runningmean(meanjulian$model,k),col="black",lwd=2)
 legend(xlabels[1],6,"mean value per census", lwd=2, bty="n")
 mtext(side=2, "log of biweekly seed production", line=2,cex=1.5)
 mtext(side=1, "time", line=2, cex=1.5)
}


#which species are masting? How can we define a masting species
masting = function(file = "nouragues results parameters per year.txt", fileresults="Nouragues-seed production with CI.pdf")
  
{
  nrg = read.delim(file)
 #I calculate the ratio for the peak of each year in relation to the previous year.
  spnames = sort(unique(nrg$species))
  pdf(file=fileresults)
  par(mar=c(5,4,3,3), las=1)
  for (i in 1:length(spnames))
  {
    
    speciesname=as.character(spnames[i]) 
    sp=nrg[nrg$species==speciesname,]
    ratio=numeric()
    label=character()  
    
     for (k in 1:length(sp)-1)
     {
      ratio[k]=log(sp$peak[k+1]+1)/log(sp$peak[k]+1) 
      label[k]=paste(sp$year[k],"-",sp$year[k+1]) 
      }
    plot(sp$peak, ylim=c(0,max(sp$CIpeak97)), bty="l", xlab="years", ylab="accumulated number of seeds",axes=F)
    axis(side=2)
    axis(side=1, at=1:dim(sp)[1],labels=as.character(sp$year))
    segments(1:dim(sp)[1], sp$peak,1:dim(sp)[1],sp$CIpeak2)
    segments(1:dim(sp)[1], sp$peak,1:dim(sp)[1],sp$CIpeak97)
    mtext(side=3, text=speciesname, font=3,cex=1.5, line=1)
   results=data.frame(species=speciesname,ratio=ratio, labels=label) 
 if (i==1) allresults=results else allresults=rbind(allresults,results)
  }
  dev.off()
  return(allresults)
  }
  #Ratios very different from one are indicating a very high inter-annual variability
  

## How many years with practically 0 seed production was found per species?
CVyears = function(file = "nouragues results parameters per year.txt") {
  
  nrg = read.delim(file)
  #nrg=nrg[-117,]
  spnames=sort(unique(nrg$sp))
  spmeans=aggregate(data.frame(peak=nrg$peak),by=list(year=nrg$year),mean)
  spsd=aggregate(data.frame(peak=nrg$peak),by=list(year=nrg$year),sd)
  condensed=merge(spmeans, spsd, by="year")
  names(condensed)=c("year","mean","sd")
  condensed$CV=condensed$sd/condensed$mean
  hist(condensed$CV, breaks=10,xlab="CV of seed production", ylab="number of species")
  barplot(condensed$CV,xlab="years",ylab="CV of seed production" ,names.arg=condensed$year, ylim=c(0,5))
}

CVpeakdays = function(file="nouragues results parameters per year with cycles.txt") {
  
  nrg=read.delim(file)
  negval=as.vector(which(nrg$peakday<0))
  nrg$peakday[negval]=nrg$peakday[negval]+365
  #lic=which(nrg$species=="Licania membranacea")
  #nrg=nrg[-lic,]
  spnames=sort(unique(nrg$sp))
  spmeans=aggregate(data.frame(peak=nrg$peakday),by=list(species=nrg$sp),mean)
  spsd=aggregate(data.frame(peak=nrg$peakday),by=list(species=nrg$sp),sd)
  condensed=merge(spmeans, spsd, by="species")
  names(condensed)=c("species","mean","sd")
  condensed$CV=condensed$sd/condensed$mean
  tt=hist(condensed$CV, breaks=10,xlab="CV of peakday",right=TRUE)
  orco=order(condensed$CV)
  par(mar=c(14,4,2,2))
  barplot( condensed$CV[orco],xlab="",ylab="CV of day of year" ,names.arg=condensed$year[orco], ylim=c(0,5),las=2)
abline(h=1)
return (condensed)
}

confidencepeakdays=function(file="nouragues results parameters per year.txt", fileresults="Nouragues peakday confidence intervals.pdf")
  
{
  nrg=read.delim(file)
  #I calculate the ratio for the peak of each year in relation to the previous year.
  spnames=sort(unique(nrg$species))
  
  #these are essaies for dealing with negative values of peakday
  #negval=which(nrg$peakday<0)
  #nrg$peakday[negval]=nrg$peakday[negval]+365
  #negvalCI2=which(nrg$CIpeakday2<0)
  #nrg$CIpeakday2[negvalCI2]=nrg$CIpeakday2[negvalCI2]+365
  #negvalCI97=which(nrg$CIpeakday97<0)
  #nrg$CIpeakday97[negvalCI97]=nrg$CIpeakday97[negvalCI97]+365
  
  pdf(file=fileresults)
  par(mar=c(5,4,3,3), las=1)
  for (i in 1:length(spnames))
  {
    
    speciesname=as.character(spnames[i]) 
    sp=nrg[nrg$species==speciesname,]
    plot(sp$peakday, ylim=c(-150,450), bty="l", xlab="years", ylab="day of year",axes=F)
    axis(side=2)
    axis(side=1, at=1:dim(sp)[1],labels=as.character(sp$year))
    segments(1:dim(sp)[1], sp$peakday,1:dim(sp)[1],sp$CIpeakday2)
    segments(1:dim(sp)[1], sp$peakday,1:dim(sp)[1],sp$CIpeakday97)
    mtext(side=3, text=speciesname, font=3,cex=1.5, line=1)
     }
  dev.off()
}


corrbioclim=function(clim=newman, biomass=biomass, est=est){
  
  biomass=biomass[-1,]
  meanclim=aggregate(data.frame(rain=clim$rain, tmin=clim$tmin, tmax= clim$tmax), by=list(month=clim$month), mean)
  meanbio=aggregate(data.frame(fruit=(biomass$fruitsb/(160*0.5)),flower=(biomass$flowersb/(160*0.5))), by=list(month=biomass$month), mean)
  meanest=aggregate(data.frame(est=est$estnumbspp), by=list(month=est$month), mean)
  
  
  cor.test(meanclim$rain, meanest$est)
  cor.test(meanclim$rain, meanbio$flower)
  cor.test(meanclim$rain, meanbio$fruit)
}


##this function removes seasonality from the local climate dataset
#anomalieslocal=nrg.removeseasonality(local=clim)
nrg.removeseasonality=function(local=clim)
{
  season=aggregate(data.frame(rain=local$rain, tmin=local$tmin, tmax=local$tmax), by=list(month=local$mon), mean, na.rm=TRUE)
  
  anomalrain=anomaltmin=anomaltmax=numeric()
  for(k in 1:dim(local)[1]) {
    
    anomalrain[k]=local$rain[k]-season$rain[season$month==local$month[k]]
    anomaltmin[k]=local$tmin[k]-season$tmin[season$month==local$month[k]]
    anomaltmax[k]=local$tmax[k]-season$tmax[season$month==local$month[k]]
    
  }
  localanomalies=data.frame(year=local$year, month=local$month,rain= anomalrain,  tmin=anomaltmin, tmax=anomaltmax)
  return(localanomalies) 
}

####### TABLE 2 ##################

## this function calculates the regression among the mean values of an ENSO index over L months and the sum of estimated parameters of peaks of seed production per year
###ensobio=ENSObiomassreg(biomass=biomass, climseries=ensoanomalies,variable="NAO", monthname="MONTH", yearname="YEAR",Linit=-12, Lfinal=12, factor=2, kval=3)  

ENSObiomassreg=function(biomass=biomass, climseries=ensoanomalies,variable="MEI", monthname="MONTH", yearname="YEAR",kval=3,Linit=-12, Lfinal=12, factor=1)    #L is the number of months used to calculate the mean of ENSO,
  
{
  #setwd( "C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  
  
  biomass$fecha=strptime(paste(biomass$Year2,"-",biomass$month,"-",1), format="%Y - %m - %d")
  
  ind=which(colnames(climseries)==variable)#this selects the climatic variable of study
  var=climseries[,ind]
  ind2=which(colnames(climseries)==monthname)
  month=climseries[,ind2]
  ind3=which(colnames(climseries)==yearname)
  year=climseries[,ind3]
  julian=climseries$julian
  clim2=data.frame(julian,year,month, var)
  
  franomal=flanomal=numeric(length=length(biomass$fruitsb))
  
  for (k in 1:length(biomass$fruitsb))
  { 
    meanbio=aggregate(data.frame(fr=biomass$fruitsb,fl=biomass$flowersb), by=list(month=biomass$month), mean)
    franomal[k]=biomass$fruitsb[k]-meanbio$fr[meanbio$month==biomass$month[k]]
    flanomal[k]=biomass$flowersb[k]-meanbio$fl[meanbio$month==biomass$month[k]]
  }
  biomass$franomal=franomal; biomass$flanomal=flanomal
  
  k=kval
  start=which(biomass$Year2==2001 & biomass$month==2)
  end=which(biomass$Year2==2012 & biomass$month==1)
  N=length(biomass$fruitsb[start:end])
  biomass=biomass[start:end,]
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))
  
  par(mar=c(4,4,2,4),cex=1.75)
  plot(biomass$franomal,type="n", xlab="time", ylab="fruit biomass anomalies (g/m?)", bty="l", axes=F, las=1, ylim=c(-2,1.5)*factor)
  axis(side=2, las=1,col="blue3")
  jan1=which(biomass$month==1)
  axis(side=1, at=c(jan1[1]-12,jan1), labels=as.character(c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")))
  lines(c(comienzo:(final+1)),runningmean(biomass$franomal,k), col="blue3", lwd=2)
  legend(end-130,-3, lty=c(1,2), c("fruit biomass",variable), bty="n",col=c("blue3", "red"), horiz=T,lwd=2)
  #lines(biomass$fecha[start:(end-2)],runningmean(biomass$flowersb[start:end],k)/(160*0.5), lty=2)
  
  startmei=which(clim2$year==2001 & clim2$month==2)
  endmei=which(clim2$year==2012 & clim2$month==1)
  par(new=T)
  plot(c(start:(end-2)),runningmean(clim2$var[startmei:endmei]),type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2, ylim=c(-2,1.5))
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, variable, line=3, cex=1.5)
  
  #plot(biomass$fecha[2:length(biomass$fecha)-2], las=1,runningmean(biomass$franomal, k=3), type="l", col="red", xlab="time", ylab="anomalies", bty="l")
  #lines(biomass$fecha[2:length(biomass$fecha)-2], runningmean(biomass$flanomal, k=3), col="blue")
  #abline(h=0)
  #legend(biomass$fecha[53],2.5, lty=c(1,1), col=c("red", "blue"), c("fruit biomass", "flower biomass"), bty="n")
  
  xval=matrix(nrow=length (biomass$fruitsb),ncol=length(c(Linit:Lfinal)))
  #colnames(xvalmean, xvalmax, yval)=as.character(c(Linit:Lfinal))
  yval=yflval=betafr=betafl=rfr=rfl=significancefr =significancefl= lag=numeric()
  
  for (l in 1:length(c(Linit:Lfinal)))
  {
    cat("L:",c(Linit:Lfinal)[l], "\n" )
    
    for(o in 1: length (biomass$fruitsb))
    {
      ###f?jate aqu?
      final=(which(clim2$year==biomass$Year2[o]&clim2$month==biomass$month[o]))-c(Linit:Lfinal)[l]
      xval[o,l]=clim2$var[final]
      yval[o]=biomass$franomal[o]
      yflval[o]=biomass$flanomal[o]
    } 
    betafr[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,1]
    rfr[l]=sqrt(summary.lm(lm(yval~xval[,l]))$adj.r.squared)
    betafl[l]=summary.lm(lm(yflval~xval[,l]))$coefficients[2,1]
    rfl[l]=sqrt(summary.lm(lm(yflval~xval[,l]))$adj.r.squared)
    significancefr[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,4]
    significancefl[l]=summary.lm(lm(yflval~xval[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  } 
  results=list()
  results$xvalues=xval
  results$yvalues=yval
  results$yvaluesflo=yflval
  results$sig=data.frame(variable,lag, betafr,significancefr, rfr,betafl, significancefl,rfl)
  return(results)
  
}
###
reglag4=function(biomass=biomass, climseries=ensoanomalies,variable="MEI", monthname="MONTH", yearname="YEAR",kval=3,Linit=-12, Lfinal=12, factor=1)    #L is the number of months used to calculate the mean of ENSO,
  
{
  #setwd( "C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  
  
  biomass$fecha=strptime(paste(biomass$Year2,"-",biomass$month,"-",1), format="%Y - %m - %d")
  
  clim2=data.frame(julian=climseries$julian,year=climseries$YEAR,month=climseries$MONTH, var=climseries$MEI)
  
  franomal=flanomal=numeric(length=length(biomass$fruitsb))
  
  meanbio=aggregate(data.frame(fr=biomass$fruitsb,fl=biomass$flowersb), by=list(month=biomass$month), mean)
  for (k in 1:length(biomass$fruitsb))
  { 
    franomal[k]=biomass$fruitsb[k]-meanbio$fr[meanbio$month==biomass$month[k]]
    flanomal[k]=biomass$flowersb[k]-meanbio$fl[meanbio$month==biomass$month[k]]
  }
  biomass$franomal=franomal; biomass$flanomal=flanomal
  
  k=kval
  start=which(biomass$Year2==2001 & biomass$month==2)
  end=which(biomass$Year2==2012 & biomass$month==1)
  N=length(biomass$fruitsb[start:end])
  biomass=biomass[start:end,]
  comienzo<-as.integer((k/2)+1)     
  final<-as.integer(N-(k/2))
  
  par(mar=c(4,4,2,4),cex=1.75)
  plot(biomass$franomal,type="n", xlab="time", ylab="fruit biomass anomalies (g/m?)", bty="l", axes=F, las=1, ylim=c(-2,1.5)*factor)
  axis(side=2, las=1,col="blue3")
  jan1=which(biomass$month==1)
  axis(side=1, at=c(jan1[1]-12,jan1), labels=as.character(c("2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")))
  lines(c(comienzo:(final+1)),runningmean(biomass$franomal,k), col="blue3", lwd=2)
  legend(end-130,-3, lty=c(1,2), c("fruit biomass",variable), bty="n",col=c("blue3", "red"), horiz=T,lwd=2)
  #lines(biomass$fecha[start:(end-2)],runningmean(biomass$flowersb[start:end],k)/(160*0.5), lty=2)
  
  startmei=which(clim2$year==2001 & clim2$month==2)
  endmei=which(clim2$year==2012 & clim2$month==1)
  par(new=T)
  plot(c(start:(end-2)),runningmean(clim2$var[startmei:endmei]),type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2, ylim=c(-2,1.5))
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, variable, line=3, cex=1.5)
  
  #plot(biomass$fecha[2:length(biomass$fecha)-2], las=1,runningmean(biomass$franomal, k=3), type="l", col="red", xlab="time", ylab="anomalies", bty="l")
  #lines(biomass$fecha[2:length(biomass$fecha)-2], runningmean(biomass$flanomal, k=3), col="blue")
  #abline(h=0)
  #legend(biomass$fecha[53],2.5, lty=c(1,1), col=c("red", "blue"), c("fruit biomass", "flower biomass"), bty="n")
  
  xval=matrix(nrow=length (biomass$fruitsb),ncol=length(c(Linit:Lfinal)))
  #colnames(xvalmean, xvalmax, yval)=as.character(c(Linit:Lfinal))
  yval=yflval=betafr=betafl=rfr=rfl=significancefr =significancefl= lag=numeric()
  
  for (l in 1:length(c(Linit:Lfinal)))
  {
    cat("L:",c(Linit:Lfinal)[l], "\n" )
    
    for(o in 1: length (biomass$fruitsb))
    {
      ###f?jate aqu?
      final=(which(clim2$year==biomass$Year2[o]&clim2$month==biomass$month[o]))-c(Linit:Lfinal)[l]
      xval[o,l]=clim2$var[final]
      yval[o]=biomass$franomal[o]
      yflval[o]=biomass$flanomal[o]
    } 
    betafr[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,1]
    rfr[l]=sqrt(summary.lm(lm(yval~xval[,l]))$adj.r.squared)
    betafl[l]=summary.lm(lm(yflval~xval[,l]))$coefficients[2,1]
    rfl[l]=sqrt(summary.lm(lm(yflval~xval[,l]))$adj.r.squared)
    significancefr[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,4]
    significancefl[l]=summary.lm(lm(yflval~xval[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  } 
  results=list()
  results$xvalues=xval
  results$yvalues=yval
  results$yvaluesflo=yflval
  results$sig=data.frame(variable,lag, betafr,significancefr, rfr,betafl, significancefl,rfl)
  return(results)
  
}

#ensoest=ENSOestreg(estfile="number total spp per month estimated.txt", climseries=ensoanomalies,variable="MEI", monthname="MONTH", yearname="YEAR",Linit=-12, Lfinal=12, kval=3)   
ENSOestreg=function(estfile="number total spp per month estimated.txt", climseries=ensoanomalies,variable="MEI", monthname="MONTH", yearname="YEAR",Linit=-12, Lfinal=12, kval=3)    #L is the number of months used to calculate the mean of ENSO,
  
{
  #setwd( "C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  est=read.delim(estfile)
  est$fecha=strptime(paste(est$year,"-",est$month,"-",1), format="%Y - %m - %d")
  
  ind=which(colnames(climseries)==variable)#this selects the climatic variable of study
  var=climseries[,ind]
  ind2=which(colnames(climseries)==monthname)
  month=climseries[,ind2]
  ind3=which(colnames(climseries)==yearname)
  year=climseries[,ind3]
  julian=climseries$julian
  clim2=data.frame(julian,year,month, var)
  
  anomal=numeric(length=length(est$estnumbspp))
  
  for (k in 1:length(est$estnumbspp))
  { 
    meanest=aggregate(data.frame(est=est$estnumbspp), by=list(month=est$month), mean)
    anomal[k]=est$estnumbspp[k]-meanest$es[meanest$month==est$month[k]]
  }
  est$anomal=anomal
  
  N=length(clim2$var[startmei:endmei])
  comienzo<-as.integer((kval/2))     
  final<-as.integer(N-(kval/2))
  par(mar=c(4,4,2,4),cex=1.75)
  plot(est$fecha[comienzo:final], las=1,runningmean(est$anomal, k=kval), type="l", col="blue", xlab="time", ylab="anomalies", bty="l", ylim=c(-20,15))
  abline(h=0)
  startmei=which(clim2$year==2001 & clim2$month==2)
  endmei=which(clim2$year==2011& clim2$month==2)
  
 
  
  par(new=T)
  plot(c(comienzo:final),runningmean(clim2$var[startmei:endmei],k=kval),type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2, ylim=c(-2,1.5))
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, variable, line=3, cex=1.5)
  
  xval=matrix(nrow=length (est$anomal),ncol=length(c(Linit:Lfinal)))
  #colnames(xvalmean, xvalmax, yval)=as.character(c(Linit:Lfinal))
  yval=beta=r=significance= lag=numeric()
  
  for (l in 1:length(c(Linit:Lfinal)))
  {
    cat("L:",c(Linit:Lfinal)[l], "\n" )
    
    for(o in 1: length (est$anomal))
    {
      final=(which(clim2$year==est$year[o]&clim2$month==est$month[o]))-c(Linit:Lfinal)[l]
      xval[o,l]=clim2$var[final]
      yval[o]=est$anomal[o]
    } 
    beta[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,1]
    r[l]=sqrt(summary.lm(lm(yval~xval[,l]))$adj.r.squared)
    significance[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  } 
  results=list()
  results$xvalues=xval
  results$yvalues=yval
  results$sig=data.frame(variable,lag, beta,significance, r)
  return(results)
  
}

#biolocal=climbiomassreg(bio=biomass, climseries=clim,variable="rain", monthname="month", yearname="year",Linit=-6, Lfinal=6) 
## this function calculates the regression among the mean values of an local climate variable over L months and the sum of estimated parameters of peaks of seed production per year
climbiomassreg=function(bio=biomass, climseries=clim,variable="rain", monthname="month", yearname="year",Linit=-6, Lfinal=6)    #L is the number of months used to calculate the mean of ENSO,
  
{
  #setwd( "C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  
  datayr=bio
  biomass$fecha=strptime(paste(biomass$Year2,"-",biomass$month,"-",1), format="%Y - %m - %d")
  biostart=which(biomass$Year2==2001 & biomass$month==7)
  bioend=which(biomass$Year2==2011 & biomass$month==1)
  biomass=biomass[biostart:bioend,]
  
  ind=which(colnames(climseries)==variable)#this selects the climatic variable of study
  var=climseries[,ind]
  ind2=which(colnames(climseries)==monthname)
  month=climseries[,ind2]
  ind3=which(colnames(climseries)==yearname)
  year=climseries[,ind3]
  julian=climseries$julian
  clim2=data.frame(julian,year,month, var)
  
  anomalclim<-numeric(length=length(clim2$var))
  
  for (z in 1:length(var))
  { 
    meanclim=aggregate(data.frame(index=clim2$var), by=list(month=clim2$month), mean)
    anomalclim[z]=clim2$var[z]-meanclim$index[meanclim$month==clim2$month[z]]
  }
  clim2$anomal=anomalclim
  plot(clim2$julian[2:length(clim2$julian)-2], las=1,runningmean(clim2$anomal, k=3), type="l", xlab="time", ylab="anomalies of monthly values of rainfall", bty="l")
  abline(h=0)
  
  
  franomal=flanomal=numeric(length=length(biomass$fruitsb))
  
  for (k in 1:length(biomass$fruitsb))
  { 
    meanbio=aggregate(data.frame(fr=biomass$fruitsb,fl=biomass$flowersb), by=list(month=biomass$month), mean)
    franomal[k]=biomass$fruitsb[k]-meanbio$fr[meanbio$month==biomass$month[k]]
    flanomal[k]=biomass$flowersb[k]-meanbio$fl[meanbio$month==biomass$month[k]]
  }
  biomass$franomal=franomal; biomass$flanomal=flanomal
  plot(biomass$fecha[2:length(biomass$fecha)-2], las=1,runningmean(biomass$franomal, k=3), type="l", col="red", xlab="time", ylab="anomalies", bty="l")
  lines(biomass$fecha[2:length(biomass$fecha)-2], runningmean(biomass$flanomal, k=3), col="blue")
  abline(h=0)
  legend(biomass$fecha[53],2.5, lty=c(1,1), col=c("red", "blue"), c("fruit biomass", "flower biomass"), bty="n")
  
  xval=matrix(nrow=length (biomass$fruitsb),ncol=length(c(Linit:Lfinal)))
  #colnames(xvalmean, xvalmax, yval)=as.character(c(Linit:Lfinal))
  yval=yflval=betafr=betafl=rfr=rfl=significancefr =significancefl= lag=numeric()
  
  for (l in 1:length(c(Linit:Lfinal)))
  {
    cat("L:",c(Linit:Lfinal)[l], "\n" )
    
    for(o in 1: length (biomass$fruitsb))
    {
      final=(which(clim2$year==biomass$Year2[o]&clim2$month==biomass$month[o]))-c(Linit:Lfinal)[l]
      xval[o,l]=clim2$anomal[final]
      yval[o]=biomass$franomal[o]
      yflval[o]=biomass$flanomal[o]
    } 
    betafr[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,1]
    rfr[l]=sqrt(summary.lm(lm(yval~xval[,l]))$adj.r.squared)
    betafl[l]=summary.lm(lm(yflval~xval[,l]))$coefficients[2,1]
    rfl[l]=sqrt(summary.lm(lm(yflval~xval[,l]))$adj.r.squared)
    significancefr[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,4]
    significancefl[l]=summary.lm(lm(yflval~xval[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  } 
  results=list()
  results$xvalues=xval
  results$yvalues=yval
  results$yvaluesflo=yflval
  results$sig=data.frame(variable,lag, betafr,significancefr, rfr,betafl, significancefl,rfl)
  return(results)
  
}


anomalieslocal=nrg.removeseasonality(local=clim)
#local=anomalieslocal; ensoanomalies=ensoanomalies;varlocal="rain";variable="MEI";k=3;Linit=-6; Lfinal=6; startyear=2001; endyear=2011; startmonth=6;endmonth=6; factor=50; unity="mm"  
#el=ensovslocal(local=anomalieslocal, ensoanomalies=ensoanomalies,varlocal="rain", variable="MEI",k=3,Linit=-6, Lfinal=6, startyear=2001, endyear=2011, startmonth=6, endmonth=6, factor=50, unity="mm")
#el=ensovslocal(local=anomalieslocal, ensoanomalies=ensoanomalies,varlocal="tmin", variable="MEI",k=3,Linit=-6, Lfinal=6, startyear=2001, endyear=2011, startmonth=6, endmonth=1, factor=0.25, unity="?C")
ensovslocal=function(local=anomalieslocal, ensoanomalies=ensoanomalies,varlocal="rain", variable="MEI",k=3,Linit=-6, Lfinal=6, startyear=2001, endyear=2011, startmonth=6, endmonth=6, factor=50, unity="mm") {
  
  iniciolocal=which(local$year==startyear&local$mon==startmonth)
  finlocal=which(local$year==endyear&local$mon==endmonth)
  local=local[iniciolocal:finlocal,]
  N=length(local$rain)
  comienzo<-floor((k/2)+1)     
  final<-floor(N-(k/2))+1
  
  ind2=which(colnames(local)==varlocal)#this selects the climatic variable of study
  var=local[,ind2]
  month=local$month
  year=local$year
  local2=data.frame(year, month, var)
  
  par(mar=c(4,4,2,4),cex=1.25)
  plot(local2$var,type="n", xlab="time", ylab=paste(varlocal,"anomalies","(",unity,")"),ylim=c(-factor*2,factor*3),bty="l", axes=F, las=1)
  #plot(local2$var,type="n", xlab="time", ylab=varlocal,bty="l", axes=F, las=1)
  axis(side=2, las=1)
  jan1=which(local2$month==1)
  axis(side=1, at=c((jan1[1]-12),jan1), labels=as.character(unique(local2$year)))
  lines(c(comienzo:final),runningmean(local2$var,k), col="blue", lwd=2)
  legend(length(comienzo:final)/2,factor*3, lty=c(1,2),lwd=2, c(varlocal, variable), col=c( "blue", "red"), bty="n", horiz=T)
  
  
  ind=which(colnames(ensoanomalies)==variable)#this selects the climatic variable of study
  var=ensoanomalies[,ind]
  year=ensoanomalies$YEAR
  month=ensoanomalies$MONTH
  clim2=data.frame(year,month, var)
  
  startmei=which(clim2$year==startyear & clim2$month==1)
  endmei=which(clim2$year==endyear & clim2$month==endmonth)
  par(new=T)
  plot(clim2$var[startmei:endmei]*factor,type="l", lty=2,ylab="",xlab="",axes=F, las=1,  lwd=2,col="red",ylim=c(-2,3)*factor)
  abline(h=0, lty=2)
  axis(side=4, las=1, at=c(-2:3)*factor, labels=c(-2:3))
  mtext(side=4, variable, line=2,cex=1.25)
  
  xval=matrix(nrow=length (local2$var),ncol=length(c(Linit:Lfinal)))
  #colnames(xvalmean, xvalmax, yval)=as.character(c(Linit:Lfinal))
  yval=beta=r=significance=lag=numeric()
  
  for (l in 1:length(c(Linit:Lfinal)))
  {
    cat("L:",c(Linit:Lfinal)[l], "\n" )
    
    for(o in 1: length (local2$var))
    {
      final=(which(clim2$year==local2$year[o]&clim2$month==local2$month[o]))-c(Linit:Lfinal)[l]
      xval[o,l]=clim2$var[final]
      yval[o]=local2$var[o]
      
    } 
    beta[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,1]
    r[l]=sqrt(summary.lm(lm(yval~xval[,l]))$adj.r.squared)
    significance[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  } 
  results=list()
  results$xvalues=xval
  results$yvalues=yval
  results$sig=data.frame(variable,lag, beta,significance, r)
  return(results)
  
  
}


####REGRESSION BETWEEN THE MODEL OF SEED PRODUCTION AND MEI####


ENSOmodreg=function(infile="nouragues results parameters per year.txt", soifile="ENSO anomalies - 2012.txt",index, Linit=1, Lfinal=12)    #L is the number of months used to calculate the mean of ENSO,
  
{
  #setwd( "C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  
  datayr=read.delim(file=infile)
  day=jmonth=newyear=numeric()
  
  for (i in 1:length(datayr$peakday))         #this part of the function calculates the month of the year corresponding to each peakdays per year
  { 
    day[i]=tojulian(pst('01/01/',datayr$year[i]))+datayr$peakday[i]
    fulldate=create.fulldate(fromjulian(day[i],dateform='%Y-%m-%d'),format='%Y-%m-%d')
    jmonth[i]=fulldate$month
    newyear[i]=fulldate$year
  } 
  datayr$month=jmonth
  datayr$nyr=newyear
  
  soi<-read.delim(soifile,as.is=T)
  #gobl=which(datayr$sp=="Gouania blanchetiana")
  #datayr=datayr[-gobl,]
  sumseed=aggregate(data.frame(peak=datayr$peak), by=list(month=datayr$month,year=datayr$nyr), sum)
  sumseed=subset(sumseed, sumseed$year>2000 & sumseed$year < 2010)  
  #nbcens=aggregate(data.frame(nbcens=datayr$julian),by=list(month=datayr$month,year=datayr$year),lengthunique )
  #sumseed$pkstd=sumseed$peak/nbcens$nbcens
  sumseed$julian=tojulian(paste("1/",sumseed$month,"/",sumseed$year), dateform="%d/ %m / %Y")
  #plot(sumseed$julian[2:(length(sumseed$julian)-2)], runningmean(sumseed$peak, k=4), type="l", xlab="dates", ylab="number of species with a peak of production") 
  
  nbspp=aggregate(datayr$sp, by=list(month=datayr$month, year=datayr$year), lengthunique)
  nbspp$julian=tojulian(paste("1/",nbspp$month,"/",nbspp$year), dateform="%d/ %m / %Y")
  #plot(nbspp$julian[2:(length(nbspp$julian)-1)], runningmean(nbspp$x, k=3), type="l", xlab="dates", ylab="number of species with a peak of production")
  indcol=which(colnames(soi)==index)
  enso<-soi[,c(1,2,indcol)]
  names(enso)=c("YEAR","MONTH","index")
  meanenso= meanmei(enso)
  
  anomal<-numeric(length=length(sumseed$peak))
  
  for (k in 1:length(anomal))
  { 
    meanpeak=aggregate(data.frame(peak=sumseed$peak), by=list(month=sumseed$month), mean)
    anomal[k]=sumseed$peak[k]-meanpeak$peak[meanpeak$month==sumseed$month[k]]
  }
  sumseed=data.frame(sumseed, anomal)
  
  #par(mar=c(4,2,2,1), oma=c(2,2,2,1))
  #plot(sumseed$julian, sumseed$anomal, xlab="date", ylab="anomalies of seed production",las=1, axes=F, type="l")
  #axis(side=1, at=sumseed$julian[c(jan1)], labels=as.character(sumseed$julian[jan1]), colnames=)
  #axis(side=2, las=1)
  #abline(h=0)
  
  xvalmean=xvalmax=xvalnormal=yval=matrix(nrow=length (sumseed$peak),ncol=length(c(Linit:Lfinal)))
  #colnames(xvalmean, xvalmax, yval)=as.character(c(Linit:Lfinal))
  betamax=betamean=betanorm=significancemean =significancemax=significancenorm= lag=numeric()
  
  for (l in 1:length(c(Linit:Lfinal)))
  {
    cat("L:",c(Linit:Lfinal)[l], "\n" )
    
    for(o in 1: length (sumseed$peak))
    {
      final=(which(enso$YEAR==sumseed$year[o]&enso$MONTH==sumseed$month[o]))-Linit
      inicio=final-c(Linit:Lfinal)[l]
      xvalmean[o,l]=mean(enso$index[inicio:final])
      xvalmax[o,l]=max(enso$index[inicio:final])
      yval[o,l]=sumseed$peak[o]
      final2=which(enso$YEAR==sumseed$year[o]&enso$MONTH==sumseed$month[o])-c(Linit:Lfinal)[l]
      xvalnormal[o,l]=enso$index[final2]
      
    } 
    betamean[l]=summary.lm(lm(log(yval[,l])~xvalmean[,l]))$coefficients[2,1]
    betamax[l]=summary.lm(lm(log(yval[,l])~xvalmax[,l]))$coefficients[2,1]
    betanorm[l]=summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$coefficients[2,1]
    significancemean[l]=summary.lm(lm(log(yval[,l])~xvalmean[,l]))$coefficients[2,4]
    significancemax[l]=summary.lm(lm(log(yval[,l])~xvalmax[,l]))$coefficients[2,4]
    significancenorm[l]=summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  } 
  results=list()
  results$xvaluesmean=xvalmean
  results$xvaluesmax=xvalmax
  results$yvalues=yval
  results$sig=data.frame(lag, betamean,betamax, betanorm,significancemean, significancemax, significancenorm)
  results$sig$index=index
  return(results)
  
}

### cross-correlation among the model of seed production and the ENSO indices

pnormmod=function(infile=NRGallspp, soifile="ENSO anomalies - 2012.txt",index, Linit=1, Lfinal=12)    #L is the number of months used to calculate the mean of ENSO,
  
{
  #setwd( "C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  
  datayr=infile
  day=jmonth=newyear=numeric()
  
  for (i in 1:length(datayr$peakday))         #this part of the function calculates the month of the year corresponding to each peakdays per year
  { 
    day[i]=tojulian(pst('01/01/',datayr$year[i]))+datayr$peakday[i]
    fulldate=create.fulldate(fromjulian(day[i],dateform='%Y-%m-%d'),format='%Y-%m-%d')
    jmonth[i]=fulldate$month
    newyear[i]=fulldate$year
  } 
  datayr$month=jmonth
  datayr$nyr=newyear
  
  soi<-read.delim(soifile,as.is=T)
  #gobl=which(datayr$sp=="Gouania blanchetiana")
  #datayr=datayr[-gobl,]
  sumseed=aggregate(data.frame(peak=datayr$peak), by=list(month=datayr$month,year=datayr$nyr), sum)
  sumseed=subset(sumseed, sumseed$year>2000 & sumseed$year < 2010)  
  #nbcens=aggregate(data.frame(nbcens=datayr$julian),by=list(month=datayr$month,year=datayr$year),lengthunique )
  #sumseed$pkstd=sumseed$peak/nbcens$nbcens
  sumseed$julian=tojulian(paste("1/",sumseed$month,"/",sumseed$year), dateform="%d/ %m / %Y")
  #plot(sumseed$julian[2:(length(sumseed$julian)-2)], runningmean(sumseed$peak, k=4), type="l", xlab="dates", ylab="number of species with a peak of production") 
  
  nbspp=aggregate(datayr$sp, by=list(month=datayr$month, year=datayr$year), lengthunique)
  nbspp$julian=tojulian(paste("1/",nbspp$month,"/",nbspp$year), dateform="%d/ %m / %Y")
  #plot(nbspp$julian[2:(length(nbspp$julian)-1)], runningmean(nbspp$x, k=3), type="l", xlab="dates", ylab="number of species with a peak of production")
  indcol=which(colnames(soi)==index)
  enso<-soi[,c(1,2,indcol)]
  names(enso)=c("YEAR","MONTH","index")
  meanenso= meanmei(enso)
  
  anomal<-numeric(length=length(sumseed$peak))
  
  for (k in 1:length(anomal))
  { 
    meanpeak=aggregate(data.frame(peak=sumseed$peak), by=list(month=sumseed$month), mean)
    anomal[k]=sumseed$peak[k]-meanpeak$peak[meanpeak$month==sumseed$month[k]]
  }
  sumseed=data.frame(sumseed, anomal)
  
  #par(mar=c(4,2,2,1), oma=c(2,2,2,1))
  #plot(sumseed$julian, sumseed$anomal, xlab="date", ylab="anomalies of seed production",las=1, axes=F, type="l")
  #axis(side=1, at=sumseed$julian[c(jan1)], labels=as.character(sumseed$julian[jan1]), colnames=)
  #axis(side=2, las=1)
  #abline(h=0)
  
  xvalmean=xvalmax=xvalnormal=yval=matrix(nrow=length (sumseed$peak),ncol=length(c(Linit:Lfinal)))
  #colnames(xvalmean, xvalmax, yval)=as.character(c(Linit:Lfinal))
  betamax=betamean=betanorm=significancemean =significancemax=significancenorm= lag=numeric()
  
  for (l in 1:length(c(Linit:Lfinal)))
  {
    cat("L:",c(Linit:Lfinal)[l], "\n" )
    
    for(o in 1: length (sumseed$peak))
    {
      final=(which(enso$YEAR==sumseed$year[o]&enso$MONTH==sumseed$month[o]))-Linit
      inicio=final-c(Linit:Lfinal)[l]
      xvalmean[o,l]=mean(enso$index[inicio:final])
      xvalmax[o,l]=max(enso$index[inicio:final])
      yval[o,l]=sumseed$peak[o]
      final2=which(enso$YEAR==sumseed$year[o]&enso$MONTH==sumseed$month[o])-c(Linit:Lfinal)[l]
      xvalnormal[o,l]=enso$index[final2]
      
    } 
    betamean[l]=summary.lm(lm(log(yval[,l])~xvalmean[,l]))$coefficients[2,1]
    betamax[l]=summary.lm(lm(log(yval[,l])~xvalmax[,l]))$coefficients[2,1]
    betanorm[l]=summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$coefficients[2,1]
    significancemean[l]=summary.lm(lm(log(yval[,l])~xvalmean[,l]))$coefficients[2,4]
    significancemax[l]=summary.lm(lm(log(yval[,l])~xvalmax[,l]))$coefficients[2,4]
    significancenorm[l]=summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  } 
  results=list()
  results$xvaluesmean=xvalmean
  results$xvaluesmax=xvalmax
  results$yvalues=yval
  results$sig=data.frame(lag, betamean,betamax, betanorm,significancemean, significancemax, significancenorm)
  results$sig$index=index
  return(results)
  
}





####REGRESSION BETWEEN THE MONTHLY VALUES WITH THE MODEL OF SEED PRODUCTION AND MEI####


## I am using for this regression the monthly values of the seed model

#file="Nouragues model all spp.txt";fileno="NRG model all spp without Dj.txt"; beginyearfile=beginyearfile;spstart=1;spend=45; soifile="local climate data Nouragues.txt";index="rain"; Linit=1; Lfinal=12; local=TRUE

##ddrain=ENSOmonthlymodreg(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", beginyearfile=beginyearfile,spstart=1,spend=45, soifile="local climate data Nouragues.txt",index="rain", Linit=1, Lfinal=12, local=TRUE,kval=3)
#dd=ENSOmonthlymodreg(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", beginyearfile=beginyearfile,spstart=1,spend=45, soifile="ENSO anomalies - 2012.txt",index="MEI", Linit=1, Lfinal=12, local=FALSE,kval=3)
ENSOmonthlymodreg=function(file="Nouragues model all spp.txt",fileno="NRG model all spp without Dj.txt", beginyearfile=beginyearfile,spstart=1,spend=45, soifile="ENSO anomalies - 2012.txt",index="MEI",Linit=1, Lfinal=12, local=FALSE,kval=3)    #L is the number of months used to calculate the mean of ENSO,
  
{
  #setwd( "C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  
  datayr=monthlyvalues(file=file,fileno=fileno, beginyearfile=beginyearfile,spstart=spstart,spend=spend)
  soi<-read.delim(soifile,as.is=T)
  #gobl=which(datayr$sp=="Gouania blanchetiana")

  sumseed=aggregate(data.frame(peak=datayr$model), by=list(month=datayr$month,year=datayr$year), sum)
  sumseed=sumseed[sumseed$year>2002 & sumseed$year < 2010,]  
  #nbcens=aggregate(data.frame(nbcens=datayr$julian),by=list(month=datayr$month,year=datayr$year),lengthunique )
  #sumseed$pkstd=sumseed$peak/nbcens$nbcens
  sumseed$julian=tojulian(paste("15/",sumseed$month,"/",sumseed$year), dateform="%d/ %m / %Y")
  #plot(sumseed$julian[2:(length(sumseed$julian)-2)], runningmean(sumseed$peak, k=4), type="l", xlab="dates", ylab="number of species with a peak of production") 
  
  nbspp=aggregate(datayr$sp, by=list(month=datayr$month, year=datayr$year), lengthunique)
  nbspp$julian=tojulian(paste("1/",nbspp$month,"/",nbspp$year), dateform="%d/ %m / %Y")
  #plot(nbspp$julian[2:(length(nbspp$julian)-1)], runningmean(nbspp$x, k=3), type="l", xlab="dates", ylab="number of species with a peak of production")
  
  
  if (local==TRUE){
    season=aggregate(data.frame(rain=soi$rain, tmin=soi$tmin, tmax=soi$tmax), by=list(month=soi$mon), mean, na.rm=TRUE)
    anomalrain=anomaltmin=anomaltmax=numeric()
    for(m in 1:dim(soi)[1]) {
      
      anomalrain[m]=soi$rain[m]-season$rain[season$month==soi$month[m]]
      anomaltmin[m]=soi$tmin[m]-season$tmin[season$month==soi$month[m]]
      anomaltmax[m]=soi$tmax[m]-season$tmax[season$month==soi$month[m]]
      
    }
    soi=data.frame(year=soi$year, month=soi$month,rain= anomalrain, tmin=anomaltmin, tmax=anomaltmax) 
    }
  
  indcol=which(colnames(soi)==index)
  enso<-soi[,c(1,2,indcol)]
  names(enso)=c("YEAR","MONTH","index")
  anomal<-numeric(length=length(sumseed$peak))
  
  for (k in 1:length(anomal))
  { 
    meanpeak=aggregate(data.frame(peak=sumseed$peak), by=list(month=sumseed$month), mean)
    anomal[k]=sumseed$peak[k]-meanpeak$peak[meanpeak$month==sumseed$month[k]]
  }
  sumseed=data.frame(sumseed, anomal)
  
  #par(mar=c(4,2,2,1), oma=c(2,2,2,1))
  #plot(sumseed$julian, sumseed$anomal, xlab="date", ylab="anomalies of seed production",las=1, axes=F, type="l")
  #axis(side=1, at=sumseed$julian[c(jan1)], labels=as.character(sumseed$julian[jan1]), colnames=)
  #axis(side=2, las=1)
  #abline(h=0)
  par(mar=c(4,4,2,4),cex=1.75)
  startmei=which(enso$YEAR==2003 & enso$MONTH==1)
  endmei=which(enso$YEAR==2009& enso$MONTH==12)
  N=length(enso$index[startmei:endmei])
  comienzo<-as.integer((kval/2))     
  final<-as.integer(N-(kval/2))
 
  plot(sumseed$julian[comienzo:final], las=1,runningmean(sumseed$anomal, k=kval), type="l", col="blue", xlab="time", ylab="anomalies of seed counts", bty="l", ylim=c(-1400,1000))
  abline(h=0)
  
  par(new=T)
  plot(c(comienzo:final),runningmean(enso$index[startmei:endmei],k=kval)*667,type="l", lty=2,ylab="",xlab="",axes=F, las=1, col="red", lwd=2, ylim=c(-2,1.5)*667)
  abline(h=0, lty=2)
  axis(side=4, las=1, col="red")
  mtext(side=4, index, line=3, cex=1.5)
  
  
  xvalmean=xvalmax=xvalnormal=yval=matrix(nrow=length (sumseed$peak),ncol=length(c(Linit:Lfinal)))
  #colnames(xvalmean, xvalmax, yval)=as.character(c(Linit:Lfinal))
  betamax=betamean=betanorm=significancemean =significancemax=significancenorm= rmean=rmax= rnorm=lag=numeric()
  
  for (l in 1:length(c(Linit:Lfinal)))
  {
    cat("L:",c(Linit:Lfinal)[l], "\n" )
    
    for(o in 1: length (sumseed$peak))
    {
      #l?nea 358 hierarchical-enso-local climate
      #final=(which(enso$YEAR==sumseed$year[o]&enso$MONTH==sumseed$month[o]))-c(Linit:Lfinal)[l]
      #yval[o,l]=sumseed$peak[o]inicio=final-c(Linit:Lfinal)[l]
      final=(which(enso$YEAR==sumseed$year[o]&enso$MONTH==sumseed$month[o]))-c(Linit:Lfinal)[l]
      inicio=which(enso$YEAR==sumseed$year[o]&enso$MONTH==sumseed$month[o])
      xvalmean[o,l]=mean(enso$index[inicio:final])
      xvalmax[o,l]=max(abs(enso$index[inicio:final]))
      yval[o,l]=sumseed$peak[o]
      final2=which(enso$YEAR==sumseed$year[o]&enso$MONTH==sumseed$month[o])-c(Linit:Lfinal)[l]
      xvalnormal[o,l]=enso$index[final2]
      
    } 
    betamean[l]=summary.lm(lm(log(yval[,l])~xvalmean[,l]))$coefficients[2,1]
    rmean[l]=sqrt(summary.lm(lm(log(yval[,l])~xvalmean[,l]))$adj.r.squared)
    betamax[l]=summary.lm(lm(log(yval[,l])~xvalmax[,l]))$coefficients[2,1]
    rmax[l]=sqrt(summary.lm(lm(log(yval[,l])~xvalmax[,l]))$adj.r.squared)
    betanorm[l]=summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$coefficients[2,1]
    rnorm[l]=sqrt(summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$adj.r.squared)
    significancemean[l]=summary.lm(lm(log(yval[,l])~xvalmean[,l]))$coefficients[2,4]
    significancemax[l]=summary.lm(lm(log(yval[,l])~xvalmax[,l]))$coefficients[2,4]
    significancenorm[l]=summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  } 
  results=list()
  results$xvaluesmean=xvalmean
  results$xvaluesmax=xvalmax
  results$yvalues=yval
  results$sig=data.frame(lag, betamean,significancemean,rmean,betamax, significancemax, rmax,betanorm, significancenorm,rnorm)
  results$sig$index=index
  return(results)
  
}

#biolocal=climbiomassreg(bio=biomass, climseries=clim,variable="rain", monthname="month", yearname="year",Linit=-6, Lfinal=6) 
## this function calculates the regression among the mean values of an local climate variable over L months and the modelled values per month of seed production

###ddrain=localmonthlymodreg(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", beginyearfile=beginyearfile,spstart=1,spend=45, variable="rain", monthname="month", yearname="year",Linit=-6, Lfinal=6)
localmonthlymodreg=function(file="Nouragues model all spp.txt",fileno="NRG model all spp without Dj.txt", beginyearfile=beginyearfile,spstart=1,spend=45,variable="rain", monthname="month", yearname="year",Linit=-6, Lfinal=6)    #L is the number of months used to calculate the mean of ENSO,
  
{
  #setwd( "C:/Brunoy/Base Datos Nouragues/Metadata Joe Wright/modeling max likelihood/hierarchical models/results hierarchical models BCI-NRG")
  
  datayr=bio
  biomass$fecha=strptime(paste(biomass$Year2,"-",biomass$month,"-",1), format="%Y - %m - %d")
  biostart=which(biomass$Year2==2001 & biomass$month==7)
  bioend=which(biomass$Year2==2011 & biomass$month==1)
  biomass=biomass[biostart:bioend,]
  
  ind=which(colnames(climseries)==variable)#this selects the climatic variable of study
  var=climseries[,ind]
  ind2=which(colnames(climseries)==monthname)
  month=climseries[,ind2]
  ind3=which(colnames(climseries)==yearname)
  year=climseries[,ind3]
  julian=climseries$julian
  clim2=data.frame(julian,year,month, var)
  
  anomalclim<-numeric(length=length(clim2$var))
  
  for (z in 1:length(var))
  { 
    meanclim=aggregate(data.frame(index=clim2$var), by=list(month=clim2$month), mean)
    anomalclim[z]=clim2$var[z]-meanclim$index[meanclim$month==clim2$month[z]]
  }
  clim2$anomal=anomalclim
  plot(clim2$julian[2:length(clim2$julian)-2], las=1,runningmean(clim2$anomal, k=3), type="l", xlab="time", ylab="anomalies of monthly values of rainfall", bty="l")
  abline(h=0)
  
  
  franomal=flanomal=numeric(length=length(biomass$fruitsb))
  
  for (k in 1:length(biomass$fruitsb))
  { 
    meanbio=aggregate(data.frame(fr=biomass$fruitsb,fl=biomass$flowersb), by=list(month=biomass$month), mean)
    franomal[k]=biomass$fruitsb[k]-meanbio$fr[meanbio$month==biomass$month[k]]
    flanomal[k]=biomass$flowersb[k]-meanbio$fl[meanbio$month==biomass$month[k]]
  }
  biomass$franomal=franomal; biomass$flanomal=flanomal
  plot(biomass$fecha[2:length(biomass$fecha)-2], las=1,runningmean(biomass$franomal, k=3), type="l", col="red", xlab="time", ylab="anomalies", bty="l")
  lines(biomass$fecha[2:length(biomass$fecha)-2], runningmean(biomass$flanomal, k=3), col="blue")
  abline(h=0)
  legend(biomass$fecha[53],2.5, lty=c(1,1), col=c("red", "blue"), c("fruit biomass", "flower biomass"), bty="n")
  
  xval=matrix(nrow=length (biomass$fruitsb),ncol=length(c(Linit:Lfinal)))
  #colnames(xvalmean, xvalmax, yval)=as.character(c(Linit:Lfinal))
  yval=yflval=betafr=betafl=rfr=rfl=significancefr =significancefl= lag=numeric()
  
  for (l in 1:length(c(Linit:Lfinal)))
  {
    cat("L:",c(Linit:Lfinal)[l], "\n" )
    
    for(o in 1: length (biomass$fruitsb))
    {
      final=(which(clim2$year==biomass$Year2[o]&clim2$month==biomass$month[o]))-c(Linit:Lfinal)[l]
      xval[o,l]=clim2$anomal[final]
      yval[o]=biomass$franomal[o]
      yflval[o]=biomass$flanomal[o]
    } 
    betafr[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,1]
    rfr[l]=sqrt(summary.lm(lm(yval~xval[,l]))$adj.r.squared)
    betafl[l]=summary.lm(lm(yflval~xval[,l]))$coefficients[2,1]
    rfl[l]=sqrt(summary.lm(lm(yflval~xval[,l]))$adj.r.squared)
    significancefr[l]=summary.lm(lm(yval~xval[,l]))$coefficients[2,4]
    significancefl[l]=summary.lm(lm(yflval~xval[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  } 
  results=list()
  results$xvalues=xval
  results$yvalues=yval
  results$yvaluesflo=yflval
  results$sig=data.frame(variable,lag, betafr,significancefr, rfr,betafl, significancefl,rfl)
  return(results)
  
}

#file="Nouragues model all spp.txt";fileno="NRG model all spp without Dj.txt"; soifile="ENSO anomalies - 2012.txt";beginyearfile=beginyearfile;spstart=1;spend=45;index="MEI";Linit=-12; Lfinal=12; local=FALSE
#file="all parameters spp.txt";fileno="all parameters spp no Dj.txt";soifile="ENSO anomalies - 2012.txt";beginyearfile=beginyearfile;spstart=1;spend=45;index="MEI";monthname="MONTH"; yearname="YEAR"; Linit=-12; Lfinal=12; local=FALSE
#file="all parameters spp.txt";fileno="all parameters spp no Dj.txt"; soifile="local climate data Nouragues.txt";beginyearfile=beginyearfile;spstart=1;spend=45;index="rain";monthname="month"; yearname="year"; Linit=-5; Lfinal=5; local=TRUE

#ddrain=sppmodelsvsENSO(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", soifile="local climate data Nouragues.txt",beginyearfile=beginyearfile,spstart=1,spend=45,index="rain",monthname="month", yearname="year", Linit=-5, Lfinal=5, local=TRUE)
#dd=sppmodelsvsENSO(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", soifile="ENSO anomalies - 2012.txt",beginyearfile=beginyearfile,spstart=1,spend=45,index="MEI",monthname="MONTH", yearname="YEAR", Linit=-12, Lfinal=12, local=FALSE)
sppmodelsvsENSO=function(file="all parameters spp.txt",fileno="all parameters spp no Dj.txt", soifile="ENSO anomalies - 2012.txt",beginyearfile=beginyearfile,spstart=1,spend=45,index="MEI", monthname="MONTH", yearname="YEAR",Linit=-12, Lfinal=12, local=FALSE)
{
  datayr=monthlyvalues(file=file,fileno=fileno, beginyearfile=beginyearfile,spstart=spstart,spend=spend)
  spnames=sort(unique(datayr$species))
  soi<-read.delim(soifile,as.is=T)
  
  if (local==TRUE){
    season=aggregate(data.frame(rain=soi$rain, tmin=soi$tmin, tmax=soi$tmax), by=list(month=soi$mon), mean, na.rm=TRUE)
    anomalrain=anomaltmin=anomaltmax=numeric()
    for(m in 1:dim(soi)[1]) {
      
      anomalrain[m]=soi$rain[m]-season$rain[season$month==soi$month[m]]
      anomaltmin[m]=soi$tmin[m]-season$tmin[season$month==soi$month[m]]
      anomaltmax[m]=soi$tmax[m]-season$tmax[season$month==soi$month[m]]
      
    }
    soi=data.frame(year=soi$year, month=soi$month,rain= anomalrain, tmin=anomaltmin, tmax=anomaltmax) 
  }
    
    ind=which(colnames(soi)==index)#this selects the climatic variable of study
    var=soi[,ind]
    ind2=which(colnames(soi)==monthname)
    month=soi[,ind2]
    ind3=which(colnames(soi)==yearname)
    year=soi[,ind3]
    #julian=soi$julian
    #enso=data.frame(julian,year,month, var)
  enso=data.frame(year,month, var)

  allresults=list()
  for (j in 1:length(spnames))
  {
    species=as.character(spnames[j])
    cat("sp:",as.character(species), "\n" )
    onesp=datayr[datayr$species==species,]
     
    anomal<-numeric(length=length(onesp$model))
    
    for (k in 1:length(anomal))
    { 
      meanpeak=aggregate(data.frame(peak=onesp$model), by=list(month=onesp$month), mean)
      anomal[k]=onesp$model[k]-meanpeak$peak[meanpeak$month==onesp$month[k]]
    }
    onesp=data.frame(onesp, anomal)
    
  xvalmean=xvalmax=xvalnormal=yval=matrix(nrow=length(onesp$anomal),ncol=length(c(Linit:Lfinal)))
  betamax=betamean=betanorm=significancemean =significancemax=significancenorm= rmean=rmax= rnorm=lag=numeric()
    
  for (l in 1:length(c(Linit:Lfinal)))
  { 
   cat("L:",c(Linit:Lfinal)[l], "\n" )
   
    for (o in 1:dim(onesp)[1])
    {
      final=(which(enso$year==onesp$year[o]&enso$month==onesp$month[o]))-c(Linit:Lfinal)[l]
      inicio=which(enso$year==onesp$year[o]&enso$month==onesp$month[o])
      xvalmean[o,l]=mean(enso$var[inicio:final])
      xvalmax[o,l]=max(abs(enso$var[inicio:final]))
      yval[o,l]=onesp$anomal[o]
      xvalnormal[o,l]=enso$var[final]
      
    } 
 
    betamean[l]=summary.lm(lm(log(yval[,l])~xvalmean[,l]))$coefficients[2,1]
    rmean[l]=sqrt(summary.lm(lm(log(yval[,l])~xvalmean[,l]))$adj.r.squared)
    betamax[l]=summary.lm(lm(log(yval[,l])~xvalmax[,l]))$coefficients[2,1]
    rmax[l]=sqrt(summary.lm(lm(log(yval[,l])~xvalmax[,l]))$adj.r.squared)
    betanorm[l]=summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$coefficients[2,1]
    rnorm[l]=sqrt(summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$adj.r.squared)
    significancemean[l]=summary.lm(lm(log(yval[,l])~xvalmean[,l]))$coefficients[2,4]
    significancemax[l]=summary.lm(lm(log(yval[,l])~xvalmax[,l]))$coefficients[2,4]
    significancenorm[l]=summary.lm(lm(log(yval[,l])~xvalnormal[,l]))$coefficients[2,4]
    lag[l]=c(Linit:Lfinal)[l]
  }  
  results=list()
  results$xvaluesmean=data.frame(species,xvalmean)  
  results$xvaluesmax=data.frame(species,xvalmax)
  results$xvaluesnorm=data.frame(species,xvalnormal)
  results$yvalues=data.frame(species,yval)
  results$sig=data.frame(species, index, lag, betamean,significancemean,betamax, significancemax,betanorm, significancenorm)

  if(j==1) allresults$xvaluesmean= results$xvaluesmean else allresults$xvaluesmean=rbind(allresults$xvaluesmean,results$xvaluesmean)
  if(j==1) allresults$xvaluesmax= results$xvaluesmax else allresults$xvaluesmax=rbind(allresults$xvaluesmax,results$xvaluesmax)
  if(j==1) allresults$xvaluesnorm= results$xvaluesnorm else allresults$xvaluesnorm=rbind( allresults$xvaluesnorm,results$xvaluesnorm)
  if(j==1) allresults$yvalues= results$yvalues else allresults$yvalues=rbind(allresults$yvalues,results$yvalues)
  if(j==1) allresults$sig= results$sig else allresults$sig=rbind( allresults$sig,results$sig)
  }  
  return(allresults) 
}  

#summaryspmoldels do the summary of all the species-specific models among a climatic variable and the specific models of seed production
#norm means that the correlation was done with the corresponding value at a lag x, mean with the mean values among all the included lags and max with the absolute maximum value
#allsppMEI=read.delim(file="species-specific correlations MEI.txt")

summaryspmodels=function(data=allsppMEI){
 
  spnames=unique(data$species)
  for (i in 1:length(spnames))
    {
   species=spnames[i]
   onesp=data[data$species==species,]
   signormal=which(onesp$significancenorm<0.05)
   posnormal=which(onesp$betanorm[signormal]>0)
   negnormal=which(onesp$betanorm[signormal]<0)
   sigmean=which(onesp$significancemean<0.05)
   posmean=which(onesp$betamean[sigmean]>0)
   negmean=which(onesp$betamean[sigmean]<0)
   sigmax=which(onesp$significancemax<0.05)
   posmax=which(onesp$betamax[sigmax]>0)
   negmax=which(onesp$betamax[sigmax]<0)
  
   results=data.frame(index=unique(onesp$index),species, norm=length(signormal),posnormal=length(posnormal), negnormal=length(negnormal), mean=length(sigmean),posmean=length(posmean), negmean=length(negmean), max=length(sigmax), posmax=length(posmax),negmax=length(negmax))
  if (i==1) allresults=results else allresults=rbind(allresults, results)
  }
  
  return(allresults)
  
}

##sums per spp of the total number of produced seeds and contribution of the total community per species
##file="all parameters spp.txt";fileno="all parameters spp no Dj.txt"; beginyearfile=beginyearfile; spstart=1; spend=45



#### STANDARD NORMAL DEVIATE #####

SND = function (NRGmonthly=NRGmonthly){
  
  nrgyear=aggregate(data.frame(model=NRGmonthly$model),by=list(year=NRGmonthly$year,sp=NRGmonthly$species),sum)
  #peaksummary=aggregate(data.frame(peak=dat$peak), by=list(cyle=dat$cycle),unique)
  spp=unique(nrgyear$sp)
  for (i in 1:length(spp)){
    
    onesp=nrgyear[nrgyear$sp==spp[i],]
    snd=(log(onesp$model+1)-mean(log(onesp$model+1)) )/sd(log(onesp$model+1)) 
    if (i ==1) sndtot=snd else sndtot=c( sndtot,snd)
  }
  nrgyear$snd=sndtot
  meansnd=aggregate(data.frame(sndmean=nrgyear$snd),by=list(year=nrgyear$year),mean)
  sdsnd=aggregate(data.frame(sdsnd=nrgyear$snd),by=list(year=nrgyear$year),sd)
  plot(meansnd$year[2:10],meansnd$sndmean[2:10],type="o",ylim=c(-1.3,1.3),xlab="years",ylab="Seed production (SND)",las=1,lwd=2)
  segments(sdsnd$year[2:10], meansnd$sndmean[2:10] - sdsnd$sdsnd[2:10], sdsnd$year[2:10], meansnd$sndmean[2:10] + sdsnd$sdsnd[2:10], lwd = 2, lend = 3, col = "black")
 
}

#SNDmodel calculates the Standard Normal Deviate of the model of seed production of each species
# and then plots the sum per month of the values of all spp

SNDmodel = function (nrgyear=NRGmonthly){
  spp=unique(nrgyear$sp)
  for (i in 1:length(spp)){    
    onesp=nrgyear[nrgyear$sp==spp[i],]
    onesp$snd=(log(onesp$model+1)-mean(log(onesp$model+1)) )/sd(log(onesp$model+1)) 
    if (i ==1) nrgdata=onesp else nrgdata=rbind(nrgdata,onesp)
  }

  meansnd=aggregate(data.frame(snd=nrgdata$snd),by=list(month=nrgdata$month,year=nrgdata$year),sum)
  #meansnd=aggregate(data.frame(snd=nrgdata$snd),by=list(year=nrgdata$year),mean)
  meansnd=meansnd[-as.vector(which(meansnd$year==2001)),]
  sdsnd=aggregate(data.frame(snd=nrgdata$snd),by=list(month=nrgdata$month,year=nrgdata$year),sd)
  sdsnd=sdsnd[-as.vector(which(sdsnd$year==2001)),]
  plot(meansnd$snd,type="o",xlab="years",ylab="Sum of seed production (SND)",las=1,lwd=2,axes=FALSE)
  axis(side=2,las=2)
  axis(side=1, at=seq(1,120,12),labels=as.character(c("2002","2003","2004","2005","2006","2007","2008","2009","2010","2011")))
  abline(lm(meansnd$snd~c(1:108)))
  #segments(c(1:108), meansnd$snd - sdsnd$snd,c(1:108), meansnd$snd+ sdsnd$snd, lwd = 1, lend = 3, col = "black")
  
}



##### CLIMATIC DATA#################

meteo <- read.delim(file="MeteoFrance-Nouragues.txt")
meteo$dates = strptime(paste("1-",meteo$Month,"-",meteo$Year), format="%d - %m - %Y")
meteo$julian = tojulian(meteo$dates,dateform = "%Y-%m-%d")
meteo$yday = meteo$dates$yday+1

rochambeau <- meteo[meteo$Station=="Rochambeau",]
stgeorges <- meteo[meteo$Station=="Saint Georges",]
saul <- meteo[meteo$Station=="Saul",]
camopi <- meteo[meteo$Station=="Camopi",]
noura <- meteo[meteo$Station=="Nouragues",]
regina<-meteo[meteo$Station=="Regina",]
roura <- meteo[meteo$Station=="Roura",]



### 1- correlation between irradiation data and rainfall or temperature

plotradiationtemp=function(meteofile="MeteoFrance-Nouragues.txt",station="Saint Georges"){
  meteo<-read.delim(meteofile)
  meteo$dates=strptime(paste("1-",meteo$Month,"-",meteo$Year), format="%d - %m - %Y")
  meteo$julian=tojulian(meteo$dates,dateform = "%Y-%m-%d")
  meteo$yday=meteo$dates$yday+1

rochambeau<-meteo[meteo$Station==station,]

par(mfrow=c(2,1),las=1,mar=c(2,3,2,3),oma=c(2,2,2,2),cex=1)
plot(rochambeau$insolat,type="n",ylab="",xlab="",las=1,bty="l",axes=F)
lines(rochambeau$insolat, lwd=2)
axis(side=2, las=1,col="black",cex.axis=0.55)
jan1=which(rochambeau$Month==1)
axis(side=1, at=c(jan1[1]-12,jan1,jan1[length(jan1)]+12),labels=as.character(c("2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")))
mtext(side=2, "monthly insolation (minutes)", line=3, cex=1,las=0)
legend(2,18500, lty=c(1,1,1), c("insolation","mean temperature","rainfall"), bty="n",col=c("black", "red","blue"), horiz=T,lwd=2,xpd=T)
par(new=T)
plot(rochambeau$TM,type="l", lty=1,ylab="",xlab="",axes=F, las=1, col="red", lwd=2)
axis(side=4,cex.axis=0.75)
mtext(side=4, "mean temperature (?T)", line=3, cex=1,las=0)
mtext(side=3, station, line=1.5, cex=1.1,las=0,font=2)
#mtext(side=3, "Saint Georges", line=2, cex=1.1,las=0)
  
plot(rochambeau$insolat,type="n",ylab="",xlab="",las=1,bty="l",axes=F)
lines(rochambeau$insolat, lwd=2)
axis(side=2, las=1,col="black",cex.axis=0.55)
jan1=which(rochambeau$Month==1)
axis(side=1, at=c(jan1[1]-12,jan1,jan1[length(jan1)]+12), labels=as.character(c("2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")),cex=1)
mtext(side=2, "monthly insolation (minutes)", line=3, cex=1,las=0)
par(new=T)
plot(rochambeau$rainfall,type="l", lty=1,ylab="",xlab="",axes=F, las=1, col="blue", lwd=2)
axis(side=4,cex.axis=0.75)
mtext(side=4, "rainfall (mm)", line=3, cex=1,las=0)

}

##anomalies radiation
radanomalies=function(meteofile="MeteoFrance-Nouragues.txt",station="Rochambeau"){
  
  meteo<-read.delim(meteofile)
  meteo$dates=strptime(paste("1-",meteo$Month,"-",meteo$Year), format="%d - %m - %Y")
  meteo$julian=tojulian(meteo$dates,dateform = "%Y-%m-%d")
  meteo$yday=meteo$dates$yday+1
  
  rochambeau<-meteo[meteo$Station==station,]
  
  radmean=aggregate(data.frame(insolat=rochambeau$insolat), by=list(month=rochambeau$Month), mean, na.rm=T)
  sdmean=aggregate(data.frame(insolat=rochambeau$insolat), by=list(month=rochambeau$Month), sd, na.rm=T)
  
  for (i in 1:nrow(rochambeau))
    
  {
    
  rochambeau$anomal[i]=  rochambeau$insolat[i]-radmean$insolat[radmean$month[(radmean$month %in% rochambeau$Month[i])]]
    
  }

return(rochambeau)
}

###plot standard values of radiation plus anomalies 

plotanomalies=function(meteofile="MeteoFrance-Nouragues.txt",station="Rochambeau"){
anomal<-radanomalies(meteofile=meteofile,station=station)

par(mfrow=c(2,1),las=1,mar=c(2,3,2,3),oma=c(2,2,2,2),cex=1)
plot(rochambeau$insolat,type="n",ylab="",xlab="",las=1,bty="l",axes=F)
lines(rochambeau$insolat, lwd=2)
axis(side=2, las=1,col="black",cex.axis=0.55)
jan1=which(rochambeau$Month==1)
axis(side=1, at=c(jan1[1]-12,jan1,jan1[length(jan1)]+12),labels=as.character(c("2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")))
mtext(side=2, "monthly insolation (minutes)", line=3, cex=1,las=0)
legend(2,18500, lty=c(1,1,1), c("insolation","mean temperature","rainfall"), bty="n",col=c("black", "red","blue"), horiz=T,lwd=2,xpd=T)
par(new=T)
plot(rochambeau$TM,type="l", lty=1,ylab="",xlab="",axes=F, las=1, col="red", lwd=2)
axis(side=4,cex.axis=0.75)
mtext(side=4, "mean temperature (?T)", line=3, cex=1,las=0)
mtext(side=3, station, line=1.5, cex=1.1,las=0,font=2)

plot(anomal$anomal,type="n",ylab="",xlab="",las=1,bty="l",axes=F)
lines(anomal$anomal, lwd=2)
axis(side=2, las=1,col="black",cex.axis=0.55)
jan1=which(rochambeau$Month==1)
axis(side=1, at=c(jan1[1]-12,jan1,jan1[length(jan1)]+12),labels=as.character(c("2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")))
}

sumanomalies=function(dataset,k=3) #this function sums k number of observations; 
  #dataset must be a vector including the variable of interest  
{
  N<-length(dataset)
  total<-N-k+1
  
  counter=1
  val=vector()
  
  for (i in 1:total)
  {
    val[i]<-sum(dataset[i:(k+i-1)])
    counter=counter+1
  }  
  return(val)
  
}

#par(mfrow=c(3,1))
#plot(sumanomalies(anomal$anomal, k=3),type="l")
#plot(sumanomalies(anomal$anomal, k=9),type="l")
#plot(sumanomalies(anomal$anomal, k=11),type="l")

plotradyr=function(meteofile="MeteoFrance-Nouragues.txt",station="Rochambeau"){
  meteo<-read.delim(meteofile)
  meteo$dates=strptime(paste("1-",meteo$Month,"-",meteo$Year), format="%d - %m - %Y")
  meteo$julian=tojulian(meteo$dates,dateform = "%Y-%m-%d")
  meteo$yday=meteo$dates$yday+1
  
  nyear=unique(rochambeau$Year)
  
  rochambeau<-meteo[meteo$Station==station,]
  radmean=aggregate(data.frame(insolat=rochambeau$insolat), by=list(month=rochambeau$Month), mean, na.rm=T)
  
  par(mar=c(4,4,3,3),mfrow=c(1,1))
  plot(rochambeau$insolat[rochambeau$Year==nyear[1]], type="n",ylim=c(0,16000),ylab="minutes of insolation",xlab="month of the year",bty="l",axes=F)
  axis(side=1,at=1:12, labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
  axis(side=2)
  
 for (i in 1:lengthunique(rochambeau$Year)) {
   lines(rochambeau$insolat[rochambeau$Year==nyear[i]],col=i)
   }
 lines(radmean$insolat,lwd=3) 
}
ccfinsolationclimate=function(meteofile="MeteoFrance-Nouragues.txt"){
  par(mfrow=c(2,1),las=1,mar=c(3,3,4,3),oma=c(2,2,2,2),cex=1)
  aa=ccf(rochambeau$insolat, rochambeau$TM,lag.max=12)
  ab=ccf(rochambeau$insolat, rochambeau$rainfall,lag.max=12)
  
}

### Analyzing raw data of climate from Nouragues (manual station)
rawclimate=function(climraw=climraw, autoraw=autoraw){
rawna=aggregate(data.frame(rain=climraw$Rainfall,tmin=climraw$TempMin, tmax=climraw$TempMax), by=list(month=climraw$Month, year=climraw$Year), lengthisna)##how many missing values are in the manual station?
rainraw=aggregate(data.frame(rain=climraw$Rainfall), by=list(month=climraw$Month, year=climraw$Year), sum)
rainmean=aggregate(data.frame(rain=rainraw$rain), by=list(month=rainraw$month), mean, na.rm=T)

##comparing manual vs. automatic station

commonraw=climraw[which(climraw$julian %in% autoraw$julian),]
commonautoraw=autoraw[which(autoraw$julian %in% climraw$julian),]
compautoman=data.frame(year=commonraw$Year,month=commonraw$Month,day=commonraw$Day,
                       julian=commonraw$julian,tminman=commonraw$TempMin, tminauto=commonautoraw$Min.AirT.,
                       tmaxman=commonraw$TempMax,tmaxauto=commonautoraw$Max.AirT.,
                       rainman=commonraw$Rainfall,rainauto=commonautoraw$Prec)

#compautoman is the common points from the manual and the automatic stations
par(mar=c(3,3,2,2),oma=c(4,4,4,4))
plot(compautoman$tminauto,compautoman$tminman, xlab="tmin auto", ylab="tmin manual")
abline(lm(compautoman$tminman~compautoman$tminauto))
summary(lm(compautoman$tminman~compautoman$tminauto))

plot(compautoman$tmaxauto, compautoman$tmaxman, xlab="tmax auto", ylab="tmax manual")
abline(lm(compautoman$tmaxman~compautoman$tmaxauto))
summary(lm(compautoman$tmaxman~compautoman$tmaxauto))

#we define a threshold of 11 mm for the differences in rainfall between the manual and automatic station
rainauto=compautoman$rainauto[-which(abs(compautoman$rainauto-compautoman$rainman)>=11)]
rainman=compautoman$rainman[-which(abs(compautoman$rainauto-compautoman$rainman)>=11)]

plot(compautoman$rainauto[-which(abs(compautoman$rainauto-compautoman$rainman)>=11)], compautoman$rainman[-which(abs(compautoman$rainauto-compautoman$rainman)>=11)], xlab="rainfall auto", ylab="rainfall manual")
abline(lm(rainman~rainauto))
summary(lm(rainman~rainauto))

##we construct a new dataframe with the estimated values of the manual station from the automatic
newman=data.frame(climraw,tminauto=autoraw$Min.AirT.[match(climraw$julian,autoraw$julian)],
                  tmaxauto=autoraw$Max.AirT.[match(climraw$julian,autoraw$julian)],
                  rainauto=autoraw$Prec[match(climraw$julian,autoraw$julian)])

##we estimate the missing values of the manual station according to the regression between automatic and manual
newman$tmine=newman$TempMin
newman$tmine[which(is.na(newman$tminauto)==FALSE&is.na(newman$TempMin==TRUE))]=1.28996+(newman$tminauto[which(is.na(newman$tminauto)==FALSE&is.na(newman$TempMin==TRUE))]*0.95469)

newman$tmaxe=newman$TempMax
newman$tmaxe[which(is.na(newman$tmaxauto)==FALSE&is.na(newman$TempMax==TRUE))]=-0.63600+(newman$tmaxauto[which(is.na(newman$tmaxauto)==FALSE&is.na(newman$TempMax==TRUE))]*1.07511)

newman$raine=newman$Rainfall
newman$raine[which(is.na(newman$rainauto)==FALSE&is.na(newman$Rainfall==TRUE))]=0.129037+(newman$rainauto[which(is.na(newman$rainauto)==FALSE&is.na(newman$Rainfall==TRUE))]*0.914345)

write.table(newman, file="Nouragues manual estimated new.txt",sep="\t",row.names=F)

###newmeans

newmanmean=aggregate(data.frame(tmin=newman$tmine, tmax=newman$tmaxe), by=list(month=newman$Month, year=newman$Year), mean, na.rm=T)
newmanrain=aggregate(data.frame(rain=newman$raine),by=list(month=newman$Month, year=newman$Year),sum)
}

### 2- Correlations between climatic stations (Nouragues vs the rest from Meteo-France)
#newman=read.delim(file="Nouragues manual estimated new.txt")

##I would need daily data from Meteo-France for doing the correlations
comparisonNrgMeteo<-function(meteo=meteo, climraw=climraw){
  
commonroch=climraw[which(climraw$julian %in% autoraw$julian),]
lm1<-lm(climraw$Rainfall~rochambeau$rainfall)

}

### STRENGTH OF THE DRY SEASON
## I need to compare the two datasets used: newman (manual + estimations from the automatic statioin) vs. "local climate data Nouragues.txt"

dryNouragues<-function(dat=newman){
  
  ## we calculate first how many days of information were lacking of the manual rainfall dataset
  dry=aggregate(data.frame(rain=dat$raine), by=list(Month=dat$Month,Year=dat$Year), sum,na.rm=T)
  for (k in 1:dim(drynolack)[1]){   
    dry$lack[k]=length(which(is.na(dat$raine[dat$Year==dry$Year[k]&dat$Month==dry$Month[k]])==T))
  }
# we establish a threshold of 4 days for keeping rainfall data of that month
 dry2=dry[-which(dry$lack>4),]
drymean=aggregate(data.frame(rain=dry2$rain),by=list(Month=dry2$Month),mean)
sdmonthrain=aggregate(data.frame(rain=dry2$rain),by=list(Month=dry2$Month),sd)
meanrain=aggregate(data.frame(rain=dry2$rain),by=list(Year=dry2$Year),sum)
dat2=read.delim(file="local climate data Nouragues.txt")
drymean2=aggregate(data.frame(rain2=dat2$rain),by=list(Month=dat2$month),mean)
merge(drymean,drymean2,by="Month")
}
