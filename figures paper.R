####FIGURES OF THE PAPER####

####DATASET####
attach('CTFSRPackage.rdata')

library(ggplot2)
library(dplyr)



####  FIGURE 1 OF THE PAPER   ######################################
### beautiful graph of mean number of species  and climatogram ###

#figure1(datfile="new",graphname="figure1.tif",est=est,biomass=biomass)
figure1 = function(graphname="figure1.tif", est = est, biomass = biomass, newman = "Nouragues manual estimated new.txt")
  #newman is the estimated dataset of local climate from the climatic station of Nouragues
  #old datafile referes to "local climate data Nouragues.txt"
{
  #dat <- read.table("local climate data Nouragues.txt", header = T) #this was for the previous version of the paper
  
    dat = read.delim(newman)
    dry = aggregate(data.frame(rain=dat$raine), by = list(Month=dat$Month,Year=dat$Year), sum,na.rm=T)
    for (k in 1:dim(dry)[1]){   
      dry$lack[k]=length(which(is.na(dat$raine[dat$Year==dry$Year[k]&dat$Month==dry$Month[k]])==T))
    }
    # we establish a threshold of 4 days for keeping rainfall data of that month
    dry2 = dry[-which(dry$lack>4),]
    monthrain = aggregate(data.frame(rain=dry2$rain),by=list(month=dry2$Month),mean)
    sdrain = aggregate(data.frame(rain=dry2$rain),by=list(month=dry2$Month),sd)
    
    means <- data.frame(aggregate(data.frame(tmin=dat[,16],tmax=dat[,17]), by = list(month = dat$Month), mean,na.rm=T),rain=monthrain$rain)
    sds  <- data.frame(aggregate(data.frame(tmin=dat[,16],tmax=dat[,17]), by = list(month = dat$Month), sd,na.rm=T),rain=sdrain$rain)
  }
  
  labels <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
  
  #est<-read.delim("number total spp per month estimated.txt",header=T)
  sppmeans = aggregate(est$estnumbspp, by = list(month = est$month), mean)
  names(sppmeans)=c("month","nbspp")
  sppsds=aggregate( est$estnumbspp, by = list(month = est$month), sd)
  names(sppsds)=c("month","nbspp")
  
  #biomass=read.delim(file="biomass all months.txt",as.is=T)
  biomass = biomass[-1,]
  meanbio = aggregate(data.frame(fruit=biomass$fruits,flower=biomass$flowers), by=list(month=biomass$month), mean)
  sdbio = aggregate(data.frame(fruit=biomass$fruits,flower=biomass$flowers), by=list(month=biomass$month), sd)
  
  #x11(height=10,width=6) 
  tiff(filename=graphname,width = 1900, height = 3000,pointsize=12, res=300)
  par(las = 1, bty = "o", tcl = 0.2, mar = c(2, 5, 1, 5), oma=c(2,3,1,1), mgp = c(0.25, 0.25, 0),cex.axis=1.5,lwd=1.5)
  par(mfrow=c(3,1))
  
  plot(meanbio$fruit,type="o",axes=FALSE, xlab="",ylab="", col="black", ylim=c(0,15) )
  segments(1:12, meanbio$fruit- sdbio$fruit, 1:12, meanbio$fruit + sdbio$fruit, lwd = 2, lend = 3, col = "black")
  axis(1, at = 1:12, labels = NA)
  axis(2, cex.axis = 2) 
  mtext(text=expression(paste("Biomass (g ",m^-2, ")", sep = " ")), 2, las = 3, line = 4, cex = 1.5) 
  text(2,15,labels="A",pos=2, offset=2, cex=2)
  par(new=TRUE)
  ji12=jitter(c(1:12))
  plot(ji12, meanbio$flower, type="o",axes=FALSE, xlab="", ylab="", col="grey70", ylim=c(0,15) ,lwd=2)
  segments(ji12, meanbio$flower- sdbio$flower, ji12, meanbio$flower + sdbio$flower, lwd = 2, lend = 3, col = "grey70")
  legend(7,15,legend=c("fruits","flowers"), col=c("black","grey70"),lty=1,lwd=3,bty="n",cex=2)
  
  print(cor.test(means$rain, meanbio$flower))
  print(cor.test(means$rain, meanbio$fruit))
  
  plot(sppmeans$nbspp,type="o",axes=FALSE, ylim=c(10,40), xlab="",ylab="",col="black")
  text(2,40,labels="B",pos=2, offset=2, cex=2)
  axis(1, at = 1:12, labels = NA)
  axis(2, cex.axis = 2)
  segments(1:12, sppmeans$nbspp - sppsds$nbspp, 1:12, sppmeans$nbspp + sppsds$nbspp, lwd = 2, lend = 3, col = "black")
  mtext("Number of fruiting species", 2, las = 3, line = 4, cex = 1.5)
  print(cor.test(means$rain,sppmeans$nbspp))
  
  #plot(1:12, means$rain, axes = F, type = "n", xlab = "", ylab = "", ylim = range(means$rain - sds$rain, means$rain+ sds$rain))
  plot(1:12, means$rain, axes = F, type = "n", xlab = "", ylab = "", ylim = c(0,800))
  text(2,800,labels="C",pos=2, offset=2, cex=2)
  polygon(c(1:12, 12:1), c(means$tmin, rev(means$tmax))*20, col = "#FF000050", border = NA)
  segments( 1:12, 0, 1:12, means$rain, lwd = 40, lend = 3, col = "blue",lend=3)
  segments(1:12, means$rain - sds$rain, 1:12, means$rain+ sds$rain, lwd = 2, lend = 3, col = "black")
  
  
  axis(1, at = 1:12, labels = labels, cex.axis = 1.7)
  axis(2, at = 100*1:5, col.axis = "blue", cex.axis = 2)
  axis(4, at = 200*1:4, labels = 10*1:4, col.axis = "red", cex.axis = 2)
  mtext("Temperature (ºC)", 4, las = 3, line = 3.5, col = "red", cex = 1.5)
  mtext("Precipitation (mm)", 2, las = 3, line = 4, col = "blue", cex = 1.5)
  
  dev.off() 
}


####FIGURE 2 OF THE PAPER ###############

#figure2 plots the hyperparameters for each species in Nouragues (in an exponential scale)
#figure2(trfile = "Nouragues results hyperparameters.txt", longnames = "total number of seeds per species.txt", filename = "figure2.tif")
figure2 = function(trfile = "Nouragues results hyperparameters.txt", longnames = "total number of seeds per species.txt", filename = "figure2.tif")
  
{
  
  #x11(height=9,width=9)
  tiff(filename=filename, height = 1200,width=2000,res=300)
  par(mar=c(3,3,1,1),oma=c(2,2,1,2))
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
  plot(1:365,oneyrpred,type='l',ylab='',xlab='',ylim=c(-0.5,max(oneyrpred+3)),lwd=1.25, cex=1, axes=F)
  axis(side=2, las=2, at = seq(0,35,5), labels = round(seq(0,35,5)/80,2))
  axis(side=1, at=c(1,32,60,92,122,152,183,214,245,275,306,336,366), labels=c("J", "F", "M", "A","M", "J", "J", "A", "S","O","N","D","J"))
  #polygon(c(248,341,341,248),c(0,0,3,3),col = "#FF000020", border = NA)
  #polygon(c(0,248,248,0),c(0,0,3,3),col = "#0000FF20", border = NA)
  polygon(c(212,334,334,212),c(0,0,30,30),col = "grey80", border = NA)
  #polygon(c(341,400,400,341),c(0,0,3,3),col = "#0000FF20", border = NA)
  #mtext(3, text="Nouragues (2001-2011)",line=2,cex=1.5)
  mtext(text = expression(paste("biweekly seedfall (seeds* ",m^-2, ")", sep = " ")), 2, las = 3, line = 3, cex = 1)
  mtext(1,text="month of the year",line=3,cex=1)
  for (i in 1:dim(tr)[1])
  {
    oneyrpred=14*exp(tr$logmu[i])*dnorm(1:365,mean=tr$mu2[i],sd=tr$bestSD[i])
    lines(1:365,oneyrpred,col=colores[i], lwd=2) 
  }
  
  lines(1:365, 14*mean(exp(tr$logmu))*dnorm(1:365,mean=mean(tr$mu2),sd=mean(tr$bestSD)),col="black",lwd=4)
  legend(10,39,lty=c(1,1),col=c("black","red","blue"),legend=c("community","biotic","abiotic"),bty="n",cex= 1.2,horiz=T,lwd=2)
  dev.off()
  
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
  
  disptest=matrix(ncol=2,nrow=2)
  colnames(disptest)=c("rainypeak","rest")
  rownames(disptest)=c("biotic","abiotic")
  disptest[1,]=c(11,13)
  disptest[2,]=c(11,10)
  biotic<-c(11,1,5,7)
  abiotic<-c(11,3,3,4)
  chisq.test(disptest) ## no signficant differences between dispersal mode of species
  

}
