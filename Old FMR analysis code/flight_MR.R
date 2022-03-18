###  Analysis of metabolic rate, somewhat taken from Michael Dillon and Abbie Reade (MDAR) scripts...
# for analysis in R of Expedata files previously exported as .csv 
#  
# These are for Stephen's honeybee experiments, flying chamber bees, 350 mL chamber, with channels:
#
#    
#    Seconds: seconds, starts with "0"
#    Barometric_Pressure:  pressure (kPa)
#    Flow_Rate: flow rate in system (mL/min)
#    Oxygen:   Oxygen %, corrected for barometric pressure
#    Water_vapor:   Water vapor (kPa) *******(mmol/mol so divide by 1000 to get molar fraction)
#    Carbon_Dioxide:   CO2 %, corrected for barometric pressure
#
### Last Edited SM 8/5/18


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#set the working directory, here, it is reinforcing the default project directory
setwd("C:/Users/jcassano/Desktop/Masters Thesis/Summer 2020 Thesis Project/FMR Analysis/FMR Raw Data")

# directory for visual checks of data
dir.checks <- "output/checktrace"
# directory where all corrected data is housed for further analysis
dir.processed <- "output/processed"


library(dplyr)
library(stringr)
library(foreach)
library(segmented)
library(grDevices)
library(tidyr)
library(zoo)


# Functions
load.expedata.file <- function(){
  file <- file.choose()
  # extract information from filename
  filename <<- unlist(strsplit(basename(file), "[ ]"))[1]
  filename <<- gsub(".csv", "", tail(strsplit(file, "/")[[1]], n=1))
  #edate    <<- strsplit(filename, "_")[[1]][2]
  # read in data to "flightfull" which will be full resting data
  flightfull   <- read.csv(file, header=T)
  # Markers .. these are saved as Ascii codes, so here we convert them back to characters
  flightfull$markers <- rawToChar(as.raw(as.character(flightfull$Marker)), multiple=T)
  # this returns decimal fraction, uncorrected for water dilution because lags have to be fixed first
  #flightfull$CO2uc <- flightfull$CO2/1000000 
  #unique(flightfull$chamber)
  return(flightfull)
  #write.csv(flightfull, "test.csv", row.names=F)
}





flightfull <- load.expedata.file()
str(flightfull)
View(flightfull)

#create new marker names
#flightfull$Markname <- ifelse(flightfull$Marker == -1, flightfull$Markname <- "none", 
#                           ifelse(flightfull$Marker == 108, flightfull$Markname <- "end", 
#                              ifelse(flightfull$Marker == 118, flightfull$Markname <- "beginbase", flightfull$Markname <- "check")))



#Timepoint determination
#This code populates the blanks in markers with NA's then fills in every NA with a previous actual value entry, 
#                                 effectively making each row fall into a timepoint position based on the markers column
flightfull$markers[1]<-"K"
flightfull$markers[flightfull$markers == ""] <- NA
flightfull$markers <- na.locf(flightfull$markers)
flightfull$markers <- toupper(flightfull$markers)

# saves for later use
# Ds <- which(rdata$markers != "") # switches between baseline and experimental by Daemon???? <<<<< Don't think I have this

seal <- which(flightfull$markers == "V") # seal, continues until 1st behavior noted
fly <- which(flightfull$markers == "F") #flying
bsln <- which(flightfull$markers == "K") #Baseline at beginning
rest <- which(flightfull$markers == "R") #rest


#first of each section
seal1 <- seal[1]
fly1 <- fly[1] #flying
bsln1 <- bsln[1] #Baseline at beginning
rest1 <- rest[1] #rest

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#lag correct 
plot(flightfull$Seconds[(seal1-20):(seal1+20)], flightfull$Carbon_Dioxide[(seal1-20):(seal1+20)])
abline(v=seal1)

#This creates a new column with "shifted" or lag-corrected readings to match the time point
lagpres <- 10
lagflow <- 10
lagwvp <- 10
lago2 <- 10
lagco2 <- 10

flightfull$co2nolag <- c(flightfull$Carbon_Dioxide[(1+lagco2):length(flightfull$Carbon_Dioxide)], rep(NA, lagco2))
flightfull$wvpnolag <- c(flightfull$Water_Vapor[(1+lagwvp):length(flightfull$Water_Vapor)], rep(NA, lagwvp))
flightfull$o2nolag <- c(flightfull$Oxygen[(1+lago2):length(flightfull$Oxygen)], rep(NA, lago2))
flightfull$flownolag <- c(flightfull$Flow_Rate[(1+lago2):length(flightfull$Flow_Rate)], rep(NA, lagflow))
flightfull$presnolag <- c(flightfull$Barometric_Pressure[(1+lagpres):length(flightfull$Barometric_Pressure)], rep(NA, lagpres))

#plot full data lag corrected
plot(flightfull$Seconds, flightfull$wvpnolag, type="l", col="blue")
par(new=T)
plot(flightfull$Seconds, flightfull$co2nolag, type="l", col="darkorange", axes=T, ann=F)
par(new=T)
plot(flightfull$Seconds, flightfull$o2nolag, type="l", col="green", axes=T, ann=F)
par(new=T)
plot(flightfull$Seconds, flightfull$flownolag, type="l", axes=T, ann=F)
abline(v=fly1-10)


## Check lag visually in PDF, this chunk graphs it ----

lagfile <- paste(filename, "lags.pdf", sep="-")

pdf(file = paste(dir.checks, lagfile, sep="/"), width=7, height=7)

#set graphical parameters
par(mfrow=c(3,1), plt=c(0.1, 0.9, 0.1, 0.9))
# plot of entire sequence, with metadata
plot(flightfull$Seconds, flightfull$wvpnolag, type="n", ann=F, axes=F)
axis(side=1,at=seq(0,1000,by=25))
abline(v=(lagco2+61), lwd=1, lty=2, col="dark orange")
abline(v=(lagwvp+61), lwd=1, lty=2, col="blue")
# grayed baseline boxes
leftx  <- chbr
boty   <- rep(min(flightfull$wvpnolag, na.rm=T), length(leftx))
rightx <- bsln
topy   <- rep(max(flightfull$wvpnolag, na.rm=T), length(rightx))
rect(leftx, boty, rightx, topy, col=rgb(0.5, 0.5, 0.5, 0.1), border=NA)
#abline(v=Ds, col="gray")
par(new=T)
plot(flightfull$Seconds, flightfull$wvpnolag, type="l", col="lightblue", axes=T, ann=F)
par(new=T)
plot(flightfull$Seconds, flightfull$co2nolag, type="l", col="darkorange", axes=T, ann=F)
top <- max(flightfull$co2nolag, na.rm=T)
#timerange <- paste(format(flightfull$Time[1], format="%m-%d-%Y %I:%M %p"), 
#format(tail(flightfull$Time, n=1), format="%I:%M %p"), sep=" - ")
#text(rightx, top, adj=c(1, 0), cex=0.8, timerange)
legend(x=rightx[1], y=top, bty="n", cex=0.8,
       legend=c(paste(lagwvp, "s WVP lag"), paste(lagco2, "s CO2 lag")),
       text.col=c("lightblue", "darkorange"))
text(chbr, max(flightfull$co2nolag, na.rm=T), adj=c(0.5, 0.5), chbr)
text(bsln, min(flightfull$co2nolag, na.rm=T), adj=c(0.5, 0.5), bsln)
# now zoom in on first switch to chamber
r <- (chbr-30):(chbr+30)
yl <- range(flightfull$wvpnolag[r], na.rm=T) + c(-.1, .1)
plot(flightfull$Seconds[r], flightfull$wvpnolag[r], type="n", axes=T, ann=F, xlim=range(r))
abline(v=(lagco2+61), lwd=1, lty=2, col="dark orange")
abline(v=(lagwvp+61), lwd=1, lty=2, col="blue")
rect(chbr, yl[1], max(r), yl[2], col=rgb(0.5, 0.5, 0.5, 0.2), border=NA)
text(chbr, yl[2]-0.1, adj=c(0.5, 0.5), chbr)
par(new=T)
plot(flightfull$Seconds[r], flightfull$wvpnolag[r], type="o", col="lightblue", lwd=0.4, axes=T, ann=F, xlim=range(r))
par(new=T)
plot(flightfull$Seconds[r], flightfull$co2nolag[r], type="o", col="darkorange", lwd=0.4, axes=T, ann=F)
abline(v=chbr, col=rgb(0.5, 0.5, 0.5, 0.5), lwd=0.4)
# now zoom in on switch to baseline
r <- (bsln-30):(bsln+30)
yl <- range(flightfull$wvpnolag[r], na.rm=T) + c(-.1, .1)
plot(flightfull$Seconds[r], flightfull$wvpnolag[r], type="n", axes=T, ann=F, xlim=range(r))
abline(v=(lagco2+661), lwd=1, lty=2, col="dark orange")
abline(v=(lagwvp+661), lwd=1, lty=2, col="blue")
rect(min(r), yl[1], bsln, yl[2], col=rgb(0.5, 0.5, 0.5, 0.2), border=NA)
text(bsln, yl[1]+0.1, adj=c(0.5, 0.5), bsln)
par(new=T)
plot(flightfull$Seconds[r], flightfull$wvpnolag[r], type="o", col="lightblue", lwd=0.4, axes=T, ann=F)
par(new=T)
plot(flightfull$Seconds[r], flightfull$co2nolag[r], type="o", col="darkorange", lwd=0.4, axes=T, ann=F)
abline(v=chbr, col=rgb(0.5, 0.5, 0.5, 0.5), lwd=0.4)

dev.off()

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

## Correct CO2 for WVP dilution ----
#  units in kPa
# correction factor: BP/(BP-WVP), where
# BP = barometric pressure, column was measured! 
# WVP = measured "H2O" column

# correction is uneccessary, see p 96 of Lighton 2008 
#flightfull$co2_wvpc <- flightfull$co2nolag * flightfull$presnolag/(flightfull$presnolag-flightfull$wvpnolag)
flightfull$co2_wvpc <- flightfull$co2nolag

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





## Drift correct and zero CO2 ----

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''
allVarBase <- data.frame(startTime=integer(0), varValBase=numeric(0))

for(x in bsln1:(seal1-5)){ 
  varNowBase <- var(restfull$co2_wvpc[x:(x+30)], na.rm = T)
  newRowBase <- c(x, varNowBase)
  allVarBase <- rbind(allVarBase, newRowBase)  
}

names(allVarBase)=c("startTime", "varValBase")
#View(allVar)
#now to find minumum variance point using the above time frame in the for loop (2 minutes) 
minVarBase <- min(allVarBase$varValBase)

#Get the estimates 
begin.base <- allVarBase$startTime[which(allVarBase$varValBase==min(allVarBase$varValBase))]
duration.base <- 30
lowVar.baseline <- begin.base:(begin.base+duration.base)
prelowvar <- mean(flightfull$co2_wvpc[lowVar.baseline], na.rm=T)

# drift correction
shiftedco2 <- flightfull$co2_wvpc - prelowvar
# subtracting out pre-baseline ("zeroing") 

flightfull$co2pct <- shiftedco2 

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''


plot(flightfull$Seconds, flightfull$co2pct, type="l")
abline(v=fly1)
plot(flightfull$Seconds[fly], flightfull$co2pct[fly], type="l")
abline(v=fly1)



#''''''''''''''''''''''''''''''''''''''''''''''''''''''




## Add VCO2 to data frame and save complete corrected trace to file ----
flightfull$VCO2  <- flightfull$co2pct*flightfull$Flow_Rate*60 # mL CO2/hr

plot(flightfull$Seconds[fly], flightfull$VCO2[fly], main=c(filename, "VCO2"), type="l", xlab = "Seconds flight", ylab = "VCO2 flight")

## Estimate metabolic rate for specified regions and plot for confirmation


##### Subset just flight points -----

# This section creates a new column telling whether the datapoint starts a new behavior, ends it, or is the same as previous
wholetrace <- flightfull$Seconds[1]:length(flightfull$Seconds)
flightfull$markz <- "problem"
for(x in wholetrace) {
ifelse (flightfull$markers[x] != flightfull$markers[x-1], flightfull$markz[x] <- "new", 
                            ifelse (flightfull$markers[x] != flightfull$markers[x+1], flightfull$markz[x] <- "end", flightfull$markz[x] <- "same"))
}

# Selects all just flight
justflight <- data.frame(subset(flightfull, flightfull$markers == "F"))

(overallflight <- c(summary(justflight$VCO2), sd=round(sd(justflight$VCO2, na.rm=T))))

(maxflightVCO2 <- overallflight[6])
(meanflightVCO2 <- overallflight[4])
(sdflightVCO2 <- overallflight[7])

# longest continuous flight subset 
subnewend <- subset(justflight, justflight$markz != "same")
#View(subnewend) # tells all the start / ends of the behaviors
contigseconds <- data.frame(subnewend$Seconds[which(subnewend$markz == "new")], subnewend$Seconds[which(subnewend$markz == "end")])
names(contigseconds) <- c("new", "end")
contigseconds$length <- contigseconds$end - contigseconds$new
startseconds <- contigseconds$new[which(contigseconds$length == max(contigseconds$length))]
endseconds <- contigseconds$end[which(contigseconds$length == max(contigseconds$length))]

60 < endseconds - startseconds #Make sure this is true, that the longest flight is greater than 60 seconds

z <- subset(justflight, justflight$Seconds >= startseconds)
useflight <- subset(z, z$Seconds <= endseconds) # This is longest continuous flight 


##### Chosen Regions MR determination------

# 1. 
#     Overall contiguous  flight

(overallcontig <- c(summary(useflight$VCO2), sd=sd(useflight$VCO2, na.rm=T)))
(meancontig <- overallcontig[4])
(sdcontig <- overallcontig[7])
(maxcontig <- overallcontig[6])
(lengthcontig <- length(useflight$Seconds))



# 2. 
#     Most stable 1 minute of flight in longest contiguous stretch
allVar <- data.frame(startTime=integer(0), time=integer(0), varVal=numeric(0))
for(x in 1:length(useflight$Seconds)) {
  varNow <- var(useflight$VCO2[x:(x+60)])
  timenow <- useflight$Seconds[x]
  newRow <- c(x, timenow, varNow)
  allVar <- rbind(allVar, newRow)  
}

names(allVar) <- c("startTime", "time", "varVal")
View(allVar)
#now to find smallest variance
minVar <- min(allVar$varVal, na.rm=T)

(begin.sec <- allVar$startTime[which(allVar$varVal == minVar)])
duration.sec <- 60 
(MRlowvar <- c(summary(useflight$VCO2[begin.sec:(begin.sec+duration.sec)]), sd = sd(useflight$VCO2[begin.sec:(begin.sec+duration.sec)]), start=allVar$time[which(allVar$varVal == minVar)]))

(meanstable1min <- MRlowvar[4])
(sdstable1min <- MRlowvar[7])
(maxstable1min <- MRlowvar[6])




# 3. 
#     Highest mean 1 minute during contiguous flight
allhighmean <- data.frame(startTime=integer(0), timenow = integer(0), meanVal=numeric(0))

#This for loop computes the VCO2 mean for the next 1 minute for each timepoint 
for(x in 1:length(useflight$Seconds)){
  meanNow <- mean(useflight$VCO2[(x):(x+(60))], na.rm = T)
  timenow <- useflight$Seconds[x]
  newRow <- c(x, timenow, meanNow)
  allhighmean <- rbind(allhighmean, newRow)  
}
names(allhighmean)=c("startTime", "time", "meanVal")
#View(allhighmean)

#Get the estimates
beginsecmean <- allhighmean$startTime[which(allhighmean$meanVal == max(allhighmean$meanVal, na.rm = T))]
durationsecmean <- 60
highmean.region <- beginsecmean:(beginsecmean+durationsecmean)
(summaryhighmean <- c(summary(useflight$VCO2[highmean.region], na.rm=T), sd=sd(useflight$VCO2[highmean.region], na.rm=T), start=beginsecmean, end=beginsecmean+durationsecmean))
(meanhighregion <- summaryhighmean[4])
(sdhighregion <- summaryhighmean[7])
(maxhighregion<- summaryhighmean[6])




#4.   Highest 120 points during longest contiguous flight
n <- 120
decendingVCO2 <- sort(useflight$VCO2, decreasing = TRUE)
high120 <- decendingVCO2[1:n]
(regionhighest120 <- c(summary(high120), sd = sd(high120, na.rm=T)))
(meanhighest120 <- regionhighest120[4])
(sdhighest120 <- regionhighest120[7])
(maxhighest120 <- regionhighest120[6])



plot(useflight$Seconds, useflight$VCO2, type="l", main = c(filename, "continuous"))




#5.   Integral ... sort of .... 
#Sum  of all points in 
# The most stable 1 minute
(intregralstable1min <- sum(useflight$VCO2[begin.sec:(begin.sec+60)]))
 
# Highest mean 1 minute
(integralhigh1min <- sum(useflight$VCO2[beginsecmean:(beginsecmean+60)]))






###### Exporting graphs and numbers to excel docs etc. ----


### EXPORT HARD MR NUMBERS

# whole continguous 
meancontig <- overallcontig[4]
sdcontig <- overallcontig[7]
maxcontig <- overallcontig[6]
lengthcontig <- length(useflight$Seconds)

# stable 1 minutes
(meanstable1min <- MRlowvar[4])
(sdstable1min <- MRlowvar[7])
(maxstable1min <- MRlowvar[6])

# highest mean 1 minute
(meanhighregion <- summaryhighmean[4])
(sdhighregion <- summaryhighmean[7])
(maxhighregion<- summaryhighmean[6])

# highest 120 points
(meanhighest120 <- regionhighest120[4])
(sdhighest120 <- regionhighest120[7])
(maxhighest120 <- regionhighest120[6])



#Combine the metrics of RMR into one row and add to the a read in data file 
this_MR_data <- c(filename, meancontig, sdcontig, maxcontig, lengthcontig, meanstable1min, sdstable1min, maxstable1min,
                  meanhighregion, sdhighregion, maxhighregion, meanhighest120, sdhighest120, maxhighest120)


FMR_all <- read.csv("/Users/stephen/Desktop/CSU/Research/Analysis/MR_Analysis/MRAnalysis/FMR_all.csv", header = TRUE, sep=",", quote="''")
FMR_all$Filename <- as.character(FMR_all$Filename)
str(FMR_all)


#Add newest processed row to previous data, overwrite the csv (then move on to the next one)
FMR_all <- rbind(FMR_all, this_MR_data)
write.csv(FMR_all, "FMR_all.csv", row.names=F)









