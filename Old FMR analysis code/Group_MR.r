###  Analysis of metabolic rate, somewhat taken from Michael Dillon and Abbie Reade (MDAR) scripts...
# for analysis in R of Expedata files previously exported as .csv 
#  
#######  '''''    Overview & File key :::::::::: -------
#
#
#
# These are for Stephen's group level honeybee experiments 
#
#    
#
#
#    Seconds : seconds, starts with "0"
#    Barometric_Pressure :  pressure (kPa)
#    O2 :   Oxygen %, corrected for barometric pressure
#    Flowset : What is the flow setting (mL/min)
#    Flow_Rate : flow rate in system (mL/min)
#    CO2 :   CO2 %, corrected for barometric pressure  ***SM noted 10/25/18 :: This is a %, so it needs to divided by 100!!
#    Temp__Aux1_ : Temperature of incurrent air
#    Temp__Aux2_ : Temperature in the styrofoam where experimental chamber is located
#    FoxBoxTemp : Temperature in the Foxbox (Basically useless info)
#    Marker : marked points through the assay
#
#
#    Marker: -1 = no marker
#            118 = "v", switched to chamber
#            108 = "l", end / switch back to baseline 
#            86 = "V" switched to chamber
#            76 = "L" end
#    

# Output of interest is summarized CO2 trace - units will be ml CO2 




### Last Edited SM 11/5/19




#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#set the working directory, here, it is reinforcing the default project directory
#setwd("/Users/stephen/Desktop/CSU/Research/Groups2019/Analysis/Data Analysis/Groups2019Analysis") #unecessary if it's the same as where you started the project 
#getwd()
{
library(dplyr)
library(stringr)
library(foreach)
library(segmented)
library(grDevices)
library(tidyr)
library(zoo)


# These directories need to be set up in an organized folder
# directory for visual checks of data
dir.checks <- "MR-output/checktrace"
# directory where all corrected data is housed for further analysis
dir.processed <- "MR-output/processed"


# Functions
load.expedata.file <- function(){
  file <- file.choose()
  # extract information from filename
  filename <<- unlist(strsplit(basename(file), "[ ]"))[1]
  filename <<- gsub(".csv", "", tail(strsplit(file, "/")[[1]], n=1))
  #group_name <<- tail(strsplit(filename, "-")[[1]], n=1)
  #edate    <<- strsplit(filename, "_")[[1]][2]
  # read in data to "restfull" which will be full resting data
  restfull   <- read.csv(file, header=T)
  restfull$markers <- rawToChar(as.raw(as.character(restfull$Marker)), multiple=T) #turn markers in ascii character codes
  # this returns decimal fraction, uncorrected for water dilution because lags have to be fixed first
  #restfull$CO2uc <- restfull$CO2/1000000 
  #unique(restfull$chamber)
  return(restfull)
  #write.csv(restfull, "test.csv", row.names=F)
}



#######  '''''    Load file and transform data into correct units & correct for drift & lag ----- 


restfull <- load.expedata.file()
}
{
str(restfull)
restfull$filename <- filename
restfull$Group_name <- group_name

restfull$markers[1]<-"S"
restfull$markers[restfull$markers == ""] <- NA
restfull$markers <- na.locf(restfull$markers)

restfull$markers <- toupper(restfull$markers) # Turns to capital letters

chbr <- which(restfull$markers == "V") #In the sample chamber
bsln <- which(restfull$markers == "S") #Baseline at beginning
endpt <- which(restfull$markers == "L") #End, typically re-baselining

#first of each section
chbr1 <- chbr[1] #In the sample chamber
bsln1 <- bsln[1] #Baseline at beginning, first point
endpt1 <- endpt[1] #End, typically re-baselining


#'''''''''''''''''''''''''''''''''''''''''''''''''''
#Check the lags
#After the visual check -- mark the appropriate numbers down below for the lag (I think it's 10 for late and 4 for early SMR's)

plot(restfull$Seconds[(chbr1-20):(chbr1+20)], restfull$CO2[(chbr1-20):(chbr1+20)], main = "CO2 lag check")
abline(v=chbr1)

plot(restfull$Seconds[(chbr1-20):(chbr1+20)], restfull$O2[(chbr1-20):(chbr1+20)], main = "O2 lag check")
abline(v=chbr1)

#plot(restfull$Seconds[(chbr1-20):(chbr1+20)], restfull$Water_Vapor[(chbr1-20):(chbr1+20)], main = "WVP lag check")
#abline(v=chbr1)

plot(restfull$Seconds[(chbr1-20):(chbr1+20)], restfull$Barometric_Pressure[(chbr1-20):(chbr1+20)], main = "BP lag check")
abline(v=chbr1)

plot(restfull$Seconds[(chbr1-20):(chbr1+30)], restfull$Flow_Rate[(chbr1-20):(chbr1+30)], main = "FR lag check")
abline(v=chbr1)

#This creates a new column with "shifted" or lag-corrected readings to match the time point
# Manually change these numbers to fit the lag time
lagco2 <- 4 
lago2 <- 5 
#lagwvp <- 5 
lagpres <- lagco2
lagflow <- 5



restfull$co2nolag <- c(restfull$CO2[(1+lagco2):length(restfull$CO2)], rep(NA, lagco2))
#restfull$wvpnolag <- c(restfull$Water_Vapor[(1+lagwvp):length(restfull$Water_Vapor)], rep(NA, lagwvp))
restfull$o2nolag <- c(restfull$O2[(1+lago2):length(restfull$O2)], rep(NA, lago2))
restfull$flownolag <- c(restfull$Flow_Rate[(1+lagflow):length(restfull$Flow_Rate)], rep(NA, lagflow))
restfull$presnolag <- c(restfull$Barometric_Pressure[(1+lagpres):length(restfull$Barometric_Pressure)], rep(NA, lagpres))

#plot full data lag corrected
#plot(restfull$Seconds, restfull$wvpnolag, type="l", col="blue", xlab = "Seconds", main = "No Lag, uncorrected for drift")
#par(new=T)
plot(restfull$Seconds, restfull$co2nolag, type="l", col="darkorange", axes=F, ann=F)
par(new=T)
plot(restfull$Seconds, restfull$o2nolag, type="l", col="green", axes=F, ann=F)
par(new=T)
legend(x="topright", legend=c("CO2", "O2"), col=c("darkorange","green"), bty = "n", lty=1, cex=.75)
par(new=T)
plot(restfull$Seconds, restfull$flownolag, type="l", axes=T, ann=F)
abline(v=chbr1-lagco2, col = "black")
abline(v=chbr1+120, col="red")




#######  '''''    Correct CO2 for WVP dilution ----
#  units in kPa
# correction factor: BP/(BP-WVP), where
# BP = barometric pressure, column was measured! 
# WVP = measured "H2O" column

#restfull$co2_wvpc <- restfull$co2nolag * restfull$presnolag/(restfull$presnolag-restfull$wvpnolag)
# This is unnecessary, see page 96 of Lighton 2008 
restfull$co2_wvpc <- restfull$co2nolag

# Units still % CO2

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





#######  '''''    Stable Baseline determination for zero correct-----
#' This gets the most stable 30 seconds for the baseline calculation
#' 
allVarBase <- data.frame(startTime=integer(0), varValBase=numeric(0))
for(x in bsln1:(chbr1-5)){ 
  varNowBase <- var(restfull$co2_wvpc[x:(x+30)], na.rm = T)
  newRowBase <- c(x, varNowBase)
  allVarBase <- rbind(allVarBase, newRowBase)  
}
names(allVarBase)=c("startTime", "varValBase")
#View(allVar)
#now to find minumum variance point using the above time frame in the for loop (2 minutes) 
minVarBase <- min(allVarBase$varValBase)

#Get the time estimates for lowest variance
begin.base <- allVarBase$startTime[which(allVarBase$varValBase==min(allVarBase$varValBase))]
duration.base <- 30
lowVar.baseline <- begin.base:(begin.base+duration.base)

prelowvar <- mean(restfull$co2_wvpc[lowVar.baseline], na.rm=T)
FiCO2 <- prelowvar

prelowvarO2 <- mean(restfull$o2nolag[lowVar.baseline], na.rm = T)
FiO2 <- prelowvarO2

# drift correction
#This is subtracting the average %C of the baseline (lowest %C variance from 30 s) from every point during the sample

driftcorrectO2 <- prelowvarO2 - restfull$o2nolag
# subtracting out baseline ("zeroing") 

restfull$co2pct <- restfull$co2_wvpc - FiCO2 # FeCO2 - FiCO2

# The below calculation gives the oxygen CONSUMED
restfull$o2pct <- FiO2 - restfull$o2nolag # FiO2 - FeO2

# units still % CO2
# units still % O2


# PLOT the drift correction in R
plot(restfull$Seconds[], restfull$co2pct[], type="l", ylab= "CO2 lag and drift corrected", xlab="Seconds", main=c(filename, " CO2 trace"))
par(new=T)
plot(restfull$Seconds[], restfull$o2pct[], type = "l", col= "red", ylab = "O2 lag and drift corrected", xlab = "Seconds", main=c(filename, " O2 trace"))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''



#######  '''''    Add VCO2 to data frame and calculate output values ---------

# Get VCO2 by multiplying by the flow rate (in ml/min)
restfull$VCO2  <- restfull$co2pct / 100 * restfull$flownolag * 60 # mL CO2/hr  ### on 10/25/18 SM discovered the mistake that multiplying by a percent is stupid, and you need to divide out 100 ... Now (7/8/19) fixed here
# units are NOW :: ml CO2 / hr and are fully corrected 
# essentially equation 10.3 from Lighton 2008 book

restfull$VO2 <- restfull$o2pct / 100 * restfull$flownolag * 60 # mL O2 / hr ### may be wrong w/o H2O correction ?? (10/30/19 SM)

restfull$RQV <- restfull$VCO2/restfull$VO2

usetime <- (chbr1+120):(chbr1+(60*25)) # make "usetime" 2 min to 10 min
  
plot(restfull$Seconds[usetime], restfull$VCO2[usetime], main=c(filename, "VCO2 10 minute trace"), type="l", xlab = "Seconds", ylab = "VCO2 ml/hr")

## Estimate metabolic rate for specific regions
#     *****Currently (7/23) I don't make any O2 corrections, and everything is based off of VCO2 without WVP correction **********


# 0. Full chamber time trace, using the "usetime" determined above
(overallVCO2 <- c(summary(restfull$VCO2[usetime]), sd=sd(restfull$VCO2[usetime], na.rm=T)))
(meanoverallVCO2 <- overallVCO2[4])
(sdoverallVCO2 <- overallVCO2[8])

(overallVO2 <- c(summary(restfull$VO2[usetime]), sd = sd(restfull$VO2[usetime], na.rm = T)))
(meanoverallVO2 <- overallVO2[4])
(sdoverallVO2 <- overallVO2[8])

(overallRQV <- c(summary(restfull$RQV[usetime]), sd = sd(restfull$RQV[usetime], na.rm = T)))
(meanoverallRQV <- overallRQV[4])
(sdoverallRQV <- overallRQV[8])


#   1.    LOWEST VARIANCE 2 MINUTES

allVar <- data.frame(startTime=integer(0), varVal=numeric(0))
beg.r <- chbr1 + 120 # 2 minutes after chamber switch
end.r <- ifelse(is.na(endpt1) == T, (chbr1 + 60*8), (endpt1-120)) # 10 minutes after chamber switch ;; beg.r : end.r = 8 minutes of data

#This for loop computes the VCO2 variance for the following two minutes for each timepoint 
for(x in beg.r : end.r){
  varNow <- var(restfull$VCO2[(x):(x+(60*2))], na.rm = T)
  newRow <- c(x, varNow)
  allVar <- rbind(allVar, newRow)  
}
names(allVar)=c("startTime", "varVal")
#View(allVar)
#now to find minumum variance point using the above time frame in the for loop (2 minutes) 
minVar=min(allVar$varVal)

#Get the estimates 
begin.seclowvar <- allVar$startTime[which(allVar$varVal==min(allVar$varVal))]
duration.seclowvar <- 60*2 
lowVar.region <- begin.seclowvar:(begin.seclowvar+duration.seclowvar)
(regionlowvar <- c(summary(restfull$VCO2[lowVar.region], na.rm=T), sd=sd(restfull$VCO2[lowVar.region], na.rm=T), start=begin.seclowvar, end=begin.seclowvar+duration.seclowvar))
(meanstable2min <- regionlowvar[4])
(sdstable2min <- regionlowvar[7])

(regionlowvarO2 <- c(summary(restfull$VO2[lowVar.region], na.rm=T), sd=sd(restfull$VO2[lowVar.region], na.rm=T), start=begin.seclowvar, end=begin.seclowvar+duration.seclowvar))
(meanstable2minO2 <- regionlowvarO2[4])
(sdstable2minO2 <- regionlowvarO2[7])

(regionlowvarRQV <- c(summary(restfull$RQV[lowVar.region], na.rm = T), sd = sd(restfull$RQV[lowVar.region], na.rm = T), start = begin.seclowvar, end = begin.seclowvar+duration.seclowvar))
(meanstable2minRQV <- regionlowvarRQV[4])
(sdstable2minRQV <- regionlowvarRQV[7])


################.
#   2.     LOWEST 2 minute mean
allMean <- data.frame(startTime=integer(0), meanVal=numeric(0))

#This for loop computes the VCO2 mean for the next two minutes for each timepoint 
for(x in beg.r:end.r){
  meanNow <- mean(restfull$VCO2[(x):(x+(60*2))], na.rm = T)
  newRow <- c(x, meanNow)
  allMean <- rbind(allMean, newRow)  
}
names(allMean)=c("startTime", "meanVal")
#View(allMean)
#now to find minumum mean point using the above time frame in the for loop (2 minutes) 
minMean=min(allMean$meanVal)

#Get the estimates
begin.seclowmean <- allMean$startTime[which(allMean$meanVal==min(allMean$meanVal))]
duration.seclowmean <- 60*2
lowMean.region <- begin.seclowmean:(begin.seclowmean+duration.seclowmean)
(regionlowmean <- c(summary(restfull$VCO2[lowMean.region], na.rm=T), sd=sd(restfull$VCO2[lowMean.region], na.rm=T), start=begin.seclowmean, end=begin.seclowmean+duration.seclowmean))
(meanlowregion <- regionlowmean[4])
(sdlowregion <- regionlowmean[7])

(regionlowmeanO2 <- c(summary(restfull$VO2[lowMean.region], na.rm=T), sd=sd(restfull$VO2[lowMean.region], na.rm=T), start=begin.seclowmean, end=begin.seclowmean+duration.seclowmean))
(meanlowregionO2 <- regionlowmeanO2[4])
(sdlowregionO2 <- regionlowmeanO2[7])

(regionlowmeanRQV<- c(summary(restfull$RQV[lowMean.region], na.rm=T), sd=sd(restfull$RQV[lowMean.region], na.rm=T), start=begin.seclowmean, end=begin.seclowmean+duration.seclowmean))
(meanlowregionRQ <- regionlowmeanRQV[4])
(sdlowregionRQ <- regionlowmeanRQV[7])



################.
#     3.     Lowest 120 seconds within the chamber, taken from the 'usetime' defined above
n <- 120

decVCO2 <- sort(restfull$VCO2[usetime])
low120 <- decVCO2[1:n]
(regionlow120 <- c(summary(low120), sd = sd(low120, na.rm=T)))
(meanlow120 <- regionlow120[4])
(sdlow120 <- regionlow120[7])

decVO2 <- sort(restfull$VO2[usetime])
low120O2 <- decVO2[1:n]
(regionlow120O2 <- c(summary(low120O2), sd = sd(low120O2, na.rm=T)))
(meanlow120O2 <- regionlow120O2[4])
(sdlow120O2 <- regionlow120O2[7])




################.
#     4.     Mean Peak  

#This makes 2 new columns, and the peaks will have BOTH "increase" AND "p" written
restfull$peak1 <- 1
restfull$peakbig <- 1
for(x in beg.r:end.r) {
  restfull$peak1[x] <- ifelse(restfull$VCO2[x] == max(restfull$VCO2[(x-10):(x+10)]), "peak", "notpeak") #****MAJOR problem, is that this assumes a period of > 20 seconds, which is not always true!
  restfull$peakbig[x] <- ifelse(restfull$VCO2[x] >= (.5*max(restfull$VCO2[beg.r:end.r])), "high", "nothigh")
  }
subpeak <- subset(restfull, peak1 == "peak")
plot(restfull$Seconds, restfull$VCO2, type="l", main = c(filename," CO2 trace peaks"), ylab="VCO2 ml/hr", xlab= "Seconds")
abline(v= subpeak$Seconds, col="red")

# Biggest peaks ( > 50% of the max)
subbigpeak <- subset(subpeak, peakbig == "high")
plot(restfull$Seconds, restfull$VCO2, type="l", main = c(filename," CO2 trace biggest peaks"), ylab="VCO2 ml/hr", xlab= "Seconds")
abline(v= subbigpeak$Seconds, col="blue")

(summarypeakVCO2 <- c(summary(subpeak$VCO2), sd = sd(subpeak$VCO2, na.rm=T)))
(peakmean <- summarypeakVCO2[4])
(peaksd <- summarypeakVCO2[7])




################.
#     5.    Mean valley
restfull$valley1 <- 1
for(x in beg.r:end.r) {
  restfull$valley1[x] <- ifelse(restfull$VCO2[x] == min(restfull$VCO2[(x-10):(x+10)]), "valley", "nonvalley") #****MAJOR problem, is that this assumes a period of > 20 seconds, which is not always true!
}
subvalley <- subset(restfull, valley1 == "valley")
plot(restfull$Seconds, restfull$VCO2, type="l", main = c(filename," CO2 trace valleys"), ylab="VCO2 ml/hr", xlab= "Seconds")
abline(v= subvalley$Seconds, col="blue")

(summaryvalleyVCO2 <- c(summary(subvalley$VCO2), sd = sd(subvalley$VCO2, na.rm = T)))
(valleymean <- summaryvalleyVCO2[4])
(valleysd <- summaryvalleyVCO2[7])



################.
#      6.         Gets the PERIOD by calculating the seconds between each peak 
subpeak$per1 <- c(1:length(subpeak$Seconds))
for(i in subpeak$per1){
  subpeak$per1[i] <- (subpeak$Seconds[i+1]-subpeak$Seconds[i])
}
(peakperiodmean <- mean(subpeak$per1, na.rm=T))
(peakperiodsd <- sd(subpeak$per1, na.rm=T))
(summarypeaks <- c(summary(subpeak$per1), sd = sd(subpeak$per1, na.rm=T)))
(peaknums <- length(subpeak$VCO2))





################.
#         7.  Most stable 4 minutes

allVar <- data.frame(startTime=integer(0), varVal=numeric(0))
beg.r <- chbr1 + 240 # 4 minutes after chamber switch
end.r <- ifelse(is.na(endpt1) == T, (chbr1 + 60*8), (endpt1-240)) # 10 minutes after chamber switch ;; beg.r : end.r = 8 minutes of data

#This for loop computes the VCO2 variance for the following two minutes for each timepoint 
for(x in beg.r : end.r){
  varNow <- var(restfull$VCO2[(x):(x+(60*2))], na.rm = T)
  newRow <- c(x, varNow)
  allVar <- rbind(allVar, newRow)  
}
names(allVar)=c("startTime", "varVal")
#View(allVar)
#now to find minumum variance point using the above time frame in the for loop (2 minutes) 
minVar=min(allVar$varVal)

#Get the estimates 
begin.seclowvar4 <- allVar$startTime[which(allVar$varVal==min(allVar$varVal))]
duration.seclowvar4 <- 60*4 
lowVar.region4 <- begin.seclowvar4:(begin.seclowvar4+duration.seclowvar4)
(region4lowvar <- c(summary(restfull$VCO2[lowVar.region4], na.rm=T), sd=sd(restfull$VCO2[lowVar.region4], na.rm=T), start=begin.seclowvar4, end=begin.seclowvar4+duration.seclowvar4))
(meanstable4min <- region4lowvar[4])
(sdstable4min <- region4lowvar[7])


(region4lowvarO2 <- c(summary(restfull$VO2[lowVar.region4], na.rm=T), sd=sd(restfull$VO2[lowVar.region4], na.rm=T), start=begin.seclowvar4, end=begin.seclowvar4+duration.seclowvar4))
(meanstable4minO2 <- region4lowvarO2[4])
(sdstable4minO2 <- region4lowvarO2[7])


#         8.        Lowest 4 minutes

allMean <- data.frame(startTime=integer(0), meanVal=numeric(0))

#This for loop computes the VCO2 mean for the next two minutes for each timepoint 
for(x in beg.r:end.r){
  meanNow <- mean(restfull$VCO2[(x):(x+(60*4))], na.rm = T)
  newRow <- c(x, meanNow)
  allMean <- rbind(allMean, newRow)  
}
names(allMean)=c("startTime", "meanVal")
#View(allMean)
#now to find minumum mean point using the above time frame in the for loop (2 minutes) 
minMean=min(allMean$meanVal)

#Get the estimates
begin.seclowmean4 <- allMean$startTime[which(allMean$meanVal==min(allMean$meanVal))]
duration.seclowmean4 <- 60*4
lowMean.region4 <- begin.seclowmean4:(begin.seclowmean4+duration.seclowmean4)
(region4lowmean <- c(summary(restfull$VCO2[lowMean.region4], na.rm=T), sd=sd(restfull$VCO2[lowMean.region4], na.rm=T), start=begin.seclowmean4, end=begin.seclowmean4+duration.seclowmean4))
(mean4lowregion <- region4lowmean[4])
(sd4lowregion <- region4lowmean[7])

(region4lowmeanO2 <- c(summary(restfull$VO2[lowMean.region4], na.rm=T), sd=sd(restfull$VO2[lowMean.region4], na.rm=T), start=begin.seclowmean4, end=begin.seclowmean4+duration.seclowmean4))
(meanlow4regionO2 <- region4lowmeanO2[4])
(sdlow4regionO2 <- region4lowmeanO2[7])


#######  '''''    Temperature values -----

#TEMP
# 2 min low var
(Tempincurlowvar2 <- c(summary(restfull$Temp__Aux1_[lowVar.region], na.rm=T), sd=sd(restfull$Temp__Aux1_[lowVar.region], na.rm=T), start=begin.seclowvar, end=begin.seclowvar+duration.seclowvar))
(meanTincurlowvar2 <- Tempincurlowvar2[4])
(Tempchbrlowvar2 <- c(summary(restfull$Temp__Aux2_[lowVar.region], na.rm=T), sd=sd(restfull$Temp__Aux2_[lowVar.region], na.rm=T), start=begin.seclowvar, end=begin.seclowvar+duration.seclowvar))
(meanTchbrlowvar2 <- Tempchbrlowvar2[4])

# 4 min low var
(Tempincurlowvar4 <- c(summary(restfull$Temp__Aux1_[lowVar.region4], na.rm=T), sd=sd(restfull$Temp__Aux1_[lowVar.region4], na.rm=T), start=begin.seclowvar4, end=begin.seclowvar4+duration.seclowvar4))
(meanTincurlowvar4 <- Tempincurlowvar4[4])
(Tempchbrlowvar4 <- c(summary(restfull$Temp__Aux2_[lowVar.region4], na.rm=T), sd=sd(restfull$Temp__Aux2_[lowVar.region4], na.rm=T), start=begin.seclowvar4, end=begin.seclowvar4+duration.seclowvar4))
(meanTchbrlowvar4 <- Tempchbrlowvar4[4])

# 2 min low mean
(Tempincurlowmean2 <- c(summary(restfull$Temp__Aux1_[lowMean.region], na.rm=T), sd=sd(restfull$Temp__Aux1_[lowMean.region], na.rm=T), start=begin.seclowmean, end=begin.seclowmean+duration.seclowmean))
(meanTincurlowmean2 <- Tempincurlowmean2[4])
(Tempchbrlowmean2 <- c(summary(restfull$Temp__Aux2_[lowMean.region], na.rm=T), sd=sd(restfull$Temp__Aux2_[lowMean.region], na.rm=T), start=begin.seclowmean, end=begin.seclowmean+duration.seclowmean))
(meanTchbrlowmean2 <- Tempchbrlowmean2[4])

# 4 min low mean
(Tempincurlowmean4 <- c(summary(restfull$Temp__Aux1_[lowMean.region4], na.rm=T), sd=sd(restfull$Temp__Aux1_[lowMean.region4], na.rm=T), start=begin.seclowmean4, end=begin.seclowmean4+duration.seclowmean4))
(meanTincurlowmean4 <- Tempincurlowmean4[4])
(Tempchbrlowmean4 <- c(summary(restfull$Temp__Aux2_[lowMean.region4], na.rm=T), sd=sd(restfull$Temp__Aux2_[lowMean.region4], na.rm=T), start=begin.seclowmean4, end=begin.seclowmean4+duration.seclowmean4))
(meanTchbrlowmean4 <- Tempchbrlowmean4[4])



#######  '''''    VdotCO2 and VdotO2 for all the same times ------
# *** Divide Fe and Fi by 100 first
restfull$FeCO2 <- restfull$co2_wvpc / 100
restfull$FiCO2 <- prelowvar / 100
  
restfull$FiO2 <- prelowvarO2 / 100
restfull$FeO2 <- restfull$o2nolag /100

# Now units are in fractional concentrations (not %'s)

restfull$VdotCO2 <- restfull$flownolag * 60 * (restfull$FeCO2 - restfull$FiCO2 + restfull$FiCO2*(restfull$FiO2 - restfull$FeO2)) / (1 + restfull$FiCO2)  # Equations in chapter 11 
# output in ml CO2 / hr

restfull$VdotO2 <- restfull$flownolag * 60 * (restfull$FiO2 - restfull$FeO2 - restfull$FiO2*(restfull$FeCO2 - restfull$FiCO2)) / (1 - restfull$FiO2)
# output in ml O2 / hr

restfull$RQVdot <- restfull$VdotCO2 / restfull$VdotO2

#######  '''''    VdotCO2 and VdotO2 summary output values -------

### 1. overall
(overallVdotCO2 <- c(summary(restfull$VdotCO2[usetime]), sd=sd(restfull$VdotCO2[usetime], na.rm=T)))
(meanoverallVdotCO2 <- overallVdotCO2[4])
(sdoverallVdotCO2 <- overallVdotCO2[8])

(overallVdotO2 <- c(summary(restfull$VdotO2[usetime]), sd = sd(restfull$VdotO2[usetime], na.rm = T)))
(meanoverallVdotO2 <- overallVdotO2[4])
(sdoverallVdotO2 <- overallVdotO2[8])

(overallRQdot <- c(summary(restfull$RQVdot[usetime]), sd = sd(restfull$RQVdot[usetime], na.rm = T)))
(meanoverallRQdot <- overallRQdot[4])
(sdoverallRQdot <- overallRQdot[8])


### 2. lowest CO2 variance 2 min

alldotVar <- data.frame(startTime=integer(0), varVal=numeric(0))
beg.r <- chbr1 + 120 # 2 minutes after chamber switch
end.r <- ifelse(is.na(endpt1) == T, (chbr1 + 60*8), (endpt1-120)) # 10 minutes after chamber switch ;; beg.r : end.r = 8 minutes of data

for(x in beg.r : end.r){
  varNow <- var(restfull$VdotCO2[(x):(x+(60*2))], na.rm = T)
  newRow <- c(x, varNow)
  alldotVar <- rbind(alldotVar, newRow)  
}
names(alldotVar)=c("startTime", "varVal")

#Get the estimates 
begin..seclowvar <- alldotVar$startTime[which(alldotVar$varVal==min(alldotVar$varVal))]
duration..seclowvar <- 60*2 
lowVar..region <- begin..seclowvar:(begin..seclowvar+duration..seclowvar)
(region.lowvar <- c(summary(restfull$VdotCO2[lowVar..region], na.rm=T), sd=sd(restfull$VdotCO2[lowVar..region], na.rm = T), start = begin..seclowvar, end = begin..seclowvar+duration..seclowvar))
(mean.stable2min <- region.lowvar[4])
(sd.stable2min <- region.lowvar[7])

(region.lowvarO2 <- c(summary(restfull$VdotO2[lowVar..region], na.rm=T), sd=sd(restfull$VdotO2[lowVar..region], na.rm = T), start = begin..seclowvar, end = begin..seclowvar+duration..seclowvar))
(mean.stable2minO2 <- region.lowvarO2[4])
(sd.stable2minO2 <- region.lowvarO2[7])

(region.lowvarRQdot <- c(summary(restfull$RQVdot[lowVar..region], na.rm=T), sd=sd(restfull$RQVdot[lowVar..region], na.rm = T), start = begin..seclowvar, end = begin..seclowvar+duration..seclowvar))
(mean.stable2minRQdot <- region.lowvarRQdot[4])
(sd.stable2minRQdot <- region.lowvarRQdot[7])


### 3. lowest CO2 mean 2 min

all.Mean <- data.frame(startTime=integer(0), meanVal=numeric(0))

#This for loop computes the VCO2 mean for the next two minutes for each timepoint 
for(x in beg.r:end.r){
  meanNow <- mean(restfull$VdotCO2[(x):(x+(60*2))], na.rm = T)
  newRow <- c(x, meanNow)
  all.Mean <- rbind(all.Mean, newRow)  
}
names(all.Mean)=c("startTime", "meanVal")

#Get the estimates
begin..seclowmean <- all.Mean$startTime[which(all.Mean$meanVal == min(all.Mean$meanVal))]
duration..seclowmean <- 60*2
lowMean..region <- begin..seclowmean:(begin..seclowmean + duration..seclowmean)
(region.lowmean <- c(summary(restfull$VdotCO2[lowMean..region], na.rm=T), sd=sd(restfull$VdotCO2[lowMean..region], na.rm=T), start=begin..seclowmean, end=begin..seclowmean+duration..seclowmean))
(mean.lowregion <- region.lowmean[4])
(sd.lowregion <- region.lowmean[7])

(region.lowmeanO2 <- c(summary(restfull$VdotO2[lowMean..region], na.rm=T), sd=sd(restfull$VdotO2[lowMean..region], na.rm=T), start=begin..seclowmean, end=begin..seclowmean+duration..seclowmean))
(mean.lowregionO2 <- region.lowmeanO2[4])
(sd.lowregionO2 <- region.lowmeanO2[7])

(region.lowmeanRQdot <- c(summary(restfull$RQVdot[lowMean..region], na.rm=T), sd=sd(restfull$RQVdot[lowMean..region], na.rm=T), start=begin..seclowmean, end=begin..seclowmean+duration..seclowmean))
(mean.lowregionRQdot <- region.lowmeanRQdot[4])
(sd.lowregionRQdot <- region.lowmeanRQdot[7])



#######  '''''    Export Images of Traces -------------

# Export graphs of whole trace + each determined region + peaks etc. 
analysis.file <- paste(filename, "traces.pdf", sep="-")
pdf(paste(dir.processed, analysis.file, sep="/"), width=10, height=7)

#plot(restfull$Seconds, restfull$wvpnolag, type="l", col="blue", xlab = "Seconds", main = paste(filename, "- No Lag, uncorrected for drift"))
#par(new=T)
plot(restfull$Seconds, restfull$co2nolag, type="l", col="darkorange", axes=F, ann=F)
par(new=T)
plot(restfull$Seconds, restfull$o2nolag, type="l", col="green", axes=F, ann=F)
par(new=T)
legend(x="topright", legend=c("WVP", "CO2", "O2"), col=c("blue", "darkorange","green"), bty = "n", lty=1, cex=.75)

plot(restfull$Seconds[], restfull$co2pct[], type="l", ylab= "CO2 lag and drift corrected", xlab="Seconds", main=c(filename, " CO2 trace"))


plot(restfull$Seconds, restfull$VCO2, type="l", main = c(filename," CO2 trace peaks"), ylab="VCO2 ml/hr", xlab= "Seconds")
abline(v= subpeak$Seconds, col="red")

plot(restfull$Seconds, restfull$VCO2, type="l", main = c(filename," CO2 trace valleys"), ylab="VCO2 ml/hr", xlab= "Seconds")
abline(v= subvalley$Seconds, col="blue")

plot(restfull$Seconds[usetime], restfull$VCO2[usetime], main=c(filename, "VCO2 10 mins"), type="l", xlab = "Seconds", ylab = "VCO2 ml/hr")

plot(restfull$Seconds[lowMean.region], restfull$VCO2[lowMean.region], type = "l", xlab=paste("Seconds (start:", begin.seclowmean,"s)"), ylab= "VCO2 ml/hr", main = paste(filename, "Low 2-min mean"))
plot(restfull$Seconds[lowVar.region], restfull$VCO2[lowVar.region], type = "l", xlab=paste("Seconds (start:",begin.seclowvar,"s)"), ylab= "VCO2 ml/hr", main = paste(filename, "Low 2-min var,"))

plot(restfull$Seconds[lowMean.region], restfull$VO2[lowMean.region], type = "l", xlab=paste("Seconds (start:", begin.seclowmean,"s)"), ylab= "VO2 ml/hr", main = paste(filename, "Low 2-min mean"))
plot(restfull$Seconds[lowVar.region], restfull$VO2[lowVar.region], type = "l", xlab=paste("Seconds (start:",begin.seclowvar,"s)"), ylab= "VO2 ml/hr", main = paste(filename, "Low 2-min var,"))

plot(restfull$Seconds[lowMean.region], restfull$RQV[lowMean.region], type = "l", xlab=paste("Seconds (start:", begin.seclowmean,"s)"), ylab= "RQ", main = paste(filename, "Low 2-min mean"))
plot(restfull$Seconds[lowVar.region], restfull$RQV[lowVar.region], type = "l", xlab=paste("Seconds (start:",begin.seclowvar,"s)"), ylab= "RQ", main = paste(filename, "Low 2-min var,"))

plot(restfull$Seconds[lowMean..region], restfull$VdotCO2[lowMean..region], type = "l", xlab=paste("Seconds (start:", begin..seclowmean,"s)"), ylab= "V.CO2 ml/hr", main = paste(filename, "Low 2-min mean"))
plot(restfull$Seconds[lowVar..region], restfull$VdotCO2[lowVar..region], type = "l", xlab=paste("Seconds (start:",begin..seclowvar,"s)"), ylab= "V.CO2 ml/hr", main = paste(filename, "Low 2-min var,"))

plot(restfull$Seconds[lowMean..region], restfull$VdotO2[lowMean..region], type = "l", xlab=paste("Seconds (start:", begin..seclowmean,"s)"), ylab= "V.O2 ml/hr", main = paste(filename, "Low 2-min mean"))
plot(restfull$Seconds[lowVar..region], restfull$VdotO2[lowVar..region], type = "l", xlab=paste("Seconds (start:",begin..seclowvar,"s)"), ylab= "V.O2 ml/hr", main = paste(filename, "Low 2-min var,"))

plot(restfull$Seconds[lowMean..region], restfull$RQV[lowMean..region], type = "l", xlab=paste("Seconds (start:", begin..seclowmean,"s)"), ylab= "RQ", main = paste(filename, "Low 2-min mean"))
plot(restfull$Seconds[lowVar..region], restfull$RQV[lowVar..region], type = "l", xlab=paste("Seconds (start:",begin..seclowvar,"s)"), ylab= "RQ", main = paste(filename, "Low 2-min var,"))

plot(restfull$Seconds[lowMean.region4], restfull$VCO2[lowMean.region4], type = "l", xlab=paste("Seconds (start:", begin.seclowmean4,"s)"), ylab= "VCO2 ml/hr", main = paste(filename, "Low 4-min mean"))
plot(restfull$Seconds[lowVar.region4], restfull$VCO2[lowVar.region4], type = "l", xlab=paste("Seconds (start:",begin.seclowvar4,"s)"), ylab= "VCO2 ml/hr", main = paste(filename, "Low 4-min var,"))

dev.off()



##### Export summary values -------

#Combine the metrics of RMR into one row and add to the a read in data file 
this_MR_data <- data.frame(group_name, filename, 
                   meanoverallVCO2, sdoverallVCO2, 
                   meanoverallVO2, sdoverallVO2,     
                      meanstable2min, sdstable2min, meanstable4min, sdstable4min, 
                      meanstable2minO2, sdstable2minO2, meanstable4minO2, sdstable4minO2,
                        meanlowregion, sdlowregion, mean4lowregion, sd4lowregion,
                        meanlowregionO2, sdlowregionO2, meanlow4regionO2, sdlow4regionO2,
                          meanlow120, sdlow120,
                          meanlow120O2, sdlow120O2,
                            peakmean, peaksd, valleymean, valleysd, peakperiodmean, peakperiodsd, peaknums,
                              meanoverallVdotCO2, sdoverallVdotCO2, meanoverallVdotO2, sdoverallVdotO2,
                              mean.stable2min, sd.stable2min, mean.stable2minO2, sd.stable2minO2,
                              mean.lowregion, sd.lowregion, mean.lowregionO2, sd.lowregionO2,
                                meanTincurlowvar2, meanTchbrlowvar2, meanTincurlowvar4, meanTchbrlowvar4, 
                                meanTincurlowmean2, meanTchbrlowmean2, meanTincurlowmean4, meanTchbrlowmean4,
                                  meanoverallRQV, sdoverallRQV, meanstable2minRQV, sdstable2minRQV, meanlowregionRQ, sdlowregionRQ,
                                    meanoverallRQdot, sdoverallRQdot, mean.stable2minRQdot, sd.stable2minRQdot, mean.lowregionRQdot, sd.lowregionRQdot
                              )

names(this_MR_data) <- c("group_name", "filename", 
                         "meanoverallVCO2", "sdoverallVCO2", 
                         "meanoverallVO2", "sdoverallVO2",     
                         "meanstable2min", "sdstable2min", "meanstable4min", "sdstable4min", 
                         "meanstable2minO2", "sdstable2minO2", "meanstable4minO2", "sdstable4minO2",
                         "meanlowregion", "sdlowregion", "mean4lowregion", "sd4lowregion",
                         "meanlowregionO2", "sdlowregionO2", "meanlow4regionO2", "sdlow4regionO2",
                         "meanlow120", "sdlow120",
                         "meanlow120O2", "sdlow120O2",
                         "peakmean", "peaksd", "valleymean", "valleysd", "peakperiodmean", "peakperiodsd", "peaknums",
                         "meanoverallVdotCO2", "sdoverallVdotCO2", "meanoverallVdotO2", "sdoverallVdotO2",
                         "mean.stable2min", "sd.stable2min", "mean.stable2minO2", "sd.stable2minO2",
                         "mean.lowregion", "sd.lowregion", "mean.lowregionO2", "sd.lowregionO2",
                         "meanTincurlowvar2", "meanTchbrlowvar2", "meanTincurlowvar4", "meanTchbrlowvar4", 
                         "meanTincurlowmean2", "meanTchbrlowmean2", "meanTincurlowmean4", "meanTchbrlowmean4",
                         "meanoverallRQV", "sdoverallRQV", "meanstable2minRQV", "sdstable2minRQV", "meanlowregionRQ", "sdlowregionRQ",
                         "meanoverallRQdot", "sdoverallRQdot", "mean.stable2minRQdot", "sd.stable2minRQdot", "mean.lowregionRQdot", "sd.lowregionRQdot"
)

#View(this_MR_data)
}
RMR_all <- read.csv("Group_MR_all.csv", header = T, sep=",") # Change this to the desired file name for your application
RMR_all$filename <- as.character(RMR_all$filename)
str(RMR_all)

#Add newest processed row to previous data, overwrite the csv (then move on to the next one)
RMR_all <- rbind(RMR_all, this_MR_data)
write.csv(RMR_all, "Group_MR_all.csv", row.names=F)


# check then clear the workspace for the next processing
#ls()
rm(list = ls())

