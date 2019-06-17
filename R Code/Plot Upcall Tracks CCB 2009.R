# Clear workspace and memory
rm(list=ls())
gc()

# Load packages
library(ggplot2)
library(dplyr)
library(gridExtra)

# Set working directory
setwd("/home/kpalmer/Desktop/CCB_Comparison/SECR")
# setwd("~/Dropbox/Kaitlin")

#==============================================
#   READ, CLEAN, SUBSET DATA
#==============================================


#________________________________
# Load ranges, capture history, and MARU locations

# File without times
CCBData_orig <- read.table(file='CCBData_export_October_5.csv', sep=',', header=T)

# CCBData_orig<-read.table(file='/home/kpalmer/Desktop/CCB_Comparison/SECR/C_Oct162018.txt',
#                          sep=',', header=T)
#CCBData_orig <- read.table(file='CCB_data_locs.csv', sep=',', header=T, na.strings = "NaN")
str(CCBData_orig)

# Load MARU locations
MARUloc <- read.table("Mass Bay MARU_locs.txt",sep="",header=T)



#________________________________
# Format data
# Retain date for easy plotting
CCBData_orig$UselessDate = CCBData_orig$date
CCBData_orig$Channel <- as.factor(CCBData_orig$Channel)

# Create date formats
mdy <- paste("2009",substr(CCBData_orig$date,1,1), substr(CCBData_orig$date,2,3), sep="-")
CCBData_orig$date <- as.POSIXct(mdy, tz="UTC") + CCBData_orig$StartTime
CCBData_orig$datum <- factor(format(CCBData_orig$date, "%Y-%m-%d"))
CCBData_orig$ndate <- as.numeric(CCBData_orig$date)
rm(mdy)

# ***Correct the Lon-Lat labels***
names(CCBData_orig)[names(CCBData_orig) %in% c("Lon","Lat")] <- c("Lon","Lat")

# Clean Captures
# Some Captures have value 2 or 4! what does it mean?
CCBData_orig[which(CCBData_orig$Capture > 1),]
# # For now, assume they are captures: replace with 1
# CCBData_orig$Capture[which(CCBData_orig$Capture > 1)] <- 1
# Also many NA's
CCBData_orig[which(!CCBData_orig$Capture%in%(0:1)),]
# For now set them all to NA
CCBData_orig$Capture[which(!CCBData_orig$Capture%in%(0:1))] <- NA

# Label all non-localized calls as NA
CCBData_orig$Lat[which(CCBData_orig$Lat==0)]=NA
# They are already NA's
# all columns of these observations are NA's


# Add an estimate for when the call was produced
CCBData_orig$CallProudcedTime = CCBData_orig$date - CCBData_orig$DistanceKm*1000/1497






#________________________________
# Subset data

# Subset the data to only localised calls, without missing noiselevel or capture value
NLsub <- subset(CCBData_orig, !is.na(Lat) & !is.na(NoiseLevel) & !is.na(Capture))

# Trim the data (remove 5% most extreme noise levels)
quantiles <- quantile(NLsub$NoiseLevel, c(0.025, .975))
NLsub <- subset(NLsub, NoiseLevel>quantiles[1] & NoiseLevel<quantiles[2])

NLsub$CallProudcedTime = as.POSIXct(NLsub$CallProudcedTime, tz = 'UTM')
NLsub$NumCallTime = unclass(NLsub$CallProudcedTime)

#________________________________
# Create the aggrigate table plotting location vs time of day
CallTimeDf = group_by(NLsub, ID)
CallTimeDf<- summarize(CallTimeDf, mean(CallProudcedTime))
colnames(CallTimeDf)[2] ='CallProudcedTime'

CallTimeDf$UselessDate = aggregate(data=NLsub, UselessDate~ID, FUN = mean)[,2]
CallTimeDf$SecondsSinceMidnight = aggregate(data=NLsub, StartTime~ID, FUN = mean)[,2]

# Get average location
CallTimeDf$Lat = aggregate(data=NLsub, Lat~ID, FUN = median)[,2]
CallTimeDf$Lon = aggregate(data=NLsub, Lon~ID, FUN = median)[,2]

# Reset the ID's
CallTimeDf=CallTimeDf[order(CallTimeDf$CallProudcedTime),]
CallTimeDf$ID = seq(1, nrow(CallTimeDf))

# Add column for hour of the day
CallTimeDf$hr = as.numeric(strftime(CallTimeDf$CallProudcedTime, "%H"))

#==============================================
#  Do some exploratory plotting
#==============================================

library(ggplot2)

ggplot(CallTimeDf, aes(Lon, Lat, col=SecondsSinceMidnight)) + 
  facet_wrap(~UselessDate)  +
  geom_point() 

ggplot(CallTimeDf[0:20,], aes(Lon,Lat, col=SecondsSinceMidnight)) + 
  geom_point() 


#==============================================
#  Figure out which calls could be connected
#==============================================
CallTimeDf$TrackNumber =NA
CallTimeDf$TrackNumber[1] = 1
TrackNumber = 1
CallTimeDf$SwimSpeedMS = NA


# Assume max swimming speed of 10mph http://www.listenforwhales.org/page.aspx?pid=451
MaxSwimSpeed = 4.4 #m/s

CurrID=1
MinSpeed =0 

while(sum(is.na(CallTimeDf$TrackNumber))>0){
  

  # Check if swim speed is less than max if so add that call to the track
 
   # Get current time and location
   Time = CallTimeDf$CallProudcedTime[CurrID]
   Loc = c(CallTimeDf$Lat[CurrID], CallTimeDf$Lon[CurrID])
   
   
   # Make a subset with only times after the initial calls
   Time_sub = subset(CallTimeDf, is.na(CallTimeDf$TrackNumber) &
                       CallProudcedTime>Time & 
                       CallProudcedTime<(Time+(60*150)))
   
   Time_sub$DeltaSec = as.numeric(Time_sub$CallProudcedTime - Time)
   
   
   # Assuming animals are unlikely to produce calls within 1 second of another call, trim
   # to remove concurrent calls
   Time_sub = subset(Time_sub, DeltaSec>2.2)
   Time_sub$DeltaXY = sqrt((Time_sub$Lat-Loc[1])^2 +(Time_sub$Lon-Loc[2])^2)
   
   
   Time_sub$SwimSpeed = Time_sub$DeltaXY/Time_sub$DeltaSec
   Time_sub = Time_sub[order(Time_sub$SwimSpeed),]
   
    
   # If there is a call within 20 minutes and the animal could have feasibly gotten there
    if (nrow(Time_sub)>1 & Time_sub$SwimSpeed[1]<MaxSwimSpeed){
        ID = Time_sub$ID[1]
        CallTimeDf$TrackNumber[ID] = TrackNumber
        CallTimeDf$SwimSpeedMS[CurrID] = Time_sub$SwimSpeed[1]
        
        # Move to the next call in the sequence
        CurrID = ID
        
    }else{

      # Run out of max times, create a new track
      TrackNumber = TrackNumber + 1
    
      # Advance the current ID
      CurrID = min(which(is.na(CallTimeDf$TrackNumber)))  
      CallTimeDf$TrackNumber[CurrID] = TrackNumber
      CallTimeDf$SwimSpeedMS[CurrID] =  0
  }

  print(as.character(TrackNumber))
}




#==============================================
#  Do some exploratory plotting of the tracks
#==============================================

# Create a dummy column to calculate tracks per day
CallTimeDf$TracksPerDay = NA

for (ii in 1:6){
  
  ids = which(CallTimeDf$UselessDate==unique(CallTimeDf$UselessDate)[ii])
  CallTimeDf$TracksPerDay[ids] = as.numeric(as.factor(CallTimeDf$TrackNumber[ids]))
  
}


# Report the track start time
CallTimeDf = group_by(CallTimeDf, TrackNumber)
CallTimeDf_new<- summarize(CallTimeDf, min(CallProudcedTime))
colnames(CallTimeDf_new)[2] ='TrackStartTime'




CallTimeDf = merge(CallTimeDf, CallTimeDf_new, by='TrackNumber', all.x=T)



for (ii in 1:6){
  data_sub = subset(CallTimeDf, UselessDate==unique(CallTimeDf$UselessDate)[ii])
  p<-ggplot(data_sub, aes(Lon,Lat, col=TracksPerDay, group=TrackStartTime)) + 
    geom_point(size=.5, color='black') +
    #geom_line(size=.5, arrow = arrow(length=unit(.2, 'cm'), angle = 15, type = "closed"))+
    scale_color_viridis()+
    ggtitle(as.Date(data_sub$CallProudcedTime[1])) 
  print(p)  
  
  
  
}
  
# Pull out a subset of the track for testing


track_sub = subset(CallTimeDf, UselessDate=='218' & 
                     Lat>4640000 & Lat<4650000 & 
                     Lon>380000 & Lon<390000)

ggplot(track_sub, aes(Lon,Lat, col=SecondsSinceMidnight)) + 
  geom_point(size=.5) +
  #geom_line(size=.5, arrow = arrow(length=unit(.2, 'cm'), angle = 15, type = "closed"))+
  scale_color_viridis()
print(p)  















# Get location change between calls
CallTimeDf$deltaLoc_m = c(0, sqrt(diff(CallTimeDf$Lat)^2 +diff(CallTimeDf$Lon)^2))



# Required maximum swimm speed
CallTimeDf$SwimSpeedKmHr = CallTimeDf$deltaLoc_m/ CallTimeDf$TimeDiff *(60*60)/1000


# Restrict data points to only those in the within the bounds of the array
NLsub$Included <- point.in.polygon(point.x=NLsub$Lon, point.y=NLsub$Lat,
                                   MARUloc$x[c(1,2,3,4,5,9,10,8,1)], MARUloc$y[c(1,2,3,4,5,9,10,8,1)])==1

dat <- subset(NLsub, Included==TRUE)
dat <- droplevels(dat)
str(dat)



