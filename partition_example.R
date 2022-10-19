
# required packages 
require(tidyverse)
require(lubridate)
require(mgcv)
require(ANOM)
require(nlme)
require(multcomp)
require(tidyr)
require(minpack.lm)

# load in data 
setwd("/Volumes/MMills1/Manuscript")
flux<- read.csv("SAFE_FluxTower_GapFilled_12-18_withCIcorrection.csv", header = TRUE)
head(flux)

flux<- flux %>%
  select(Year, 
         Hour, 
         DoY, 
         Rg, # global radiation in Wm-2 
         VPD, # Vapour pressure defict in hPa
         NEE_orig) # observed NEE values, not gapfilled

# global radiation to be converted to PAR
flux$PAR<- flux$Rg/0.5143

# filter data to 08:00 - 19:00 
data.day<- subset(flux, flux$Hour<19.5)
data.day<- subset(data.day, data.day$Hour>7.5)
summary(data.day$Hour)

# remove any NA values in observed NEE
data.day <- data.day %>%
  drop_na(NEE_orig)

#######################  light saturdated VPD 

# VPD0 threshold of 10hPa (Körner, 1995)
# Ch. Körner, “Leaf diffusive conductances in the major vegetation types of the globe”
# in Ecophysiology of Photosynthesis, 100th Ed., E. D. Schulze, M. M. M.Caldwell, Eds. 
# (Springer Study Edition, 1995), pp. 463–490.
VPD.0<- 10

# create a subset of light saturated conditions where PAR >1200 
data.all.lightsaturated <- subset(data.day, PAR>1200) 
data.all.lightsaturated.VPDlimited <- subset(data.all.lightsaturated, VPD>VPD.0) 

m.VPD.all.lightsaturated <- nlsLM(NEE_orig ~ -(Bm*exp(-k*(VPD-VPD.0))),
                                  data = data.all.lightsaturated.VPDlimited, start = list(k=0.4, Bm = 20),
                                  na.action=na.exclude, control=nls.lm.control(maxiter=1024,maxfev=1024))

summary(m.VPD.all.lightsaturated)
k <- as.numeric(coef(m.VPD.all.lightsaturated)[1])

####################### 
# now establish constant value for \alpha\, which is the slope of the linear 
# regression between NEE and PAR under low light conditions (<200 μmol) (Xu et al., 2019).
# J. Xu, et al., A general non-rectangular hyperbola equation for photosynthetic light response 
# curve of rice at various leaf ages. Sci Rep 9, 1–8 (2019)

# I am running this as a subset using only data from 2018 
# this is optional 
f18.day<- subset(data.day, data.day$Year==2018)

# be aware the code does not run if there are large gaps 
# so i am subsetting it to a shorter period 
tapply(f18.day$NEE_orig, f18.day$DoY, length)
f18.day.A<- subset(f18.day, f18.day$DoY>0 & f18.day$DoY<61)


#check the slope 
f18.day.lowPAR <- subset(f18.day, PAR<=200)
plot(NEE_orig ~ PAR, data=f18.day.lowPAR)
lm.lowPAR.f18 <- lm(NEE_orig ~ PAR, data=f18.day.lowPAR)
summary(lm.lowPAR.f18)
as.numeric(coef(lm.lowPAR.f18)[2])

# Am represents \alpha\ which is canopy light utilization
Am.18 <- as.numeric(coef(lm.lowPAR.f18)[2]) * -1

###########################################################
estimate.Bm <- vector("numeric") # the maximum CO2 uptake rate of the canopy at light saturation
estimate.Ym <- vector("numeric") # ecosystem respiration 
estimate.Am <- vector("numeric") # anopy light utilization which represents the initial slope of the light–response curve
estimate.Bm.SE <- vector("numeric")
estimate.Ym.SE <- vector("numeric")
estimate.Am.SE <- vector("numeric")
Data.7day.length <- vector("numeric")
StartDay <- vector("numeric")
EndDay <- vector("numeric")

DoY.seq <- seq(from=min(f18.day.A$DoY), to=max(f18.day.A$DoY), by=1)


# set the fixed values:
# Fixed values
VPD.0 <- 10
k <- as.numeric(coef(m.VPD.all.lightsaturated)[1])
Am <- as.numeric(coef(lm.lowPAR.f18)[2]) * -1


for(i in seq_along(DoY.seq)){
  
  StartDay[i] <- DoY.seq[i]-3  
  EndDay[i] <- StartDay[i] + 7
  
  Data.7day <- f18.day.A[f18.day.A$DoY>=StartDay[i] & f18.day.A$DoY<=EndDay[i] , ]
  
  lrc.f18.day.A <- nlsLM(NEE_orig ~ ifelse(VPD > VPD.0, (-(Bm*exp(-k*(VPD-VPD.0)) + Ym) * (1-exp((-Am*PAR)/(Bm *exp(-k*(VPD-VPD.0)) + Ym))) + Ym), (-(Bm + Ym) * (1-exp((-Am*PAR)/(Bm + Ym))) + Ym)),        
                         data = Data.7day, start = list(Bm = 5, Ym = 0),
                         na.action=na.exclude, control=nls.lm.control(maxiter=1024,maxfev=1024))
  
  estimate.Bm[i] <- summary(lrc.f18.day.A)$coefficients[[1]]
  estimate.Ym[i] <- summary(lrc.f18.day.A)$coefficients[[2]]
  estimate.Bm.SE[i] <- summary(lrc.f18.day.A)$coefficients[[3]]
  estimate.Ym.SE[i] <- summary(lrc.f18.day.A)$coefficients[[4]]
  
}

# make a table and save the results:
f18.day.A.part<- cbind(DoY.seq, estimate.Ym, estimate.Ym.SE, estimate.Bm, estimate.Bm.SE)
f18.day.A.part<- as.data.frame(f18.day.A.part)


