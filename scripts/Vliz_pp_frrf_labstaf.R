title: "Estimating integrated production by merging FRRF(Vliz edition) or Labstaf data with CTD data"
author: "Jens H. Dujardin"
date: "18-July-2024"

### Libraries ###

library(FME)
library(deSolve)
library(dtWad)
library(dtBioG)
library(tidyverse)
library(ggsci)

### Loading Labstaf data###

DIR <- "../raw_data"
ff <- list.files(DIR, pattern=".csv")
LABSlist <- NULL
1
for (f in ff)
  LABSlist <- rbind(LABSlist,
                    data.frame(file=f,
                               read.csv2(paste(DIR,f, sep="/"))[,1:22]))
cf <- colnames(LABSlist)
cf[cf %in% c("F.","Fm.")] <-
  c("Fo","Fm")
colnames(LABSlist) <- cf
LABSlist[,-(1:2)] <- lapply(LABSlist[,-(1:2)], FUN=as.numeric)
DIRS <- c("../raw_data/45O", "../raw_data/48F", "../raw_data/49O",
          "../raw_data/50B", "../raw_data/51O", "../raw_data/53M",
          "../raw_data/54O", "../raw_data/56F", "../raw_data/58B",
          "../raw_data/59B", "../raw_data/60B")
rLABs <- function(fn, dir){
  specs <- strsplit(fn, "m.txt")[[1]]
  specs <- strsplit(specs, "_")[[1]]
  specs <- specs[c(1, length(specs))]
  names(specs) <- c("station", "depth")
  fr <- readFRRF(file=fn, dir=dir,
                 specs=specs ,
                 txt = "delim")
  fr
}
Flist <- NULL
for (d in DIRS){
  ff <- list.files(d, pattern=".txt")
  for (f in ff)
    Flist <- rbind(Flist,rLABs(f, d))
}
FF <- Flist[,colnames(LABSlist)]
LabSTAFlist <- rbind(LABSlist, FF)
LabSTAFlist$statNr <- substr(LabSTAFlist$station, 1, nchar(LabSTAFlist$station)-1)
LabSTAFlist$longitude <- LabSTAFlist$latitude <- NA
LabSTAFlist[,3:24] <- lapply(LabSTAFlist[,3:24], FUN=as.numeric)

### Add Fblanc ###

FRRFsettings <- read.csv("../raw_data/settings/FRRFsettings.csv")
Fblanc <- aggregate(FRRFsettings$Fblanc, by=list(FRRFsettings$nr), FUN=mean)
Fblanc <- rbind(Fblanc, c(3, 0.1)) # Fblanc is unknown for station 3
colnames(Fblanc) <- c("statNr", "Fblanc")
LabSTAFlist <- merge(LabSTAFlist, Fblanc, all.x=TRUE)

### add latitude and longitude ###

for (i in 1:nrow(stations)){
  ii <- which (LabSTAFlist$station == stations$station[i])
  if (length(ii)){
    LabSTAFlist$latitude[ii] <- stations$latitude[i]
    LabSTAFlist$longitude[ii] <- stations$longitude[i]
  }
}

### Processing the LabSTAF data ###

LabSTAFlist$JVPII <- LabSTAFlist$JVPII*3600/1000

profile <- unique(LabSTAFlist[,c("station", "depth")])
LSlist <- NULL
for (i in 1: nrow(profile)){
  sub <- subset(LabSTAFlist,
                subset=station == profile[i,1] &
                  depth == profile[i,2] )
  std <- standardizeFRRF(sub,
                         Fblanc = sub$Fblanc,
                         convJVPII = 3.6)
  LSlist <- rbind(LSlist, std)
}
# negative JVPII values are set to NA
LSlist$JVPII[LSlist$JVPII <= 0 ] <- NA

save(file = "../processed_data/LSlist.rda", LSlist)
save(file = "../processed_data/LabSTAFlist.rda", LabSTAFlist)
write.csv(file = "../processed_data/LSlist.csv", LSlist)

### Fitting the PI curves ###

par(mfrow=c(6,5), mar=c(4,4,4,0), las=1)
# unique combinations of station and measurement depth
Uin <- unique(data.frame(station=LSlist$station,
                         depth =LSlist$depth))
# Fit all PI curves and put in one list
FIT <- NULL
Pmax <-1000 # ignore measurements above 1000 uEinst/m2/s (i.e. no)
minRows <- 4
Uin$used <- TRUE
for (i in 1:nrow(Uin)) {
  frrf <- subset(LSlist,
                 subset = station==Uin$station[i] &
                   depth ==Uin$depth[i])
  if(Uin$station[i] == "53M" & Uin$depth[i]==5){
    nf <- nrow(frrf)
    NAremove <- frrf[(nf-1):nf,]
    frrf <- frrf[1:(nf-2),]
  } else NAremove <- list(E=NA, JVPII=NA)
  frrf <- frrf[! is.na(frrf$JVPII), ]
  if (nrow(frrf) > minRows){
    STAT <- paste(Uin$station[i], Uin$depth[i]) # to label the plot
    fitone <- fitPI(I = frrf$E, response = frrf$JVPII, model="EP")
    xlim <- range(c(frrf$E, NAremove$E), na.rm=TRUE)
    ylim <- range(c(frrf$JVPII, frrf$JVPII_uc, NAremove$JVPII), na.rm=TRUE)
    plot(fitone, main=STAT, ylim=ylim, pch=18, xlim=xlim)
    points(frrf$E, frrf$JVPII_uc, pch=18, col="grey")
    if (!is.na(NAremove$E[1])) points(NAremove$E, NAremove$JVPII, col="red")
    FIT <- rbind(FIT,
                 data.frame(
                   station = Uin$station[i],
                   depth = Uin$depth[i],
                   alpha = fitone$par["alpha"],
                   eopt = fitone$par["eopt"],
                   ps = fitone$par["ps"],
                   rsq = r.squared(fitone))
    )
  } else Uin$used[i] <- FALSE
}

FIT <- as.data.frame(FIT)
FIT$statnr <- as.integer(substr(FIT$station, 1, nchar(FIT$station)-1))
fs <- unique(FIT$station)
FIT$chl_ctd <- NA
chl_c <- NULL
for (s in fs){
  ii <- which(FIT$station==s)
  ctd <- subset(CTDall, subset=station==s)
  if (nrow(ctd) > 0){
    for (i in ii){
      D <- FIT$depth[i]
      Dmin <- max(0, D-1)
      ctd_s <- subset(ctd, subset=
                        Depth >= Dmin & Depth <= D+1 &
                        ! is.na(Chlorophyll) & Chlorophyll > 0)
      chl_c <- mean(ctd_s$Chlorophyll, na.rm=TRUE)
      FIT$chl_ctd[i] <- chl_c
    }
  }
}

### Here is a list of the depth strata per station whose data were either successfully fitted or ignored ###
USE <- cbind(aggregate(Uin$used, by=list(Uin$station), FUN=sum),
             aggregate(1-Uin$used, by=list(Uin$station), FUN=sum)[,2])
names(USE) <- c("station", "used", "ignored")
USE <- merge(stations[, 1:2], USE)
USE <- USE[order(USE$statNr),]
knitr::kable(USE, row.names = FALSE,
             caption="Depth levels per station that were fitted (used) or had no good JVPII values (ignored)")

### Resulting parameters ###
fit_Labstaf <- FIT
fit_Labstaf <- merge(fit_Labstaf,
                     stations[,c("station", "statNr", "latitude", "longitude", "Fjord")],
                     by="station")
knitr::kable(fit_Labstaf[,c("station", "depth", "alpha", "eopt", "ps", "chl_ctd", "rsq")],
             digits=c(0,0,4,0,2, 2, 2),
             caption="Photosynthesis parameters, derived from the LabSTAF measurements")


fit_Labstaf$station <- factor(fit_Labstaf$station,
                              levels = c("2O", "3M" ,
                                         "4O", "5O", "6B", "7F", "8O", "9B",
                                         "10O", "11M", "12O", "13B", "14O", "15F",
                                         "15O", "16O", "17O", "19O", "45O", "48F",
                                         "49O", "50B", "51O", "52B", "53B", "54O",
                                         "55O", "56F", "57O", "58B", "59B", "60B" ))
fit_Labstaf <- fit_Labstaf[order(as.integer(fit_Labstaf$statNr)),]
save(file="../processed_data/fit_Labstaf.rda", fit_Labstaf)
write.csv(file="../processed_data/fit_Labstaf.csv", fit_Labstaf, row.names = FALSE)

### Plotting the measurement positions ###
J0 <- 19553 # Julian day for the first day of the cruise
par(mfrow=c(1,2), las=1)
asp <- 1/cos((mean(LSlist$latitude) * pi)/180)
with(LSlist,
     plot(longitude, latitude, asp=asp, pch=16,
          main="LABSTAF station data", xlim=c(-46.7, -45.3),
          xlab="Longitude",
          ylab="Latitude"))
lines2D(LATLONcoast$longitude, LATLONcoast$latitude,
        type="l", col="grey", add=TRUE)

LL <- aggregate(LSlist$depth, by=list(LSlist$station),
                FUN=function(x) length(unique(x)))
names(LL) <- c("station", "numDepths")
LL <- merge(LL, stations[, c(1:5, ncol(stations))], by="station")
LL <- LL[order(as.numeric(LL$statNr)),]
knitr::kable(LL, digits=c(0,0,0,3,3,0), row.names=FALSE,
             caption="position of stations and number of Labstaf measurment depths")

### About the standardization ###
par(mfrow=c(2,2))
f_ori <- subset(LabSTAFlist, subset = station== "54O" & depth == 40)
frrf <- standardizeFRRF(f_ori, Fblanc=f_ori$Fblanc)
f_ori2 <- subset(LabSTAFlist, subset = station== "54O" & depth == 1)
frrf2 <- standardizeFRRF(f_ori2, Fblanc=f_ori2$Fblanc)
with(frrf, matplot(x=E, y=cbind(JVPII, JVPII_uc), pch=16:17,
                   ylab="JVPII", main="station 54O, 40m"))
legend("right", legend=c("Fblanc=0", "Fblanc=0.132"), col=1:2, pch=16:17)
with(frrf2, matplot(x=E, y=cbind(JVPII, JVPII_uc), pch=16:17,
                    ylab="JVPII", main="station 54O, 1m"))
legend("right", legend=c("Fblanc=0", "Fblanc=0.132"), col=1:2, pch=16:17)
ZZ <- standardizeFRRF(f_ori,
                      Fblanc = 0,
                      aLHII_0 = 0.007194)
ZZb <- standardizeFRRF(f_ori,
                       Fblanc = 0)
Ztrue <- standardizeFRRF(f_ori,
                         Fblanc = f_ori$Fblanc)
data.frame(corrected=Ztrue$JVPII, corrFb0=ZZ$JVPII, uncorrected = ZZ$JVPII_uc) # difference

(aLHII_0_true <- attributes(Ztrue)$aLHII_0) # compare with 0.007194 !!!

(aLHII_0_LS <- attributes(ZZ)$aLHII_0) # compare with 0.007194 !!!

(aLHII_0_b0 <- attributes(ZZb)$aLHII_0) # compare with 0.007194 !!!

plot(ZZ$E, ZZ$a_LHII, ylim=c(0.0001, 0.03), log="y", pch=17, col=1)
points(Ztrue$E, Ztrue$a_LHII, pch=16, col=2)
abline(h=c(aLHII_0_true,aLHII_0_LS,aLHII_0_b0), col=1:3, lty=2)
legend("right", legend=c("Fblanc=0", "Fblanc=0.132"), col=1:2, pch=16:17)
ZZ2 <- standardizeFRRF(f_ori2,
                       Fblanc = 0,
                       aLHII_0 = 0.0132)
ZZb2 <- standardizeFRRF(f_ori2,
                        Fblanc = 0)
Ztrue2 <- standardizeFRRF(f_ori2,
                          Fblanc = f_ori2$Fblanc)
aLHII_0_true2 <- attributes(Ztrue2)$aLHII_0 # compare with 0.0132 !!!
aLHII_0_LS2 <- attributes(ZZ2)$aLHII_0 # compare with 0.0132 !!!
aLHII_0_b02 <- attributes(ZZb2)$aLHII_0 # compare with 0.0132 !!!

plot(ZZ2$E, ZZ2$a_LHII, ylim=c(0.0001, 0.03), log="y", pch=17, col=1)
points(Ztrue2$E, Ztrue2$a_LHII, pch=16, col=2)
abline(h=c(aLHII_0_true2,aLHII_0_LS2,aLHII_0_b02), col=1:3, lty=2)
legend("right", legend=c("Fblanc=0", "Fblanc=0.132"), col=1:2, pch=16:17)




