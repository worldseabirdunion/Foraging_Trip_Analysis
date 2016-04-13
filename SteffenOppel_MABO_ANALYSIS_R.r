############################################################################################################
#######  COMPARISON OF MABO FORAGING TRIPS FROM ASCENSION AND ST HELENA  ###################################
############################################################################################################
## this analysis is based on "A:\RSPB\UKOT\MABO_COMPARISON_ASI_StHx64.r"
## reduced on 26 June 2014, see above script for basic analyses of trips etc.
## marine environmental data and maps were outsourced to "A:\RSPB\Marine\MARINE_CHARTS.r"

## PLEASE CITE:
##Oppel, S., Beard, A., Fox, D., Mackley, E., Leat, E., Henry, L., Clingham, E., Fowler, N., Sim, J., 
##Sommerfeld, J., Weber, N., Weber, S., Bolton, M., 2015. Foraging distribution of a tropical seabird 
##supports Ashmole's hypothesis of population regulation. Behavioral Ecology and Sociobiology 69: 915-926.



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(CircStats)
library(fields)
library(circular)
library(marmap)
library(boot)
require(maps)
require(mapdata)
require(adehabitat)
#require(adehabitatHR)
#require(adehabitatLT)
require(gpclib)
require(foreign)
require(maptools)
require(geosphere)
require(sp)
require(rgdal)
require(rgeos)
library(raster)
library(RHmm)
library(BBMM)
library(move)
library(rptR) 
library(plotrix)
library(lme4)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM R WORKSPACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### this work space includes all summary tables and raw data, last saved on 26 June 2014
### reduced from previous very long version with command below

# load("C:\\STEFFEN\\RSPB\\UKOT\\StHelena\\Science\\Birds\\seabirds\\MABO_analysis3.RData")
# save.image("C:\\STEFFEN\\RSPB\\UKOT\\StHelena/Science/Birds/seabirds\\MABO_analysis3.RData")
# rm(list=ls()[c(1,2,6,9,13,14,16,17,18,20,21,22,23,24,25,26,27,28,30,31:67,69:70,72:77,79:86,95:108)])
# save.image("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\Ascension\\MABO_track_summaries.RData")


#load("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_press\\Ascension\\MABO_track_summaries.RData")
load("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\Ascension\\MABO_track_summaries.RData")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT THE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#pdf("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\Ascension\\Fig1.pdf", width=8, height=7)
#jpeg("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\Ascension\\Fig1.jpg", width=8, height=7, units="in", quality = 75, res=200)

#pdf("A:\\MANUSCRIPTS\\submitted\\Ascension\\Fig1.pdf", width=8, height=7)
#jpeg("A:\\MANUSCRIPTS\\submitted\\Ascension\\Fig1.jpg", width=8, height=7, units="in", quality = 75, res=200)

### remove the outlier track from St Helena:
outlier<-tracks[tracks$ID==204,]
tracks<-tracks[!tracks$ID==204,]


par(mar=c(5,6,0,0), oma=c(0,0,0,0))
plot(Latitude~Longitude, data=tracks[tracks$island=="Ascension",], type='l',col='gray10', xlim=c(-18,-2), ylim=c(-18,-4), main="", frame=F, axes=F, cex.lab=1.8, cex.axis=1.8, mgp=c(3.4,0.7,0), asp=1)
par(new=T)
plot(Latitude~Longitude, data=tracks[tracks$island=="StHelena",], type='l',col='gray10', xlim=c(-18,-2), ylim=c(-18,-4), frame=F, axes=F, cex.lab=1.8, cex.axis=1.8, xlab="", ylab="", asp=1)
par(new=T)
plot(Latitude~Longitude, data=outlier, type='l',col='gray50', xlim=c(-18,-2), ylim=c(-18,-4), frame=F, axes=F, cex.lab=1.8, cex.axis=1.8, xlab="", ylab="", asp=1)

map("worldHires", add=T, fill=T, col="grey91", asp=1)

draw.circle(-5.734380,-16.003095,radius=1,border='gray27',lty=2,lwd=1)
draw.circle(-14.301331,-7.946484,radius=1,border='gray27',lty=2,lwd=1)
draw.circle(-5.734380,-16.003095,radius=3,border='gray27',lty=1,lwd=1)
draw.circle(-14.301331,-7.946484,radius=3,border='gray27',lty=1,lwd=1)
axis(1,at=seq(-18,-4,2), labels=T, cex=2, cex.lab=1.8, cex.axis=1.8)
axis(2,at=seq(-18,-4,2), labels=T, cex=2, cex.lab=1.8, cex.axis=1.8, las=1, mgp=c(3.4,0.7,0))
text(-8,-8,"Ascension Island", cex=1.8)
text(-11,-16,"St Helena", cex=1.8)
#dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARY TABLE AND NUMBERS FOR PAPER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aggregate(n_trips~island+season, data=deployments, FUN=sum)

summary<-aggregate(max_dist~Stage+island, data=trip_distances, FUN=mean)
summary$dist_sd<-aggregate(max_dist~Stage+island, data=trip_distances, FUN=sd)[,3]
summary$time<-aggregate(time~Stage+island, data=trip_distances, FUN=mean)[,3]
summary$time_sd<-aggregate(time~Stage+island, data=trip_distances, FUN=sd)[,3]
summary$length<-aggregate(total_dist~Stage+island, data=trip_distances, FUN=mean)[,3]
summary$length_sd<-aggregate(total_dist~Stage+island, data=trip_distances, FUN=sd)[,3]
summary$n_trips<-aggregate(sum~Stage+island, data=trip_distances, FUN=sum)[,3]
summary$trips_per_day<-aggregate(trips_per_day~breeding_status+island, data=deployments, FUN=mean)[,3]
summary$n_ind<-aggregate(count~breeding_status+island, data=deployments[deployments$n_trips>0,], FUN=sum)[,3]		## manually extracted from aggregate(ID~Stage+island, data=trip_distances, FUN=unique)[,3]
summary

papertable<-summary[,c(2,1,11,9,10)]
papertable$max_dist<-paste(round(summary$max_dist,0)," + ", round(summary$dist_sd,0))
papertable$total_dist<-paste(round(summary$length,0)," + ", round(summary$length_sd,0))
papertable$time<-paste(round(summary$time,1)," + ", round(summary$time_sd,1))
#write.table(papertable, "C:\\STEFFEN\\MANUSCRIPTS\\submitted\\Ascension\\Table1.csv", sep=",", row.names=F)


hist(trip_distances$max_dist[trip_distances$island=="StHelena"], breaks=100)
hist(trip_distances$max_dist[trip_distances$island=="StHelena"], breaks=100)

hist(trip_distances$time[trip_distances$island=="StHelena"], breaks=10, xlim=c(0,100), ylim=c(0,80))
par(new=T)
hist(trip_distances$time[trip_distances$island=="Ascension"], breaks=10, xlim=c(0,100), ylim=c(0,80), col='red')


mean(deployments$n_trips)
aggregate(n_trips~island, data=deployments, FUN=mean)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT TRIP TIME AGAINST DISTANCE (FIG 2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#pdf("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\Ascension\\Fig2.pdf", width=10, height=7)
#jpeg("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\Ascension\\Fig2.jpg", width=10, height=7, units="in", quality = 75, res=200)

#pdf("A:\\MANUSCRIPTS\\submitted\\Ascension\\Fig2.pdf", width=10, height=7)
#jpeg("A:\\MANUSCRIPTS\\submitted\\Ascension\\Fig2.jpg", width=10, height=7, units="in", quality = 75, res=200)

trip_distances<-trip_distances[!trip_distances$ID==204,]

par(mar=c(5,5,0,0), oma=c(0,0,0,0))
plot(max_dist~time, data=trip_distances[trip_distances$island=="StHelena",], pch=21,bg="white",xlim=c(0,90), ylim=c(0,350), las=1, cex=1.7, cex.lab=1.5, cex.axis=1.5, axes=F, xlab="", ylab="")
par(new=T)
plot(max_dist~time, data=trip_distances[trip_distances$island=="Ascension",], pch=16, xlim=c(0,90), ylim=c(0,350), las=1, cex=1.7, cex.lab=1.5, cex.axis=1.5, axes=F, xlab="Trip duration (hrs)", ylab="Maximum distance from colony (km)")
axis(1,at=seq(0,90,10), labels=T,cex.axis=1.3)
axis(2,at=seq(0,350,50), cex.axis=1.3, las=1, mgp=c(4,0.5,0), tck = -0.009)
legend(0,350,pch=c(16,1), legend=c("Ascension", "St Helena"), cex=1.5, bty = "n")
abline(v=20, lty=2, col='gray50')

#dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DO INDIVIDUALS ON ST HELENA HAVE MORE TRIPS PER DAY?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(glm(trips_per_day~breeding_status+island, data=deployments))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DO TRIP CHARACTERISTICS DIFFER BETWEEN ISLANDS (after accounting for stage and sex)?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### MAXIMUM DISTANCE FROM COLONY ###

max.dist.m0<-lmer(max_dist~Stage+(1|season/ID), data=trip_distances)
max.dist.m1<-lmer(max_dist~Stage+island+(1|season/ID), data=trip_distances)
#dredge(max.dist.m1)
summary(max.dist.m1)
#max.dist.m1@fixef[3]/max(summary[,3])
max.dist.m1@beta[4]/max(summary[,3])

#Likelihood ratio test to check whether m1 is better than m0#
anova(max.dist.m0,max.dist.m1, test = "Chisq")




### TOTAL TRIP DISTANCE ###

tot.dist.m0<-lmer(total_dist~Stage+(1|season/ID), data=trip_distances)
tot.dist.m1<-lmer(total_dist~Stage+island+(1|season/ID), data=trip_distances)
dredge(tot.dist.m1)
summary(tot.dist.m1)
#tot.dist.m1@fixef[3]/max(summary[,7])
tot.dist.m1@beta[4]/max(summary[,7])

#Likelihood ratio test to check whether m1 is better than m0#
anova(tot.dist.m0,tot.dist.m1, test = "Chisq")



### TOTAL TRIP DURATION ###

time.m0<-lmer(time~Stage+(1|season/ID), data=trip_distances)
time.m1<-lmer(time~Stage+island+(1|season/ID), data=trip_distances)
#dredge(time.m1)
summary(time.m1)
#time.m1@fixef[3]/max(summary[,5])			### for old version of lme4
time.m1@beta[4]/max(summary[,5])

#Likelihood ratio test to check whether m1 is better than m0#
anova(time.m0,time.m1, test = "Chisq")



### TRIP DIRECTIONS #################################

aov.circular(circular(trip_distances[,10],type = "directions", units = "degrees"),trip_distances$island)
mean(circular(trip_distances[trip_distances$island=="Ascension",10],type = "directions", units = "degrees", modulo = "2pi"))
mean(circular(trip_distances[trip_distances$island=="StHelena",10],type = "directions", units = "degrees", modulo = "2pi"))

median(circular(trip_distances[trip_distances$island=="Ascension",10],type = "directions", units = "degrees", modulo = "2pi"))
median(circular(trip_distances[trip_distances$island=="StHelena",10],type = "directions", units = "degrees", modulo = "2pi"))

plot(circular(trip_distances[trip_distances$island=="Ascension",10],type = "directions", units = "degrees", modulo = "2pi"),pch=16, cex=1.2,units = "degrees", template =  "geographics", zero = 0, rotation = "clock", shrink=1.5, tol=0, tcl.text=-0.4, stack=T)
par(new=T)
plot(circular(trip_distances[trip_distances$island=="StHelena",10],type = "directions", units = "degrees", modulo = "2pi"),pch=1, cex=1.2,units = "degrees", template =  "geographics", zero = 0, rotation = "clock", shrink=1.5, tol=0, tcl.text=-0.4, stack=T)
legend(-1.5,1.5,pch=c(16,1), legend=c("ASI", "STH"), cex=1.7, bty = "n")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ESTIMATE REPEATABILITIES IN DISTANCE, TIME, AND DIRECTION FOR INDIVIDUALS AND SPECIES
# PREDICTION: St Helena could be more consistent than Ascension
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### sample size of birds with >1 trip
aggregate(count~island, data=deployments[deployments$n_trips>1,], FUN=sum)


repeatabilities<-data.frame(island=rep(c("Ascension","StHelena"), each=4), trip_char=rep(c("time","max_dist","tot_dist","direction"),2), r=0, lcl=0,ucl=0,p=0)


r1<-rpt(trip_distances$time[trip_distances$island=="StHelena"], trip_distances$ID[trip_distances$island=="StHelena"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="time",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_distances$time[trip_distances$island=="Ascension"], trip_distances$ID[trip_distances$island=="Ascension"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="time",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_distances$max_dist[trip_distances$island=="StHelena"], trip_distances$ID[trip_distances$island=="StHelena"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="max_dist",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_distances$max_dist[trip_distances$island=="Ascension"], trip_distances$ID[trip_distances$island=="Ascension"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="max_dist",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_distances$total_dist[trip_distances$island=="StHelena"], trip_distances$ID[trip_distances$island=="StHelena"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="tot_dist",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_distances$total_dist[trip_distances$island=="Ascension"], trip_distances$ID[trip_distances$island=="Ascension"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="tot_dist",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])



## MANUAL CALCULATION OF REPEATABILITY FOLLOWING LESSELS AND BOAG 1987

#### UPDATED REVIEW 7 March 2015 to add CI and p-value ####
#### SE estimation from Sam Patrick

n_trips<-aggregate(count~ID, data=trip_distances[trip_distances$island=="Ascension",], FUN=sum)
tab<-table(n_trips$count)
freq<-as.numeric(dimnames(tab)[[1]])
sum_ni<-sum(tab*freq)
sum_ni2<-sum(tab*freq*freq)

n0<-(1/(length(unique(trip_distances$ID[trip_distances$island=="Ascension"]))-1))*(sum(n_trips$count)-(sum(n_trips$count^2)/sum(n_trips$count)))			## mean sample size
test<-aov.circular(circular(trip_distances[trip_distances$island=="Ascension",10],type = "directions", units = "degrees", modulo = "2pi"),trip_distances$ID[trip_distances$island=="Ascension"])
among_var<-(test$MS[1]-test$MS[2])/n0
within_var<-test$MS[2]
rep<- among_var / (within_var + among_var)
SE<-  sqrt((2*(sum_ni-1)*(1-rep)*(1-rep)*(1+(n0-1)*rep)*(1+(n0-1)*rep))/(n0*n0*(sum_ni-dim(n_trips)[1])*(dim(n_trips)[1]-1)))
confint<- SE * 3.92
z <- rep/SE
2*pnorm(-abs(z))

repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="direction",3]<- rep
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="direction",4]<- rep-0.5*confint
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="direction",5]<- rep+0.5*confint
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="direction",6]<- 2*pnorm(-abs(z))





n_trips<-aggregate(count~ID, data=trip_distances[trip_distances$island=="StHelena",], FUN=sum)
tab<-table(n_trips$count)
freq<-as.numeric(dimnames(tab)[[1]])
sum_ni<-sum(tab*freq)
sum_ni2<-sum(tab*freq*freq)


n0<-(1/(length(unique(trip_distances$ID[trip_distances$island=="StHelena"]))-1))*(sum(n_trips$count)-(sum(n_trips$count^2)/sum(n_trips$count)))			## mean sample size
test<-aov.circular(circular(trip_distances[trip_distances$island=="StHelena",10],type = "directions", units = "degrees", modulo = "2pi"),trip_distances$ID[trip_distances$island=="StHelena"])
among_var<-(test$MS[1]-test$MS[2])/n0
within_var<-test$MS[2]
rep<- among_var / (within_var + among_var)
SE<-  sqrt((2*(sum_ni-1)*(1-rep)*(1-rep)*(1+(n0-1)*rep)*(1+(n0-1)*rep))/(n0*n0*(sum_ni-dim(n_trips)[1])*(dim(n_trips)[1]-1)))
confint<- SE * 3.92
z <- rep/SE
2*pnorm(-abs(z))

repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="direction",3]<- rep
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="direction",4]<- rep-0.5*confint
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="direction",5]<- rep+0.5*confint
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="direction",6]<- 2*pnorm(-abs(z))


repeatabilities
write.table(repeatabilities, "C:\\STEFFEN\\MANUSCRIPTS\\submitted\\Ascension\\Repeatabilities.csv", sep=",", row.names=F)
write.table(repeatabilities, "A:\\MANUSCRIPTS\\submitted\\Ascension\\Repeatabilities.csv", sep=",", row.names=F)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ESTIMATE REPEATABILITIES WHEN SUBSAMPLING ST HELENA DATA FOR SIMILAR NUMBER OF TRIPS PER IND
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rep_sim<-data.frame()
nsim<-100
ntrips<-aggregate(n_trips~island, data=deployments, FUN=mean)[1,2]
sdtrips<-aggregate(n_trips~island, data=deployments, FUN=sd)[1,2]

for (i in 1:nsim){

xA<-trip_distances[trip_distances$island=="Ascension",]
xS<-trip_distances[trip_distances$island=="StHelena",]
xSsel<-data.frame()

for (l in unique(xS$ID)){

xSind<-subset(xS, ID==l)
N<-as.integer(rnorm(1,ntrips,sdtrips))
N<-ifelse(N>dim(xSind)[1],dim(xSind)[1],N)
picktrips<-sample(xSind$trip,ifelse(N<1,1,N))
xSsel<-rbind(xSsel,xSind[xSind$trip %in% picktrips,])
}

trip_dist_sim<-rbind(xA, xSsel)

repeatabilities<-data.frame(island=rep(c("Ascension","StHelena"), each=4), trip_char=rep(c("time","max_dist","tot_dist","direction"),2), r=0, lcl=0,ucl=0,p=0, simul=i)


r1<-rpt(trip_dist_sim$time[trip_dist_sim$island=="StHelena"], trip_dist_sim$ID[trip_dist_sim$island=="StHelena"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="time",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_dist_sim$time[trip_dist_sim$island=="Ascension"], trip_dist_sim$ID[trip_dist_sim$island=="Ascension"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="time",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_dist_sim$max_dist[trip_dist_sim$island=="StHelena"], trip_dist_sim$ID[trip_dist_sim$island=="StHelena"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="max_dist",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_dist_sim$max_dist[trip_dist_sim$island=="Ascension"], trip_dist_sim$ID[trip_dist_sim$island=="Ascension"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="max_dist",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_dist_sim$total_dist[trip_dist_sim$island=="StHelena"], trip_dist_sim$ID[trip_dist_sim$island=="StHelena"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="tot_dist",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])

r1<-rpt(trip_dist_sim$total_dist[trip_dist_sim$island=="Ascension"], trip_dist_sim$ID[trip_dist_sim$island=="Ascension"], datatype="Gaussian", method="ANOVA")
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="tot_dist",3:6]<-c(r1$R, r1$CI.R[1],r1$CI.R[2], r1$P[2])



## MANUAL CALCULATION OF REPEATABILITY FOLLOWING LESSELS AND BOAG 1987

n_trips<-aggregate(count~ID, data=trip_dist_sim[trip_dist_sim$island=="Ascension",], FUN=sum)
tab<-table(n_trips$count)
freq<-as.numeric(dimnames(tab)[[1]])
sum_ni<-sum(tab*freq)
sum_ni2<-sum(tab*freq*freq)

n0<-(1/(length(unique(trip_dist_sim$ID[trip_dist_sim$island=="Ascension"]))-1))*(sum(n_trips$count)-(sum(n_trips$count^2)/sum(n_trips$count)))			## mean sample size
test<-aov.circular(circular(trip_dist_sim[trip_dist_sim$island=="Ascension",10],type = "directions", units = "degrees", modulo = "2pi"),trip_dist_sim$ID[trip_dist_sim$island=="Ascension"])
among_var<-(test$MS[1]-test$MS[2])/n0
within_var<-test$MS[2]
rep<- among_var / (within_var + among_var)
SE<-  sqrt((2*(sum_ni-1)*(1-rep)*(1-rep)*(1+(n0-1)*rep)*(1+(n0-1)*rep))/(n0*n0*(sum_ni-dim(n_trips)[1])*(dim(n_trips)[1]-1)))
confint<- SE * 3.92
z <- rep/SE
2*pnorm(-abs(z))

repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="direction",3]<- rep
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="direction",4]<- rep-0.5*confint
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="direction",5]<- rep+0.5*confint
repeatabilities[repeatabilities$island=="Ascension" & repeatabilities$trip_char=="direction",6]<- 2*pnorm(-abs(z))



n_trips<-aggregate(count~ID, data=trip_dist_sim[trip_dist_sim$island=="StHelena",], FUN=sum)
tab<-table(n_trips$count)
freq<-as.numeric(dimnames(tab)[[1]])
sum_ni<-sum(tab*freq)
sum_ni2<-sum(tab*freq*freq)


n0<-(1/(length(unique(trip_dist_sim$ID[trip_dist_sim$island=="StHelena"]))-1))*(sum(n_trips$count)-(sum(n_trips$count^2)/sum(n_trips$count)))			## mean sample size
test<-aov.circular(circular(trip_dist_sim[trip_dist_sim$island=="StHelena",10],type = "directions", units = "degrees", modulo = "2pi"),trip_dist_sim$ID[trip_dist_sim$island=="StHelena"])
among_var<-(test$MS[1]-test$MS[2])/n0
within_var<-test$MS[2]
rep<- among_var / (within_var + among_var)
SE<-  sqrt((2*(sum_ni-1)*(1-rep)*(1-rep)*(1+(n0-1)*rep)*(1+(n0-1)*rep))/(n0*n0*(sum_ni-dim(n_trips)[1])*(dim(n_trips)[1]-1)))
confint<- SE * 3.92
z <- rep/SE
2*pnorm(-abs(z))

repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="direction",3]<- rep
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="direction",4]<- rep-0.5*confint
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="direction",5]<- rep+0.5*confint
repeatabilities[repeatabilities$island=="StHelena" & repeatabilities$trip_char=="direction",6]<- 2*pnorm(-abs(z))


rep_sim<-rbind(rep_sim,repeatabilities)

} # end of 100 simulations



##### SUMMARIZE THE SIMULATIONS ####

sim_summary<-aggregate(r ~ trip_char+island, data = rep_sim, FUN=mean)
sim_summary$lcl<-aggregate(lcl ~ trip_char+island, data = rep_sim, FUN=mean)[,3]
sim_summary$ucl<-aggregate(ucl ~ trip_char+island, data = rep_sim, FUN=mean)[,3]
sim_summary$p<-aggregate(p ~ trip_char+island, data = rep_sim, FUN=mean)[,3]
sim_summary




write.table(rep_sim, "C:\\STEFFEN\\MANUSCRIPTS\\submitted\\Ascension\\SimulatedRepeatabilities.csv", sep=",", row.names=F)
write.table(rep_sim, "A:\\MANUSCRIPTS\\submitted\\Ascension\\SimulatedRepeatabilities.csv", sep=",", row.names=F)
write.table(sim_summary, "C:\\STEFFEN\\MANUSCRIPTS\\submitted\\Ascension\\SimulRepeat_SUMMARY.csv", sep=",", row.names=F)
write.table(sim_summary, "A:\\MANUSCRIPTS\\submitted\\Ascension\\SimulRepeat_SUMMARY.csv", sep=",", row.names=F)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARISE THE TIME SPENT IN EACH BEHAVIOUR AND LENGTH OF FORAGING BOUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

trip_distances$time_foraging<-trip_distances$prop_foraging*trip_distances$time
trip_distances$time_travelling<-trip_distances$prop_travelling*trip_distances$time
trip_distances$time_resting<-trip_distances$prop_resting*trip_distances$time

fbouts<-aggregate(forage_bouts~island+Stage, data=trip_distances, FUN=mean)
fbouts$sd<-aggregate(forage_bouts~island+Stage, data=trip_distances, FUN=sd)[,3]
fbouts$forage<-aggregate(prop_foraging~island+Stage, data=trip_distances, FUN=mean)[,3]
fbouts$foragesd<-aggregate(prop_foraging~island+Stage, data=trip_distances, FUN=sd)[,3]
fbouts$travel<-aggregate(prop_travelling~island+Stage, data=trip_distances, FUN=mean)[,3]
fbouts$travelsd<-aggregate(prop_travelling~island+Stage, data=trip_distances, FUN=sd)[,3]
fbouts$rest<-aggregate(prop_resting~island+Stage, data=trip_distances, FUN=mean)[,3]
fbouts$restsd<-aggregate(prop_resting~island+Stage, data=trip_distances, FUN=sd)[,3]

fbouts$foragetime<-aggregate(time_foraging~island+Stage, data=trip_distances, FUN=mean)[,3]
fbouts$foragetimesd<-aggregate(time_foraging~island+Stage, data=trip_distances, FUN=sd)[,3]
fbouts$traveltime<-aggregate(time_travelling~island+Stage, data=trip_distances, FUN=mean)[,3]
fbouts$traveltimesd<-aggregate(time_travelling~island+Stage, data=trip_distances, FUN=sd)[,3]
fbouts$resttime<-aggregate(time_resting~island+Stage, data=trip_distances, FUN=mean)[,3]
fbouts$resttimesd<-aggregate(time_resting~island+Stage, data=trip_distances, FUN=sd)[,3]

fbouts$ON_trips<-aggregate(ON_trip~island+Stage, data=trip_distances, FUN=sum)[,3]
fbouts$bout_length<-aggregate(foragebout_length~island+Stage, data=trip_distances, FUN=mean)[,3]
fbouts$bout_length_sd<-aggregate(foragebout_length~island+Stage, data=trip_distances, FUN=sd)[,3]
fbouts
#write.table(fbouts, "clipboard", sep="\t", row.names=F)


hist(trip_distances$forage_bouts[trip_distances$island=="StHelena"], breaks=c(0,5,10,15,20,25,35,50,70,100,125,150,200), main="", xlab="N of foraging bouts", ylab="% of trips")
hist(trip_distances$forage_bouts[trip_distances$island=="Ascension"], breaks=c(0,5,10,15,20,25,35,50,70,100,125,150,200), col='grey57', axes=F, main="", xlab="", ylab="", add=T)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REVISION 1: calculate energy expenditure for each individual
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


trip_distances<-trip_distances[!trip_distances$ID==204,]
trip_distances$time_foraging<-trip_distances$prop_foraging*trip_distances$time
trip_distances$time_travelling<-trip_distances$prop_travelling*trip_distances$time
trip_distances$time_resting<-trip_distances$prop_resting*trip_distances$time
 

### get body mass data - not available from St Helena!!

### Load data from database
library(RODBC)
setwd("A:/RSPB/UKOT/Ascension/Controlled_Area/Data")
setwd("C:/STEFFEN/RSPB/UKOT/Ascension/Controlled_Area/Data")
db <- odbcConnectAccess2007('Ascension_seabird_tracking.accdb')
ASImass <- sqlQuery(db, "SELECT * FROM weight_change")
odbcClose(db)

### because we do not have individ. mass from St Hel, we use a sex-specific mean
MABOmass<-aggregate(weight~sex, ASImass[ASImass$species=='MABO',], FUN=mean)
deployments$mass<-ifelse(deployments$sex=="male",MABOmass$weight[2]/1000,MABOmass$weight[1]/1000)
aggregate(count~sex+island, deployments, FUN=sum)
deployments$mass[is.na(deployments$mass)]<-MABOmass$weight[1]/1000			### assume that unsexed birds on St Helena were females


enexp<-aggregate(time_foraging~ID, data=trip_distances, FUN=sum)
enexp$time_travelling<-aggregate(time_travelling~ID, data=trip_distances, FUN=sum)[,2]
enexp$time_resting<-aggregate(time_resting~ID, data=trip_distances, FUN=sum)[,2]

enexp<-merge(deployments[,c(1,3,5,9,10,11,12,16,17)],enexp, by="ID")
enexp$time_deployed<-as.numeric(enexp$end-enexp$start)*24
enexp$time_active<-enexp$time_travelling+enexp$time_foraging
enexp$time_rest<-enexp$time_deployed-enexp$time_active
head(enexp)
enexp$sex[is.na(enexp$sex)]<-"female"			### assume that unsexed birds on St Helena were females




### CONVERT ACTIVITY TIMES INTO ENERGY EXPENDITURE

## Ballance 1995: FMR = 1224 kJ/day for smaller RF Booby
## Birt-Friesen 1989: FMR = 4865 kJ/day for Northern gannets; rest = 144 kJ/h, flight = 349 kJ/h
## Weimerskirch 2008: extrapolated FMR = 2300 kJ/day based on regression of the former two studies and body mass
## Boyd et al. 2014: provides equations below for Peruvian Boobies
## therefore we calulate the offset and the hourly energy expenditure for resting and flying based on Boyd et al. 2014

MABO_flight_rate=10^(1.86+(0.748*log10(enexp$mass)))
MABO_rest_rate=10^(1.45+(0.737*log10(enexp$mass)))

enexp$EE_active<-enexp$time_active*MABO_flight_rate
enexp$EE_rest<-enexp$time_rest*MABO_rest_rate
enexp$total_EE<-enexp$EE_active+enexp$EE_rest
enexp$DEE<-(enexp$total_EE/enexp$time_deployed)*24


dee.m0<-lmer(DEE~breeding_status+sex+(1|season), data=enexp)
dee.m1<-lmer(DEE~breeding_status+sex+island+(1|season), data=enexp)

#Likelihood ratio test to check whether m1 is better than m0#
anova(dee.m0,dee.m1, test = "Chisq")

plotdat<-aggregate(DEE~island+sex+breeding_status, enexp, FUN=mean)
plotdat$SD<-aggregate(DEE~island+sex+breeding_status, enexp, FUN=sd)[,4]
plotdat
plotdat$SD[is.na(plotdat$SD)]<-mean(plotdat$SD, na.rm = T)


pdf("MABO_energy_expend.pdf", width=8, height=7)
par(mar=c(3,6,0,0))
errbar(c(0.9,1.1,1.9,2.1), plotdat$DEE, plotdat$DEE+0.5*plotdat$SD, plotdat$DEE-0.5*plotdat$SD, pch=16, col=c(1,2,1,2), xlim=c(0,3), ylim=c(1000,1400), ylab="Daily energy expenditure (kJ)", xlab="", axes=F, cex=1.5, cex.lab=1.3, mgp=c(4,0.5,0))
legend('topright', pch=16, col=c(1,2), legend=c("Ascension", "St Helena"), bty="n", cex=1.3)
axis(1, at=seq(0,3,1), labels=c("","chick rearing","incubation",""),cex.axis=1.3)
axis(2, at=seq(1000, 1400,100), labels=T, cex.axis=1.3, las=1)
dev.off()


pdf("Fig3.pdf", width=10, height=7)
par(mar=c(3,6.5,0,0))
errbar(c(0.7,0.8,1.2,1.3,1.7,1.8,2.2,2.3), plotdat$DEE, plotdat$DEE+0.5*plotdat$SD, plotdat$DEE-0.5*plotdat$SD, pch=c(16,1,17,2,16,1,17,2,16,1,17,2), xlim=c(0,3), ylim=c(1000,1400), ylab="Daily energy expenditure (kJ)", xlab="", axes=F, cex=1.5, cex.lab=1.7, mgp=c(4.5,0.5,0))
legend('topright', pch=c(16,1,17,2), legend=c("Ascension female", "St Helena female","Ascension male", "St Helena male"), bty="n", cex=1.7)
axis(1, at=seq(0,3,1), labels=c("","chick rearing","incubation",""),cex.axis=1.7)
axis(2, at=seq(1000, 1400,100), labels=T, cex.axis=1.7, las=1)
dev.off()


### CONFERENCE PRESENTATION ####
jpeg("A:\\MANUSCRIPTS\\in_press\\Ascension\\Fig3.jpg", width=10, height=7, units="in", quality = 100, res=600)

par(mar=c(3,6.5,0,0))
errbar(c(0.7,0.8,1.2,1.3,1.7,1.8,2.2,2.3), plotdat$DEE, plotdat$DEE+0.5*plotdat$SD, plotdat$DEE-0.5*plotdat$SD, pch=c(16,16,17,17,16,16,17,17,16,16,17,17), col=c('darkred','blue','darkred','blue','darkred','blue','darkred','blue','darkred','blue','darkred','blue'), xlim=c(0.5,2.5), ylim=c(1000,1400), ylab="Daily energy expenditure (kJ)", xlab="", axes=F, cex=1.7, cex.lab=1.7, mgp=c(4.5,0.5,0))
legend('topright', pch=c(16,16), legend=c("Ascension", "St Helena"), col=c('darkred','blue'), bty="n", cex=1.7)
axis(1, at=c(0.5,1,2,2.5), labels=c("","chick rearing","incubation",""),cex.axis=1.7, tck=0)
axis(1, at=c(0.75,1.25,1.75,2.25), labels=c("females","males","females","males"),cex.axis=1.5, tck=0, mgp=c(-2,-1.5,0))
abline(v=1.5, lty=2, col='darkgray')
axis(2, at=seq(1000, 1400,100), labels=T, cex.axis=1.7, las=1)
dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMPARE THE NUMBER OF FORAGING INCIDENTS PER TRIP BETWEEN ISLANDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### DOES NUMBER OF FORAGING BOUTS DIFFER BETWEEN ISLANDS? ###

m1<-lmer(forage_bouts~island+Stage+ (1|season/ID), data=trip_distances)
m0<-lmer(forage_bouts~Stage+ (1|season/ID), data=trip_distances)
anova(m1,m0)



### DOES LENGTH OF FORAGING BOUTS DIFFER BETWEEN ISLANDS? ###

m1<-lmer(foragebout_length~island+Stage+ (1|season/ID), data=trip_distances)
m0<-lmer(foragebout_length~Stage+ (1|season/ID), data=trip_distances)
anova(m1,m0)


### DOES PROPORTION OF TIME SPENT IN EACH BEHAVIOUR DIFFER BETWEEN ISLANDS? ###

m1<-lmer(prop_foraging~island+Stage+ (1|season/ID), data=trip_distances)
m0<-lmer(prop_foraging~Stage+ (1|season/ID), data=trip_distances)
anova(m1,m0)

m1<-lmer(prop_travelling~island+Stage+ (1|season/ID), data=trip_distances)
m0<-lmer(prop_travelling~Stage+ (1|season/ID), data=trip_distances)
anova(m1,m0)

m1<-lmer(prop_resting~island+Stage+ (1|season/ID), data=trip_distances)
m0<-lmer(prop_resting~Stage+ (1|season/ID), data=trip_distances)
anova(m1,m0)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RELATE THE LENGTH OF FORAGING BOUTS TO DISTANCE FROM COLONY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(bouts)
par(mfrow=c(2,2))
plot(ColDist~length, data=bouts[bouts$behav=="forage" & bouts$direction=="outward" ,], main="foraging")
plot(ColDist~length, data=bouts[bouts$behav=="travel" & bouts$direction=="outward" ,], main="travelling")
plot(ColDist~length, data=bouts[bouts$behav=="other" & bouts$direction=="outward" ,], main="nocturnal")
plot(ColDist~length, data=bouts[bouts$behav=="resting" & bouts$direction=="outward" ,], main="resting")

outbouts<-subset(bouts, subset=direction=="outward")
m1<-lmer(length~ColDist+island+(1|trip_id), data=outbouts[outbouts$behav=="forage",])
m0<-lmer(length~island+(1|trip_id), data=outbouts[outbouts$behav=="forage",])
anova(m1,m0)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE SIMPLE 95% KERNEL DENSITIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(adehabitatHR)
tracks2<-tracks2[order(tracks2$ID,tracks2$DateTime),]
head(tracks2)
inHR<-SpatialPointsDataFrame(coords=SpatialPoints(tracks2[tracks2$behav=="forage",8:7]),data=tracks2[tracks2$behav=="forage",c(4,1)])
inHR@data$ID<-NULL

### AREA OF HOME RANGE ####
#kernel.area(xy=tracks2[tracks2$behav=="forage",8:7], id=tracks2[,4],h = 6000, grid = 500, same4all = T, kern = "bivnorm", levels = c(50,95),unin = "m", unout = "km2", extent = 1.0)
#ud<-kernelUD(xy=tracks2[tracks2$behav=="forage",8:7], id=tracks2[,4],h = 6000, grid = 500, same4all = T, kern = "bivnorm",extent = 1.0)
ud<-kernelUD(xy=inHR,h = (ScaleOutASI/2)*1000, grid = 500, same4all = F, kern = "bivnorm",extent = 1.0)
image(ud) ## Note that the contours correspond to values of probability density
kernel.area(ud, percent = c(50,95),unin = "m", unout = "km2")

## EEZ of Ascension: 441658 
## EEZ of St Helena: 444916

79084.73/441658
19403.768/444916


## Plotting of the 95 percent home range

ver <- getverticeshr(ud, 95)
plot(ver, col=ver@data$id)

ver50 <- getverticeshr(ud, 50)
plot(ver50, col=ver@data$id)
axis(1,at=c(100000,300000,500000), labels=c('0', '200', '400 km'), cex=2, cex.lab=1.8, cex.axis=1.8)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARISE SAMPLE SIZES FOR PAPER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aggregate(count~island+season, data=deployments, FUN=sum)

