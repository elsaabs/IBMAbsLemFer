library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(dplyr)
library("RColorBrewer")
library(scales)

########################## PROBA MUTATION = 0
############# POPULATIONS DYNAMICS
#load 20 files "scoresTotaux_ID_Njob.txt"
nom.var <- c("temps","Ztot","Ctot","Dtot","Mtot","totalNbrTrait","phi","nbrMphi")
h1 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_60.txt',sep= "",dec='.',header=F, col.names=nom.var)
h2 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_61.txt',sep= "",dec='.',header=F, col.names=nom.var)
h3 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_62.txt',sep= "",dec='.',header=F, col.names=nom.var)
h4 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_63.txt',sep= "",dec='.',header=F, col.names=nom.var)
h5 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_64.txt',sep= "",dec='.',header=F, col.names=nom.var)
h6 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_65.txt',sep= "",dec='.',header=F, col.names=nom.var)
h7 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_66.txt',sep= "",dec='.',header=F, col.names=nom.var)
h8 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_67.txt',sep= "",dec='.',header=F, col.names=nom.var)
h9 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_68.txt',sep= "",dec='.',header=F, col.names=nom.var)
h10 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_69.txt',sep= "",dec='.',header=F, col.names=nom.var)
h11 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_70.txt',sep= "",dec='.',header=F, col.names=nom.var)
h12 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_71.txt',sep= "",dec='.',header=F, col.names=nom.var)
h13 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_72.txt',sep= "",dec='.',header=F, col.names=nom.var)
h14 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_73.txt',sep= "",dec='.',header=F, col.names=nom.var)
h15 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_74.txt',sep= "",dec='.',header=F, col.names=nom.var)
h16 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_75.txt',sep= "",dec='.',header=F, col.names=nom.var)
h17 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_76.txt',sep= "",dec='.',header=F, col.names=nom.var)
h18 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_77.txt',sep= "",dec='.',header=F, col.names=nom.var)
h19 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_78.txt',sep= "",dec='.',header=F, col.names=nom.var)
h20 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1966_79.txt',sep= "",dec='.',header=F, col.names=nom.var)
mycols <- colors()[c(7,8,26,30,31,32,33,37,42,47,51,52,56,62,75,84,86,92,100,107)]
plot(h1$temps,h1$Mtot, 
     xlim=c(0,1.2e+7),
     ylim=c(min(h1$Mtot,h2$Mtot,h3$Mtot,h4$Mtot,h5$Mtot,
                h6$Mtot,h7$Mtot,h8$Mtot,h9$Mtot,h10$Mtot,
                h11$Mtot,h12$Mtot,h13$Mtot,h14$Mtot,h15$Mtot,
                h16$Mtot,h17$Mtot,h18$Mtot,h19$Mtot,h20$Mtot),
            28), 
     type='l',col=mycols[1],xlab = "Time (hours)", ylab = "Cell number, M", cex=1.5)
lines(h2$temps,h2$Mtot, col=mycols[2])
lines(h3$temps,h3$Mtot, col=mycols[3])
lines(h4$temps,h4$Mtot, col=mycols[4])
lines(h5$temps,h5$Mtot, col=mycols[5])
lines(h6$temps,h6$Mtot, col=mycols[6])
lines(h7$temps,h7$Mtot, col=mycols[7])
lines(h8$temps,h8$Mtot, col=mycols[8])
lines(h9$temps,h9$Mtot, col=mycols[9])
lines(h10$temps,h10$Mtot, col=mycols[10])
lines(h11$temps,h11$Mtot, col=mycols[11])
lines(h12$temps,h12$Mtot, col=mycols[12])
lines(h13$temps,h13$Mtot, col=mycols[13])
lines(h14$temps,h14$Mtot, col=mycols[14])
lines(h15$temps,h15$Mtot, col=mycols[15])
lines(h16$temps,h16$Mtot, col=mycols[16])
lines(h17$temps,h17$Mtot, col=mycols[17])
lines(h18$temps,h18$Mtot, col=mycols[18])
lines(h19$temps,h19$Mtot, col=mycols[19])
lines(h20$temps,h20$Mtot, col=mycols[20])
title("Mtot")

########################## PROBA MUTATION POSITIVE
############# POPULATIONS DYNAMICS
#load 10 files "scoresTotaux_ID_Njob.txt"
nom.var <- c("temps","Ztot","Ctot","Dtot","Mtot","totalNbrTrait","phi","nbrMphi")
h1 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_30.txt',sep= "",dec='.',header=F, col.names=nom.var)
h2 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_31.txt',sep= "",dec='.',header=F, col.names=nom.var)
h3 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_32.txt',sep= "",dec='.',header=F, col.names=nom.var)
h4 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_33.txt',sep= "",dec='.',header=F, col.names=nom.var)
h5 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_34.txt',sep= "",dec='.',header=F, col.names=nom.var)
h6 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_35.txt',sep= "",dec='.',header=F, col.names=nom.var)
h7 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_36.txt',sep= "",dec='.',header=F, col.names=nom.var)
h8 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_37.txt',sep= "",dec='.',header=F, col.names=nom.var)
h9 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_38.txt',sep= "",dec='.',header=F, col.names=nom.var)
h10 <- read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_39.txt',sep= "",dec='.',header=F, col.names=nom.var)
mycols <- colors()[c(7,8,26,30,31,32,33,37,42,47)]
plot(h1$temps,h1$Mtot, 
     xlim=c(0,1.2e+7),
     ylim=c(min(h1$Mtot,h2$Mtot,h3$Mtot,h4$Mtot,h5$Mtot,
                h6$Mtot,h7$Mtot,h8$Mtot,h9$Mtot,h10$Mtot),
            max(h1$Mtot,h2$Mtot,h3$Mtot,h4$Mtot,h5$Mtot,
                h6$Mtot,h7$Mtot,h8$Mtot,h9$Mtot,h10$Mtot)), 
     type='l',col=mycols[1],xlab = "Time (hours)", ylab = "Cell number, M", cex=1.5)
lines(h2$temps,h2$Mtot, col=mycols[2])
lines(h3$temps,h3$Mtot, col=mycols[3])
lines(h4$temps,h4$Mtot, col=mycols[4])
lines(h5$temps,h5$Mtot, col=mycols[5])
lines(h6$temps,h6$Mtot, col=mycols[6])
lines(h7$temps,h7$Mtot, col=mycols[7])
lines(h8$temps,h8$Mtot, col=mycols[8])
lines(h9$temps,h9$Mtot, col=mycols[9])
lines(h10$temps,h10$Mtot, col=mycols[10])
title("Mtot")

############# TRAIT DYNAMICS
#load 1 file "scoresTotaux_ID_Njob.txt"
nom.var <- c("temps","Ztot","Ctot","Dtot","Mtot","totalNbrTrait","phi","nbrMphi")
f<-read.table('~/Desktop/non-spatial-PDMP/plots-Figure-1/scoresTotaux_1959_39.txt',sep= "",dec='.',header=F, col.names=nom.var)
ggplot(f,aes(x=temps,y=phi,color=log(nbrMphi),size=nbrMphi)) + geom_point(shape=19) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) +
  labs(x = "Time (hours)", y = expression(paste("Exoenzyme allocation trait ", (phi))))+
  scale_colour_gradient2(low="black", mid="red", high="orange", 
                         midpoint=(log(min(f$nbrMphi))+log(max(f$nbrMphi)))/2, 
                         guide="legend",limits=c(log(min(f$nbrMphi)),log(max(f$nbrMphi)))) +
  scale_size_continuous(limits=c(min(f$nbrMphi),max(f$nbrMphi)),range=c(1,8)) +
  theme(legend.direction = "horizontal", legend.position = c(0.8, 0.8)) +
  labs(size="Cell Number", colour="Log(Cell Number)") +
  ylim(min(f$phi),max(f$phi))




