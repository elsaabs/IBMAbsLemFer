library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(dplyr)
library("RColorBrewer")
library(scales)

########################## 20 GRIDS IN FIGURE 2
############ CASES
#load file "scoresCases_ID_Njob.txt" (example with simulation ID=1989, Njob=0)
f<-read.table('~//Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresCases_1989_0.txt',sep= "",dec='.',header=F)
#Side length of lattice (10x10 --> 10)
tailleReseau = 10
#N=number of microsites in lattice
N = tailleReseau*tailleReseau
#Outputs for t=5x10^5. Change "tpsaobserver" from 501 to 0, 201 or 301 for the 3 first time steps in Figure 3.
tpsaobserver = 501
{
  x=rep(1:tailleReseau,tailleReseau)
  y=rep(1:tailleReseau,each = tailleReseau)
  z=f[tpsaobserver,]
  zZ=as.data.frame(matrix(c(x,y,t(z[2:(N+1)])),ncol=3, nrow=N))
  zC=as.data.frame(matrix(c(x,y,t(z[(N+2):(2*N+1)])),ncol=3, nrow=N))
  zD=as.data.frame(matrix(c(x,y,t(z[(2*N+2):(3*N+1)])),ncol=3, nrow=N))
  zM1=as.data.frame(matrix(c(x,y,t(z[(3*N+2):(4*N+1)])),ncol=3, nrow=N))
  zM2=as.data.frame(matrix(c(x,y,t(z[(4*N+2):(5*N+1)])),ncol=3, nrow=N))
  
  
  gM1 <- ggplot(zM1,aes(zM1$V1,zM1$V2,color=zM1$V3))+geom_point(size = 8, shape=15)+
    scale_colour_gradient(low="rosybrown1", high="red4", guide="legend", name="Number per\ncase",
                          limits=c(1,max(f[,(3*N+2):(5*N+1)]))) +
    ggtitle("Bacteria type 1 (M1)") + scale_y_continuous(name = "") +
    scale_x_continuous(name = "") + labs(x = "t=0", y = "Mmut")
  
  gM2 <- ggplot(zM2,aes(zM2$V1,zM2$V2,color=zM2$V3))+geom_point(size = 8, shape=15)+
    scale_colour_gradient(low="paleturquoise", high="navy",guide="legend", name="Number per\ncase",
                          limits=c(1, max(f[,(3*N+2):(5*N+1)]))) +
    ggtitle("Bacteria type 2 (M2)") + scale_y_continuous(name = "") +
    scale_x_continuous(name = "") 
  
  gZ <- ggplot(zZ,aes(zZ$V1,zZ$V2,color=zZ$V3))+geom_point(size = 8, shape=15)+
    scale_colour_gradient(low="yellow", high="red", guide="legend", name="Biomass\nper case (mg)",
                          limits=c(min(f[,2:(N+1)]), max(f[,2:(N+1)]))) +
    ggtitle("Enzyme (Z)") + scale_y_continuous(name = "") +
    scale_x_continuous(name = "")
  
  gC <- ggplot(zC,aes(zC$V1,zC$V2,color=zC$V3))+geom_point(size = 8, shape=15)+
    scale_colour_gradient(low="yellow", high="red", guide="legend", name="Biomass\nper case (mg)",
                          limits=c(min(f[,(N+2):(2*N+1)]), max(f[,(N+2):(2*N+1)]))) +
    ggtitle("SOC (C)") + scale_y_continuous(name = "") +
    scale_x_continuous(name = "")
  
  gD <- ggplot(zD,aes(zD$V1,zD$V2,color=zD$V3))+geom_point(size = 8, shape=15)+
    scale_colour_gradient(low="yellow", high="red", guide="legend", name="Biomass\nper case (mg)",
                          limits=c(min(f[,(2*N+2):(3*N+1)]), max(f[,(2*N+2):(3*N+1)]))) +
    ggtitle("DOC (D)") + scale_y_continuous(name = "") +
    scale_x_continuous(name = "")
  
  
  grid.arrange(gM1,gM2,gZ,gC,gD,ncol=3,top="step 501 - t = 5e+5")
}




########################## 20 GRIDS IN FIGURE 2
############ TOTAUX (example with 20 simulations: ID=1989, Njob=0 to 19)
#load 20 files "scoresTotaux_ID_Njob.txt"
nom.var <- c("temps", "Ztot", "Ctot", "Dtot", "accessDOCtype1", "accessDOCtype1Rel", "accessDOCtype2", "accessDOCtype2Rel", "M1tot", "M2tot")
h1 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_0.txt',sep= "",dec='.',header=F, col.names=nom.var)
h2 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_1.txt',sep= "",dec='.',header=F, col.names=nom.var)
h3 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_2.txt',sep= "",dec='.',header=F, col.names=nom.var)
h4 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_3.txt',sep= "",dec='.',header=F, col.names=nom.var)
h5 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_4.txt',sep= "",dec='.',header=F, col.names=nom.var)
h6 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_5.txt',sep= "",dec='.',header=F, col.names=nom.var)
h7 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_6.txt',sep= "",dec='.',header=F, col.names=nom.var)
h8 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_7.txt',sep= "",dec='.',header=F, col.names=nom.var)
h9 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_8.txt',sep= "",dec='.',header=F, col.names=nom.var)
h10 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_9.txt',sep= "",dec='.',header=F, col.names=nom.var)
h11 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_10.txt',sep= "",dec='.',header=F, col.names=nom.var)
h12 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_11.txt',sep= "",dec='.',header=F, col.names=nom.var)
h13 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_12.txt',sep= "",dec='.',header=F, col.names=nom.var)
h14 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_13.txt',sep= "",dec='.',header=F, col.names=nom.var)
h15 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_14.txt',sep= "",dec='.',header=F, col.names=nom.var)
h16 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_15.txt',sep= "",dec='.',header=F, col.names=nom.var)
h17 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_16.txt',sep= "",dec='.',header=F, col.names=nom.var)
h18 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_17.txt',sep= "",dec='.',header=F, col.names=nom.var)
h19 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_18.txt',sep= "",dec='.',header=F, col.names=nom.var)
h20 <- read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-2/scoresTotaux_1989_19.txt',sep= "",dec='.',header=F, col.names=nom.var)
#change unit for Z,C,D from number of molecules to density (mass in mg)
h1$Ztot <- h1$Ztot*1e-16
h2$Ztot <- h2$Ztot*1e-16
h3$Ztot <- h3$Ztot*1e-16
h4$Ztot <- h4$Ztot*1e-16
h5$Ztot <- h5$Ztot*1e-16
h6$Ztot <- h6$Ztot*1e-16
h7$Ztot <- h7$Ztot*1e-16
h8$Ztot <- h8$Ztot*1e-16
h9$Ztot <- h9$Ztot*1e-16
h10$Ztot <- h10$Ztot*1e-16
h11$Ztot <- h11$Ztot*1e-16
h12$Ztot <- h12$Ztot*1e-16
h13$Ztot <- h13$Ztot*1e-16
h14$Ztot <- h14$Ztot*1e-16
h15$Ztot <- h15$Ztot*1e-16
h16$Ztot <- h16$Ztot*1e-16
h17$Ztot <- h17$Ztot*1e-16
h18$Ztot <- h18$Ztot*1e-16
h19$Ztot <- h19$Ztot*1e-16
h20$Ztot <- h20$Ztot*1e-16
h1$Ctot <- h1$Ctot*1e-16
h2$Ctot <- h2$Ctot*1e-16
h3$Ctot <- h3$Ctot*1e-16
h4$Ctot <- h4$Ctot*1e-16
h5$Ctot <- h5$Ctot*1e-16
h6$Ctot <- h6$Ctot*1e-16
h7$Ctot <- h7$Ctot*1e-16
h8$Ctot <- h8$Ctot*1e-16
h9$Ctot <- h9$Ctot*1e-16
h10$Ctot <- h10$Ctot*1e-16
h11$Ctot <- h11$Ctot*1e-16
h12$Ctot <- h12$Ctot*1e-16
h13$Ctot <- h13$Ctot*1e-16
h14$Ctot <- h14$Ctot*1e-16
h15$Ctot <- h15$Ctot*1e-16
h16$Ctot <- h16$Ctot*1e-16
h17$Ctot <- h17$Ctot*1e-16
h18$Ctot <- h18$Ctot*1e-16
h19$Ctot <- h19$Ctot*1e-16
h20$Ctot <- h20$Ctot*1e-16
h1$Dtot <- h1$Dtot*1e-19
h2$Dtot <- h2$Dtot*1e-19
h3$Dtot <- h3$Dtot*1e-19
h4$Dtot <- h4$Dtot*1e-19
h5$Dtot <- h5$Dtot*1e-19
h6$Dtot <- h6$Dtot*1e-19
h7$Dtot <- h7$Dtot*1e-19
h8$Dtot <- h8$Dtot*1e-19
h9$Dtot <- h9$Dtot*1e-19
h10$Dtot <- h10$Dtot*1e-19
h11$Dtot <- h11$Dtot*1e-19
h12$Dtot <- h12$Dtot*1e-19
h13$Dtot <- h13$Dtot*1e-19
h14$Dtot <- h14$Dtot*1e-19
h15$Dtot <- h15$Dtot*1e-19
h16$Dtot <- h16$Dtot*1e-19
h17$Dtot <- h17$Dtot*1e-19
h18$Dtot <- h18$Dtot*1e-19
h19$Dtot <- h19$Dtot*1e-19
h20$Dtot <- h20$Dtot*1e-19
h1$accessDOCtype1 <- h1$accessDOCtype1*1e-19
h2$accessDOCtype1 <- h2$accessDOCtype1*1e-19
h3$accessDOCtype1 <- h3$accessDOCtype1*1e-19
h4$accessDOCtype1 <- h4$accessDOCtype1*1e-19
h5$accessDOCtype1 <- h5$accessDOCtype1*1e-19
h6$accessDOCtype1 <- h6$accessDOCtype1*1e-19
h7$accessDOCtype1 <- h7$accessDOCtype1*1e-19
h8$accessDOCtype1 <- h8$accessDOCtype1*1e-19
h9$accessDOCtype1 <- h9$accessDOCtype1*1e-19
h10$accessDOCtype1 <- h10$accessDOCtype1*1e-19
h11$accessDOCtype1 <- h11$accessDOCtype1*1e-19
h12$accessDOCtype1 <- h12$accessDOCtype1*1e-19
h13$accessDOCtype1 <- h13$accessDOCtype1*1e-19
h14$accessDOCtype1 <- h14$accessDOCtype1*1e-19
h15$accessDOCtype1 <- h15$accessDOCtype1*1e-19
h16$accessDOCtype1 <- h16$accessDOCtype1*1e-19
h17$accessDOCtype1 <- h17$accessDOCtype1*1e-19
h18$accessDOCtype1 <- h18$accessDOCtype1*1e-19
h19$accessDOCtype1 <- h19$accessDOCtype1*1e-19
h20$accessDOCtype1 <- h20$accessDOCtype1*1e-19
h1$accessDOCtype1Rel <- h1$accessDOCtype1Rel*1e-19
h2$accessDOCtype1Rel <- h2$accessDOCtype1Rel*1e-19
h3$accessDOCtype1Rel <- h3$accessDOCtype1Rel*1e-19
h4$accessDOCtype1Rel <- h4$accessDOCtype1Rel*1e-19
h5$accessDOCtype1Rel <- h5$accessDOCtype1Rel*1e-19
h6$accessDOCtype1Rel <- h6$accessDOCtype1Rel*1e-19
h7$accessDOCtype1Rel <- h7$accessDOCtype1Rel*1e-19
h8$accessDOCtype1Rel <- h8$accessDOCtype1Rel*1e-19
h9$accessDOCtype1Rel <- h9$accessDOCtype1Rel*1e-19
h10$accessDOCtype1Rel <- h10$accessDOCtype1Rel*1e-19
h11$accessDOCtype1Rel <- h11$accessDOCtype1Rel*1e-19
h12$accessDOCtype1Rel <- h12$accessDOCtype1Rel*1e-19
h13$accessDOCtype1Rel <- h13$accessDOCtype1Rel*1e-19
h14$accessDOCtype1Rel <- h14$accessDOCtype1Rel*1e-19
h15$accessDOCtype1Rel <- h15$accessDOCtype1Rel*1e-19
h16$accessDOCtype1Rel <- h16$accessDOCtype1Rel*1e-19
h17$accessDOCtype1Rel <- h17$accessDOCtype1Rel*1e-19
h18$accessDOCtype1Rel <- h18$accessDOCtype1Rel*1e-19
h19$accessDOCtype1Rel <- h19$accessDOCtype1Rel*1e-19
h20$accessDOCtype1Rel <- h20$accessDOCtype1Rel*1e-19
h1$accessDOCtype2 <- h1$accessDOCtype2*1e-19
h2$accessDOCtype2 <- h2$accessDOCtype2*1e-19
h3$accessDOCtype2 <- h3$accessDOCtype2*1e-19
h4$accessDOCtype2 <- h4$accessDOCtype2*1e-19
h5$accessDOCtype2 <- h5$accessDOCtype2*1e-19
h6$accessDOCtype2 <- h6$accessDOCtype2*1e-19
h7$accessDOCtype2 <- h7$accessDOCtype2*1e-19
h8$accessDOCtype2 <- h8$accessDOCtype2*1e-19
h9$accessDOCtype2 <- h9$accessDOCtype2*1e-19
h10$accessDOCtype2 <- h10$accessDOCtype2*1e-19
h11$accessDOCtype2 <- h11$accessDOCtype2*1e-19
h12$accessDOCtype2 <- h12$accessDOCtype2*1e-19
h13$accessDOCtype2 <- h13$accessDOCtype2*1e-19
h14$accessDOCtype2 <- h14$accessDOCtype2*1e-19
h15$accessDOCtype2 <- h15$accessDOCtype2*1e-19
h16$accessDOCtype2 <- h16$accessDOCtype2*1e-19
h17$accessDOCtype2 <- h17$accessDOCtype2*1e-19
h18$accessDOCtype2 <- h18$accessDOCtype2*1e-19
h19$accessDOCtype2 <- h19$accessDOCtype2*1e-19
h20$accessDOCtype2 <- h20$accessDOCtype2*1e-19
h1$accessDOCtype2Rel <- h1$accessDOCtype2Rel*1e-19
h2$accessDOCtype2Rel <- h2$accessDOCtype2Rel*1e-19
h3$accessDOCtype2Rel <- h3$accessDOCtype2Rel*1e-19
h4$accessDOCtype2Rel <- h4$accessDOCtype2Rel*1e-19
h5$accessDOCtype2Rel <- h5$accessDOCtype2Rel*1e-19
h6$accessDOCtype2Rel <- h6$accessDOCtype2Rel*1e-19
h7$accessDOCtype2Rel <- h7$accessDOCtype2Rel*1e-19
h8$accessDOCtype2Rel <- h8$accessDOCtype2Rel*1e-19
h9$accessDOCtype2Rel <- h9$accessDOCtype2Rel*1e-19
h10$accessDOCtype2Rel <- h10$accessDOCtype2Rel*1e-19
h11$accessDOCtype2Rel <- h11$accessDOCtype2Rel*1e-19
h12$accessDOCtype2Rel <- h12$accessDOCtype2Rel*1e-19
h13$accessDOCtype2Rel <- h13$accessDOCtype2Rel*1e-19
h14$accessDOCtype2Rel <- h14$accessDOCtype2Rel*1e-19
h15$accessDOCtype2Rel <- h15$accessDOCtype2Rel*1e-19
h16$accessDOCtype2Rel <- h16$accessDOCtype2Rel*1e-19
h17$accessDOCtype2Rel <- h17$accessDOCtype2Rel*1e-19
h18$accessDOCtype2Rel <- h18$accessDOCtype2Rel*1e-19
h19$accessDOCtype2Rel <- h19$accessDOCtype2Rel*1e-19
h20$accessDOCtype2Rel <- h20$accessDOCtype2Rel*1e-19
k1 <- h1[-nrow(h1),]
k2 <- h2[-nrow(h2),]
k3 <- h3[-nrow(h3),]
k4 <- h4[-nrow(h4),]
k5 <- h5[-nrow(h5),]
k6 <- h6[-nrow(h6),]
k7 <- h7[-nrow(h7),]
k8 <- h8[-nrow(h8),]
k9 <- h9[-nrow(h9),]
k10 <- h10[-nrow(h10),]
k11 <- h11[-nrow(h11),]
k12 <- h12[-nrow(h12),]
k13 <- h13[-nrow(h13),]
k14 <- h14[-nrow(h14),]
k15 <- h15[-nrow(h15),]
k16 <- h16[-nrow(h16),]
k17 <- h17[-nrow(h17),]
k18 <- h18[-nrow(h18),]
k19 <- h19[-nrow(h19),]
k20 <- h20[-nrow(h20),]
h <- rbind(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20)
h <- cbind(h, nbrRep=1)
h <- cbind(h, tempsround=floor(h$temps))
h <- cbind(h, carreM1tot=h$M1tot*h$M1tot)
h <- cbind(h, carreM2tot=h$M2tot*h$M2tot)
g <- aggregate(list(sumZtot = h$Ztot, sumCtot = h$Ctot, sumDtot = h$Dtot, 
                    sumaccessDOCtype1 = h$accessDOCtype1, sumaccessDOCtype1Rel = h$accessDOCtype1Rel, 
                    sumaccessDOCtype2 = h$accessDOCtype2, sumaccessDOCtype2Rel = h$accessDOCtype2Rel, 
                    sumM1tot= h$M1tot, sumM2tot = h$M2tot, sumcarreM1 = h$carreM1tot,
                    sumcarreM2 = h$carreM2tot, sumnbrRep = h$nbrRep), 
               by = list(tempsround = h$tempsround), sum)
g$meanZtot <- g$sumZtot / g$sumnbrRep
g$meanCtot <- g$sumCtot / g$sumnbrRep
g$meanDtot <- g$sumDtot / g$sumnbrRep
g$meanaccessDOCtype1 <- g$sumaccessDOCtype1 / g$sumnbrRep
g$meanaccessDOCtype1Rel <- g$sumaccessDOCtype1Rel / g$sumnbrRep
g$meanaccessDOCtype2 <- g$sumaccessDOCtype2 / g$sumnbrRep
g$meanaccessDOCtype2Rel<- g$sumaccessDOCtype2Rel / g$sumnbrRep
g$meanM1tot<- g$sumM1tot / g$sumnbrRep
g$meanM2tot <- g$sumM2tot / g$sumnbrRep
g$sdM1tot <- sqrt(g$sumcarreM1 / g$sumnbrRep - g$meanM1tot*g$meanM1tot)
g$sdM2tot <- sqrt(g$sumcarreM2 / g$sumnbrRep - g$meanM2tot*g$meanM2tot)
g$varM1tot <- g$sumcarreM1 / g$sumnbrRep - g$meanM1tot*g$meanM1tot
g$varM2tot <- g$sumcarreM2 / g$sumnbrRep - g$meanM2tot*g$meanM2tot
#plots MEAN (last right right column in figure 2)
par(mfrow = c(2,3))
plot(g$tempsround,g$meanM1tot, 
     xlim=c(0,1e+6),
     ylim=c(min(h1$M1tot,h2$M1tot,h3$M1tot,h4$M1tot,h5$M1tot,
                h6$M1tot,h7$M1tot,h8$M1tot,h9$M1tot,h10$M1tot,
                h11$M1tot,h12$M1tot,h13$M1tot,h14$M1tot,h15$M1tot,
                h16$M1tot,h17$M1tot,h18$M1tot,h19$M1tot,h20$M1tot),
            max(h1$M1tot,h2$M1tot,h3$M1tot,h4$M1tot,h5$M1tot,
                h6$M1tot,h7$M1tot,h8$M1tot,h9$M1tot,h10$M1tot,
                h11$M1tot,h12$M1tot,h13$M1tot,h14$M1tot,h15$M1tot,
                h16$M1tot,h17$M1tot,h18$M1tot,h19$M1tot,h20$M1tot)), 
     type='l',col="red",xlab = "Time", ylab = "meanM1tot")
title("meanM1tot")
#min and max for range around mean:
#lines(cbind(g$tempsround,NA,g$tempsround,NA),cbind(g$meanM1tot-g$sdM1tot,NA,g$meanM1tot+g$sdM1tot,NA))
plot(g$tempsround,g$meanM2tot, 
     xlim=c(0,1e+6),
     ylim=c(min(h1$M2tot,h2$M2tot,h3$M2tot,h4$M2tot,h5$M2tot,
                h6$M2tot,h7$M2tot,h8$M2tot,h9$M2tot,h10$M2tot,
                h11$M2tot,h12$M2tot,h13$M2tot,h14$M2tot,h15$M2tot,
                h16$M2tot,h17$M2tot,h18$M2tot,h19$M2tot,h20$M2tot),
            max(h1$M2tot,h2$M2tot,h3$M2tot,h4$M2tot,h5$M2tot,
                h6$M2tot,h7$M2tot,h8$M2tot,h9$M2tot,h10$M2tot,
                h11$M2tot,h12$M2tot,h13$M2tot,h14$M2tot,h15$M2tot,
                h16$M2tot,h17$M2tot,h18$M2tot,h19$M2tot,h20$M2tot)), 
     type='l',col="red",xlab = "Time", ylab = "meanM2tot")
title("meanM2tot")
#min and max for range around mean:
#lines(cbind(g$tempsround,NA,g$tempsround,NA),cbind(g$meanM2tot-g$sdM2tot,NA,g$meanM2tot+g$sdM2tot,NA))
plot(g$tempsround,g$meanZtot, 
     xlim=c(0,1e+6),
     ylim=c(min(h1$Ztot,h2$Ztot,h3$Ztot,h4$Ztot,h5$Ztot,
                h6$Ztot,h7$Ztot,h8$Ztot,h9$Ztot,h10$Ztot,
                h11$Ztot,h12$Ztot,h13$Ztot,h14$Ztot,h15$Ztot,
                h16$Ztot,h17$Ztot,h18$Ztot,h19$Ztot,h20$Ztot),
            max(h1$Ztot,h2$Ztot,h3$Ztot,h4$Ztot,h5$Ztot,
                h6$Ztot,h7$Ztot,h8$Ztot,h9$Ztot,h10$Ztot,
                h11$Ztot,h12$Ztot,h13$Ztot,h14$Ztot,h15$Ztot,
                h16$Ztot,h17$Ztot,h18$Ztot,h19$Ztot,h20$Ztot)), 
     type='l',col="red",xlab = "Time", ylab = "meanZtot")
title("meanZtot")
plot(g$tempsround,g$meanCtot, 
     xlim=c(0,1e+6),
     ylim=c(min(h1$Ctot,h2$Ctot,h3$Ctot,h4$Ctot,h5$Ctot,
                h6$Ctot,h7$Ctot,h8$Ctot,h9$Ctot,h10$Ctot,
                h11$Ctot,h12$Ctot,h13$Ctot,h14$Ctot,h15$Ctot,
                h16$Ctot,h17$Ctot,h18$Ctot,h19$Ctot,h20$Ctot),
            max(h1$Ctot,h2$Ctot,h3$Ctot,h4$Ctot,h5$Ctot,
                h6$Ctot,h7$Ctot,h8$Ctot,h9$Ctot,h10$Ctot,
                h11$Ctot,h12$Ctot,h13$Ctot,h14$Ctot,h15$Ctot,
                h16$Ctot,h17$Ctot,h18$Ctot,h19$Ctot,h20$Ctot)), 
     type='l',col="red",xlab = "Time", ylab = "meanCtot")
title("meanCtot")
plot(g$tempsround,g$meanDtot, 
     xlim=c(0,1e+6),
     ylim=c(min(h1$Dtot,h2$Dtot,h3$Dtot,h4$Dtot,h5$Dtot,
                h6$Dtot,h7$Dtot,h8$Dtot,h9$Dtot,h10$Dtot,
                h11$Dtot,h12$Dtot,h13$Dtot,h14$Dtot,h15$Dtot,
                h16$Dtot,h17$Dtot,h18$Dtot,h19$Dtot,h20$Dtot),
            max(h1$Dtot,h2$Dtot,h3$Dtot,h4$Dtot,h5$Dtot,
                h6$Dtot,h7$Dtot,h8$Dtot,h9$Dtot,h10$Dtot,
                h11$Dtot,h12$Dtot,h13$Dtot,h14$Dtot,h15$Dtot,
                h16$Dtot,h17$Dtot,h18$Dtot,h19$Dtot,h20$Dtot)), 
     type='l',col="red",xlab = "Time", ylab = "meanDtot")
title("meanDtot")
#plots 1 job (Njob=0) dynamics (before last right right column in figure 2)
par(mfrow = c(2,3))
plot(h1$temps,h1$M1tot, 
    xlim=c(0,5e+5),
    ylim=c(min(h1$M1tot,h2$M1tot,h3$M1tot,h4$M1tot,h5$M1tot,
               h6$M1tot,h7$M1tot,h8$M1tot,h9$M1tot,h10$M1tot,
               h11$M1tot,h12$M1tot,h13$M1tot,h14$M1tot,h15$M1tot,
               h16$M1tot,h17$M1tot,h18$M1tot,h19$M1tot,h20$M1tot),
           max(h1$M1tot,h2$M1tot,h3$M1tot,h4$M1tot,h5$M1tot,
               h6$M1tot,h7$M1tot,h8$M1tot,h9$M1tot,h10$M1tot,
               h11$M1tot,h12$M1tot,h13$M1tot,h14$M1tot,h15$M1tot,
               h16$M1tot,h17$M1tot,h18$M1tot,h19$M1tot,h20$M1tot)), 
    type='l',col="black",xlab = "Time")
plot(h1$temps,h1$M2tot, 
     xlim=c(0,5e+5),
     ylim=c(min(h1$M2tot,h2$M2tot,h3$M2tot,h4$M2tot,h5$M2tot,
                h6$M2tot,h7$M2tot,h8$M2tot,h9$M2tot,h10$M2tot,
                h11$M2tot,h12$M2tot,h13$M2tot,h14$M2tot,h15$M2tot,
                h16$M2tot,h17$M2tot,h18$M2tot,h19$M2tot,h20$M2tot),
            max(h1$M2tot,h2$M2tot,h3$M2tot,h4$M2tot,h5$M2tot,
                h6$M2tot,h7$M2tot,h8$M2tot,h9$M2tot,h10$M2tot,
                h11$M2tot,h12$M2tot,h13$M2tot,h14$M2tot,h15$M2tot,
                h16$M2tot,h17$M2tot,h18$M2tot,h19$M2tot,h20$M2tot)),
     type='l',col="black",xlab = "Time")
plot(h1$temps,h1$Ztot, 
     xlim=c(0,5e+5),
     ylim=c(min(h1$Ztot,h2$Ztot,h3$Ztot,h4$Ztot,h5$Ztot,
                h6$Ztot,h7$Ztot,h8$Ztot,h9$Ztot,h10$Ztot,
                h11$Ztot,h12$Ztot,h13$Ztot,h14$Ztot,h15$Ztot,
                h16$Ztot,h17$Ztot,h18$Ztot,h19$Ztot,h20$Ztot),
            max(h1$Ztot,h2$Ztot,h3$Ztot,h4$Ztot,h5$Ztot,
                h6$Ztot,h7$Ztot,h8$Ztot,h9$Ztot,h10$Ztot,
                h11$Ztot,h12$Ztot,h13$Ztot,h14$Ztot,h15$Ztot,
                h16$Ztot,h17$Ztot,h18$Ztot,h19$Ztot,h20$Ztot)),
     type='l',col="black",xlab = "Time")
plot(h1$temps,h1$Ctot, 
     xlim=c(0,5e+5),
     ylim=c(min(h1$Ctot,h2$Ctot,h3$Ctot,h4$Ctot,h5$Ctot,
                h6$Ctot,h7$Ctot,h8$Ctot,h9$Ctot,h10$Ctot,
                h11$Ctot,h12$Ctot,h13$Ctot,h14$Ctot,h15$Ctot,
                h16$Ctot,h17$Ctot,h18$Ctot,h19$Ctot,h20$Ctot),
            max(h1$Ctot,h2$Ctot,h3$Ctot,h4$Ctot,h5$Ctot,
                h6$Ctot,h7$Ctot,h8$Ctot,h9$Ctot,h10$Ctot,
                h11$Ctot,h12$Ctot,h13$Ctot,h14$Ctot,h15$Ctot,
                h16$Ctot,h17$Ctot,h18$Ctot,h19$Ctot,h20$Ctot)),
     type='l',col="black",xlab = "Time")
plot(h1$temps,h1$Dtot, 
     xlim=c(0,5e+5),
     ylim=c(min(h1$Dtot,h2$Dtot,h3$Dtot,h4$Dtot,h5$Dtot,
                h6$Dtot,h7$Dtot,h8$Dtot,h9$Dtot,h10$Dtot,
                h11$Dtot,h12$Dtot,h13$Dtot,h14$Dtot,h15$Dtot,
                h16$Dtot,h17$Dtot,h18$Dtot,h19$Dtot,h20$Dtot),
            max(h1$Dtot,h2$Dtot,h3$Dtot,h4$Dtot,h5$Dtot,
                h6$Dtot,h7$Dtot,h8$Dtot,h9$Dtot,h10$Dtot,
                h11$Dtot,h12$Dtot,h13$Dtot,h14$Dtot,h15$Dtot,
                h16$Dtot,h17$Dtot,h18$Dtot,h19$Dtot,h20$Dtot)),
     type='l',col="black",xlab = "Time")


