library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(dplyr)
library("RColorBrewer")
library(scales)
library(extrafont)
font_import()
loadfonts()




########################## FIGURE 4
############# COMMUNS
#load files "scoresCommuns_ID.txt"
nom.var = c("IDcluster","IDjob","timeMax","T0Mean","phi","diffD",
            "VmU","KmU","dM","dZ","gM","gZ","p","VmD","IC","ID","lC","lD",
            "timeIndex","popTotFin","meanZ","txDecomp","varZ","meanC","varC",
            "meanD","varD","meanM","varM","meanAccessD","varAccessD")
#With evolution
fcom1 = read.table('~/Desktop/spatial-PDMP-one-strain/plots-Figure-4/scoresCommuns_3707.txt',sep= "",dec='.',header=F, col.names=nom.var)
#Mean, min and max of decomposition rate and carbon stock out of the 20 simulations
fcom1 = cbind(fcom1, nbrRep=1)
reptab1 = aggregate(list(sumnbrRep1 = fcom1$nbrRep), 
                   by = list(phi = fcom1$phi, diffD = fcom1$diffD), sum, na.rm=T)
mtab1 = aggregate(list(meanDecompEvol = fcom1$txDecomp, meanCEvol = fcom1$meanC), 
                 by = list(phi = fcom1$phi, Diff = fcom1$diffD), mean, na.rm=T)
mintab1 = aggregate(list(mintxDecomp1 = fcom1$txDecomp, minmeanC1 = fcom1$meanC), 
                   by = list(phi = fcom1$phi, Diff = fcom1$diffD), min, na.rm=T)
maxtab1 = aggregate(list(maxtxDecomp1 = fcom1$txDecomp, maxmeanC1 = fcom1$meanC), 
                   by = list(phi = fcom1$phi, Diff = fcom1$diffD), max, na.rm=T)
gcom1 = cbind(mtab1,minDecompEvol=mintab1$mintxDecomp1,maxDecompEvol=maxtab1$maxtxDecomp1,minCEvol=mintab1$minmeanC1,maxCEvol=maxtab1$maxmeanC1,sumnbrRepEvol=reptab1$sumnbrRep1)
#Without evolution
fcom2 = read.table('~/Desktop/spatial-PDMP-one-strain/plots-Figure-4/scoresCommuns_3708.txt',sep= "",dec='.',header=F, col.names=nom.var)
#Mean, min and max of decomposition rate and carbon stock out of the 20 simulations
fcom2 = cbind(fcom2, nbrRep=1)
reptab2 = aggregate(list(sumnbrRep2 = fcom2$nbrRep), 
                    by = list(phi = fcom2$phi, diffD = fcom2$diffD), sum, na.rm=T)
mtab2 = aggregate(list(meanDecompEcos = fcom2$txDecomp, meanCEcos = fcom2$meanC), 
                  by = list(phi = fcom2$phi, Diff = fcom2$diffD), mean, na.rm=T)
mintab2 = aggregate(list(mintxDecomp2 = fcom2$txDecomp, minmeanC2 = fcom2$meanC), 
                    by = list(phi = fcom2$phi, Diff = fcom2$diffD), min, na.rm=T)
maxtab2 = aggregate(list(maxtxDecomp2 = fcom2$txDecomp, maxmeanC2 = fcom2$meanC), 
                    by = list(phi = fcom2$phi, Diff = fcom2$diffD), max, na.rm=T)
gcom2 = cbind(mtab2,minDecompEcos=mintab2$mintxDecomp2,maxDecompEcos=maxtab2$maxtxDecomp2,minCEcos=mintab2$minmeanC2,maxCEcos=maxtab2$maxmeanC2,sumnbrRepEcos=reptab2$sumnbrRep2)
#Figure 5a
phiESS = c(0.2,0.175,0.175,0.125,0.125,0.125,0.125,0.125,0.125,
           0.1,0.1,0.1,0.1,0.1,0.1,0.075,0.075,0.075,0.075,0.075)
diffESS = seq(5e-7, 1e-5, by=5e-7)
dffig5a = data.frame(diffESS, phiESS)
#pdf("fig5a.pdf", family="Arial Narrow", width=6, height=4)
plot(diffESS, rep(0.075, 20), type='l', lwd = 5, col=alpha(rgb(0,0,1), 0.5), frame = FALSE,
     xlim=range(c(1e-7, 1e-5)),
     ylim=c(0.05,0.2),
     xlab = expression(paste("Diffusion Rate ", (cm^{2} ~ h^{-1}))),
     ylab = expression(paste("Exoenzyme allocation ESS, ", phi,"*")))
lines(diffESS, phiESS, type='p', pch = 19, col = "red", 
      lty = 2, lwd = 1)
#dev.off()
#Figure 5b
#pdf("fig5b.pdf", family="Arial Narrow", width=6, height=4)
plot(gcom2$Diff, gcom2$meanDecompEcos,
     xlim=range(c(1e-7, 1e-5)),
     ylim=range(c(gcom1$minDecompEvol, 0.0002)),
     type='l', lwd = 5, col=alpha(rgb(0,0,1), 0.5), frame = FALSE, pch=19,
     xlab=expression(paste("Diffusion Rate ", (cm^{2} ~ h^{-1}))), 
     ylab=expression(paste("Decomposition Rate ", (h^{-1}))))
lines(gcom1$Diff, gcom1$meanDecompEvol, type='p', pch = 19, col = "red")
# hack: we draw arrows but with very special "arrowheads"
arrows(gcom1$Diff, gcom1$minDecompEvol, gcom1$Diff, gcom1$maxDecompEvol, 
       length=0.05, angle=90, code=3)
#dev.off()
#Figure 5c
#("fig5c.pdf", family="Arial Narrow", width=6, height=4)
plot(gcom2$Diff, gcom2$meanCEcos,
     xlim=range(c(1e-7, 1e-5)),
     ylim=range(c(gcom1$minCEvol, gcom1$maxCEvol)),
     type='l', lwd = 5, col=alpha(rgb(0,0,1), 0.5), frame = FALSE, pch=19,
     xlab=expression(paste("Diffusion Rate ", (cm^{2} ~ h^{-1}))), 
     ylab=expression(paste("Soil C Stock ", (mg))))
lines(gcom1$Diff, gcom1$meanCEvol, type='p', pch = 19, col = "red")
# hack: we draw arrows but with very special "arrowheads"
arrows(gcom1$Diff, gcom1$minCEvol, gcom1$Diff, gcom1$maxCEvol, 
       length=0.05, angle=90, code=3)
#dev.off()



