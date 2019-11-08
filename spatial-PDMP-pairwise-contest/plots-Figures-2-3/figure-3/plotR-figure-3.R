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
library(Cairo)



########################## FIGURE 3
############# COMMUNS
#load files "scoresCommuns_ID.txt"
nom.var = c("IDcluster","IDjob","timeMax","phi1","phi2","diffD","VmU","KmU","dM","dZ",
            "gM","gZ","p","VmD","IC","ID","lC","lD","fracInit2","fracFin2","timeIndex",
            "popTotFin","txAccr","txAccrCond")
#Following files contain most values of Diff presented in Figure 3 (a few are missing, e.g. 1e-5)
fcom1 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2762.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom2 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2763.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom3 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2766.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom4 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2767.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom5 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2770.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom6 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2775.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom7 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2781.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom8 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2783.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom9 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2845.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom10 = read.table('~/Desktop/spatial-PDMP-pairwise-contest/plots-Figures-2-3/figure-3/scoresCommuns_2863.txt',sep= "",dec='.',header=F, col.names=nom.var)
fcom = rbind(fcom1,fcom2,fcom3,fcom4,fcom5,fcom6,fcom7,fcom8,fcom9,fcom10)
#fcom = read.table('~/Dropbox/Helene_Elsa_ecoevo_model/programmes/tests_Elsa/all_Diff_sur_R_mars2019/scoresCommuns_2767.txt',sep= "",dec='.',header=F, col.names=nom.var)
#NA in simulations where microbes went extinct
fcom$txAccrCond[fcom$txAccrCond == -1000] = NA
#Multiply growth rate by timeMax
fcom$txAccrScale = fcom$timeMax*fcom$txAccrCond
#To check number of simulations per test
fcom = cbind(fcom, nbrRep=1)
#New table per test with mean, standard deviation, number of simulations
mtab = aggregate(list(meantxAccrScale = fcom$txAccrScale), 
                 by = list(phi1 = fcom$phi1, phi2 = fcom$phi2, diffD = fcom$diffD), mean, na.rm=T)
sdtab = aggregate(list(sdtxAccrScale = fcom$txAccrScale), 
                  by = list(phi1 = fcom$phi1, phi2 = fcom$phi2, diffD = fcom$diffD), sd, na.rm=T)
reptab = aggregate(list(sumnbrRep = fcom$nbrRep), 
                   by = list(phi1 = fcom$phi1, phi2 = fcom$phi2, diffD = fcom$diffD), sum, na.rm=T)
gcom = cbind(mtab, sdtxAccrScale = sdtab$sdtxAccrScale, sumnbrRep = reptab$sumnbrRep)
#Add column "Number of survivors"
fcom = cbind(fcom, surv=round(fcom$fracFin2,2))
fcom$surv[fcom$surv != 0.00] = 1
hcom = aggregate(list(nbrSurv = fcom$surv), 
                 by = list(phi1 = fcom$phi1, phi2 = fcom$phi2, diffD = fcom$diffD), sum, na.rm=T)
kcom = cbind(gcom, hcom$nbrSurv)
#Add column "estimation of survival probability" (Bernoulli)
kcom$psurv = kcom$`hcom$nbrSurv` / kcom$sumnbrRep
#Indicate "producer" (higher phi) or "cheater" (lower phi)
kcom$mutantType = ifelse(kcom$phi1 < kcom$phi2, "producer", "cheater")
#Remove negative fitness
kcom$meantxAccrScale2 = ifelse(kcom$meantxAccrScale < 0, 0, kcom$meantxAccrScale)
kcom$sdtxAccrScale2 = ifelse(kcom$meantxAccrScale < 0, 0, kcom$sdtxAccrScale)
#Add column growth rate x survival probability
kcom$txaccXpsurv = kcom$meantxAccrScale2*kcom$psurv
#Subset per value of Diff  (example here for Diff=5x10^-7)
kcomdiff = subset(kcom,kcom$diffD==5e-7 & kcom$phi1<0.3 & kcom$phi2<0.3)
#Single - producer
kcomprod = subset(kcomdiff,kcomdiff$mutantType=='producer')
kcomcheat = subset(kcomdiff,kcomdiff$mutantType=='cheater')
#PLOT
#cairo_pdf("fig4k.pdf", family="Arial Narrow", width=5, height=4)
barplot(kcomprod$meantxAccrScale2*kcomprod$psurv,kcomprod$phi2, 
        width=0.5, 
        col=rgb(0,205/255,205/255,1/2), border=NA,
        names.arg = NULL, 
        legend.text = NULL,
        space=0,
        ylim=c(0,3.5))
barplot(kcomcheat$meantxAccrScale2*kcomcheat$psurv,kcomcheat$phi2, axes=F,
        width=0.5, 
        col=rgb(1,193/255,193/255,1/2), border=NA, 
        names.arg = c("0.05-0.1","0.1-0.15","0.15-0.2","0.2-0.25"),
        legend.text = NULL,
        space=0,
        ylim=c(0,3.5),
        add=T,
        xlab=expression(paste("Exoenzyme allocation traits (\u03C6) of competing strains")),
        ylab="Fitness x TimeMax x Proba Survival") 
#dev.off()
