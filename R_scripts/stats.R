#   This script includes the statistical analyses in the paper "Relationships among 
#       powered flight, metabolic rate, body mass, genome size, and the retrotransposon 
#       complement of volant birds" (Y. Ji and J. A. DeWoody)

### Summary
### 1. read in data
### 2. stats 
###    2.1 genome size measurements
###    2.2 linear models (a) - (f), including measurements of phylogenetic signals.
###    2.3 Mann-Whitney tests (g) - (j), including power analyses using simulation.
###    2.4 Levene tests and associated power analyses using simulation.
###    2.5 D(OG) 
### 3. plotting
###    3.1 correlations (a) - (f)
###    3.2 box plots (g) - (j)
###    3.3 box plot D(OG)

###    Yanzhu Ji (yanzhuji20@gmail.com), last modified 10/1/2016

########### script starts here ###########
###  1. read in ave.48, and add the four variables:
###     simple genome size ("simple.gs"), body mass ("bm"), "mean ratio", and FY/GD ("dummy")
library(phytools)
setwd("/Users/yanzhuji/R-scripts/R-chp3/") # folder on Mac load("stats.RData") # for the convenience of rerunning data; will delete in scripts for uploading.

##   read the phylogeny, and root it.
ave.tree=read.tree("TENT+outgroup_abb.ExaML.tre")# load("ave.rooted.RData") # the same...
fastMRCA(ave.tree, "Human", "Green_sea_turtle") # node 58 is the outgroup.
ave.rooted = root(ave.tree, node=58, resolve.root=TRUE) 
plot(ave.rooted, cex=0.7) # double check

##   load csv file with species names, abbreviations, sequencing assembly data 
##     from Zhang et al. 2014 Science paper, and AGSD genome size (agsd.gs)
ave.48=read.csv("assembly-size-48.csv", sep=",", row.names=3) # spotted there is a typo in "opiHoa"

##   run "automated_genomesize.R" to get the "simple.gs.RData".
load("simple.gs.RData")
ave.48$simple.gs=simple.gs[row.names(ave.48),]$simple.gs/1e9
ave.48$log.sgs=log(ave.48$simple.gs) # "sgs": simple.gs

##   body mass (bm)
bm=read.csv("CRC-selected-species.csv", sep=",", row.names=4)
ave.48$bm=bm[row.names(ave.48),]$mean
ave.48$log.bm=log(ave.48$bm)
length(na.omit(ave.48$bm)) # 46; 2 species do not have body mass records

##   basal metabolic rates (bmr)
bmr=read.csv("metabolic_Makarieva-2008-PNAS.csv", header=TRUE, row.names=3) # added to RData 11/9/2015
ave.48$bmr=bmr[row.names(ave.48),]$qWkg
ave.48$log.bmr=log(ave.48$bmr)

##  scaled length of CR1
load("cmean.ratio.RData") # cmean.ratio: consensus-scaled mean ratio (CR1 divided by consensus)
ave.48$cmean.ratio=cmean.ratio[row.names(ave.48)]

gd.list=c("strCam", "pygAde", "aptFor", "galGal", "melGal", 
          "tinGut", "carCri") 
fl.list=c("acaChl", "anaPla", "apaVit",           "balReg",
          "bucRhi", "calAnn", "capCar",           "catAur",
          "chaPel", "chaVoc", "chlMac", "colLiv", "colStr",
          "corBra", "cucCan", "egrGar", "eurHel", "falPer",
          "fulGla",           "gavSte", "geoFor", "halAlb",
          "halLeu", "lepDis", "manVit",           "melUnd",
          "merNub", "mesUni", "nesNot", "nipNip", "opiHoa",
          "pelCri", "phaCar", "phaLep", "phoRub", "picPub",
          "podCri", "pteGut",           "taeGut", "tauEry",
          "tytAlb"
)
length(gd.list) # 7
length(fl.list) # 41

View(ave.48)

###  2. stats
###  2.1 Genome size estimations
cor.test(ave.48$simple.gs, ave.48$agsd.gs, alternative = "g") 
# r2 = 0.72^2 = 0.52, p = 0.001426, df = 12

cor.test(ave.48$assembly.size, ave.48$agsd.gs, alternative = "g")  
# r2 = 0.48^2 = 0.23, p = 0.03627, df = 13

cor.test(ave.48$assembly.size, ave.48$simple.gs, alternative = "g") 
# r = 0.43^2 = 0.19, p = 0.001229, df = 45

## 2.2 use simple regression models to test pairwise correlations for hypotheses (a) - (f) 
# (a) body mass vs. bmr
# phylogenetic signal of residuals
bmr.bm.sp=row.names(ave.48[ (!is.na(ave.48$bmr)==TRUE) & (!is.na(ave.48$log.bm)==TRUE), ]) # 13 species
bmr.bm.pruned=drop.tip(ave.rooted, setdiff(ave.rooted$tip, bmr.bm.sp))
bmr.bm.df=ave.48[bmr.bm.pruned$tip,]
nrow(bmr.bm.df)
phylosig(bmr.bm.pruned, as.numeric(residuals(lm(bmr.bm.df$log.bmr~bmr.bm.df$log.bm))), test=T, method="lambda") 
phylosig(bmr.bm.pruned, as.numeric(residuals(lm(bmr.bm.df$log.bmr~bmr.bm.df$log.bm))), test=T) 

summary(lm(ave.48$log.bmr~ave.48$log.bm))
# adjusted R2 = 0.7298, p = 0.000123/2 < 0.0001

#d pwr.r.test(n=13, power=0.8, sig.level = 0.05)
# r = 0.69, passed!


# (b) genome size vs. TE length
gs.cmr.sp=row.names(ave.48[ (!is.na(ave.48$simple.gs)==TRUE) & (!is.na(ave.48$cmean.ratio)==TRUE), ]) 
gs.cmr.pruned=drop.tip(ave.rooted, setdiff(ave.rooted$tip, gs.cmr.sp))
gs.cmr.df=ave.48[gs.cmr.pruned$tip,]
nrow(gs.cmr.df) # 46
phylosig(gs.cmr.pruned, as.numeric(residuals(lm(gs.cmr.df$log.sgs~gs.cmr.df$cmean.ratio))), test=T, method="lambda") 
phylosig(gs.cmr.pruned, as.numeric(residuals(lm(gs.cmr.df$log.sgs~gs.cmr.df$cmean.ratio))), test=T) 

summary(lm(ave.48$log.sgs~ave.48$cmean.ratio))
# adjusted r2 = 0.05, p = 0.0659/2 = 0.03

#d pwr.r.test(n=46, sig.level=0.05, power=0.8)
# r = 0.40, not passed...

# (c) genome size vs bmr
gs.bmr.sp=row.names(ave.48[ (!is.na(ave.48$simple.gs)==TRUE) & (!is.na(ave.48$bmr)==TRUE), ]) # 13 species
gs.bmr.pruned=drop.tip(ave.rooted, setdiff(ave.rooted$tip, gs.bmr.sp))
gs.bmr.df=ave.48[gs.bmr.pruned$tip,]
nrow(gs.bmr.df)
phylosig(gs.bmr.pruned, as.numeric(residuals(lm(gs.bmr.df$log.sgs~gs.bmr.df$log.bmr))), test=T, method="lambda") 
phylosig(gs.bmr.pruned, as.numeric(residuals(lm(gs.bmr.df$log.sgs~gs.bmr.df$log.bmr))), test=T) 
# good without correction. 

summary(lm(ave.48$log.sgs~ave.48$log.bmr))
# adjusted r2 = 0.10, p = 0.1498/2 = 0.07
#d pwr.r.test(n=13, sig.level=0.05, power=0.8)
# r = 0.69, not passed...

# (d) cmean.ratio vs. BMR
cmr.bmr.sp=row.names(ave.48[ (!is.na(ave.48$cmean.ratio)==TRUE) & (!is.na(ave.48$bmr)==TRUE), ]) # 13 species
cmr.bmr.pruned=drop.tip(ave.rooted, setdiff(ave.rooted$tip, cmr.bmr.sp))
cmr.bmr.df=ave.48[cmr.bmr.pruned$tip,]
nrow(cmr.bmr.df)
phylosig(cmr.bmr.pruned, as.numeric(residuals(lm(cmr.bmr.df$cmean.ratio~cmr.bmr.df$log.bmr))), test=T, method="lambda") 
phylosig(cmr.bmr.pruned, as.numeric(residuals(lm(cmr.bmr.df$cmean.ratio~cmr.bmr.df$log.bmr))), test=T) 
# good without correction.

summary(lm(ave.48$cmean.ratio~ave.48$log.bmr))
# adjusted R-squared: 0.4676, p-value = 0.00596/2 = 0.003,

#d pwr.r.test(n=13, sig.level=0.05, power=0.8)
# power test: r = 0.69, nearly passed...(r = 0.68)

# (e) genome size vs body mass
gs.bm.sp=row.names(ave.48[ (!is.na(ave.48$simple.gs)==TRUE) & (!is.na(ave.48$log.bm)==TRUE), ]) # 45 species
gs.bm.pruned=drop.tip(ave.rooted, setdiff(ave.rooted$tip, gs.bm.sp))
gs.bm.df=ave.48[gs.bm.pruned$tip,]
nrow(gs.bm.df)
phylosig(gs.bm.pruned, as.numeric(residuals(lm(gs.bm.df$log.sgs~gs.bm.df$log.bm))), test=T, method="lambda") 
phylosig(gs.bm.pruned, as.numeric(residuals(lm(gs.bm.df$log.sgs~gs.bm.df$log.bm))), test=T) 
# good without correction.

summary(lm(ave.48$log.sgs~ave.48$log.bm))
# adjusted r2 = 0.18, p = 0.002/2 = 0.001

# pwr.r.test(n=45, power=0.8, sig.level=0.05)
# r = 0.40 (r2 = 0.16), passed!


# (f) cmean.ratio vs. body mass
cmr.bm.sp=row.names(ave.48[ (!is.na(ave.48$cmean.ratio)==TRUE) & (!is.na(ave.48$log.bm)==TRUE), ]) # 45 species
cmr.bm.pruned=drop.tip(ave.rooted, setdiff(ave.rooted$tip, cmr.bm.sp))
cmr.bm.df=ave.48[cmr.bm.pruned$tip,]
nrow(cmr.bm.df)
phylosig(cmr.bm.pruned, as.numeric(residuals(lm(cmr.bm.df$cmean.ratio~cmr.bm.df$log.bm))), test=T, method="lambda") 
phylosig(cmr.bm.pruned, as.numeric(residuals(lm(cmr.bm.df$cmean.ratio~cmr.bm.df$log.bm))), test=T)
# good without correction.

summary(lm(ave.48$cmean.ratio~ave.48$log.bm))
# adjusted r2 = 0.4493, p << 0.001

#d pwr.r.test(n=45, sig.level=0.05, power=0.8)
#d r = 0.4493^0.5 = 0.67, r(min) = 0.40, passed!

### 2.3 (g) - (j), use Mann-Whitney tests to detect differences between FY and GD birds
##      (g) FY/GD vs bmr
wilcox.test(ave.48[fl.list,]$log.bmr, ave.48[gd.list,]$log.bmr, alternative="g")
# W=30, n1 = 3, n2 = 10, p = 0.003

## power by simulation
bmr.mw.p.all = vector(length=10000)
for ( i in c(1:10000) ) {
  fl.bmr.sample = rnorm(length(na.omit(ave.48[fl.list,]$log.bmr)), mean=mean(ave.48[fl.list,]$log.bmr, na.rm=T), sd=sd(ave.48[fl.list,]$log.bmr, na.rm=T))
  gd.bmr.sample = rnorm(length(na.omit(ave.48[gd.list,]$log.bmr)), mean=mean(ave.48[gd.list,]$log.bmr, na.rm=T), sd=sd(ave.48[gd.list,]$log.bmr, na.rm=T))
  group = as.factor(c(rep(1, length(fl.bmr.sample)), rep(2, length(gd.bmr.sample))))
  
  bmr.mw.p.all[i] = wilcox.test(fl.bmr.sample, gd.bmr.sample, alternative="g")$p.value
}
mean(bmr.mw.p.all<=0.05)
# 0.81

## (h) FY/GD vs body mass: n=45, n=13 respectively
wilcox.test(ave.48[fl.list,]$log.bm, ave.48[gd.list,]$log.bm, alternative="l")
length(na.omit(ave.48[fl.list,]$log.bm)) # 39
length(na.omit(ave.48[gd.list,]$log.bm)) # 7
# w = 40, p = 0.00097

## power by simulation
bm.mw.p.all = vector(length=10000)
for ( i in c(1:10000) ) {
  fl.bm.sample = rnorm(length(na.omit(ave.48[fl.list,]$log.bm)), mean=mean(ave.48[fl.list,]$log.bm, na.rm=T), sd=sd(ave.48[fl.list,]$log.bm, na.rm=T))
  gd.bm.sample = rnorm(length(na.omit(ave.48[gd.list,]$log.bm)), mean=mean(ave.48[gd.list,]$log.bm, na.rm=T), sd=sd(ave.48[gd.list,]$log.bm, na.rm=T))
  group = as.factor(c(rep(1, length(fl.bm.sample)), rep(2, length(gd.bm.sample))))
  
  bm.mw.p.all[i] = wilcox.test(fl.bm.sample, gd.bm.sample, alternative="l")$p.value
}
mean(bm.mw.p.all<=0.05)
# power = 0.96

# (i). genome size vs. FY/GD
wilcox.test(ave.48[fl.list,]$log.sgs, ave.48[gd.list,]$log.sgs, alternative="l")
#  W = 104, p = 0.15
length(na.omit(ave.48[fl.list,]$log.sgs)) # 40
length(na.omit(ave.48[gd.list,]$log.sgs)) # 7

# power by simulation
sgs.mw.p.all = vector(length=10000)
for ( i in c(1:10000) ) {
  fl.sgs.sample = rnorm(length(na.omit(ave.48[fl.list,]$log.sgs)), mean=mean(ave.48[fl.list,]$log.sgs, na.rm=T), sd=sd(ave.48[fl.list,]$log.sgs, na.rm=T))
  gd.sgs.sample = rnorm(length(na.omit(ave.48[gd.list,]$log.sgs)), mean=mean(ave.48[gd.list,]$log.sgs, na.rm=T), sd=sd(ave.48[gd.list,]$log.sgs, na.rm=T))
  group = as.factor(c(rep(1, length(fl.sgs.sample)), rep(2, length(gd.sgs.sample))))
  
  sgs.mw.p.all[i] = wilcox.test(fl.sgs.sample, gd.sgs.sample, alternative="l")$p.value
}

mean(sgs.mw.p.all<=0.05)
# power: 0.18...

# (j). cmean.ratio vs. FY/GD
wilcox.test(ave.48[fl.list,]$cmean.ratio, ave.48[gd.list,]$cmean.ratio, alternative="l")
# W = 45, p = 0.0053
length(na.omit(ave.48[fl.list,]$cmean.ratio)) # 41
length(na.omit(ave.48[gd.list,]$cmean.ratio)) # 6

# power by simulation
cmr.mw.p.all = vector(length=10000)
for ( i in c(1:10000) ) {
  fl.cmr.sample = rnorm(length(na.omit(ave.48[fl.list,]$cmean.ratio)), mean=mean(ave.48[fl.list,]$cmean.ratio, na.rm=T), sd=sd(ave.48[fl.list,]$cmean.ratio, na.rm=T))
  gd.cmr.sample = rnorm(length(na.omit(ave.48[gd.list,]$cmean.ratio)), mean=mean(ave.48[gd.list,]$cmean.ratio, na.rm=T), sd=sd(ave.48[gd.list,]$cmean.ratio, na.rm=T))
  group = as.factor(c(rep(1, length(fl.cmr.sample)), rep(2, length(gd.cmr.sample))))
  
  cmr.mw.p.all[i] = wilcox.test(fl.cmr.sample, gd.cmr.sample, alternative="l")$p.value
}

mean(cmr.mw.p.all<=0.05)
# power: 0.89

### 2.4 Levene's tests and assotiated power analysis using simulation
##  first test if the variances are the same using current data
ave.48[gd.list,"dummy"]=1
ave.48[fl.list,"dummy"]=0

##  CR1 length
leveneTest(ave.48$cmean.ratio, as.factor(ave.48$dummy))
# p = 0.15

## genome size
leveneTest(ave.48$log.sgs, as.factor(ave.48$dummy))
# p = 0.40

## next, use simulation to calculate power for Levene's test 
## genome size
sgs.l.p.all = vector(length=10000)
for ( i in c(1:10000) ) {
  fl.sgs.sample = rnorm(length(na.omit(ave.48[fl.list,]$log.sgs)), mean=mean(ave.48[fl.list,]$log.sgs, na.rm=T), sd=sd(ave.48[fl.list,]$log.sgs, na.rm=T))
  gd.sgs.sample = rnorm(length(na.omit(ave.48[gd.list,]$log.sgs)), mean=mean(ave.48[gd.list,]$log.sgs, na.rm=T), sd=sd(ave.48[gd.list,]$log.sgs, na.rm=T))
  group = as.factor(c(rep(1, length(fl.sgs.sample)), rep(2, length(gd.sgs.sample))))
  sgs.l.p.all[i] = leveneTest(c(fl.sgs.sample, gd.sgs.sample), group)$Pr[1]
}

mean(sgs.l.p.all)
# power = 0.46

## CR1 length
cmr.l.p.all = vector(length=10000)
for ( i in c(1:10000) ) {
  fl.cmr.sample = rnorm(length(na.omit(ave.48[fl.list,]$cmean.ratio)), mean=mean(ave.48[fl.list,]$cmean.ratio, na.rm=T), sd=sd(ave.48[fl.list,]$cmean.ratio, na.rm=T))
  gd.cmr.sample = rnorm(length(na.omit(ave.48[gd.list,]$cmean.ratio)), mean=mean(ave.48[gd.list,]$cmean.ratio, na.rm=T), sd=sd(ave.48[gd.list,]$cmean.ratio, na.rm=T))
  group = as.factor(c(rep(1, length(fl.cmr.sample)), rep(2, length(gd.cmr.sample))))
  cmr.l.p.all[i] = leveneTest(c(fl.cmr.sample, gd.cmr.sample), group)$Pr[1]
}
mean(cmr.l.p.all)
# power = 0.25

### 2.5 D(OG) comparison between species pairs within FY and GD using Mann-Whitney test
setwd("lengths")
my.files = list.files(pattern="*_trimmed_length.txt")
total.gd=data.frame(#"sp1"=character(), "sp2"=character(), 
  "gd1"=numeric(), "gd2"=numeric(), 
  "gd1.gs"=numeric(), "gd2.gs"=numeric(), stringsAsFactors=FALSE)

pair=rbind(combn(gd.list,2))

for ( i in 1:ncol(pair)) {
  gd1=pair[1,i]
  gd2=pair[2,i]
  
  pairwise=data.frame("gd1"=numeric(), "gd2"=numeric(), stringsAsFactors=FALSE) # 
  for (file in my.files){
    locus.info=read.table(file, row.names=1)
    if (!is.na(locus.info[gd1,]) & !is.na(locus.info[gd2,]) ){ #
      pairwise[file,]=c(locus.info[gd1,], locus.info[gd2,] )
    }  
  }
  # save combined lengths to total.gd
  total.gd[paste(gd1, gd2, sep=","),]=c(#gd1, gd2, 
    sum( pairwise$gd1 ),
    sum( pairwise$gd2 ),
    ave.48[gd1,]$simple.gs,
    ave.48[gd2,]$simple.gs)
  # stats on correlation between genome size and combined lengths in this two species...PIC needed? perhaps not..
  gd1.gs = ave.48[gd1,]$simple.gs
  gd2.gs = ave.48[gd2,]$simple.gs
}

# now total.fl
# sample.fl.list=sample(fl.list, length(gd.list) )

total.fl=data.frame(#"sp1"=character(), "sp2"=character(), 
  "fl1"=numeric(), "fl2"=numeric(), 
  "fl1.gs"=numeric(), "fl2.gs"=numeric(), stringsAsFactors=FALSE)

# or,
pair=rbind(combn(fl.list,2))

for ( i in 1:ncol(pair)) {
  fl1=pair[1,i]
  fl2=pair[2,i]
  
  pairwise=data.frame("fl1"=numeric(), "fl2"=numeric(), stringsAsFactors=FALSE) # 
  for (file in my.files){
    locus.info=read.table(file, row.names=1)
    if (!is.na(locus.info[fl1,]) & !is.na(locus.info[fl2,]) ){ #
      pairwise[file,]=c(locus.info[fl1,], locus.info[fl2,] )
    }  
  }
  # save combined lengths to total.fl
  total.fl[paste(fl1, fl2, sep=","),]=c(#fl1, fl2, 
    sum( pairwise$fl1 ),
    sum( pairwise$fl2 ),
    ave.48[fl1,]$simple.gs,
    ave.48[fl2,]$simple.gs)
  # stats on correlation between genome size and combined lengths in this two species...PIC needed? perhaps not..
  fl1.gs = ave.48[fl1,]$simple.gs
  fl2.gs = ave.48[fl2,]$simple.gs
}

wilcox.test(abs(total.fl$fl1-total.fl$fl2)/((total.fl$fl1+total.fl$fl2)/2),
            abs(total.gd$gd1-total.gd$gd2)/((total.gd$gd1+total.gd$gd2)/2),
            alternative="g")

save.image("stats.RData")
### 3. plotting
##  3.1 correlations (a) - (f)
# (a) body mass vs. bmr
plot(ave.48$log.bmr~ave.48$log.bm, xlab="log (body mass)", ylab="log (BMR)", cex.lab=1.5,
     pch=21, bg="gold", cex=1.2)
points(ave.48[gd.list,]$log.bm, ave.48[gd.list,]$log.bmr, pch=19, cex=1.2, col="dodgerblue3")
abline(lm(ave.48$log.bmr~ave.48$log.bm), col="red", lwd=2)
# optional
text(ave.48$log.bm, ave.48$log.bmr, labels=row.names(ave.48))

# (b) genome size vs. TE length
plot(ave.48$log.sgs~ave.48$cmean.ratio, xlab="scaled length of CR1", ylab="log (genome size)", cex.lab=1.5,
     pch=21, bg="gold", cex=1.2)
points(ave.48[gd.list,]$cmean.ratio, ave.48[gd.list,]$log.sgs, pch=19, cex=1.2, col="dodgerblue3" )
abline(lm(ave.48$log.sgs~ave.48$cmean.ratio), col="red", lwd=2)

# (c) genome size vs bmr
plot(ave.48$log.sgs~ave.48$log.bmr, xlab="log (BMR)", ylab="log (genome size)", cex.lab=1.5,
     pch=21, bg="gold", cex=1.2)
points(ave.48[gd.list,]$log.bmr, ave.48[gd.list,]$log.sgs, pch=19, cex=1.2, col="dodgerblue3"  )
abline(lm(ave.48$log.sgs~ave.48$log.bmr), col="red", lwd=2)

# (d) cmean.ratio vs. BMR
plot(ave.48$cmean.ratio~ave.48$log.bmr, xlab="log (BMR)", ylab="scaled length of CR1", cex.lab=1.5,
     pch=21, bg="gold", cex=1.2)
points(ave.48[gd.list,]$log.bmr, ave.48[gd.list,]$cmean.ratio, pch=19, cex=1.2, col="dodgerblue3" )
abline(lm(ave.48$cmean.ratio~ave.48$log.bmr), col="red", lwd=2)

# (e) genome size vs body mass
plot(ave.48$log.sgs~ave.48$log.bm, xlab="log (body mass)", ylab="log (genome size)", cex.lab=1.5,
     pch=21, bg="gold", cex=1.2)
points(ave.48[gd.list,]$log.bm, ave.48[gd.list,]$log.sgs, pch=19, cex=1.2, col="dodgerblue3" )
abline(lm(ave.48$log.sgs~ave.48$log.bm), col="red", lwd=2)

# (f) cmean.ratio vs. body mass
plot(ave.48$cmean.ratio~ave.48$log.bm, xlab="log (body mass)", ylab="scaled length of CR1", cex.lab=1.5,
     pch=21, bg="gold", cex=1.2)
points(ave.48[gd.list,]$log.bm, ave.48[gd.list,]$cmean.ratio, pch=19, cex=1.2, col="dodgerblue3" )
abline(lm(ave.48$cmean.ratio~ave.48$log.bm), col="red", lwd=2)

## 3.2 plot Mann-Whitney tests (g) - (j), using box plots
ave.48[gd.list,"dummy"]=1
ave.48[fl.list,"dummy"]=0

# (g) FY/GD vs bmr
boxplot(ave.48$log.bmr~ave.48$dummy,
        names=c("FY (n=10)", "GD (n=3)"), boxwex=0.4, ylab="log (BMR)", col=c("gold", "dodgerblue3"), cex.lab=1.2)

nrow(ave.48[!is.na(ave.48$bmr)==TRUE,])
row.names(ave.48[!is.na(ave.48$bmr)==TRUE,][gd.list,]) 

# (h) FY/GD vs body mass: n=45, n=13 respectively
boxplot(ave.48$log.bm~ave.48$dummy,
        names=c("FY (n=39)", "GD (n=7)"), boxwex=0.4, ylab="log (body mass)", col=c("gold", "dodgerblue3"), cex.lab=1.2)

nrow(ave.48[!is.na(ave.48$bm)==TRUE,]) # total = 46
row.names(ave.48[!is.na(ave.48$bm)==TRUE,][gd.list,])  # no.gd = 7


# (i). genome size vs. FY/GD
boxplot(ave.48$log.sgs~ave.48$dummy,
        names=c("FY (n=40)", "GD (n=7)"), boxwex=0.4, ylab="log (genome size)", col=c("gold", "dodgerblue3"), cex.lab=1.2)

nrow(ave.48[!is.na(ave.48$log.sgs)==TRUE,]) # total = 47
row.names(ave.48[!is.na(ave.48$log.sgs)==TRUE,][gd.list,]) # no.gd = 7

# (j). cmean.ratio vs. FY/GD
boxplot(ave.48$cmean.ratio~ave.48$dummy,
        names=c("FY (n=41)", "GD (n=6)"), boxwex=0.4, ylab="scaled length of CR1", col=c("gold", "dodgerblue3"), cex.lab=1.2)

nrow(ave.48[!is.na(ave.48$cmean.ratio)==TRUE,]) # total = 47
row.names(ave.48[!is.na(ave.48$cmean.ratio)==TRUE,][gd.list,])  # no.gd = 6

##  3.3 plot variation comparison between gd-gd and fl-fl

boxplot(abs(total.fl$fl1-total.fl$fl2)/((total.fl$fl1+total.fl$fl2)/2),
        abs(total.gd$gd1-total.gd$gd2)/((total.gd$gd1+total.gd$gd2)/2), 
        names=c("FY","GD"), boxwex=0.4,
        ylab="pairwise difference of CR1", 
        col=c("gold", "dodgerblue3"))
