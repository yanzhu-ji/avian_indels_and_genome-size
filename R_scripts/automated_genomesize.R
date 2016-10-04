# This script is aimed to automatically calculate genome size based on ksa-estimated coverage

# 0. estimate k-mer coverage using Kmerspectrumanalyzer; record estimated coverage in file "ksa-coverage_47-species.list".
# 1. read file "ksa-coverage_47-species.list"
# 2. read species_m31_U5.histo files
# 3. estimate where should be cutted off based on coverage (usually 1/4-1/3 of coverage, so define a region and get the lowest value, then throw lower values off)

#   NOTE: this was an earlier measurement for k-mer genome size that was eventually abandoned due to various reasons 
#         (for example, we did not evaluate whether the modified genome size estimation accurately 
#         reflects k-mer frequencies of eukaryotes).  However, the results can be very useful in 
#         terms of measuring the simpler version of k-mer genome size efficiently.  Of course this 
#         could be resolved with proper R codes, but we didn't investigate additional time on it.
#         Instead we just used the results of KSA to facilitate estimating simple k-mer genome sizes
#         automatically.

#   Yanzhu Ji, yanzhuji20@gmail.com, last modified 10.1.2016

library(phytools)
setwd("/Users/yanzhuji/R-scripts/R-chp3/kmer-c-value/")

species.list=c("acaChl", "anaPla", "apaVit", "aptFor", "balReg",
               "bucRhi", "calAnn", "capCar", "carCri", "catAur",
               "chaPel", "chaVoc", "chlMac", "colLiv", "colStr",
               "corBra", "cucCan", "egrGar", "eurHel", "falPer",
               "fulGla", "galGal", "gavSte", "geoFor", "halAlb",
               "halLeu", "lepDis", "manVit", "melGal", "melUnd",
               "merNub", "mesUni", "nesNot", "nipNip", "opiHoa",
               "pelCri", "phaCar", "phaLep", "phoRub", "picPub",
               "podCri", "pteGut", "pygAde", "strCam", "tauEry",
               "tinGut", "tytAlb"
)

length(species.list) # 47, taeGut is excluded because insufficient amount of short reads


### kmer genomesize with coverage estimated by modified version of Kmerspectrumanalyzer (KSA)
ksa.coverage=read.table("ksa-coverage_47-species.list", row.names = 1)
colnames(ksa.coverage) = "ksa.coverage"

ksa.gs = data.frame( error=numeric(0), genomesize=numeric(0) )

for (i in species.list){
  infile = paste("histo/", i, "_m31_U5.histo", sep="")
  #  print (infile)
  i.histo = read.table(file=infile) # , row.names = 1
  colnames(i.histo) = c("cov", "freq")
  #  plot (i.histo[5:100,])
  i.ksapeak = ksa.coverage[i,]
  i.error = i.histo[i.histo$freq == min(i.histo[i.histo$cov <= 0.35 * i.ksapeak, ]$freq), ]$cov
  i.genomesize = sum(as.numeric(i.histo[i.error:nrow(i.histo), 1] * i.histo[i.error : nrow(i.histo), 2]))/i.ksapeak
  
  i.estimates = data.frame(i.genomesize, i.error, row.names=i)
  
  ksa.gs = rbind(ksa.gs, i.estimates)
}

ksa.gs$cov = ksa.coverage[row.names(ksa.gs),]

colnames(ksa.gs)=c("ksa.gs", "error.cov", "ksa.cov")

### simple way of kmer coverage

simple.gs = data.frame( error=numeric(0), simple.gs=numeric(0), simple.cov=numeric(0) )

for (i in species.list){
  infile = paste("histo/", i, "_m31_U5.histo", sep="")
#  infile = ("histo/tinGut_m31_U5.histo") #test
  #  print (infile)
  i.histo = read.table(file=infile) # , row.names = 1
  colnames(i.histo) = c("cov", "freq")
  #  plot (i.histo[5:100,])
  ksa.cov = ksa.coverage[i,]
  i.mainpeak = i.histo[i.histo$freq == max(i.histo[i.histo$cov > 0.5 * ksa.cov , ]$freq), ]$cov # used ksa coverage to help delimit the location of the mainpeak
  i.error = i.histo[i.histo$freq == min(i.histo[i.histo$cov <= 0.35 * ksa.cov, ]$freq), ]$cov
  i.genomesize = sum(as.numeric(i.histo[i.error:nrow(i.histo), 1] * i.histo[i.error : nrow(i.histo), 2]))/i.mainpeak
  
  i.estimates = data.frame(i.genomesize, i.error, i.mainpeak, row.names=i)
  
  simple.gs = rbind(simple.gs, i.estimates)
}

colnames(simple.gs) = c("simple.gs", "error.cov", "simple.cov")

### plot the k-mer spectrum to check if every min and max has been correctly identified

pdf(file='k-mer_spectrum.pdf', width=8, height=40)
par(mfrow=c(18, 3),mar=c(2,2,2,2))
for (i in species.list){
  #test i = "tinGut"
  infile = paste("histo/", i, "_m31_U5.histo", sep="")
  #  print (infile)
  i.histo = read.table(file=infile) # , row.names = 1
  if ( i == "catAur"){ # except for "catAur", all other plots have lowest x-axis value of (error_peak - 2) but catAur has (error_peak - 1) to avoid flatterned k-mer spectrum..
    plot (i.histo[(simple.gs[i, "error.cov"]-1):150,], main=i, xlab="k-mer frequency", ylab="number of k-mers",
          pch=20, cex=0.8, xaxp=c(0, 150, 3), xlim=c(0, 150))
    points(simple.gs[i, "error.cov"], i.histo[(simple.gs[i, "error.cov"]),]$V2, col="red")
    points(simple.gs[i, "simple.cov"], i.histo[(simple.gs[i, "simple.cov"]),]$V2, col="green")
  }else{
    plot (i.histo[(simple.gs[i, "error.cov"]-2):150,], main=i, xlab="k-mer frequency", ylab="number of k-mers",
          pch=20, cex=0.8, xaxp=c(0, 150, 3), xlim=c(0, 150))
    points(simple.gs[i, "error.cov"], i.histo[(simple.gs[i, "error.cov"]),]$V2, col="red")
    points(simple.gs[i, "simple.cov"], i.histo[(simple.gs[i, "simple.cov"]),]$V2, col="green")
  }
}
dev.off()

### three errors:
i="catAur"
simple.gs["catAur",2]=4 # used to be 3
infile = paste("histo/", i, "_m31_U5.histo", sep="")
i.histo = read.table(file=infile)
simple.gs["catAur", "simple.gs"] = 
  sum(as.numeric(i.histo[4:nrow(i.histo), 1]* i.histo[4 : nrow(i.histo), 2]))/9

i="tinGut"
simple.gs["tinGut",3]=71 # used to be 38
infile = paste("histo/", i, "_m31_U5.histo", sep="")
i.histo = read.table(file=infile)
simple.gs["tinGut", "simple.gs"]=
  sum(as.numeric(i.histo[15:nrow(i.histo), 1]* i.histo[15 : nrow(i.histo), 2]))/71
 # 1096625944

i="colLiv"
simple.gs["colLiv", "error.cov"]=34
infile = paste("histo/", i, "_m31_U5.histo", sep="")
i.histo = read.table(file=infile)
simple.gs["colLiv", "simple.gs"]=
  sum(as.numeric(i.histo[34:nrow(i.histo), 1] * i.histo[34 : nrow(i.histo), 2]))/71

### re-run the plotting part to double check.

### save the RData
save(simple.gs, file="../simple.gs.RData")

##############################################################
