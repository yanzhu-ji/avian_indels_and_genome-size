# this script is used to measure scaled length of CR1 (aka consensus-scaled mean ratio, 
#   abbreviated as cmean.ratio or cmr)

setwd("/Users/yanzhuji/R-scripts/R-chp3/lengths/")
length.matrix=data.frame(row.names=row.names(ave.48))

my.files = list.files(pattern="*_trimmed_length.txt")

for (file in my.files){
  # file="KL217251.1_864634_866046_trimmed_length.txt"
  locus.info=read.table(file, row.names=1)
  length.matrix[,file]=locus.info[row.names(length.matrix),] 
}

write.table(length.matrix, file="CR1-lengths.csv", sep=",") 

### scale lengths of loci by consensus lengths
# ratio compared to consensus
consensus.length = read.table("consensus-length.txt", row.names=1)

cratio.matrix=matrix(nrow=48, ncol=50)

for ( n in 1:ncol(length.matrix)){
  cratio.matrix[,n] = length.matrix[,n]/consensus.length[colnames(length.matrix)[n],"V2"]
}

write.table(cratio.matrix, file="scaled.length.csv", sep=",")

cmean.ratio=vector()

for ( m in 1:nrow(length.matrix)){
  cmean.ratio[m] = mean(cratio.matrix[m,], na.rm=T)
}

names(cmean.ratio)=row.names(length.matrix)
save(cmean.ratio, file="../cmean.ratio.RData")
