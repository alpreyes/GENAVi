
final.counts <- read.csv("final.counts.csv",header = TRUE) ##read in data going to view
count.matrix <- final.counts[,c(1,8:25)]

##filter rows with only zeroes or only one count (little to no info on gene expression)
dim(count.matrix)
count.matrix.filtered <- count.matrix[rowSums(count.matrix[,2:19]) > 1,-1] ##remove gene names column for now, to figure out norm techniques
dim(count.matrix.filtered)

cmf.rownorm <- data.frame()
#colnames(count.matrix.filtered.rnorm) <- colnames(count.matrix.filtered)

for(i in 1:dim(count.matrix.filtered)[1]) ##using cmf bc cant have sd of all-zero-row
{
  cmf.rownorm <- rbind(cmf.rownorm, (as.numeric(count.matrix.filtered[i,]) - mean(as.numeric(count.matrix.filtered[i,])))/sd(as.numeric(count.matrix.filtered[i,])))
}

colnames(cmf.rownorm) <- colnames(count.matrix.filtered)

##need some way to pull out gene names from filtered subset matrices
write.csv(cmf.rownorm,"cmf.rownorm.csv",row.names = FALSE)

##trying to make row norm matrix without gene filtering step

final.rownorm <- data.frame()

for(i in 1:dim(count.matrix)[1])
{
  final.rownorm <- rbind(final.rownorm, (as.numeric(count.matrix[i,-1]) - mean(as.numeric(count.matrix[i,-1]))) / sd(as.numeric(count.matrix[i,-1])))
}

##this command works but has NaN for rows with all zeroes bc divide by sd of zero

##add other info columns before writing final version of table
final.rownorm <- cbind(final.counts[,1:7],final.rownorm)

colnames(final.rownorm) <- colnames(final.counts)

##now write tables to include in app
write.csv(final.rownorm, "final.rownorm.csv",row.names = FALSE)

rm(final.rownorm)
rm(cmf.rownorm)

#final.rownorm <- read.csv("final.rownorm.csv",header = TRUE)
#cmf.rownorm <- read.csv("cmf.rownorm.csv",header = TRUE)

####### for publication/rmarkdown put all making data table code in this script #####



