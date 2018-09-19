
# # source("https://bioconductor.org/biocLite.R")
# # biocLite("rhdf5")
# library(rhdf5)
# 
# h5f <- H5Fopen("/Users/reyesa1/Desktop/human_matrix.h5")
# 
# h5f
# 
# #h5f_test <- head(h5f$data)
# 
# h5f$data
# 
# h5read(h5f$data)

library(readr)

test_large_counts <- read_csv("/Users/reyesa1/Desktop/test_large_counts/CountTable_withDups.csv") ##clarify with becem: how this featurecounts matrix was generated

samples_test_large_counts <- read_csv("/Users/reyesa1/Desktop/test_large_counts/SamplesTable.csv") ##see that order of rows in samples table is same as order of cols in count matrix
##--> use same index vectors -1 to subset rows from sample table

## some samples/columns have NAs for group info so can't do DEA, only sample from cols with all info
cols_samples_test_large_counts <- samples_test_large_counts[,c(2,3,22,23)] ##these cols used as metadata fro DEA
dim(cols_samples_test_large_counts)

which(is.na(cols_samples_test_large_counts), arr.ind = TRUE)[,1] ##exclude these cols from full data table (+1 for gene id col)
test_large_counts <- test_large_counts[,-(which(is.na(cols_samples_test_large_counts), arr.ind = TRUE)[,1]+1)] ##### use this as full data table --> write to a new file
dim(test_large_counts) ## reduced from 985 cols to 979 cols, 6 excluded had NA in sample info

cols_samples_test_large_counts <- na.omit(cols_samples_test_large_counts) ### use this as full samples metadata --> write to a new file
dim(cols_samples_test_large_counts)



## get column indeces to pull x random columns from large data set to run through GENAVi (time and see when it crashes)
set.seed(1234)
col_index_25  <- sample(2:979, size = 25, replace = FALSE)
col_index_30  <- sample(2:979, size = 30, replace = FALSE)
col_index_50  <- sample(2:979, size = 50, replace = FALSE)
col_index_100 <- sample(2:979, size = 100, replace = FALSE)
col_index_105 <- sample(2:979, size = 105, replace = FALSE)
col_index_110 <- sample(2:979, size = 110, replace = FALSE)
col_index_125 <- sample(2:979, size = 125, replace = FALSE)
col_index_150 <- sample(2:979, size = 150, replace = FALSE)
col_index_175 <- sample(2:979, size = 175, replace = FALSE)
col_index_190 <- sample(2:979, size = 190, replace = FALSE)
col_index_200 <- sample(2:979, size = 200, replace = FALSE)
col_index_500 <- sample(2:979, size = 500, replace = FALSE)
col_index_800 <- sample(2:979, size = 800, replace = FALSE)

test_large_counts_25  <- test_large_counts[,append(c(1),col_index_25)]
test_large_counts_30  <- test_large_counts[,append(c(1),col_index_30)]
test_large_counts_50  <- test_large_counts[,append(c(1),col_index_50)]
test_large_counts_100 <- test_large_counts[,append(c(1),col_index_100)]
test_large_counts_105 <- test_large_counts[,append(c(1),col_index_105)]
test_large_counts_110 <- test_large_counts[,append(c(1),col_index_110)]


test_large_counts_125 <- test_large_counts[,append(c(1),col_index_125)]
test_large_counts_150 <- test_large_counts[,append(c(1),col_index_150)]
test_large_counts_175 <- test_large_counts[,append(c(1),col_index_175)]
test_large_counts_190 <- test_large_counts[,append(c(1),col_index_190)]
test_large_counts_200 <- test_large_counts[,append(c(1),col_index_200)]
test_large_counts_500 <- test_large_counts[,append(c(1),col_index_500)]
test_large_counts_800 <- test_large_counts[,append(c(1),col_index_800)]
test_large_counts_978 <- test_large_counts

write.csv(test_large_counts_25, file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_25.csv", row.names = FALSE)
write.csv(test_large_counts_30, file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_30.csv", row.names = FALSE)
write.csv(test_large_counts_50, file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_50.csv", row.names = FALSE)
write.csv(test_large_counts_100,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_100.csv",row.names = FALSE)
write.csv(test_large_counts_105,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_105.csv",row.names = FALSE)
write.csv(test_large_counts_110,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_110.csv",row.names = FALSE)


write.csv(test_large_counts_125,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_125.csv",row.names = FALSE)
write.csv(test_large_counts_150,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_150.csv",row.names = FALSE)
write.csv(test_large_counts_175,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_175.csv",row.names = FALSE)
write.csv(test_large_counts_190,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_190.csv",row.names = FALSE)
write.csv(test_large_counts_200,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_200.csv",row.names = FALSE)
write.csv(test_large_counts_500,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_500.csv",row.names = FALSE)
write.csv(test_large_counts_800,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_800.csv",row.names = FALSE)
write.csv(test_large_counts_978,file = "/Users/reyesa1/Desktop/test_large_counts/test_large_counts_978.csv",row.names = FALSE)


samples_test_large_counts_25  <- cols_samples_test_large_counts[col_index_25-1, ]
samples_test_large_counts_30  <- cols_samples_test_large_counts[col_index_30-1, ]
samples_test_large_counts_50  <- cols_samples_test_large_counts[col_index_50-1, ]
samples_test_large_counts_100 <- cols_samples_test_large_counts[col_index_100-1,]
samples_test_large_counts_105 <- cols_samples_test_large_counts[col_index_105-1,]
samples_test_large_counts_110 <- cols_samples_test_large_counts[col_index_110-1,]


samples_test_large_counts_125 <- cols_samples_test_large_counts[col_index_125-1,]
samples_test_large_counts_150 <- cols_samples_test_large_counts[col_index_150-1,]
samples_test_large_counts_175 <- cols_samples_test_large_counts[col_index_175-1,]
samples_test_large_counts_190 <- cols_samples_test_large_counts[col_index_190-1,]
samples_test_large_counts_200 <- cols_samples_test_large_counts[col_index_200-1,]
samples_test_large_counts_500 <- cols_samples_test_large_counts[col_index_500-1,]
samples_test_large_counts_800 <- cols_samples_test_large_counts[col_index_800-1,]
samples_test_large_counts_978 <- cols_samples_test_large_counts

write.csv(samples_test_large_counts_25,  file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_25.csv",  row.names = FALSE)
write.csv(samples_test_large_counts_30,  file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_30.csv",  row.names = FALSE)
write.csv(samples_test_large_counts_50,  file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_50.csv",  row.names = FALSE)
write.csv(samples_test_large_counts_100, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_100.csv", row.names = FALSE)
write.csv(samples_test_large_counts_105, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_105.csv", row.names = FALSE)
write.csv(samples_test_large_counts_110, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_110.csv", row.names = FALSE)


write.csv(samples_test_large_counts_125, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_125.csv", row.names = FALSE)
write.csv(samples_test_large_counts_150, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_150.csv", row.names = FALSE)
write.csv(samples_test_large_counts_175, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_175.csv", row.names = FALSE)
write.csv(samples_test_large_counts_190, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_190.csv", row.names = FALSE)
write.csv(samples_test_large_counts_200, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_200.csv", row.names = FALSE)
write.csv(samples_test_large_counts_500, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_500.csv", row.names = FALSE)
write.csv(samples_test_large_counts_800, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_800.csv", row.names = FALSE)
write.csv(samples_test_large_counts_978, file = "/Users/reyesa1/Desktop/test_large_counts/samples_test_large_counts_978.csv", row.names = FALSE)

## how to get groups/number ot groups same for all rounds of DEA?
