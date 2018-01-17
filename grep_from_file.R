
library(stringr)

toMatch <- c("red", "blue", "yellow")

animal <- c("dog", "bird", "fish", "cat", "snake", "otter", "tiger")
color  <- c("green", "purple", "red", "yellow", "red", "blue", "yellow")

mydata <- cbind(animal,color)
mydata <- as.data.frame(mydata)

matches <- unique(grep(paste(toMatch, collapse = "|"), mydata$color, value = TRUE))
matches <- grep(paste(toMatch, collapse = "|"), mydata$color, value = TRUE) ##see how these two lines are different

matched.indices <- grep(paste(toMatch, collapse = "|"), mydata$color) ##grep gives indices automatically

matched.animals <- mydata[matched.indices,1]

test_gene_list <- read.table("test_gene_list.txt",header = FALSE)

paste(test_gene_list$V1, collapse = "|")

grep(paste(test_gene_list$V1, collapse = "|"), final.counts$Genename) ##need to find a way to do exact match...use match()

match(paste(test_gene_list$V1, collapse = "|"), final.counts$Genename)

match(test_gene_list$V1,final.counts$Genename) ##yup use match instead of grep
