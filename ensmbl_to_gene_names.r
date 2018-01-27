
gencode.rt <- rtracklayer::import("/Users/reyesa1/Desktop/SG04202017_merged_featurecounts/gencode.v24.annotation.gtf", ##change this to the annotation file used in alignment and quantification
                                   format = "gtf", genome = "hg38") ##change this to hg38
library(GenomicRanges)
gencode.table <- data.frame(mcols(gencode.rt), stringsAsFactors=FALSE)
gencode.table <- gencode.table[gencode.table[, "type"] %in% "gene", ]

#h.gene.names <- data.frame(ENSEMBL=rownames(m.gene.rpkm), stringsAsFactors=FALSE) ##don't need this?
#h.gene.names <- dplyr::left_join(m.gene.names, gencode.table, by=c("ENSEMBL" = "gene_id")) ##don't need this?

## actually making vector of gene names according to ensembl ids
gene.names <- c()

for(i in 1:dim(counts3)[1]) #### change counts3 to whatever counts table you are using to generate new vector of corresponding gene names
{
  gene.names[i] <- gencode.table[grep(counts3[i,1],gencode.table[,5]),8] ## change counts3 to whatever counts table you are using to generate new vector of corresponding gene names
}

test <- cbind(gene.names,counts3)
