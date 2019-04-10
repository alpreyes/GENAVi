library(TCGAbiolinks)
library(SummarizedExperiment)
proj <- "TCGA-BRCA"
dir.create(proj,showWarnings = FALSE)
query <- GDCquery(project = proj,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query,files.per.chunk = 100)
data <- GDCprepare(query = query,
                   save = FALSE)
data <- data[,data$subtype_BRCA_Subtype_PAM50 %in% c("Basal","LumA")]
data <- data[,!is.na(data$race)]
data <- data[,data$race == "black or african american"]
samples <- cbind(rownames(colData(data)), colData(data))
readr::write_csv(as.data.frame(samples),path = file.path(proj,"metadata.csv"))
readr::write_csv(cbind(data.frame("GeneID" = rownames(assay(data))),as.data.frame(assay(data))),path = file.path(proj,"genavi_raw_counts_hg38_TCGA.csv"))
