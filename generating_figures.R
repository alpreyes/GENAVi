
##### have to confirm with tiago: versions of featurecounts matrix (different cols btwn ensembleid.csv and all_cell_lines_ordered.csv) #####
##### may not need all_cell_lines_ordered.csv in github/ local repo

all_cell_lines_ordered <- read_csv("all_cell_lines_ordered.csv", col_types = readr::cols(), col_names = TRUE)
raw_counts_all_cols <- read_csv("raw_counts.csv", col_types = readr::cols(), col_names = TRUE)
vst_all_cols <- read_csv("vst.csv", col_types = readr::cols(), col_names = TRUE)

pairwise_DEA_genelist <- read.table("P_10_6_genelist_each_pairwise.txt", header = FALSE, stringsAsFactors = FALSE)
test_unique <- unique(pairwise_DEA_genelist[,1])


all_cell_lines_ordered_DEA_genes <- all_cell_lines_ordered[grep(paste(test_unique, collapse = "|"), all_cell_lines_ordered$Genename),]
raw_counts_all_cols_DEA_genes <- raw_counts_all_cols[grep(paste(test_unique, collapse = "|"), raw_counts_all_cols$Genename),]
vst_all_cols_DEA_genes <- vst_all_cols[grep(paste(test_unique, collapse = "|"), vst_all_cols$Genename),]


##### number of genes that can be found in either version of data table does NOT match number of genes in (unique) list
##### but consensus DEA gene list came from the data table
##### confirm with tiago and michelle



#####just try using tables from app to make heatmaps for figure 1 (resolve data table versions later)

### first vst transform data table
# all_cell_lines_ordered_DEA_genes_vst <- cbind(all_cell_lines_ordered_DEA_genes[,c(1:7)], vst(as.matrix(all_cell_lines_ordered_DEA_genes[,-c(1:7)])))
# test_vst <- vst(as.matrix(all_cell_lines_ordered_DEA_genes[,-c(1:7)]))


##first do expr heatmap

heatmap_expr <- main_heatmap(as.matrix(vst_all_cols_DEA_genes[,-c(1:8)]), name = "Expression", colors = custom_pal_blues) %>%
  add_col_labels(ticktext = colnames(vst_all_cols_DEA_genes[,-c(1:8)])) %>%
  add_row_labels(ticktext = vst_all_cols_DEA_genes$Genename) %>%
  add_col_dendro(hclust(dist(t(as.matrix(vst_all_cols_DEA_genes[,-c(1:8)])))), reorder = TRUE) %>%
  add_row_dendro(hclust(dist(t(as.matrix(vst_all_cols_DEA_genes[,-c(1:8)])))), reorder = TRUE, side = "right") ###looks weird...

heatmap_cor <- main_heatmap(as.matrix(cor(vst_all_cols_DEA_genes[,-c(1:8)], method = "pearson")), name = "Correlation", colors = custom_pal_blues) %>%
  add_col_labels(ticktext = colnames(vst_all_cols_DEA_genes[,-c(1:8)])) %>%
  add_row_labels(ticktext = colnames(vst_all_cols_DEA_genes[,-c(1:8)])) %>%
  add_col_dendro(hclust(as.dist(1-cor(vst_all_cols_DEA_genes[,-c(1:8)], method = "pearson"))), reorder = TRUE) %>%
  add_row_dendro(hclust(as.dist(1-cor(vst_all_cols_DEA_genes[,-c(1:8)], method = "pearson"))), reorder = TRUE, side = "right")
  

