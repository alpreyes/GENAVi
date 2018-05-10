
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

# all_cell_lines_ordered_DEA_genes_vst <- cbind(all_cell_lines_ordered_DEA_genes[,c(1:7)], vst(as.matrix(all_cell_lines_ordered_DEA_genes[,-c(1:7)])))
# test_vst <- vst(as.matrix(all_cell_lines_ordered_DEA_genes[,-c(1:7)]))


##first do expr heatmap
heatmap_expr <- main_heatmap(as.matrix(vst_all_cols_DEA_genes[,-c(1:8)]), name = "Expression", colors = custom_pal_blues) %>%
  add_col_labels(ticktext = colnames(vst_all_cols_DEA_genes[,-c(1:8)])) %>%
  add_row_labels(ticktext = vst_all_cols_DEA_genes$Genename) %>%
  add_col_dendro(hclust(dist(t(as.matrix(vst_all_cols_DEA_genes[,-c(1:8)])))), reorder = TRUE) %>%
  add_row_dendro(hclust(dist(t(as.matrix(vst_all_cols_DEA_genes[,-c(1:8)])))), reorder = TRUE, side = "right") ###looks weird...

##now do cor heatmap
heatmap_cor <- main_heatmap(as.matrix(cor(vst_all_cols_DEA_genes[,-c(1:8)], method = "pearson")), name = "Correlation", colors = custom_pal_blues) %>%
  add_col_labels(ticktext = colnames(vst_all_cols_DEA_genes[,-c(1:8)])) %>%
  add_row_labels(ticktext = colnames(vst_all_cols_DEA_genes[,-c(1:8)])) %>%
  add_col_dendro(hclust(as.dist(1-cor(vst_all_cols_DEA_genes[,-c(1:8)], method = "pearson"))), reorder = TRUE) %>%
  add_row_dendro(hclust(as.dist(1-cor(vst_all_cols_DEA_genes[,-c(1:8)], method = "pearson"))), reorder = TRUE, side = "right")
  
##as of 5/9/18 volcano plot in draft figure panel was generated in GENAVi, now try to generate here to control dimensions and resolution






########### michelle's original pairwise DEA results gene list used pval instead of padj.....so now use new correct gene list

pairwise_DEA_results <- read_csv("PrecursorNormal_vs_CancerLines.csv", col_types = readr::cols(), col_names = TRUE) ###contains a bunch of empty rows at the end, just use first 66 rows
pairwise_DEA_results <- pairwise_DEA_results[c(1:66),]

## as of 5/9/18 biomaRt and GENCODE version issue in GENAVi causing not all gene names to be found, so must use BOTH name1 and name2 cols for pairwise_DEA_gene_list used for search and figures

pairwise_DEA_genelist <- append(pairwise_DEA_results[[1]],pairwise_DEA_results[[2]])
##now remove NA's from name2 column
pairwise_DEA_genelist <- as.character(na.omit(pairwise_DEA_genelist))

###have DEA gene search list, now search downloaded data table from GENAVi, still have to resolve data table versions
vst_all_cols_DEA_genes <- vst_all_cols[grep(paste(pairwise_DEA_genelist, collapse = "|"), vst_all_cols$Genename),]
x <- unique(pairwise_DEA_genelist) ##still can't find all the gene names.....why?



##first do expr heatmap
heatmap_expr <- main_heatmap(as.matrix(vst_all_cols_DEA_genes[,-c(1:8)]), name = "Expression", colors = custom_pal_blues) %>%
  add_col_labels(ticktext = colnames(vst_all_cols_DEA_genes[,-c(1:8)])) %>%
  add_row_labels(ticktext = vst_all_cols_DEA_genes$Genename, font = list(size = 6)) %>%
  add_col_dendro(hclust(dist(t(as.matrix(vst_all_cols_DEA_genes[,-c(1:8)])))), reorder = TRUE) %>%
  add_row_dendro(hclust(dist(t(as.matrix(vst_all_cols_DEA_genes[,-c(1:8)])))), reorder = TRUE, side = "right") ###looks weird...

##now do cor heatmap
heatmap_cor <- main_heatmap(as.matrix(cor(vst_all_cols_DEA_genes[,-c(1:8)], method = "pearson")), name = "Correlation", colors = custom_pal_blues) %>%
  add_col_labels(ticktext = colnames(vst_all_cols_DEA_genes[,-c(1:8)])) %>%
  add_row_labels(ticktext = colnames(vst_all_cols_DEA_genes[,-c(1:8)])) %>%
  add_col_dendro(hclust(as.dist(1-cor(vst_all_cols_DEA_genes[,-c(1:8)], method = "pearson"))), reorder = TRUE) %>%
  add_row_dendro(hclust(as.dist(1-cor(vst_all_cols_DEA_genes[,-c(1:8)], method = "pearson"))), reorder = TRUE, side = "right")


write.table(pairwise_DEA_genelist, "pairwise_DEA_genelist.txt", row.names = FALSE) ##had to manipulate file to input into GENAVi

test_filter <- all_cell_lines_ordered[rowSums(all_cell_lines_ordered[,-c(1:7)]),-c(1:7)] ### not same as ensembleid.csv, still diff num of rows
test_filter <- na.omit(test_filter) ### still diff num of rows

test_filter




