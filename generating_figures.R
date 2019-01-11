
##### have to confirm with tiago: versions of featurecounts matrix (different cols btwn ensembleid.csv and all_cell_lines_ordered.csv) #####
##### may not need all_cell_lines_ordered.csv in github/ local repo

all_cell_lines_ordered <- read_csv("all_cell_lines_ordered.csv", col_types = readr::cols(), col_names = TRUE)
raw_counts_all_cols <- read_csv("test/raw_counts.csv", col_types = readr::cols(), col_names = TRUE)
vst_all_cols <- read_csv("test/vst.csv", col_types = readr::cols(), col_names = TRUE)

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

test_filter <- all_cell_lines_ordered[rowSums(all_cell_lines_ordered[,-c(1:7)]) > 1,-c(1:7)] ### not same as ensembleid.csv, still diff num of rows
test_filter <- na.omit(test_filter) ### still diff num of rows

View(test_filter)

rm(test_filter)

test_filter <- all_cell_lines_ordered[rowSums(all_cell_lines_ordered[,-c(1:7)]) >1,]
test_filter <- na.omit(test_filter)



#----------------------------------------------------------------------------
#
#     remake figs w/ new counts matrix --> new DE gene list
#
#----------------------------------------------------------------------------
library(tidyverse)
library(DESeq2)
new_fig_matrix_vst <- read_table2("./draft_figures/updated_count_matrix_top_DE_genes.txt")

heatmap_expr <- main_heatmap(as.matrix(new_fig_matrix_vst[,-c(1:8)]), name = "Expression", colors = custom_pal_blues) %>%
  add_col_labels(ticktext = colnames(new_fig_matrix_vst[,-c(1:8)]), font = list(size = 17)) %>%
  add_row_labels(ticktext = new_fig_matrix_vst$Genename, font = list(size = 13)) %>%
  add_col_dendro(hclust(dist(t(as.matrix(new_fig_matrix_vst[,-c(1:8)])))), reorder = TRUE) %>% ### this messes up the col labels
  add_row_dendro(hclust(dist((as.matrix(new_fig_matrix_vst[,-c(1:8)])))), reorder = TRUE, side = "right") ### messes up right side
####### diff than heatmap generated in GENAVi #######
heatmap_expr
heatmap_expr %>% save_iheatmap("genavi_expr_heatmap.pdf")

grid_params <- setup_colorbar_grid(y_spacing = .45)
heatmap_cor <- main_heatmap(as.matrix(cor(new_fig_matrix_vst[,-c(1:8)], method = "pearson")), name = "Correlation", colors = custom_pal_blues, colorbar_grid = grid_params) %>%
  add_col_labels(ticktext = colnames(new_fig_matrix_vst[,-c(1:8)]), font = list(size = 17)) %>%
  add_row_labels(ticktext = colnames(new_fig_matrix_vst[,-c(1:8)]), font = list(size = 17)) %>%
  add_col_dendro(hclust(as.dist(1-cor(new_fig_matrix_vst[,-c(1:8)], method = "pearson"))), reorder = TRUE) %>%
  add_row_dendro(hclust(as.dist(1-cor(new_fig_matrix_vst[,-c(1:8)], method = "pearson"))), reorder = TRUE, side = "right")
heatmap_cor
heatmap_cor %>% save_iheatmap("genavi_cor_heatmap.pdf")

# library(plotly)
# plotly_IMAGE(heatmap_cor, format = "pdf", out_file = "test_genavi_cor_heatmap.pdf")

#----------------------------------------------------------------------------
#
#                         address reviwer 2 comment 4
#
#----------------------------------------------------------------------------
comment_table <- new_fig_matrix_vst[,c(2,9:30)] 
# row.names(comment_table)
# row.names(comment_table) <- comment_table$Genename[1:18]
# row.names(comment_table)
View(comment_table)
comment_table <- comment_table[,-1]
View(comment_table)
row.names(comment_table) 
row.names(comment_table) <- new_fig_matrix_vst$Genename
row.names(comment_table)
View(comment_table)
comment_table <- t(comment_table)
View(comment_table)

hmc_arc <- main_heatmap(as.matrix(cor(comment_table, method = "pearson")), name = "Correlation", colors = custom_pal_blues) %>% # heatmap cor address reviewer comments
  add_col_labels(ticktext = colnames(comment_table)) %>%
  add_row_labels(ticktext = colnames(comment_table)) %>%
  add_col_dendro(hclust(as.dist(1-cor(comment_table, method = "pearson"))), reorder = TRUE) %>%
  add_row_dendro(hclust(as.dist(1-cor(comment_table, method = "pearson"))), reorder = TRUE, side = "right")


hmc_arc ### gene-to-gene correaltion heatmap ACROSS SAMPLES






#---------------------------------------------------------------------------------------------
#
#     make norm gene expr dis for each cell line (new table) and highlight BRCA 1 and 2
#
#---------------------------------------------------------------------------------------------
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(DT)
library(apeglm)
library(readr)
library(edgeR)

data <- read_csv("Cell_Line_RNA_seq_2017_and_2018_FAIL_SAMPLES_REMOVED_featurecounts_matrix.csv")

data_vst <- cbind(data[,1:7], vst(as.matrix(data[,8:29])))

hist(data_vst$IOSE11)

ggplot(data_vst, aes(x = IOSE11)) + 
  geom_histogram() +
  geom_vline(aes(xintercept = data_vst$IOSE11[46251]), color="red", linetype="dashed", size = 1) +
  geom_vline(aes(xintercept = data_vst$IOSE11[36671]), color="blue",linetype="dashed", size = 1)

# macl_vst <- rowMeans(data_vst[,8:29]) ### means across cell lines vst
# ggplot(as.data.frame(macl_vst), aes(x=macl_vst)) + geom_histogram()
# rm(macl_vst)

data_vst$macl <- rowMeans(data_vst[,8:29])

ggplot(data_vst, aes(x=macl)) +
  geom_histogram() +
  geom_vline(aes(xintercept = data_vst$macl[46251]), color="red",  linetype="dashed", size = 1) +
  geom_vline(aes(xintercept = data_vst$macl[36671]), color="blue", linetype="dashed", size = 1)

ggplot(data_vst, aes(x=macl)) +
  geom_violin()

table_brca_vst <- data_vst[c(36671,46251),c(1,2,8:29)]
write.csv(table_brca_vst, file = "/Users/reyesa1/Downloads/table_brca_vst.csv")



range(data_vst$macl)

