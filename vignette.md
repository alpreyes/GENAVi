---
title: "GENAVi"
output:
  html_document:
    toc: true
    toc_depth: 2
    number_sections: true
    theme: united
    highlight: tango
---

# Introduction

In order to improve access to normalization, visualization and analysis of RNA-Seq data we have combined a number of R packages using the Shiny Web App format to create GENAVi: Gene Expression Normalization Analysis and Visualization. This application allows users to browse a provided dataset; a panel of 20 cell lines commonly used in breast and ovarian cancer research, or upload their own dataset as a featureCounts matrix (PMID: 24227677). featureCounts uses read counting for summarization of all reads that overlap any exon for each gene, followed by read annotation with RefSeq or Ensembl.

The application is separated into three tabs:

1. Gene Expression: Users can search or browse a table of provided data or upload their own data, and apply different normalization methods. Genes can be selected in the table using the single field search bar and then selected within the table, or a list of gene symbols can be uploaded and then will be selected within the table. Selected genes will be used for the plotting functions on the next tab. The genes and normalization method selected on this tab will be retained for plotting and Differential Analysis functions.
2. Visualization: Users can select between plotting a histogram of a single gene (Expression plot), or heatmaps (Cluster plots) plotting the Euclidean distance between samples based on either the selected genes from the data table or all genes (which can be toggled between from the Cluster by what genes menu.
3. Differential Expression Analysis: Users can upload a metadata file providing information on sample groups to be compared, run differential analysis based on the DESeq2 package from R, and plot the resulting P values in a customizable volcano plot.

# Uploading User Data/Querying Displayed Data

## Functions

On the Data Expression tab of the app users can upload a featureCounts table as a .csv file via the Data upload option. Large matrices can take several minutes to upload and be processed by the normalization functions, and progress of this upload and normalization can be tracked by the blue progress bar on the bottom right of the app. Once uploaded the users’ data is displayed replacing the provided data table. Desired genes can be selected by searching in the search field at the top right of the data table and clicking on individual genes within the displayed table. Additionally, the user can also upload a list of gene symbols as a .txt file. Selected genes are automatically highlighted and moved to the top of the data table. The user can also download a data table containing only the selected genes by clicking the “Download” option in the top left of the data table. The displayed data table can be sorted by any of the columns using the bi-directional arrows next to the column title. 

## Format of input gene list

In addition to being able to search for and select individual genes within the displayed data table, the user can select multiple genes by uploading a list of gene symbols or Ensembl ID’s as a .txt file. This text file should be formatted such that each line in the file is a gene symbol or Ensemble ID. Once the input gene list is uploaded, the corresponding genes within the displayed data table are automatically selected and can be used in the visualization tab of the application.


## Format of input count matrix

In addition to being able to query the default data set- a panel of 20 breast and ovarian cancer cell lines, the user can also upload their own gene expression data set in the form of count a matrix produced by the featurecounts function in the subread package.
the addgeneinfo() function reads in count matrix in csv format. it checks if the first column is an ENSEMBLE ID or gene name. If the first column of the inputted count matrix is the proper format, addgeneinfo() adds the genes metadata from gencode (by default uses hg38 build)

# Transformation and Normalization of Gene Expression Data (Tab 1)

## Background

i. Variability in RNA-seq count data usually increases with the mean, indicating that genes that are more highly expressed can often also have a greater variation; this type of data is called heteroskedastic. Expression data in a raw count format causes cluster analysis such as PCA, hierarchical clustering, or k-means analysis to be driven by genes or features with the highest expression in the form of high raw count values. This is because in its raw count format, these highly expressed genes/features are also the most variable. To improve cluster analysis we need to convert this expression data to a homoskedastic form wherein genes/features have the same variance even for a higher range of the mean. One method to impose homoskedasticity is the log2 transformation (log2(count + pseudocount)). While this method does prevent the most highly expressed genes/features from overpowering cluster analysis it can also give undue weight to the genes/features with the lowest counts depending on the choice of pseudocount. To understand this, consider this example; with the choice of 1 as the pseudocount in the case of log2 transforming a count value of 2 (2 reads mapped to given gene/feature in given sample), the resulting expression value is 1.584. Whereas with the same choice of pseudocount in the case of log2 transforming a count value of 200, the resulting expression value is 7.651. In this comparison, a difference in raw count expression that is two orders of magnitude is compressed to a difference of 6.067 under the log2 transform.
ii. Another problem that arises with RNA-seq count data is overdispersion- the phenomenon wherein the observed variance within a dataset is greater than the expected variance under a chosen statistical model. Previously, the Poisson distribution has been used to model RNA-seq count data. However, a property of the Poisson distribution is that the mean and variance are equal which prevents the model from accurately modeling RNA-seq count data which is often heteroskedastic. A more appropriate statistical model is the negative binomial model used in the edgeR package as the DESeq2 package. This is an extension of the Poisson distribution that has different parameters allowing for separate modeling of mean and variance. Both edgeR and DESeq2 assume that the expected value of counts for a gene in a given sample can be modeled using this distribution.
iii. The user can toggle between transformation options by choosing from the “Select Transform” dropdown menu on left of the screen.


## raw counts

i. The default version of expression data displayed in GENAVi is raw count data. This is the measure of expression produced by featurecounts. These raw count measures are integer values counting the number of reads mapped to a feature in the reference genome to which fastq files were aligned in a specific sample. Raw count data must be normalized to account for differences in sequencing depth between samples so that the user can make meaningful comparisons of gene expression levels across different samples.

## Row normalized

i. The “row normalized” transformation option is applied to raw count data similar to a Z-score normalization. For each gene in the displayed dataset, the mean expression and standard deviation across all samples/columns are calculated. The raw count expression values for each gene/feature are then transformed by subtracting the mean and scaling by the standard deviation across all samples. A custom function was implemented to perform row normalization on each gene of the displayed data table.
ii. When opening GENAVi and/or uploading your own count matrix, the application acquires metadata matching the featurecounts table (default or uploaded) and calculates the row normalized transformation and saves it as an object to avoid calculating it each time the “row normalized” option is selected.
iii. This row normalization rescales the displayed data so that the gene expression is reported in units of standard deviation away from the mean expression.

## Log-counts-per-million: logCPM

i. The logCPM transformation implemented in GENAVi through the cpm() function comes from the edgeR package.
ii. This transformation is calculated in two steps. First, counts per million reads are calculated from raw counts by scaling the number of reads overlapping a feature in a specific sample by the total number of reads in that sample (library size). This is performed for every feature in every sample. These CPM values are then log2 transformed with the log=TRUE argument in the cpm() function.
iii. When opening GENAVi and/or uploading your own count matrix, the application acquires metadata matching the featurecounts table (default or uploaded) and calculates the logCPM transform and saves it as an object to avoid calculating it each time the option is selected (implemented in the cpm() function in edgeR).
iv. The logCPM transformation allows the comparison of gene expression across separate samples by accounting for differences in sequencing depth which allows for meaningful comparisons of gene expression across samples by accounting for between-sample variability.

## Variance stabilizing transformation- vst

i. The variance stabilizing transformation, implemented in the DESeq2 package, attempts to address the issue of overdispersion of RNA-seq count data.
ii. The VST transformation produces transformed expression measures similar to those produced by the log2 transform when used on genes with high counts.
iii. Similar to the edgeR package, the VST transformation also assumes a negative binomial model to allow for a change in variance based on the dynamic range of the mean.
iv. For entire dataset, VST transformed counts approaches homoskedasticity and can be used directly for clustering. 
v. VST is much faster to compute than rlog and is not as sensitive to counts with high values unlike the rlog function. So it is more appropriate to use VST to transform larger datasets (up to hundreds of samples/columns) instead of the rlog transformation for example.
vi. When opening GENAVi and/or uploading your own count matrix, the application acquires metadata matching the featurecounts table (default or uploaded) and calculates the vst transform and saves it as an object to avoid calculating it each time the option is selected (implemented in the rlog function from DESeq2)

## Regularized-logarithm transformation- rlog (add warning about number of samples because compute time)

i. To address this issue, GENAVi uses the regularized logarithm (rlog) transform implemented in DESeq2.
ii. The rlog transform performs similar to the log2 transform when applied to genes/features with high raw count expression. When rlog is applied to lower count values, the resulting transformed values are more shrunken towards the mean value of a gene’s expression across all samples in data. Rlog transformed values are calculated by converting the raw count data to the log2 scale and then defining a size factor for each sample which accounts for differences in sequencing depth between samples. This approach is referred to as shrinkage.
iii. This shrinkage approach can affect raw count expression of genes/features to a varying degree. Therefore the rank or order of genes based on expression may not be conserved.
iv. Rlog transformed counts approach homoskedasticity and can be used directly for clustering. 
v. Because the rlog function calculates a shrinkage term for each sample it takes the longest to compute (in GENAVi), so it is best used for small data sets (fewer samples/columns).
vi. When opening GENAVi and/or uploading your own count matrix, the application acquires metadata matching the featurecounts table (default or uploaded) and calculates the rlog transform and saves it as an object to avoid calculating it each time the option is selected (implemented in the rlog function from DESeq2).

# Visualization (tab 2)

## Barplot

i. If the user wants to view the expression of a single gene/feature across all samples, GENAVi displays a simple barplot. The user can search the displayed data set under any transformation for a particular gene of interest and select it to be displayed in the barplot in the Visualization tab.
ii. If more than one gene is selected from the data table, the simple barplot in the Visualization tab is disabled and an interactive heatmap is displayed instead.

## Expression heatmap

GENAVi uses the R package iheatmapr to generate interactive heatmaps to visualize gene expression information. All capabilities of manipulating these heatmaps that iheatmapr provides are available in GENAVi. The user can toggle over cells in the displayed heat map to see row and column identification information, subset the heatmap to zoom in on cells of interest, and download the heatmap by clicking on the camera icon in the top right corner.

i. When more than one gene/feature from the displayed data is selected, an interactive heatmap is displayed in place of the barplot to visualize, cluster, and compare gene/ expression across samples.
    
ii. To generate the gene expression heatmap, the user must first select a version of the data table by selecting from the “Select Transform” dropdown menu on the left side of the application. See the above section on properties and advantages of specific transformations. It is worth noting that visualizing the “raw counts” version of an RNA-seq data set may make the heatmap color scale difficult to interpret because of the magnitude of the range of gene expression in raw counts. Next, the user can select genes/features (either by searching for and clicking on individual genes or uploading a gene list test file). Once the desired genes/features are selected in the displayed data table, the user can navigate to the Visualization tab and select the “Expression plots” subtab to see the resulting expression heatmap. 
    
iii. The expression heatmap displays the expression matrix of the user-selected genes across all samples of the displayed dataset under the selected version of the data table.

iv. dendrogram/hierarchical clustering of samples
    1. The user may see that the order of the samples (columns) and selected genes/features in the expression heatmap is not the same as the order in the displayed data table or the order of selection. 
    2. This is due to the add_col_dendro() and add_row_dendro() options in the iheatmapr package. These options order the columns and rows of the expression heatmap by similarity of expression of the selected genes/rows. This similarity is measured by calculating the Euclidian distance between the samples of the gene expression matrix by applying the dist() function directly to the gene expression information.

## Correlation heatmap (Clustering heatmap)

i. In addition to viewing the expression matrix of selected genes, the user can also view how the expression of selected genes affects the similarity between samples measured through correlation. By selecting “Clustering plots” subtab, the user can view a heatmap representing the the Pearson correlation matrix calculated from the gene expression matrix. The correlation matrix is calculated using the cor() function with the option method=”pearson” option which produces a matrix containing the pairwise pearson correlations between the columns of the displayed data table. When viewing the euclidian distance heatmap, the user can select to either cluster samples based on the entire expression data set by selecting “All genes” from the dropdown menu in the top left of the tab. The user can also choose “Selected genes” to cluster samples using the gene expression information only from the selected genes/features.

ii. dendrogram/hierarchical clustering of samples
  1. You will notice that the order of the samples (columns) in the clustering heatmap is not the same as the order in the displayed data table. 
  2. This is due to the add_col_dendro() and add_row_dendro() options in the iheatmapr package. Unlike the expression heatmap where these options order columns and rows of the data by similarity of expression of the selected genes/rows by calculating the Euclidian distance between the samples. based on the Euclidian distance calculated directly from the gene expression matrix, the correlation heatmap orders columns and rows by calculating the Euclidian distance between samples from the correlation matrix.


# Differential Expression Analysis
