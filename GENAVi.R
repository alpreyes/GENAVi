
library(shiny)
library(tidyr) ##use this for barplot
library(ggplot2) ##use this for barplot
library(DESeq2)
library(iheatmapr)

final.counts <- read.csv("final.counts.csv",header = TRUE) ##read in data going to view
count.matrix <- final.counts[,c(1,8:25)] ##subset data for functions
#final.RPKM   <- read.csv("final.RPKM.csv",header = TRUE) ##read in data going to view
#RPKM.matrix  <- final.RPKM[,c(1,8:25)] ##subset data for functions
meh <- gather(count.matrix,"cell_line", "value", 2:19) ##use these for barplot but may be deleting soon, will have to change name
meh2 <- data.frame() ##will have to change name to something that makes sense in rmarkdown, need this for making reactive barplot

final.rownorm <- read.csv("final.rownorm.csv",header = TRUE)
cmf.rownorm <- read.csv("cmf.rownorm.csv",header = TRUE) ##will need to have gene names and other first data columns to id genes
cmf.rownorm <- as.matrix(cmf.rownorm)
rownorm.matrix <- final.rownorm[,c(1,8:25)]
meh.rownorm <- data.frame()
meh.rownorm <- gather(rownorm.matrix,"cell_line","value",2:19) ##now use this for barplot
##cmf data tables need gene names and other identifiers, will fix this when reading in data tables (will calc manually instead of reading in objects)

#counts.rlog <- rlog(as.matrix(count.matrix[,2:19]))
#write.csv(counts.rlog, "counts.rlog.csv",row.names = FALSE)
#rm(counts.rlog)
counts.rlog <- read.csv("counts.rlog.csv",header = TRUE)

final.rlog <- cbind(final.counts[,1:7],counts.rlog) ##might change to just reading in final.rlog matrix with read.csv (after writing final versions of tables)
rlog.matrix <- final.rlog[,c(1,8:25)]
meh.rlog   <- data.frame()
meh.rlog <- gather(rlog.matrix,"cell_line","value",2:19)

#counts.vst <- vst(as.matrix(count.matrix[,2:19]))
#write.csv(counts.vst, "counts.vst.csv",row.names = FALSE)
#rm(counts.vst)
counts.vst <- read.csv("counts.vst.csv",header = TRUE)

final.vst <- cbind(final.counts[,1:7],counts.vst)
vst.matrix  <- final.vst[,c(1,8:25)]
meh.vst    <- data.frame()
meh.vst <- gather(vst.matrix,"cell_line","value",2:19)


#matrix_expr <- matrix()

transforms <- c("raw counts", "row normalized", "rlog", "vst") ##adding "t-score" -> row normalized, ask michelle about normalizing with FTSEC lines 
cell.line.clusters <- c("-no selection-","filtered genes", "selected genes") ##do this for cell line cluster heatmaps, changed all to filtered, should "all" be an option?

###objects for cell line clustering
count.matrix.filtered <- count.matrix[rowSums(count.matrix[,-1]) > 1,] ##treat this as if it were "cmf"
cmf.rlog <- rlog(as.matrix(count.matrix.filtered[,-1]),blind = FALSE)
cmf.vst  <- vst(as.matrix(count.matrix.filtered[,-1]),blind = FALSE)
####will have to add one for zscore when zscore is finally implemented####

ui <- fluidPage(tabsetPanel(
  tabPanel("Expression", ##changing from tab 1, but still using tab1 in the other parts of code
           h2('GENAVi'), 
           #selectInput("select", "Select Columns", names(final.counts), multiple = TRUE),
           selectInput("select_tab1", "Select Transform", transforms, multiple = FALSE), ##need individual selectInputs for each tab
           fileInput("input_gene_list_tab1", "Input Gene List (Optional)", multiple = FALSE, accept = NULL, width = NULL, buttonLabel = "Browse", placeholder = "No file selected"), ##what is accept parameter???
           #actionButton("but_sortSelectedFirst_tab1", "Selected Rows First"), ##do this to put selected rows at top of data table, trying it out
           DT::dataTableOutput('tbl.tab1'),
           verbatimTextOutput("warning_message_tab1", placeholder = FALSE), ##trying this out, 
           #numericInput("num","Input row index",1,min = 1,max = 60554),
           #verbatimTextOutput("warning_message"), ##may need to take this out for  now
           #DT::dataTableOutput('tbl.tab1.selected'), ##doing this to keep track of selected rows, take this out
           plotOutput("barplot"), 
           iheatmaprOutput("heatmap_expr") ##expr for expression
  ),
  tabPanel("Clustering", ##changing from tab 2, but still usibg tab2 in other parts of code
           h2("GENAVi"), 
           selectInput("select_tab2", "Select Transform", transforms, multiple = FALSE),
           fileInput("input_gene_list_tab2", "Input Gene List (Optional)", multiple = FALSE, accept = NULL, width = NULL, buttonLabel = "Browse", placeholder = "No file selected"),
           #actionButton("but_sortSelectedFirst_tab2", "Selected Rows First"), ##repeat in tab2
           DT::dataTableOutput('tbl.tab2'),
           #DT::dataTableOutput('tbl.tab2.selected'),
           selectInput("select_clus", "Cluster by what genes", cell.line.clusters, multiple = FALSE),
           verbatimTextOutput("warning_message_tab2", placeholder = FALSE), ##do this to guide making heatmap in tab2
           iheatmaprOutput("heatmap_clus") ##clus for clusters
  )
)
)



#### sample ui code for separate tabs

#ui <- fluidPage(tabsetPanel(
#  tabPanel("tab 1", 
#           h2('Ovarian Cancer Cell Lines'), 
#           #selectInput("select", "Select Columns", names(final.counts), multiple = TRUE),
#           selectInput("select", "Select Transform", transforms, multiple = FALSE),
#           DT::dataTableOutput('tbl'),
#           #numericInput("num","Input row index",1,min = 1,max = 60554),
#           #verbatimTextOutput("value"),
#           plotOutput("barplot")
#  ),
#  tabPanel("tab 2",
#           iheatmaprOutput("heat") 
#  ))
#)

####

server <- function(input,output) 
{
  #output$tbl = DT::renderDataTable({
  # columns = names(final.counts)
  #if(!is.null(input$select)) {
  # columns = input$select
  #}
  #final.counts[,columns,drop=FALSE]
  #})
  
  # if(input$select == "raw counts") {output$tbl = DT::renderDataTable({final.counts})}
  #if(input$select == "rlog") {output$tbl = DT::renderDataTable({final.rlog})}
  #if(input$select == "vst") {output$tbl = DT::renderDataTable({final.vst})}
  
  # output$tbl = DT::renderDataTable({final.counts})
  
  #output$tbl = DT::renderDataTable({final.rlog})
  
  #table.counts <- final.counts #reactive({(final.counts)}) ##might not need to do this step at all, just change "table"s in later commands
  #table.rlog   <- final.rlog #reactive({(final.rlog)})
  #table.vst    <- final.vst #reactive({(final.vst)})
  
  ###might have to change how to render selected data table section bc rows dont stay selected
  #output$tbl.tab1 <-  DT::renderDataTable({
  # if(input$select_tab1 == "rlog") tbl.tab1 <- final.rlog#table.rlog #DT::datatable(table.rlog)
  #  if(input$select_tab1 == "vst") tbl.tab1 <- final.vst#table.vst #DT::datatable(table.vst)
  # if(input$select_tab1 == "raw counts") tbl.tab1 <- final.counts#table.counts #DT::datatable(table.counts)
  #  if(input$select_tab1 == "row normalized") tbl.tab1 <- final.rownorm
  # tbl.tab1
  #})
  
  output$tbl.tab1 <-  DT::renderDataTable({
    
    if(input$select_tab1 == "rlog") tbl.tab1 <- final.rlog#table.rlog #DT::datatable(table.rlog)
    if(input$select_tab1 == "vst") tbl.tab1 <- final.vst#table.vst #DT::datatable(table.vst)
    if(input$select_tab1 == "raw counts") tbl.tab1 <- final.counts#table.counts #DT::datatable(table.counts)
    if(input$select_tab1 == "row normalized") tbl.tab1 <- final.rownorm
    
    selected_rows <- input$tbl.tab1_rows_selected
    
    status <- factor("Unselected",levels = c("Unselected","Selected"))
    tbl.tab1 <- cbind(status,tbl.tab1)
    tbl.tab1$status[selected_rows] <- 'Selected'
    
    tbl.tab1 %>% 
      DT::datatable(extensions = c('Buttons',"FixedHeader"),
                    class = 'cell-border stripe',
                    selection = list(mode = "multiple", target= 'row', selected = selected_rows),
                    options = list(dom = 'Blfrtip',
                                   order = c(1,"desc"),
                                   deferRender = TRUE,
                                   paging = T,
                                   buttons =
                                     list('colvis', list(
                                       extend = 'collection',
                                       buttons = list(list(extend='csv',
                                                           filename = "Cell"),
                                                      list(extend='excel',
                                                           filename = "Cell"),
                                                      list(extend='pdf',
                                                           title = "",
                                                           filename= "Cell")),
                                       text = 'Download'
                                     )),
                                   fixedHeader = TRUE,
                                   pageLength = 20,
                                   scrollX = TRUE,
                                   lengthMenu = list(c(10, 20, -1), c('10', '20', 'All')) ##might have to change this
                    ),
                    filter = 'top'
      )
  }) ##works to get selected rows on top but fucks up if select more than one at a time.....worry about it later
  
  #output$tbl.tab1 <-  reactiveValues({ ##trying this to get reordering rows to work...
  #if(input$select_tab1 == "rlog") tbl.tab1 <- final.rlog#table.rlog #DT::datatable(table.rlog)
  #if(input$select_tab1 == "vst") tbl.tab1 <- final.vst#table.vst #DT::datatable(table.vst)
  #if(input$select_tab1 == "raw counts") tbl.tab1 <- final.counts#table.counts #DT::datatable(table.counts)
  #if(input$select_tab1 == "row normalized") tbl.tab1 <- final.rownorm
  #tbl.tab1
  #})
  
  ##display selected rows separately in another data table object
  #output$tbl.tab1.selected <- DT::renderDataTable({
  #if(input$select_tab1 == "rlog") {tbl.tab1.selected <- final.rlog[c(input$tbl.tab1_rows_selected),]}
  #if(input$select_tab1 == "vst") {tbl.tab1.selected <- final.vst[c(input$tbl.tab1_rows_selected),]}
  #if(input$select_tab1 == "raw counts") {tbl.tab1.selected <- final.counts[c(input$tbl.tab1_rows_selected),]}
  #if(input$select_tab1 == "row normalized") tbl.tab1.selected <- final.rownorm[c(input$tbl.tab1_rows_selected),]
  #tbl.tab1.selected
  #})
  
  #output$tbl.tab2 <-  DT::renderDataTable({
  #if(input$select_tab2 == "rlog") tbl.tab2 <- final.rlog#table.rlog #DT::datatable(table.rlog)
  #if(input$select_tab2 == "vst") tbl.tab2 <- final.vst#table.vst #DT::datatable(table.vst)
  #if(input$select_tab2 == "raw counts") tbl.tab2 <- final.counts#table.counts #DT::datatable(table.counts)
  #if(input$select_tab2 == "row normalized") tbl.tab2 <- final.rownorm
  #tbl.tab2
  #})
  
  output$tbl.tab2 <-  DT::renderDataTable({
    
    if(input$select_tab2 == "rlog") tbl.tab2 <- final.rlog#table.rlog #DT::datatable(table.rlog)
    if(input$select_tab2 == "vst") tbl.tab2 <- final.vst#table.vst #DT::datatable(table.vst)
    if(input$select_tab2 == "raw counts") tbl.tab2 <- final.counts#table.counts #DT::datatable(table.counts)
    if(input$select_tab2 == "row normalized") tbl.tab2 <- final.rownorm
    
    selected_rows <- input$tbl.tab2_rows_selected ###may need to label this for tab1 and tab2...nope don't need to
    
    status <- factor("Unselected",levels = c("Unselected","Selected"))
    tbl.tab2 <- cbind(status,tbl.tab2)
    tbl.tab2$status[selected_rows] <- 'Selected'
    
    tbl.tab2 %>% 
      DT::datatable(extensions = c('Buttons',"FixedHeader"),
                    class = 'cell-border stripe',
                    selection = list(mode = "multiple", target= 'row', selected = selected_rows),
                    options = list(dom = 'Blfrtip',
                                   order = c(1,"desc"),
                                   deferRender = TRUE,
                                   paging = T,
                                   buttons =
                                     list('colvis', list(
                                       extend = 'collection',
                                       buttons = list(list(extend='csv',
                                                           filename = "Cell"),
                                                      list(extend='excel',
                                                           filename = "Cell"),
                                                      list(extend='pdf',
                                                           title = "",
                                                           filename= "Cell")),
                                       text = 'Download'
                                     )),
                                   fixedHeader = TRUE,
                                   pageLength = 20,
                                   scrollX = TRUE,
                                   lengthMenu = list(c(10, 20, -1), c('10', '20', 'All'))
                    ),
                    filter = 'top'
      )
  })
  
  #output$tbl.tab2.selected <- DT::renderDataTable({
  # if(input$select_tab2 == "rlog") {tbl.tab2.selected <- final.rlog[c(input$tbl.tab2_rows_selected),]}
  #if(input$select_tab2 == "vst") {tbl.tab2.selected <- final.vst[c(input$tbl.tab2_rows_selected),]}
  #if(input$select_tab2 == "raw counts") {tbl.tab2.selected <- final.counts[c(input$tbl.tab2_rows_selected),]}
  #if(input$select_tab2 == "row normalized") tbl.tab2.selected <- final.rownorm[c(input$tbl.tab2_rows_selected),]
  #tbl.tab2.selected
  #})
  
  
  
  #output$tbl = DT::renderDataTable({if(input$select == "vst") final.vst})
  
  ##print error message(s): if number of selected genes is NOT equal to exactly one
  
  ### tbl.tab1 object not found in app???
  #warning_message <- renderPrint({
  # if(nrow(tbl.tab1[input$tbl.tab1_rows_selected,]) != 1) expr = "Barplot is only displayed with one selected gene"
  #})
  #warning_message() ##this breaks the app???
  
  #warning_message <- "Barplot is only displayed when exactly one gene is selected" ##putting in top of code doesn't work either
  
  #output$warning_message <- renderText({ ##making this conditional did not work
  # if(nrow(matrix_expr[c(input$tbl.tab1_rows_selected),]) != 1) warning_message <- "Barplot is only displayed when exactly one gene is selected"
  #warning_message
  #})
  
  
  ################### making selected rows first in tabl.tab1...currently breaks the whole app ###################
  
  #test_data <- reactiveValues(table_data = tbl.tab1) ##trying this out, may not need this bc i render the data object above, possibly dont need this
  
  #observeEvent(input$but_sortSelectedFirst_tab1, ##adding just this part within the ({}) works but breaks when you click the button
  #            {
  #             ## selected row info is stored in var input$tbl.tab1_rows_selected
  #            selected_rows <- input$tbl.tab1_rows_selected
  #           ##alc new row order w/ selected rows on top
  #          row_order <- order(seq_along(output$tbl.tab1[[1]]) %in% selected_rows, decreasing = TRUE)
  #       })
  
  #output$tbl.tab1 <- output$tbl.tab1[row_order,] ##trying this out, including this line breaks the app
  
  #proxy <- DT::dataTableProxy('tbl.tab1')
  #DT::replaceData(proxy, output$tbl.tab1)
  ## make sure to select rows again
  #DT::selectRows(proxy, seq_along(selected_rows))
  
  #output$tbl.tab1 <- DT::renderDataTable(isolate(output$tbl.tab1))
  
  ### start manipulating input file to grep for and select genes ###
  #output$input_gene_list_tab1 <- DT::dataTableProxy({ ##this breaks the app
  
  # input_gene_list_tab1 <- input$input_gene_list_tab1 ##see if this breaks the app
  
  #})
  
  output$warning_message_tab1 <- renderText({
    if(length(input$tbl.tab1_rows_selected) > 1)
    {
      return("Barplot is only displayed when exactly one gene is selected")
    }
    else
    {
      return(NULL)
    }
  })
  
  
  output$barplot <- renderPlot({
    
    if(is.null(input$tbl.tab1_rows_selected)) {return(NULL)} ##may need to put this in heatmap section and in tab2
    if(length(input$tbl.tab1_rows_selected) > 1) {return(NULL)}
    
    if(input$select_tab1 == "raw counts") ##try adding condition here...we'll see
    {
      for(i in seq(input$tbl.tab1_rows_selected,1089972,60554)) ##1089972 is 18 repetitions of 60554 change to vars instead of numbers later, will have to change for filetered data???
      {
        meh2 <- rbind(meh2,meh[i,]) ##will have to change name of meh and meh2
      }
      barplot <- ggplot(meh2, aes(x=cell_line, y=value)) + geom_bar(stat = "identity")
    }  
    
    if(input$select_tab1 == "rlog")
    {
      for(i in seq(input$tbl.tab1_rows_selected,1089972,60554))
      {
        meh2 <- rbind(meh2,meh.rlog[i,])
      }
      barplot <- ggplot(meh2, aes(x=cell_line, y=value)) + geom_bar(stat = "identity")
    }
    
    if(input$select_tab1 == "vst")
    {
      for(i in seq(input$tbl.tab1_rows_selected,1089972,60554))
      {
        meh2 <- rbind(meh2,meh.vst[i,])
      }
      barplot <- ggplot(meh2, aes(x=cell_line, y=value)) + geom_bar(stat = "identity")
    }
    
    if(input$select_tab1 == "row normalized")
    {
      for(i in seq(input$tbl.tab1_rows_selected,1089972,60554))
      {
        meh2 <- rbind(meh2,meh.rownorm[i,])
      }
      barplot <- ggplot(meh2, aes(x=cell_line, y=value)) + geom_bar(stat = "identity") ##this is weird bc full data table of rownorm has NA's, how to address???
    }
    
    barplot
    
  })
  
  output$heatmap_expr <- renderIheatmap({ ##expr for expression other barplot will have clus
    
    if(is.null(input$tbl.tab1_rows_selected)) {return(NULL)} ##necessary???
    
    if(input$select_tab1 == "raw counts") matrix_expr <- count.matrix ##first instance of matrix_expr...why cant app find it???
    if(input$select_tab1 == "rlog") matrix_expr <- rlog.matrix
    if(input$select_tab1 == "vst") matrix_expr <- vst.matrix
    if(input$select_tab1 == "row normalized") matrix_expr <- rownorm.matrix
    
    ##may need to change order of cell lines from default alphabetic to histotype specific???...do that with dendro???
    heatmap_expr <- main_heatmap(as.matrix(matrix_expr[c(input$tbl.tab1_rows_selected),-1])) %>%
      add_col_labels(ticktext = colnames(matrix_expr[,-1])) %>%
      add_row_labels(ticktext = matrix_expr[input$tbl.tab1_rows_selected,1]) %>% ##trying to add dendro
      add_col_dendro(hclust(dist(t(matrix_expr[input$tbl.tab1_rows_selected,-1])))) ##may have to take out -1 to avoid losing 1st data col
    if(nrow(matrix_expr[input$tbl.tab1_rows_selected,] ) > 1) ##currently still trying to cluster genes selected
    {heatmap_expr <- heatmap_expr %>% add_row_dendro(hclust(dist((matrix_expr[c(input$tbl.tab1_rows_selected),]))), reorder = TRUE, side = "right")} ##taking out t() works but still has to be there...see DESeq2 workflow
    #add_row_dendro(hclust(as.matrix(matrix_expr[c(input$tbl.tab1_rows_selected),-1]))) %>% ##trying row dendro first to see possible clustering???
    #add_col_dendro(hclust(as.matrix(matrix_expr[c(input$tbl.tab1_rows_selected),-1])))
    
    #if(input$select_tab1 == "raw counts")
    #{
    # heatmap_expr <- main_heatmap(as.matrix(count.matrix[c(input$tbl.tab1_rows_selected),-1])) %>%
    #  add_col_labels(ticktext = colnames(count.matrix[,-1])) %>%
    # add_row_labels(ticktext = count.matrix[input$tbl.tab1_rows_selected,1])
    #}
    
    #if(input$select_tab1 == "rlog")
    #{
    # heatmap_expr <- main_heatmap(as.matrix(rlog.matrix[c(input$tbl.tab1_rows_selected),-1])) %>%
    #  add_col_labels(ticktext = colnames(rlog.matrix[,-1])) %>% 
    # add_row_labels(ticktext = rlog.matrix[input$tbl.tab1_rows_selected,1])
    #}
    
    #if(input$selec_tab1 == "vst")
    #{
    # heatmap_expr <- main_heatmap(as.matrix(vst.matrix[c(input$tbl.tab1_rows_selected),-1])) %>%
    #  add_col_labels(ticktext = colnames(vst.matrix[,-1])) %>%
    # add_row_labels(ticktext = vst.matrix[input$tbl.tab1_rows_selected,1])
    #}
    
    heatmap_expr  
    
  })
  
  ##now working on heatmap on tab2 for clustering cell lines based on input genes
  
  ##will probably have to put tab2 warning message here
  output$warning_message_tab2 <- renderText({
    if(length(input$tbl.tab2_rows_selected) < 2 && input$select_clus == "selected genes")
    {
      return("Selected gene heatmap will be displayed when 2 or more genes are selected")
    }
    #else
    #{
    # return(NULL) ##might not need this, this might be messing up filtered genes heatmap
    #}
  })
  
  
  output$heatmap_clus <- renderIheatmap({
    
    ##will need to do something to show filtering step (taking out genes with all-zero-rows), or is filtering step necessary????? review DESeq2 vignette
    ##for filtering steps will have to put all object manipulation in top of code
    
    if(input$select_tab2 == "raw counts") matrix_clus <- count.matrix.filtered[,-1] ##clus for cluster, treat as cmf, do -1 to take out gene names
    if(input$select_tab2 == "rlog") matrix_clus <- cmf.rlog
    if(input$select_tab2 == "vst") matrix_clus <- cmf.vst
    if(input$select_tab2 == "row normalized") matrix_clus <- cmf.rownorm
    
    #replace above command with this based on select input
    if(input$select_clus == "-no selection-") {return(NULL)} ##commenting it out still has filtered hm show automatically
    #if(is.null(input$tbl.tab2_rows_selected)) {return(NULL)} ##necessary???
    
    if(length(input$tbl.tab2_rows_selected) < 2 && input$select_clus == "selected genes") {return(NULL)} ##try this out to see how affects heatmaps
    
    ##BT549 disapears from list of cell lines???
    ##how to make this heatmap show by default/automatically
    ##doesn't work with raw counts
    if(input$select_clus == "filtered genes")
    {
      #dend.clus <- hclust(dist(t(matrix_clus))) ##try not creating it as an object
      
      heatmap_clus <- main_heatmap(as.matrix(dist(t(matrix_clus)))) %>%
        add_col_labels(ticktext = colnames(matrix_clus[,-1])) %>%
        add_row_labels(ticktext = colnames(matrix_clus[,-1])) %>%
        add_col_dendro(hclust(dist(t(matrix_clus))), reorder = TRUE) %>%
        add_row_dendro(hclust(dist(t(matrix_clus))), reorder = TRUE, side = "right")
    }
    
    ##selecting one works with raw counts but not other transforms??? is it calculating distances right???
    if(input$select_clus == "selected genes")
    {
      #if(is.null(input$tbl.tab2_rows_selected)) {return(NULL)} ##might need to take this out (but its in tiagos code???)
      #dend.clus <- hclust(dist(t(matrix_clus))) ##try not creating it as an object ##dont need the object?
      heatmap_clus <- main_heatmap(as.matrix(dist(t(matrix_clus[c(input$tbl.tab2_rows_selected),-1])))) %>% ##partially working,
        add_col_labels(ticktext = colnames(matrix_clus[,-1])) %>%
        add_row_labels(ticktext = colnames(matrix_clus[,-1])) %>% ##works when not using add dendro, but calculates dist wrong?
        add_col_dendro(hclust(dist(t(matrix_clus[c(input$tbl.tab2_rows_selected),-1]))), reorder = TRUE) %>% ##add_dendro not working...save for later, try taking out t(matrix[]), but put back in later if it doesnt work
        add_row_dendro(hclust(dist(t(matrix_clus[c(input$tbl.tab2_rows_selected),-1]))), reorder = TRUE, side = "right") ##try taking out t(matrix[]), but put back in later if it doesnt work
      
      ##line ~218 doesn't stop this from showing, it's the dendrogram fcts fault
      
      ##checking row selection is working right
      #testprint <- renderPrint({input$tbl.tab2_rows_selected})
      #testprint()
      
      
    }
    
    heatmap_clus
    
  })
  
}

shinyApp(ui = ui, server = server)


