library(shiny)
library(tidyr) 
library(ggplot2) 
library(DESeq2)
library(edgeR)
library(iheatmapr)
library(tidyverse)
library(readr)
library(dplyr)
library(plotly)

options(shiny.maxRequestSize=1024^3) # Max file upload 1GB 
all_cell_lines <- read_csv("all_cell_lines.csv", col_names = TRUE) ##read in data going to view

transforms <- c("raw counts", "row normalized", "logCPM", "vst") ##adding "t-score" -> row normalized, ask michelle about normalizing with FTSEC lines, have to add log(CPM+1)
sortby <- c("-no selection-","mean", "standard deviation")
cell.line.clusters <- c("-no selection-","filtered genes", "selected genes") ##do this for cell line cluster heatmaps, changed all to filtered, should "all" be an option?

# This will be used to parse the text areas input
# possibilities of separation , ; \n
parse.textarea.input <- function(text){
  sep <- NULL
  if(grepl(";",text)) sep <- ";"
  if(grepl(",",text)) sep <- ","
  if(grepl("\n",text)) sep <- "\n"
  if(is.null(sep)) {
    text <- text
  } else {
    text <- unlist(stringr::str_split(text,sep))
  }
  return (text)
}

createTable <- function(df,selected_rows,tableType = "TCGAbiolinks"){
  DT::datatable(df,
                extensions = c('Buttons',"FixedHeader"),
                class = 'cell-border stripe',
                selection = list(mode = "multiple", target= 'row', selected = selected_rows),
                options = list(dom = 'Blfrtip',
                               columnDefs = list(
                                 list(visible=FALSE, targets=c(0))
                               ),
                               order = c(0,"desc"),
                               deferRender = TRUE,
                               paging = T,
                               buttons =
                                 list('colvis', list(
                                   extend = 'collection',
                                   buttons = list(list(extend='csv',
                                                       filename = tableType),
                                                  list(extend='excel',
                                                       filename = tableType),
                                                  list(extend='pdf',
                                                       title = "",
                                                       filename= tableType)),
                                   text = 'Download'
                                 )),
                               fixedHeader = TRUE,
                               pageLength = 20,
                               scrollX = TRUE,
                               lengthMenu = list(c(10, 20, -1), c('10', '20', 'All'))
                ),
                rownames = FALSE,
                filter   = 'top'
  )
}

rownorm <- function(counts.filtered)
{
  rownorm.tbl <-  (counts.filtered - rowMeans(counts.filtered,na.rm = TRUE)) / apply(counts.filtered,1,sd)
  colnames(rownorm.tbl) <- colnames(counts.filtered)
  rownorm.tbl
}

ui <- fluidPage(tabsetPanel(
  
  
  tabPanel("Expression", ##changing from tab 1, but still using tab1 in the other parts of code
           h2('GENAVi'), 
           #selectInput("select", "Select Columns", names(final.counts), multiple = TRUE),
           selectInput("select_tab1", "Select Transform", transforms, multiple = FALSE), ##need individual selectInputs for each tab
           fileInput("input_gene_list_tab1", "Input Gene List (Optional)", multiple = FALSE, accept = NULL, width = NULL, buttonLabel = "Browse", placeholder = "No file selected"), ##how to increase max upload size
           #textAreaInput(inputId = "geneList",label = "Gene list filter: separate gene names by , or ; or newline",value =  "", width = "100%"),
           #actionButton("but_sortSelectedFirst_tab1", "Selected Rows First"), ##do this to put selected rows at top of data table, trying it out
           selectInput("select_sort_tab1", "Sort Table By", sortby, multiple = FALSE),
           DT::dataTableOutput('tbl.tab1'), ##dont think i need to change this to calc/render data tables live
           verbatimTextOutput("warning_message_tab1", placeholder = FALSE), ##NEED TO ADD ERROR MESSAGE FOR HEATMAP
           #verbatimTextOutput("warning_message_tab1_hm",placeholder = FALSE), ## hm for heatmap
           #verbatimTextOutput("test_message",placeholder = FALSE), ###testing to see if reactive fct is being called
           #numericInput("num","Input row index",1,min = 1,max = 60554),
           #verbatimTextOutput("warning_message"), ##may need to take this out for  now
           #DT::dataTableOutput('tbl.tab1.selected'), ##doing this to keep track of selected rows, take this out
           plotlyOutput("barplot", width = "auto"),
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
  ),
  # App title ----
  tabPanel("Upload file",
           # App title ----
           titlePanel("Uploading Files"),
           # Sidebar layout with input and output definitions ----
           sidebarLayout(
             
             # Sidebar panel for inputs ----
             sidebarPanel(
               
               # Input: Select a file ----
               fileInput("rawcounts", "Choose CSV File",
                         multiple = TRUE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv"))
             ),
             # Main panel for displaying outputs ----
             mainPanel(
               # Output: Data file ----
               DT::dataTableOutput("contents")
             )
           )
  )
)
)


server <- function(input,output) 
{
  output$contents <-  DT::renderDataTable({
    data <- readData()
    if(!is.null(data)) data %>% DT::datatable(
      class = 'cell-border stripe',
      rownames = FALSE,
      filter   = 'top'
    )
  })
  
  readData <- reactive({
    ret <- NULL
    inFile <- input$rawcounts
    if (!is.null(inFile))  {
      ret <-  read_csv(inFile$datapath)
    }
    ret
  })
  ### reactive fct that calcs the transforms and saves them so it doesnt take too long each time
  getNormalizedData <- reactive({
    if (!is.null(readData())) all_cell_lines <- readData()
    tbl.tab1 <- all_cell_lines[rowSums(all_cell_lines[,8:ncol(all_cell_lines)]) > 1,] ##filtering step, actually change the object
    data <- as.matrix(tbl.tab1[,8:ncol(tbl.tab1)])
    metadata <- tbl.tab1[,1:7]
    
    # normalization: rlog takes a lot of time (hours for a big matrix)
    raw      <- cbind(metadata, data) ##might might have to take out blind option???
    #rlog     <- cbind(metadata, rlog(data))
    vst      <- cbind(metadata, vst(data))
    rownorm  <- cbind(metadata, rownorm(data))
    cpm      <- cbind(metadata, cpm(data))
    ret      <- list(vst,rownorm,raw,cpm)
    names(ret) <- c("vst","rownorm","raw","cpm")
    return(ret)
  })
  getTab1 <- reactive({
    data <- getNormalizedData()
    select <- input$select_tab1
    if(select == "raw counts") tbl.tab1 <- data$raw #table.counts #DT::datatable(table.counts)
    if(select == "rlog")  tbl.tab1 <- data$rlog
    if(select == "vst")  tbl.tab1 <- data$vst
    if(select == "row normalized")  tbl.tab1 <- data$rownorm
    if(select == "logCPM")  tbl.tab1 <- data$cpm
    tbl.tab1
  })
  getTab2 <- reactive({
    data <- getNormalizedData()
    select <- input$select_tab2
    if(select == "raw counts") tbl.tab1 <- data$raw #table.counts #DT::datatable(table.counts)
    if(select == "rlog")  tbl.tab1 <- data$rlog
    if(select == "vst")  tbl.tab1 <- data$vst
    if(select == "row normalized")  tbl.tab1 <- data$rownorm
    if(select == "logCPM")  tbl.tab1 <- data$cpm
    tbl.tab1
  })
  output$tbl.tab1 <-  DT::renderDataTable({
    tbl.tab1 <- getTab1()
    
    ######### sorting by mean and sd ##################### ....fucks up the select sorting thing...
    #if(input$select_sort_tab1 == "-no selection-") {return(tbl.tab1)}
    #if(input$select_sort_tab1 == "mean")
    #{
    # tbl.tab1 <- tbl.tab1[order(apply(tbl.tab1[,8:dim(tbl.tab1)[2]],1,mean), decreasing = TRUE),]
    #}
    #if(input$select_sort_tab1 == "standard deviation")
    #{
    # tbl.tab1 <- tbl.tab1[order(apply(tbl.tab1[,8:dim(tbl.tab1)[2]],1,sd), decreasing = TRUE),]
    #}
    
    ######## this section sorts the table so that selected rows are first #######
    ####### ordering rows like this makes the selection wonky in the figures, other rows than what you select are being displayed
    ####### try taking this out to see how tables are rendered...or not rendered???
    selected_rows <- input$tbl.tab1_rows_selected
    
    
    inFile <- input$input_gene_list_tab1
    if (!is.null(inFile)) {
      geneList <- read_lines(inFile$datapath)
      selected_rows <- unique(c(selected_rows,which(tbl.tab1[,1] %in% geneList)))
    }    
    #  Parse textarea
    #text.samples <- isolate({input$geneList})
    #if(!is.null(text.samples)){
    #  geneList <- parse.textarea.input(text.samples)
    #  selected_rows <- unique(c(selected_rows,which(tbl.tab1[,1] %in% geneList)))
    #}
    
    status <- factor("Unselected",levels = c("Unselected","Selected"))
    tbl.tab1 <- cbind(status,tbl.tab1)
    tbl.tab1$status[selected_rows] <- 'Selected'
    
    tbl.tab1 %>% createTable(selected_rows)
    
    ## try adding the genes list to match() here, see if it breaks the app
    ##gene_list_tab1 <- input$input_gene_list_tab1 ##create object here see if it breaks the app.....yup breaks the app
    
  }) ##works to get selected rows on top but fucks up if select more than one at a time.....worry about it later
  
  output$tbl.tab2 <-  DT::renderDataTable({ ###need to update this to match for tab1 whatever we discover works
    tbl.tab2 <- getTab2()
    
    selected_rows <- input$tbl.tab2_rows_selected ###may need to label this for tab1 and tab2...nope don't need to
    
    inFile <- input$input_gene_list_tab2
    if (!is.null(inFile)) {
      geneList <- read_lines(inFile$datapath)
      selected_rows <- unique(c(selected_rows,which(tbl.tab2[,1] %in% geneList)))
    }    
    
    status <- factor("Unselected",levels = c("Unselected","Selected"))
    tbl.tab2 <- cbind(status,tbl.tab2)
    tbl.tab2$status[selected_rows] <- 'Selected'
    
    tbl.tab2 %>% createTable(selected_rows)
  })
  
  output$warning_message_tab1 <- renderText({
    if(length(input$tbl.tab1_rows_selected) > 1) return("Barplot is only displayed when exactly one gene is selected")
    return(NULL)
  })
  
  
  output$barplot <- renderPlotly({
    
    if(is.null(input$tbl.tab1_rows_selected)) return(NULL) ##may need to put this in heatmap section and in tab2
    if(length(input$tbl.tab1_rows_selected) > 1) return(NULL)
    tbl.tab1 <- getTab1()
    tbl.tab1 <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected)
    p <- as.data.frame(t(tbl.tab1[,8:ncol(tbl.tab1)]))
    colnames(p) <- "value"
    p$cell_line <- rownames(p)
    barplot <- ggplot(p, aes(x=cell_line, y=value)) + geom_bar(stat = "identity")
    ggplotly(barplot)
  })
  
  #output$warning_message_tab1_hm <- renderText({ ##see how this is displayed
  # if(length(input$tbl.tab1_rows_selected) < 2)
  #{
  # return("Heatmap is only displayed when two or more genes are selected")
  #}
  #else
  #{
  # return(NULL)
  #}
  #})
  
  output$heatmap_expr <- renderIheatmap({ ###### heatmap is under construction too...raw counts doesnt work...need to get saving obj code to work
    
    #if(is.null(input$tbl.tab1_rows_selected)) {return(NULL)} ##necessary???
    if(length(input$tbl.tab1_rows_selected) < 2) return(NULL)
    
    tbl.tab1 <- getTab1()
    # Columns 1 to 7: Genename  Geneid Chr   Start   End Strand Length 
    geneNames <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% pull(1)
    matrix_expr <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% select(8:ncol(tbl.tab1)) ##raw counts heat map not working...has to do with rows???
    
    ##may need to change order of cell lines from default alphabetic to histotype specific???...do that with dendro???
    heatmap_expr <- main_heatmap(as.matrix(matrix_expr)) %>%
      add_col_labels(ticktext = colnames(matrix_expr)) %>%
      add_row_labels(ticktext = geneNames) %>% ##trying to add dendro
      add_col_dendro(hclust(dist(t(as.matrix(matrix_expr))))) ##may have to take out -1 to avoid losing 1st data col
    
    if(nrow(matrix_expr) > 1) ##currently still trying to cluster genes selected
    {
      heatmap_expr <- heatmap_expr %>% add_row_dendro(hclust(dist((matrix_expr))), reorder = TRUE, side = "right")
    } ##taking out t() works but still has to be there...see DESeq2 workflow
    print(heatmap_expr)  ## currently rlog visualization takes too long
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
    
    tbl.tab2 <- getTab2()
    matrix_clus <- tbl.tab2[,c(1,8:ncol(tbl.tab2))] ### trying this out
    
    #replace above command with this based on select input
    if(input$select_clus == "-no selection-") return(NULL) ##commenting it out still has filtered hm show automatically
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
        add_col_dendro(hclust(dist(t(matrix_clus[,-1]))), reorder = TRUE) %>%
        add_row_dendro(hclust(dist(t(matrix_clus[,-1]))), reorder = TRUE, side = "right")
    }
    
    ##selecting one works with raw counts but not other transforms??? is it calculating distances right???
    if(input$select_clus == "selected genes")
    {
      selected_rows <- input$tbl.tab2_rows_selected
      inFile <- input$input_gene_list_tab2
      if (!is.null(inFile)) {
        geneList <- read_lines(inFile$datapath)
        selected_rows <- unique(c(selected_rows,which(matrix_clus[,1] %in% geneList)))
      }    
      
      #if(is.null(input$tbl.tab2_rows_selected)) {return(NULL)} ##might need to take this out (but its in tiagos code???)
      #dend.clus <- hclust(dist(t(matrix_clus))) ##try not creating it as an object ##dont need the object?
      heatmap_clus <- main_heatmap(as.matrix(dist(t(matrix_clus[selected_rows,-1])))) %>% ##partially working,
        add_col_labels(ticktext = colnames(matrix_clus[,-1])) %>%
        add_row_labels(ticktext = colnames(matrix_clus[,-1])) %>% ##works when not using add dendro, but calculates dist wrong?
        add_col_dendro(hclust(dist(t(matrix_clus[selected_rows,-1]))), reorder = TRUE) %>% ##add_dendro not working...save for later, try taking out t(matrix[]), but put back in later if it doesnt work
        add_row_dendro(hclust(dist(t(matrix_clus[selected_rows,-1]))), reorder = TRUE, side = "right") ##try taking out t(matrix[]), but put back in later if it doesnt work
    }
    heatmap_clus
  })
  
}
shinyApp(ui = ui, server = server)