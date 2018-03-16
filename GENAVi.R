source("aux_functions.R")$value

ui <- fluidPage(theme = shinytheme("spacelab"),
                useShinyjs(),
                tabsetPanel(
                  tabPanel("Data Expression",  ##changing from tab 1, but still using tab1 in the other parts of code
                           icon = icon("table"),
                           h1('GENAVi'), 
                           sidebarPanel(
                             h2('Table'), 
                             selectInput("select_tab1", "Select Transform", transforms, multiple = FALSE), ##need individual selectInputs for each tab
                             fileInput("input_gene_list_tab1", "Input Gene Symbol List (Optional)", multiple = FALSE, accept = NULL, width = NULL, buttonLabel = "Browse", placeholder = "No file selected"), ##how to increase max upload size
                             #textAreaInput(inputId = "geneList",label = "Gene list filter: separate gene names by , or ; or newline",value =  "", width = "100%"),
                             #actionButton("but_sortSelectedFirst_tab1", "Selected Rows First"), ##do this to put selected rows at top of data table, trying it out
                             selectInput("select_sort_tab1", "Sort Table By", sortby, multiple = FALSE),
                             tags$hr(),
                             h2('Data upload'), 
                             # Input: Select a file ----
                             fileInput("rawcounts", "Choose CSV File",
                                       multiple = TRUE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
                             tags$div(
                               HTML(paste(help_text))
                             )
                           ),
                           mainPanel(
                             DT::dataTableOutput('tbl.tab1') ##dont think i need to change this to calc/render data tables live
                           )
                  ),
                  tabPanel("Plots", ##changing from tab 2, but still usibg tab2 in other parts of code
                           icon = icon("object-group"),
                           h2("GENAVi"), 
                           
                           div(id = "expression_plots",
                               h2('Counts Barplot'), 
                               plotlyOutput("barplot", width = "auto")
                           ),
                           h2('Counts Heatmap'),
                           iheatmaprOutput("heatmap_expr",height = "auto"),
                           div(id = "cluster_plots",
                               h2('Correlation Heatmap'), 
                               selectInput("select_clus", "Cluster by what genes", cell.line.clusters, multiple = FALSE),
                               iheatmaprOutput("heatmap_clus") 
                           )
                           
                  )
                )
                
)

server <- function(input,output) 
{
  
  output$contents <-  DT::renderDataTable({
    data <- getNormalizedData()$raw
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
      withProgress(message = 'Reading the data',
                   detail = "This may take a while", value = 0, {
                     ret <-  read_csv(inFile$datapath)
                     setProgress(1, detail = paste("Completed"))
                   }
      )
    }
    ret
  })
  ### reactive fct that calcs the transforms and saves them so it doesnt take too long each time
  getNormalizedData <- reactive({
    if (!is.null(readData())) all_cell_lines <- readData()
    
    # Add gene metadata information
    withProgress(message = 'Adding gene metadata',
                 detail = "This may take a while", value = 0, {
                   all_cell_lines <- addgeneinfo(all_cell_lines) 
                 }
    )
    
    tbl.tab1 <- all_cell_lines[rowSums(all_cell_lines[,7:ncol(all_cell_lines)]) > 1,] ##filtering step, actually change the object
    data <- as.matrix(tbl.tab1[,7:ncol(tbl.tab1)])
    metadata <- tbl.tab1[,1:6]
    
    withProgress(message = 'Normalizing data',
                 detail = "This may take a while", value = 0, {
                   # normalization: rlog takes a lot of time (hours for a big matrix)
                   raw      <- cbind(metadata, data) ##might might have to take out blind option???
                   #rlog     <- cbind(metadata, rlog(data))
                   vst      <- cbind(metadata, vst(data))
                   setProgress(0.2, detail = paste("vst completed"))
                   rownorm  <- cbind(metadata, rownorm(data))
                   setProgress(0.5, detail = paste("rownorm completed"))
                   cpm      <- cbind(metadata, cpm(data))
                   setProgress(0.7, detail = paste("cpm completed"))
                   ret      <- list(vst,rownorm,raw,cpm)
                   names(ret) <- c("vst","rownorm","raw","cpm")
                   setProgress(1, detail = paste("Completed"))
                 }
    )
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
      selected_rows <- unique(c(selected_rows,which(tbl.tab1$Symbol %in% geneList)))
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
  
  
  observeEvent(input$tbl.tab1_rows_selected, {
    if(length(input$tbl.tab1_rows_selected) == 1){
      shinyjs::show(id = "expression_plots", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    } else {
      shinyjs::hide(id = "expression_plots", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    }
  })
  
  output$barplot <- renderPlotly({
    
    if(is.null(input$tbl.tab1_rows_selected)) return(NULL) ##may need to put this in heatmap section and in tab2
    if(length(input$tbl.tab1_rows_selected) > 1) return(NULL)
    tbl.tab1 <- getTab1()
    tbl.tab1 <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected)
    p <- as.data.frame(t(tbl.tab1[,7:ncol(tbl.tab1)]))
    colnames(p) <- "value"
    p$cell_line <- rownames(p)
    barplot <- ggplot(p, aes(x=cell_line, y=value)) + geom_bar(stat = "identity") +  theme_bw()
    ggplotly(barplot)
  })
  
  
  output$heatmap_expr <- renderIheatmap({ ###### heatmap is under construction too...raw counts doesnt work...need to get saving obj code to work
    
    #if(is.null(input$tbl.tab1_rows_selected)) {return(NULL)} ##necessary???
    if(length(input$tbl.tab1_rows_selected) < 2) return(NULL)
    
    tbl.tab1 <- getTab1()
    # Columns 1 to 6: Genename  Geneid Chr   Start   End Strand  
    geneNames <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% pull("Symbol")
    matrix_expr <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% select(7:ncol(tbl.tab1)) 
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
  
  
  output$heatmap_clus <- renderIheatmap({
    
    tbl.tab2 <- getTab1()
    matrix_clus <- tbl.tab2[,c(1,7:ncol(tbl.tab2))] ### trying this out
    
    #replace above command with this based on select input
    if(input$select_clus == "-no selection-") return(NULL) ##commenting it out still has filtered hm show automatically
    #if(is.null(input$tbl.tab2_rows_selected)) {return(NULL)} ##necessary???
    
    if(length(input$tbl.tab1_rows_selected) < 2 && input$select_clus == "selected genes") {return(NULL)} ##try this out to see how affects heatmaps
    
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
      selected_rows <- input$tbl.tab1_rows_selected
      inFile <- input$input_gene_list_tab1
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