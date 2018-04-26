source("aux_functions.R")$value

ui <- fluidPage(title = "GENAVi",
                theme = shinytheme("spacelab"),
                tags$head(tags$style(
                  HTML('
                       #sidebar {
                       background-color: #ffffff;
                       }
                       
                       body, label, input, button, select { 
                       font-family: "Arial";
                       }')
                )),
                titlePanel("GENAVi"),
                useShinyjs(),
                tabsetPanel( #type = "pills",
                  tabPanel("Data Expression",  ##changing from tab 1, but still using tab1 in the other parts of code
                           icon = icon("table"),
                           sidebarPanel(id="sidebar",
                                        h3('Table'), 
                                        selectInput("select_tab1", "Select Transform", transforms, multiple = FALSE), ##need individual selectInputs for each tab
                                        fileInput("input_gene_list_tab1", "Input Gene Symbol List (Optional)", multiple = FALSE, accept = NULL, width = NULL, buttonLabel = "Browse", placeholder = "No file selected"), ##how to increase max upload size
                                        #textAreaInput(inputId = "geneList",label = "Gene list filter: separate gene names by , or ; or newline",value =  "", width = "100%"),
                                        #actionButton("but_sortSelectedFirst_tab1", "Selected Rows First"), ##do this to put selected rows at top of data table, trying it out
                                        #selectInput("select_sort_tab1", "Sort Table By", sortby, multiple = FALSE),
                                        tags$hr(),
                                        h3('Data upload'), 
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
                           #icon = icon("object-group"),
                           icon = icon("image"),
                           tabsetPanel(type = "pills",
                                       tabPanel("Expression plots",
                                                icon = icon("bar-chart-o"),
                                                bsAlert("genemessage"),
                                                hidden(
                                                  div(id = "expression_plots",
                                                      h3('Expression Barplot'), 
                                                      plotlyOutput("barplot", width = "auto")
                                                  )),
                                                hidden(
                                                  div(id = "expression_heatmap",
                                                      h3('Expression Heatmap'),
                                                      iheatmaprOutput("heatmap_expr",height = "auto")
                                                  )
                                                )
                                       ),
                                       tabPanel("Clustering plots",
                                                icon = icon("object-group"),
                                                div(id = "cluster_plots",
                                                    h3('Euclidian Distance Heatmap'), 
                                                    selectInput("select_clus", "Cluster by what genes", cell.line.clusters, multiple = FALSE),
                                                    bsAlert("genemessage2"),
                                                    iheatmaprOutput("heatmap_clus",height = "800px")
                                                )
                                       )
                           )
                  ),
                  tabPanel("Differential Expression Analysis", 
                           icon = icon("flask"),
                           sidebarPanel(id="sidebar",
                                        h3('Metata upload'), 
                                        # Input: Select a file ----
                                        downloadButton('downloadData', 'Download example metadata file'),
                                        fileInput("metadata", "Choose CSV File",
                                                  multiple = TRUE,
                                                  accept = c("text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv")),
                                        tags$div(
                                          HTML(paste(help_text2))
                                        ),
                                        tags$hr(),
                                        h3('DEA - DESeq2'), 
                                        selectInput("condition", "Select condition column for DEA", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                        selectInput("covariates", 
                                                    label = "Select covariates for DEA",
                                                    choices =  NULL, 
                                                    multiple = FALSE), ##need individual selectInputs for each tab
                                        verbatimTextOutput("formulatext"),
                                        selectInput("reference", "Select reference level for DEA", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                        actionButton("dea", "Perform DEA"),
                                        tags$hr(),
                                        h3('DEA - Select Results'), 
                                        selectInput("deaSelect", "Select results", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                        checkboxInput(inputId="lfc", label = "Perform Log fold change shrinkage", value = FALSE, width = NULL),
                                        tags$hr(),
                                        h3('Volcano plot'), 
                                        numericInput("log2FoldChange", "log2FoldChange  cut-off:", 1, min = 0, max = 10, step = 0.1),
                                        numericInput("padj", "P adjusted cut-off:", 0.01, min = 0, max = 1,step = 0.1),
                                        actionButton("volcanoplotBt", "Plot volcano plot")
                                        
                           ),
                           mainPanel(
                             bsAlert("deamessage"),
                             tabsetPanel(type = "pills",
                                         id = "DEA",
                                         tabPanel("Metadata",
                                                  tags$hr(),
                                                  DT::dataTableOutput('metadata.tbl')
                                         ), 
                                         tabPanel("DEA results",
                                                  tags$hr(),
                                                  DT::dataTableOutput('dea.results') 
                                         ),
                                         tabPanel("Volcano plot",
                                                  tags$hr(),
                                                  plotlyOutput('volcanoplot') 
                                         )
                             )
                           )
                  )
                )
)

server <- function(input,output,session) 
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
  
  
  readMetaData <- reactive({
    ret <- NULL
    inFile <- input$metadata
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
    if (!is.null(readData())) {
      all_cell_lines <- readData()
      
      # Add gene metadata information
      withProgress(message = 'Adding gene metadata',
                   detail = "This may take a while", value = 0, {
                     # We will check if metadata was added
                     res <- getEndGeneInfo(all_cell_lines)
                     all_cell_lines <- res$data
                     ngene <- res$ngene
                   }
      )
      
      tbl.tab1 <- all_cell_lines[rowSums(all_cell_lines[,(ngene + 1):ncol(all_cell_lines)]) > 1,] ##filtering step, actually change the object
      data <- as.matrix(tbl.tab1[,(ngene + 1):ncol(tbl.tab1)])
      metadata <- tbl.tab1[,1:(ngene)]
      
      withProgress(message = 'Normalizing data',
                   detail = "This may take a while", value = 0, {
                     # normalization: rlog takes a lot of time (hours for a big matrix)
                     raw      <- cbind(metadata, data) ##might might have to take out blind option???
                     setProgress(0.1, detail = paste("Starting VST"))
                     vst      <- cbind(metadata, vst(data))
                     setProgress(0.2, detail = paste("VST completed, starting rownorm"))
                     rownorm  <- cbind(metadata, rownorm(data))
                     setProgress(0.5, detail = paste("rownorm completed, starting CPM"))
                     cpm      <- cbind(metadata, cpm(data))
                     setProgress(0.7, detail = paste("CPM completed, starting rlog"))
                     rlog     <- cbind(metadata, rlog(data))
                     ret      <- list(vst,rownorm,raw,cpm,rlog)
                     names(ret) <- c("vst","rownorm","raw","cpm","rlog")
                     setProgress(1, detail = paste("Completed"))
                   }
      )
    } else {
      ret <- get(load("genavi.rda"))
    }
    return(ret)
  })
  getTab1 <- reactive({
    data <- getNormalizedData()
    select <- input$select_tab1
    if(select == "raw counts") tbl.tab1 <- data$raw #table.counts #DT::datatable(table.counts)
    if(select == "rlog")  tbl.tab1 <- data$rlog ##took this out from app??
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
    
    tbl.tab1 %>% createTable(selected_rows,tableType = isolate({input$select_tab1}))
    
    ## try adding the genes list to match() here, see if it breaks the app
    ##gene_list_tab1 <- input$input_gene_list_tab1 ##create object here see if it breaks the app.....yup breaks the app
    
  }) ##works to get selected rows on top but fucks up if select more than one at a time.....worry about it later
  
  sel <- reactive({!is.null(input$tbl.tab1_rows_selected)})  
  
  observeEvent(sel(), {
    closeAlert(session, "geneAlert2")
    
    if(!sel()){
      createAlert(session, "genemessage", "geneAlert2", title = "Missing data", style =  "danger",
                  content = paste0("Please select genes in Data expression tab"),
                  append = FALSE)
      shinyjs::hide(id = "expression_heatmap", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::hide(id = "expression_plots", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    } else if(length(input$tbl.tab1_rows_selected) == 1){
      shinyjs::show(id = "expression_plots", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::hide(id = "expression_heatmap", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    } else {
      shinyjs::show(id = "expression_heatmap", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
      shinyjs::hide(id = "expression_plots", anim = FALSE, animType = "slide", time = 0.5,selector = NULL)
    }
  })
  
  output$barplot <- renderPlotly({
    
    if(is.null(input$tbl.tab1_rows_selected)) {
      createAlert(session, "genemessage", "geneAlert2", title = "Missing data", style =  "danger",
                  content = paste0("Please select genes in Data expression tab"),
                  append = FALSE)
      return(NULL) 
    }
    if(length(input$tbl.tab1_rows_selected) == 0) {
      createAlert(session, "genemessage", "geneAlert2", title = "Missing data", style =  "danger",
                  content = paste0("Please select genes in Data expression tab"),
                  append = FALSE)
      return(NULL) 
    }
    ##may need to put this in heatmap section and in tab2
    if(length(input$tbl.tab1_rows_selected) > 1) return(NULL)
    tbl.tab1 <- getTab1()
    
    # Columns 1 to 7: Genename  Geneid Chr   Start   End Strand Length 
    res <- getEndGeneInfo(tbl.tab1)
    ngene <- res$ngene

    tbl.tab1 <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected)
    p <- as.data.frame(t(tbl.tab1[,(res$ngene+1):ncol(tbl.tab1)]))
    colnames(p) <- "value"
    p$cell_line <- rownames(p)
    barplot <- ggplot(p, aes(x=cell_line, y=value)) + 
      geom_bar(stat = "identity") +  
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
      
    ggplotly(barplot)
  })
  
  #---------------------------------
  # Heatmap plot tab
  #---------------------------------
  
  output$heatmap_expr <- renderIheatmap({ 
    
    # we can't do a heatmap with only one gene
    if(length(input$tbl.tab1_rows_selected) < 2) return(NULL)
    
    tbl.tab1 <- getTab1()
    
    
    if("Symbol" %in% colnames(tbl.tab1)){
      geneNames <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% pull("Symbol")
    } else if("Genename" %in% colnames(tbl.tab1))  {
      geneNames <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% pull("Genename")
    }
    
    # Columns 1 to 7: Genename  Geneid Chr   Start   End Strand Length 
    res <- getEndGeneInfo(tbl.tab1)
    ngene <- res$ngene

    matrix_expr <- tbl.tab1 %>% slice(input$tbl.tab1_rows_selected) %>% select((res$ngene+1):ncol(tbl.tab1)) 
    ##may need to change order of cell lines from default alphabetic to histotype specific???...do that with dendro???
    heatmap_expr <- main_heatmap(as.matrix(matrix_expr)) %>%
      add_col_labels(ticktext = colnames(matrix_expr)) %>%
      add_row_labels(ticktext = geneNames) %>% ##trying to add dendro
      add_col_dendro(hclust(dist(t(as.matrix(matrix_expr)))), reorder = TRUE) ##may have to take out -1 to avoid losing 1st data col
    
    if(nrow(matrix_expr) > 1) ##currently still trying to cluster genes selected
    {
      heatmap_expr <- heatmap_expr %>% add_row_dendro(hclust(dist((matrix_expr))), reorder = TRUE, side = "right")
    } ##taking out t() works but still has to be there...see DESeq2 workflow
    print(heatmap_expr)  ## currently rlog visualization takes too long
  })
  
  
  output$heatmap_clus <- renderIheatmap({
    closeAlert(session, "geneAlert")
    tbl.tab2 <- getTab1()
    
    # Columns 1 to 7: Genename  Geneid Chr   Start   End Strand Length 
    res <- getEndGeneInfo(tbl.tab2)
    ngene <- res$ngene

    matrix_clus <- tbl.tab2[,c(1,(res$ngene+1):ncol(tbl.tab2))] 
    
    #replace above command with this based on select input
    if(input$select_clus == "-no selection-") return(NULL) ##commenting it out still has filtered hm show automatically
    #if(is.null(input$tbl.tab2_rows_selected)) {return(NULL)} ##necessary???
    
    ##BT549 disapears from list of cell lines???
    ##how to make this heatmap show by default/automatically
    ##doesn't work with raw counts
    if(input$select_clus == "All genes")
    {
      #dend.clus <- hclust(dist(t(matrix_clus))) ##try not creating it as an object
      
      heatmap_clus <- main_heatmap(as.matrix(dist(t(matrix_clus[,-1])))) %>%
        add_col_labels(ticktext = colnames(matrix_clus[,-1])) %>%
        add_row_labels(ticktext = colnames(matrix_clus[,-1])) %>%
        add_col_dendro(hclust(dist(t(matrix_clus[,-1]))), reorder = TRUE) %>%
        add_row_dendro(hclust(dist(t(matrix_clus[,-1]))), reorder = TRUE, side = "right")
    } else { # selected genes
      selected_rows <- input$tbl.tab1_rows_selected
      if(length(selected_rows) < 1) {
        createAlert(session, "genemessage2", "geneAlert", title = "Missing data", style =  "danger",
                    content = paste0("Please select genes in Data expression tab"),
                    append = FALSE)
        return(NULL)
      }
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
  
  #------------------------------------------
  # DEA - differential expression analysis
  #------------------------------------------
  output$downloadData <- downloadHandler(
    filename = function() {
      file  <-   paste0("Genavi-metadata.csv")
    },
    content = function(con) {
      if(!is.null(readData())) all_cell_lines <- readData()
      data <- all_cell_lines
      file  <-  paste0("Genavi-metadata.csv")
      samples <- colnames(data)[-1]
      metadata <- data.frame(samples,groups = "Control")
      write_csv(metadata, con)
    }
  )
  
  observe({
    metadata <- readMetaData()
    if(!is.null(metadata)){
      updateSelectizeInput(session, 'condition', choices =  colnames(metadata)[-1], server = TRUE)
    }
    if(!is.null(metadata)){
      updateSelectizeInput(session, 'covariates', choices =  c(" ",colnames(metadata)[-1]), server = TRUE)
    }
    #if(is.null(getDataType(as.logical(input$tcgaDatabase),input$tcgaDataCategoryFilter))) {
    #  shinyjs::hide("tcgaDataTypeFilter")
    #} else {
    #  shinyjs::show("tcgaDataTypeFilter")
    #}
  })
  
  
  get.DEA.results <- reactive({
    closeAlert(session, "deaAlert")
    input$dea
    metadata <- readMetaData()
    if(is.null(metadata)) {
      createAlert(session, "deamessage", "deaAlert", title = "Missing metadata", style =  "danger",
                  content = paste0("Please upload metadata file"),
                  append = FALSE)
      return(NULL)
    }
    if(!is.null(readData())) all_cell_lines <- readData()
    
    res <- getEndGeneInfo(all_cell_lines)
    all_cell_lines <- res$data
    ngene <- res$ngene
    
    genes <- all_cell_lines$Symbol
    cts <- as.matrix(all_cell_lines[,(ngene + 1):ncol(all_cell_lines)])
    rownames(cts) <-  genes
    
    # Read aux values required for analysis (condition, covariates and reference value)
    cond <- isolate(input$condition)
    cov <- isolate(input$covariates)
    ref <-  isolate(input$reference)
    
    if(is.null(cond))   {
      createAlert(session, "deamessage", "deaAlert", title = "Missing metadata", style =  "danger",
                  content = paste0("Please select condition file"),
                  append = FALSE)
      return(NULL)
    } 
    form <- getFormula()
    if(is.null(form)){
      createAlert(session, "deamessage", "deaAlert", title = "Missing formula", style =  "danger",
                  content = paste0("Please select condition column"),
                  append = FALSE)
      return(NULL)
    }
    if(str_length(ref) == 0)   {
      createAlert(session, "deamessage", "deaAlert", title = "Missing reference level", style =  "danger",
                  content = paste0("Please select reference level"),
                  append = FALSE)
      return(NULL)
    } 
    
    withProgress(message = 'DESeq2 Analysis',
                 detail = "Creating input file", value = 0, {
                   dds <-  tryCatch({
                     dds <- DESeqDataSetFromMatrix(countData = cts,
                                                   colData = metadata,
                                                   design = form)
                     setProgress(0.2, detail = paste("Performing DEA"))
                     keep <- rowSums(counts(dds)) >= 10
                     dds[[cond]] <- relevel(dds[[cond]], ref = ref)
                     dds <- dds[keep,]
                     dds <- DESeq(dds)
                     return(dds)
                   }, error = function(e){
                     createAlert(session, "deamessage", "deaAlert", 
                                 title = "Error in DEA", style =  "danger",
                                 content = paste0(e),
                                 append = FALSE)
                     
                     return(NULL)
                   })
                 }
    )
    return(dds)
  })  
  
  getFormula <- reactive({
    form <- NULL
    cond <- input$condition
    cov <- input$covariates
    if(str_length(cond) > 0 & (str_length(cov) == 0 | str_length(cov) == 1 & cov == " " )) {
      form <- as.formula(paste0("~ ", cond))
    } else if(str_length(cov) > 0) {
      form <- as.formula(paste0("~ ",paste(cov,collapse = "+")," + ", cond))
    }
    return(form)
  })
  
  output$formulatext <- renderText({
    f <- getFormula()
    if(!is.null(f)) return(as.character(f))
    return("")
  })
  
  output$metadata.tbl <-  DT::renderDataTable({
    metadata <- readMetaData()
    if(is.null(metadata)) {
      return(NULL)
    }
    metadata  %>% createTable2(show.rownames=F)
  })
  
  observeEvent(input$condition, {
    metadata <- readMetaData()
    if(!is.null(metadata)) {
      updateSelectizeInput(session, 'reference', choices =  as.character(unique(metadata %>% pull(input$condition))), server = TRUE)
    }
  })
  
  observeEvent(input$dea, {
    updateTabsetPanel(session, inputId="DEA", selected = "DEA results")
    if(!is.null(get.DEA.results())) updateSelectizeInput(session, 'deaSelect', choices =  resultsNames(get.DEA.results()), server = TRUE)
    output$dea.results <-  DT::renderDataTable({
      res <- get.DEA.results()
      if(is.null(res)) return(NULL)
      deaSelect <- input$deaSelect
      if(str_length(deaSelect) == 0) {
        tbl <-  as.data.frame(results(res))
      } else {
        if(input$lfc) {
          tbl <-  as.data.frame(lfcShrink(res, coef = deaSelect))
        } else {
          tbl <-  as.data.frame(results(res, name = deaSelect))
        }
      }
      tbl %>% createTable2(show.rownames=T)
    })
  })
  
  
  #---------------------------------
  # Volcano plot tab
  #---------------------------------
  observeEvent(input$volcanoplotBt, {
    updateTabsetPanel(session, inputId="DEA", selected = "Volcano plot")
    output$volcanoplot <- renderPlotly({
      res <- get.DEA.results()
      
      if(is.null(res)) return(NULL)
      deaSelect <- input$deaSelect
      if(str_length(deaSelect) == 0) {
        dea <-  as.data.frame(results(res))
      } else {
        if(input$lfc) {
          dea <-  as.data.frame(lfcShrink(res, coef = deaSelect))
        } else {
          dea <-  as.data.frame(results(res, name = deaSelect))
        }
      }
      x.cut <- isolate({input$log2FoldChange})
      y.cut <- isolate({input$padj})
      
      dea$group <- "Not Significant"
      dea[which(dea$padj < y.cut & dea$log2FoldChange < -x.cut ),"group"] <- "Downregulated"
      dea[which(dea$padj < y.cut & dea$log2FoldChange > x.cut ),"group"] <- "Upregulated"
      
      
      f <- list(
        family = "Courier New, monospace",
        size = 18,
        color = "#7f7f7f"
      )
      x <- list(
        title = "log2FoldChange",
        titlefont = f
      )
      y <- list(
        title = "-log10(p-value adjusted)",
        titlefont = f
      )
      
      p <- plot_ly(data = dea, 
                   x = dea$log2FoldChange, 
                   y = -log10(dea$padj), 
                   text = rownames(dea), 
                   mode = "markers", 
                   color = dea$group) %>% 
        layout(title ="Volcano Plot") %>%
        layout(xaxis = x, yaxis = y)  %>%
        layout(shapes=list(list(type='line', 
                                x0 = x.cut, 
                                x1 = x.cut, 
                                y0 = 0, 
                                y1 = max(-log10(dea$padj),na.rm = T), 
                                line=list(dash='dot', width=1)),
                           list(type='line', 
                                x0 = -x.cut, 
                                x1 = -x.cut, 
                                y0 = 0, 
                                y1 = max(-log10(dea$padj),na.rm = T), 
                                line =list(dash='dot', width=1)),
                           list(type ='line', 
                                x0 = min(dea$log2FoldChange), 
                                x1 = max(dea$log2FoldChange), 
                                y0 =  -log10(y.cut), 
                                y1 =  -log10(y.cut), 
                                line = list(dash='dot', width=1))
        ) 
        )
      return(p)
    })
  })
}

shinyApp(ui = ui, server = server)
