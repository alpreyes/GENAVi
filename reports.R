
# ---------------------------------
# Creating reports with code
# We will create for each section a HTML with the code to 
# used in Genavi
# - Required reports:
#    - DEA
#    - Normalization
#    - PCA/t-sne
#    - Volcano/Heatmap plots
# ---------------------------------
output$reportNorm <- downloadHandler(
  
  # For PDF output, change this to "report.pdf"
  filename = "report_normalization_genavi.html",
  content = function(file) {
    
    data <- readData()
    if(is.null(data)) {
      sendSweetAlert(
        session = session,
        title = "Opps...",
        text = "Missing data input",
        type = "error"
      )
      return(NULL)
    }
    # Set up parameters to pass to Rmd document
    params <- list( 
      file = isolate({input$rawcounts}),
      mouse.genes =  GRCm38,
      human.genes =  hg38,
      norm.method = isolate({input$select_tab1})
    )
    
    shiny::withProgress(value = 0,
                        message = 'Rendering Normalization report.',
                        detail =  'This might take several minutes.',{
                          rmarkdown::render(input = "report/normalization.Rmd", 
                                            params = params,
                                            output_file = file,
                                            envir = new.env(parent = globalenv()))
                        })
    
  }
)

output$reportPCA <- downloadHandler(
  
  # For PDF output, change this to "report.pdf"
  filename = "report_PCA_genavi.html",
  content = function(file) {
    
    tbl.tab1 <- getTab1()
    
    # Columns 1 to 7: Genename  Geneid Chr   Start   End Strand Length 
    res <- getEndGeneInfo(tbl.tab1)
    ngene <- res$ngene
    m <- tbl.tab1 %>% dplyr::select((res$ngene + 1):ncol(tbl.tab1)) %>% as.matrix  
    rownames(m) <- tbl.tab1 %>% pull(1)
    # Set up parameters to pass to Rmd document
    params <- list(matrix = m)
    
    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the document
    # from the code in this app).
    shiny::withProgress(value = 0,
                        message = 'Rendering PCA report.',
                        detail =  'This might take several minutes.',{
                          rmarkdown::render(input = "report/pca.Rmd", 
                                            params = params, 
                                            output_file = file,
                                            envir = new.env(parent = globalenv()))
                        })
    
  }
)

output$reportDEA <- downloadHandler(
  
  # For PDF output, change this to "report.pdf"
  filename = "report_DEA_genavi.html",
  content = function(file) {
    shiny::withProgress(value = 0,message = 'Rendering DEA report.',
                        detail =  'This might take several minutes.',{
                          
                          metadata <- readMetaData()
                          
                          if(is.null(metadata)){
                            sendSweetAlert(
                              session = session,
                              title = "Opps...",
                              text = "Missing metadata input",
                              type = "error"
                            )
                            return(NULL)
                          }
                          
                          if(!is.null(readData())) all_cell_lines <- readData()
                          res <- getEndGeneInfo(all_cell_lines)
                          matrix <- res$data
                          ngene <- res$ngene
                          
                          if("Symbol" %in% colnames(matrix)){
                            genes <- matrix %>% pull("Symbol")
                          } else if("Genename" %in% colnames(matrix))  {
                            genes <- matrix %>% pull("Genename")
                          } else {
                            genes <- matrix %>% pull(1)
                          }
                          
                          cts <- as.matrix(matrix[,(ngene + 1):ncol(matrix)])
                          rownames(cts) <-  genes
                          rownames(cts) <- matrix %>% pull(1)
                          
                          # Read aux values required for analysis (condition, covariates and reference value)
                          cond <- isolate(input$condition)
                          cov <- isolate(input$covariates)
                          ref <-  isolate(input$reference)
                          lfc <-  isolate(input$lfc)
                          log2FoldChange <- isolate(input$log2FoldChange)
                          lfc <-  isolate(input$lfc)
                          deaSelect <-  isolate(input$deaSelect)
                          padj <-  isolate(input$padj)
                          
                          # Set up parameters to pass to Rmd document
                          params <- list( 
                            log2FoldChange = log2FoldChange,
                            padj = padj,
                            metadata = metadata,
                            deaSelect = deaSelect,
                            lfc =  lfc,
                            condition = cond,
                            covariates = cov,
                            raw_cts = cts,
                            reference = ref)
                          
                          # Knit the document, passing in the `params` list, and eval it in a
                          # child of the global environment (this isolates the code in the document
                          # from the code in this app).
                          rmarkdown::render(input = "report/DEA.Rmd", 
                                            params = params, 
                                            output_file = file,
                                            envir = new.env(parent = globalenv()))
                        })
    
  }
)


output$reportEA <- downloadHandler(
  
  # For PDF output, change this to "report.pdf"
  filename = "report_EA_genavi.html",
  content = function(file) {
    shiny::withProgress(value = 0,message = 'Rendering EA report.',
                        detail =  'This might take several minutes.',{
                          
                          data <- readDEA()
                          if(is.null(data)) {
                            sendSweetAlert(
                              session = session,
                              title = "Opps...",
                              text = "Missing data input",
                              type = "error"
                            )
                            return(NULL)
                          }
                          # Set up parameters to pass to Rmd document
                          params <- list( 
                            deaanalysistype = isolate({input$deaanalysistype}),
                            ea_plottype = isolate({input$ea_plottype}),
                            gsea_gene_sets = isolate({input$gsea_gene_sets}),
                            enrichmentfdr = isolate({input$enrichmentfdr}),
                            msigdbtype = isolate({input$msigdbtype}),
                            gotype = isolate({input$gotype}),
                            deaanalysisselect = isolate({input$deaanalysisselect}),
                            ea_nb_categories = isolate({input$ea_nb_categories}),
                            earankingmethod = isolate({input$earankingmethod}),
                            ea_subsetfdr = isolate({input$ea_subsetfdr}),
                            ea_subsetlc = isolate({input$ea_subsetlc}),
                            ea_subsettype = isolate({input$ea_subsettype}),
                            data = data)
                          
                          # Knit the document, passing in the `params` list, and eval it in a
                          # child of the global environment (this isolates the code in the document
                          # from the code in this app).
                          rmarkdown::render(input = "report/ea.Rmd", 
                                            params = params, 
                                            output_file = file,
                                            envir = new.env(parent = globalenv()))
                        })
    
  }
)
