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
                       }
                       
                       .btn-file {
                        background-color:#5B81AE; 
                        border-color: #5B81AE; 
                        background: #5B81AE;
                       }
                      
         
                      .bttn-bordered.bttn-sm {
                          width: 200px;
                          text-align: left;
                          margin-bottom : 0px;
                          margin-top : 20px;
                       }
                       '
                  )
                )),
                titlePanel("GENAVi"),
                useShinyjs(),
                tabsetPanel( #type = "pills",
                  tabPanel("Gene Expression",  ##changing from tab 1, but still using tab1 in the other parts of code
                           icon = icon("table"),
                           column(2,
                                  #sidebarPanel(id="sidebar",width = 3,
                                  dropdown(label = "Data upload",
                                           icon = icon("file-excel"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           fileInput("rawcounts", "Choose CSV File",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           tags$div(
                                             HTML(paste(help_text))
                                           )
                                  ), 
                                  dropdown(label = "Table selection",
                                           icon = icon("table"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           selectInput("select_tab1", "Select Transform", transforms, multiple = FALSE)
                                  ), ##need individual selectInputs for each tab
                                  dropdown(label = "Download data",
                                           icon = icon("download"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton("downloadNormalizedData", "Download normalized files",class = "btn-primary")
                                  ),
                                  dropdown(label = "Generate report",
                                           icon = icon("file-code"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton("reportNorm", "Download report",class = "btn-primary")
                                  ),
                                
                                  dropdown(label = "Gene selection",
                                           icon = icon("mouse-pointer"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           fileInput("input_gene_list_tab1", "Input Gene Symbol List (Optional)", 
                                                     multiple = FALSE, 
                                                     accept = NULL, 
                                                     width = NULL, 
                                                     buttonLabel = "Browse", 
                                                     placeholder = "No file selected"), ##how to increase max upload size
                                           textAreaInput(inputId = "input_gene_list_area",
                                                         label = "Gene list filter: separate gene names by , or ; or newline",
                                                         value =  "", 
                                                         width = "100%"),
                                           actionButton("input_gene_list_but", 
                                                        "Select Rows",
                                                        width = "100%", 
                                                        class = "btn-primary"), ##do this to put selected rows at top of data table, trying it out
                                           actionButton("select_most_variable", 
                                                        "Select 1000 genes of highest variance",
                                                        width = "100%", 
                                                        class = "btn-primary"), ##do this to put selected rows at top of data table, trying it out
                                           actionButton("unselect_all", 
                                                        "Deselect all genes",
                                                        width = "100%", 
                                                        class = "btn-primary")
                                  ) ##do this to put selected rows at top of data table, trying it out
                                  #selectInput("select_sort_tab1", "Sort Table By", sortby, multiple = FALSE),
                           ),
                           column(10,
                                  DT::dataTableOutput('tbl.tab1') 
                           )
                  ),
                  tabPanel("Visualization", ##changing from tab 2, but still usibg tab2 in other parts of code
                           #icon = icon("object-group"),
                           icon = icon("image"),
                           tabsetPanel(type = "pills",
                                       tabPanel("Expression plots",
                                                icon = icon("bar-chart-o"),
                                                bsAlert("genemessage"),
                                                hidden(
                                                  div(id = "expression_plots",
                                                      h3('Expression Barplot'), 
                                                      plotlyOutput("barplot", width = "auto") %>% withSpinner(type = 6)
                                                  )),
                                                hidden(
                                                  div(id = "expression_heatmap",
                                                      h3('Expression Heatmap'),
                                                      #selectInput("select_z_score", 
                                                      #            label = "Standardized scores?", 
                                                      #            choices = c("No","Rows z-score", "Columns z-score"), 
                                                      #            multiple = FALSE),
                                                      iheatmaprOutput("heatmap_expr",height = "auto") %>% withSpinner(type = 6)
                                                  )
                                                )
                                       ),
                                       tabPanel("Clustering plots",
                                                icon = icon("object-group"),
                                                div(id = "cluster_plots",
                                                    column(2,
                                                           h3('Correlation Heatmap'), 
                                                           selectInput("select_clus_type", 
                                                                       label = "Cluster correlation", 
                                                                       choices = c("Sample","Genes"), 
                                                                       multiple = FALSE),
                                                           selectInput("select_clus", "Cluster by what genes", 
                                                                       c("All genes", "Selected genes"), 
                                                                       multiple = FALSE)
                                                    ),
                                                    column(9,
                                                           iheatmaprOutput("heatmap_clus",height = "800px") %>% withSpinner(type = 6)
                                                    )
                                                )
                                       ), tabPanel("PCA plots",
                                                   icon = icon("object-group"),
                                                   div(id = "pca_plots",
                                                       bsAlert("genemessage3"),
                                                       column(2,
                                                              selectInput("select_pca_type", 
                                                                          label = "PCA genes", 
                                                                          choices = c("Top 1000 variable genes", "All genes", "Selected genes"), 
                                                                          multiple = FALSE),
                                                              selectInput("pca_dimensions", 
                                                                          label = "Number of dimensions", 
                                                                          choices = c("2D", "3D"), 
                                                                          multiple = FALSE),
                                                              selectInput("pcacolor", "Color samples by", NULL, multiple = FALSE), 
                                                              downloadButton("reportPCA", "Generate report",class = "btn-primary")),
                                                       column(6,
                                                              plotlyOutput("pca_plot",height = "600",width = "600") %>% withSpinner(type = 6)
                                                       )
                                                   )
                                       )
                           )
                  ),
                  tabPanel("Differential Expression Analysis", 
                           icon = icon("flask"),
                           column(2,
                                  #sidebarPanel(id="sidebar",width = 3,
                                  dropdown(label = "Metadata upload",
                                           icon = icon("file-excel"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeft,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           # Input: Select a file ----
                                           downloadButton('downloadData', 'Download example metadata file',class = "btn-primary"),
                                           fileInput("metadata", "Choose CSV File",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           tags$div(
                                             HTML(paste(help_text2))
                                           )
                                  ),
                                  dropdown(label = "DE analysis",
                                           icon = icon("flask"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           
                                           selectInput("condition", "Select condition column for DEA", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                           selectInput("covariates", 
                                                       label = "Select covariates for DEA",
                                                       choices =  NULL, 
                                                       multiple = TRUE), ##need individual selectInputs for each tab
                                           verbatimTextOutput("formulatext"),
                                           selectInput("reference", "Select reference level for DEA", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                           actionButton("dea", "Perform DEA", class = "btn-primary")
                                  ),
                                  dropdown(label = "Select Results",
                                           icon = icon("table"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           
                                           selectInput("deaSelect", "Select results", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                           checkboxInput(inputId="lfc", label = "Perform Log fold change shrinkage", value = FALSE, width = NULL),
                                           downloadButton("downloadDEAFiles", "Download DEA Results",class = "btn-primary")
                                  ),
                                  dropdown(label = "Volcano plot",
                                           icon = icon("chart-bar"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           numericInput("log2FoldChange", "log2FoldChange  cut-off:", value = 0, min = 0, max = 10, step = 0.1),
                                           numericInput("padj", "P adjusted cut-off:", 0.01, min = 0, max = 1,step = 0.1)
                                  ),
                                  dropdown(label = "Generate report",
                                           icon = icon("file-code"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton("reportDEA", "Download Report",class = "btn-primary")
                                  )
                           ),
                           column(10,
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
                                                       plotlyOutput('volcanoplot') %>% withSpinner(type = 6)
                                              )
                                  )
                           )
                  ),
                  tabPanel("Enrichment analysis", 
                           icon = icon("flask"),
                           column(2,
                                  #sidebarPanel(id="sidebar",width = 3,
                                  dropdown(label = "DEA results upload      ",
                                           icon = icon("file-excel"),
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           size =  "sm",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           downloadButton('downloadExampleDEAData', 'Download example DEA file',class = "btn-primary"),
                                           fileInput("deafile", "Choose CSV File",
                                                     multiple = TRUE,
                                                     accept = c("text/csv",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv"))
                                  ),
                                  dropdown(label = "Enrichment Analysis",
                                           icon = icon("flask"),
                                           style = "bordered", 
                                           size =  "sm",
                                           status = "primary", 
                                           width = "300px",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           selectInput("deaanalysistype", 
                                                       "Select the type of analysis", 
                                                       c("ORA (over representation analysis)" = "ORA",
                                                         "GSEA (gene set enrichment analysis)" = "GSEA"), 
                                                       multiple = FALSE),
                                           selectInput("deaanalysisselect", 
                                                       "Select the analysis", 
                                                       c("WikiPathways analysis",
                                                         "MSigDb analysis",
                                                         "Gene Ontology Analysis",
                                                         "KEGG Analysis",
                                                         "Disease Ontology Analysis"), 
                                                       multiple = FALSE),
                                           selectInput("msigdbtype", 
                                                       "Select collection for Molecular Signatures Database", 
                                                       c("All human gene sets" = "All",
                                                         "H: hallmark gene sets" = "H",
                                                         "C1: positional gene sets" = "C1",
                                                         "C2: curated gene sets" = "C2",
                                                         "C3: motif gene sets" = "C3",
                                                         "C4: computational gene sets" = "C4",
                                                         "C5: GO gene sets" = "C5",
                                                         "C6: oncogenic signatures" = "C6",
                                                         "C7: immunologic signatures" = "C7"), 
                                                       multiple = FALSE),
                                           selectInput("gotype", 
                                                       "Select collection for Molecular Signatures Database", 
                                                       c("Molecular Function"="MF",
                                                         "Cellular Component"="CC",
                                                         "Biological Process" = "BP"), 
                                                       multiple = FALSE),
                                           numericInput("enrichmentfdr", 
                                                        "P-value cut-off:", 
                                                        value = 0.05, 
                                                        min = 0, 
                                                        max = 1, 
                                                        step = 0.05),
                                           div(id = "eaorasectui",
                                               tags$hr(),
                                               h3('ORA - selecting genes'), 
                                               numericInput("ea_subsetfdr", "P-adj cut-off", value = 0.05, min = 0, max = 1, step = 0.05),
                                               numericInput("ea_subsetlc", "LogFC cut-off", value = 1, min = 0, max = 3, step = 1),
                                               selectInput("ea_subsettype", 
                                                           "Gene status", 
                                                           c("Upregulated",
                                                             "Downregulated"), 
                                                           multiple = FALSE)
                                           ),
                                           div(id = "eagsearankingui",
                                               tags$hr(),
                                               h3('GSEA - ranking method'), 
                                               selectInput("earankingmethod", 
                                                           "Select the ranking method", 
                                                           c("log Fold Change",
                                                             "-log10(P-value) * sig(log2FC)",
                                                             "-log10(P-value) * log2FC"), 
                                                           multiple = FALSE)
                                           ),
                                           
                                           actionButton("enrichementbt", "Perform analysis", class = "btn-primary")
                                  ),
                                  dropdown(label = "Plot options",
                                           icon = icon("image"),
                                           size =  "sm",
                                           style = "bordered", 
                                           status = "primary", 
                                           width = "300px",
                                           animate = animateOptions(
                                             enter = animations$fading_entrances$fadeInLeftBig,
                                             exit = animations$fading_exits$fadeOutLeft
                                           ),
                                           
                                           selectInput("ea_plottype", 
                                                       "Plot type", 
                                                       c("Dot plot",
                                                         "Ridgeline",
                                                         "Running score and preranked list",
                                                         "Ranked list of genes"), 
                                                       multiple = FALSE),
                                           numericInput("ea_nb_categories", "Number of categories", value = 10, min = 2, max = 30, step = 1),
                                           selectInput("gsea_gene_sets", "Plot gene sets", NULL, multiple = TRUE)
                                  ),
                                  dropdown(
                                    label = "Export image",
                                    icon = icon("save"),
                                    size =  "sm",
                                    style = "bordered", 
                                    status = "primary", 
                                    width = "300px",
                                    animate = animateOptions(
                                      enter = animations$fading_entrances$fadeInLeftBig,
                                      exit = animations$fading_exits$fadeOutLeft
                                    ),
                                    tooltip = tooltipOptions(title = "Export image"),
                                    textInput("enrichementPlot.filename", label = "Filename", value = "enrichement_plot.pdf"),
                                    bsTooltip("enrichementPlot.filename", "Filename (pdf, png, svg)", "left"),
                                    numericInput("ea_width", "Figure width (in)", value = 10, min = 5, max = 30, step = 1),
                                    numericInput("ea_height", "Figure height (in)", value = 10, min = 5, max = 30, step = 1),
                                    downloadButton('saveenrichementpicture', 'Export figure',class = "btn-primary")
                                  ),
                                  dropdown(
                                    label = "Generate report",
                                    size =  "sm",
                                    icon = icon("file-code"),
                                    style = "bordered", 
                                    status = "primary", 
                                    width = "300px",
                                    animate = animateOptions(
                                      enter = animations$fading_entrances$fadeInLeftBig,
                                      exit = animations$fading_exits$fadeOutLeft
                                    ),
                                    downloadButton('reportEA', 'Download HTML report',class = "btn-primary"))
                                  ,
                                  dropdown(
                                    label = "Help material",
                                    icon = icon("book"),
                                    style = "bordered", 
                                    size =  "sm",
                                    status = "primary", 
                                    width = "300px",
                                    animate = animateOptions(
                                      enter = animations$fading_entrances$fadeInLeftBig,
                                      exit = animations$fading_exits$fadeOutLeft
                                    ),
                                    shiny::actionButton(inputId='ab1', 
                                                        label = "Learn More", 
                                                        icon = icon("th"), 
                                                        onclick = "window.open('https://guangchuangyu.github.io/pathway-analysis-workshop/', '_blank')", 
                                                        class = "btn-primary"),
                                    shiny::actionButton(inputId = 'ab1', 
                                                        label="MSigDB Collections", 
                                                        icon = icon("th"), 
                                                        onclick = "window.open('http://software.broadinstitute.org/gsea/msigdb/collection_details.jsp', '_blank')", 
                                                        class = "btn-primary"))
                           ),
                           column(8,
                                  
                                  tabsetPanel(type = "pills",
                                              tabPanel("Plots",
                                                       jqui_resizable(
                                                         plotOutput("plotenrichment", height = "600")
                                                       ) #%>% withSpinner(type = 6)
                                              ),
                                              tabPanel("Table",
                                                       DT::dataTableOutput('tbl.analysis') # %>% withSpinner(type = 6)
                                                       )
                                  )
                           )
                  ),
                  tabPanel("Help", ##changing from tab 2, but still usibg tab2 in other parts of code
                           #icon = icon("object-group"),
                           icon = icon("book"),
                           tabsetPanel(type = "pills",
                                       tabPanel("Vignette", 
                                                icon = icon("book"),
                                                includeMarkdown("Genavi.Rmd")
                                       ),
                                       tabPanel("Tutorial",
                                                icon = icon("book"),
                                                includeHTML("GENAVi_Tutorial.html")
                                       ),
                                       tabPanel("References",
                                                icon = icon("book"),
                                                includeMarkdown("References.Rmd")
                                       )
                           )
                  )
                )
)

