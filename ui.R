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
                  tabPanel("Gene Expression",  ##changing from tab 1, but still using tab1 in the other parts of code
                           icon = icon("table"),
                           sidebarPanel(id="sidebar",
                                        h3('Table'), 
                                        selectInput("select_tab1", "Select Transform", transforms, multiple = FALSE), ##need individual selectInputs for each tab
                                        fileInput("input_gene_list_tab1", "Input Gene Symbol List (Optional)", multiple = FALSE, accept = NULL, width = NULL, buttonLabel = "Browse", placeholder = "No file selected"), ##how to increase max upload size
                                        #textAreaInput(inputId = "geneList",label = "Gene list filter: separate gene names by , or ; or newline",value =  "", width = "100%"),
                                        #actionButton("but_sortSelectedFirst_tab1", "Selected Rows First"), ##do this to put selected rows at top of data table, trying it out
                                        #selectInput("select_sort_tab1", "Sort Table By", sortby, multiple = FALSE),
                                        downloadButton("downloadNormalizedData", "Download normalized files"),
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
                             bsAlert("tab1message"),
                             DT::dataTableOutput('tbl.tab1') ##dont think i need to change this to calc/render data tables live
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
                                                    h3('Correlation Heatmap'), 
                                                    selectInput("select_clus_type", 
                                                                label = "Cluster correlation", 
                                                                choices = c("Sample","Genes"), 
                                                                multiple = FALSE),
                                                    selectInput("select_clus", "Cluster by what genes", c("All genes", "Selected genes"), multiple = FALSE),
                                                    bsAlert("genemessage2"),
                                                    iheatmaprOutput("heatmap_clus",height = "800px")
                                                )
                                       )
                           )
                  ),
                  tabPanel("Differential Expression Analysis", 
                           icon = icon("flask"),
                           sidebarPanel(id="sidebar",
                                        h3('Metadata upload'), 
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
                                                    multiple = TRUE), ##need individual selectInputs for each tab
                                        verbatimTextOutput("formulatext"),
                                        selectInput("reference", "Select reference level for DEA", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                        actionButton("dea", "Perform DEA"),
                                        tags$hr(),
                                        h3('DEA - Select Results'), 
                                        selectInput("deaSelect", "Select results", NULL, multiple = FALSE), ##need individual selectInputs for each tab
                                        checkboxInput(inputId="lfc", label = "Perform Log fold change shrinkage", value = FALSE, width = NULL),
                                        downloadButton("downloadDEAFiles", "Download DEA Results"),
                                        tags$hr(),
                                        h3('Volcano plot'), 
                                        numericInput("log2FoldChange", "log2FoldChange  cut-off:", value = 0, min = 0, max = 10, step = 0.1),
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
                  ),
                  tabPanel("Vignette", 
                           icon = icon("book"),
                           includeMarkdown("Genavi.Rmd")
                  ),
                  tabPanel("Tutorial",
                           icon = icon("book"),
                           includeHTML("GENAVi_Tutorial.html"))
                )
)

