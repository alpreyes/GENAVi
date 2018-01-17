library(shiny)
library(DT)

shinyApp(
  ui = fluidPage(
    fluidRow(
      column(2, actionButton("but_sortSelectedFirst", "Selected Rows First")),
      column(10, DT::dataTableOutput('foo'))
    )
  ),
  server = function(input, output, session) {
    data <- reactiveValues(table_data = iris)
    
    observeEvent(input$but_sortSelectedFirst, {
      # the selected rows are in the variable input$tableid_rows_selected
      selected_rows <- input$foo_rows_selected
      
      # calculate new row order
      row_order <- order(
        seq_along(cars[[1]]) %in% selected_rows,
        decreasing = TRUE
      )
      data$table_data <- data$table_data[row_order, ]
      
      proxy <- DT::dataTableProxy('foo')
      DT::replaceData(proxy, data$table_data)
      # make sure to select the rows again
      DT::selectRows(proxy, seq_along(selected_rows))
    })
    
    output$foo <- DT::renderDataTable(isolate(data$table_data))
  }
)
