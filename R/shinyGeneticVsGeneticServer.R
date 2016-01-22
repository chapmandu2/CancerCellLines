#' shinyGeneticVsGeneticServer
#'
#' Create a shiny server for the GeneticsvsGenetic shiny app
#'
#' @param input shiny input
#' @param output shiny output
#' @param con SQLite connection object
#'
#' @return a shiny server
shinyGeneticVsGeneticServer <- function(input, output, con) {

  #get the cell lines for the given tissues
  proc_cls <- reactive({
    getTissueCellLines(con, input$tissue)
  })


  proc_data <- reactive({
    makeGeneticVsGeneticDataFrame(con,
                                 cell_lines=proc_cls(),
                                 gene1 = input$gene_name1,
                                 gene2 = input$gene_name2,
                                 data_type1 = input$data_type1,
                                 data_type2 = input$data_type2)
  })

  #make a reactive tissue selection UI
  output$tissuesUI <- renderUI({

    tissues.df <- src_sqlite(con@dbname) %>%
      tbl('ccle_sampleinfo') %>%
      transmute(tissue=Site_primary) %>%
      collect %>%
      distinct() %>%
      arrange(tissue)


    tissues <- tissues.df$tissue

    selectInput("tissue", label = h3("Select a tissue type"),
                choices = tissues,
                selected = tissues,
                multiple = TRUE,
                selectize = FALSE,
                size=10)

  })


  output$resultsUI <- renderUI({

      if (input$output_option == 1) {
        mainPanel(plotOutput("plot1", width=input$plot_width, height=input$plot_height),
                  downloadButton('downloadData', 'Download Data'))
      } else if (input$output_option == 2) {
        mainPanel(plotOutput("plot2", width=input$plot_width, height=input$plot_height),
                  downloadButton('downloadData', 'Download Data'))
      } else {
        mainPanel(DT::dataTableOutput('df'),
                  downloadButton('downloadData', 'Download Data'))
      }


  })

  output$text1 <- renderText({
    sprintf("Tissue is length: %s", length(input$tissue))
  })

  output$df <- DT::renderDataTable({
    proc_data()
    #get_shiny_cell_lines(con, input$tissue) %>% as.data.frame
  }, filter='top')

  output$plot1 <- renderPlot({

    plotGeneticVsGeneticPoint(proc_data(),  facet_option=input$facet_option)

  })

  output$plot2 <- renderPlot({

    plotGeneticVsGeneticHist(proc_data(),  facet_option=input$facet_option, label_option=input$label_option)

  })

  output$downloadData <- downloadHandler(
    filename = sprintf("%s_%s_vs_%s_%s_data.txt", input$gene_name1, input$data_type1, input$gene_name2, input$data_type2),
    content = function(file) {
      write.table (proc_data(), file = file, sep = '\t', row.names = FALSE, col.names=TRUE, na='')
    }
  )

}
