#' shinyRespVsRespServer
#'
#' Create a shiny server for the RespVsResp shiny app
#'
#' @param input shiny input
#' @param output shiny output
#' @param con SQLite connection object
#' @param drug_df drug data frame
#'
#' @return a shiny server
shinyRespVsRespServer <- function(input, output, con, drug_df=NULL) {

  #get the cell lines for the given tissues
  proc_cls <- reactive({

    if (is.null(drug_df)) {
      intersect(getTissueCellLines(con, input$tissue, input$tissue_option), proc_resp_data()$unified_id)
    } else {
      intersect(getTissueCellLines(con, input$tissue, input$tissue_option), proc_resp_data()$unified_id)
    }

  })

  #get custom and database response data into a common format for reference of what cell lines/drugs are present
  proc_resp_data <- reactive({

    if (is.null(drug_df)) {
      #if no drug data frame provided just use CCLE
      ccle_drugs.df <- src_sqlite(con@dbname) %>% tbl("ccle_drug_data") %>% select(Compound, CCLE_name) %>% distinct %>% collect
      resp_data <- getDrugData_CCLE(con, unique(ccle_drugs.df$Compound), unique(ccle_drugs.df$CCLE_name)) %>%
        transmute(unified_id, resp_id=assayed_id, resp_value=value)
      return(resp_data)
    } else {
      #if drug data provided then use that
      resp_data <- getDrugData_custom(drug_df, unique(drug_df$compound_id), unique(drug_df$unified_id) ) %>%
        transmute(unified_id, resp_id=assayed_id, resp_value=value)
      return(resp_data)
    }
  })


  proc_data <- reactive({

    makeRespVsRespDataFrame(con,
                               cell_lines=proc_cls(),
                               drugs=input$respid,
                               resp_type=ifelse(is.null(drug_df),'ccle', 'custom'),
                               tissue_info = input$tissue_option,
                               drug_df = drug_df)
  })

  proc_data_wide <- reactive({

    proc_data() %>% makeWideFromRespVsRespDataFrame()

  })

  #make a reactive tissue selection UI
  output$tissuesUI <- renderUI({

    tissues.df <- getTissueInfo(con, input$tissue_option)
    resp_cls <- unique(proc_resp_data()$unified_id)

    tissues.df <- tissues.df %>%
      dplyr::filter(unified_id %in% resp_cls) %>%
      dplyr::select(tissue) %>%
      distinct %>% arrange(tissue)

    tissues <- tissues.df$tissue

    selectInput("tissue", label = h3("Select a tissue type"),
                choices = tissues,
                selected = tissues,
                multiple = TRUE,
                selectize = FALSE,
                size=10)

  })

  #make a reactive response variable selection UI
  output$respUI <- renderUI({

    resp.df <- proc_resp_data() %>% select(resp_id) %>%
      distinct() %>% arrange(resp_id)

    resps <- resp.df$resp_id

    selectInput("respid", label = h3("Select response variables"),
                choices = resps,
                selected = resps[1:2],
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
    } else if (input$output_option == 3) {
      mainPanel(plotOutput("plot3", width=input$plot_width, height=input$plot_height),
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
    proc_data_wide()
  }, filter='top')

  output$plot1 <- renderPlot({

    plotRespVsRespPairs(proc_data())

  })

  output$plot2 <- renderPlot({

    plotRespVsRespDensity(proc_data())

  })

  output$plot3 <- renderPlot({

    plotRespVsRespWaterfall(proc_data())

  })

  output$downloadData <- downloadHandler(
    filename = sprintf("%s_data.txt",  paste(sample(c(0:9, LETTERS, letters), 10, replace=TRUE), collapse='')),
    content = function(file) {
      write.table (proc_data_wide(), file = file, sep = '\t', row.names = FALSE, col.names=TRUE, na='')
    }
  )


}
