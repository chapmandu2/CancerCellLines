shinyRespVsGeneticServer <- function(input, output, con, drug_df=NULL) {

  #get the cell lines for the given tissues
  proc_cls <- reactive({

    if (is.null(drug_df)) {
      intersect(getTissueCellLines(con, input$tissue), proc_resp_data()$CCLE_name)
    } else {
      intersect(getTissueCellLines(con, input$tissue), proc_resp_data()$CCLE_name)
    }

  })

  #get custom and database response data into a common format for reference of what cell lines/drugs are present
  proc_resp_data <- reactive({

    if (is.null(drug_df)) {
      #if no drug data frame provided just use CCLE
      ccle_drugs.df <- src_sqlite(con@dbname) %>% tbl("ccle_drug_data") %>% select(Compound, CCLE_name) %>% distinct %>% collect
      resp_data <- getDrugData_CCLE(con, unique(ccle_drugs.df$Compound), unique(ccle_drugs.df$CCLE_name)) %>%
        transmute(CCLE_name, resp_id=ID, resp_value=value)
      return(resp_data)
    } else {
      #if drug data provided then use that
      resp_data <- getDrugData_custom(drug_df, unique(drug_df$compound_id), unique(drug_df$unified_id) ) %>%
        transmute(CCLE_name, resp_id=ID, resp_value=value)
      return(resp_data)
    }
  })


  proc_data <- reactive({
    makeRespVsGeneticDataFrame(con, gene=input$geneid,
                               cell_lines=proc_cls(),
                               drug=input$respid,
                               data_types = input$data_type,
                               drug_df = drug_df)
  })

  #make a reactive tissue selection UI
  output$tissuesUI <- renderUI({

    tissues.df <- src_sqlite(con@dbname) %>%
      tbl('ccle_sampleinfo') %>%
      filter(CCLE_name %in% proc_resp_data()$CCLE_name) %>%
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

  #make a reactive response variable selection UI
  output$respUI <- renderUI({

    resp.df <- proc_resp_data() %>% select(resp_id) %>%
      distinct() %>% arrange(resp_id)

    resps <- resp.df$resp_id

    selectInput("respid", label = h3("Select a response variable"),
                choices = resps,
                selected = resps[1])

  })

  output$resultsUI <- renderUI({
    mainPanel(#textOutput('text1'),
      if (input$output_option == 1) {
        plotOutput("plot1")
      } else if (input$output_option == 2) {
        plotOutput("plot2")
      } else {
        tableOutput('df')
      }

    )
  })

  output$text1 <- renderText({
    sprintf("Tissue is length: %s", length(input$tissue))
  })

  output$df <- renderTable({
    proc_data()
    #get_shiny_cell_lines(con, input$tissue) %>% as.data.frame
  })

  output$plot1 <- renderPlot({

    plotRespVsGeneticPoint(proc_data(), data_type=input$data_type, facet_option=input$facet_option)

  })

  output$plot2 <- renderPlot({

    plotRespVsGeneticHist(proc_data(), input$data_type, input$facet_option)

  })

}
