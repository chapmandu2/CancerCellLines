#functions to create a shiny app

#' Get data for shiny
#'
#' This function creates a \code{data.frame} suitable for the shiny app
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param drugs A vector of compound identifiers
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{data.frame} containing the drug response data for the requested compounds and cell lines
#' @export
make_shiny_df <- function(con, gene, cell_lines, drug, data_types=c('affy', 'cn', 'hybcap', 'cosmicclp'), drug_df=NULL) {

  #gene <- 'TP53'
  #drug <- 'GSKMI00000714_pIC50'
  #cell_lines <- cls.df$CCLE_name
  #data_types <- c('affy', 'hybcap')
  #drug_df <- dnmt1_data

  gene_data <- make_tall_df(con,
                            gene,
                            cell_lines,
                            drug,
                            data_types) %>%
                transmute(CCLE_name, ID, feature_type=Type, feature_name=paste(ID, Type, sep="_"), feature_value=value)

  if (is.null(drug_df)) {
    stop('Need a drug dataframe')
  } else {
    resp_data <- getDrugData_custom(drug_df, drug, cell_lines ) %>%
      transmute(CCLE_name, resp_id=ID, resp_value=value)
  }

  cls.df <- src_sqlite(con@dbname) %>%
    tbl('ccle_sampleinfo') %>%
    transmute(CCLE_name, tissue=Site_primary) %>%
    collect

  plot_data <- gene_data %>%
    inner_join(resp_data, by='CCLE_name') %>%
    inner_join(cls.df, by='CCLE_name')

}

#' Convert tissue to cell lines
#'
#' This function creates a \code{vector} of cell lines from a tissue
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{vector} of cell lines for the provided tissue
#' @export
get_shiny_cell_lines <- function (con, tissue_id) {

  cls.df <- src_sqlite(con@dbname) %>%
    tbl('ccle_sampleinfo') %>%
    transmute(CCLE_name, tissue=Site_primary) %>%
    collect %>%
    filter(tissue %in% tissue_id)


  if(nrow(cls.df) == 0)  {
    warning('No cell lines found for that tissue')
    return(NULL)
  } else {
    return(cls.df$CCLE_name)
  }


}

plotAppUI <- function () {

  fluidPage(
    sidebarLayout(
      sidebarPanel(
        helpText("blah blah"),
        textInput("geneid", label = h3("Enter Gene Name"),
                  value = "TP53"),
        uiOutput('respUI'),
        radioButtons("data_type", label = h3("Select a genomic data type"),
                     choices = list("Affy" = 'affy', "CCLE hybcap" = 'hybcap'),
                     selected = 'affy'),
        uiOutput('tissuesUI'),
       # h3('Other options'),
       # checkboxInput("log2", label = "Log2 scale", value = TRUE)
        submitButton('Submit')
      ),
      mainPanel(textOutput('text1'),
                plotOutput("plot1"),
                tableOutput('df')

      )
    )
  )
}

plotAppServer <- function(input, output, con, drug_df) {

  proc_cls <- reactive({
    intersect(get_shiny_cell_lines(con, input$tissue), drug_df$unified_id)
  })

  proc_data <- reactive({
    make_shiny_df(con, gene=input$geneid,
                  cell_lines=proc_cls(),
                  drug=input$respid,
                  data_types = input$data_type,
                  drug_df = drug_df)
  })

  #make a reactive tissue selection UI
  output$tissuesUI <- renderUI({

    tissues.df <- src_sqlite(con@dbname) %>%
      tbl('ccle_sampleinfo') %>%
      filter(CCLE_name %in% drug_df$unified_id) %>%
      transmute(tissue=Site_primary) %>%
      collect %>%
      distinct() %>%
      arrange(tissue)


    tissues <- tissues.df$tissue

    checkboxGroupInput("tissue", label = h3("Select a tissue type"),
                 choices = tissues,
                 selected = tissues[1])

  })

  #make a reactive response variable selection UI
  output$respUI <- renderUI({

    resp.df <- drug_df %>% transmute(resp_id=paste(compound_id, endpoint, sep='_')) %>%
      distinct() %>% arrange(resp_id)

    resps <- resp.df$resp_id

    radioButtons("respid", label = h3("Select a response variable"),
                 choices = resps,
                 selected = resps[1])

  })

    output$text1 <- renderText({
      sprintf("Tissue is length: %s", length(input$tissue))
    })

  output$df <- renderTable({
    proc_data()
    #get_shiny_cell_lines(con, input$tissue) %>% as.data.frame
  })

  output$plot1 <- renderPlot({
    if(input$data_type == 'hybcap') {
      ggplot(proc_data(), aes(x=as.factor(feature_value), y=resp_value) ) + geom_boxplot() + facet_wrap(~tissue)
    } else {
      ggplot(proc_data(), aes(x=feature_value, y=resp_value) ) + geom_point() + stat_smooth(method = 'lm') + facet_wrap(~tissue)
    }

  })

}

#' Add title
#'
#' add description
#'
#' @param drug_df A data frame containing the drug data
#' @return returns
#' @export
run_shiny_app <- function(con, drug_df) {

  shinyApp(
    ui = plotAppUI(),
    server = function(input, output) {
      plotAppServer(input, output, con, drug_df)
    }
  )
}

