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
makeRespVsGeneticDataFrame <- function(con, gene, cell_lines, drug, data_types=c('affy', 'cn', 'hybcap', 'cosmicclp'), drug_df=NULL) {

  #gene <- 'TP53'
  #drug <- 'GSKMI00000714_pIC50'
  #cell_lines <- cls.df$CCLE_name
  #data_types <- c('affy', 'hybcap')
  #drug_df <- dnmt1_data

  gene_data <- makeTallDataFrame(con,
                            gene,
                            cell_lines,
                            drug,
                            data_types) %>%
                transmute(CCLE_name, ID, feature_type=Type, feature_name=paste(ID, Type, sep="_"), feature_value=value)

  if (is.null(drug_df)) {
    #if no drug data frame provided just use CCLE
    ccle_drugs.df <- src_sqlite(con@dbname) %>% tbl("ccle_drug_data") %>% select(Compound) %>% distinct %>% collect
    resp_data <- getDrugData_CCLE(con, drug, cell_lines) %>%
      transmute(CCLE_name, resp_id=ID, resp_value=value)
  } else {
    #if drug data provided then use that
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

  return(plot_data)
}

#' Convert tissue to cell lines
#'
#' This function creates a \code{vector} of cell lines from a tissue
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{vector} of cell lines for the provided tissue
#' @export
getTissueCellLines <- function (con, tissue_id) {

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

shinyRespVsGeneticUI <- function () {

  fluidPage(
    sidebarLayout(
      sidebarPanel(
        helpText("Correlate genomic and response data"),
        textInput("geneid", label = h3("Enter Gene Name"),
                  value = "TP53"),
        uiOutput('respUI'),
        selectInput("data_type", label = h3("Select a genomic data type"),
                     choices = list("Affy" = 'affy', "CCLE hybcap" = 'hybcap', "Cosmic CLP" = 'cosmicclp'),
                     selected = 'affy'),
        selectInput("output_option", label = h3("Generate"),
                     choices = list('Plot 1' = 1, 'Plot 2'=2, 'Data'=3),
                     selected = 1),
        h3("Plot Options:"),
        checkboxInput("facet_option", label = "Facet by tissue"),
        uiOutput('tissuesUI')
       # submitButton('Submit')
      ),
      uiOutput('resultsUI')


    )
  )
}

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

    plot_data <- proc_data() %>% arrange(desc(resp_value))

    if(input$data_type %in% c('hybcap', 'cosmicclp')) {
      p <- ggplot(plot_data, aes(x=as.factor(feature_value), y=resp_value) ) +
              geom_boxplot(outlier.size=0, aes(colour=as.factor(feature_value))) +
              geom_point(aes(fill=as.factor(feature_value)), shape=21, size=rel(2), position = position_jitter(width=0.1)) +
              xlab(unique(plot_data$feature_name)) + ylab(unique(plot_data$resp_id)) +
              theme_bw()

      if (input$facet_option) {
        p <- p + facet_wrap(~tissue)
      }

      return(p)

    } else {
      p <- ggplot(plot_data, aes(x=feature_value, y=resp_value) ) +
              geom_point(aes(fill=scale(feature_value)), shape=21, size=rel(2)) +
              stat_smooth(method = 'lm') +
              scale_fill_gradient2(low='blue', mid='white', high='red') +
              xlab(unique(plot_data$feature_name)) + ylab(unique(plot_data$resp_id)) +
              theme_bw() + theme(legend.position='none')

      if (input$facet_option) {
        p <- p + facet_wrap(~tissue)
      }
      return(p)

      }

  })

  output$plot2 <- renderPlot({

    plot_data <- proc_data() %>% filter(!is.na(resp_value)) %>% arrange(desc(resp_value))

    if(input$data_type %in% c('hybcap', 'cosmicclp')) {
      p <- ggplot(plot_data, aes(x=CCLE_name, y=resp_value, fill=as.factor(feature_value)) ) +
        geom_bar(stat='identity') +
        scale_x_discrete(limits=plot_data$CCLE_name) +
        ylab(unique(plot_data$resp_id)) +
        theme_bw() + theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1))

      if(nrow(plot_data) > 30 | input$facet_option) {
        p <- p + theme(axis.text.x=element_text(size=0))  #remove cell line id's if more than 30 cell lines
      }

      if (input$facet_option) {
        p <- p + facet_wrap(~tissue)
      }

      return(p)

    } else {
      p <- ggplot(plot_data, aes(x=CCLE_name, y=resp_value, fill=scale(feature_value)) ) +
        geom_bar(stat='identity') +
        scale_x_discrete(limits=plot_data$CCLE_name) +
        scale_fill_gradient2(low='blue', mid='white', high='red') +
        ylab(unique(plot_data$resp_id)) +
        theme_bw() + theme(legend.position='none', axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1))

      if(nrow(plot_data) > 30 | input$facet_option) {
        p <- p + theme(axis.text.x=element_text(size=0))  #remove cell line id's if more than 30 cell lines
      }

      if (input$facet_option) {
        p <- p + facet_wrap(~tissue)
      }

      return(p)
    }

  })

}

#' Runs the Shiny visualisation
#'
#' Takes a user supplied data frame of drug response data and combines it with user defined gene feature data
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param drug_df A data frame containing the drug data
#' @return Launches an interactive Shiny application
#' @export
shinyRespVsGeneticApp <- function(con, drug_df=NULL) {

  shinyApp(
    ui = shinyRespVsGeneticUI(),
    server = function(input, output) {
      shinyRespVsGeneticServer(input, output, con, drug_df)
    }
  )
}

