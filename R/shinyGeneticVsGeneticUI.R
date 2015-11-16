shinyGeneticVsGeneticUI <- function () {

  fluidPage(
    sidebarLayout(
      sidebarPanel(
        helpText("Correlate different types of genetic data from one or more genes"),
        h3("Gene 1 Details:"),
        textInput("gene_name1", label = "Enter Gene Name",
                  value = "SMARCA4"),
        selectInput("data_type1", label = "Select a genomic data type",
                    choices = list("Affy" = 'affy', "CCLE hybcap" = 'hybcap', "Cosmic CLP" = 'cosmicclp', "Copy Number" = 'cn'),
                    selected = 'affy'),
        h3("Gene 2 Details:"),
        textInput("gene_name2", label = "Enter Gene Name",
                  value = "SMARCA4"),
        selectInput("data_type2", label = "Select a genomic data type",
                    choices = list("Affy" = 'affy', "CCLE hybcap" = 'hybcap', "Cosmic CLP" = 'cosmicclp', "Copy Number" = 'cn'),
                    selected = 'hybcap'),
        selectInput("output_option", label = h3("Generate"),
                    choices = list('Scatter/Boxplot' = 1, 'By Cell Line'=2, 'View data'=3),
                    selected = 1),
        h3("Plot Options:"),
        checkboxInput("facet_option", label = "Facet by tissue"),
        checkboxInput("label_option", label = "Label cell lines", value = FALSE),
        sliderInput('plot_width', label = 'Control plot width', min=200, max=2000, value=600),
        sliderInput('plot_height', label = 'Control plot height', min=200, max=2000, value=400),
        uiOutput('tissuesUI')
        # submitButton('Submit')
      ),
      uiOutput('resultsUI')


    )
  )
}
