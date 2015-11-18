#' Waterfall plot for selected response variables
#'
#' Generate a set of waterfall plots for the selected response variables
#'
#' @param df A data frame from the \code{makeRespVsRespDataFrame} function
#' @return A ggplot2 object
#' @export
plotRespVsRespWaterfall <- function(df) {

  #warning if more than one ID
  if (length(unique(df$ID)) > 1) {
    selected_id <- unique(df$ID)[1]
    df <- df %>% filter(ID == selected_id)
    warning(sprintf('More than one response provided, using %s', selected_id))
  } else {
      selected_id <- unique(df$ID)
    }

  #sort plot data
  plot_data <- df %>% arrange(tissue, desc(value)) %>% filter(!is.na(value))

  #make a custom pallete for cell lines
  require(RColorBrewer)
  mypal <- brewer.pal(8,'Set2')
  set.seed(813379)
  mypal <- sample(colorRampPalette(mypal)(length(unique(plot_data$tissue))))

  #do waterfall plot ordered by potency
  histo_plot <- ggplot(plot_data, aes(x=CCLE_name, y=value)) + geom_bar(stat='identity', aes(fill = tissue)) +
    scale_x_discrete(limits=rev(plot_data$CCLE_name)) +
    #scale_y_log10(breaks=c(10,100,1000,10000)) +
    scale_fill_manual(values=mypal) +
    ggtitle(sprintf('Waterfall plot for %s',selected_id)) +
    theme_bw() +
    theme(#axis.text.x = element_text(size=rel(1.2), angle=90, hjust=1, vjust=0),
      axis.title.y = element_text(size=rel(1)),
      axis.title.x = element_text(size=0),
      axis.text.x = element_text(size=rel(0.6), angle=90, hjust=1, vjust=0),
      axis.text.y = element_text(size=rel(0.9)),
      axis.ticks.x = element_line(size=0),
      #    plot.margin=unit(c(0.5,5.5,0,1,0.7), "lines"),
      #    legend.key.size = unit(0.8, 'lines'),
      legend.position = 'right',
      legend.title = element_text(size=0))
  histo_plot

}
