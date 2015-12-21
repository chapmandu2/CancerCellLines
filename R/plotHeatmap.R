#' Make a heatmap from a tall data frame
#'
#' This function creates a heatmap from the output of the makeTallDataFrame function.
#'
#' @param tall_df A data frame from makeTallDataFrame
#' @return A ggplot2 object
#' @export
#' @import ggplot2
plotHeatmap <- function(tall_df, order_feature=NULL) {

  #scale the plot data by feature name
  plot_data <- tall_df %>% dplyr::mutate(fn=paste(ID, Type, sep='_')) %>%
                           dplyr::group_by(fn) %>%
                           dplyr::mutate(zscore=ifelse(Type %in% c('cosmicclp', 'hybcap'), value, scale(value))) %>%
                           dplyr::ungroup()

  #work out min, mean and max plus a suitable offset for each type to create a unified scale
  scale_info <- plot_data %>% dplyr::group_by(Type) %>%
                              dplyr::summarise(min=min(zscore, na.rm=TRUE), max=max(zscore, na.rm=TRUE), mean=0) %>%
                              dplyr::ungroup() %>%
                              dplyr::mutate(type_offset=20*(row_number()-1),
                                     min=min+type_offset,
                                     mean=mean+type_offset,
                                     max=max+type_offset)

  #calculate the zscore offset
  plot_data <- plot_data %>%
    dplyr::inner_join(scale_info %>% dplyr::select(Type, type_offset), by='Type') %>%
    dplyr::mutate(zscore_offset=zscore + type_offset) %>%
    dplyr::select(-type_offset)

  #combine the calculated scale info with colour information to get the final scale
  scale_info_tall <- scale_info %>% dplyr::select(-type_offset) %>% tidyr::gather(Level, value, -Type, convert=TRUE)
  scale_colours <- readr::read_tsv(system.file("extdata", "scale_colours.txt", package = "CancerCellLines"))
  scale_info_final <- scale_info_tall %>% dplyr::inner_join(scale_colours, by=c('Type', 'Level')) %>% dplyr::arrange(value)

  #order for x and y axis
  ordered_feature_names <- plot_data %>% dplyr::select(Type, ID, fn) %>%
    dplyr::distinct() %>% dplyr::mutate(isresp=Type!='resp') %>%
    dplyr::arrange(isresp, Type, ID)

  if (is.null(order_feature)) {
    ordered_feature_names <- ordered_feature_names$fn
    order_feature <- ordered_feature_names[1]
  } else if (order_feature %in% ordered_feature_names$fn) {
    message('Using user specified feature to order cell lines')
    ordered_feature_names <- ordered_feature_names %>% filter(fn != order_feature)
    ordered_feature_names <- c(order_feature, ordered_feature_names$fn)
  } else {
    stop(sprintf('%s not a valid feature name', order_feature))
  }

  ordered_cell_lines <- plot_data %>% dplyr::filter(fn == order_feature & !is.na(value)) %>% dplyr::arrange(desc(value))

  #make the plot!
  p <- ggplot(plot_data, aes(x=CCLE_name, y=fn)) +
    geom_tile(aes(fill = zscore_offset), linetype=0 ) +
    scale_fill_gradientn(colours = scale_info_final$Colour, values = scales::rescale(scale_info_final$value)) +
    scale_y_discrete(limits=ordered_feature_names) + #this orders the y axis as we want it
    scale_x_discrete(limits=ordered_cell_lines$CCLE_name) + #this orders the x axis as we want it
    coord_flip() +
    theme_bw(base_size = 9) +
    labs(x = "", y = "") +
    theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1),
          axis.text.y = element_text(size=rel(1))) +
    theme(panel.background=element_rect(fill="gray90", color='white'), legend.position = "none")
  p
  return(p)

}
