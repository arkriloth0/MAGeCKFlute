#' Violin plot
#'
#' Violin plot showing the distribution of numeric vectors with the same length.
#'
#' @docType methods
#' @name ViolinView
#' @rdname ViolinView
#' @aliases violinview
#'
#' @param dat A data frame.
#' @param samples A character vector, specifying the columns in the \code{dat} for plotting.
#' @param main A character, specifying title.
#' @param ylab A character, specifying title of y-axis.
#' @param filename A character, specifying a file name to create on disk.
#' Set filename to be "NULL", if don't want to save the figure.
#' @param width Numeric, specifying width of figure.
#' @param height Numeric, specifying height of figure.
#' @param units The units of figure size, one of "in", "cm", "mm", "px".
#' @param ... Other available parameters in function 'ggsave'.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{DensityView}}
#'
#' @examples
#' file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/mle.gene_summary.txt")
#' dd = ReadBeta(file3)
#' ViolinView(dd[, -1])
#'
#'
#' @importFrom reshape2 melt scales label_wrap
#' @import ggplot2 
#'
#' @export
#'
ViolinView <- function(dat, samples = NULL, main = NULL, ylab = "Score",
                       filename=NULL, width=5, height=4, units="in", ...){
  requireNamespace("ggplot2", quietly=TRUE) || stop("need ggplot2 package")
  requireNamespace("reshape2", quietly=TRUE) || stop("need reshape2 package")
  requireNamespace("scales", quietly=TRUE) || stop("need scales package")
  if(!is.null(samples) && length(samples)>1){
    dat = dat[, samples]
  }
  
  dd1 = reshape2::melt(dat, id.vars=NULL)
  if(!"variable" %in% colnames(dd1)){
    dd1$variable = colnames(dat)
  }
  dd1$Sample = factor(dd1$variable)
  levels(dd1$Sample) = snakecase::to_parsed_case(levels(dd1$Sample), sep_out = " ", parsing_option = 0)
  ## Plotting
  p = ggplot(data=dd1, aes(x=Sample,y=value,color=Sample)) +
    geom_violin() +
    geom_boxplot(width=.1, outlier.colour=NA) +
    # theme_bw(base_size = 14) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none") +
    labs(x=NULL, y=ylab, title=main)  +
    scale_x_discrete(labels = scales::label_wrap(5))
  if(!is.null(filename)){
    ggsave(plot=p,filename=filename, units = units, width=width, height=height, ...)
  }
  return(p)
}