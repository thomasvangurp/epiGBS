#' Analysis of genome-scale DNA methylation data with RnBeads
#'
#' RnBeads facilitates comprehensive analysis of various types of DNA methylation data at the genome scale. It extends
#' previous approaches for such analysis by high throughput capabilities, as well as presenting results in a
#' comprehensive, highly interpretable fashion. 
#'
#' The complete analysis can be performed by calling the function \code{\link{rnb.run.analysis}}. 
#' 
#' @references Yassen Assenov*, Fabian Mueller*, Pavlo Lutsik*, Joern Walter, Thomas Lengauer and Christoph Bock (2014) Compehensive Analysis of DNA Methylation Data with RnBeads, Nature Methods, in press.
#' @import methods MASS cluster RColorBrewer fields ggplot2 matrixStats IRanges GenomicRanges methylumi ff gridExtra limma
#' @importFrom gplots colorpanel
#' @importFrom gplots heatmap.2
#' @importFrom plyr rbind.fill
#' @docType package
#' @name RnBeads
NULL
