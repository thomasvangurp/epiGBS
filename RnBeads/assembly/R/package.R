## D A T A S E T S #####################################################################################################

#' MM9 - Annotation tables
#'
#' Scaffold of annotation tables for MM9. This structure is automatically loaded upon initialization of the annotation,
#' that is, by the first valid call to any of the following functions: \code{\link{rnb.get.assemblies}},
#' \code{\link{rnb.get.chromosomes}}, \code{\link{rnb.get.annotation}}, \code{\link{rnb.set.annotation}},
#' \code{\link{rnb.get.mapping}}, \code{\link{rnb.annotation.size}}. Adding an annotation amounts to attaching its
#' table(s) and mapping structures to this scaffold.
#'
#' @docType data
#' @keywords datasets
#' @name mm9
#' @format \code{list} of four elements - \code{"regions"}, \code{"sites"}, \code{"controls"} and \code{"mappings"}.
#'         These elements are described below.
#'         \describe{
#'           \item{\code{"regions"}}{\code{list} of \code{NULL}s; the names of the elements correspond to the built-in
#'                region annotation tables. Once the default annotations are loaded, the attribute \code{"builtin"} is
#'                a \code{logical} vector storing, for each region annotation, whether it is the default (built-in) or
#'                custom.}
#'           \item{\code{"sites"}}{\code{list} of \code{NULL}s; the names of the elements correspond to the site
#'                annotation tables.}
#'           \item{\code{"mappings"}}{\code{list} of \code{NULL}s; the names of the elements correspond to the built-in
#'                region annotation tables.}
#'         }
#' @author Yassen Assenov
NULL
