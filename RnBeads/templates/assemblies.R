########################################################################################################################
## annotations.R
## created: 2012-08-16
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Collection of helper constants and functions related to the management of probe and region annotations.
########################################################################################################################

#' RnBeads Annotation Tables
#'
#' RnBeads uses sets of annotation tables and mappings (from regions to sites) for each of the supported genomes. The
#' structures for one assembly are stored in a separate dedicated data package. Currently, the following assemblies are
#' supported:
#' \describe{
#'   \item{\code{"hg19"}}{through the package \pkg{RnBeads.hg19}}
#'   \item{\code{"mm10"}}{through the package \pkg{RnBeads.mm10}}
#'   \item{\code{"mm9"}}{through the package \pkg{RnBeads.mm9}}
#'   \item{\code{"rn5"}}{through the package \pkg{RnBeads.rn5}}
#' }
#' 
#' @details
#' The assembly-specific structures are automatically loaded upon initialization of the annotation, that is, by the
#' first valid call to any of the following functions: \code{\link{rnb.get.chromosomes}},
#' \code{\link{rnb.get.annotation}}, \code{\link{rnb.set.annotation}}, \code{\link{rnb.get.mapping}},
#' \code{\link{rnb.annotation.size}}. Adding an annotation amounts to attaching its table(s) and mapping structures to
#' the scaffold.
#'
#' @docType data
#' @keywords datasets
#' @name RnBeads.data
#' @aliases hg19 mm10 mm9 rn5
#' @format \code{list} of four elements - \code{"regions"}, \code{"sites"}, \code{"controls"} and \code{"mappings"}.
#'         These elements are described below.
#'         \describe{
#'           \item{\code{"regions"}}{\code{list} of \code{NULL}s; the names of the elements correspond to the built-in
#'                region annotation tables. Once the default annotations are loaded, the attribute \code{"builtin"} is
#'                a \code{logical} vector storing, for each region annotation, whether it is the default (built-in) or
#'                custom.}
#'           \item{\code{"sites"}}{\code{list} of \code{NULL}s; the names of the elements correspond to the site and
#'                probe annotation tables.}
#'           \item{\code{"controls"}}{\code{list} of \code{NULL}s; the names of the elements correspond to the control
#'                probe annotation tables. The attribute \code{"sites"} is a \code{character} vector pointing to the
#'                site annotation that encompasses the respective control probes.}
#'           \item{\code{"mappings"}}{\code{list} of \code{NULL}s; the names of the elements correspond to the built-in
#'                region annotation tables.}
#'         }
#' @author Yassen Assenov
NULL

## G L O B A L S #######################################################################################################

## Environment to contain all probe, site and region annotation tables.
##
##  hg19
##      $regions
##          $tiling             GRangesList
##          $genes              GRangesList
##          $promoters          GRangesList
##          $cpgislands         GRangesList
##      $sites
##          $CpG                GRangesList
##          $probes450          GRangesList
##      $controls
##          $controls450        data.frame
##      $mappings
##          $tilinig
##              $CpG            list of IRanges
##              $probes450      list of IRanges
##          $genes
##              $CpG            list of IRanges
##              $probes450      list of IRanges
##          $promoters
##              $CpG            list of IRanges
##              $probes450      list of IRanges
##          $cpgislands
##              $CpG            list of IRanges
##              $probes450      list of IRanges
##      $lengths                int[ <chromosomes> , <annotations> ]
.rnb.annotations <- new.env()

## Chromosomes supported by the annotation packages
##%(chromosomes)s
CHROMOSOMES.L2S <- list("hg19" = c(1:22, "X", "Y"), "mm9" = c(1:19, "X", "Y"), "mm10" = c(1:19, "X", "Y"),
	"rn5" = c(1:20, "X")
	##%(assembly_table)s
	)
CHROMOSOMES.S2L <- lapply(CHROMOSOMES.L2S, function(x) { paste0("chr", x) })
CHROMOSOMES <- CHROMOSOMES.S2L
for (assembly.name in names(CHROMOSOMES)) {
	names(CHROMOSOMES.S2L[[assembly.name]]) <- CHROMOSOMES.L2S[[assembly.name]]
	names(CHROMOSOMES[[assembly.name]]) <- names(CHROMOSOMES.L2S[[assembly.name]]) <- CHROMOSOMES[[assembly.name]]
}
rm(assembly.name)

## Control probe types
HM450.CONTROL.TARGETS <- c(
	"bisulfite conversion I" = "BISULFITE CONVERSION I",
	"bisulfite conversion II" = "BISULFITE CONVERSION II",
	"extension" = "EXTENSION",
	"hybridization" = "HYBRIDIZATION",
	"negative control" = "NEGATIVE",
	"non-polymorphic" = "NON-POLYMORPHIC",
	"norm A" = "NORM_A",
	"norm C" = "NORM_C",
	"norm G" = "NORM_G",
	"norm T" = "NORM_T",
	"specificity I" = "SPECIFICITY I",
	"specificity II" = "SPECIFICITY II",
	"staining" = "STAINING",
	"target removal" = "TARGET REMOVAL")


HM27.CONTROL.TARGETS<-c(
		"bisulfite conversion" = "Bisulfite conversion",
		"extension" = "Extension",
		"hybridization" = "Hybridization",
		"negative control" = "Negative",
		"SNP" = "Genotyping",
		"non-polymorphic" = "Non-Polymorphic",
		"norm Grn" = "Normalization-Green",
		"norm Red" = "Normalization-Red",
		"specificity" = "Specificity",
		"staining" = "Staining",
		"target removal" = "Target Removal",
		"pACYC174" = "pACYC174",
		"pUC19" = "pUC19",
		"phiX174" = "phiX174"
		)

## Sample-independent control probe types (subset of CONTROL.TARGETS)
CONTROL.TARGETS.SAMPLE.INDEPENDENT <- c("STAINING", "HYBRIDIZATION", "TARGET REMOVAL", "EXTENSION")

## Genotyping probes on the 27k microarray

HM27.CY3.SNP.PROBES<-c(
		"rs798149",
		"rs2959823",
		"rs2235751",
		"rs2125573",
		"rs2804694"
		)
		
HM27.CY5.SNP.PROBES<-c(
		"rs1941955",
		"rs845016",
		"rs866884",
		"rs739259",
		"rs1416770",
		"rs1019916",
		"rs2521373",
		"rs10457834",
		"rs6546473",
		"rs5931272",
		"rs264581"
		)

## F U N C T I O N S ###################################################################################################

#' get.genome.data
#'
#' Gets the specified genome.
#'
#' @param assembly Genome assembly of interest. Currently the only supported genomes are \code{"hg19"}, \code{"mm9"},
#'                 \code{"mm10"} and \code{"rn5"}.
#' @return Sequence data object for the specified assembly.
#'
#' @author Yassen Assenov
#' @noRd
get.genome.data <- function(assembly) {
  if (assembly == "hg19") {
    suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg19))
    genome.data <- Hsapiens
  } else if (assembly == "mm9") {
    suppressPackageStartupMessages(require(BSgenome.Mmusculus.UCSC.mm9))
    genome.data <- Mmusculus
  } else if (assembly == "mm10") {
    suppressPackageStartupMessages(require(BSgenome.Mmusculus.UCSC.mm10))
    genome.data <- Mmusculus
  } else if (assembly == "rn5") {
    suppressPackageStartupMessages(require(BSgenome.Rnorvegicus.UCSC.rn5))
    genome.data <- Rnorvegicus
##%(assembly_package)s
  else {
    stop("unsupported assembly")
  }
  return(genome.data)
}