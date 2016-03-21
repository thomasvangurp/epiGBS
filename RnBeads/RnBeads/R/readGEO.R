########################################################################################################################
## readGEO.R
## created: 2012-03-27
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Loading HumanMethylation450K datasets from GEO.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

## Constructs the intersection of two or more sets
##
## @param sets Non-empty list of vectors of the same type. Each vector corresponds to a set of elements.
## @return Vector containing all values that appear in each of the given sets.
intersect.elements <- function(sets) {
	if (length(sets) == 1) {
		return(sets)
	}
	result <- sets[[1]]
	for (i in 2:length(sets)) {
		result <- intersect(result, sets[[i]])
	}
	return(result)
}

########################################################################################################################

#' read.geo.parse.characteristics_ch1
#'
#' Parses the sample information in all of the \code{characteristics_ch1} columns from a phenoData data frame as
#' obtained from \code{getGEO}.
#'
#' @param phenoData Parsed phenotypic data frame object as output by \code{getGEO}.
#' @return Phenotypic data frame with parsed sample information instead of \code{characteristics_ch1}.
#'
#' @seealso \code{\link[GEOquery]{getGEO}} in package \pkg{GEOquery}
#' @author Fabian Mueller
read.geo.parse.characteristics_ch1 <- function(phenoData) {
	#phenoData <- result[[1]]@phenoData@data
	ch1.cols <- which(grepl("characteristics_ch1",colnames(phenoData)))
	ch1.df <- phenoData[ch1.cols]
	ch1.df <- data.frame(lapply(ch1.df, as.character), stringsAsFactors=FALSE)
	#keys are the identifier before the ":" (excluding whitespaces)
	keys <- unique(na.omit(as.vector(apply(ch1.df,c(1,2),FUN=function(x){
		unlist(strsplit(x,":", fixed = TRUE))[1]
	}))))
	key.df <- data.frame(matrix(NA,nrow=nrow(ch1.df),ncol=length(keys)))
	colnames(key.df) <- keys
	ch1.df.res <- ch1.df
	for (kk in keys){
		for (ss in 1:nrow(ch1.df)){
			fitting.cols <- which(grepl(kk,ch1.df[ss,],fixed=TRUE))
			if (length(fitting.cols)>0){
				cur.split <- strsplit(ch1.df[ss,fitting.cols[1]],paste("\\s*",kk,"\\s*:\\s*",sep=""))[[1]]
				key.df[ss,kk] <- cur.split[length(cur.split)]
				ch1.df.res[ss,fitting.cols[1]] <- "" #set the corresponding entry to blanc in the residual table as it is now assigned to the key table
			} else {
				key.df[ss,kk] <- NA
			}
		}
	}
#	#get the maximum number of ch1 characteristics NOT included in the key table
#	n.unparsed <- apply(ch1.df.res,1,FUN=function(x){
#		sum(!is.na(x) & x!="")
#	})
#	max.res <- max(n.unparsed)
#	ch1.df.unparsed <- data.frame(matrix(NA,nrow=nrow(ch1.df.res),ncol=max.res))
#	for (i in 1:length(n.unparsed)){
#		if (n.unparsed[i]>0){
#			cur.unparsed <- ch1.df.res[ss,]
#			ch1.df.unparsed[i,1:n.unparsed[i]] <-  cur.unparsed[!is.na(cur.unparsed) & cur.unparsed!=""]
#		}
#	}
	res <- cbind(phenoData[-ch1.cols],key.df)
	return(res)
}

########################################################################################################################

#' read.geo
#'
#' Imports Infinium 450K data series from the Gene Expression Omnibus.
#'
#' @param accession                 Character string representing the GEO series for download and parsing. It must start
#'                                  with \code{"GSE"}.
#' @param filename                  File name of a previously downloaded GEO series matrix file or its gzipped
#'                                  representation (in which case the filename must end in \code{".gz"}). Other file
#'                                  formats, such as SOFT files, are not supported. Exactly one of \code{accession} or
#'                                  \code{filename} must be specified.
#' @param verbose                   Flag indicating if messages should be created informing about the progress. If the
#'                                  logger is initialized prior to calling this function, the informative messages are
#'                                  sent to the logger. Warnings and errors are not affected by this parameters, the
#'                                  function always outputs them.
#' @param destdir                   The destination directory for any downloads. Defaults to the
#'                                  (architecture-dependent) temporary directory. Keep in mind that GEO series can be
#'                                  demanding in terms of storage space.
#' @param parse.characteristics_ch1 Flag indicating if additional sample characteristics (if such exist) are to be
#'                                  parsed from the downloaded series matrix file(s).
#' @return \code{\linkS4class{RnBeadSet}} object with phenotypic and beta value information; \code{NULL} if the given
#'         series contain no Infinium450K samples.
#'
#' @seealso \code{getGEO} in package \code{GEOquery}
#' @author Yassen Assenov
#' @export
read.geo <- function(accession = NULL, filename = NULL, verbose = logger.isinitialized(), destdir = tempdir(),
	parse.characteristics_ch1 = TRUE) {

	rnb.require("GEOquery")
	if (verbose) {
		rnb.logger.start("Loading GEO Data Series")
	}

	## Download or load the data series from a file
	if (is.null(filename)) {
		if (is.null(accession)) {
			rnb.error("None of accession and filename specified")
		}
		if (!grepl("^GSE", accession)) {
			rnb.error("Invalid accession identifier")
		}
		result <- GEOquery::getGEO(GEO = accession, destdir = destdir)
	} else if (!is.null(accession)) {
		rnb.error("Only one of accession and filename must be specified")
	} else if (!grepl("_series_matrix\\.txt(\\.gz)?$", filename)) {
		rnb.error("Only series matrix files are supported")
	} else {
		accession <- filename
		result <- GEOquery::getGEO(filename = filename, destdir = destdir)
		result <- list(result)
	}
	if (verbose) {
		rnb.info(c("Loaded", length(result), "data series from", accession))
	}
	betas.all <- list()
	phenotype.all <- list()
	i <- 0

	## Extract all Infinium450K samples from the loaded series
	for (result.series in result) {

		## Validate platform
		pheno.data <- result.series@phenoData@data
		sinds.accepted <- which(pheno.data$platform_id %in% c("GPL13534", "GPL8490")) # GEO accession for Infinium450K
		if (length(sinds.accepted) == 0) {
			next
		}

		## Extract beta values and phenotype data
		i <- i + 1
		betas.all[[i]] <- exprs(result.series)[, sinds.accepted]
		if (length(sinds.accepted) == 1) {
			betas.all[[i]] <- t(betas.all[[i]])
			colnames(betas.all[[i]]) <- rownames(pheno.data)[sinds.accepted]
		}

		phenotype.all[[i]] <- pheno.data[sinds.accepted, ]
		if (verbose) {
			rnb.info(c("Extracted", length(sinds.accepted), "samples from series", i))
		}
		if (parse.characteristics_ch1){
			phenotype.all[[i]] <- read.geo.parse.characteristics_ch1(phenotype.all[[i]])
			if (verbose) {
				rnb.info("Parsed characteristics_ch1 sample information")
			}
		}
	}
	rm(result.series, pheno.data, sinds.accepted)

	## Convert to a single beta value matrix and single phenotype data.frame
	if (i == 0) {
		if (verbose) {
			rnb.logger.completed()
		}
		return(NULL)
	}
	if (i == 1) {
		betas <- betas.all[[1]]
		phenotype <- phenotype.all[[1]]
	} else { # i > 1
		## Combine beta value matrices
		rnames <- intersect.elements(lapply(betas.all, rownames))
		## TODO: Handle the case when length(rnames) == 0
		betas.all <- lapply(betas.all, function(x) { x[rnames, ] })
		cnames.all <- lapply(betas.all, colnames)
		cnames <- unlist(cnames.all, use.names = FALSE)
		betas <- matrix(as.double(NA), nrow = length(rnames), ncol = length(cnames), dimnames = list(rnames, cnames))
		for (j in 1:length(betas.all)) {
			betas[, cnames.all[[j]]] <- betas.all[[j]]
		}
		## Combine the data frames of pheontype data
		cnames <- intersect.elements(lapply(phenotype.all, colnames))
		## TODO: Handle the case when length(cnames) == 0
		phenotype <- phenotype.all[[1]][, cnames]
		for (j in 2:length(phenotype.all)) {
			phenotype <- rbind(phenotype, phenotype.all[[j]][, cnames])
		}
	}
	beta.range <- tryCatch(range(betas, na.rm = TRUE), warning = function(wr) { NULL })
	if (is.null(beta.range)) {
		## All beta values are NA or NaN; accept it
	} else if (!(0 <= beta.range[1] && beta.range[2] <= 1)) {
		betas<-rnb.mval2beta(betas)		
		rnb.warning("Loaded values are not in the range [0,1], treated as M-values and converted to beta-values")
	}	
	
	if (verbose) {
		rnb.logger.completed()
	}
	return(new("RnBeadSet", phenotype, betas, p.values = NULL, bead.counts = NULL))
}