########################################################################################################################
## dataset_generator.R
## created: 2012-11-27
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Routines for generating random data to be used in unit tests.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

## rnb.generate.betas
##
## Generates a matrix of random values in the range [0, 1].
##
## @param site.names   \code{character} vector to use for row names of the matrix.
## @param sample.names \code{character} vector to use for column names of the matrix.
## @return Matrix of beta values with dimension names \code{site.names} and \code{sample.names}. There are no missing
##         values (\code{NA}s) in the matrix.
##
## @author Yassen Assenov
rnb.generate.betas <- function(site.names, sample.names = LETTERS[1:12]) {
	N <- length(site.names) * length(sample.names)
	N.bump <- as.integer(0.03 * N) # maximum number of values to direct towards the "bump" at 0.5
	values.bump <- rnorm(N.bump, mean = 0.5, sd = 0.15)
	values.bump <- values.bump[0 <= values.bump & values.bump <= 1]
	values.peaks <- rbeta(N - length(values.bump), 0.5, 0.5, ncp = 0)
	values <- sample(c(values.bump, values.peaks))
	matrix(values, nrow = length(site.names), ncol = length(sample.names), dimnames = list(site.names, sample.names))
}

########################################################################################################################

## rnb.generate.intensities
##
## Generates a list of two matrices with intensity values corresponding to a given beta-value matrix. The matrix should
## have row and column names
##
## @param beta.values   \code{numeric} matrix of beta-values
## @param average.intensity \code{numeric} vector of length 1 giving mean total intensity
## @param sd.intensity	\code{numeric} vector of length 1 giving standard deviation of the total intensity
##
## @return List of two matrices with dimension equal to \code{beta.values}. There are no missing
##         values (\code{NA}s) in the matrix.
##
## @author Yassen Assenov
rnb.generate.intensities <- function(beta.values, average.intensity=20000, sd.intensity=2000) {
	
	nsites<-nrow(beta.values)
	nsamples<-ncol(beta.values)
	
	site.total<-rnorm(nsites, average.intensity, sd.intensity)
	
	total<-t(sapply(site.total, rnorm, n=nsamples, sd=sd.intensity/100))
	colnames(total)<-colnames(beta.values)
	rownames(total)<-rownames(beta.values)
	
	methylated<-beta.values * total
	unmethylated<-(1-beta.values) * total
			
	list(M=methylated, U=unmethylated)
}

########################################################################################################################

## rnb.generate.pvalues
##
## Generates a matrix of random p-values, the vast majority of which are very small.
##
## @param site.names   \code{character} vector to use for row names of the matrix.
## @param sample.names \code{character} vector to use for column names of the matrix.
## @return Matrix of p-values values with dimension names \code{site.names} and \code{sample.names}. There are no
##         missing values (\code{NA}s) in the matrix.
##
## @author Yassen Assenov
rnb.generate.pvalues <- function(site.names, sample.names = LETTERS[1:12]) {
	N <- length(site.names) * length(sample.names)
	values <- rexp(as.integer(N * 0.1), rate = 12) / 12
	values <- sample(c(rep(min(values) / 2, N - length(values)), values))
	matrix(values, nrow = length(site.names), ncol = length(sample.names), dimnames = list(site.names, sample.names))
}

########################################################################################################################

## rnb.generate.annotation
##
## Generates an annotation table for 12 samples coming from the same array.
##
## @return Sample annotation table in the form of a \code{data.frame}.
##
## @author Yassen Assenov
rnb.generate.annotation <- function() {
	sentrix.id <- paste("1", paste(sample.int(10, size = 9, replace = TRUE) - 1L, collapse = ""), sep = "")
	sentrix.id <- as.integer(sentrix.id)
	sentrix.pos <- paste("R0", rep(1:6, each = 2), "C0", rep(1:2, 6), sep = "")
	tissue <- "lung"
	disease.states <- c("healthy", "primary tumor")
	disease <- factor(disease.states[c(1, 1, 2, 2, 2, 1, 2, NA, NA, 2, 2, 2)])
	diagnoses <- c("2012-11-26", "2012-11-26", "2012-11-26", "2012-11-26", "2012-11-27", "2012-11-27", NA, NA,
		"2012-12-01", "2012-12-01", "2012-12-01", "2012-12-02")
	result <- data.frame(
		"Sentrix_ID" = sentrix.id,
		"Sentrix_Position" = sentrix.pos,
		"Tissue" = tissue,
		"Disease" = disease,
		"Diagnose Date" = diagnoses,
		"barcode" = paste(sentrix.id, sentrix.pos, sep = "_"), check.names = FALSE, stringsAsFactors = FALSE)
	rownames(result) <- LETTERS[1:nrow(result)]
	return(result)
}

########################################################################################################################

## rnb.generate.infinium450.dataset
##
## Generates an object of type \code{RnBeadSet} containing approximately 4800 probes and 12 samples.
##
## @return The newly initialized dataset.
##
## @author Yassen Assenov
rnb.generate.infinium450.dataset <- function() {
	## Select randomly 1% of the probe names
	annot <- rnb.annotation2data.frame(rnb.get.annotation("probes450"))
	site.names <- rownames(annot)[sort(sample(1:nrow(annot), size = round(nrow(annot) * 0.01)))]
	## Generate data matrices and sample annotation table
	beta.matrix <- rnb.generate.betas(site.names)
	pval.matrix <- rnb.generate.pvalues(site.names)
	phenotypes <- rnb.generate.annotation()
	## Initialize an object
	new("RnBeadSet", pheno = phenotypes, betas = beta.matrix, p.values = pval.matrix, bead.counts = NULL)
}

## rnb.generate.infinium450.raw.dataset
##
## Generates an object of type \code{RnBeadRawSet} containing approximately 4800 probes and 12 samples.
##
## @return The newly initialized dataset.
##
## @author Yassen Assenov
rnb.generate.infinium450.raw.dataset <- function() {
	## Select randomly 1% of the probe names
	annot <- rnb.annotation2data.frame(rnb.get.annotation("probes450"))
	site.names <- rownames(annot)[sort(sample(1:nrow(annot), size = round(nrow(annot) * 0.01)))]
	## Generate data matrices and sample annotation table
	beta.matrix <- rnb.generate.betas(site.names)
	intensities <- rnb.generate.intensities(beta.matrix)
	pval.matrix <- rnb.generate.pvalues(site.names)
	phenotypes <- rnb.generate.annotation()
	## Initialize an object
	new("RnBeadRawSet", pheno = phenotypes, M=intensities$M, U=intensities$U, p.values = pval.matrix, bead.counts.M = NULL, bead.counts.U = NULL)
}