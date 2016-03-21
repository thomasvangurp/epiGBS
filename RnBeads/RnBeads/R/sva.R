## get.components.sva
##
## Computes Surrogate Variables for a given dataset and target
##
## @param df.data an MxN data.frame containing the values for M sites on N samples
## @param df.pheno an NxK annotation table containing sample annotations
## @param pheno.col the column name in df.pheno containing the target variable for SVA
## @param conf.cols a vector of confounding variable column names
## @param numSVmethod Method for estimating the number of surrogate variables. passed directly to \code{sva}
## @return A matrix containing the computed Surrogate Variables. Rows correspond to samples, columns correspond to the SVs
##
## @author Fabian Mueller
get.components.sva <- function(df.data,df.pheno,pheno.col,conf.cols=NULL,numSVmethod="leek"){
	require(sva)
	if (length(conf.cols)>0){
		logger.info(paste0("Adjusting for known covariates: ",paste(conf.cols,collapse=",")))
	}
	n.samples <- nrow(df.pheno)
	sel.samples <- !is.na(df.pheno[,pheno.col]) #remove sample with no primary variable from dataset
	df.data <- df.data[,sel.samples] 
	df.pheno <- df.pheno[sel.samples,]
	fff.txt <- paste0("~1+",paste(pheno.col,collapse="+"))
	mod <- model.matrix(as.formula(fff.txt),data=df.pheno)
	fff0.txt <- "~1"
	if (length(conf.cols)>0){
		fff0.txt <- paste0("~0+",paste(conf.cols,collapse="+"))
	}
	mod0 <- model.matrix(as.formula(fff0.txt),data=df.pheno)
	svobj <- NULL
	tryCatch(
		svobj <- sva(df.data,mod,mod0,numSVmethod=numSVmethod),
		error = function(ee) {
			logger.warning(c("Could not compute SVA components:",ee$message))
	})
	valid.res <- !is.null(svobj)
	if (valid.res && is.null(ncol(svobj$sv))){
		logger.info("No surragate variables found")
		valid.res <- FALSE
	}
	if (valid.res){
		res <- matrix(NA,n.samples,ncol(svobj$sv))
		res[sel.samples,] <- svobj$sv
	} else {
		res <- matrix(NA,n.samples,0)
	}
	return(res)
}

## get.components.isva
##
## Computes Independent Surrogate Variables for a given dataset and target
##
## @param df.data an MxN data.frame containing the values for M sites on N samples
## @param df.pheno an NxK annotation table containing sample annotations
## @param pheno.col the column name in df.pheno containing the target variable for SVA
## @param conf.cols a vector of confounding variable column names
## @return A matrix containing the computed Surrogate Variables. Rows correspond to samples, columns correspond to the SVs
##
## @author Fabian Mueller
get.components.isva <- function(df.data,df.pheno,pheno.cols,conf.cols=NULL){
	require(isva)
	if (length(conf.cols)>0){
		conf.is.factor <- apply(df.pheno[,conf.cols],2,FUN=function(x){
			is.factor(x) || is.logical(x)
		})
		isvobj <- DoISVA(df.data,df.pheno[,pheno.cols],df.pheno[,conf.cols],conf.is.factor)
	} else {
		# logger.start("Computing components: ISVA")
		isvobj <- DoISVA(df.data,df.pheno[,pheno.cols])
		# logger.completed()
	}
	return(isvobj$isv)
}

## get.pheno.tab
##
## helper function thart retrieves a list of columns from the sample annotation table
## that will be used in the batch effects and SVA steps of RnBeads
##
## @param rnb.set an \code{RnBSet} object
## @return A list of inferred phenotype data
##
## @author Yassen Assenov, Fabian Mueller
get.pheno.tab <- function(rnb.set){
	pheno.table <- pheno(rnb.set)
	if (is.null(pheno.table)) {
		stop("invalid value for rnb.set; missing phenotype data")
	}
	
	## Get a list of comparable traits
	predefined.columns <- rnb.getOption("exploratory.columns")
	if (!is.null(predefined.columns)) {
		if (is.character(predefined.columns)) {
			predefined.columns <- intersect(predefined.columns, colnames(pheno.table))
		} else { # is.integer(columns)
			predefined.columns <- intersect(predefined.columns, 1:ncol(pheno.table))
		}
		pheno.table <- pheno.table[, predefined.columns]
	}
	traits <- lapply(pheno.table, function(x) {
				if (length(unique(na.omit(x))) < 2) {
					return(NULL)
				}
				if (is.character(x)) {
					x <- as.factor(x)
				}
				if (is.logical(x)) {
					x <- as.factor(x)
				}
				if (is.factor(x)) {
					y <- na.omit(x)
					if (anyDuplicated(y) == 0) {
						return(NULL)
					}
				}
				return(x)
			}
	)
	traits <- traits[sapply(traits, is.null) == FALSE]
	return(traits)
}

## compute.sva.assoc
##
## Helper function to compute the associations of Surrogate Variables (SVs)
## with principal components and sample annotations
## implemented in analogy to the pca analysis by Yassen
##
## @param meth.vals an MxN data.frame containing the values for M sites on N samples. Should not contain missing values.
## @param pheno.list a list of sample annotations as returned by \code{get.pheno.tab}
## @param sva.tables a named list of computed SVA tables as returned by \code{get.components.sva}. It should contain one element
##					 per target variable
## @return An object of class \code{SvaAssoc}: basically a list whose only element is named \code{sva}. This element is
## 		   again a list containing containing the association information for each target variable SVA.
##         The association information in turn consists of association with the Principal components (element \code{pca}) and
##         with the sample annation (element \code{traits}). These elements contain matrices with correlation coefficients,
##         the name of the association test, the resulting p-value from an association test and names of errors that occurred.
##
## @author Fabian Mueller
compute.sva.assoc <- function(meth.vals,pheno.list,sva.tables){
	result <- list()
	logger.start("Computing components: PCA")
	pcaobj <- prcomp(t(meth.vals), center = TRUE, scale. = FALSE)
	pca.tab <- pcaobj$x
	logger.completed()
	n.pc <- 3
	n.ph <- length(pheno.list)
	result[["sva"]] <- list()
	for (cc in names(sva.tables)){
		result[["sva"]][[cc]] <- list()
		logger.start(c("Target variable:",cc))
		n.sv <- ncol(sva.tables[[cc]])
		if (n.sv > 0){
			sv.tab <- sva.tables[[cc]]
			names.svs <- paste0("sv",1:n.sv)
			# names.svs <- 1:n.sv
			names.pcs <- paste0("pc",1:n.pc)
			# names.pcs <- 1:n.pc
			names.traits <- names(pheno.list)
			cor.tab.pca  <- matrix(nrow = n.pc, ncol = n.sv, dimnames = list("Principal component" = names.pcs, "Surrogate Variable" = names.svs))
			cor.tab.ph   <- matrix(nrow = n.ph, ncol = n.sv, dimnames = list("Sample Trait" = names.traits, "Surrogate Variable" = names.svs))
			test.tab.ph  <- matrix(nrow = n.ph, ncol = n.sv, dimnames = list("Sample Trait" = names.traits, "Surrogate Variable" = names.svs))
			pval.tab.ph  <- matrix(nrow = n.ph, ncol = n.sv, dimnames = list("Sample Trait" = names.traits, "Surrogate Variable" = names.svs))
			err.tab.ph   <- matrix(nrow = n.ph, ncol = n.sv, dimnames = list("Sample Trait" = names.traits, "Surrogate Variable" = names.svs))
			for (j in 1:n.sv){
				for (i in 1:n.pc){
					t.result <- test.traits(sv.tab[,j],pca.tab[,i])
					cor.tab.pca[i,j] <- t.result[["correlation"]]
				}
				for (i in 1:n.ph){
					t.result <- test.traits(sv.tab[,j],pheno.list[[i]])
					cor.tab.ph[i,j]  <- t.result[["correlation"]]
					test.tab.ph[i,j] <- t.result[["test"]]
					pval.tab.ph[i,j] <- t.result[["pvalue"]]
					err.tab.ph[i,j]  <- t.result[["error"]]
				}
			}
		} else {
			cor.tab.pca <- NA
			cor.tab.ph  <- NA
			test.tab.ph <- NA
			pval.tab.ph <- NA
			err.tab.ph  <- NA
		}
		result[["sva"]][[cc]][["pca"]] <- list(
					"correlations" = cor.tab.pca)
		result[["sva"]][[cc]][["traits"]] <- list(
					"failures" = err.tab.ph,
					"tests" = test.tab.ph,
					"correlations" = cor.tab.ph,
					"pvalues" = pval.tab.ph)

		logger.completed()
	}
	class(result) <- "SvaAssoc"
	return(result)
}

#' rnb.execute.sva
#'
#' Conduct Surrogate Variable Analysis (SVA) on the beta values of an RnBSet for given target variables
#'
#' @param rnb.set The \code{RnBSet} object on which the SVA should be conducted
#' @param cmp.cols a vector of sample annotation column names which will be the targets of the SVA.
#' @param columns.adj Column names in the table of phenotypic information to be used for confounder adjustment.
#' @param assoc a flag indicating whether association information with principal components and other sample annotation should
#'				be returned
#' @param numSVmethod method to estimate the number of surrogate variables. Passed to \code{sva}.
#' @return An object of class \code{SvaResult}: basically a list containing the following elements:
#' \describe{
#'   \item{\code{num.components}}{a vector storing the number of detected SVs for each target variable}
#'   \item{\code{sva.performed}}{a vector storing whether SVA was performed on a target variable and whether more than 0 SVs were found}
#'   \item{\code{targets}}{a vector storing the names of the target variables}
#'   \item{\code{components}}{a list storing for each target variable a matrox containing the sample-wise SVs as rows}
#'   \item{\code{assoc}}{a special object containing association information of SVs with principal components and sample annotations
#'	 typically only used \code{rnb.section.sva}.}
#' }
#' @export
#' @author Fabian Mueller
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' sva.obj <- rnb.execute.sva(rnb.set.example,c("Sample_Group","Treatment"),numSVmethod="be")
#' sva.obj$sva.performed
#' sva.obj$num.components
#' rnb.set.mod <- set.covariates.sva(rnb.set.example, sva.obj)
#' has.covariates.sva(rnb.set.example,"Sample_Group")
#' has.covariates.sva(rnb.set.mod,"Sample_Group")
#' has.covariates.sva(rnb.set.mod,"Treatment")
#' }
rnb.execute.sva <- function(rnb.set, cmp.cols=rnb.getOption("inference.targets.sva"),
		columns.adj=rnb.getOption("covariate.adjustment.columns"), assoc=TRUE,
		numSVmethod=rnb.getOption("inference.sva.num.method")){
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!(is.character(cmp.cols) || is.null(cmp.cols))){
		stop("invalid value for cmp.cols")
	}
	result <- list()
	logger.start("Conducting SVA analysis")

	if (is.null(cmp.cols)){
		cmp.cols <- rnb.getOption("inference.targets.sva")
	}

	ph <- pheno(rnb.set)
	ph.cols <- colnames(ph)
	
	if (is.character(cmp.cols)){
		cmp.cols.ph <- intersect(ph.cols,cmp.cols)
		if (length(cmp.cols.ph)<length(cmp.cols)){
			logger.warning(c("The following columns are not specified in the sample annotation table:",
							 paste(setdiff(cmp.cols,cmp.cols.ph),collapse=","),
							 "--> ignoring them in the SVA adjustment"))
		}
		cmp.cols <- cmp.cols.ph
	}
	logger.info(c("Considering the following target variables:",paste(cmp.cols,collapse=",")))
	logger.start("Preparing data")
		mm <- meth(rnb.set)
		nonzero.var.rows <- rowVars(mm,na.rm=TRUE) > 0
		mm.vao <- mm[nonzero.var.rows,]
		logger.info(c("Removed",nrow(mm)-sum(nonzero.var.rows),"sites without variation"))
		mm.nao <- na.omit(mm.vao)
		logger.info(c("Removed",nrow(mm.vao)-nrow(mm.nao),"sites containing missing values"))
	logger.completed()

	logger.start("Computing components: SVA")
		sva.tables <- list()
		for (cc in cmp.cols){
			logger.start(c("Target variable:",cc))
			adj.cols <- setdiff(columns.adj,cc)
			sv <- get.components.sva(mm.nao,ph,cc, conf.cols=adj.cols, numSVmethod=numSVmethod)
			sva.tables[[cc]] <- sv
			logger.completed()
		}
		num.components <- sapply(sva.tables,ncol)
		names(num.components) <- names(sva.tables)
		computed.sva <- num.components > 0
		result[["num.components"]] <- num.components
		result[["sva.performed"]] <- computed.sva
		logger.info(c("Found SVs for ",sum(computed.sva),"comparisons"))

	logger.completed()
	sva.assoc <- list()
	if (assoc && any(computed.sva)){
		logger.start("Computing association of SVA with PCA and sample annotations")
			pheno.list <- get.pheno.tab(rnb.set)
			sva.assoc <- compute.sva.assoc(mm.nao,pheno.list,sva.tables)
		logger.completed()
	}

	result[["targets"]] <- cmp.cols
	result[["components"]] <- sva.tables
	result[["assoc"]] <- sva.assoc
	class(result) <- "SvaResult"
	logger.completed()
	return(result)
}

#' rnb.section.sva
#'
#' Adds a section on Surrogate Variable Analysis (SVA) to an RnBeads Report
#' @param report An object of class \code{\linkS4class{Report}}.
#' @param sva.obj An object of class \code{SvaResult} as returned by \code{rnb.execute.sva}.
#' @return The modified report
#' @author Fabian Mueller
#' @noRd
rnb.section.sva <- function(report, sva.obj) {
	if (!inherits(sva.obj,"SvaResult")){
		stop("invalid value for sva.obj. Expected SvaResult!")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	do.it <- any(sva.obj$sva.performed)

	txt <- c("Surrogate variables were tested for association ",
			"with principal components and with sample traits using the specified target variables. The table below lists the ",
			"target variables along with the number of computed Surrogate Variables (SVs)")
	if (!do.it) {
		txt <- "Surrogate variable analysis did not result in any surrogate variables."
	}
	report <- rnb.add.section(report, "Surrogate Variable Analysis", txt, level = 1)
	
	if (do.it) {
		no.sva.inds <- which(sva.obj$num.components<1)
		targets <- sva.obj$targets
		names(targets) <- paste0("t",1:length(targets))
		tbl <- data.frame(targets,sva.obj$num.components[targets])
		colnames(tbl) <- c("Target","#SV")
		rnb.add.table(report, tbl)

		targets <- targets[sva.obj$sva.performed]
		setting.names <- list("Target" = targets)

		data.dir.abs <- rnb.get.directory(report, "data", TRUE)
		data.dir.rel <- rnb.get.directory(report, "data")

		txt <- c("Surrogate variables were tested for association ",
				 "with sample traits using the tests outlined in the table below.")
		report <- rnb.add.section(report, "Association of surrogate variables with sample traits", txt, level = 2)
		## Display summary table of tests performed
		i <- which(sva.obj$sva.performed)[1]
		assoc.tables <- sva.obj$assoc$sva[[i]][["traits"]]
		tbl <- matrix(c(rownames(assoc.tables$tests), assoc.tables$tests[,1]), ncol = 2)
		colnames(tbl) <- c("Trait", "Test")
		rnb.add.table(report, tbl)

		append.table <- function(tbl, fname.plot, fname.tab, heatmap.fun, header.r) {
			mytbl <- cbind(rownames(tbl), as.data.frame(tbl, check.names = FALSE))
			colnames(mytbl)[1] <- paste0(header.r, " \\ SV") 
			fname.tab.abs <- file.path(data.dir.abs,fname.tab)
			fname.tab.rel <- paste(data.dir.rel, fname.tab, sep = "/")
			utils::write.csv(mytbl, file = fname.tab.abs, na = "", row.names = FALSE)
			sprintf("<a href=\"%s\">csv</a>", fname.tab.rel)
			hmap <- heatmap.fun(report, tbl, fname.plot)
			list(link = sprintf("<a href=\"%s\">csv</a>", fname.tab.rel), plot = hmap$plot, description = hmap$description)
		}

		tab.header.r <- "Trait"
		#Traits: correlations
		file.tbl <- matrix(character(), nrow = 0, ncol = 2, dimnames = list(NULL, c("Target", "File Name")))
		rplots <- list()
		for (i in 1:length(targets)){
			tt <- targets[i]
			tn <- names(targets)[i]
			tbl <- sva.obj$assoc$sva[[tt]][["traits"]]$correlations
			if (!all(is.na(tbl))) {
				#reassign names to make plot layout better
				names.traits <- rownames(tbl)
				tbl <- unname(tbl)
				dimnames(tbl) = list(names.traits,"Surrogate Variable"=1:ncol(tbl))
				fname.plot <- paste0("heatmap_cor_sva_trait_",tn)
				fname.tab  <- paste0("cor_sva_trait_",tn,".csv")
				result <- append.table(tbl, fname.plot, fname.tab, plot.heatmap.pc.correlations, tab.header.r)
				file.tbl <- rbind(file.tbl, c(targets[tn], result$link))
				rplots[[targets[tn]]] <- result$plot
				description <- result$description
			}
		}
		tbl.correlations.present <- (length(rplots) != 0)
		if (tbl.correlations.present) {
			txt <- c("The following table shows the corralation between surrogate variables and sample annotations:")
			rnb.add.paragraph(report, txt)
			report <- rnb.add.figure(report, description, rplots, setting.names)
			colnames(file.tbl) <- c("Target", "Table file")
			rnb.add.exported.tables(report, file.tbl)
		}
		#Traits: p-values
		file.tbl <- matrix(character(), nrow = 0, ncol = 2, dimnames = list(NULL, c("Target", "File Name")))
		rplots <- list()
		for (i in 1:length(targets)){
			tt <- targets[i]
			tn <- names(targets)[i]
			tbl <- sva.obj$assoc$sva[[tt]][["traits"]]$pvalues
			if (!all(is.na(tbl))) {
				#reassign names to make plot layout better
				names.traits <- rownames(tbl)
				tbl <- unname(tbl)
				dimnames(tbl) = list(names.traits,"Surrogate Variable"=1:ncol(tbl))
				fname.plot <- paste0("heatmap_pval_sva_trait_",tn)
				fname.tab  <- paste0("pval_sva_trait_",tn,".csv")
				result <- append.table(tbl, fname.plot, fname.tab, plot.heatmap.pc.pvalues, tab.header.r)
				file.tbl <- rbind(file.tbl, c(targets[tn], result$link))
				rplots[[targets[tn]]] <- result$plot
				description <- result$description
			}
		}
		tbl.pvalues.present <- (length(rplots) != 0)
		if (tbl.pvalues.present) {
			txt <- c("The following table shows p-values resulting from the statistical tests for association of ",
				"surrogate variables and sample annotations:")
			rnb.add.paragraph(report, txt)
			report <- rnb.add.figure(report, description, rplots, setting.names)
			colnames(file.tbl) <- c("Target", "Table file")
			rnb.add.exported.tables(report, file.tbl)
		}

		txt <- c("Surrogate variables were tested for association ",
				 "with principal components using correlation coefficients.")
		report <- rnb.add.section(report, "Association of surrogate variables with principal components", txt, level=2)
		tab.header.r <- "PC"
		#PCA: correlations
		file.tbl <- matrix(character(), nrow = 0, ncol = 2, dimnames = list(NULL, c("Target", "File Name")))
		rplots <- list()
		for (i in 1:length(targets)){
			tt <- targets[i]
			tn <- names(targets)[i]
			tbl <- sva.obj$assoc$sva[[tt]][["pca"]]$correlations
			if (!all(is.na(tbl))) {
				#reassign names to make plot layout better
				tbl <- unname(tbl)
				dimnames(tbl) = list("Principal Component"=1:nrow(tbl),"Surrogate Variable"=1:ncol(tbl))
				fname.plot <- paste0("heatmap_cor_sva_pca_",tn)
				fname.tab  <- paste0("cor_sva_pca_",tn,".csv")
				result <- append.table(tbl, fname.plot, fname.tab, plot.heatmap.pc.correlations, tab.header.r)
				file.tbl <- rbind(file.tbl, c(targets[tn], result$link))
				rplots[[targets[tn]]] <- result$plot
				description <- result$description
			}
		}
		tbl.correlations.present <- (length(rplots) != 0)
		if (tbl.correlations.present) {
			txt <- c("The following table shows the corralation between surrogate variables and principal components:")
			rnb.add.paragraph(report, txt)
			report <- rnb.add.figure(report, description, rplots, setting.names)
			colnames(file.tbl) <- c("Target", "Table file")
			rnb.add.exported.tables(report, file.tbl)
		}
	}
	return(report)
}

#' set.covariates.sva
#'
#' Adds the results of Surrogate Variable Analysis (SVA) to an RnBSet
#' @param rnb.set The \code{RnBSet} object to which the results should be added
#' @param sva.obj An object of class \code{SvaResult} as returned by \code{rnb.execute.sva}.
#' @return The modified \code{RnBSet}. Note that the association information will not be stored.
#' @export
#' @author Fabian Mueller
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' sva.obj <- rnb.execute.sva(rnb.set.example,c("Sample_Group","Treatment"),numSVmethod="be")
#' sva.obj$sva.performed
#' sva.obj$num.components
#' rnb.set.mod <- set.covariates.sva(rnb.set.example, sva.obj)
#' has.covariates.sva(rnb.set.example,"Sample_Group")
#' has.covariates.sva(rnb.set.mod,"Sample_Group")
#' }
set.covariates.sva <- function(rnb.set, sva.obj) {
	sva.obj[["assoc"]] <- NULL # the association tables are not stored in the rnb.set
	rnb.set@inferred.covariates[["sva"]] <- sva.obj
	return(rnb.set)
}

#' get.covariates.sva
#'
#' Retrieves an NxK table of Surrogate variables stored in an RnBSet for a given target variable
#' @param rnb.set \code{RnBSet} object
#' @param target target variable. Must be in \code{pheno(rnb.set)} and belong to target variables for which the
#'		  SVs have already been computed and stored in the RnBSet.
#' @return an NxK table of K Surrogate variables stored for N samples of the \code{rnb.set}. \code{NULL}
#'		   if the components have not been computed or added to \code{rnb.set}.
#' @export
#' @author Fabian Mueller
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' sva.obj <- rnb.execute.sva(rnb.set.example,c("Sample_Group","Treatment"),numSVmethod="be")
#' sva.obj$sva.performed
#' sva.obj$num.components
#' rnb.set.mod <- set.covariates.sva(rnb.set.example, sva.obj)
#' get.covariates.sva(rnb.set.mod,"Sample_Group")
#' }
get.covariates.sva <- function(rnb.set, target) {
	if (!.hasSlot(rnb.set,"inferred.covariates")) { #.hasSlot ensure backwards compatibility
		return(NULL)
	}
	if (!is.element(target,names(rnb.set@inferred.covariates[["sva"]]$components))) {
		return(NULL)
	}
	res <- rnb.set@inferred.covariates[["sva"]]$components[[target]]
	return(res)
}

#' has.covariates.sva
#'
#' Returns whether Surrogate Variables have been computed and added to the \code{rnb.set} for a given target variable
#' @param rnb.set \code{RnBSet} object
#' @param target target variable. Must be in \code{pheno(rnb.set)} and belong to target variables for which the
#'		  SVs have already been computed and stored in the RnBSet.
#' @return \code{logical(1)}
#' @export
#' @author Fabian Mueller
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' sva.obj <- rnb.execute.sva(rnb.set.example,c("Sample_Group","Treatment"),numSVmethod="be")
#' sva.obj$sva.performed
#' sva.obj$num.components
#' rnb.set.mod <- set.covariates.sva(rnb.set.example, sva.obj)
#' has.covariates.sva(rnb.set.example,"Sample_Group")
#' has.covariates.sva(rnb.set.mod,"Sample_Group")
#' has.covariates.sva(rnb.set.mod,"Treatment")
#' }
has.covariates.sva <- function(rnb.set, target) {
	if (!.hasSlot(rnb.set,"inferred.covariates")) { #.hasSlot ensure backwards compatibility
		return(FALSE)
	}
	res <- is.element(target,names(rnb.set@inferred.covariates[["sva"]]$components))
	if (res) res <- ncol(rnb.set@inferred.covariates[["sva"]]$components[[target]]) > 0
	return(res)
}

#' rnb.step.sva
#'
#' Performs Surrogate Variable Analysis (SVA) and adds a corresponding section to an RnBeads Report
#' @param rnb.set The \code{RnBSet} object on which the SVA should be conducted
#' @param report An object of class \code{\linkS4class{Report}}.
#' @param ... arguments passed to \code{rnb.execute.sva}
#' @return A list containing the modified \code{RnBSet} and the modified report
#' @author Fabian Mueller
#' @noRd
rnb.step.sva <- function(rnb.set, report, ...){
	if (!inherits(rnb.set, "RnBSet")) {
		stop("invalid value for rnb.set")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	logger.start("Computing Surrogate Variables and Associations")
	sva.obj <- rnb.execute.sva(rnb.set,...)
	logger.completed()
	logger.start("Adding Report Section")
	report <- rnb.section.sva(report, sva.obj)
	logger.completed()

	#modify the rnb.set
	rnb.set <- set.covariates.sva(rnb.set, sva.obj)
	if (any(sva.obj[["sva.performed"]])){
		logger.info("Added SVA results to the RnBSet")
	}

	#result: modified rnb.set and modified report
	res <- list(rnb.set=rnb.set,report=report)
	return(res)
}
