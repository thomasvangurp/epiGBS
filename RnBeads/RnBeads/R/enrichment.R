########################################################################################################################
## differentialMethylation.R
## created: 2012-09-03
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Methods for enrichment analysis.
########################################################################################################################

#' performGOenrichment.diffMeth.entrez
#'
#' performs Gene Ontology (GO) enrichment analysis for a list of Entrez identifiers
#' @param gids gene ids to test (entrez IDs)
#' @param uids ids to test against (universe)
#' @param ontology which ontology should be used (see \code{GOHyperGParams} from the \code{GOstats} package for details)
#' @param assembly Genome to be used. One of the following: hg19, mm9, mm10 or rn5
#' @param ... arguments passed on to the parameters of \code{GOHyperGParams} from the \code{GOstats} package
#' @return a \code{GOHyperGresult} object (see the \code{GOstats} package for furthert details)
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' dmt <- get.table(dm,get.comparisons(dm)[1],"promoters")
#' annot <- annotation(rnb.set.example,"promoters")
#' all.promoters <- annot$entrezID
#' #get the hypermethylated promoters
#' hyper.promoters <- annot$entrezID[dmt[,"mean.mean.diff"]>0]
#' result <- performGOenrichment.diffMeth.entrez(hyper.promoters,all.promoters,"BP",assembly="hg19")
#' }
performGOenrichment.diffMeth.entrez <- function(gids,uids,ontology,assembly="hg19",...){
	suppressPackageStartupMessages(require(GOstats))
#	#testing
#	dmt <- diffmeth$region[[1]][["promoters"]]
#	annot <- rnb.get.annotations(type = "promoters")
#	uids <- unique(unlist(sapply(annot[rownames(dmt),]$entrezid,FUN=function(idstr){
#		unlist(strsplit(idstr, ";", fixed = TRUE))
#	})))
#	dmrs <- dmt[dmt[,"combinedRank"] < 1000,]
#	gids <- unique(unlist(sapply(annot[rownames(dmrs),]$entrezid,FUN=function(idstr){
#		unlist(strsplit(idstr, ";", fixed = TRUE))
#	})))
#	#/testing
	if (assembly == "hg19"){
		ass <- "org.Hs.eg.db"
	} else if (is.element(assembly,c("mm9","mm10"))){
		ass <- "org.Mm.eg.db"
	} else if (is.element(assembly,c("rn5"))){
		ass <- "org.Rn.eg.db"
	} else {
		stop("Unsupported assembly")
	}
	params <- new("GOHyperGParams",annotation=ass,geneIds = gids, universeGeneIds = uids, 
			ontology = ontology,conditional = TRUE, testDirection = "over",...)
	hgOver <- tryCatch({
			GOstats::hyperGTest(params)
		}, error = function(ee) {
			if (ee$message=="The genes you are testing do not have any corresponding GO terms for the ontology you are searching."){
				logger.info("Could not conduct enrichment analysis as associated genes are not in GO database.")
				NULL
			} else {
				stop(ee)
			}
		}
	)
	return(hgOver)
}

#' performEnrichment.diffMeth
#'
#' performs Geno Ontology (GO) enrichment analysis for a given differential methylation table.
#' @author Fabian Mueller
#' @param rnbSet RnBSet object for which dirrential methylation was computed
#' @param diffmeth RnBDiffMeth object. See \code{\link{RnBDiffMeth-class}} for details.
#' @param ontologies GO ontologies to use for enrichment analysis
#' @param rank.cuts.region Cutoffs for combined ranking that are used to determine differentially methylated regions
#' @param add.auto.rank.cut flag indicating whether an automatically computed cut-off should also be considered.
#' @param rerank For deterimining differential methylation: should the ranks be ranked again or should the absolute ranks be used.
#' @param verbose Enable for detailed status report
#' @param ... arguments passed on to the parameters of \code{GOHyperGParams} from the \code{GOstats} package
#' @return a DiffMeth.enrich object (S3) containing the following attributes
#' \item{region}{Enrichment information for differential methylation on the region level. See \code{GOHyperGresult} from the \code{GOstats} package for furthert details}
#' @export
#' @examples
#' \dontrun{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' res <- performEnrichment.diffMeth(rnb.set.example,dm)
#' }
performEnrichment.diffMeth <- function(rnbSet,diffmeth,ontologies=c("BP","MF"),rank.cuts.region=c(100,500,1000),add.auto.rank.cut=TRUE,rerank=TRUE,verbose=TRUE,...){
	gene.id.col <- "entrezID"
	logger.start("Differential Methylation Enrichment Analysis")
	if (!is.element(assembly(rnbSet),c("hg19","mm9","mm10","rn5"))){
		logger.warning(c("Enrichment analysis currently not supported for genome assembly:",assembly(rnbSet),"--> skipping enrichment analysis"))
		logger.completed()
		return(NULL)
	}
	comps <- get.comparisons(diffmeth)
	region.types <- get.region.types(diffmeth)
	dm.enrich <- list(probe=list(),region=list())
	for (cc in comps){
		logger.start(c("Comparison: ",cc))
		dm.enrich$probe[[cc]] <- list()
		dm.enrich$region[[cc]] <- list()
		for (oo in ontologies){
			logger.start(c("Ontology: ",oo))
			logger.start("Region Level")
			dm.en.region.list.list <- list()
			for (rt in region.types){
				logger.start(c("Region Type:",rt))
				#rt <- "promoters"
				annot <- annotation(rnbSet,type = rt)
				if (!(gene.id.col %in% colnames(annot))){
					logger.info(c("Not annotated with",gene.id.col,"--> Skipped"))
					logger.completed()
					next
				}
#				dmt <- diffmeth$region[[cc]][[rt]]
				dmt <- get.table(diffmeth,cc,rt,undump=TRUE,return.data.frame=TRUE)
				rrs <- dmt[,"combinedRank"]
				rrs.hyper <- rrs
				rrs.hypo <- rrs
				rrs.hyper[dmt[,"mean.mean.diff"] <= 0] <- NA
				rrs.hypo[dmt[,"mean.mean.diff"] >= 0] <- NA
				#automatically select rank cutoff
				rc.auto <- 0L
				if (add.auto.rank.cut){
					rc.auto <- as.integer(auto.select.rank.cut(dmt$comb.p.adj.fdr,dmt$combinedRank))
				}

				if (rerank){
					rrs <- rank(rrs,na.last="keep",ties.method="min")
					rrs.hyper <- rank(rrs.hyper,na.last="keep",ties.method="min")
					rrs.hypo <- rank(rrs.hypo,na.last="keep",ties.method="min")
					if (add.auto.rank.cut && rc.auto > 0L){
						rc.auto <- na.omit(rrs[dmt[,"combinedRank"]==rc.auto])[1] #arbitrarily select the first rerank if the rank cut threshold is occupied by multiple ranks
					}
				}
				
				gene.ids.chr <- as.character(annot[,gene.id.col])
				uids <- unique(unlist(sapply(gene.ids.chr,FUN=function(idstr){
					unlist(strsplit(idstr, ";", fixed = TRUE))
				})))
				uids <- na.omit(uids)
				dm.en.region.list <- list()
				rank.cuts <- rank.cuts.region
				if (add.auto.rank.cut){
					rank.cuts <- c(rank.cuts,"autoSelect")
				}
				for (i in 1:length(rank.cuts)){
					rc <- rank.cuts[i]
					if (rc == "autoSelect") {
						rc <- rc.auto
						if (verbose) logger.info(c("Rank cutoff:",rc,"(auto-select)"))
					} else {
						rc <- as.integer(rc)
						if (verbose) logger.info(c("Rank cutoff:",rc))
					}
					is.dmr.hyper <- rrs.hyper <= rc
					is.dmr.hypo <- rrs.hypo <= rc
					gids.hyper <- unique(unlist(sapply(gene.ids.chr[is.dmr.hyper],FUN=function(idstr){
						unlist(strsplit(idstr, ";", fixed = TRUE))
					})))
					gids.hyper <- na.omit(gids.hyper)
					if (length(gids.hyper)>0 && length(uids)>0){
						dm.en.region.hyper <- performGOenrichment.diffMeth.entrez(gids.hyper,uids,ontology=oo,assembly=assembly(rnbSet),...)
					} else {
						dm.en.region.hyper <- NULL
					}
					
					gids.hypo <- unique(unlist(sapply(gene.ids.chr[is.dmr.hypo],FUN=function(idstr){
						unlist(strsplit(idstr, ";", fixed = TRUE))
					})))
					gids.hypo <- na.omit(gids.hypo)
					if (length(gids.hypo)>0 && length(uids)>0){
						dm.en.region.hypo <- performGOenrichment.diffMeth.entrez(gids.hypo,uids,ontology=oo,assembly=assembly(rnbSet),...)
					} else {
						dm.en.region.hypo <- NULL
					}
					dm.en.region.list <- c(dm.en.region.list,list(list(hyper=dm.en.region.hyper,hypo=dm.en.region.hypo)))
				}
				names(dm.en.region.list) <- paste("rankCut_",rank.cuts,sep="")
				dm.en.region.list.list <- c(dm.en.region.list.list,list(dm.en.region.list))
				names(dm.en.region.list.list)[length(dm.en.region.list.list)] <- rt
				logger.completed()
			}
			dm.enrich$region[[cc]][[oo]] <- dm.en.region.list.list
			logger.completed()
			logger.completed()
		}
		logger.completed()
	}
	class(dm.enrich) <- "DiffMeth.enrich"
	logger.completed()
	return(dm.enrich)
}
