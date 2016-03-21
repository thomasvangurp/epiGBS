########################################################################################################################
## CNV.R
## created: 2013-06-11
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Copy number variation based on HumanMethylation450 data.
########################################################################################################################

## utility functions

rnb.get.cnv.annotations<-function(platform="probes450"){
	load_env<-new.env()
	load(system.file("extdata/cnv.reference.data.RData", package="RnBeads"), 
			envir=load_env)
	get("cnv.reference.data", envir=load_env)
}

rnb.get.chrom.sizes<-function(assembly="hg19"){
	read.delim(system.file(paste0("extdata/chromSizes/", assembly, ".chrom.sizes"), package="RnBeads"), 
			header=F)
}

#######################################################################################################################

getGLADProfiles<-function(rnb.set, cgs){
	
	if(!inherits(rnb.set,"RnBeadRawSet")){
		stop("RnBeadRawSet object expected")
	}
	
	cnv.reference.data<-rnb.get.cnv.annotations("probes450")
	
	annot<-rnb.annotation2data.frame(rnb.get.annotation("probes450"), add.names=T)
	annot<-annot[rownames(annot) %in% cgs,]
	
	I<-M(rnb.set)+unmethylated(rnb.set)
	I<-I[rownames(annot),]
	
	sample.names<-samples(rnb.set)
	
	primary.data<-lapply(sample.names, function(cl){
				data.frame(
						Clone=rownames(annot),
						PosOrder=1:nrow(annot), 
						LogRatio=log2(I[,cl]/cnv.reference.data$ref.intensity),
						PosBase=annot[,"Start"],
						Chromosome=annot[,"Chromosome"],
						BAC=rownames(annot))
			})
	names(primary.data)<-sample.names
		
	cgh.profiles<-lapply(primary.data, function(df) as.profileCGH(df, value=3))
	names(cgh.profiles)<-sample.names
	
	glad.profiles<-lapply(cgh.profiles, function(cl) daglad(cl,
						#param=0.3,
						#smoothfunc="aws",
						#base=T,
						#smoothfunc="haarseg",
						#genomestep=T,
						#qlambda=0.85,
						mediancenter=TRUE, 
				#lambdabreak=25,
				#bandwidth=1,
				#lambdacluster=2,
				#lambdaclusterGen=2,
				#breaksFdrQ=0.001
				#OnlySmoothing=TRUE
				))
	names(glad.profiles)<-sample.names
	glad.profiles
	
}

#######################################################################################################################

rnb.plot.GLAD.profile<-function(glad.profile, label, sample.names = NA, writeToFile = FALSE, numeric.names = FALSE, ...){
	
	if (!parameter.is.flag(writeToFile)) {
		stop("invalid value for writeToFile; expected TRUE or FALSE")
	}
	cytoband_env<-new.env()
	data(cytoband, package="GLAD", envir = cytoband_env)
	
	if(writeToFile){
		plot.file<-createReportPlot(paste('GLADProfilePlot',  
						ifelse(!is.na(sample.names) && numeric.names, 
								match(sample, sample.names), gsub("[ |_]", ".", sample)) , sep="_"), ...)
	}
		
	GLAD::plotProfile(glad.profile,cytoband=get("cytoband", envir = cytoband_env), Bkp=T, Smoothing="Smoothing", main=paste(label))
		
	if(writeToFile){
		off(plot.file)
	}
	
	if(writeToFile){
		return(plot.file)
	}
			
}

#######################################################################################################################

getCGCounts<-function(cnv.profiles, rnb.set){
	
	cnv.reference.data<-rnb.get.cnv.annotations("probes450")
	chrom.data<-rnb.get.chrom.sizes("hg19")
	
	si<-cnv.reference.data$seq.info.hg19
	annot<-annotation(rnb.set)
	
	results<-sapply(cnv.profiles, function(profile){
		
		chr.basic<-GRanges(seqnames=ChrNumeric(chrom.data[,1]), 
				IRanges(start=rep(1,nrow(chrom.data)), end=chrom.data[,2]), 
				strand=(rep("*", nrow(chrom.data))), regid=rep("native", nrow(chrom.data)),
				seqinfo=si)		
				
		sample.cgh<-as.data.frame(profile)
		sample.bkps<-profile$BkpInfo
					
		sample.ranges<-GRanges(seqnames=sample.bkps$Chromosome, 
				IRanges(start=sample.bkps$PosBase,width=1), strand=rep("*", nrow(sample.bkps)),
				seqinfo=si, regid=sample.bkps[,1])
		sample.ranges<-c(chr.basic, sample.ranges)
		sample.ranges<-disjoin(sample.ranges)
		start(sample.ranges[which(width(sample.ranges)==1)+1L])<-start(sample.ranges[which(width(sample.ranges)==1)+1L])-1
		sample.ranges<-sample.ranges[width(sample.ranges)>1]
					
		sample.prof.ranges<-GRanges(seqnames=ChrNumeric(sample.cgh$Chromosome), 
				IRanges(start=sample.cgh$PosBase, width=1), strand="*", seqinfo=si)
		sample.olaps<-findOverlaps(sample.prof.ranges, sample.ranges)
		
		mean.zonegnl<-tapply(sample.cgh[queryHits(sample.olaps),"ZoneGNL"], subjectHits(sample.olaps), mean)
		Level<-round(mean.zonegnl)+2L # expecting 2n, can be very wrong if the sample is on average kn, k>2
		sample.counts<-cbind(as.data.frame(sample.ranges),Level)	
		
		sample.counts$chromosome<-rnb.get.chromosomes()[sample.counts$seqnames]
		sample.counts$seqnames<-NULL
		sample.counts<-sample.counts[,c(6,1,2,4,5)]			
		
		cg.ranges<-GRanges(seqnames=ChrNumeric(annot$Chromosome), 
				IRanges(start=annot$Start, width=1), strand="*", seqinfo=si)
		cg.olaps<-findOverlaps(cg.ranges, sample.ranges)
		#cpg.gnl<-tapply(sample.cgh[queryHits(sample.olaps),"ZoneGNL"], 
		#	subjectHits(sample.olaps), function(cgs) rep(subjectHits(sample.olaps)[queryHits]) )
		
		cnv.vector<-sample.counts[subjectHits(cg.olaps),"Level"]
		#names(cnv.vector)<-rownames(annot)[queryHits(cg.olaps)]
			
	})
	
	colnames(results)<-names(cnv.profiles)
	rownames(results)<-rownames(annot)
	results
	
}

#######################################################################################################################
#' rnb.execute.cnv
#'
#' Copy number variation calling using GLAD
#'
#' @param object \code{\linkS4class{RnBeadRawSet}} object
#' 
#' @return a \code{list} with elements \code{glad.profiles} and \code{cg.counts}
#'
#' @author Pavlo Lutsik
#' @noRd

rnb.execute.cnv<-function(object){
	
	cnv.reference.data<-rnb.get.cnv.annotations()
	
	glad.profiles <- getGLADProfiles(object, cnv.reference.data$bac.array.cgs)
	cg.counts <- getCGCounts(glad.profiles, object)
	
	list(glad.profiles=glad.profiles, cg.counts=cg.counts)
	
}

#######################################################################################################################

#' rnb.section.cnv
#'
#' Adds CNV section to quality control report
#'
#' @param report analysis report
#' @param cnv.data a \code{list} output by \code{rnb.execute.cnv}
#' @param cnv.profile.plot add genome-wide plots of 
#' @param qc.snpbarplot add bar plots of beta-values at SNP probes, one bar plot par sample
#' 
#' @return report object
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.section.cnv<-function(report, cnv.data){
	
	logger.start("Copy Number Variation Section")
	
	report <- rnb.add.section(report, "Copy number variation analysis", 
			"Visualization and analysis of the copy number variations based on HumanMethylation450 data")
	
	#ids<-annotation(object, add.names=TRUE)[,"ID"]
	
	if(rnb.getOption("cnv.profile.plot")){
		
		txt <- "CNV profile plot visualizes the results of CNV analysis using the GLAD package"
		rnb.add.paragraph(report, c("<b>CNV profile</b>. ", txt))
		report<-add.profile.plots(report, cnv.data$cnv.profiles)
		logger.status("Added CNV profile plots")
		
	}
	
	fname <- "cnv_per_cpg.csv"
	fname.full <- file.path(rnb.get.directory(report, "data", absolute = TRUE), fname)
	utils::write.csv(cnv.data$cg.counts, file = fname.full, na = "", row.names = FALSE)
	
	txt<-c("Based on the segments with the identical copy number, one can calculate copy number for each CpG. ",
			"The table of computed copy numbers is available as a <a href=\"",
			rnb.get.directory(report, "data"), "/", fname, "\">comma-separated file</a> accompanying this report.")
	rnb.add.paragraph(report, c("<b>Copy gains per CpG</b>. ", txt))
	
	logger.status("Added CpG counts")
	
	logger.completed()
	return(report)
	
}
#######################################################################################################################

#' rnb.step.cnv
#'
#' Performs copy nubmer calling from the HumanMethylation450 intenstity data
#'
#' @param intensities	\code{list} with methylated and unmethylated probe intensities
#' @param report  		Report on quality control to contain the generated sections. This must be an object of type
#'                		\code{\linkS4class{Report}}.
#' @return The modified report.
#'
#' @author Pavlo Lutsik
#' @noRd
rnb.step.cnv<-function(rnb.set, report){
	
	suppressPackageStartupMessages(require(GLAD))
	
	if (!is.list(intensities)){
		stop("invalid value for intensities")
	}
	if (!inherits(report, "Report")) {
		stop("invalid value for report")
	}
	
	logger.start("Copy Number Variation")
	
	cnv.data<-rnb.execute.cnv(rnb.set)
	
	report<-rnb.section.cnv(report, cnv.data)
	
	logger.completed()
	
	report
	
}

#######################################################################################################################
#######################################################################################################################

add.profile.plots<-function(report, cnv.profiles){
	
	descr<-"Profiles visualize copy number variation across the genome of HumanMethylation450 profiled samples."
	
	ids<-names(cnv.profiles)
	
	cplots<-lapply(ids, function(id) {
				rnb.plot.GLAD.profile(glad.profile=cnv.profiles[[id]], label=id, sample.names=ids, 
						report=report, writeToFile=T, numeric.names=T, width=8, height=7, low.png=100, high.png=300)
			})
	
	names(cplots)<-1:length(ids)
	
	sn<-list("Sample labels" = ids)
	names(sn[[1]])<-1:length(ids)
	
	report<-rnb.add.figure(report, description=descr, report.plots=cplots, setting.names=sn)
	report
}
#######################################################################################################################