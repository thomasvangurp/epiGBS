################################################################################
# RnBClusterRun
################################################################################

#' RnBClusterRun Class
#'
#' A class for configuring and running RnBeads on a scientific compute cluster.
#' 
#'
#' @section Slots:
#' \describe{
#'   \item{\code{architecture}}{A \code{\linkS4class{ClusterArchitecture}} object managing the settings for a scientific compute cluster}
#'   \item{\code{modules}}{A vector of pipeline modules} 
#'   \item{\code{module.res.req}}{Stores the resource requirements for each module. A list containing named vectors for the resources} 
#'   \item{\code{module.num.cores}}{Stores the number of cores for each module}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{setModuleResourceRequirements,RnBClusterRun,character,character-method}}}{Sets the resource requirements for the different pipeline modules}
#'   \item{\code{\link{setModuleNumCores,RnBClusterRun,integer,character-method}}}{Sets the number of cores used by the different pipeline modules}
#'   \item{\code{\link{getModuleNumCores,RnBClusterRun-method}}}{Gets the number of cores used by the different pipeline modules}
#'   \item{\code{\link{run,RnBClusterRun-method}}}{Submit the pipeline modules to the cluster}
#' }
#'
#' @name RnBClusterRun-class
#' @rdname RnBClusterRun-class
#' @author Fabian Mueller
#' @exportClass RnBClusterRun
setClass("RnBClusterRun",
	slots = list(
		architecture = "ClusterArchitecture",
		modules = "character",
		module.res.req = "list",
		module.num.cores = "integer"
	)
)

#' initialize.RnBClusterRun
#'
#' Initialize an RnBClusterRun object
#' 
#' @param architecture A \code{\linkS4class{ClusterArchitecture}} object managing the settings for a scientific compute cluster.
#'
#' @author Fabian Mueller
#' @docType methods
setMethod("initialize","RnBClusterRun",
	function(
		.Object,
		architecture
	) {
		.Object@architecture <- architecture
		.Object@modules <- c (
			import        = "data import",
			qc            = "quality control",
			preprocessing = "preprocessing",
			tnt           = "tracks and tables",
			inference     = "covariate inference",
			exploratory   = "exploratory analysis",
			differential  = "differential methylation",
			differential_chunk   = "differential methylation (chunk)",
			differential_wrapup  = "differential methylation (wrapup)",
			wrapup        = "wrapup"
		)
		module.res.req <- rep(list(character(0)),length(.Object@modules))
		names(module.res.req) <- names(.Object@modules)

		module.num.cores <- rep(1L,length(.Object@modules))
		names(module.num.cores) <- names(.Object@modules)

		.Object
	}
)


if (!isGeneric("setModuleResourceRequirements")) setGeneric("setModuleResourceRequirements", function(object,resources,modules) standardGeneric("setModuleResourceRequirements"))
#' setModuleResourceRequirements-methods
#'
#' Specifies resource requirements for the different pipeline modules
#'
#' @param object \code{\linkS4class{RnBClusterRun}} object
#' @param resources A NAMED character vector containing the resource reuirements as value and the resource name as name
#' @param modules vector of applicable pipeline modules. Can be \code{"all"} to specify all modules
#' @return The modified object
#'
#' @rdname setModuleResourceRequirements-RnBClusterRun-methods
#' @docType methods
#' @aliases setModuleResourceRequirements
#' @aliases setModuleResourceRequirements,RnBClusterRun-method
#' @author Fabian Mueller
#' @export
setMethod("setModuleResourceRequirements",
	signature(
		object="RnBClusterRun",
		resources="character",
		modules="character"
	),
	function(
		object,
		resources,
		modules="all"
	) {
		if (length(modules)<1 || !all(modules %in% c(names(object@modules),"all"))){
			stop("Invalid modules specified")
		}
		if (length(resources) < 1){
			stop("Invalid resources parameter. Has to have at least length 1.")
		}
		if (any(is.na(names(resources))) || is.null(names(resources))){
			stop("Invalid resource requirement specification. Need names for all requirements")
		}
		if (any(modules=="all")){
			modules = names(object@modules)
		}
		for (mm in modules){
			object@module.res.req[[mm]][names(resources)] <- resources
		}
		return(object)
	}
)

if (!isGeneric("setModuleNumCores")) setGeneric("setModuleNumCores", function(object,num.cores,modules) standardGeneric("setModuleNumCores"))
#' setModuleNumCores-methods
#'
#' Specifies the number of cores used by the different pipeline modules
#'
#' @param object \code{\linkS4class{RnBClusterRun}} object
#' @param num.cores an integer specifying the number of cores to be used
#' @param modules vector of applicable pipeline modules. Can be \code{"all"} to specify all modules
#' @return The modified object
#'
#' @rdname setModuleNumCores-RnBClusterRun-methods
#' @docType methods
#' @aliases setModuleNumCores
#' @aliases setModuleNumCores,RnBClusterRun-method
#' @author Fabian Mueller
#' @export
setMethod("setModuleNumCores",
	signature(
		object="RnBClusterRun",
		num.cores="integer",
		modules="character"
	),
	function(
		object,
		num.cores,
		modules="all"
	) {
		if (length(modules)<1 || !all(modules %in% c(names(object@modules),"all"))){
			stop("Invalid modules specified")
		}
		if (any(num.cores<1)){
			stop("invalid number of cores specified. Must be 1 or more.")
		}
		num.cores <- num.cores[1]
		if (any(modules=="all")){
			modules = names(object@modules)
		}
		for (mm in modules){
			object@module.num.cores[mm] <- num.cores
		}
		return(object)
	}
)

if (!isGeneric("getModuleNumCores")) setGeneric("getModuleNumCores", function(object) standardGeneric("getModuleNumCores"))
#' getModuleNumCores-methods
#'
#' Retrieves the number of cores used by each module
#'
#' @param object \code{\linkS4class{RnBClusterRun}} object
#' @return A named vector containing the number of cores for each module
#'
#' @rdname getModuleNumCores-RnBClusterRun-methods
#' @docType methods
#' @aliases getModuleNumCores
#' @aliases getModuleNumCores,RnBClusterRun-method
#' @author Fabian Mueller
#' @export
setMethod("getModuleNumCores",
	signature(
		object="RnBClusterRun"
	),
	function(
		object
	) {
		return(object@module.num.cores)
	}
)

if (!isGeneric("run")) setGeneric("run", function(rnb.cr,...) standardGeneric("run"))
#' run-methods
#'
#' Runs the analysis by submitting jobs for each module to the compute cluster
#'
#' @param rnb.cr \code{\linkS4class{RnBClusterRun}} object
#' @param analysis.id analysis id. used for naming submitted jobs and log files
#' @param config.xml XML file specifying the analysis options and parameter settings
#' @param split.differential flag indicating whether to split the differnetial methylation module
#'        into seperate jobs according to sample annotation column and region type.
#' @param dry.run Prevent the actual job submission. Rather only write to a shell script file
#' @param long.cmd.thres commands that are longer than this number will be encapsulated in shell scripts
#' 		  rather than being submitted as direct command
#' @return Nothing of importance
#'
#' @rdname run-RnBClusterRun-methods
#' @docType methods
#' @aliases run
#' @aliases run,RnBClusterRun-method
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' #specify the xml file for your analysis
#' xml.file <- "MY_ANALYSIS_SETTINGS.XML"
#' #set the cluster architecture specific to your environment
#' arch <- new("ClusterArchitectureSGE")
#' rnb.cr <- new("RnBClusterRun",arch)
#' #set up the cluster so that 32GB of memory are required (SGE resource is called "mem_free")
#' rnb.cr <- setModuleResourceRequirements(rnb.cr,c(mem_free="32G"),"all")
#' #set up the cluster to use 4 cores on each node for all modules
#' rnb.cr <- setModuleNumCores(rnb.cr,4L,"all")
#' #set up the cluster to use 2 cores for the exploratory analysis module
#' rnb.cr <- setModuleNumCores(rnb.cr,2L,"exploratory")
#' #run the actual analysis (remove dry.run=TRUE, to really submit the jobs)
#' run(rnb.cr, "rnbeads_analysis", xml.file, dry.run=TRUE)
#' }
setMethod("run",
	signature(
		rnb.cr="RnBClusterRun"
	),
	function(
		rnb.cr,
		analysis.id,
		config.xml,
		split.differential=TRUE,
		dry.run=FALSE,
		long.cmd.thres=1024L
	) {
		arch = rnb.cr@architecture
		r.exec <- getExecutable(arch,"Rscript")

		rnb.settings <- rnb.xml2options(config.xml,return.full.structure=TRUE)
		if (!("dir.reports" %in% names(rnb.settings$analysis.params))) {
			stop("Invalid analysis parameters: 'dir.reports' must be specified")
		}
		dir.reports <- rnb.settings$analysis.params[["dir.reports"]]
		data.type   <- rnb.settings$analysis.params[["data.type"]]
		if (!rnb.initialize.reports(dir.reports)) {
			stop(paste("Could not initialize reports in", dir.reports, "; make sure this path does not exist."))
		}
		cluster.dir <- file.path(dir.reports,"cluster_run")
		create.path(cluster.dir)
		log.dir <- cluster.dir #file.path(cluster.dir,"logs")
		#create.path(log.dir)
		logger.start("RnBeads Cluster Run",fname=c(file.path(log.dir,"cluster_run.log"),NA))

		if ("preanalysis.script" %in% names(rnb.settings)){
			source(rnb.settings$preanalysis.script)
		}

		do.call(rnb.options, rnb.settings$options)
		analysis.options <- rnb.options()

		## Set options
		if (length(analysis.options) != 0) {
			do.call(rnb.options, analysis.options)
		}

		logger.start("Saving settings")
		settings.file <- file.path(cluster.dir,"clusterRunSettings.RData")
		save(rnb.cr,rnb.settings,analysis.options,file=settings.file)
		file.copy(config.xml, file.path(cluster.dir,"options.xml"))
		logger.completed()

		#Can we submit shell scripts instead of binary commands for long commands?
		shell.script.for.long.commands <- is.element("sub.binary",arch@getSubCmdTokens.optional.args) && long.cmd.thres > 0L
		shell.script.dir <- cluster.dir
		submit.job <- function(name,cmd.tokens,...){
			r.cmd <- paste(cmd.tokens,collapse=" ")
			cmd <- getSubCmdStr(arch, cmd.tokens, ...)
			#make sure the command is not too long. else, wrap it in a shell script
			if (shell.script.for.long.commands){
				if (nchar(r.cmd)>long.cmd.thres){
					shell.script.file <- file.path(shell.script.dir,paste0(name,".sh"))
					fileConn<-file(shell.script.file)
					writeLines(c("#!/bin/sh",r.cmd), fileConn)
					close(fileConn)
					Sys.chmod(shell.script.file, mode = "0755")
					cmd.tokens.shell <- c(shell.script.file)
					cmd <- getSubCmdStr(arch, cmd.tokens.shell, sub.binary=FALSE, ...)
				}
			}

			#actually submit
			if (!dry.run){
				system(cmd)
			}
			logger.info(c("Command used:",cmd))
			return(cmd)
		}

		cmds.submit <- c()
		deps.wrapup <- c()

		mm <- "import"
		logger.start(c("Running:",mm))
			log.file <- file.path(log.dir,paste0(mm,".log"))
			jid <- paste(analysis.id,mm,sep="_")
			res.req <- rnb.cr@module.res.req[[mm]]
			script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")

			cmd.tokens <- c(r.exec,script.file, "-x",config.xml, "-o",cluster.dir, "-c",getModuleNumCores(rnb.cr)[mm])
			cmds.submit[mm] <- submit.job(mm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req)
			deps.wrapup <- c(deps.wrapup,jid)
		logger.completed()
		jid.import <- jid

		rnb.set.file <- file.path(cluster.dir,paste0("import","_RnBSet"))
		depend.jobs <- jid.import

		mm <- "preprocessing"
		logger.start(c("Running:",mm))
			log.file <- file.path(log.dir,paste0(mm,".log"))
			jid <- paste(analysis.id,mm,sep="_")
			res.req <- rnb.cr@module.res.req[[mm]]
			script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")
			
			cmd.tokens <- c(
				r.exec,script.file,
				"-x",config.xml,
				"-s",rnb.set.file,
				"-o",cluster.dir,
				"-c",getModuleNumCores(rnb.cr)[mm]
			)
			cmds.submit[mm] <- submit.job(mm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req, depend.jobs=depend.jobs)
			deps.wrapup <- c(deps.wrapup,jid)
		logger.completed()
		jid.preprocessing <- jid

		if (rnb.getOption("qc")){
			mm <- "qc"
			logger.start(c("Running:",mm))
				log.file <- file.path(log.dir,paste0(mm,".log"))
				jid <- paste(analysis.id,mm,sep="_")
				res.req <- rnb.cr@module.res.req[[mm]]
				script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")
				cmd.tokens <- c(
					r.exec,script.file,
					"-x",config.xml,
					"-s",rnb.set.file,
					"-c",getModuleNumCores(rnb.cr)[mm]
				)
				cmds.submit[mm] <- submit.job(mm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req, depend.jobs=depend.jobs)
				deps.wrapup <- c(deps.wrapup,jid)
			logger.completed()
		}

		rnb.set.file <- file.path(cluster.dir,paste0("preprocessing","_RnBSet"))
		depend.jobs <- jid.preprocessing

		if (rnb.getOption("export.to.bed") || (length(rnb.getOption("export.to.trackhub")) > 0)){
			mm <- "tnt"
			logger.start(c("Running:","Tracks and Tables (tnt)"))
				log.file <- file.path(log.dir,paste0(mm,".log"))
				jid <- paste(analysis.id,mm,sep="_")
				res.req <- rnb.cr@module.res.req[[mm]]
				script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")
				
				cmd.tokens <- c(
					r.exec,script.file,
					"-x",config.xml,
					"-s",rnb.set.file,
					"-c",getModuleNumCores(rnb.cr)[mm]
				)
				cmds.submit[mm] <- submit.job(mm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req, depend.jobs=depend.jobs)
				deps.wrapup <- c(deps.wrapup,jid)
			logger.completed()
		}

		if (rnb.getOption("inference")){
			mm <- "inference"
			logger.start(c("Running:",mm))
				log.file <- file.path(log.dir,paste0(mm,".log"))
				jid <- paste(analysis.id,mm,sep="_")
				res.req <- rnb.cr@module.res.req[[mm]]
				script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")
				cmd.tokens <- c(
					r.exec,script.file,
					"-x",config.xml,
					"-s",rnb.set.file,
					"-o",cluster.dir,
					"-c",getModuleNumCores(rnb.cr)[mm]
				)
				cmds.submit[mm] <- submit.job(mm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req, depend.jobs=depend.jobs)
				deps.wrapup <- c(deps.wrapup,jid)
			logger.completed()

			jid.inference <- jid
			rnb.set.file <- file.path(cluster.dir,paste0("inference","_RnBSet"))
			depend.jobs <- jid.inference
		}

		if (rnb.getOption("exploratory")){
			mm <- "exploratory"
			logger.start(c("Running:",mm))
				log.file <- file.path(log.dir,paste0(mm,".log"))
				jid <- paste(analysis.id,mm,sep="_")
				res.req <- rnb.cr@module.res.req[[mm]]
				script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")
				cmd.tokens <- c(
					r.exec,script.file,
					"-x",config.xml,
					"-s",rnb.set.file,
					"-c",getModuleNumCores(rnb.cr)[mm]
				)
				cmds.submit[mm] <- submit.job(mm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req, depend.jobs=depend.jobs)
				deps.wrapup <- c(deps.wrapup,jid)
			logger.completed()
		}

		if (rnb.getOption("differential")){
			logger.start(c("Running:","differential"))
				cmp.cols <- rnb.getOption("differential.comparison.columns")
				if ((rnb.getOption("differential.permutations") != 0L) && split.differential){
					logger.warning(c("Option 'differential.permutations' is currently not compatible with splitting jobs for differential methylation on a cluster",
								"--> proceed without splitting the differential methylation job"))
					split.differential <- FALSE
				}
				if (is.null(cmp.cols) && split.differential){
					logger.warning(c("differential.comparison.columns=NULL is currently not compatible with splitting differential methylation task into chunks","--> proceed without splitting the differential methylation job"))
					split.differential <- FALSE
				}
				if (split.differential){
					logger.info(c("Splitting up into chunks"))
					logger.start(c("Running chunks"))
						mm <- "differential_chunk"
						res.req <- rnb.cr@module.res.req[[mm]]
						script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")

						reg.types <- rnb.region.types.for.analysis(rnb.getOption("assembly"))
						if(rnb.getOption("region.subsegments")>1L){
							reg.types <- c(reg.types,paste(reg.types,"subsegments",sep="."))
						}
						if (length(reg.types)>0){
							logger.info(c("Using region types:",paste(reg.types,collapse=",")))
						} else {
							logger.info(c("Skipping region level analysis"))
						}
						#split by pheno column and region type
						do.regions <- length(reg.types)>0
						combinations <- data.frame(Var1=cmp.cols,Var2=rep(NA,length(cmp.cols)))
						if (do.regions){
							combinations <- expand.grid(cmp.cols,reg.types)
						}

						#submit the jobs for each chunk
						chunk.ids <- c()
						depend.jobs.wrapup <- c()
						for (i in 1:nrow(combinations)){
							pp <- cmp.cols[combinations[i,1]]
							rr <- reg.types[combinations[i,2]]
							chunk.id.cur <-  paste0(pp,"_")
							if (do.regions){
								chunk.id.cur <- paste(pp,rr,sep="_")
							}

							mmm <- paste0(mm,"_",chunk.id.cur)
							log.file <- file.path(log.dir,paste0(mmm,".log"))
							jid <- paste(analysis.id,mmm,sep="_")
							chunk.ids <- c(chunk.ids,chunk.id.cur)
							depend.jobs.wrapup <- c(depend.jobs.wrapup,jid)

							r.tokens <- NULL
							if (do.regions){
								r.tokens <- c("-r",rr)
							}
							cmd.tokens <- c(
								r.exec,script.file,
								"-x",config.xml,
								"-s",rnb.set.file,
								"-i",chunk.id.cur,
								"-o",cluster.dir,
								"-p",pp,
								r.tokens,
								"-c",getModuleNumCores(rnb.cr)[mm]
							)
							cmds.submit[mmm] <- submit.job(mmm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req, depend.jobs=depend.jobs)
						}
					logger.completed()
					logger.start(c("Running wrapup"))
						mm <- "differential_wrapup"
						log.file <- file.path(log.dir,paste0(mm,".log"))
						jid <- paste(analysis.id,mm,sep="_")
						res.req <- rnb.cr@module.res.req[[mm]]
						script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")

						chunk.tokens <- as.vector(rbind(rep("-d",length(chunk.ids)),chunk.ids))

						cmd.tokens <- c(
							r.exec,script.file,
							"-x",config.xml,
							"-s",rnb.set.file,
							"-o",cluster.dir,
							chunk.tokens,
							"-c",getModuleNumCores(rnb.cr)[mm]
						)
						cmds.submit[mm] <- submit.job(mm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req, depend.jobs=depend.jobs.wrapup)
						deps.wrapup <- c(deps.wrapup,jid)
					logger.completed()

				} else {
					mm <- "differential"
					log.file <- file.path(log.dir,paste0(mm,".log"))
					jid <- paste(analysis.id,mm,sep="_")
					res.req <- rnb.cr@module.res.req[[mm]]
					script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")
					cmd.tokens <- c(
						r.exec,script.file,
						"-x",config.xml,
						"-s",rnb.set.file,
						"-o",cluster.dir,
						"-c",getModuleNumCores(rnb.cr)[mm]
					)
					cmds.submit[mm] <- submit.job(mm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req, depend.jobs=depend.jobs)
					deps.wrapup <- c(deps.wrapup,jid)
				}
			logger.completed()
		}

		mm <- "wrapup"
		logger.start(c("Running:",mm))
			log.file <- file.path(log.dir,paste0(mm,".log"))
			jid <- paste(analysis.id,mm,sep="_")
			res.req <- rnb.cr@module.res.req[[mm]]
			script.file <- system.file(file.path("extdata","Rscript",paste0("rscript_",mm,".R")), package = "RnBeads")
			cmd.tokens <- c(
				r.exec,script.file,
				"-x",config.xml,
				"-c",getModuleNumCores(rnb.cr)[mm]
			)
			cmds.submit[mm] <- submit.job(mm, cmd.tokens, log=log.file, job.name=jid, res.req=res.req, depend.jobs=deps.wrapup)
		logger.completed()

		#save the submission commands to a shell script
		cmd.shellscript.file <- file(file.path(cluster.dir,"submitted_jobs.sh"))
		writeLines(c("#!/bin/sh",cmds.submit), cmd.shellscript.file)
		close(cmd.shellscript.file)

		logger.completed()
		invisible(TRUE)
	}
)

################################################################################
# Helper and utility functions
################################################################################
#' logger.machine.name
#'
#' log the machine name the analysis is run on
#' @author Fabian Mueller
#' @export
logger.machine.name <- function(){
	uname <- Sys.info()[["nodename"]]
	if (!is.null(uname)) logger.info(c("Machine name:",uname))
}

#' combine.diffMeth.objs
#'
#' combine differential methylation objects (output from \code{rnb.run.differential}).
#' To be more precise, the \code{diffmeth} and \code{dm.enrich} are merged.
#' individual objects that are merged are assumed to belong to the same analysis
#' and vary only in their indexing of region types and comparisons
#' @param obj.list a list containing outputs from \code{rnb.run.differential}
#' @author Fabian Mueller
combine.diffMeth.objs <- function(obj.list){
	if (length(obj.list) < 1){
		return(NULL)
	}
	if (length(obj.list) == 1){
		return(obj.list[[1]])
	}
	ontologies <- c("BP","MF")
	ontol.list.empty <- rep(list(list()),length(ontologies))
	names(ontol.list.empty) <- ontologies
	dm.enrich.comb <- list(probe=list(),region=list())
	
	diffmeth <- obj.list[[1]]$diffmeth
	for (i in 1:length(obj.list)){
		dm <- obj.list[[i]]$diffmeth
		if (i > 1){
			diffmeth <- join.diffMeth(diffmeth,dm)
		}
		
		#merge dm.enrich parts
		new.comps <- get.comparisons(dm)
		if (i > 1){
			new.comps <- new.comps[!(new.comps %in% get.comparisons(diffmeth))]
		}
		n.new.comps <- length(new.comps)
		if (!is.null(obj.list[[i]]$dm.enrich)) {
			dmer <- obj.list[[i]]$dm.enrich$region
			#add empty lists for new comparisons
			new.comp.list.empty <- rep(list(ontol.list.empty),n.new.comps)
			names(new.comp.list.empty) <- new.comps
			dm.enrich.comb$region <- c(dm.enrich.comb$region,new.comp.list.empty)
			#fill in the empty entries with the new entries
			for (cc in names(dmer)){
				for (oo in names(dmer[[cc]])){
					for (rr in names(dmer[[cc]][[oo]])){
						if (!is.element(rr,names(dm.enrich.comb$region[[cc]][[oo]]))){
							dm.enrich.comb$region[[cc]][[oo]][[rr]] <- dmer[[cc]][[oo]][[rr]]
						}
					}
				}
			}
		}
	}
	
	#set emtpy enrichment analysis to NULL object
	if (length(dm.enrich.comb$region)==0){
		dm.enrich.comb <- NULL
	} else {
		class(dm.enrich.comb) <- "DiffMeth.enrich"
	}

	res <- list(report=NULL,diffmeth=diffmeth,dm.enrich=dm.enrich.comb)
	return(res)
}
