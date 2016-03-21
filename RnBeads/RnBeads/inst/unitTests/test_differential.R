########################################################################################################################
## test_differential.R
## creator: Fabian Müller
## ---------------------------------------------------------------------------------------------------------------------
## Unit tests for differential methylation functions
########################################################################################################################

test_rowWelchP <- function(){
	require(RnBeads.hg19)
	data(small.example.object)
	meth.mat <- meth(rnb.set.example)
	sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
	p.vals <- rowWelchP(meth.mat,sample.groups[[1]],sample.groups[[2]])
	passes <- is.numeric(p.vals) && all(p.vals>=0) && all(p.vals<=1)
	checkTrue(passes)
}

test_limmaP_1 <- function(){
	require(RnBeads.hg19)
	data(small.example.object)
	meth.mat <- meth(rnb.set.example)
	sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
	p.vals <- limmaP(meth.mat,sample.groups[[1]],sample.groups[[2]])
	passes <- is.numeric(p.vals) && all(p.vals>=0) && all(p.vals<=1)
	checkTrue(passes)
}

test_limmaP_2 <- function() {
	require(RnBeads.hg19)
	data(small.example.object)
	X <- meth(rnb.set.example)
	tbl <- data.frame(
		stype = c("a", "b", "a", "b", "a", "b", "b", "b", "b", "a", "b", "a"),
		gender = c("f", "m", "f", "f", "f", "f", "m", "m", "f", "f", "m", "m"),
		group = c("1", "2", "1", "1", "1", "1", "2", "2", "1", "1", "2", "2"))

	p.vals <- limmaP(X, tbl$stype == "a", tbl$stype == "b", tbl[, -1])
	checkTrue(is.numeric(p.vals) && all(p.vals>=0) && all(p.vals<=1))
}

test_computeDiffTab.extended.site <- function(){
	require(RnBeads.hg19)
	data(small.example.object)
	meth.mat <- meth(rnb.set.example)
	sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
	dm <- computeDiffTab.extended.site(meth.mat,sample.groups[[1]],sample.groups[[2]])
	checkTrue(dim(dm)[1]>1)
	checkEquals(dim(dm)[2],15)
}

test_computeDiffTab.default.region <- function(){
	require(RnBeads.hg19)
	data(small.example.object)
	meth.mat <- meth(rnb.set.example)
	sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
	dm.sites <- computeDiffTab.extended.site(meth.mat,sample.groups[[1]],sample.groups[[2]])
	map.regions.to.sites <- regionMapping(rnb.set.example,"promoters")
	dm.promoters <- computeDiffTab.default.region(dm.sites,map.regions.to.sites)
	checkTrue(dim(dm.promoters)[1]>1)
	checkEquals(dim(dm.promoters)[2],9)
}

test_combineTestPvalsMeth <- function(){
	p.vals <- 10^-c(0,1,5)
	res <- combineTestPvalsMeth(p.vals)
	checkEquals(res,0.00011, tolerance=1e-5,scale=1)
}

test_get.adjustment.variables <- function(){
	require(RnBeads.hg19)
	data(small.example.object)
	sample.groups <- rnb.sample.groups(rnb.set.example)[[1]]
	res <- get.adjustment.variables(rnb.set.example,sample.groups[[1]],sample.groups[[2]],"Cell_Line")
	checkTrue(is.data.frame(res))
	checkEquals(colnames(res),"Cell_Line")
}

test_get.comparison.info <- function(){
	require(RnBeads.hg19)
	data(small.example.object)
	cmp.info <- get.comparison.info(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
	checkTrue(is.list(cmp.info))
	checkTrue(length(cmp.info)>0)
	req.names <- c("comparison","pheno.colname","group.names","group.inds","paired","adj.sva","adj.celltype","adjustment.table","region.types")
	checkTrue(all(req.names %in% names(cmp.info[[1]])))
}

test_rnb.execute.computeDiffMeth <- function(){
	require(RnBeads.hg19)
	data(small.example.object)
	dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
	checkEquals(as.character(class(dm)),"RnBDiffMeth",checkNames=FALSE)
}

test_class_RnBDiffMeth <- function(){
	require(RnBeads.hg19)
	data(small.example.object)
	#DiffMeth example objects
	dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
	dm1 <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group"),region.types=c("genes","tiling"))
	dm1b <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group"),region.types=c("genes","tiling"))
	dm2 <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"),region.types=c("promoters"))
	dm3 <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Treatment"),region.types=c("genes","tiling"))
	checkEquals(as.character(class(dm)), "RnBDiffMeth")
	checkEquals(as.character(class(dm1)),"RnBDiffMeth")
	checkEquals(as.character(class(dm2)),"RnBDiffMeth")
	checkEquals(as.character(class(dm3)),"RnBDiffMeth")
	logger.start(c("Testing function:","get.region.types"))
		res <- get.region.types(dm)
		checkEquals(res,rnb.region.types.for.analysis(rnb.set.example),checkNames=FALSE)
	logger.completed()
	logger.start(c("Testing function:","get.comparisons"))
		res <- get.comparisons(dm)
		checkTrue(length(res)>1)
	logger.completed()
	logger.start(c("Testing function:","get.comparison.grouplabels"))
		res <- get.comparison.grouplabels(dm)
		checkTrue(nrow(res)>1)
		checkEquals(ncol(res),2)
	logger.completed()
	logger.start(c("Testing function:","get.site.test.method"))
		res <- get.site.test.method(dm)
		checkEquals(res,rnb.getOption("differential.site.test.method"))
	logger.completed()
	logger.start(c("Testing function:","get.table"))
		res <- get.table(dm,get.comparisons(dm)[1],"promoters",return.data.frame=TRUE)
		checkTrue(nrow(res)>1)
		checkEquals(ncol(res),14)
	logger.completed()
	logger.start(c("Testing function:","addDiffMethTable"))
		sample.groups <- rnb.sample.groups(rnb.set.example,"Sample_Group")[[1]]
		dmt.sites <- computeDiffTab.extended.site(meth(rnb.set.example),sample.groups[[1]],sample.groups[[2]])
		map.regions.to.sites <- regionMapping(rnb.set.example,"promoters")
		dmt.promoters <- computeDiffTab.default.region(dmt.sites,map.regions.to.sites)
		cmp.name <- get.comparisons(dm1)[1]
		grp.labs <- get.comparison.grouplabels(dm1)[1,]
		#add the promoter level differential methylation table
		dm.add <- addDiffMethTable(dm1,dmt.promoters,cmp.name,"promoters",grp.labs)
		checkEquals(get.region.types(dm.add),c("genes","tiling","promoters"),checkNames=FALSE)
	logger.completed()
	logger.start(c("Testing functions:","join.diffMeth, is.valid"))
		dm.join1 <- join.diffMeth(dm1b,dm2)
		#the following joint object is invalid, because some region type - comparison combinations are missing
		checkTrue(!is.valid(dm.join1))
		dm.join2 <- join.diffMeth(dm.join1,dm3)
		#after joining the missing information, the new object is valid
		checkTrue(is.valid(dm.join2))
	logger.completed()
	logger.start("Destructor")
		destroy(dm.join2)
		destroy(dm.add)
		# destroy(dm1) #already included in the destructor of dm.add
		# destroy(dm1b) #already included in the destructor of dm.join2
		# destroy(dm.join2) #already included in the destructor of dm.join2
		destroy(dm2)
		destroy(dm3)
	logger.completed()
}

test_differential <- function(){
	require(RUnit)
	logger.start(fname=NA)
	logger.start(c("Unit testing:","differential"))
		logger.start(c("Testing function:","rowWelchP"))
			test_rowWelchP()
		logger.completed()
		logger.start(c("Testing function:","limmaP"))
			test_rowWelchP()
		logger.completed()
		logger.start(c("Testing function:","computeDiffTab.extended.site"))
			test_computeDiffTab.extended.site()
		logger.completed()
		logger.start(c("Testing function:","computeDiffTab.default.region"))
			test_computeDiffTab.default.region()
		logger.completed()
		logger.start(c("Testing function:","combineTestPvalsMeth"))
			test_combineTestPvalsMeth()
		logger.completed()
		logger.start(c("Testing function:","get.adjustment.variables"))
			test_get.adjustment.variables()
		logger.completed()
		logger.start(c("Testing function:","get.comparison.info"))
			test_get.comparison.info()
		logger.completed()
		logger.start(c("Testing function:","rnb.execute.computeDiffMeth"))
			test_rnb.execute.computeDiffMeth()
		logger.completed()
		logger.start(c("Testing class:","RnBDiffMeth"))
			test_class_RnBDiffMeth()
		logger.completed()
	logger.completed()
}

test_differential()
