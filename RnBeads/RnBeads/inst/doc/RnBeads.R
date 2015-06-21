### R code from vignette source 'RnBeads.Rnw'

###################################################
### code chunk number 1: RnBeads.Rnw:104-107
###################################################
suppressPackageStartupMessages(library(RnBeads))
#for knitr: avoid tidy because it messes up line breaks
# opts_chunk$set(tidy=FALSE)


###################################################
### code chunk number 2: RnBeads.Rnw:117-118 (eval = FALSE)
###################################################
## source("http://rnbeads.mpi-inf.mpg.de/install.R")


###################################################
### code chunk number 3: RnBeads.Rnw:164-173 (eval = FALSE)
###################################################
## library(RnBeads)
## # Directory where your data is located
## data.dir <- "~/RnBeads/data/Ziller2011_PLoSGen_450K"
## idat.dir <- file.path(data.dir, "idat")
## sample.annotation <- file.path(data.dir, "sample_annotation.csv")
## # Directory where the output should be written to
## analysis.dir <- "~/RnBeads/analysis"
## # Directory where the report files should be written to
## report.dir <- file.path(analysis.dir, "reports")


###################################################
### code chunk number 4: RnBeads.Rnw:183-185 (eval = FALSE)
###################################################
## rnb.options(filtering.sex.chromosomes.removal=TRUE,
##     identifiers.column="Sample_ID")


###################################################
### code chunk number 5: RnBeads.Rnw:214-216 (eval = FALSE)
###################################################
## rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation,
##     data.dir=idat.dir, data.type="infinium.idat.dir")


###################################################
### code chunk number 6: RnBeads.Rnw:263-264 (eval = FALSE)
###################################################
## rnb.options(differential=FALSE)


###################################################
### code chunk number 7: RnBeads.Rnw:278-280 (eval = FALSE)
###################################################
## report.dir <- file.path(analysis.dir, "reports_details")
## rnb.initialize.reports(report.dir)


###################################################
### code chunk number 8: RnBeads.Rnw:284-285 (eval = FALSE)
###################################################
## logger.start(fname=NA)


###################################################
### code chunk number 9: RnBeads.Rnw:306-307 (eval = FALSE)
###################################################
## data.source <- c(idat.dir, sample.annotation)


###################################################
### code chunk number 10: RnBeads.Rnw:312-315 (eval = FALSE)
###################################################
## result <- rnb.run.import(data.source=data.source,
##     data.type="infinium.idat.dir", dir.reports=report.dir)
## rnb.set <- result$rnb.set


###################################################
### code chunk number 11: RnBeads.Rnw:323-325 (eval = FALSE)
###################################################
## rnb.set <- rnb.execute.import(data.source=data.source,
## 	data.type="infinium.idat.dir")


###################################################
### code chunk number 12: RnBeads.Rnw:336-349 (eval = FALSE)
###################################################
## rnb.set
## summary(pheno(rnb.set))
## dim(M(rnb.set))
## summary(M(rnb.set))
## dim(meth(rnb.set))
## summary(meth(rnb.set))
## dim(meth(rnb.set, type="promoters"))
## summary(meth(rnb.set, type="promoters"))
## summary(dpval(rnb.set))
## str(qc(rnb.set))
## samples(rnb.set)
## # Visualize the bimodal distribution of beta values in sample 5
## hist(meth(rnb.set)[, 5])


###################################################
### code chunk number 13: RnBeads.Rnw:353-367 (eval = FALSE)
###################################################
## # Remove samples
## rnb.set2 <- remove.samples(rnb.set, samples(rnb.set)[c(2, 6, 10)])
## setdiff(samples(rnb.set), samples(rnb.set2))
## 
## # Remove probes with a p-value > 0.01 in any sample
## any.bad.p.val <- apply(dpval(rnb.set)>0.01, 1, any)
## rnb.set3 <- remove.sites(rnb.set, any.bad.p.val)
## nrow(meth(rnb.set3))
## 
## # Add a sample annotation column indicating whether the given
## # sample represents an iPS cell line
## is.hiPSC <- pheno(rnb.set)[, "Sample_Group"]=="hiPSC"
## rnb.set <- addPheno(rnb.set, is.hiPSC, "is_hiPSC")
## pheno(rnb.set)


###################################################
### code chunk number 14: RnBeads.Rnw:411-413 (eval = FALSE)
###################################################
## rnb.set.geo <- rnb.execute.import(data.source="GSE31848",
## 	data.type="infinium.GEO")


###################################################
### code chunk number 15: RnBeads.Rnw:469-470 (eval = FALSE)
###################################################
## rnb.run.qc(rnb.set, report.dir)


###################################################
### code chunk number 16: RnBeads.Rnw:482-483 (eval = FALSE)
###################################################
## rnb.plot.control.boxplot(rnb.set, "BISULFITE CONVERSION I")


###################################################
### code chunk number 17: RnBeads.Rnw:487-488 (eval = FALSE)
###################################################
## rnb.plot.control.barplot(rnb.set, "BISULFITE CONVERSION I.2")


###################################################
### code chunk number 18: RnBeads.Rnw:491-492 (eval = FALSE)
###################################################
## rnb.plot.negative.boxplot(rnb.set)


###################################################
### code chunk number 19: RnBeads.Rnw:519-521 (eval = FALSE)
###################################################
## annot450 <- rnb.annotation2data.frame(rnb.get.annotation("probes450"))
## annot450[grep("rs", rownames(annot450)), ]


###################################################
### code chunk number 20: RnBeads.Rnw:535-536 (eval = FALSE)
###################################################
## rnb.plot.snp.heatmap(rnb.set.unfiltered)


###################################################
### code chunk number 21: RnBeads.Rnw:543-544 (eval = FALSE)
###################################################
## rnb.plot.snp.barplot(rnb.set.unfiltered, samples(rnb.set)[1])


###################################################
### code chunk number 22: RnBeads.Rnw:559-560 (eval = FALSE)
###################################################
## rnb.options(qc.snp.boxplot=TRUE)


###################################################
### code chunk number 23: RnBeads.Rnw:586-589 (eval = FALSE)
###################################################
## rnb.set.unfiltered <- rnb.set
## result <- rnb.run.preprocessing(rnb.set.unfiltered, dir.reports=report.dir)
## rnb.set <- result$rnb.set


###################################################
### code chunk number 24: RnBeads.Rnw:604-628 (eval = FALSE)
###################################################
## nrow(meth(rnb.set.unfiltered)) # the number of sites in the unfiltered object
## # Remove probes outside of CpG context
## rnb.set.filtered <- rnb.execute.context.removal(rnb.set.unfiltered)$dataset
## nrow(meth(rnb.set.filtered)) # the number of CpG sites in the unfiltered object
## # SNP filtering allowing no SNPs in the probe sequence
## rnb.set.filtered <-
##     rnb.execute.snp.removal(rnb.set.filtered, snp="any")$dataset
## # the number of CpG sites in the unfiltered object
## # that do not contain a SNP
## nrow(meth(rnb.set.filtered)) 
## 
## # Remove CpGs on sex chromosomes
## rnb.set.filtered <- rnb.execute.sex.removal(rnb.set.filtered)$dataset
## nrow(meth(rnb.set.filtered))
## # Remove probes and samples based on a greedy approach
## rnb.set.filtered <- rnb.execute.greedycut(rnb.set.filtered)$dataset
## nrow(meth(rnb.set.filtered))
## # Remove probes containing NA for beta values
## rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered)$dataset
## nrow(meth(rnb.set.filtered))
## # Remove probes for which the beta values have low standard deviation
## rnb.set.filtered <-
##     rnb.execute.variability.removal(rnb.set.filtered, 0.005)$dataset
## nrow(meth(rnb.set.filtered))


###################################################
### code chunk number 25: RnBeads.Rnw:679-681 (eval = FALSE)
###################################################
## rnb.set.norm <- rnb.execute.normalization(rnb.set.unfiltered, method="illumina",
##     bgcorr.method="methylumi.noob")


###################################################
### code chunk number 26: RnBeads.Rnw:691-692 (eval = FALSE)
###################################################
## rnb.options(inference=TRUE)


###################################################
### code chunk number 27: RnBeads.Rnw:696-697 (eval = FALSE)
###################################################
## rnb.set <- rnb.run.inference(rnb.set, report.dir)$rnb.set


###################################################
### code chunk number 28: RnBeads.Rnw:706-711 (eval = FALSE)
###################################################
## rnb.options(
## 	inference.targets.sva=c("Sample_Group"),
## 	inference.sva.num.method="be"
## )
## rnb.set.new <- rnb.run.inference(rnb.set, report.dir)$rnb.set


###################################################
### code chunk number 29: RnBeads.Rnw:715-717 (eval = FALSE)
###################################################
## sva.obj <- rnb.execute.sva(rnb.set, cmp.cols="Sample_Group", numSVmethod="be")
## rnb.set <- rnb.set.add.covariates.sva(rnb.set, sva.obj)


###################################################
### code chunk number 30: RnBeads.Rnw:728-729 (eval = FALSE)
###################################################
## rnb.options(inference.reference.methylome.column="CellType")


###################################################
### code chunk number 31: RnBeads.Rnw:754-756 (eval = FALSE)
###################################################
## ct.obj<-rnb.execute.ct.estimation(rnb.set, cell.type.column="CellType",
##     test.max.markers=10000, top.markers=500)


###################################################
### code chunk number 32: RnBeads.Rnw:763-764 (eval = FALSE)
###################################################
## ct.obj$contributions


###################################################
### code chunk number 33: RnBeads.Rnw:767-768 (eval = FALSE)
###################################################
## ct.obj$contributions.constr


###################################################
### code chunk number 34: RnBeads.Rnw:777-778 (eval = FALSE)
###################################################
## rnb.run.exploratory(rnb.set, report.dir)


###################################################
### code chunk number 35: RnBeads.Rnw:785-791 (eval = FALSE)
###################################################
## dred.sites <- rnb.execute.dreduction(rnb.set)
## dred.promoters <- rnb.execute.dreduction(rnb.set, target="promoters")
## dred <- list(sites=dred.sites, promoters=dred.promoters)
## sample.colors <- ifelse(pheno(rnb.set)$Sample_Group=="hESC", "red", "blue")
## plot(dred[[1]]$mds$euclidean, col=sample.colors,
##     xlab="Dimension 1", ylab="Dimension 2")


###################################################
### code chunk number 36: RnBeads.Rnw:802-804 (eval = FALSE)
###################################################
## assoc <- rnb.execute.batcheffects(rnb.set, pcoordinates=dred)
## str(assoc)


###################################################
### code chunk number 37: RnBeads.Rnw:808-809 (eval = FALSE)
###################################################
## rnb.options(exploratory.columns=c("Sample_Group", "Passage_No"))


###################################################
### code chunk number 38: RnBeads.Rnw:812-813 (eval = FALSE)
###################################################
## assoc.qc <- rnb.execute.batch.qc(rnb.set, pcoordinates=dred)


###################################################
### code chunk number 39: RnBeads.Rnw:817-822 (eval = FALSE)
###################################################
## probe.types <- as.character(annotation(rnb.set)[, "Design"])
## pdf(file="deviationPlot_probeType.pdf", width=11)
## deviation.plot.beta(meth(rnb.set), probe.types,
##     c.legend=c("I"="blue", "II"="red"))
## dev.off()


###################################################
### code chunk number 40: RnBeads.Rnw:836-839 (eval = FALSE)
###################################################
## clusterings.sites <- rnb.execute.clustering(rnb.set, region.type="sites")
## clusterings.promoters <-
##     rnb.execute.clustering(rnb.set, region.type="promoters")


###################################################
### code chunk number 41: RnBeads.Rnw:844-856 (eval = FALSE)
###################################################
## # Get the methylation values
## X <- meth(rnb.set, type="promoters")[1:100, ]
## # Convert the clustering result to a dendrogram
## # Index 7 holds euclidean distance, average linkage clustering
## cresult <- clusterings.promoters[[7]]@result
## attr(cresult, "class") <- "hclust"
## cresult <- as.dendrogram(cresult)
## # Save the heatmap as pdf file
## pdf(file="promoter_heatmap.pdf")
## heatmap.2(X, Rowv=TRUE, Colv=cresult, dendrogram="both",
##     scale="none", trace="none")
## dev.off()


###################################################
### code chunk number 42: RnBeads.Rnw:870-871 (eval = FALSE)
###################################################
## rnb.run.differential(rnb.set, report.dir)


###################################################
### code chunk number 43: RnBeads.Rnw:874-882 (eval = FALSE)
###################################################
## # Specify the sample annotation table columns for which
## # differential methylation is to be computed
## cmp.cols <- "Sample_Group"
## # Specify the region types
## reg.types <- c("genes", "promoters")
## # Conduct the analysis
## diffmeth <- rnb.execute.computeDiffMeth(rnb.set, cmp.cols,
##     region.types=reg.types)


###################################################
### code chunk number 44: RnBeads.Rnw:888-889 (eval = FALSE)
###################################################
## str(diffmeth)


###################################################
### code chunk number 45: RnBeads.Rnw:892-898 (eval = FALSE)
###################################################
## comparison <- get.comparisons(diffmeth)[1]
## tab.sites <- get.table(diffmeth, comparison, "sites", return.data.frame=TRUE)
## str(tab.sites)
## tab.promoters <- get.table(diffmeth, comparison, "promoters",
## 	return.data.frame=TRUE)
## str(tab.promoters)


###################################################
### code chunk number 46: RnBeads.Rnw:929-932 (eval = FALSE)
###################################################
## dmps <- tab.sites[order(tab.sites[, "combinedRank"]), ]
## summary(dmps[1:100, ])
## summary(dmps[1:1000, ])


###################################################
### code chunk number 47: RnBeads.Rnw:936-938 (eval = FALSE)
###################################################
## which.diffmeth <- abs(dmps[, "mean.diff"])>0.2 & dmps$diffmeth.p.adj.fdr<0.05
## dmps.cutoff <- dmps[which.diffmeth, ]


###################################################
### code chunk number 48: RnBeads.Rnw:943-946 (eval = FALSE)
###################################################
## dmrs <- get.table(diffmeth, comparison, "promoters")
## plot(dmrs[, "mean.mean.diff"], -log10(dmrs[, "comb.p.val"]),
##     xlab="mean difference", ylab="-log10(combined p-value)", col="blue")


###################################################
### code chunk number 49: RnBeads.Rnw:976-978 (eval = FALSE)
###################################################
## rnb.options("differential.comparison.columns"=c("diseaseState"),
##     "columns.pairing"=c("diseaseState"="individual"))


###################################################
### code chunk number 50: RnBeads.Rnw:982-984 (eval = FALSE)
###################################################
## rnb.options("differential.comparison.columns"=NULL,
##     "columns.pairing"=NULL)


###################################################
### code chunk number 51: RnBeads.Rnw:990-993 (eval = FALSE)
###################################################
## rnb.options("covariate.adjustment.columns"=c("Passage_No"))
## diffmeth.adj <- rnb.execute.computeDiffMeth(rnb.set, cmp.cols,
##     region.types=reg.types)


###################################################
### code chunk number 52: RnBeads.Rnw:999-1001 (eval = FALSE)
###################################################
## rnb.options(export.to.csv=TRUE)
## rnb.run.tnt(rnb.set, report.dir)


###################################################
### code chunk number 53: RnBeads.Rnw:1110-1112 (eval = FALSE)
###################################################
## rnb.run.analysis(dir.reports=report.dir, data.source=data.source,
##     data.type="bs.bed.dir")


###################################################
### code chunk number 54: RnBeads.Rnw:1142-1152 (eval = FALSE)
###################################################
## # Remove CpGs on sex chromosomes
## rnb.set.filtered <- rnb.execute.sex.removal(rnb.set.unfiltered)$dataset
## # Remove sites that have an exceptionally high coverage
## rnb.set.filtered <-
##     rnb.execute.highCoverage.removal(rnb.set.filtered)$dataset
## # Remove sites containing NA for beta values
## rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered)$dataset
## # Remove sites for which the beta values have low standard deviation
## rnb.set.filtered <-
##     rnb.execute.variability.removal(rnb.set.filtered, 0.005)$dataset


###################################################
### code chunk number 55: RnBeads.Rnw:1165-1166 (eval = FALSE)
###################################################
## rnb.options()


###################################################
### code chunk number 56: RnBeads.Rnw:1169-1170 (eval = FALSE)
###################################################
## rnb.getOption("analyze.sites")


###################################################
### code chunk number 57: RnBeads.Rnw:1173-1175 (eval = FALSE)
###################################################
## rnb.options(filtering.sex.chromosomes.removal=TRUE,
##     colors.category=rainbow(8))


###################################################
### code chunk number 58: RnBeads.Rnw:1180-1181 (eval = FALSE)
###################################################
## ?rnb.options


###################################################
### code chunk number 59: RnBeads.Rnw:1192-1195 (eval = FALSE)
###################################################
## library(RnBeads)
## rnb.get.annotation(type="CpG") # all CpGs in the human genome
## rnb.get.annotation(type="probes450") # all Infinium 450k probes


###################################################
### code chunk number 60: RnBeads.Rnw:1199-1201 (eval = FALSE)
###################################################
## probes.450k <- rnb.annotation2data.frame(rnb.get.annotation(type="probes450"))
## head(probes.450k)


###################################################
### code chunk number 61: RnBeads.Rnw:1205-1206
###################################################
rnb.region.types()


###################################################
### code chunk number 62: RnBeads.Rnw:1210-1211 (eval = FALSE)
###################################################
## rnb.get.annotation("promoters",assembly="mm9")


###################################################
### code chunk number 63: RnBeads.Rnw:1215-1217 (eval = FALSE)
###################################################
## head(annotation(rnb.set, type="sites"))
## head(annotation(rnb.set, type="genes"))


###################################################
### code chunk number 64: RnBeads.Rnw:1221-1224 (eval = FALSE)
###################################################
## aa <- annotation(rnb.set, type="promoters")
## annotated.dmrs <- data.frame(aa, dmrs)
## head(annotated.dmrs)


###################################################
### code chunk number 65: RnBeads.Rnw:1234-1264 (eval = FALSE)
###################################################
## # Retrieve the chromHMM state segmentation from UCSC
## library(rtracklayer)
## mySession <- browserSession()
## genome(mySession) <- "hg19"
## tab.chromHMM.h1 <- getTable(ucscTableQuery(mySession,
##     track="wgEncodeBroadHmm", table="wgEncodeBroadHmmH1hescHMM"))
## 
## # Filter for enhancer states
## tab.enhancers <- tab.chromHMM.h1[grep("Enhancer", tab.chromHMM.h1$name), ]
## 
## # Select the interesting parts of the table and rename columns
## tab.enhancers <- tab.enhancers[, c("chrom", "chromStart", "chromEnd", "name")]
## colnames(tab.enhancers) <- c("chromosome", "start", "end", "name")
## 
## # Create RnBeads annotation by providing a data.frame
## rnb.set.annotation("enhancersH1EscChromHMM", tab.enhancers, assembly="hg19")
## 
## # Set the options to include the enhancer annotation
## rnb.options(region.types=c(rnb.getOption("region.types"),
##     "enhancersH1EscChromHMM"))
## # Parse the input again, this time with the enhancer annotation added
## rnb.set.enh <-
##     rnb.execute.import(data.source=data.source, data.type="idat.dir")
## rnb.set.enh
## 
## # Annotation and methylation levels of enhancer regions in this object
## annot.enh <- annotation(rnb.set.enh, "enhancersH1EscChromHMM")
## head(annot.enh)
## meth.enh <- meth(rnb.set.enh, "enhancersH1EscChromHMM")
## head(meth.enh)


###################################################
### code chunk number 66: RnBeads.Rnw:1268-1277 (eval = FALSE)
###################################################
## annotation.filename <- file.path(analysis.dir,
##     "annotation_hg19_enhancersH1EscChromHMM.RData")
## # Save the enhancer annotation to disk
## rnb.save.annotation(annotation.filename, "enhancersH1EscChromHMM",
##     assembly="hg19")
## # Load the enhancer annotation as a duplicate
## rnb.load.annotation(annotation.filename, "enhancersH1EscChromHMM_duplicate")
## # Check that the annotation has been successfully loaded
## rnb.region.types()


###################################################
### code chunk number 67: RnBeads.Rnw:1283-1288
###################################################
logger.start(fname=NA)
parallel.isEnabled()
num.cores <- 2
parallel.setup(num.cores)
parallel.isEnabled()


###################################################
### code chunk number 68: RnBeads.Rnw:1293-1294
###################################################
if (parallel.isEnabled()) parallel.getNumWorkers()


###################################################
### code chunk number 69: RnBeads.Rnw:1298-1299
###################################################
parallel.disable()


###################################################
### code chunk number 70: RnBeads.Rnw:1315-1316 (eval = FALSE)
###################################################
## rnb.options(disk.dump.big.matrices=TRUE)


###################################################
### code chunk number 71: RnBeads.Rnw:1320-1321 (eval = FALSE)
###################################################
## options(fftempdir="/path/to/a/very/large/directory/")


###################################################
### code chunk number 72: RnBeads.Rnw:1326-1327 (eval = FALSE)
###################################################
## tempdir()


###################################################
### code chunk number 73: RnBeads.Rnw:1333-1334 (eval = FALSE)
###################################################
## rnb.options(logging.disk=TRUE)


###################################################
### code chunk number 74: RnBeads.Rnw:1337-1338 (eval = FALSE)
###################################################
## rnb.options(enforce.memory.management=TRUE)


###################################################
### code chunk number 75: RnBeads.Rnw:1350-1351 (eval = FALSE)
###################################################
## rnb.set.disk <- rnb.execute.import(data.source=data.source, data.type="idat.dir")


###################################################
### code chunk number 76: RnBeads.Rnw:1358-1360 (eval = FALSE)
###################################################
## save.dir <- file.path(report.dir, "analysis")
## save.rnb.set(rnb.set.disk, path=save.dir, archive=TRUE)


###################################################
### code chunk number 77: RnBeads.Rnw:1365-1366 (eval = FALSE)
###################################################
## rns <- load.rnb.set(path=paste0(save.dir, ".zip"))


###################################################
### code chunk number 78: RnBeads.Rnw:1370-1371 (eval = FALSE)
###################################################
## destroy(rnb.set.disk)


###################################################
### code chunk number 79: RnBeads.Rnw:1393-1395 (eval = FALSE)
###################################################
## rnb.plot.betadistribution.sampleGroups(meth(rnb.set),
##     rnb.sample.groups(rnb.set)[["Sample_Group"]], "ESC or iPSC")


###################################################
### code chunk number 80: RnBeads.Rnw:1408-1412 (eval = FALSE)
###################################################
## theme_set(theme_bw())
## rnb.options(colors.category = c("red", "blue", "grey"))
## print(rnb.plot.dreduction(rnb.set, dimensions = c(2, 5),
##     point.colors="Predicted Gender"))


###################################################
### code chunk number 81: RnBeads.Rnw:1423-1429 (eval = FALSE)
###################################################
## # Coordinates around the HOXD3 locus
## chrom <- "chr2"
## start <- 177010000
## end <- 177040000
## sample.grouping <- rnb.sample.groups(rnb.set)[[1]]
## rnb.plot.locus.profile(rnb.set, chrom, start, end, grps=sample.grouping)


###################################################
### code chunk number 82: RnBeads.Rnw:1445-1481 (eval = FALSE)
###################################################
## # set the target comparison column and region types for 
## # differential methylation analysis
## cmp.cols <- "Sample_Group"
## cmp.name <- "hESC vs. hiPSC (based on Sample_Group)"
## reg.types <- c("cpgislands", "promoters")
## # if you have not done so yet, compute the SVs and add them to the RnBSet object
## sva.obj <- rnb.execute.sva(rnb.set,cmp.cols=cmp.cols,numSVmethod="be")
## rnb.set.sva <- rnb.set.add.covariates.sva(rnb.set, sva.obj)
## 
## # compute differential methylation tables: unadjusted and SVA adjusted
## diffmeth.base <- rnb.execute.computeDiffMeth(
## 	rnb.set.sva, cmp.cols, region.types=reg.types,
## 	adjust.sva=FALSE
## )
## diffmeth.sva <- rnb.execute.computeDiffMeth(
## 	rnb.set.sva, cmp.cols, region.types=reg.types,
## 	adjust.sva=TRUE, pheno.cols.adjust.sva=cmp.cols
## )
## dm.tab.base <- get.table(diffmeth.base, cmp.name, "sites",
## 	return.data.frame=TRUE)
## dm.tab.sva  <- get.table(diffmeth.sva,  cmp.name, "sites",
## 	return.data.frame=TRUE)
## # compute quantiles of -log10 p-values and prepare a data.frame for plotting
## p.val.perc.base <- quantile(-log10(dm.tab.base$diffmeth.p.val),
## 	probs = seq(0, 1, 0.01))
## p.val.perc.sva  <- quantile(-log10(dm.tab.sva$diffmeth.p.val),
## 	probs = seq(0, 1, 0.01))
## df <- data.frame(
## 	neg.log10.p.val.base=p.val.perc.base,
## 	neg.log10.p.val.sva=p.val.perc.sva,
## 	quantile=seq(0, 1, 0.01)
## )
## # plot using ggplot2
## library(ggplot2)
## ggplot(df,aes(x=neg.log10.p.val.base,y=neg.log10.p.val.sva,color=quantile)) + 
## geom_point() + geom_abline(intercept=0, slope=1) + xlim(0, 5) + ylim(0, 5)


###################################################
### code chunk number 83: RnBeads.Rnw:1497-1506 (eval = FALSE)
###################################################
## comparison <- "hESC vs. hiPSC (based on Sample_Group)"
## dframe <- get.table(diffmeth, comparison, "sites", return.data.frame=TRUE)
## # define the probes with an FDR corrected p-value below 0.05
## # as differentially methylated
## isDMP <- dframe[,"diffmeth.p.adj.fdr"] < 0.05
## 
## create.densityScatter(dframe, is.special=isDMP, sparse.points=0.001,
## 	add.text.cor=TRUE) + labs(x="mean beta: hESC", y="mean beta: hiPSC") + 
## 	coord_fixed()


###################################################
### code chunk number 84: RnBeads.Rnw:1515-1516 (eval = FALSE)
###################################################
## rnb.options(differential.site.test.method="refFreeEWAS")


###################################################
### code chunk number 85: RnBeads.Rnw:1527-1528 (eval = FALSE)
###################################################
## rnb.options(export.to.ewasher=TRUE)


###################################################
### code chunk number 86: RnBeads.Rnw:1539-1555 (eval = FALSE)
###################################################
## library(RnBeads)
## # specify the xml file for your analysis
## xml.file <- "MY_ANALYSIS_SETTINGS.XML"
## # set the cluster architecture specific to your environment. We use SGE here
## arch <- new("ClusterArchitectureSGE")
## # initialize the object containing the RnBeads specific cluster configuration
## rnb.cr <- new("RnBClusterRun",arch)
## # set up the cluster so that 32GB of memory are required
## # (SGE resource is called "mem_free")
## rnb.cr <- setModuleResourceRequirements(rnb.cr,c(mem_free="32G"),"all")
## # set up the cluster to use 4 cores on each node for all modules
## rnb.cr <- setModuleNumCores(rnb.cr,4L,"all")
## # set up the cluster to use 2 cores for the exploratory analysis module
## rnb.cr <- setModuleNumCores(rnb.cr,2L,"exploratory")
## # run the actual analysis
## run(rnb.cr, "rnbeads_analysis", xml.file)


###################################################
### code chunk number 87: RnBeads.Rnw:1565-1599 (eval = FALSE)
###################################################
## # Define the new class, extending RnBeads ClusterArchitecture class
## setClass("ClusterArchitectureMyEnv",
## 	contains = "ClusterArchitecture"
## )
## #define the getSubCmdTokens method
## setMethod("getSubCmdTokens",
##   signature(
##     object="ClusterArchitectureMyEnv"
##   ),
##   function(
##     object,
##     cmd.tokens,
##     log,
##     job.name = "",
##     res.req = character(0),
##     depend.jobs = character(0)
##   ) {
##     dependency.token <- NULL
##     if (length(depend.jobs)>0){
##       dependency.token <- c( "--waitForJobs", paste(depend.jobs,collapse=","))
##     }
##     job.name.token <- NULL
##     if (nchar(job.name)>0) {
##       job.name.token <- c("--jobName",job.name)
##     }
##     result <- c(
##       "iSubmit", "--myParamter", "5",
##       dependency.token, job.name.token,
##       "--logFile",log,
##       paste0("'",paste(cmd.tokens,collapse=" "),"'")
##     )
##     return(result)
##   }
## )


###################################################
### code chunk number 88: RnBeads.Rnw:1606-1607 (eval = FALSE)
###################################################
## ?`getSubCmdTokens,ClusterArchitectureSGE-method`


###################################################
### code chunk number 89: RnBeads.Rnw:1622-1626 (eval = FALSE)
###################################################
## report.dir <- "RnBeads_report"
## rnb.initialize.reports(report.dir)
## report <- createReport(file.path(report.dir, "myreport.html"), "An Eye Opener",
##     page.title="Fascinating Report", authors=c("Me", "You"))


###################################################
### code chunk number 90: RnBeads.Rnw:1629-1630 (eval = FALSE)
###################################################
## str(report)


###################################################
### code chunk number 91: RnBeads.Rnw:1642-1667 (eval = FALSE)
###################################################
## stext <- "Here is some text for our awesome report"
## report <- rnb.add.section(report, "Adding a text section", stext)
## lorem1 <- c("Lorem ipsum dolor sit amet, consectetur adipiscing elit.
##             Maecenas vestibulum placerat lobortis. Ut viverra fringilla
##             urna at rutrum. In hac habitasse platea dictumst. Vestibulum
##             adipiscing rutrum libero et interdum. Etiam sed odio ac nibh
##             ullamcorper congue. Proin ac ipsum elit. Ut porta lorem sed
##             lorem molestie pharetra.",
##             "Vestibulum ante ipsum primis in faucibus orci luctus et
##             ultrices posuere cubilia Curae; Cras ac augue eu turpis
##             dignissim aliquam. Vivamus at arcu ligula, vel scelerisque nisi.
##             Vivamus ac lorem libero, quis venenatis metus. Fusce et lectus
##             at lectus vestibulum faucibus ac id sem.")
## report <- rnb.add.section(report, "A subsection", lorem1, level=2)
## lorem2 <- "Nunc congue molestie fringilla. Aliquam erat volutpat.
##            Integer consequat turpis nec dolor pulvinar et vulputate magna
##            adipiscing. Curabitur purus dolor, porttitor vel dapibus quis,
##            eleifend at lacus. Cras at mauris est, quis aliquam libero.
##            Nulla facilisi. Nam facilisis placerat aliquam. Morbi in odio
##            non ligula mollis rhoncus et et erat. Maecenas ut dui nisl.
##            Mauris consequat cursus augue quis euismod."
## rnb.add.paragraph(report, lorem2)
## rnb.add.paragraph(report, "TODO: Add content here", paragraph.class="task")
## rnb.add.paragraph(report, "To be or not to be, that is the question",
##     paragraph.class="note")


###################################################
### code chunk number 92: RnBeads.Rnw:1671-1680 (eval = FALSE)
###################################################
## report <- rnb.add.section(report, "Lists and Tables", "")
## ll <- lapply(LETTERS[1:10], function(x) { paste(rep(x, 3), collapse="") })
## rnb.add.list(report, ll, type="o")
## tt <- matrix(sapply(LETTERS[1:24],
##     function(x) { paste(rep(x, 3), collapse="") }), ncol=4)
## colnames(tt) <- paste("Col", 1:4)
## rownames(tt) <- paste("Row", 1:6)
## rnb.add.table(report, tt, row.names=TRUE, first.col.header=TRUE,
##     tcaption="A table")


###################################################
### code chunk number 93: RnBeads.Rnw:1684-1694 (eval = FALSE)
###################################################
## stext <- c('<p>Some German umlauts:
##                 <ul>
##                     <li>&Auml;</li>
##                     <li>&Ouml;</li>
##                     <li>&Uuml;</li>
##                 </ul>
##             </p>',
##             '<p> A link: <a href="http://rnbeads.mpi-inf.mpg.de/">
##             RnBeads website</a></p>')
## report <- rnb.add.section(report, "HTML code", stext)


###################################################
### code chunk number 94: RnBeads.Rnw:1699-1708 (eval = FALSE)
###################################################
## report <- rnb.add.section(report, "Plots", "")
## report.plot <- createReportPlot("hist_normal", report, create.pdf=TRUE,
##     high.png=200)
## hist(rnorm(1000), breaks=50, col="blue")
## off(report.plot)
## desc <- 'A histogramm of samples drawn from the normal distribution
##          Click <a href="myreport_pdfs/hist_normal.pdf">
##          here</a> for the pdf version.'
## report <- rnb.add.figure(report, desc, report.plot)


###################################################
### code chunk number 95: RnBeads.Rnw:1714-1739 (eval = FALSE)
###################################################
## # Plot a sine curve
## plotSine <- function(a, b, fname) {
##     report.plot <- createReportPlot(fname, report, create.pdf=FALSE,
##         high.png=200)
##     curve(sin(a*x)+b, -2*pi, 2*pi, ylim=c(-2, 2), col="blue", lwd=2)
##     return(off(report.plot))
## }
## 
## # Set parameters for different sine curves
## period.params <- c(a05=0.5, a1=1, a2=2)
## intercept.params <- c(im1=-1, i0=0, i1=1)
## plot.list <- list()
## for (aa in names(period.params)){
##     for (bb in names(intercept.params)){
##         fname <- paste("sinePlot", aa, bb, sep="_")
##         current.plot <- plotSine(period.params[aa], intercept.params[bb],
##             fname)
##         plot.list <- c(plot.list, current.plot)
##     }
## }
## setting.names <- list('period parameter'=period.params,
##     'intercept parameter'=intercept.params)
## description <- 'Sine for various parameters'
## report <- rnb.add.figure(report, description, plot.list, setting.names,
##      selected.image=5)


###################################################
### code chunk number 96: RnBeads.Rnw:1744-1745 (eval = FALSE)
###################################################
## off(report)


###################################################
### code chunk number 97: RnBeads.Rnw:1754-1755
###################################################
logger.start("Logging tutorial", fname=NA)


###################################################
### code chunk number 98: RnBeads.Rnw:1759-1760
###################################################
logger.isinitialized()


###################################################
### code chunk number 99: RnBeads.Rnw:1764-1765
###################################################
logger.start("Some tedious task")


###################################################
### code chunk number 100: RnBeads.Rnw:1769-1770
###################################################
logger.completed()


###################################################
### code chunk number 101: RnBeads.Rnw:1773-1791
###################################################
f <- function(){
    logger.start("Another tedious task")
    Sys.sleep(2)
        logger.start("Some tedious subtask")
        Sys.sleep(2)
        logger.completed()
        logger.start("Another, more tedious subtask")
        Sys.sleep(2)
            logger.start("Some tedious subsubtask")
            Sys.sleep(2)
            logger.completed()
            logger.start("Some even more tedious subsubtask")
            Sys.sleep(3)
            logger.completed()
        logger.completed()
    logger.completed()
}
f()


###################################################
### code chunk number 102: RnBeads.Rnw:1802-1824
###################################################
draw.lotto <- function(count=6) {
	urn <- 1:49
	primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
	logger.start(c("Drawing", count, "numbers from", length(urn)))
	for (i in 1:count) {
		if (i > 9) {
			msg <- "Too many numbers drawn"
			logger.error(msg, terminate=FALSE)
			stop(msg)
		}
		drawn <- sample(urn, 1)
		urn <- setdiff(urn, drawn)
		logger.info(c("The next number is", drawn))
		if (drawn %in% primes) {
			logger.warning("A prime number was drawn")
		}
	}
	logger.completed()
}
draw.lotto()
try(draw.lotto(15))
# draw.lotto(15)


###################################################
### code chunk number 103: RnBeads.Rnw:1888-1894 (eval = FALSE)
###################################################
## install.packages(c("doParallel", "fields", "foreach", "gplots", "matrixStats",
##     "mclust", "mgcv", "RPMM", "wordcloud", "XML"))
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("biomaRt", "GenomicRanges", "GEOquery", "ggbio", "GOstats", "Gviz",
##     "IlluminaHumanMethylation450k.db", "IlluminaHumanMethylation450kmanifest",
## 	"methylumi", "rtracklayer"))


###################################################
### code chunk number 104: RnBeads.Rnw:1905-1906 (eval = FALSE)
###################################################
## library(RnBeads)


###################################################
### code chunk number 105: RnBeads.Rnw:1909-1910 (eval = FALSE)
###################################################
## install.packages("XYZ")


###################################################
### code chunk number 106: RnBeads.Rnw:1913-1915 (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("XYZ")


###################################################
### code chunk number 107: RnBeads.Rnw:1930-1932 (eval = FALSE)
###################################################
## install.packages("[PATH_OF_DOWNLOADED_DATA_ARCHIVE]", repos=NULL)
## install.packages("[PATH_OF_DOWNLOADED_MAIN_ARCHIVE]", repos=NULL)


