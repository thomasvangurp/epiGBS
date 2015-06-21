test_report_2 <- function() {

	report.dir <- "RnBeads_report_testing"

	rnb.initialize.reports(report.dir)
	report <- createReport(file.path(report.dir, "myreport.html"), "An Eye Opener", page.title = "Fascinating Report", authors = c("Me","You"))

	stext <- "Here is some text for our awesome report"
	report <- rnb.add.section(report, "Adding a text section", stext)
	lorem1 <- c("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Maecenas vestibulum placerat lobortis.
				Ut viverra fringilla urna at rutrum. In hac habitasse platea dictumst.
				Vestibulum adipiscing rutrum libero et interdum. Etiam sed odio ac nibh ullamcorper congue.
				Proin ac ipsum elit. Ut porta lorem sed lorem molestie pharetra.",
				"Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia Curae;
				Cras ac augue eu turpis dignissim aliquam. Vivamus at arcu ligula, vel scelerisque nisi.
				Vivamus ac lorem libero, quis venenatis metus. Fusce et lectus at lectus vestibulum
				faucibus ac id sem.")
	report <- rnb.add.section(report,"A subsection",lorem1,level=2)
	lorem2 <- "Nunc congue molestie fringilla. Aliquam erat volutpat.
	Integer consequat turpis nec dolor pulvinar et vulputate magna adipiscing.
	Curabitur purus dolor, porttitor vel dapibus quis, eleifend at lacus.
	Cras at mauris est, quis aliquam libero. Nulla facilisi.
	Nam facilisis placerat aliquam. Morbi in odio non ligula mollis rhoncus
	et et erat. Maecenas ut dui nisl. Mauris consequat cursus augue quis euismod."
	rnb.add.paragraph(report,lorem2)

	rnb.add.paragraph(report,"TODO: Add more content here",paragraph.class="task")
	rnb.add.paragraph(report,"To be, or not to be, that is the question",paragraph.class="note")
	ll <- lapply(LETTERS[1:10],FUN=function(x){paste(rep(x,3),collapse="")})
	rnb.add.list(report,ll,type="u")
	rnb.add.list(report,ll,type="o")
	tt <- matrix(sapply(LETTERS[1:24],FUN=function(x){paste(rep(x,3),collapse="")}),ncol=4)
	colnames(tt) <- paste("Col",1:4)
	rownames(tt) <- paste("Row",1:6)
	rnb.add.table(report,tt,row.names = TRUE, first.col.header = TRUE,
			  tcaption = "A table")

	stext <- c("<p>Some German umlauts:
			   <ul>
			   <li>&Auml;</li>
			   <li>&Ouml;</li>
			   <li>&Uuml;</li>
			   </ul>
			   </p>",
			   '<p> A link: <a href="http://rnbeads.bioinf.mpi-inf.mpg.de/">RnBeads website</a></p>')
	report <- rnb.add.section(report, "HTML code", stext)

	report <- rnb.add.section(report, 'Plots', "")
	report.plot <- createReportPlot("hist_normal",report,create.pdf=TRUE,high.png=200)
	hist(rnorm(1000),breaks=50,col="blue")
	off(report.plot)
	desc <- 'A histogramm of samples drawn from the normal distribution
	Click <a href="myreport_pdfs/hist_normal.pdf">
	here</a> for the pdf version.'
	report <- rnb.add.figure(report, desc, report.plot)

	plotSine <- function(a,b,fname){
		report.plot <- createReportPlot(fname,report,create.pdf=FALSE,high.png=200)
		curve(sin(a*x)+b, -2*pi, 2*pi,ylim=c(-2,2),col="blue",lwd=2)
		off(report.plot)
		return(report.plot)
	}

	period.params <- c(a05=0.5,a1=1,a2=2)
	intercept.params <- c(im1=-1,i0=0,i1=1)
	plot.list <- list()
	for (aa in names(period.params)){
		for (bb in names(intercept.params)){
			fname <- paste("sinePlot",aa,bb,sep="_")
			current.plot <- plotSine(period.params[aa],intercept.params[bb],fname)
			plot.list <- c(plot.list,list(current.plot))
		}
		
	}
	setting.names <- list('period parameter' = period.params,
						  'intercept parameter' = intercept.params)
	description <- 'Sine for various parameters'
	report <- rnb.add.figure(report, description, plot.list, setting.names, selected.image=5)

	off(report)

	unlink(report.dir,recursive=TRUE)
}
