% Generated by roxygen2 (4.0.2): do not edit by hand
\docType{methods}
\name{summarized.regions,RnBSet-method}
\alias{summarized.regions}
\alias{summarized.regions,RnBSet-method}
\title{summarized.regions-methods}
\usage{
\S4method{summarized.regions}{RnBSet}(object)
}
\arguments{
\item{object}{Methylation dataset of interest.}
}
\value{
\code{character} vector listing all genomic annotations summarized in the given dataset. If the dataset
        contains methylation in sites only, an empty vector is returned.
}
\description{
Gets the genomic annotations for which methylation data is present in the \code{RnBSet} object.
}
\examples{
\dontrun{
library(RnBeads.hg19)
data(small.example.object)
summarized.regions(rnb.set.example)
}
}
\author{
Yassen Assenov
}
\seealso{
\code{\link[=summarize.regions,RnBSet-method]{summarize.regions}} for calculating region-wise methylation in a dataset;
         \code{\link{rnb.set.annotation}} for adding or replacing a region annotation table
}

