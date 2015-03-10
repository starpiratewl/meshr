\name{geneid.cummeRbund}
\alias{geneid.cummeRbund}
\docType{data}
\title{
Test data of significant differentially expressed genes used in cummeRbund package.
}
\description{
This RNA-Seq data were taken from three samples, "iPS", "hESC", and  "Fibroblasts". We first create two objects of gene sets, i.e., selected and universal genes, by comparing significantly regulated genes between iPS and hESC under the significance level of 0.05 by getSig method in \pkg{cummeRbund} package. 303 genes were finally choosed and 104 of them were differentially expressed.
}
\usage{data(geneid.cummeRbund)}
\source{
\url{http://www.bioconductor.org/packages/release/bioc/vignettes/cummeRbund/inst/doc/cummeRbund-manual.pdf}
}
\seealso{
\code{\link{sig.geneid.cummeRbund}}.
}
\examples{
data(geneid.cummeRbund)
names(geneid.cummeRbund)
}
\keyword{datasets}