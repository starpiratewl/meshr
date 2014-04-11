\name{sig.geneid.cummeRbund}
\alias{sig.geneid.cummeRbund}
\docType{data}
\title{
Test data of significant differentially expressed genes used in cummeRbund package.
}
\description{
This RNA-Seq data were taken from three samples, "iPS", "hESC", and  "Fibroblasts". We first create two objects of gene sets, i.e., selected and universal genes, by comparing significantly regulated genes between iPS and hESC under the significance level of 0.05 by getSig method in \pkg{cummeRbund} package. 303 genes were finally choosed and 104 of them were differentially expressed.
}
\usage{data(sig.geneid.cummeRbund)}
\source{
\url{http://www.bioconductor.org/packages/release/bioc/vignettes/cummeRbund/inst/doc/cummeRbund-manual.pdf}
}
\seealso{
\code{\link{geneid.cummeRbund}}.
}
\examples{
data(sig.geneid.cummeRbund)
names(sig.geneid.cummeRbund)

## This data is also available by following scripts.

# library(cummeRbund)
# library(org.Hs.eg.db)
# cuff <- readCufflinks(dir = system.file("extdata", package = "cummeRbund"))

# gene.symbols <- annotation(genes(cuff))[,4]
# mySigGeneIds <- getSig(cuff,x='hESC',y='iPS',alpha=0.05,level='genes')
# mySigGenes <- getGenes(cuff,mySigGeneIds)

# sig.gene.symbols <- annotation(mySigGenes)[,4]
# gene.symbols <- gene.symbols[!is.na(gene.symbols)]
# sig.gene.symbols <- sig.gene.symbols[!is.na(sig.gene.symbols)]

# geneid <- select(org.Hs.eg.db, keys=gene.symbols, keytype="SYMBOL", columns="ENTREZID")
# sig.geneid <- select(org.Hs.eg.db, keys=sig.gene.symbols, keytype="SYMBOL", columns="ENTREZID")

# na.index1 <- which(is.na(geneid[,2]))
# for (i in na.index1){
#   s <- unlist(strsplit(as.character(geneid[i,][1]), ","))[1]
#   sym <- get(s, org.Hs.egALIAS2EG)[1]
#   geneid[i,2] <- as.integer(sym)
# }
# na.index2 <- which(is.na(sig.geneid[,2]))
# for (i in na.index2){
#   s <- unlist(strsplit(as.character(sig.geneid[i,][1]), ","))[1]
#   sym <- get(s, org.Hs.egALIAS2EG)[1]
#   sig.geneid[i,2] <- as.integer(sym)
# }

# geneid <- geneid[!duplicated(geneid[,2]), ]
# sig.geneid <- sig.geneid[!duplicated(sig.geneid[,2]), ]
}
\keyword{datasets}