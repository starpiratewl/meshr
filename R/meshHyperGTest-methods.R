setMethod("meshHyperGTest", signature(p="MeSHHyperGParams"),
          function(p) .meshHyperGTestInternal(p) )

.meshHyperGTestInternal <- function(p) {
  ##
  ## MeSH enrichment analysis 
  ##

  ## Map gene ID to MeSH ID through annotation data
  my.keytype <- c("GENEID")
  my.cols <- c("GENEID", "MESHID")
  my.geneids <- as.data.frame(p@geneIds)
  names(my.geneids) <- "GENEID"
  universe.geneids <- as.data.frame(p@universeGeneIds)
  names(universe.geneids) <- "GENEID"

  ## Retive data of specific database
  selectedDatabase <- select(eval(parse(text=p@annotation)), keys = p@database, columns = c("GENEID", "MESHID","MESHCATEGORY"), keytype = "SOURCEDB")
  selectedDatabase <- selectedDatabase[which(selectedDatabase[,3] == p@category), 1:2]

  ## Error against impossible category-database combination was choosed
  if(nrow(selectedDatabase) == 0){
    stop("Impossible MeSH category - database combination was choosed.")
  }

  selected.mesh <- merge(my.geneids, selectedDatabase, "GENEID")
  universe.mesh <- merge(universe.geneids, selectedDatabase, "GENEID")

  selected.table <- table(selected.mesh[,2])
  universe.table <- table(universe.mesh[,2])

  ## Hypergeometric test 
  pvals <- array()
  for (i in 1:length(selected.table)) {
    numWdrawn <- selected.table[i]
    mesh.index <- which(names(selected.table[i])==names(universe.table))
    numW <- universe.table[mesh.index]
    numB <- length(p@universeGeneIds) - numW
    numDrawn <- length(p@geneIds)
    pvals[i] <- phyper(numWdrawn - 1L, numW, numB, numDrawn, lower.tail=FALSE)
  }
  
  ## Multiple testing correction  
  stats <- switch(p@pAdjust, BH = {
    p.adjust(pvals, "BH")
  }, QV = {
    suppressWarnings(fdrtool(pvals, statistic="pvalue", plot=FALSE, verbose=FALSE)$qval)
  }, lFDR = {
    suppressWarnings(fdrtool(pvals, statistic="pvalue", plot=FALSE, verbose=FALSE)$lfdr)
  }, none = pvals)

  ## Retrieve full name of MeSH category 
  mesh.cat <- p@category

  switch(mesh.cat,
    "A" = {
      mesh.full.cat <- "Anatomy"
    }, "B" = {
      mesh.full.cat <- "Organisms"
    }, "C" = {
      mesh.full.cat <- "Diseases"
    }, "D" = {
      mesh.full.cat <- "Chemicals and Drugs"
    }, "E" = {
      mesh.full.cat <- "Analytical, Diagnostic and Therapeutic Techniques and Equipment"
    }, "F" = {
      mesh.full.cat <- "Psychiatry and Psychology"
    }, "G" = {
      mesh.full.cat <- "Phenomena and Processes"
    }, "H" = {
      mesh.full.cat <- "Disciplines and Occupations"
    }, "I" = {
      mesh.full.cat <- "Anthropology, Education, Sociology and Social Phenomena"
    }, "J" = {
      mesh.full.cat <- "Technology and Food and Beverages"
    }, "K" = {
      mesh.full.cat <- "Humanities"
    }, "L" = {
      mesh.full.cat <- "Information Science"
    }, "M" = {
      mesh.full.cat <- "Persons"
    }, "N" = {
      mesh.full.cat <- "Health Care"
    }, "V" = {
      mesh.full.cat <- "Publication Type"
    }, "Z" = {
      mesh.full.cat <- "Geographical Locations"
    }
  )

  ## Choose siginificantly enriched MeSH terms 
  mesh.list <- names(selected.table[which(stats < p@pvalueCutoff)])
  ## Create a data.frame
  tmp.df <- data.frame("MESHID" = mesh.list, "PVALUE" = stats[which(stats < p@pvalueCutoff)], stringsAsFactors=FALSE)
  ## Mapping  MeSH ID to MeSH term via MeSH.db
  mesh.df <- select(MeSH.db, keys=mesh.list, columns=c("MESHID", "MESHTERM", "CATEGORY"), keytype="MESHID")
  # remove categories not specified 
  mesh.df <- mesh.df[mesh.cat==mesh.df[,3],][,c(1,2)]
  ## merge mesh.df and tmp.df
  mesh.df <- merge(mesh.df, tmp.df, by = "MESHID", all=TRUE)
  ## calculate sorted index
  sort.index <- unlist(sort(mesh.df[,3], index.return=TRUE)[2])
  ## sort mesh.df based on p-values
  mesh.df <- mesh.df[sort.index,]
  # remove MeSH terms appearing multiple times within same category 
  mesh.df <- mesh.df[!duplicated(mesh.df[,1]),]
  if (nrow(mesh.df) != 0){
    rownames(mesh.df) <- 1:nrow(mesh.df)
  }
  
  new("MeSHHyperGResult",
      meshCategory=mesh.full.cat,
      meshAnnotation=p@annotation, 
      meshDatabase=p@database, 
      meshIds=mesh.df[,1],
      meshTerms=mesh.df[,2],
      pvalues=mesh.df[,3])
}

