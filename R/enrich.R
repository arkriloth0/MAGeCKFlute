enrich <- function(geneList,
                   listType = NULL,
                   annoType = "Entrez",
                   organism = "hsa",
                   TERM2GENE = NULL,
                   TERM2NAME = NULL,
                   pvalueCutoff = 1,
                   limit = c(2, 100),
                   universe=NULL,
                   gmtpath = NULL,
                   verbose = TRUE,
                   ...){
  ## Gene ID conversion
  gene = names(geneList)
  
  if(!is.null(listType)){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, listType, annoType, organism = organism)
    idx = duplicated(gene)|is.na(gene)
    geneList = geneList[!idx]; names(geneList) = gene[!idx]
  }
  if(!is.null(universe)){
    if(!is.null(listType)){
      universe = TransGeneID(universe, listType, annoType, organism = organism)
    }
    universe = universe[!is.na(universe)]
  }else{
    universe = unique(c(gene, TERM2GENE[[2]]))
  }
  gene = names(geneList)
  
  ## Enrichment analysis
  len = length(unique(intersect(gene, TERM2GENE[[2]])))
  if(verbose) message("\t", len, " genes are mapped ...")
  enrichedRes = clusterProfiler::enricher(gene, universe = universe,
                                          minGSSize = 0, maxGSSize = max(limit),
                                          TERM2NAME = TERM2NAME, pvalueCutoff = pvalueCutoff,
                                          TERM2GENE = TERM2GENE)
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    res = enrichedRes@result[enrichedRes@result$p.adjust<=pvalueCutoff, ]
    res = res[order(res$pvalue), ]
    enrichedRes@result = res
  }
  ## Add enriched gene symbols into enrichedRes table
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    allsymbol = TransGeneID(gene, annoType, "Symbol", organism = organism)
    geneID = strsplit(enrichedRes@result$geneID, split = "/")
    geneName = lapply(geneID, function(gid){
      SYMBOL = allsymbol[gid]; paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
    enrichedRes@result$NES = as.vector(sapply(enrichedRes@result$geneID, function(x){
      enrichGenes = unlist(strsplit(x, split = "/"))
      NES = mean(geneList[enrichGenes]) * length(enrichGenes)^0.6
      return(NES)
    }))
    observed = sapply(strsplit(enrichedRes@result$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    expected = sapply(strsplit(enrichedRes@result$BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    enrichedRes@result$strength = log10(observed/expected)
    idx = c("ID", "Description", "NES", "strength", "pvalue", "p.adjust",
            "GeneRatio", "BgRatio", "geneID", "geneName", "Count")
    enrichedRes@result = enrichedRes@result[, idx]
  }
  return(enrichedRes)
}