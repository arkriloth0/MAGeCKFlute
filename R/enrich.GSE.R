#' Gene set enrichment analysis
#'
#' A universal gene set enrichment analysis tools
#'
#' @docType methods
#' @name enrich.GSE
#' @rdname enrich.GSE
#' @aliases enrichGSE
#'
#' @param geneList A order ranked numeric vector with geneid as names
#' @param keytype "Entrez", "Ensembl", or "Symbol"
#' @param type Molecular signatures for testing, available datasets include
#' Pathway (KEGG, REACTOME, C2_CP_PID, C2_CP_BIOCARTA, C2_CP_WIKIPATHWAYS, C2_CP_KEGG_MEDICUS), GO (GOBP, GOCC, GOMF),
#' MSIGDB (C1, C2 (C2_CP (C2_CP_PID, C2_CP_BIOCARTA, C2_CP_WIKIPATHWAYS, C2_CP_KEGG_MEDICUS), C2_CGP),
#' C3 (C3_MIR_MIRDB, C3_TFT_GTRD), C4 (C4_CGN, C4_CM, c4_3CA), C5 (C5_GO_BP, C5_GO_CC, C5_GO_MF, C5_HPO), C6, C7 (C7_IMMUNESIGDB, C7_VAX), C8, H)
#' and Complex (CORUM). Any combination of them are also accessible
#' (e.g. 'GOBP+GOMF+KEGG+REACTOME')
#' @param organism 'hsa' or 'mmu'
#' @param pvalueCutoff FDR cutoff
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets for enrichent analysis
#' @param gmtpath The path to customized gmt file
#' @param by One of 'fgsea' or 'DOSE'
#' @param verbose Boolean
#' @param ... Other parameter
#'
#' @return An enrichResult instance
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{enrich.HGT}}
#' @seealso \code{\link{enrich.ORT}}
#' @seealso \code{\link{EnrichAnalyzer}}
#'
#' @examples
#' data(geneList, package = "DOSE")
#' \dontrun{
#'     enrichRes = enrich.GSE(geneList, keytype = "entrez")
#'     head(slot(enrichRes, "result"))
#' }
#' @export

enrich.GSE <- function(geneList,
                       keytype = "Symbol",
                       type = "GOBP",
                       organism = 'hsa',
                       pvalueCutoff = 1,
                       limit = c(2, 100),
                       gmtpath = NULL,
                       # nPerm = 2000,
                       by = "fgsea",
                       verbose = TRUE,
                       ...){
  requireNamespace("clusterProfiler", quietly=TRUE) || stop("Please install the clusterProfiler package")
  requireNamespace("DOSE", quietly=TRUE) || stop("Please install the DOSE package")
  # type[type=="GOBP"] <- "C5_BP"
  # type[type=="GOCC"] <- "C5_CC"
  # type[type=="GOMF"] <- "C5_MF"
  geneList = sort(geneList, decreasing = TRUE)

  ## Prepare gene set annotation
  gene2path = gsGetter(gmtpath, type, limit, organism)
  idx = duplicated(gene2path$PathwayID)
  pathways = data.frame(PathwayID = gene2path$PathwayID[!idx],
                        PathwayName = gene2path$PathwayName[!idx],
                        stringsAsFactors = FALSE)

  ## Gene ID conversion
  if(keytype != "Entrez"){
    allsymbol = names(geneList)
    gene = TransGeneID(allsymbol, keytype, "Entrez", organism = organism)
    idx = is.na(gene) | duplicated(gene)
    geneList = geneList[!idx]
    names(geneList) = gene[!idx]
  }

  ## Enrichment analysis
  len = length(unique(intersect(names(geneList), gene2path$Gene)))
  if(verbose) message("\t", len, " genes are mapped ...")
  enrichedRes = clusterProfiler::GSEA(geneList = geneList, pvalueCutoff = pvalueCutoff,
                     minGSSize = min(limit), maxGSSize = max(limit),
                     TERM2NAME = pathways, by = by,
                     TERM2GENE = gene2path[,c("PathwayID","Gene")],
                     verbose = verbose, ...)
  ## Add enriched gene symbols into enrichedRes table
  if(!is.null(enrichedRes) && nrow(enrichedRes@result)>0){
    colnames(enrichedRes@result)[11] = "geneID"
    # enrichedRes@result = EnrichedFilter(enrichedRes@result)
    geneID = strsplit(enrichedRes@result$geneID, "/")
    allsymbol = TransGeneID(names(geneList), "Entrez", "Symbol",
                            organism = organism)
    geneName = lapply(geneID, function(gid){
      SYMBOL = allsymbol[gid]; paste(SYMBOL, collapse = "/")
    })
    enrichedRes@result$geneName = unlist(geneName)
    enrichedRes@result$Count = unlist(lapply(geneID, length))
    cnames = c("ID", "Description", "NES", "pvalue", "p.adjust",
               "geneID", "geneName", "Count")
    res = enrichedRes@result[, cnames]
    res = res[order(res$pvalue), ]
    enrichedRes@result = res
  }
  return(enrichedRes)
}
