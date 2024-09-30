#' Enrichment analysis for Positive and Negative selection genes
#'
#' Do enrichment analysis for selected genes, in which positive selection and negative selection
#' are termed as Positive and Negative
#'
#' @docType methods
#' @name EnrichAB
#' @rdname EnrichAB
#'
#' @param data A data frame.
#' @param enrich_method One of "ORT" (Over-Representing Test) and "HGT" (HyperGemetric test).
#' @param top An integer, specifying the number of pathways to show.
#' @param limit A two-length vector, specifying the min and max size of pathways
#' for enrichent analysis.
#' @param filename Suffix of output file name.
#' @param out.dir Path to save plot to (combined with filename).
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param units The units of figure size, one of "in", "cm", "mm", "px".
#' @param verbose Boolean
#' @param ... Other available parameters in ggsave.
#'
#' @return A list containing enrichment results for each group genes. This list contains eight
#' items, which contain subitems of \code{gridPlot} and \code{enrichRes}.
#'
#' @author Wubing Zhang
#'
EnrichAB <- function(data,
                     enrich_method = "HGT",
                     top = 10,
                     limit = c(2, 100),
                     filename = NULL, out.dir = ".",
                     width = 6.5, height = 4, units = "in", 
                     verbose = TRUE, ...){

  requireNamespace("clusterProfiler", quietly=TRUE) || stop("Need clusterProfiler package")
  if(verbose) message(Sys.time(), " # Enrichment analysis of Positive and Negative selected genes ...")
  gg = data

  ##=====enrichment for Positive======
  idx1 = gg$group=="top"; genes = gg$GeneID[idx1]
  geneList = gg$Diff[idx1]; names(geneList) = genes
  enrichA = EnrichAnalyzer(geneList = geneList, universe = gg$GeneID, keytype = "entrez",
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           pvalueCutoff = 1, limit = limit)
  if(!is.null(enrichA) && nrow(enrichA@result)>0){
    keggA = enrichA@result[grepl("KEGG", enrichA@result$ID), ]
    gobpA = enrichA@result[grepl("^GO", enrichA@result$ID), ]
    reactomeA = enrichA@result[grepl("REACTOME", enrichA@result$ID), ]
    complexA = enrichA@result[grepl("CORUM", enrichA@result$ID), ]
    keggA = list(enrichRes = keggA, gridPlot = EnrichedView(keggA, top = top, bottom = 0)
                 + labs(title = "KEGG: Positive") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
    gobpA = list(enrichRes = gobpA, gridPlot = EnrichedView(gobpA, top = top, bottom = 0)
               + labs(title = "GOBP: Positive") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
    reactomeA = list(enrichRes = reactomeA, gridPlot = EnrichedView(reactomeA, top = top, bottom = 0)
                     + labs(title = "REACTOME: Positive") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
    complexA = list(enrichRes = complexA, gridPlot = EnrichedView(complexA, top = top, bottom = 0)
                    + labs(title = "Complex: Positive") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
  }else{
    keggA = gobpA = reactomeA = complexA = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  idx2 = gg$group=="bottom"; genes = gg$GeneID[idx2]
  geneList = gg$Diff[idx2]; names(geneList) = genes
  enrichB = EnrichAnalyzer(geneList = geneList, universe = gg$GeneID,
                           method = enrich_method, type = "KEGG+REACTOME+GOBP+Complex",
                           pvalueCutoff = 1, limit = limit, keytype = "entrez", verbose = verbose)
  if(!is.null(enrichB) && nrow(enrichB@result)>0){
    keggB = enrichB@result[grepl("KEGG", enrichB@result$ID), ]
    gobpB = enrichB@result[grepl("^GO", enrichB@result$ID), ]
    reactomeB = enrichB@result[grepl("REACTOME", enrichB@result$ID), ]
    complexB = enrichB@result[grepl("CORUM", enrichB@result$ID), ]
    keggB = list(enrichRes = keggB, gridPlot = EnrichedView(keggB, top = 0, bottom = top)
                 + labs(title = "KEGG: Negative") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
    gobpB = list(enrichRes = gobpB, gridPlot = EnrichedView(gobpB, top = 0, bottom = top)
               + labs(title = "GOBP: Negative") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
    reactomeB = list(enrichRes = reactomeB, gridPlot = EnrichedView(reactomeB, top = 0, bottom = top)
                     + labs(title = "REACTOME: Negative") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
    complexB = list(enrichRes = complexB, gridPlot = EnrichedView(complexB, top = 0, bottom = top)
                    + labs(title = "Complex: Negative") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
  }else{
    keggB = gobpB = reactomeB = complexB = list(enrichRes = NULL, gridPlot = noEnrichPlot())
  }

  if(!is.null(filename)){
    ## Save Positive enrichment results
    if(!is.null(enrichA) && nrow(enrichA@result)>0){
      write.table(keggA$enrichRes, file.path(out.dir,paste0("GroupA_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(keggA$gridPlot, filename=file.path(out.dir,paste0("GroupA_kegg_", filename,".png")),
             units = units, width=6.5, height=4, ...)
      write.table(reactomeA$enrichRes, file.path(out.dir,paste0("GroupA_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactomeA$gridPlot, filename=file.path(out.dir,paste0("GroupA_reactome_", filename,".png")),
             units = units, width=6.5, height=4, ...)
      write.table(gobpA$enrichRes, file.path(out.dir,paste0("GroupA_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobpA$gridPlot, filename=file.path(out.dir,paste0("GroupA_gobp_", filename,".png")),
             units = units, width=6.5, height=4, ...)
      write.table(complexA$enrichRes, file.path(out.dir,paste0("GroupA_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complexA$gridPlot, filename=file.path(out.dir,paste0("GroupA_complex_", filename,".png")),
             units = units, width=6.5, height=4, ...)
    }
    if(!is.null(enrichB) && nrow(enrichB@result)>0){
      write.table(keggB$enrichRes, file.path(out.dir,paste0("GroupB_kegg_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(keggB$gridPlot, filename=file.path(out.dir,paste0("GroupB_kegg_", filename,".png")),
             units = units, width=6.5, height=4, ...)
      write.table(reactomeB$enrichRes, file.path(out.dir,paste0("GroupB_reactome_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(reactomeB$gridPlot, filename=file.path(out.dir,paste0("GroupB_reactome_", filename,".png")),
             units = units, width=6.5, height=4, ...)
      write.table(gobpB$enrichRes, file.path(out.dir,paste0("GroupB_gobp_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(gobpB$gridPlot, filename=file.path(out.dir,paste0("GroupB_gobp_", filename,".png")),
             units = units, width=6.5, height=4, ...)
      write.table(complexB$enrichRes, file.path(out.dir,paste0("GroupB_complex_",filename,".txt")),
                  sep="\t", row.names = FALSE, col.names = TRUE,quote=FALSE)
      ggsave(complexB$gridPlot, filename=file.path(out.dir,paste0("GroupB_complex_", filename,".png")),
             units = units, width=6.5, height=4, ...)
    }
  }
  return(list(keggA=keggA, gobpA=gobpA, reactomeA=reactomeA, complexA=complexA,
              keggB=keggB, gobpB=gobpB, reactomeB=reactomeB, complexB=complexB))
}
