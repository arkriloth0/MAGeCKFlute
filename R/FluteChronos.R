#' Downstream analysis based on MAGeCK-MLE result
#'
#' Integrative analysis pipeline using the gene summary table in MAGeCK MLE results
#'
#' @docType methods
#' @name FluteChronos
#' @rdname FluteChronos
#' @aliases flutechronos
#'
#' @param gene_summary A data frame containing a column labeled `Gene` containing HGNC symbols or EntrezIDs, 
#' with other columns labeled with sample names containing scaled Chronos gene effect scores. 
#' @param treatname A character vector, specifying the names of treatment samples (column names in `gene_summary`).
#' @param ctrlname A character vector, specifying the names of control samples (column names in `gene_summary`).
#' If there is no controls in your CRISPR screen, you can specify "Depmap" as ctrlname and set `incorporateDepmap=TRUE`.
#' @param lab_treat A character, for labelling the treatment condition on graphs.
#' @param lab_ctrl A character, for labelling the control condition on graphs.
#' @param reps Boolean, indicating whether `gene_summary` contains data of replicates in separate columns denoted by a suffix of "_R1", "_R2", "_R3", etc.
#' @param keytype "Entrez" or "Symbol".
#' @param organism "hsa" or "mmu".
#' @param incorporateDepmap Boolean, indicating whether incorporate Depmap data into analysis.
#' @param cell_lines A character vector, specifying the cell lines in Depmap to be considered.
#' @param lineages A character vector, specifying the lineages in Depmap to be considered.
#' @param omitEssential Boolean, indicating whether omit common essential genes from the downstream analysis.
#' @param EG A character vector of symbols of common essential genes to be omitted. 
#' If `Null` (default), common essential genes built into MAGeCKFlute package will be used.
#' @param FalsePositive A character vector of symbols of unexpressed genes to be omitted. Default: `NULL`.
#' @param top An integer, specifying the number of top selected genes to be labeled
#' in rank figure and the number of top pathways to be shown.
#' @param toplabels A character vector, specifying interested genes to be labeled in rank figure.
#' @param repDiff_cutoff Numeric, Replicates with difference in gene effect scores greater than the cutoff will be filtered out before visualisation and enrichment analyses. 
#' @param scale_cutoff Boolean or numeric, specifying how many standard deviation will be used as cutoff.
#' @param fdr_cutoff Boolean or numeric, specifying the fdr threshold for identifying differentially dependent genes.
#' @param limit A two-length vector, specifying the minimal and maximal size of gene sets for enrichent analysis.
#' @param enrich_method One of "ORT"(Over-Representing Test) and "HGT"(HyperGemetric test).
#' @param proj A character, indicating the prefix of output file name, which can't contain special characters.
#' @param width The width of summary pdf in inches.
#' @param height The height of summary pdf in inches.
#' @param figwidth The width of png figures.
#' @param units The units of png figure sizes, one of "in", "cm", "mm", "px".
#' @param dpi dpi of png figures, default = 1200.
#' @param outdir Output directory on disk.
#' @param pathview.top Integer, specifying the number of pathways for pathview visualization.
#' @param verbose Boolean
#'
#' @author Wubing Zhang & Allan Lui
#'
#' @return All of the pipeline results is output into the \code{out.dir}/MAGeCKFlute_\code{proj},
#' which includes a pdf file and many folders. The pdf file 'FluteMLE_\code{proj}.pdf' is the
#' summary of pipeline results. For each section in this pipeline, figures and useful data are
#' outputed to corresponding subfolders.
#' \itemize{
#'   \item {QC}: {Quality control}
#'   \item {Selection}: {Positive selection and negative selection.}
#'   \item {Enrichment}: {Enrichment analysis for positive and negative selection genes.}
#'   \item {PathwayView}: {Pathway view for top enriched pathways.}
#' }
#'
#' @details MAGeCK-MLE can be used to analyze screen data from multi-conditioned experiments. MAGeCK-MLE
#' also normalizes the data across multiple samples, making them comparable to each other. The most important
#' ouput of MAGeCK MLE is `gene_summary` file, which includes the beta scores of multiple conditions and
#' the associated statistics. The `beta score`  for each gene describes how the gene is selected: a positive
#' beta score indicates a positive selection, and a negative beta score indicates a negative selection.
#'
#' The downstream analysis includes identifying essential, non-essential, and target-associated genes,
#' and performing biological functional category analysis and pathway enrichment analysis of these genes.
#' The function also visualizes genes in the context of pathways to benefit users exploring screening data.
#'
#'
#' @seealso \code{\link{FluteRRA}}
#'
#' @examples
#' file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#' "testdata/mle.gene_summary.txt")
#' \dontrun{
#'   # functional analysis for MAGeCK MLE results
#'   FluteMLE(file3, treatname = "Pmel1", ctrlname = "Pmel1_Ctrl", proj = "Pmel1")
#' }
#'
#' @import data.table dplyr ggplot2 stats grDevices utils gridExtra grid 
#' @export

FluteChronos <- function(gene_summary, treatname, ctrlname = "Depmap",
                         lab_treat = "Treatment", lab_ctrl = "Control", reps = F,
                         keytype = "Symbol", organism = "hsa", # Input dataset
                         incorporateDepmap = FALSE,
                         cell_lines = NA, lineages = "All",
                         omitEssential = TRUE, EG = NULL,
                         FalsePositive = NULL,
                         top = 10, toplabels = NA,
                         repDiff_cutoff = 0.5, scale_cutoff = 2, fdr_cutoff = 0.1,
                         limit = c(2,500),
                         enrich_method = "ORT", proj = NA,
                         width = 10, height = 7, figwidth = 13.5, units = "cm",  dpi = 1200,
                         outdir = ".", pathview.top = 4, verbose = TRUE){
  ## Prepare the running environment ##
  {
    message(Sys.time(), " # Create output dir and pdf file...")
    if(omitEssential | !is.null(FalsePositive)){
      filter = "filter"
    } else {filter = "noFilter"}
    outdir = file.path(outdir, paste0("MAGeCKFlute_Chronos_", proj))
    dir.create(file.path(outdir), showWarnings = FALSE)
    output_pdf = paste0("FluteChronos_", proj, "_", filter, ".pdf")
    pdf(file.path(outdir, output_pdf), width = width, height = height)
  }
  
  ## Beta Score Preparation ##
  {
    beta = gene_summary
    if(tolower(keytype)!="entrez"){
      beta$EntrezID = TransGeneID(beta$Gene, keytype, "Entrez", organism = organism)
    }else{
      beta$EntrezID = beta$Gene
    }
    if(tolower(keytype)!="symbol"){
      beta$Symbol = TransGeneID(beta$Gene, keytype, "Symbol", organism = organism)
    }else{
      beta$Symbol = beta$Gene
    }
    if(organism != "hsa"){
      beta$HumanGene = TransGeneID(beta$Symbol, "Symbol", "Symbol",
                                   fromOrg = organism, toOrg = "hsa")
      beta$GeneID = TransGeneID(beta$HumanGene, "Symbol", "Entrez")
    }else{
      beta$HumanGene = beta$Symbol
      beta$GeneID = beta$EntrezID
    }
    
    message(Sys.time(), " # Transform id to official human gene name ...")
    idx1 = is.na(beta$EntrezID)
    idx2 = !is.na(beta$EntrezID) & duplicated(beta$EntrezID)
    idx = idx1|idx2
    if(sum(idx1)>0) message(sum(idx1), " genes fail to convert into Entrez IDs: ",
                            paste0(beta$Gene[idx1], collapse = ", "))
    if(sum(idx2)>0) message(sum(idx2), " genes have duplicate Entrez IDs: ",
                            paste0(beta$Gene[idx2], collapse = ", "))
    
    dd = beta[!idx, ]
    if(incorporateDepmap | "Depmap"%in%ctrlname)
      dd = IncorporateDepmap(dd, symbol = "HumanGene", cell_lines = cell_lines,
                             lineages = lineages)
    if(!all(c(ctrlname, treatname) %in% colnames(dd)))
      stop("Sample name doesn't match !!!")
  }
  ## Distribution of beta scores ##
  {
    ## All genes ##
    outputDir1 = file.path(outdir, paste0("QC_",filter))
    dir.create(outputDir1, showWarnings = FALSE)
    idx_distr = c(ctrlname, treatname)
    if (length(idx_distr) <= 4) {
      P1 = ViolinView(dd[, idx_distr], ylab = "Gene effect score", main = "All genes",
                      filename = paste0(outputDir1,"/", "ViolinView_all_", filter, ".png"),
                      width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
    } else  {
      message(Sys.time(), " # >4 reps to be graphed in ViolinView, figure will be double specified width")
      P1 = ViolinView(dd[, idx_distr], ylab = "Gene effect score", main = "All genes") +
        theme(axis.text.x = element_blank())
      ggsave(paste0(outputDir1,"/", "ViolinView_all_", filter, ".png"), 
             P1 + theme(legend.position = "right"),
             width = figwidth*4/3, height = figwidth*3/4, units = units, dpi = dpi)
    }
    # P1 = ViolinView(dd[, idx_distr], ylab = "Gene effect score", main = "All genes",
    #                 filename = paste0(outputDir1,"/", "ViolinView_all_", filter, ".png"))
    P2 = DensityView(dd[, idx_distr], xlab = "Gene effect score", main = "All genes",
                     filename = paste0(outputDir1, "/", "DensityView_all_", filter, ".png"),
                     width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
    P3 = ConsistencyView(dd, ctrlname, treatname, lab_treat = lab_treat, lab_ctrl = lab_ctrl, main = "All genes",
                         filename = paste0(outputDir1, "/", "Consistency_all_", filter, ".png"),
                         width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
    P4 = MAView(dd, ctrlname, treatname, main = "All genes",
                filename = paste0(outputDir1, "/", "MAView_all_", filter, ".png"),
                width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
    gridExtra::grid.arrange(P1, P2, P3, P4, ncol = 2)
    
    ## Essential genes ##
    Zuber_Essential = readRDS(file.path(system.file("extdata", package = "MAGeCKFlute"),
                                        "Zuber_Essential.rds"))
    if(is.null(EG))
      idx = toupper(dd$HumanGene) %in% toupper(Zuber_Essential$GeneSymbol)
    else
      idx = toupper(dd$Gene) %in% toupper(EG)
    
    if(sum(idx) > 10){
      if (length(idx_distr) <= 4) {
        P1 = ViolinView(dd[idx, idx_distr], ylab = "Gene effect score", main = "Essential genes",
                        filename = paste0(outputDir1, "/", "ViolinView_posctrl_", filter, ".png"),
                        width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
      } else  {
        P1 = ViolinView(dd[idx, idx_distr], ylab = "Gene effect score", main = "Essential genes") +
          theme(axis.text.x = element_blank())
        ggsave(paste0(outputDir1,"/", "ViolinView_posctrl_", filter, ".png"), 
               P1 + theme(legend.position = "right"),
               width = figwidth*4/3, height = figwidth*3/4, units = units, dpi = dpi)
      }
      # P1 = ViolinView(dd[idx, idx_distr], ylab = "Gene effect score", main = "Essential genes",
      #                 filename = paste0(outputDir1, "/", "ViolinView_posctrl_", filter, ".png"))
      P2 = DensityView(dd[idx, idx_distr], xlab = "Gene effect score", main = "Essential genes",
                       filename = paste0(outputDir1, "/", "DensityView_posctrl_", filter, ".png"),
                       width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
      P3 = ConsistencyView(dd[idx, ], ctrlname, treatname, lab_treat = lab_treat, lab_ctrl = lab_ctrl, main = "Essential genes",
                           filename = paste0(outputDir1, "/", "Consistency_posctrl_", filter, ".png"),
                           width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
      P4 = MAView(dd[idx, ], ctrlname, treatname, main = "Essential genes",
                  filename = paste0(outputDir1, "/", "MAView_posctrl_", filter, ".png"),
                  width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
      gridExtra::grid.arrange(P1, P2, P3, P4, ncol = 2)
    }
  }
  
  
  
  ## Combine replicates (and hypothesis testing) ##
  {
    dd$Control = Biobase::rowMedians(as.matrix(dd[, ctrlname, drop = FALSE]))
    dd$Treatment = Biobase::rowMedians(as.matrix(dd[, treatname, drop = FALSE]))
    dd$Diff = dd$Treatment - dd$Control
    x_cut = c(-CutoffCalling(dd$Control, scale_cutoff),
              CutoffCalling(dd$Control, scale_cutoff))
    y_cut = c(-CutoffCalling(dd$Treatment, scale_cutoff),
              CutoffCalling(dd$Treatment, scale_cutoff))
    intercept = c(-CutoffCalling(dd$Diff, scale_cutoff),
                  CutoffCalling(dd$Diff, scale_cutoff))
    
    ## Hypothesis testing ##
    {
      FDR = NULL
      if (reps == T) {
        condition_treat <- str_split_i(treatname, "_(?!.*_)",1) %>% unique()
        condition_ctrl <- str_split_i(ctrlname, "_(?!.*_)",1) %>% unique()
        
        rep <- beta %>% 
          dplyr::select(-any_of(c("EntrezID","Symbol","HumanGene","GeneID"))) %>%
          pivot_longer(cols = !Gene, 
                       names_to = c("condition","rep"), 
                       names_sep = "_(?!.*_)", 
                       values_to = "geneEffect") %>%
          dplyr::filter(condition %in% c(condition_treat,condition_ctrl)) %>%
          as.data.table()
        

        if (length(unique(rep$condition)) >= 2 & length(unique(rep$rep)) == 2) {
          message(Sys.time(), 
                  " # >=2 conditions & 2 replicates per condition detected, performing significance testing by t-test comparing difference between selected conditions to difference between replicates")
          
          # calculate deviation from median of replicate gene effect scores for each gene and condition
          median <- rep %>% 
            group_by(Gene, condition) %>%
            summarise(median = median(geneEffect))
          rep.wide <- rep %>% 
            pivot_wider(names_from = rep, values_from = geneEffect) %>%
            as.data.table()
          meddev <- rep.wide %>%
            left_join(median) %>%
            mutate(across(!Gene & !condition & !median, ~ .x - median)) %>%
            dplyr::select(!median) %>% 
            pivot_longer(cols = !Gene & !condition,
                         names_to = "rep",
                         values_to = "medDev") %>%
            nest(medDev = medDev, .by = Gene)
          
          dd <- dd %>% left_join(meddev, by = "Gene") %>% group_by(Gene) %>%
            mutate(p.value = t.test(medDev[[1]]$medDev, mu = Diff)$p.value) %>%
            ungroup() %>%
            mutate(FDR = p.adjust(p.value, method = "fdr"),
                   logFDR = -log10(FDR),
                   stat = Diff*logFDR) %>%
            dplyr::select(!medDev)
          
          # identify and filter genes with >0.5 difference between reps
          rep.filter <- rep.wide %>%
            dplyr::filter(condition %in% c(condition_treat,condition_ctrl)) %>%
            mutate(diff = unlist(abs(rep.wide[,3]-rep.wide[,4]))) %>%
            dplyr::filter(condition %in% c(condition_treat,condition_ctrl)) %>%
            group_by(Gene) %>%
            summarise(maxRepDiff = max(diff))
          dd <- dd %>% left_join(rep.filter, by = "Gene") %>%
            mutate(filter.repDiff = maxRepDiff >= repDiff_cutoff)
          
          FDR = "FDR_"
          
        } else if (length(unique(rep$rep)) >= 3) {
          message(Sys.time(), 
                  paste0(" # ", length(unique(rep$rep))," replicates & ", length(unique(rep$condition))," conditions detected, performing significance testing by t-test comparing selected conditions"))
          
          test <- rep %>%
            mutate(condition = case_when(condition %in% condition_treat ~ "treat",
                                         condition %in% condition_ctrl ~ "ctrl",
                                         .default = NA)) %>%
            dplyr::filter(!is.na(condition)) %>%
            pivot_wider(id_cols = Gene, names_from = condition, values_from = geneEffect, values_fn = list) %>%
            group_by(Gene) %>%
            summarise(p.value = t.test(treat[[1]],ctrl[[1]])$p.value) %>% 
            mutate(FDR = p.adjust(p.value, method = "fdr"),
                   logFDR = -log10(FDR), 
                   stat = Diff*logFDR)
          
          dd <- dd %>% left_join(test, by = "Gene") 
          
          # identify and filter genes with >0.5 difference between reps
          rep.filter <- rep %>% 
            dplyr::filter(condition %in% c(condition_treat,condition_ctrl)) %>%
            group_by(Gene, condition) %>%
            summarise(RepDiff = max(geneEffect) - min(geneEffect)) %>%
            group_by(Gene) %>%
            summarise(maxRepDiff = max(RepDiff))
          dd <- dd %>% left_join(rep.filter, by = "Gene") %>%
            mutate(filter.repDiff = maxRepDiff >= repDiff_cutoff)
          
          FDR = "FDR_"
          
        } else {
          message(Sys.time(), 
                  paste0(" # ", length(unique(rep$rep))," replicates & ", length(unique(rep$condition))," conditions detected, 
                insufficient for significance testing"))
        }
      } 
    }
    
    dd$Rank = rank(dd$Diff)
    write.table(dd, paste0(outputDir1, "/", proj, "_processed_data.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    if(omitEssential){
      if(!is.null(EG)){
        dd = dplyr::filter(dd, !(HumanGene %in% EG) )
      } else {
        dd = OmitCommonEssential(dd, symbol = "HumanGene", dependency = -0.4)
      }
    }
    if(!is.null(FalsePositive)) {
      dd = dplyr::filter(dd, !(HumanGene %in% FalsePositive) )
    }
    if(omitEssential | !is.null(FalsePositive)){
      write.table(dd, paste0(outputDir1, "/", proj, "_cleaned_data.txt"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
    dd <- dd %>% dplyr::filter(!filter.repDiff)
    dd$RandomIndex = sample(1:nrow(dd), nrow(dd))
    dd$Gene = dd$Symbol
  }
  
  ### Diff values ###
  {
    # ## Differentially dependent genes ##
    # {
    #   outputDir2 = file.path(outdir, paste0("Selection_",filter))
    #   dir.create(outputDir2, showWarnings = FALSE)
    #   
    #   p1 = ScatterView(dd, "Control", "Treatment", groups = c("top", "bottom"),
    #                    groupnames = c("GroupA", "GroupB"), intercept = intercept,
    #                    xlab = lab_ctrl, ylab = lab_treat)
    #   ggsave(paste0(outputDir2, "/", "ScatterView_TreatvsCtrl_", filter, ".png"),
    #          p1, width = 4, height = 3)
    #   write.table(p1$data, paste0(outputDir2, "/", "Data_ScatterView_TreatvsCtrl_", filter, ".txt"),
    #               sep = "\t", row.names = FALSE, quote = FALSE)
    #   p2 = ScatterView(dd, x = "Rank", y = "Diff", label = "Symbol",
    #                    groups = c("top", "bottom"), groupnames = c("GroupA", "GroupB"),
    #                    top = top, y_cut = intercept)
    #   ggsave(paste0(outputDir2, "/", "RankView_Treat-Ctrl_", filter, ".png"),
    #          p2, width = 3, height = 5)
    #   p3 = ScatterView(dd[dd$Diff>0, ], x = "RandomIndex", y = "Diff", label = "Symbol",
    #                    y_cut = intercept, groups = "top", groupnames = c("GroupA"), top = top)
    #   ggsave(paste0(outputDir2, "/", "ScatterView_Treat-Ctrl_Positive_", filter, ".png"),
    #          p3, width = 4, height = 3)
    #   p4 = ScatterView(dd[dd$Diff<0, ], x = "RandomIndex", y = "Diff", label = "Symbol",
    #                    y_cut = intercept, groups = "bottom", groupnames = c("GroupB"), top = top)
    #   ggsave(paste0(outputDir2, "/", "ScatterView_Treat-Ctrl_Negative_", filter, ".png"),
    #          p4, width = 4, height = 3)
    #   
    #   gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
    # }
    # 
    # ## Enrichment analysis of negative and positive selected genes ##
    # {
    #   outputDir3 = file.path(outdir, paste0("/Enrichment_", filter))
    #   outputDir4 = file.path(outdir, paste0("/PathwayView_", filter))
    #   dir.create(outputDir3, showWarnings=FALSE)
    #   dir.create(outputDir4, showWarnings=FALSE)
    #   
    #   E1 = EnrichAB(p1$data, enrich_method = enrich_method,
    #                 organism = organism, limit = limit, top = top,
    #                 filename = proj, out.dir = outputDir3)
    #   # EnrichedView
    #   gridExtra::grid.arrange(E1$keggA$gridPlot, E1$reactomeA$gridPlot, E1$gobpA$gridPlot, E1$complexA$gridPlot, ncol = 2)
    #   gridExtra::grid.arrange(E1$keggB$gridPlot, E1$reactomeB$gridPlot, E1$gobpB$gridPlot, E1$complexB$gridPlot, ncol = 2)
    #   
    #   # Pathway view for top 4 pathway
    #   if(!is.null(E1$keggA$enrichRes) && nrow(E1$keggA$enrichRes)>0)
    #     arrangePathview(dd, gsub("KEGG_", "", E1$keggA$enrichRes$ID),
    #                     top = pathview.top, ncol = 2, title="Group A",
    #                     organism=organism, output=outputDir4)
    #   if(!is.null(E1$keggB$enrichRes) && nrow(E1$keggB$enrichRes)>0)
    #     arrangePathview(dd, gsub("KEGG_", "", E1$keggB$enrichRes$ID),
    #                     top = pathview.top, ncol = 2,
    #                     title="Group B", organism=organism,
    #                     output=outputDir4)
    # }
  }
  
  
  ### FDR ###
  {
    ## Differentially dependent genes ##
    {
      outputDir2 = file.path(outdir, paste0("Selection_",FDR,filter))
      dir.create(outputDir2, showWarnings = FALSE)
      omit = NULL
      
      if (!is.null(FDR)) {
        p0 = ScatterView(dd, x = "Diff", y = "logFDR", label = "Symbol",
                         group_col = c("#8da0cb", "#e78ac3"),
                         groupnames = c("GroupA", "GroupB"),
                         x_cut = intercept, y_cut = -log10(fdr_cutoff),
                         model = "volcano", top = top, toplabels = toplabels,
                         display_cut = TRUE, ylab = expression(paste("- log"[10]," FDR")), 
                         max.overlaps = 8)
        ggsave(paste0(outputDir2, "/", "VolcanoView_", filter, ".png"), p0, 
               width = figwidth*4/3, height = figwidth, units = units, dpi = dpi)
        gridExtra::grid.arrange(p0, ncol = 1)
        p0$data <- p0$data %>% as.data.table()
        p0$data[group == "topleft"]$group = "bottom"
        p0$data[group == "topright"]$group = "top"
        write.table(p0$data, paste0(outputDir2, "/", "Data_VolcanoView_", filter, ".txt"),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        dd <- dd %>% as.data.table()
        omit =  dd[logFDR < -log10(fdr_cutoff)]$Symbol
      }
      
      p1 = ScatterView(dd, "Control", "Treatment", groups = c("top", "bottom"),
                       groupnames = c("GroupA", "GroupB"), intercept = intercept,
                       omit = omit,
                       xlab = lab_ctrl, ylab = lab_treat)
      ggsave(paste0(outputDir2, "/", "ScatterView_TreatvsCtrl_", FDR, filter, ".png"),
             p1, width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
      if (is.null(FDR)) {
        write.table(p1$data, paste0(outputDir2, "/", "Data_ScatterView_TreatvsCtrl_", filter, ".txt"),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        p0 = p1
      }
      p2 = ScatterView(dd, x = "Rank", y = "Diff", label = "Symbol",
                       groups = c("top", "bottom"), groupnames = c("GroupA", "GroupB"),
                       y_cut = intercept, top = top, toplabels = toplabels, omit = omit, 
                       max.overlaps = 10, force = 10, max.time = 1, box.padding = 0.6)
      ggsave(paste0(outputDir2, "/", "RankView_Treat-Ctrl_", FDR, filter, ".png"),
             p2, width = figwidth*2/3, height = figwidth, units = units, dpi = dpi)
      p3 = ScatterView(dd[dd$Diff>0, ], x = "RandomIndex", y = "Diff", label = "Symbol",
                       y_cut = intercept, groups = "top", groupnames = c("GroupA"), top = top, omit = omit) + 
        theme_bw(base_size = 14) + 
        theme(panel.grid = element_blank(), legend.position = "none")
      ggsave(paste0(outputDir2, "/", "ScatterView_Treat-Ctrl_Positive_", FDR, filter, ".png"),
             p3, width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
      p4 = ScatterView(dd[dd$Diff<0, ], x = "RandomIndex", y = "Diff", label = "Symbol",
                       y_cut = intercept, groups = "bottom", groupnames = c("GroupB"), top = top, omit = omit) + 
        theme_bw(base_size = 14) + 
        theme(panel.grid = element_blank(), legend.position = "none")
      ggsave(paste0(outputDir2, "/", "ScatterView_Treat-Ctrl_Negative_", FDR, filter, ".png"),
             p4, width = figwidth, height = figwidth*3/4, units = units, dpi = dpi)
      
      gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
      
      
      # p1 = ScatterView(dd, "Control", "Treatment", groups = c("top", "bottom"),
      #                  groupnames = c("GroupA", "GroupB"), intercept = intercept_fdr, size = "logFDR",
      #                  toplabels = dd %>% dplyr::filter(Rank <= top | Rank > nrow(.)-top) %>% .$Symbol %>% c(toplabels),
      #                  label = "Symbol",
      #                  xlab = ctrlname, ylab = treatname
      # )
      # ggsave(paste0(outputDir5, "/", "ScatterView_TreatvsCtrl_FDR_", filter, ".png"),
      #        p1 + theme(legend.position = "right") +
      #          guides(colour = guide_legend(override.aes = list(size = 5))),
      #        width = 8, height = 6)
      # write.table(p1$data, paste0(outputDir5, "/", "Data_ScatterView_TreatvsCtrl_FDR_", filter, ".txt"),
      #             sep = "\t", row.names = FALSE, quote = FALSE)
      # p2 = ScatterView(dd, x = "Rank", y = "zscore_diff", label = "Symbol",
      #                  groups = c("top", "bottom"), groupnames = c("GroupA", "GroupB"),
      #                  top = top, y_cut = c(-intercept_zscore,intercept_zscore), ylab = "Diff Z-Score")
      # ggsave(paste0(outputDir5, "/", "RankView_Treat-Ctrl_Zscore_", filter, ".png"),
      #        p2, width = 3, height = 5)
      
      
    }
    
    ## Enrichment analysis of negative and positive selected genes ##
    {
      outputDir3 = file.path(outdir, paste0("/Enrichment_", FDR, filter))
      outputDir4 = file.path(outdir, paste0("/PathwayView_", FDR, filter))
      dir.create(outputDir3, showWarnings=FALSE)
      dir.create(outputDir4, showWarnings=FALSE)
      
      E1 = EnrichAB(p0$data, enrich_method = enrich_method,
                    # organism = organism, 
                    limit = limit, top = top,
                    filename = proj, out.dir = outputDir3,
                    # width = figwidth, height = figwidth*3/4, units = units, dpi = dpi
      )
      # EnrichedView
      gridExtra::grid.arrange(E1$keggA$gridPlot, E1$reactomeA$gridPlot, E1$gobpA$gridPlot, E1$complexA$gridPlot, ncol = 2)
      gridExtra::grid.arrange(E1$keggB$gridPlot, E1$reactomeB$gridPlot, E1$gobpB$gridPlot, E1$complexB$gridPlot, ncol = 2)
      
      # Pathway view for top 4 pathway
      if(!is.null(E1$keggA$enrichRes) && nrow(E1$keggA$enrichRes)>0)
        arrangePathview(dd, gsub("KEGG_", "", E1$keggA$enrichRes$ID),
                        top = pathview.top, ncol = 2, title="Group A",
                        organism=organism, output=outputDir4)
      if(!is.null(E1$keggB$enrichRes) && nrow(E1$keggB$enrichRes)>0)
        arrangePathview(dd, gsub("KEGG_", "", E1$keggB$enrichRes$ID),
                        top = pathview.top, ncol = 2,
                        title="Group B", organism=organism,
                        output=outputDir4)
    }
  }
  
  ## GSEA analysis ##
  {
    # dd = dd[!(is.na(dd$HumanGene)|duplicated(dd$HumanGene)), ]
    universe = dd$HumanGene
    geneList = dd$Diff; names(geneList) = dd$HumanGene
    enrichRes = EnrichAnalyzer(geneList=geneList, universe=universe,
                               organism="hsa", pvalueCutoff=1,
                               method = "GSEA", limit = limit,
                               type = "KEGG+REACTOME+GOBP+Complex+H")
    kegg = enrichRes@result[grepl("KEGG", enrichRes@result$ID), ]
    gobp = enrichRes@result[grepl("^GO", enrichRes@result$ID), ]
    reactome = enrichRes@result[grepl("REACTOME", enrichRes@result$ID), ]
    complex = enrichRes@result[grepl("CPX|CORUM", enrichRes@result$ID), ]
    hallmark = enrichRes@result[grepl("HALLMARK", enrichRes@result$ID), ]
    hallmark = list(enrichRes = hallmark, gridPlot = EnrichedView(hallmark, x = "LogFDR", rank_by = "NES")
                    + labs(title = "GSEA - Hallmark") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
    
    gsea.pos = enrichRes@result[enrichRes@result$NES>0, ]
    if(!is.null(gsea.pos) && nrow(gsea.pos)>0){
      keggA = gsea.pos[grepl("KEGG", gsea.pos$ID), ]
      gobpA = gsea.pos[grepl("^GO", gsea.pos$ID), ]
      reactomeA = gsea.pos[grepl("REACTOME", gsea.pos$ID), ]
      complexA = gsea.pos[grepl("CPX|CORUM", gsea.pos$ID), ]
      keggA = list(enrichRes = keggA, gridPlot = EnrichedView(keggA, top = top, bottom = 0)
                   + labs(title = "GSEA - KEGG: Positive") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
      reactomeA = list(enrichRes = reactomeA, gridPlot = EnrichedView(reactomeA, top = top, bottom = 0)
                       + labs(title = "GSEA - REACTOME: Positive") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
      gobpA = list(enrichRes = gobpA, gridPlot = EnrichedView(gobpA, top = top, bottom = 0)
                   + labs(title = "GSEA - GOBP: Positive") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
      complexA = list(enrichRes = complexA, gridPlot = EnrichedView(complexA, top = top, bottom = 0)
                      + labs(title = "GSEA - Complex: Positive") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
    }else{
      keggA = gobpA = reactomeA = complexA = list(enrichRes = NULL, gridPlot = noEnrichPlot())
    }
    gsea.neg = enrichRes@result[enrichRes@result$NES<0, ]
    if(!is.null(gsea.neg) && nrow(gsea.neg)>0){
      keggB = gsea.neg[grepl("KEGG", gsea.neg$ID), ]
      gobpB = gsea.neg[grepl("^GO", gsea.neg$ID), ]
      reactomeB = gsea.neg[grepl("REACTOME", gsea.neg$ID), ]
      complexB = gsea.neg[grepl("CPX|CORUM", gsea.neg$ID), ]
      keggB = list(enrichRes = keggB,
                   gridPlot = EnrichedView(keggB, top = 0, bottom = top)
                   + labs(title = "GSEA - KEGG: Negative") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
      reactomeB = list(enrichRes = reactomeB,
                       gridPlot = EnrichedView(reactomeB, top = 0, bottom = top)
                       + labs(title = "GSEA - REACTOME: Negative") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
      gobpB = list(enrichRes = gobpB,
                   gridPlot = EnrichedView(gobpB, top = 0, bottom = top)
                   + labs(title = "GSEA - GOBP: Negative") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
      complexB = list(enrichRes = complexB,
                      gridPlot = EnrichedView(complexB, top = 0, bottom = top)
                      + labs(title = "GSEA - Complex: Negative") + theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()))
    }else{
      keggB = gobpB = reactomeB = complexB = list(enrichRes = NULL, gridPlot = noEnrichPlot())
    }
    gridExtra::grid.arrange(hallmark$gridPlot, ncol = 1)
    gridExtra::grid.arrange(keggA$gridPlot, reactomeA$gridPlot, gobpA$gridPlot, complexA$gridPlot, ncol = 2)
    gridExtra::grid.arrange(keggB$gridPlot, reactomeB$gridPlot, gobpB$gridPlot, complexB$gridPlot, ncol = 2)
    
    ## Save enrichment results ##
    write.table(kegg, file.path(outputDir3, "GSEA_kegg.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    write.table(gobp, file.path(outputDir3, "GSEA_gobp.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    write.table(reactome, file.path(outputDir3, "GSEA_reactome.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    write.table(complex, file.path(outputDir3, "GSEA_complex.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    write.table(hallmark$enrichRes, file.path(outputDir3, "GSEA_hallmark.txt"),
                sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
    
    ggsave(hallmark$gridPlot, filename=file.path(outputDir3, "GSEA_Positive_hallmark.png"),
           units = "in", width=13, height=8)
    if(!is.null(gsea.pos) && nrow(gsea.pos)>0){
      ggsave(keggA$gridPlot, filename=file.path(outputDir3, "GSEA_Positive_kegg.png"),
             units = "in", width=6.5, height=4)
      ggsave(reactomeA$gridPlot, filename=file.path(outputDir3, "GSEA_Positive_reactome.png"),
             units = "in", width=6.5, height=4)
      ggsave(gobpA$gridPlot, filename=file.path(outputDir3, "GSEA_Positive_gobp.png"),
             units = "in", width=6.5, height=4)
      ggsave(complexA$gridPlot, filename=file.path(outputDir3, "GSEA_Positive_complex.png"),
             units = "in", width=6.5, height=4)
    }
    if(!is.null(gsea.neg) && nrow(gsea.neg)>0){
      ggsave(keggB$gridPlot, filename=file.path(outputDir3, "GSEA_Negative_kegg.png"),
             units = "in", width=6.5, height=4)
      ggsave(reactomeB$gridPlot, filename=file.path(outputDir3, "GSEA_Negative_reactome.png"),
             units = "in", width=6.5, height=4)
      ggsave(gobpB$gridPlot, filename=file.path(outputDir3, "GSEA_Negative_gobp.png"),
             units = "in", width=6.5, height=4)
      ggsave(complexB$gridPlot, filename=file.path(outputDir3, "GSEA_Negative_complex.png"),
             units = "in", width=6.5, height=4)
    }
  }
  
  ## Nine-squares ##
  {
    p1 = ScatterView(dd, x = "Control", y = "Treatment", label = "Symbol",
                     groups = c("midleft", "topcenter", "midright", "bottomcenter"),
                     groupnames = c("Group1", "Group2", "Group3", "Group4"),
                     top = top, toplabels = toplabels, omit = omit, display_cut = TRUE,
                     x_cut = x_cut, y_cut = y_cut, intercept = intercept,
                     xlab = lab_ctrl, ylab = lab_treat)
    ggsave(paste0(outputDir2, "/", "SquareView_", FDR, filter, ".png"), p1, 
           width = figwidth*4/3, height = figwidth, units = units, dpi = dpi)
    write.table(p1$data, paste0(outputDir2, "/", proj, "_squareview_data_", FDR, filter, ".txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    gridExtra::grid.arrange(p1, ncol = 1)
  }
  
  ## Nine-Square grouped gene enrichment ##
  {
    E1 = EnrichSquare(p1$data, id = "GeneID", keytype = "entrez",
                      x = "Control", y = "Treatment", top = top,
                      enrich_method = enrich_method, limit = limit,
                      filename=proj, out.dir=outputDir3)
    # EnrichView
    gridExtra::grid.arrange(E1$kegg1$gridPlot, E1$reactome1$gridPlot, E1$gobp1$gridPlot, E1$complex1$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg2$gridPlot, E1$reactome2$gridPlot, E1$gobp2$gridPlot, E1$complex2$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg3$gridPlot, E1$reactome3$gridPlot, E1$gobp3$gridPlot, E1$complex3$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg4$gridPlot, E1$reactome4$gridPlot, E1$gobp4$gridPlot, E1$complex4$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg12$gridPlot, E1$reactome12$gridPlot, E1$gobp12$gridPlot, E1$complex12$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg13$gridPlot, E1$reactome13$gridPlot, E1$gobp13$gridPlot, E1$complex13$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg24$gridPlot, E1$reactome24$gridPlot, E1$gobp24$gridPlot, E1$complex24$gridPlot, ncol = 2)
    gridExtra::grid.arrange(E1$kegg34$gridPlot, E1$reactome34$gridPlot, E1$gobp34$gridPlot, E1$complex34$gridPlot, ncol = 2)
    
    # PathwayView
    if(!is.null(E1$kegg1$enrichRes) && nrow(E1$kegg1$enrichRes)>0)
      arrangePathview(dd, gsub("KEGG_", "", E1$kegg1$enrichRes$ID), ncol = 2,
                      top = pathview.top, title = "Group 1", organism=organism,
                      output=outputDir4)
    if(!is.null(E1$kegg2$enrichRes) && nrow(E1$kegg2$enrichRes)>0)
      arrangePathview(dd, gsub("KEGG_", "", E1$kegg2$enrichRes$ID), ncol = 2,
                      top = pathview.top, title = "Group 2",
                      organism=organism,output=outputDir4)
    if(!is.null(E1$kegg3$enrichRes) && nrow(E1$kegg3$enrichRes)>0)
      arrangePathview(dd, gsub("KEGG_", "", E1$kegg3$enrichRes$ID), ncol = 2,
                      top = pathview.top, title = "Group 3",
                      organism=organism, output=outputDir4)
    if(!is.null(E1$kegg4$enrichRes) && nrow(E1$kegg4$enrichRes)>0)
      arrangePathview(dd, gsub("KEGG_", "", E1$kegg4$enrichRes$ID), ncol = 2,
                      title = "Group 4", organism = organism,
                      top = pathview.top, output=outputDir4)
    if(!is.null(E1$kegg12$enrichRes) && nrow(E1$kegg12$enrichRes)>0)
      arrangePathview(dd, gsub("KEGG_", "", E1$kegg12$enrichRes$ID), ncol = 2,
                      title = "Group 1 & Group 2", organism=organism,
                      top = pathview.top, output=outputDir4)
    if(!is.null(E1$kegg13$enrichRes) && nrow(E1$kegg13$enrichRes)>0)
      arrangePathview(dd, gsub("KEGG_", "", E1$kegg13$enrichRes$ID), ncol = 2,
                      title = "Group 1 & Group 3", organism=organism,
                      top = pathview.top, output=outputDir4)
    if(!is.null(E1$kegg24$enrichRes) && nrow(E1$kegg24$enrichRes)>0)
      arrangePathview(dd, gsub("KEGG_", "", E1$kegg24$enrichRes$ID), ncol = 2,
                      title = "Group 2 & Group 4", organism=organism,
                      top = pathview.top, output=outputDir4)
    if(!is.null(E1$kegg34$enrichRes) && nrow(E1$kegg34$enrichRes)>0)
      arrangePathview(dd, gsub("KEGG_", "", E1$kegg34$enrichRes$ID), ncol = 2,
                      title = "Group 3 & Group 4", organism=organism,
                      top = pathview.top, output=outputDir4)
  }
  dev.off()
}

#' Blank figure
#'
#' @docType methods
#' @name noEnrichPlot
#' @rdname noEnrichPlot
#' @param main The title of figure.
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#' @author Wubing Zhang
noEnrichPlot = function(main = "No enriched terms"){
  p1 = ggplot()
  p1 = p1 + geom_text(aes(x=0, y=0, label="No enriched terms"), size=6)
  p1 = p1 + labs(title=main)
  p1 = p1 + theme(plot.title = element_text(size=12))
  p1 = p1 + theme_bw(base_size = 14)
  p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
  p1
}