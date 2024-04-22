
#### Function for DESeq2 and cedar and TOAST ####
##   DESeq2
DESeq2_imple_part <- function(nsample_each_group, est_CT_prop, RNAseq_final_count,
                             ncell_type, gene_CT_DE_connect){

  reorder.ind <- c(paste0("control",seq_len(nsample_each_group) ),
                   paste0("case",seq_len(nsample_each_group) ))
  est_CT_prop <- est_CT_prop[reorder.ind, ]
  RNAseq_final_count <- RNAseq_final_count[,reorder.ind]
  coldata <- cbind(disease = factor(as.numeric(gl(2,nsample_each_group))-1), est_CT_prop)
  cell_types <- paste0("Celltype", seq_len(ncell_type))
  main_effects_str <- paste("~0 +", paste(cell_types, collapse = " + "))
  interaction_effects_str <- paste(paste("disease:", cell_types, sep = ""), collapse = " + ")
  formula_str <- paste(main_effects_str, interaction_effects_str, sep = " + ")
  formula <- as.formula(formula_str)

  dds.bl <- DESeq2::DESeqDataSetFromMatrix(countData = RNAseq_final_count,
                       colData = coldata, design= formula)
  dds.bl <- DESeq2::DESeq(dds.bl)
  celltype_results_list <- list()

  # Loop through each cell type to get DESeq2 results
  for (cell_type in cell_types) {
    result_name <- paste(cell_type, "disease", sep = ".")
    res <- as.data.frame(DESeq2::results(dds.bl, name = result_name))
    res$padj <- ifelse(is.na(res$padj), 1, res$padj) # Adjust p-values
    res$gene_id <- rownames(res)
    # Store the results with the cell type name as the adjusted p-values column
    celltype_results_list[[cell_type]] <- res %>%
      dplyr::select("gene_id", "padj") %>%
      dplyr::rename(!!cell_type := "padj")
  }

  gene.CT.FDR.bl <- gene_CT_DE_connect %>%
    dplyr::select("gene_id", "DE_CT", "LFC", "strata")

  # Use Reduce to merge all cell type results into gene.CT.FDR.bl
  gene.CT.FDR.bl <- Reduce(function(x, y) merge(x = x, y = y, by = "gene_id", all.x = TRUE),
                           c(celltype_results_list,list(gene.CT.FDR.bl)) )
  return(gene.CT.FDR.bl)
}

##   cedar
cedar_imple_part <- function(nsample_each_group,est_CT_prop,RNAseq_final_count,
                             ncell_type,gene_CT_DE_connect){
  sim.design <- as.data.frame(factor(rep(c(1,2),nsample_each_group)))
  colnames(sim.design) <- "disease"
  res_cedar <- TOAST::cedar(Y_raw = as.matrix(RNAseq_final_count), prop = est_CT_prop, design.1 = sim.design,
                            factor.to.test = 'disease',cutoff.tree = c('pval',0.01),
                            cutoff.prior.prob = c('pval',0.01))
  # Pull out model outputs
  gene.CT.FDR <- as.data.frame(1- res_cedar$tree_res$full$pp)
  gene.CT.FDR$gene_id <- (rownames(gene.CT.FDR))
  gene.CT.FDR <- gene.CT.FDR[,c(seq_len(ncell_type)+1,seq_len(ncell_type) )]

  gene_CT_DE_connect <- gene_CT_DE_connect %>%
    dplyr::select("gene_id", "DE_CT", "LFC","strata")


  gene.CT.FDR <- merge(x = gene.CT.FDR, y = gene_CT_DE_connect, by = "gene_id", all.x =TRUE)
  return(gene.CT.FDR)
}

##   TOAST
TOAST_imple_part <- function(nsample_each_group,est_CT_prop,RNAseq_final_count,
                             ncell_type,gene_CT_DE_connect){

sim_res_TOAST_strata <- run_TOAST(nsample_each_group, est_CT_prop,
                                  RNAseq_final_count)
TOAST_out <- summary_TOAST(sim_res_TOAST_strata, ncell_type)
gene_CT_FDR <- merge(TOAST_out, gene_CT_DE_connect, by = "gene_id")
return(gene_CT_FDR)
}


