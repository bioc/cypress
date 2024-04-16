###### note ##########
# TOASDT DE methods need: disease group, est_CT_prop ,RNAseq_final_count
# TOAST::findRefinx TRIM GENE

# OUTPUT: sim_res_TOAST_strata <- TOAST::csTest
# ct_TDR_bio<-ct_tdr_bio, TDR_bio <- tdr_bio ,
# ct_PWR_bio<-ct_pwr,  PWR_bio<-ct_pwr_m,  PWR_strata_bio<-POWER_strata_bio,
# PWR_strata_ct_bio<-POWER_strata_ct
# ct_FDC_bio<-FDC_bio_ct FDC_bio<-FDC_bio_m

# DESeq2: > colnames(gene.CT.FDR.bl)
# [1] "gene_id"   "DE.CT"     "strata"    "Celltype1" "Celltype2" "Celltype3" "Celltype4"
# [8] "Celltype5" "Celltype6"
# cedar: > colnames(gene.CT.FDR)
# [1] "gene_id"   "Celltype1" "Celltype2" "Celltype3" "Celltype4" "Celltype5" "Celltype6"
# [8] "DE.CT"     "strata"
# TOAST: > colnames(gene_CT_FDR)
# [1] "gene_id"   "Celltype1" "Celltype2" "Celltype3" "Celltype4" "Celltype5" "Celltype6"
# [8] "DE.CT"     "strata"

###### transfer the data ##########

# se_object <- SummarizedExperiment(
#   assays = SimpleList(counts = as(RNAseq_final_count, "DataFrame")),
#   colData =  data.frame(cbind(disease=rep(c(1,2),100),
#                               sample_id=rownames(sample_CT_prop), sample_CT_prop) )
# )
#
# str(se_object)





###### Function for DESeq2 and cedar and TOAST##########

#   DESeq2
DESeq2_imple_part <- function(nsample_each_group,est_CT_prop,RNAseq_final_count,
                             ncell_type,gene_CT_DE_connect){

  reorder.ind <- c(paste0("control",seq_len(nsample_each_group) ),paste0("case",seq_len(nsample_each_group) ))
  est_CT_prop <- est_CT_prop[reorder.ind,]
  # should order the matrix for DES2q2
  RNAseq_final_count <- RNAseq_final_count[,reorder.ind]
  # create coldata for model fitting
  coldata <- cbind(disease = factor(as.numeric(gl(2,nsample_each_group))-1), est_CT_prop)
  # Organize the formular used for DEseq2
  cell_types <- paste0("Celltype", seq_len(ncell_type))
  main_effects_str <- paste("~0 +", paste(cell_types, collapse = " + "))
  interaction_effects_str <- paste(paste("disease:", cell_types, sep = ""), collapse = " + ")
  formula_str <- paste(main_effects_str, interaction_effects_str, sep = " + ")
  formula <- as.formula(formula_str)

  # Fitting DESeq2
  dds.bl <- DESeq2::DESeqDataSetFromMatrix(countData = RNAseq_final_count,
                                           colData = coldata,
                                           design= formula)
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
    dplyr::select("gene_id", "DE_CT","LFC", "strata")

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


  gene.CT.FDR <- merge(x = gene.CT.FDR, y = gene_CT_DE_connect, by = "gene_id",all.x =TRUE)
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
#
#
#
# #
# i=8
# n_sim = 30; n_gene = 30000; DE_pct = 0.05;
# ncell_type = 6; ss_group_set = c(10,20,50);
# lfc_set = c(0, 0.5, 1, 1.5);
# sim_param = ibd_prop_param ; DEmethod="TOAST";
# lfc_target = 0.5; fdr_thred = 0.1; BPPARAM=BiocParallel::bpparam()
# scenarios <- expand.grid(ss_group = ss_group_set,
#                          lfc = lfc_set)
# lfc_mean <- scenarios[i,'lfc']
# nsample_each_group <- scenarios[i,'ss_group']
# lfc_sd <- 0.5
# top_g <- seq(50, 600, by = 100)
# health_lmean_m <- sim_param@health_lmean_m
# health_lmean_d <- sim_param@health_lmean_d
# lod_m <- sim_param@lod_m
# lod_d <- sim_param@lod_d
# health_alpha <- sim_param@health_alpha
# case_alpha <- sim_param@case_alpha
#
# exp_strata <- c(-Inf, 10, 20, 40, 80,
#                 160, 320, 640, 1280, Inf)
#
# n_scenarios <- dim(scenarios)[1]
# n_strata <- length(exp_strata)-1
#
# BiocParallel::bplapply(seq_len(n_sim),csRNA_seq_sim,
#                        n_gene,DE_pct,ncell_type,lfc_mean,lfc_sd,
#                        health_lmean_m,health_lmean_d,lod_m,lod_d,
#                        nsample_each_group,lfc_target,fdr_thred,top_g,
#                        health_alpha,case_alpha,n_strata,DEmethod,
#                        BPPARAM=BPPARAM)
#
# csRNA_seq_sim(i=8,n_gene=30000,DE_pct=0.05,ncell_type=6,
#    lfc_mean=lfc_mean,lfc_sd=lfc_sd,health_lmean_m=health_lmean_m,
#     health_lmean_d=health_lmean_d,lod_m=lod_m,lod_d=lod_d,
#               nsample_each_group=nsample_each_group,lfc_target=lfc_target,
#    fdr_thred=fdr_thred,top_g=top_g,
#               health_alpha=health_alpha,case_alpha=case_alpha,n_strata=n_strata,DEmethod=DEmethod
# )
