## ----setup, include = FALSE---------------------------------------------------
# knitr::opts_chunk$set(collapse = TRUE, comment = ">", dev = "svg",
#                       warning = F, message = F, fig.align = "center",
#                       rows.print=12)
#knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_

if (Sys.info()[["sysname"]] == "Darwin") {
  knitr::opts_chunk$set(collapse = TRUE, comment = ">",
                        warning = F, message = F, fig.align = "center",
                        rows.print = 12, dev = "png", dpi = 96)
} else {
  knitr::opts_chunk$set(collapse = TRUE, comment = ">",
                        warning = F, message = F, fig.align = "center",
                        rows.print = 12, dev = "png", dpi = 96)
}

#options(tibble.print_min = 4L, tibble.print_max = 4L)


# Load libraries required
require(snakeATAC.utility)

## ----citation, message=FALSE, warning=FALSE-----------------------------------
citation("snakeATAC.utility")

## ----eval=FALSE---------------------------------------------------------------
# res_subset <-
#   subset_BINDetect(bindetect_dir = bindetect_dir, # path
#                    bed = roi_bed, # GRanges, file_path or data.frame/data.table
#                    condition1 = "groupA",
#                    condition2 = "groupB",
#                    background = "all", # genome-wide null (default)
#                    summary_stat = "mean",
#                    bound_in_any = FALSE,
#                    min_sites = 10,
#                    top_quantile = 0.05)
# 
# head(res_subset) # one row per TF, sorted by padj then |subset_change|
# sum(res_subset$top) # TFs in the volcano tails

## ----background-outside, eval=FALSE-------------------------------------------
# res_specific <- subset_bindetect_by_bed(
#   bindetect_dir = bindetect_dir,
#   bed = roi_bed,
#   condition1 = "groupA",
#   condition2 = "groupB",
#   background = "outside")

## ----bound-filter, eval=FALSE-------------------------------------------------
# res_bound <-
#   subset_BINDetect(bindetect_dir = bindetect_dir,
#                    bed = roi_bed,
#                    condition1 = "groupA",
#                    condition2 = "groupB",
#                    bound_in_any = TRUE)

## ----top, eval=FALSE----------------------------------------------------------
# top_TFs <- res_subset[res_subset$top, c("TF", "subset_change", "padj")]
# top_TFs

## -----------------------------------------------------------------------------
jaspar_file <- "H14CORE_jaspar_format.txt"
annotation_file <- "H14CORE_annotation.jsonl"

## ----eval=F-------------------------------------------------------------------
# # By gene symbol (case-insensitive; a motif is kept if any of its genes matches):
# motifs <- filter.hocomoco.jaspar(jaspar_file = jaspar_file,
#                                  genes = c("CTCF", "GATA1", "AR"),
#                                  annotation_file = annotation_file)

## ----eval=F-------------------------------------------------------------------
# best <- filter.hocomoco.jaspar(jaspar_file = jaspar_file,
#                                genes = c("CTCF", "GATA1", "AR"),
#                                subtypes = 0,
#                                qualities = c("A", "B"),
#                                annotation_file = annotation_file)

## ----eval=F-------------------------------------------------------------------
# picked <- filter.hocomoco.jaspar(jaspar_file = jaspar_file,
#                                  motif_ids = c("CTCF.H14CORE.0.PSM.A", "GATA1.H14CORE.0.PSM.A"))

## ----eval=F-------------------------------------------------------------------
# filter.hocomoco.jaspar(jaspar_file = jaspar_file,
#                        genes = "CTCF",
#                        annotation_file = annotation_file,
#                        names_as = "gene",
#                        out_file = "ctcf_motifs.jaspar")

## ----eval=F-------------------------------------------------------------------
# corr <- file.path("06_Overall_quality_and_info", "Sample_comparisons",
#                   "Sample_correlation",
#                   "Correlation_matrix_on_BigWigs_wholeGenome_spearmanMethod.txt")
# 
# correlation.heatmap(corr)

## ----eval=F-------------------------------------------------------------------
# # Tighten the colour range, move the dendrogram, use Ward clustering, hide labels:
# correlation.heatmap(correlation.matrix = corr,
#                     correlation.scale.limits = c(0.8, 1),
#                     dendrogram.position = "top",
#                     clustering.method = "ward.D2",
#                     display.values = FALSE)
# 
# # Drop samples and also return the clustering objects:
# out <- correlation.heatmap(correlation.matrix = corr,
#                            excluded.samples = c("sampleX", "sampleY"),
#                            return.cluster.object = TRUE)
# out$heatmap
# out$cluster

## ----eval=F-------------------------------------------------------------------
# npz <- file.path("06_Overall_quality_and_info", "Sample_comparisons",
#                  "Sample_correlation", "multiBigWigSummary_matrix_allSamples.npz")
# 
# mat <- read.multibigwigsummary(npz)   # columns = samples, rows = genomic bins
# dim(mat)
# colnames(mat)

## ----eval=F-------------------------------------------------------------------
# # Default PCA on PC1 / PC2:
# plot.PCA(mat)
# 
# # Drop samples, look at PC2 vs PC3, and set custom (named) colours:
# plot.PCA(matrix = mat,
#          excluded.samples = "sampleX",
#          PCs = c(2, 3),
#          color.list = c(WT_rep1 = "#1b9e77", WT_rep2 = "#1b9e77",
#                         KO_rep1 = "#d95f02", KO_rep2 = "#d95f02"),
#          point.size = 4)

## ----eval=F-------------------------------------------------------------------
# bench <- benchmark.summary("path/to/snakeATAC_output/benchmarks")
# bench

