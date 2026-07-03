#' @title title
#'
#' @description
#' Subset TOBIAS BINDetect results to regions overlapping a BED file.
#' Re-summarises a single BINDetect condition-comparison, restricting the
#' analysis to the transcription-factor binding sites (TFBS) that overlap a
#' user-supplied set of regions. It reads the per-TF \code{<TF>_overview.txt}
#' files produced by \code{TOBIAS BINDetect}, keeps the sites overlapping
#' \code{bed}, and recomputes per-TF summary statistics (mean/median footprint
#' scores, bound counts, and a differential-binding effect) together with a
#' significance test against a background distribution of per-site log2 fold
#' changes.
#'
#' @section Important statistical note:
#' In \code{TOBIAS BINDetect}, the \code{*_change} and \code{*_pvalue} columns of
#' \code{bindetect_results.txt} are \emph{not} self-contained per-TF quantities:
#' the between-condition score normalisation, the bound/unbound threshold and the
#' null model for the p-values are all estimated from the \emph{whole} peak set.
#' Restricting to a sub-region and rebuilding the background from that subset
#' therefore skews the distribution (see loosolab/TOBIAS issue #48). This
#' function instead keeps the per-site log2 fold changes exactly as BINDetect
#' computed them (they are already normalised between conditions) and, by
#' default, tests each TF's subset against a \emph{genome-wide} background
#' (\code{background = "all"}). If your regions are large (>500) and you want a
#' properly re-centred volcano, prefer re-running BINDetect with your BED as
#' \code{--peaks} (see \code{\link{rerun_bindetect_on_bed}}).
#'
#' @param bindetect_dir Path to the BINDetect output directory for one
#'   comparison, i.e. the folder that contains \code{bindetect_results.txt} and
#'   one sub-directory per TF (e.g.
#'   \code{05b_Differential_TF_binding_TOBIAS/C_BINDetect_merged_BAMs/groupA.vs.groupB}).
#' @param bed Path to a BED file (or a \code{data.frame}/\code{GRanges}) with the
#'   regions of interest to which the analysis should be restricted.
#' @param condition1,condition2 Character scalars with the two condition names,
#'   exactly as passed to \code{--cond_names} (used to locate the
#'   \code{<cond>_score}, \code{<cond>_bound} and
#'   \code{<cond1>_<cond2>_log2fc} columns). If \code{NULL} (default) they are
#'   inferred from the header of \code{bindetect_results.txt}.
#' @param bound_in_any Logical. If \code{TRUE}, keep only sites predicted bound
#'   in at least one condition before summarising. Note this biases the log2fc
#'   distribution (it removes unbound-unbound sites, which sit near zero) and
#'   therefore your top-quantile selection; exposed as an option rather than
#'   hard-coded. Default \code{FALSE}.
#' @param background Which sites form the null for the significance test:
#'   \code{"all"} (default) uses the pooled per-site log2fc of every TF
#'   genome-wide (mirrors the "true background = full peak set" logic); "outside"
#'   uses, per TF, that TF's sites \emph{not} overlapping \code{bed} (answers the
#'   different question "is this TF's differential binding specific to my
#'   regions?").
#' @param summary_stat Central tendency for the per-TF effect, one of
#'   \code{"mean"} or \code{"median"} of the subset per-site log2fc. Default
#'   \code{"mean"}.
#' @param min_sites Minimum number of subset sites required to test a TF; TFs
#'   below this get \code{NA} p-values. Default \code{10}.
#' @param top_quantile Tail fraction used to flag "top" TFs, mirroring the
#'   BINDetect volcano (effect in the lower/upper \code{top_quantile} tails
#'   and/or \code{-log10(padj)} above the upper tail). Default \code{0.05}.
#' @param padjust_method Passed to \code{\link[stats]{p.adjust}}. Default
#'   \code{"BH"}.
#'
#' @return A \code{data.frame}, one row per TF, with columns: \code{TF},
#'   \code{n_sites_total}, \code{n_sites_subset}, \code{<cond>_mean_score},
#'   \code{<cond>_bound}, \code{subset_change} (the chosen central tendency of
#'   the subset log2fc; negative = higher in \code{condition2}), \code{pvalue},
#'   \code{padj} and \code{top} (logical). Sorted by \code{padj} then
#'   \code{abs(subset_change)}.
#'
#' @examples
#' \dontrun{
#' res = subset_bindetect_by_bed(
#'   bindetect_dir = "05b_Differential_TF_binding_TOBIAS/C_BINDetect_merged_BAMs/groupA.vs.groupB",
#'   bed = "my_regions.bed",
#'   condition1 = "groupA",
#'   condition2 = "groupB")
#' head(res[res$top, ])
#' }
#'
#' @importFrom data.table fread rbindlist
#' @import GenomicRanges
#' @import IRanges
#' @import stats
#'
#' @export


subset_BINDetect = function(bindetect_dir,
                            bed,
                            condition1 = NULL,
                            condition2 = NULL,
                            bound_in_any = FALSE,
                            background = c("all", "outside"),
                            summary_stat = c("mean", "median"),
                            min_sites = 10L,
                            top_quantile = 0.05,
                            padjust_method = "BH") {

  background = match.arg(background)
  summary_stat = match.arg(summary_stat)
  stat_fun = if (summary_stat == "mean") {function(x) mean(x, na.rm = TRUE)} else {function(x) stats::median(x, na.rm = TRUE)}


  ## resolve condition names from bindetect_results.txt header if needed ----
  results_txt = file.path(bindetect_dir, "bindetect_results.txt")
  if (is.null(condition1) || is.null(condition2)) {

    if (!file.exists(results_txt))
      stop("Cannot infer conditions: ", results_txt, " not found. ",
           "Pass condition1/condition2 explicitly.", call. = FALSE)
    hdr = strsplit(readLines(results_txt, n = 1L), "\t")[[1]]
    bound = sub("_bound$", "", grep("_bound$", hdr, value = TRUE))

    if (length(bound) < 2L)
      stop("Could not detect two conditions from the header of ", results_txt,
           call. = FALSE)
    condition1 = bound[1]; condition2 = bound[2]
  }

  log2fc_col = paste0(condition1, "_", condition2, "_log2fc")
  s1 = paste0(condition1, "_score"); s2 = paste0(condition2, "_score")
  b1 = paste0(condition1, "_bound"); b2 = paste0(condition2, "_bound")

  ## locate the per-TF overview files ---------------------------------------
  overview_files = list.files(bindetect_dir,
                              pattern = "_overview\\.txt$",
                              full.names = TRUE, recursive = TRUE)

  if (length(overview_files) == 0L)
    stop("No '*_overview.txt' files found under ", bindetect_dir,
         ". Was BINDetect run without --skip-overview?", call. = FALSE)

  ## read all sites into one table ------------------------------------------
  need = c("TFBS_chr", "TFBS_start", "TFBS_end", s1, s2, b1, b2, log2fc_col)

  read_one = function(f) {
    dt = data.table::fread(f, sep = "\t", header = TRUE, showProgress = FALSE)
    miss = setdiff(need, names(dt))
    if (length(miss))
      stop("File ", basename(f), " is missing column(s): ",
           paste(miss, collapse = ", "), call. = FALSE)
    dt = dt[, need, with = FALSE]
    dt$TF = sub("_overview\\.txt$", "", basename(f))
    dt
  }

  sites = data.table::rbindlist(lapply(overview_files, read_one))

  ## overlap sites with the regions of interest -----------------------------
  ## TOBIAS coordinates are BED-style (0-based, half-open); convert start to
  ## 1-based for GRanges. The ROI BED is treated identically, so the overlap is
  ## internally consistent.
  roi = .as_granges(bed)

  tf_gr = GenomicRanges::GRanges(seqnames = sites$TFBS_chr,
                                 ranges = IRanges::IRanges(start = sites$TFBS_start + 1L, end = sites$TFBS_end))

  sites$in_subset = IRanges::overlapsAny(tf_gr, roi)

  ## background vector of per-site log2fc ------------------------------------
  bg_all = sites[[log2fc_col]]
  bg_all = bg_all[is.finite(bg_all)]

  ## per-TF summarisation ----------------------------------------------------
  split_idx = split(seq_len(nrow(sites)), sites$TF)

  rows = lapply(names(split_idx),
                function(tf) {idx = split_idx[[tf]]
                sub = idx[sites$in_subset[idx]]
                if (bound_in_any)
                  sub = sub[(sites[[b1]][sub] > 0) | (sites[[b2]][sub] > 0)]

                l2fc_sub = sites[[log2fc_col]][sub]; l2fc_sub = l2fc_sub[is.finite(l2fc_sub)]

                ## choose the null
                bg = if (background == "all") bg_all else {
                  out = idx[!sites$in_subset[idx]]
                  v = sites[[log2fc_col]][out]; v[is.finite(v)]
                }

                pval = NA_real_

                if (length(l2fc_sub) >= min_sites && length(bg) >= min_sites)
                  pval = tryCatch(
                    stats::wilcox.test(l2fc_sub, bg)$p.value,
                    error = function(e) NA_real_)

                data.frame(
                  TF = tf,
                  n_sites_total = length(idx),
                  n_sites_subset = length(sub),
                  c1_mean_score = mean(sites[[s1]][sub], na.rm = TRUE),
                  c2_mean_score = mean(sites[[s2]][sub], na.rm = TRUE),
                  c1_bound = sum(sites[[b1]][sub] > 0, na.rm = TRUE),
                  c2_bound = sum(sites[[b2]][sub] > 0, na.rm = TRUE),
                  subset_change = if (length(l2fc_sub)) stat_fun(l2fc_sub) else NA_real_,
                  pvalue = pval,
                  stringsAsFactors = FALSE)
                })

  out = do.call(rbind, rows)

  ## rename the generic score/bound columns to carry the condition names
  names(out)[names(out) == "c1_mean_score"] = paste0(condition1, "_mean_score")
  names(out)[names(out) == "c2_mean_score"] = paste0(condition2, "_mean_score")
  names(out)[names(out) == "c1_bound"] = paste0(condition1, "_bound")
  names(out)[names(out) == "c2_bound"] = paste0(condition2, "_bound")

  ## multiple testing + top-quantile flag (mirrors BINDetect volcano) --------
  out$padj = stats::p.adjust(out$pvalue, method = padjust_method)

  ch = out$subset_change
  lo = stats::quantile(ch, top_quantile,     na.rm = TRUE)
  hi = stats::quantile(ch, 1 - top_quantile, na.rm = TRUE)
  nlp = -log10(out$padj)
  nlp_hi = stats::quantile(nlp[is.finite(nlp)], 1 - top_quantile, na.rm = TRUE)
  out$top = (!is.na(ch) & (ch <= lo | ch >= hi)) |
    (is.finite(nlp) & nlp >= nlp_hi)

  return(out[order(out$padj, -abs(out$subset_change)), ])
}


#' @keywords internal
#' @importFrom methods is
#' @importFrom data.table fread as.data.table
#' @import IRanges
#' @import package GenomicRanges
#' @noRd


.as_granges = function(bed) {
  if (methods::is(bed, "GRanges")) return(bed)
  if (is.character(bed) && length(bed) == 1L) {
    df = data.table::fread(bed, header = FALSE, sep = "\t", showProgress = FALSE)
  } else if (is.data.frame(bed)) {
    df = data.table::as.data.table(bed)
  } else {
    stop("'bed' must be a file path, data.frame, or GRanges.", call. = FALSE)
  }
  GenomicRanges::GRanges(
    seqnames = df[[1]],
    ranges = IRanges::IRanges(start = df[[2]] + 1L, end = df[[3]]))
}
