#' @title filter.hocomoco.jaspar
#'
#' @description
#' Filter HOCOMOCO (v12+) motifs from a JASPAR-format file
#' Reads a HOCOMOCO JASPAR-format motif file and returns the subset matching the
#' requested transcription factors and/or motif-quality criteria. Motifs can be
#' selected by gene symbol (bridged through the HOCOMOCO annotation, because a
#' motif-ID prefix is a UniProt mnemonic and \emph{not} always the HGNC symbol --
#' e.g. \code{ANDR} -> \code{AR}, \code{AP2A} -> \code{TFAP2A}, \code{SRBP2} ->
#' \code{SREBF2}) and/or by exact HOCOMOCO motif IDs, and can be further
#' restricted by model subtype, quality grade, and experiment type. Kept motifs
#' are returned as a list of \code{universalmotif} objects and, if \code{out_file}
#' is supplied, written out with \code{universalmotif::write_jaspar()}.
#'
#' @section HOCOMOCO v14 motif-ID anatomy:
#' A v14 ID such as \code{ANDR.H14CORE.0.P.B} splits (from the right) into
#' \emph{quality} (\code{B}; grades A-D, A best), \emph{experiment} (\code{P} =
#' ChIP-seq, \code{S} = HT-SELEX, \code{M} = Methyl-HT-SELEX, or combinations
#' such as \code{PSM}), \emph{subtype} (\code{0}; 0 is the primary model) and
#' \emph{collection} (\code{H14CORE}), with the TF mnemonic in front. The
#' familiar v11 shorthand "0.A" therefore corresponds here to
#' \code{subtypes = 0, qualities = "A"} -- the experiment letters sit in between,
#' so "0.A" is not a literal substring of the v14 ID.
#'
#' @param jaspar_file Path to the HOCOMOCO JASPAR-format file (e.g. \code{"H14CORE_jaspar_format.txt"}).
#' @param genes Character vector of gene symbols to keep (case-insensitive). Requires \code{annotation_file}. A motif is kept if \emph{any} of its associated genes matches.
#' @param motif_ids Character vector of exact HOCOMOCO motif IDs to keep.
#' @param subtypes Character/numeric vector of model subtypes to keep (e.g. \code{0}). \code{NULL} (default) keeps all.
#' @param qualities Character vector of quality grades to keep (e.g. \code{c("A", "B")}). \code{NULL} (default) keeps all.
#' @param experiments Character vector of experiment codes to keep, matched exactly on the field (e.g. \code{"P"}, \code{"PSM"}). \code{NULL} keeps all.
#' @param annotation_file Path to the HOCOMOCO annotation JSONL (e.g. \code{"H14CORE_annotation.jsonl"}). Required when \code{genes} is given or \code{names_as = "gene"}.
#' @param names_as Name the returned/exported motifs by \code{"motif"} ID (default) or by \code{"gene"} symbol. The motif ID is always retained in the \code{altname} slot.
#' @param gene_sep Separator used to join symbols when a motif maps to several genes and \code{names_as = "gene"} (default \code{"::"}).
#' @param out_file Optional output path; if supplied the filtered motifs are written there in JASPAR format.
#' @param overwrite Passed to \code{universalmotif::write_jaspar()}. Default: \code{TRUE}.
#'
#' @return A list of \code{universalmotif} objects (the kept motifs), returned invisibly when \code{out_file} is written and visibly otherwise.
#'
#' @examples
#' \dontrun{
#' jf <- "H14CORE_jaspar_format.txt"
#' af <- "H14CORE_annotation.jsonl"
#'
#' # All motifs for a set of genes, primary model at top quality only:
#' filter.hocomoco.jaspar(jf, genes = c("CTCF", "GATA1", "AR"),
#'                        subtypes = 0, qualities = c("A", "B"),
#'                        annotation_file = af)
#'
#' # Specific motif IDs, renamed to gene symbols, exported to disk:
#' filter.hocomoco.jaspar(jf, motif_ids = "CTCF.H14CORE.0.PSM.A",
#'                        annotation_file = af, names_as = "gene",
#'                        out_file = "selected.jaspar")
#' }
#'
#' @importFrom jsonlite fromJSON
#' @importFrom universalmotif read_jaspar write_jaspar
#' @importFrom methods is
#'
#' @export filter.hocomoco.jaspar



filter.hocomoco.jaspar =
  function(jaspar_file,
           genes = NULL,
           motif_ids = NULL,
           subtypes = NULL,
           qualities = NULL,
           experiments = NULL,
           annotation_file = NULL,
           names_as = c("motif", "gene"),
           gene_sep = "::",
           out_file = NULL,
           overwrite = TRUE) {

    names_as = match.arg(names_as)
    if (!requireNamespace("universalmotif", quietly = TRUE))
      stop("Package 'universalmotif' is required. Install it from Bioconductor.")
    if ((!is.null(genes) || names_as == "gene") && is.null(annotation_file))
      stop("annotation_file is required when 'genes' is supplied or names_as = 'gene' ",
           "(motif-ID prefixes are UniProt mnemonics, not reliable gene symbols).")
    if (!file.exists(jaspar_file)) stop("jaspar_file not found: ", jaspar_file)

    ## --- read motifs + authoritative IDs (aligned by file order) ---
    motifs = universalmotif::read_jaspar(jaspar_file)
    if (methods::is(motifs, "universalmotif")) motifs = list(motifs)
    ids = .hoco_jaspar_ids(jaspar_file)
    if (length(ids) != length(motifs)) {
      warning("Header count (", length(ids), ") != motif count (", length(motifs),
              "); falling back to motif 'name' slots for IDs.")
      ids = vapply(motifs, function(m) as.character(m["name"]), character(1))
    }

    ## --- gene mapping (only if needed) ---
    id2gene = if (!is.null(annotation_file)) .hoco_id2genes(annotation_file) else NULL
    genes_of = function(id) { g = if (is.null(id2gene)) NULL else id2gene[[id]]; if (is.null(g)) character(0) else g }
    gl     = lapply(ids, genes_of)
    parsed = lapply(ids, .hoco_parse_id)

    ## --- filters ---
    keep = rep(TRUE, length(ids))
    if (!is.null(subtypes))
      keep = keep & vapply(parsed, function(p) !is.na(p$subtype) &&
                             p$subtype %in% as.character(subtypes), logical(1))
    if (!is.null(qualities))
      keep = keep & vapply(parsed, function(p) !is.na(p$quality) &&
                             toupper(p$quality) %in% toupper(qualities), logical(1))
    if (!is.null(experiments))
      keep = keep & vapply(parsed, function(p) !is.na(p$experiment) &&
                             p$experiment %in% experiments, logical(1))
    if (!is.null(genes) || !is.null(motif_ids)) {              # union of identity filters
      idm = rep(FALSE, length(ids))
      if (!is.null(genes)) {
        gU = toupper(genes)
        idm = idm | vapply(gl, function(g) any(toupper(g) %in% gU), logical(1))
      }
      if (!is.null(motif_ids)) idm = idm | (ids %in% motif_ids)
      keep = keep & idm
    }

    idx  = which(keep)
    kept = motifs[idx]
    kid  = ids[idx]
    kgl  = gl[idx]

    if (length(kept) == 0L) {
      warning("No motifs matched the requested filters.")
      return(invisible(list()))
    }

    ## --- set output names (motif ID always kept in altname) ---
    n_nogene = 0L
    for (i in seq_along(kept)) {
      kept[[i]]["altname"] = kid[i]
      if (names_as == "gene") {
        if (length(kgl[[i]]))
          kept[[i]]["name"] = paste(sort(unique(kgl[[i]])), collapse = gene_sep)
        else { kept[[i]]["name"] = kid[i]; n_nogene = n_nogene + 1L }
      } else {
        kept[[i]]["name"] = kid[i]
      }
    }
    if (n_nogene > 0L)
      warning(n_nogene, " kept motif(s) had no gene in the annotation; kept the ID as name.")

    message(sprintf("Kept %d / %d motifs (%d distinct genes).",
                    length(kept), length(motifs),
                    length(unique(unlist(kgl)))))

    ## --- optional export ---
    if (!is.null(out_file)) {
      universalmotif::write_jaspar(kept, file = out_file, overwrite = overwrite)
      message("Wrote ", out_file)
      return(invisible(kept))
    }
    kept
  }

## ------------------------- internal helpers -------------------------

# Motif IDs from the ">" headers of a JASPAR file, in file order.
.hoco_jaspar_ids = function(jaspar_file) {
  lines = readLines(jaspar_file, warn = FALSE)
  hdr   = lines[grepl("^>", lines)]
  vapply(hdr, function(h) strsplit(trimws(sub("^>", "", h)), "\\s+")[[1]][1],
         character(1), USE.NAMES = FALSE)
}

# Split a HOCOMOCO v12+ ID from the right (robust to dots in the TF field).
.hoco_parse_id = function(id) {
  p = strsplit(id, ".", fixed = TRUE)[[1]]
  n = length(p)
  if (n < 5L) return(list(tf = NA, collection = NA, subtype = NA,
                          experiment = NA, quality = NA))
  list(tf = paste(p[seq_len(n - 4L)], collapse = "."), collection = p[n - 3L],
       subtype = p[n - 2L], experiment = p[n - 1L], quality = p[n])
}

# motif ID from an annotation record
.hoco_id_from_rec = function(rec) {
  for (k in c("name", "full_name", "motif", "id"))
    if (!is.null(rec[[k]]) && is.character(rec[[k]])) return(rec[[k]][1])
  NA_character_
}

# gene symbol(s) from an annotation record (prefer a HUMAN branch)
.hoco_genes_from_rec = function(rec) {
  flat = unlist(rec)
  if (length(flat) == 0) return(character(0))
  nm = names(flat); if (is.null(nm)) nm = rep("", length(flat))
  is_sym = grepl("gene_symbol|gene_name|hgnc|(^|\\.)symbol[0-9]*$", nm, ignore.case = TRUE)
  human  = is_sym & grepl("human", nm, ignore.case = TRUE)
  vals   = if (any(human)) flat[human] else flat[is_sym]
  vals   = trimws(as.character(vals))
  unique(vals[nzchar(vals)])
}

# motif_id -> character vector of gene symbols, from the annotation JSONL
.hoco_id2genes = function(annotation_file) {
  if (!file.exists(annotation_file)) stop("annotation_file not found: ", annotation_file)
  lines = readLines(annotation_file, warn = FALSE)
  lines = lines[nzchar(trimws(lines))]
  map = list()
  for (ln in lines) {
    rec = jsonlite::fromJSON(ln, simplifyVector = TRUE)
    id  = .hoco_id_from_rec(rec)
    if (!is.na(id)) map[[id]] = .hoco_genes_from_rec(rec)
  }
  map
}
