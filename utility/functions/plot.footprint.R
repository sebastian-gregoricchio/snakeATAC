plot.footprint =
  function(matrix.gz,
           factor = NULL,
           colors = NULL,
           normalize.zero.one = FALSE,
           white.list = NULL) {

    # Load libraries
    require(dplyr)
    require(ggplot2)


    # read matrix
    matrix = as.data.frame(data.table::fread(matrix.gz, skip = 1))


    # Import metadata
    metadata = data.table::fread(matrix.gz, nrows = 1, stringsAsFactors = F, sep = "@", h = F)$V2
    metadata = gsub(x = metadata, pattern = "[{]|[}]", replacement = "")
    metadata = gsub(x = metadata, pattern = "\\s", replacement = "_")
    metadata = gsub(x = metadata, pattern = "[\"],\"", replacement = ",")
    metadata = gsub(x = metadata, pattern = "\":\"", replacement = ":")
    metadata = gsub(x = metadata, pattern = "\":", replacement = ":")
    metadata = gsub(x = metadata, pattern = ",\"", replacement = "@")
    metadata = gsub(x = metadata, pattern = "\"", replacement = "")
    metadata = gsub(x = metadata, pattern = "[,]missing", replacement = "@missing")
    metadata = gsub(x = metadata, pattern = "[,]sort", replacement = "@sort")
    metadata = gsub(x = metadata, pattern = "[,]unscaled", replacement = "@unscaled")

    metadata = unlist(lapply(strsplit(x = metadata, split = "@")[[1]], function(x)(strsplit(x, ":"))))

    metadata = data.frame("parameters" = metadata[seq(1, length(metadata),2)],
                          "values" = metadata[seq(2, length(metadata),2)],
                          stringsAsFactors = F)

    # Collect info
    samples = unlist(strsplit(gsub("[[]|[]]","",metadata[21,2]), ","))
    upstream = as.numeric(unique(unlist(strsplit(gsub("[[]|[]]","",metadata[1,2]), ","))))
    downstream = as.numeric(unique(unlist(strsplit(gsub("[[]|[]]","",metadata[2,2]), ","))))
    boundaries = as.numeric(unique(unlist(strsplit(gsub("[[]|[]]","",metadata[22,2]), ","))))

    range = (-upstream):downstream
    range = range[range != 0]


    # Re-shape table
    matrix.reshaped = data.frame()

    for (i in 1:length(samples)) {
      start.col = 7 + boundaries[i]
      end.col = 7 + boundaries[i+1] - 1

      mat = matrix[,c(1:6,start.col:end.col)]
      colnames(mat)[7:ncol(mat)] = range

      mat$sample = samples[i]
      matrix.reshaped = rbind(matrix.reshaped, mat)
    }

    matrix.reshaped =
      reshape2::melt(data = matrix.reshaped,
                     id.vars = c("V1", "V2", "V3", "V4", "V5", "V6", "sample"),
                     variable.name = "position",
                     value.name = "score") %>%
      dplyr::mutate(position = as.character(as.character(position)))


    ### Define the TF name
    if (is.null(factor)) {
      factor = gsub("_[A-Z]*$","", gsub("_[^_]*$","",matrix.reshaped$V4[1]))
    }

    ### Define colors
    if (is.null(colors) | (length(colors) < length(samples))) {
      colors = rainbow(length(samples))
    }

    ### Filter for the white.list (if required)
    if (!is.null(white.list)) {
      if ("data.frame" %in% class(white.list)) {
        white.ls = white.list[,1:3]
        colnames(white.ls) = c("seqnames", "start", "end")
        white.ls$seqnames = gsub("chr", "", white.ls$seqnames, ignore.case = T)
        white.ls = GenomicRanges::makeGRangesFromDataFrame(white.ls)
      } else if ("character" %in% class(white.list)) {
        white.ls = data.table::fread(white.list, data.table = F)[,1:3]
        colnames(white.ls) = c("seqnames", "start", "end")
        white.ls$seqnames = gsub("chr", "", white.ls$seqnames, ignore.case = T)
        white.ls = GenomicRanges::makeGRangesFromDataFrame(white.ls)
      } else if ("GRanges" %in% class(white.list)) {
        white.ls = as.data.frame(white.list)
        white.ls$seqnames = gsub("chr", "", white.ls$seqnames, ignore.case = T)
        white.ls = GenomicRanges::makeGRangesFromDataFrame(white.ls)
      } else {return(warning("The 'white.list' must be a bed in Granges, data.frame format or the path to a bed file."))}
    }


    ## Extract genomic ranges from matrix
    if (!is.null(white.list)) {
      matrix.bed = unique(matrix[,1:3])
      colnames(matrix.bed) = c("seqnames", "start", "end")
      matrix.bed$seqnames = gsub("chr", "", matrix.bed$seqnames, ignore.case = T)
      matrix.bed = dplyr::mutate(unique(matrix.bed), name = paste(seqnames, start, end, sep="_"))
      matrix.bed_gr = GenomicRanges::makeGRangesFromDataFrame(matrix.bed, keep.extra.columns = T) }


    ## Filter matrix for white.list
    if (!is.null(white.list)) {
      keep.region = IRanges::countOverlaps(query = matrix.bed_gr, subject = white.ls) > 0
      matrix.bed.filtered = matrix.bed[keep.region,]

      matrix.reshaped =
        dplyr::inner_join(matrix.reshaped %>% dplyr::mutate(V1 = gsub("chr","",V1, ignore.case = T)),
                          matrix.bed.filtered[,1:3],
                          by = c("V1" = "seqnames", "V2" = "start", "V3" = "end"))
    }


    # normalize within 0-1 if required
    if (normalize.zero.one == TRUE) {
      matrix.reshaped$score = (matrix.reshaped$score - min(matrix.reshaped$score, na.rm = T))
      matrix.reshaped$score = matrix.reshaped$score / max(matrix.reshaped$score, na.rm = T)
    }


    # generate stats
    mat.plot =
      matrix.reshaped %>%
      dplyr::group_by(sample, position) %>%
      dplyr::summarise(n = n(),
                       mean = mean(score, na.rm = T),
                       sd = sd(score, na.rm = T),
                       .groups = "keep") %>%
      dplyr::mutate(sem = sd / sqrt(n)) %>%
      dplyr::mutate(position = as.numeric(position),
                    sample = factor(sample, levels = samples))



    # Make plot
    plot =
      ggplot(data = mat.plot,
             aes(x = position,
                 y = mean,
                 ymin = mean - sem,
                 ymax = mean + sem,
                 color = sample,
                 fill = sample)) +
      geom_ribbon(alpha = 0.15, color = NA) +
      geom_line() +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      ylab("mean Footprint score \u00b1 SEM") +
      xlab("Distance from motif center [bp]") +
      ggtitle(factor) +
      #ylim(c(0,1)) +
      theme_classic() +
      theme(axis.text = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5, color = "black"),
            axis.ticks = element_line(color = "black"),
            aspect.ratio = 1)

    return(plot)
    # pdf(file = plot.file, height = 7, width = 7)
    # print(plot)
    # invisible(dev.off())
  }

