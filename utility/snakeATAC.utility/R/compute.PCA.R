#' @title plot.PCA
#'
#' @description Computes and plots PCA starting from a multiBigwigSummary file.
#'
#' @param matrix Matrix where each column is a sample and each row a genomic bin/region. A multiBigwigSummary file can be read using \link{read.multibigwigsummary}.
#' @param excluded.samples String vector indicating the list of sample IDs to exclude from the plot. Default: \code{NULL}.
#' @param scale Logical value to indicate whether the data should be scaled. Default: \code{TRUE}.
#' @param center Logical value to indicate whether the data should be centered. Default: \code{TRUE}.
#' @param PCs Numeric vector of length 2, indicating which principal components (PCs) should be plotted. Default: \code{c(1,2)}.
#' @param point.size Numeric value indicating the dot size. Default: \code{3}.
#' @param point.alpha Numeric value between 0-1 indicating the dot transparency (1 full, 0 invisible). Default: \code{0.75}.
#' @param color.list (Named) String vector indicating any R-supported color to assign to each sample. Default: \code{NULL} (automatic).
#' @param shape.list Vector indicating any R-supported shape to assign to each sample. Default: \code{NULL} (automatic).
#'
#' @import dplyr
#' @import ggplot2
#'
#' @return A ggplot2 object.
#'
#' @export plot.PCA

plot.PCA =
  function(matrix,
           excluded.samples = NULL,
           scale = TRUE,
           center = TRUE,
           PCs = c(1,2),
           point.size = 3,
           point.alpha = 0.75,
           color.list = NULL) {

    ### Exclude samples if required
    if (!is.null(excluded.samples)) {
      if (FALSE %in% (excluded.samples %in% colnames(matrix))) {
        stop(paste0("At least one of the 'excluded samples' is not present in the ID list.\n",
             "Missing 'excluded.samples': \n", paste0(excluded.samples[!(excluded.samples %in% colnames(matrix))], collapse = ", "), ".\n",
             "IDs available: \n", paste0(sort(colnames(matrix)), collapse = ", "), "."))
      }
      summary.matrix = matrix[, !(colnames(matrix) %in% excluded.samples)]
    } else {
      summary.matrix = matrix
    }


    ### Perform PCA analyses
    pca = prcomp(t(summary.matrix),
                 scale. = scale,
                 center = center,
                 retx = TRUE)

    pca.summary = data.frame(summary(pca)$importance)


    ### generate PCA plot
    ## filter the table for the PCA of interest
    pca.tb =
      as.data.frame(pca$x) %>%
      dplyr::select(dplyr::all_of(paste0("PC", PCs[1:2]))) %>%
      dplyr::mutate(ID = rownames(pca$x))

    colnames(pca.tb) = c("x","y","ID")


    ### Add formatting
    pca.tb = dplyr::mutate(pca.tb, color = ID, shape = "Global")
    show.legend.shape = FALSE


    if (is.null(color.list)) {
      color.list = rainbow(n = nrow(pca.tb))
      names(color.list) = pca.tb$ID
    }

    if (is.null(shape.list)) {
      shape.list = c("Global" = 1)
    }


    ## generate the PCA plot
    pca.plot =
      ggplot(data = pca.tb,
             aes(x = x,
                 y = y,
                 color = color,
                 shape = shape)) +
      geom_hline(yintercept = 0, linetype = 2, color = "gray") +
      geom_vline(xintercept = 0, linetype = 2, color = "gray") +
      geom_point(stroke = NA, size = point.size, alpha = point.alpha) +
      scale_color_manual(values = color.list, name = NULL) +
      #scale_shape_manual(values = shape.list, name = NULL) +
      theme_classic() +
      xlab(paste0("PC", PCs[1], " (", round(pca.summary[2,PCs[1]]*100,1), "%)")) +
      ylab(paste0("PC", PCs[2], " (", round(pca.summary[2,PCs[2]]*100,1), "%)")) +
      theme(axis.text = element_text(colour = "black"),
            axis.line = element_blank(),
            axis.ticks = element_line(color = "black"),
            panel.background = element_rect(fill = NA, color = "black"),
            aspect.ratio = 1)

    if (show.legend.shape == FALSE) {
      pca.plot = pca.plot + guides(shape = "none")
    }


    return(pca.plot)
  }
