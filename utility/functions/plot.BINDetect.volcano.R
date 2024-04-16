plot.BINDetect.volcano =
  function(results,
           motif.pattern = "_(MOUSE|HUMAN)[.]H[0-9]*MO",
           right.color = "indianred",
           left.color = "steelblue",
           ns.color = "gray30",
           points.transparency = 0.5,
           labels.min.segment.length = 0,
           add.signif.lables = TRUE,
           labels.repel.force = 10,
           labels.max.overlaps = 100,
           extra.labels = "none") {

    # Libs
    require(ggplot2)
    require(dplyr)
    if (add.signif.lables == TRUE){require(ggrepel)}

    # Read table if required
    if ("character" %in% class(results)) {
      results = dplyr::mutate(read.delim(results, h=T), motif_id = gsub(motif.pattern, "",motif_id))
    }

    # Get conditions names
    condition.A = gsub("_bound", "", colnames(results)[7])
    condition.B = gsub("_bound", "", colnames(results)[9])
    title = paste(condition.A, "vs", condition.B)

    # Standardize table col names
    colnames(results)[10] = "change"
    colnames(results)[11] = "pvalue"
    colnames(results)[12] = "signif"

    # Add Factor column for points colors
    results = mutate(results,
                     diff.status = factor(ifelse(signif == "True",
                                                 yes = ifelse(sign(change) == 1,
                                                              yes = paste("Higher in", condition.A),
                                                              no = paste("Higher in", condition.B)),
                                                 no = "NS"),
                                          levels = c(paste("Higher in", condition.B),
                                                     paste("Higher in", condition.A),
                                                     "NS")))

    # Define colors
    colors = c(left.color, right.color, ns.color)
    names(colors) = c(paste("Higher in", condition.B),
                      paste("Higher in", condition.A),
                      "NS")


    # Generate the base plot
    plot =
      ggplot() +
      geom_point(data = results,
                 aes(x = change,
                     y = -log10(pvalue),
                     color = diff.status,
                     size = abs(change)),
                 stroke = NA,
                 alpha = points.transparency) +
      scale_color_manual(values = colors, name = "Differential binding status", drop = FALSE) +
      scale_size_continuous(name = "|Differential binding score|") +
      ggtitle(title) +
      ylab("-log~10~(*P*-value)") +
      xlab(paste0("Differential binding score [", condition.A, " - ", condition.B, "]")) +
      theme_classic() +
      theme(axis.title.y = ggtext::element_markdown(color = "black"),
            axis.text = ggtext::element_markdown(color = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black"))

    # make plot simmetric
    # get the x-range max
    x.max = max(abs(ggplot_build(plot)$layout$panel_params[[1]]$x.range))

    plot = plot + xlim(c(-1,1) * x.max)



    if (add.signif.lables == TRUE) {
      plot =
        plot +
        geom_text_repel(data = results %>% filter(diff.status != "NS" | motif_id %in% extra.labels),
                        mapping = aes(x = change,
                                      y = -log10(pvalue),
                                      label = motif_id,
                                      color = diff.status),
                        show.legend = F,
                        min.segment.length = labels.min.segment.length,
                        force = labels.repel.force,
                        max.overlaps = labels.max.overlaps)
    }

    return(plot)
  }

