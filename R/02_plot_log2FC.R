#' Plot log2 fold change
#'
#' This method plots the log2 fold change for each metabolite.
#'
#' @param log2FC_df fold change data frame
#' @param signif_colors signif_colors
#' @param class_colors class_colors
#' @param polarity_file polarity_file
#'
#' @return ggplot object
#'
#' @import ggplot2
#' @import ggrepel
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' print(1)
#' }

# classes <- c('Acylcarnitines', 'Alkaloids', 'Amine Oxides', 'Aminoacids',
#              'Aminoacids Related', 'Bile Acids', 'Biogenic Amines', 'Carboxylic Acids',
#              'Ceramides', 'Cholesterol Esters', 'Cresols', 'Diacylglycerols',
#              'Dihydroceramides', 'Fatty Acids', 'Glycerophospholipids',
#              'Glycosylceramides', 'Hormones', 'Indoles Derivatives',
#              'Nucleobases Related', 'Sphingolipids', 'Sugars', 'Triacylglycerols',
#              'Vitamins & Cofactors')
# class_colors <- c("#a6cee3", "#1f78b4", "#ffff99", "#b2df8a", "#33a02c", "#fb9a99",
#                   "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#8dd3c7",
#                   "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69",
#                   "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
# polarity <- c('FIA', 'LC', 'LC', 'LC', 'LC', 'LC', 'LC', 'LC', 'FIA', 'FIA',
#               'LC', 'FIA','FIA', 'LC', 'FIA',  'FIA', 'LC', 'LC', 'LC', 'FIA',
#               'FIA', 'FIA', 'LC')
# polarity_df <- data.frame(Class = classes,
#                           Polarity = factor(polarity, levels = c('LC', 'FIA'))) %>%
#   arrange(Polarity)
#
# polarity_file <- "inst/extdata/polarity.csv"

# signif_colors=c("#5F5F5F"=1, "#FEBF6E"=0.1, "#EE5C42"=0.05, "#8B1A1A"=0.01)
# class_colors="MetAlyzer"
# polarity_file="MxPQuant500"

plot_log2FC <- function(log2FC_df,
                        signif_colors=c("#5F5F5F"=1,
                                        "#FEBF6E"=0.1,
                                        "#EE5C42"=0.05,
                                        "#8B1A1A"=0.01),
                        class_colors="MetAlyzer",
                        polarity_file="MxPQuant500") {

  ## Background: Load polarity data
  if (polarity_file == "MxPQuant500") {
    polarity_file <- system.file("extdata", "polarity.csv", package = "MetAlyzer")
  }
  polarity_df <- utils::read.csv(polarity_file) %>%
    select(.data$Class,
           .data$Polarity) %>%
    mutate(Class = factor(.data$Class),
           Polarity = factor(.data$Polarity, levels = c('LC', 'FIA'))) %>%
    arrange(.data$Polarity)

  ## Background: Set class colors
  if (class_colors == "MetAlyzer") {
    class_colors <- readRDS(system.file("extdata", "metalyzer_colors.RDS", package = "MetAlyzer"))
  }
  names(class_colors) <- levels(polarity_df$Class)

  ## Background: Define LC and FIA classes with color
  lc_polarity_df <- filter(polarity_df,
                           .data$Polarity == 'LC',
                           .data$Class %in% log2FC_df$Class)
  lc_colors <- class_colors[which(names(class_colors) %in% lc_polarity_df$Class)]
  fia_polarity_df <- filter(polarity_df,
                            .data$Polarity == 'FIA',
                            .data$Class %in% log2FC_df$Class)
  fia_colors <- class_colors[which(names(class_colors) %in% fia_polarity_df$Class)]

  ## Data: Replace NAs
  log2FC_df$log2FC[is.na(log2FC_df$log2FC)] <- 0
  log2FC_df$qval[is.na(log2FC_df$qval)] <- 1

  ## Data: Add color to data based on significance
  log2FC_df$signif_color <- sapply(log2FC_df$qval, function(q_val) {
    for (t in signif_colors) {
      if (q_val <= t) {
        color <- names(signif_colors)[which(signif_colors == t)]
      }
    }
    return(color)
  })

  ## Data: Add pseudo x-value to data as a order of metabolites
  ordered_classes <- c(names(lc_colors), names(fia_colors))
  p_data <- lapply(ordered_classes, function(class) {
    log2FC_df %>%
      filter(.data$Class == class) %>%
      bind_rows(data.frame(Class = rep(NA, 5)))
  }) %>%
    bind_rows()
  p_data <- bind_rows(data.frame(Class = rep(NA, 5)), p_data)
  p_data$x <- seq(nrow(p_data))
  p_data <- filter(p_data, !is.na(.data$Class))

  ## Data: Determine lables
  signif_p_data <- filter(p_data, .data$signif_color != names(signif_colors)[1])
  lables <- sapply(p_data$Metabolite, function(m) {
    m <- as.character(m)
    label <- if_else(as.character(m) %in% signif_p_data$Metabolite, m, "")
    return(label)
  })

  ## Legend: Significance color
  signif_colors <- sort(signif_colors, decreasing = TRUE)
  signif_labels <- list()
  for (i in seq_along(signif_colors)) {
    t <- signif_colors[i]
    names(t) <- NULL
    if (i < length(signif_colors)) {
      t2 <- signif_colors[i+1]
      names(t2) <- NULL
      label <- bquote(.(t) ~ "\u2265 q-value >" ~ .(t2))
    } else {
      label <- bquote(.(t) ~ "\u2265 q-value")
    }
    signif_labels[[i]] <- label
  }

  ## Legend: Manage breaks and values for background rects
  len_diff <- length(lc_colors) - length(fia_colors)
  if (len_diff != 0) {
    blank_names <- sapply(1:abs(len_diff), function(i) {
      paste(rep(' ', i), collapse = '')
    })
    extension <- rep("white", abs(len_diff))
    names(extension) <- blank_names
    if (len_diff > 0) {
      # more classes from lc than fia
      # -> extend fia colors
      fia_colors <- c(fia_colors, extension)
    } else if (len_diff < 0) {
      # more classes from fia than lc
      # -> extend lc colors
      lc_colors <- c(lc_colors, extension)
    }
  }
  breaks <- c('LC:', names(lc_colors), 'FIA:', names(fia_colors))
  values <- c('white', lc_colors, 'white', fia_colors)
  names(values) <- NULL

  ## Background: Create data for background rects
  rects_df <- p_data %>%
    group_by(.data$Class) %>%
    summarise(Start = min(.data$x)-1,
              End = max(.data$x)+1,
              Color = class_colors[unique(.data$Class)])
  rects_df$Class <- factor(rects_df$Class, levels = breaks)

  ## Background: Determine border line between last LC and first FIA class
  lc_fia_border <- p_data %>%
    filter(.data$Class %in% names(lc_colors)) %>%
    select(.data$x) %>%
    max()

  ## Plot graph
  p_fc <- ggplot(p_data,
                 aes(x = .data$x,
                     y = .data$log2FC,
                     color = .data$signif_color,
                     label = lables)) +
    geom_rect(data = rects_df,
              inherit.aes = FALSE,
              aes(xmin = .data$Start, xmax = .data$End,
                  ymin = -Inf, ymax = Inf,
                  fill = .data$Class),
              show.legend = TRUE,
              alpha = 0.4) +
    geom_vline(xintercept = 0, size = 0.5, color = 'black') +
    geom_vline(xintercept = lc_fia_border+3, size = 0.5, color = 'black', linetype="dotted") +
    geom_hline(yintercept = 0, size = 0.5, color = 'black') +
    geom_point(size = 0.5) +
    scale_color_manual(paste0('Significance\n(linear model fit with FDR correction)'),
                       labels = signif_labels,
                       values = names(signif_colors),
                       guide = guide_legend(order=1)) +
    scale_fill_manual('Classes',
                      breaks = breaks,
                      values = values,
                      drop = FALSE,
                      guide = guide_legend(override.aes = list(alpha = 0.5),
                                           order=2, ncol = 2)) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          legend.key = element_rect(fill = 'white'),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line('#ECECEC'),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_line('#ECECEC'),
          panel.background = element_blank()) +
    labs(x = 'Metabolites') +
    geom_label_repel(size = 2, color = 'black',
                     box.padding = 0.6, # additional padding around each text label
                     point.padding = 0, # additional padding around each point
                     min.segment.length = 0,
                     max.overlaps = Inf,
                     force = 10)
  return(p_fc)
}
