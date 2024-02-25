#' Plot log2 fold change
#'
#' This method plots the log2 fold change for each metabolite.
#'
#' @param metalyzer_se A Metalyzer object
#' @param signif_colors signif_colors
#' @param hide_labels_for vector of Metabolites or Classes for which no labels
#' are printed
#' @param class_colors class_colors
#' @param polarity_file polarity_file
#' @param vulcano boolean value to plot a vulcano plot
#'
#' @return ggplot object
#'
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import SummarizedExperiment
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_mutation_data_xl())
#' metalyzer_se <- filterMetabolites(
#'   metalyzer_se,
#'   drop_metabolites = "Metabolism Indicators"
#' )
#' metalyzer_se <- renameMetaData(
#'   metalyzer_se,
#'   Mutant_Control = "Sample Description"
#' )
#' metalyzer_se <- calculate_log2FC(
#'   metalyzer_se,
#'   categorical = "Mutant_Control",
#'   impute_perc_of_min = 0.2,
#'   impute_NA = TRUE
#' )
#' 
#' # p_vulcano <- plot_log2FC(metalyzer_se, vulcano=TRUE)
#' # p_fc <- plot_log2FC(metalyzer_se, vulcano=FALSE)
plot_log2FC <- function(metalyzer_se,
                        signif_colors=c("#5F5F5F"=1,
                                        "#FEBF6E"=0.1,
                                        "#EE5C42"=0.05,
                                        "#8B1A1A"=0.01),
                        hide_labels_for=c(),
                        class_colors="MetAlyzer",
                        polarity_file="MxPQuant500",
                        vulcano=FALSE) {
  log2FC_df <- metalyzer_se@metadata$log2FC

  ## Background: Load polarity data
  if (polarity_file == "MxPQuant500") {
    polarity_file <- polarity()
  }
  polarity_df <- utils::read.csv(polarity_file) %>%
    select(.data$Class,
           .data$Polarity) %>%
    mutate(Class = factor(.data$Class),
           Polarity = factor(.data$Polarity, levels = c('LC', 'FIA'))) %>%
    arrange(.data$Polarity)

  ## Background: Set class colors
  if (class_colors == "MetAlyzer") {
    class_colors <- metalyzer_colors()
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

  if (isFALSE(vulcano)) {
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

    ## Data: Determine labels
    signif_p_data <- filter(p_data, .data$signif_color != names(signif_colors)[1])
    if (length(hide_labels_for) > 0) {
      signif_p_data$Metabolite[which(signif_p_data$Metabolite %in% hide_labels_for)] <- NA
      signif_p_data$Metabolite[which(signif_p_data$Class %in% hide_labels_for)] <- NA
    }
    labels <- sapply(p_data$Metabolite, function(m) {
      m <- as.character(m)
      label <- if_else(m %in% signif_p_data$Metabolite, m, "")
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

  } else {

    ## Data: only color classes that are significantly differentially expressed
    log2FC_df$Class <- as.character(log2FC_df$Class)
    log2FC_df$Class[log2FC_df$qval > 0.05] <- NA
    log2FC_df$Class[abs(log2FC_df$log2FC) < log2(1.5)] <- NA

    ## Data: Determine labels
    log2FC_df$labels <- as.character(log2FC_df$Metabolite)
    log2FC_df$labels[which(is.na(log2FC_df$Class))] <- ""
    if (length(hide_labels_for) > 0) {
      log2FC_df$labels[which(log2FC_df$Metabolite %in% hide_labels_for)] <- ""
      log2FC_df$labels[which(log2FC_df$Class %in% hide_labels_for)] <- ""
    }

    ## Update lc_colors and fia_colors
    lc_colors <- lc_colors[which(names(lc_colors) %in% log2FC_df$Class)]
    fia_colors <- fia_colors[which(names(fia_colors) %in% log2FC_df$Class)]
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

  if (isFALSE(vulcano)) {
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

    # Create y-axis limits for the rectangles
    ylims <- c(min(log2FC_df$log2FC) - 0.75, max(log2FC_df$log2FC) + 0,75)

    ## Plot graph
    p_fc <- ggplot(p_data,
                   aes(x = .data$x,
                       y = .data$log2FC,
                       color = .data$signif_color,
                       label = labels)) +
      geom_rect(data = rects_df,
                inherit.aes = FALSE,
                aes(xmin = .data$Start, xmax = .data$End,
                    ymin = ylims[1], ymax = ylims[2],
                    fill = .data$Class),
                show.legend = TRUE,
                alpha = 0.4) +
      geom_vline(xintercept = 0, linewidth = 0.5, color = 'black') +
      geom_vline(xintercept = lc_fia_border+3, linewidth = 0.5, color = 'black',
                 linetype="dotted") +
      geom_hline(yintercept = 0, linewidth = 0.5, color = 'black') +
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
  } else {
    p_fc <- ggplot(log2FC_df,
                   aes(x = .data$log2FC,
                       y = -log10(.data$qval),
                       color = .data$Class,
                       label = labels)) +
      geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="black",
                 linetype="dashed") +
      geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed") +
      geom_point(size = 1) +
      scale_color_manual('Classes',
                        breaks = breaks,
                        values = values,
                        drop = FALSE,
                        guide = guide_legend(override.aes = list(size = 2),
                                             order=2, ncol = 2)) +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            legend.key = element_rect(fill = 'white')) +
      labs(x = 'log2(FC)', y = "-log10(p)") +
      geom_label_repel(size = 2, color = 'black',
                       box.padding = 0.6, # additional padding around each text label
                       point.padding = 0, # additional padding around each point
                       min.segment.length = 0,
                       max.overlaps = Inf,
                       force = 10)
  }
  return(p_fc)
}

