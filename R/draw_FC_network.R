#' Draw metabolic pathway
#'
#' This function ...
#'
#' @param log2FC_df log2FC_df
#' @param output output
#' @param q_value q_value
#' @param pathway_file pathway_file
#' @param figsize figsize
#' @param colbar_width colbar_width
#' @param bg_color bg_color
#' @param node_size node_size
#' @param font_size font_size
#' @param font_color font_color
#' @param edge_color edge_color
#' @param ann_font_size ann_font_size
#' @param ann_font_color ann_font_color
#' @param colors colors
#' @param path_to_conda path_to_conda
#'
#' @import dplyr
#' @importFrom utils install.packages
#' @importFrom utils installed.packages
#' @importFrom stringr str_trim
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' print(1)
#' }

# log2FC_df <- readRDS("../change_test_df.RDS") %>%
#   filter(Tissue == "Zebrafish Liver")
# output <- "../pathway.png"
# pathway_file=fpath <- "inst/extdata/pathway.xlsx"
#
# draw_network(log2FC_df, output, pathway_file=pathway_file)

draw_FC_network <- function(
    log2FC_df,
    output,
    q_value=0.05,
    figsize=c(40, 36),
    dpi=300,
    bg_color="#FFFFFF", # white
    na_color="#D3DDDC", # lightgrey
    node_size=500,
    font_size=15,
    font_color="#000000", # black
    edge_color="#757575", # grey
    ann_font_size=25,
    ann_font_color="#046C9A",  # blue
    font_family="sans-serif",
    pathway_file=system.file("extdata", "pathway.xlsx", package = "MetAlyzer"),
    path_to_conda=NULL) {

  ## Check for reticulate installation
  installed_packages <- utils::installed.packages()[, "Package"]
  if (! "reticulate" %in% installed_packages) {
    cat('Installing package "reticulate":\n')
    utils::install.packages("reticulate")
    cat("\n")
  }

  cat("Configuring Python...  ")
  ## Check for conda environment
  if (!is.null(path_to_conda)) {
    options(reticulate.conda_binary = path_to_conda)
  }
  envs <- tryCatch(reticulate::conda_list(),
                   error = function(e) {
                     cat("\nNo installation of conda could be found.",
                         paste0("Installing miniconda to ",
                                reticulate::miniconda_path()), "\n")
                     reticulate::install_miniconda(
                       path = reticulate::miniconda_path(),
                       update = TRUE,
                       force = FALSE
                     )
                   },
                   finally = reticulate::conda_list()
  )
  if (!"MetAlyzer_networkx_env" %in% envs$name) {
    reticulate::conda_create(
      envname = "MetAlyzer_networkx_env",
      packages = c("networkx", "matplotlib", "pandas", "openpyxl", "regex")
    )
    envs <- reticulate::conda_list()
  }

  ## Initialize Python
  python_dir <- envs$python[envs$name == "MetAlyzer_networkx_env"]
  tryCatch(
    reticulate::use_python(python_dir),
    error = function(e) {
      # cat("\nPython version has already been initialized.\n")
      ""
    }
  )

  ## Source python script for network visualization
  reticulate::source_python(system.file("python", "draw_network.py",
                                        package = "MetAlyzer"))
  cat("finished!\n")

  ## Read pathway file
  nodes_df <- openxlsx::read.xlsx(pathway_file, sheet="Nodes") %>%
    select(.data$Label,
           .data$x,
           .data$y,
           .data$Shape,
           .data$Metabolites)
  nodes_df$Label <- stringr::str_trim(nodes_df$Label)

  edge_df <- openxlsx::read.xlsx(pathway_file, sheet="Edges") %>%
    select(.data$Node1,
           .data$Node2,
           .data$Rad)
  edge_df$Node1 <- stringr::str_trim(edge_df$Node1)
  edge_df$Node2 <- stringr::str_trim(edge_df$Node2)

  annotation_df <- openxlsx::read.xlsx(pathway_file, sheet="Annotations") %>%
    select(.data$Annotation,
           .data$x,
           .data$y)
  annotation_df$Annotation <- stringr::str_trim(annotation_df$Annotation)

  ## Remove invalid nodes, edges or annotations (= NA)
  invalid_nodes <- which(rowSums(is.na(nodes_df[,1:4])) > 0)
  n_invalid_nodes <- length(invalid_nodes)

  invalid_edges <- which(rowSums(is.na(edge_df)) > 0)
  n_invalid_edges <- length(invalid_edges)

  invalid_annotation <- which(rowSums(is.na(annotation_df)) > 0)
  n_invalid_annotation <- length(invalid_annotation)

  if (n_invalid_nodes > 0) {
    cat("Dropping invalid nodes (", n_invalid_nodes, "):\n  - ",
        paste(nodes_df$Label[invalid_nodes], collapse = "\n  - "), "\n",
        sep = "")
    nodes_df <- nodes_df[-invalid_nodes,]
  }
  if (n_invalid_edges) {
    cat("Dropping invalid edges (", n_invalid_edges, "):\n  - ",
        paste(sapply(invalid_edges, function(i) {
          paste(edge_df[i, 1:2], collapse = " <-> ")
        }), collapse = "\n  - "), "\n", sep = "")
    edge_df <- edge_df[-invalid_edges,]
  }
  if (n_invalid_annotation) {
    cat("Dropping invalid annotations (", n_invalid_annotation, "):\n  - ",
        paste(annotation_df$Annotation[invalid_annotation], collapse = "\n  - "),
        "\n", sep = "")
    annotation_df <- annotation_df[-invalid_annotation,]
  }

  ## Print nodes with no assigned metabolite
  not_assigned_nodes <- which(is.na(nodes_df$Metabolites))
  n_not_assigned_nodes <- length(not_assigned_nodes)

  if (n_not_assigned_nodes > 0) {
    cat("Nodes which are not assigned to any metabolite (", n_not_assigned_nodes, "):\n  - ",
        paste(nodes_df$Label[not_assigned_nodes], collapse = "\n  - "), "\n",
        sep = "")
  }

  ## Remove not assignable edges
  mismatched_edges <- which(sapply(1:nrow(edge_df), function(i) {
    any(!edge_df$Node1[i] %in% nodes_df$Label,
        !edge_df$Node2[i] %in% nodes_df$Label)
  }))
  n_mismatched_edges <- length(mismatched_edges)

  if (n_mismatched_edges > 0) {
    cat("Edges with at least one missing node (", n_mismatched_edges, "):\n  - ",
        paste(sapply(mismatched_edges, function(i) {
          paste(edge_df[i, 1:2], collapse = " <-> ")
        }), collapse = "\n  - "), "\n", sep = "")
    edge_df <- edge_df[-mismatched_edges,]
  }

  ## Print metabolites that could not been found in input data
  network_metabos <- nodes_df$Metabolites[!is.na(nodes_df$Metabolites)]
  network_metabos <- stringr::str_trim(unlist(strsplit(network_metabos, ";")))
  network_metabos <- unique(network_metabos)
  unfound_metabos <- which(!network_metabos %in% levels(log2FC_df$Metabolite))
  n_unfound_metabos <- length(unfound_metabos)

  if (n_unfound_metabos > 0) {
    cat("Metabolites which could not be found in input (Metabolite levels) (",
        n_unfound_metabos, "):\n  - ",
        paste(network_metabos[unfound_metabos], collapse = "\n  - "), "\n",
        sep = "")
  }

  ## Print un-mapped metabolites
  unmapped_metabos <- which(!levels(log2FC_df$Metabolite) %in% network_metabos)
  n_unmapped_metabos <- length(unmapped_metabos)

  if (n_unmapped_metabos > 0) {
    cat("Metabolites which are not included in the network (",
        n_unmapped_metabos, "):\n  - ",
        paste(levels(log2FC_df$Metabolite)[unmapped_metabos],
              collapse = "\n  - "), "\n", sep = "")
  }

  ## Add log2FC to nodes_df
  signif_df <- filter(log2FC_df,
                      !is.na(.data$log2FC),
                      !is.na(.data$qval),
                      .data$qval <= q_value)

  nodes_df$log2FC <- sapply(strsplit(nodes_df$Metabolites, ";"), function(m_vec) {
    if (length(m_vec) > 1) {
      # Nodes with more than 1 metabolite assigned
      tmp_df <- filter(signif_df, .data$Metabolite %in% m_vec)
      if (nrow(tmp_df) > 0) {
        # At least one of the metabolites is significantly change
        # -> take the mean log2 fold change
        l2fc <- sum(tmp_df$log2FC) / length(m_vec)
      } else {
        if (any(tmp_df$Metabolite %in% levels(log2FC_df$Metabolite))) {
          # At least one metabolite was measured but none are significantly changed
          l2fc <- 0
        } else {
          # None of the metabolites were measured
          l2fc <- NA
        }
      }
    } else {
      # Nodes with 0 or 1 metabolite assigned
      if (m_vec %in% signif_df$Metabolite) {
        # Metabolite is significantly changed
        l2fc <- signif_df$log2FC[which(signif_df$Metabolite == m_vec)]
      } else if (m_vec %in% levels(log2FC_df$Metabolite)) {
        # Metabolite was measured but is not significantly changed
        l2fc <- 0
      } else {
        # Metabolite was not measured
        l2fc <- NA
      }
    }
    return(l2fc)
  })

  ## Add node color based on log2FC
  colors <- list(
    "3 ≤ FC      "= "#FF0000",
    "1.5 ≤ FC < 3   "= "#FF6666",
    "0.5 ≤ FC < 1.5"= "#FF9999",
    "No fold change"= "#A0A0A0",
    "Not measured"= na_color,
    "-1.5 < FC ≤ -0.5"= "#00CCCC",
    "   -3 < FC ≤ -1.5"= "#3399FF",
    "       FC ≤ -3"= "#0066CC"
  )

  nodes_df$Color <- ""
  nodes_df$Color[is.na(nodes_df$log2FC)] <- colors["Not measured"]
  nodes_df$Color[nodes_df$log2FC >= 3] <- colors["3 ≤ FC      "]
  nodes_df$Color[nodes_df$log2FC < 3] <- colors["1.5 ≤ FC < 3   "]
  nodes_df$Color[nodes_df$log2FC < 1.5] <- colors["0.5 ≤ FC < 1.5"]
  nodes_df$Color[nodes_df$log2FC < 0.5] <- colors["No fold change"]
  nodes_df$Color[nodes_df$log2FC <= -0.5] <- colors["-1.5 < FC ≤ -0.5"]
  nodes_df$Color[nodes_df$log2FC <= -1.5] <- colors["   -3 < FC ≤ -1.5"]
  nodes_df$Color[nodes_df$log2FC <= -3] <- colors["       FC ≤ -3"]
  nodes_df$Color <- as.character(nodes_df$Color)

  ## Add nodes of legend to nodes_df
  legend_df <- data.frame(matrix(nrow = length(colors), ncol = ncol(nodes_df)))
  colnames(legend_df) <- colnames(nodes_df)
  legend_df$Label <- names(colors)
  legend_df$x <- max(nodes_df$x, na.rm=TRUE) * 1.1
  legend_df$Shape <- "o"
  legend_df$Color <- as.character(colors)
  span_y <- (max(nodes_df$y, na.rm=TRUE) - min(nodes_df$y, na.rm=TRUE))
  legend_y_center <- span_y / 2
  legend_y_step <- span_y * 0.05
  legend_ys <- legend_y_center + seq(from=3.5 * legend_y_step,
                                     to=-3.5 * legend_y_step,
                                     by=-legend_y_step)
  legend_df$y <- round(legend_ys)
  # legend_df$Metabolites <- ""
  # legend_df$log2FC <- 0

  # nodes_legend_df <- rbind(nodes_df, legend_df)


  ## Draw network
  cat("Drawing network...  ")
  draw_network_py(
    nodes_df=rbind(nodes_df, legend_df),
    edge_df=edge_df,
    annotation_df=annotation_df,
    output=output,
    f_width=figsize[1],
    f_height=figsize[2],
    dpi=dpi,
    bg_color=bg_color,
    node_size=node_size,
    font_size=font_size,
    font_color=font_color,
    font_family=font_family,
    edge_color=edge_color,
    ann_font_size=ann_font_size,
    ann_font_color=ann_font_color
  )
  cat("finished!\n")
}
