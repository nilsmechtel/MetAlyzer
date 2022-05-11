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
    log2FC_df, output,
    q_value=0.05,
    pathway_file=system.file("extdata", "pathway.xlsx", package = "MetAlyzer"),
    figsize=c(35, 35),
    colbar_width=2,
    bg_color="white",  # matplotlib
    # unassigned_color="dimgray",  # matplotlib
    node_size=500,
    font_size=15,
    font_color="black",  # matplotlib
    edge_color="grey",  # matplotlib
    ann_font_size=25,
    ann_font_color="blue",  # matplotlib
    colors=c("royalblue", "grey80", "tomato2"), # https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
    path_to_conda=NULL) {

  ## Check for reticulate installation
  installed_packages <- utils::installed.packages()[, "Package"]
  if (! "reticulate" %in% installed_packages) {
    cat("Installing package 'reticulate':\n")
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
  }

  ## Initialize Python
  dummy <- tryCatch(
    reticulate::use_condaenv(condaenv = "MetAlyzer_networkx_env",
                             required = TRUE),
    error = function(e) {
      cat("\nPython version has already been initialized.\n")
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
    cat("Not assigned nodes (", n_not_assigned_nodes, "):\n  - ",
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
    cat("Not assignable edges (", n_mismatched_edges, "):\n  - ",
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

  ## Print unmapped metabolites
  unmapped_metabos <- which(!levels(log2FC_df$Metabolite) %in% network_metabos)
  n_unmapped_metabos <- length(unmapped_metabos)

  if (n_unmapped_metabos > 0) {
    cat("Metabolites which are not included in the network (",
        n_unmapped_metabos, "):\n  - ",
        paste(levels(log2FC_df$Metabolite)[unmapped_metabos],
              collapse = "\n  - "), "\n", sep = "")
  }

  ## Add log2FC to nodes_df
  log2FC_df <- filter(log2FC_df,
                      !is.na(.data$log2FC),
                      !is.na(.data$qval),
                      .data$qval <= q_value)

  nodes_df$log2FC <- sapply(strsplit(nodes_df$Metabolites, ";"), function(m_vec) {
    if (length(m_vec) > 1) {
      tmp_df <- filter(log2FC_df, .data$Metabolite %in% m_vec)
      if (nrow(tmp_df) == 0) {
        value <- 0
      } else {
        value <- sum(tmp_df$log2FC) / length(m_vec)
      }
    } else if (m_vec %in% log2FC_df$Metabolite) {
      value <- log2FC_df$log2FC[which(log2FC_df$Metabolite == m_vec)]
    } else {
      value <- 0
    }
    return(value)
  })
  # nodes_df$log2FC[is.na(nodes_df$Metabolites)] <- NA

  ## Set last parameters
  f_width <- figsize[1]
  f_height <- figsize[2]
  palette <- grDevices::colorRampPalette(colors = colors)
  cmap <- palette(1000)

  cat("Drawing network...  ")
  draw_network_py(nodes_df, edge_df, annotation_df, output, f_width, f_height,
                  colbar_width, cmap, bg_color, "black", node_size,
                  font_size, font_color, edge_color, ann_font_size, ann_font_color)
  cat("finished!\n")
}
