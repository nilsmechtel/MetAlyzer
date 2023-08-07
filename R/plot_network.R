# Load the required packages
library(dplyr)
library(ggplot2)
library(ggrepel)


pathway_file <- "/Users/nilsm/Library/CloudStorage/OneDrive-bwedu/Uni/HiWi/HiWi COS/Package - MetAlyzer/MetAlyzer/inst/extdata/pathway.xlsx" # MetAlyzer::pathway()

read_named_region <- function(file_path, named_region) {
  full_sheet <- openxlsx::read.xlsx(
    file_path,
    sheet = 1,
    colNames = FALSE,
    skipEmptyRows = FALSE,
    skipEmptyCols = FALSE,
  )
  full_sheet[nrow(full_sheet) + 1, ] <- NA
  header <- colnames(openxlsx::read.xlsx(
    file_path,
    namedRegion = named_region
  ))
  coordinates <- lapply(header, function(col_name) {
    data.frame(which(full_sheet == col_name, arr.ind = TRUE))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(row, col) %>%
    dplyr::group_by(row) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::filter(n == length(header))

  header_row <- unique(coordinates$row)
  first_row <- header_row + 1
  cols <- coordinates$col
  df <- full_sheet[
    first_row:nrow(full_sheet),
    cols
  ]
  colnames(df) <- header
  last_row <- min(which(rowSums(is.na(df)) == length(header))) - 1
  df <- df[1:last_row, ]

  if ("x" %in% header) {
    df$x <- as.numeric(df$x)
  }
  if ("y" %in% header) {
    df$y <- as.numeric(df$y)
  }
  if ("Radius" %in% header) {
    df$Radius <- as.numeric(df$Radius)
  }
  rownames(df) <- NULL
  return(df)
}


## Read network nodes, edges and annotations
pathways <- read_named_region(pathway_file, "Pathways_Header")
pathways$Label <- stringr::str_trim(pathways$Label)
invalid_annotations <- which(
  is.na(pathways$Label) |
  duplicated(pathways$Label) |
  is.na(pathways$x) |
  is.na(pathways$y) |
  is.na(pathways$Color)
)
if (length(invalid_annotations) > 0) {
  # print warning and remove
  cat("Warning: Removing", length(invalid_annotations), "invalid pathways.\n")
  pathways <- pathways[-invalid_annotations, ]
}
rownames(pathways) <- pathways$Label

nodes <- read_named_region(pathway_file, "Metabolites_Header")
nodes$Label <- stringr::str_trim(nodes$Label)
nodes$Pathway[is.na(nodes$Pathway)] <- ""
invalid_nodes <- which(
  is.na(nodes$Label) |
  duplicated(nodes$Label) |
  is.na(nodes$x) |
  is.na(nodes$y) |
  !nodes$Pathway %in% c(pathways$Label, "")
)
if (length(invalid_nodes) > 0) {
  # print warning and remove
  cat("Warning: Removing", length(invalid_nodes), "invalid nodes.\n")
  nodes <- nodes[-invalid_nodes, ]
}
rownames(nodes) <- nodes$Label

edges <- read_named_region(pathway_file, "Connections_Header")
edges$Node1 <- stringr::str_trim(edges$Node1)
edges$Node2 <- stringr::str_trim(edges$Node2)
invalid_edges <- which(
  !edges$Node1 %in% nodes$Label |
  !edges$Node2 %in% nodes$Label |
  edges$Node1 == edges$Node2
)
if (length(invalid_edges) > 0) {
  # print warning and remove
  cat("Warning: Removing", length(invalid_edges), "invalid connections.\n")
  edges <- edges[-invalid_edges, ]
}


edges$x_start <- nodes[edges$Node1, "x"]
edges$y_start <- nodes[edges$Node1, "y"]
edges$x_end <- nodes[edges$Node2, "x"]
edges$y_end <- nodes[edges$Node2, "y"]
edges$Pathway <- sapply(rownames(edges), function(rowname) {
  from <- edges[rowname, "Node1"]
  to <- edges[rowname, "Node2"]
  from_pathway <- nodes[from, "Pathway"]
  to_pathway <- nodes[to, "Pathway"]
  color <- "grey"
  if (from_pathway == to_pathway & from_pathway != "") {
    color <- pathways[from_pathway, "Color"]
  }
  return(color)
})


## Add log2FC to nodes
nodes$FC_thresh <- floor(runif(rownames(nodes), -3, 4)) # !

## Draw network

# Create a plot of the network using ggplot2 and ggrepel
label_size <- 3
area_size <- 5
edge_size <- 1
annotation_size <- 6

network <- ggplot()
for (radius in unique(edges$Radius)) {
  tmp_edges <- filter(edges, Radius == radius)

  network <- network +
    # Add the round area behind the edges
    geom_curve(
      data = tmp_edges,
      aes(
        x = x_start,
        y = y_start,
        xend = x_end,
        yend = y_end,
        color = Pathway
      ),
      # color = "lightblue",
      size = area_size,
      alpha = 0.3,
      curvature = radius,
      show.legend = FALSE
    ) +
    # Or add the edges as curved lines
    geom_curve(
      data = tmp_edges,
      aes(
        x = nodes[Node1, "x"],
        y = nodes[Node1, "y"],
        xend = nodes[Node2, "x"],
        yend = nodes[Node2, "y"]
      ),
      color = "grey", size = edge_size, curvature = radius
    )
}
network <- network +
  # Add labels at the position of the nodes
  geom_label(
    data = nodes,
    aes(
      x = x,
      y = y,
      label = Label,
      fill = FC_thresh
    ),
    size = label_size,
    color = "white"
  ) +
  # Add annotations
  geom_text(
    data = pathways,
    aes(
      x = x,
      y = y,
      label = Label,
      color = Color
    ),
    size = annotation_size,
    show.legend = FALSE
  ) +
  # # Set the x and y axis limits
  # xlim(0, 10) +
  # ylim(0, 10) +
  theme_void() +
  # Add a title and remove the x and y axis labels
  ggtitle("Example Network Plot with Colored Area Behind Curved Edges") +
  theme(plot.title = element_text(hjust = 0.5))
network

ggsave("network.png", network, width = 15, height = 10, bg = "white")
