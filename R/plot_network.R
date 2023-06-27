# Load the required packages
library(dplyr)
library(ggplot2)
library(ggrepel)

pathway_file <- pathway()

# Prepare network input ---------------------------------------------------

## Read pathway file
nodes <- openxlsx::read.xlsx(pathway_file, sheet="Nodes") %>%
  select(Label,
         x,
         y,
         Shape,
         Metabolites)
nodes$Label <- stringr::str_trim(nodes$Label)
colnames(nodes)[1] <- "ID"
colnames(nodes)[2] <- "Position_x"
colnames(nodes)[3] <- "Position_y"
nodes$Position_x <- nodes$Position_x / 10
nodes$Position_y <- nodes$Position_y / 10
nodes$Pathway <- NA
nodes$Pathway[1:30] <- "Bile Acids"

edges <- openxlsx::read.xlsx(pathway_file, sheet="Edges") %>%
  select(Node1,
         Node2,
         Rad)
edges$Node1 <- stringr::str_trim(edges$Node1)
edges$Node2 <- stringr::str_trim(edges$Node2)
colnames(edges)[1] <- "from_ID"
colnames(edges)[2] <- "to_ID"
colnames(edges)[3] <- "Curvature"
edges$from <- sapply(edges$from_ID, function(id) {
  which(nodes$ID == id)
})
edges$to <- sapply(edges$to_ID, function(id) {
  which(nodes$ID == id)
})

pathways <- openxlsx::read.xlsx(pathway_file, sheet="Annotations") %>%
  select(Annotation,
         x,
         y)
pathways$Annotation <- stringr::str_trim(pathways$Annotation)
colnames(pathways)[1] <- "Pathway"
colnames(pathways)[2] <- "Position_x"
colnames(pathways)[3] <- "Position_y"
pathways$Position_x <- pathways$Position_x / 10
pathways$Position_y <- pathways$Position_y / 10
pathways$Color <- rownames(pathways)

nodes <- filter(nodes,
                !is.na(Position_x),
                !is.na(Position_y))
nodes$FC_thresh <- floor(runif(rownames(nodes), -3, 4)) #!

nodes$Pathway[is.na(nodes$Pathway)] <- ""
edges$Pathway <- sapply(seq_along(rownames(edges)), function(i) {
  from <- edges[i, "from"]
  to <- edges[i, "to"]
  from_pathway <- nodes[from, "Pathway"]
  to_pathway <- nodes[to, "Pathway"]
  color <- ifelse(all(c(from_pathway == to_pathway),
                      from_pathway %in% pathways$Pathway),
                  pathways$Color[which(pathways$Pathway == from_pathway)],
                  "grey")
                  return(color)
})

# Draw network ------------------------------------------------------------


# Create a plot of the network using ggplot2 and ggrepel
label_size <- 3
area_size <- 5
edge_size <- 1
annotation_size <- 6

network <- ggplot()
for (curvature in unique(edges$Curvature)) {
  tmp_edges <- filter(edges, Curvature == curvature)
  network <- network +
    # Add the round area behind the edges
    geom_curve(data = tmp_edges,
               aes(x = nodes[from, "Position_x"],
                   y = nodes[from, "Position_y"],
                   xend = nodes[to, "Position_x"],
                   yend = nodes[to, "Position_y"],
                   color = Pathway),
               # color = "lightblue",
               size = area_size,
               alpha=0.3,
               curvature = curvature,
               show.legend = FALSE) +
    # Or add the edges as curved lines
    geom_curve(data = tmp_edges,
               aes(x = nodes[from, "Position_x"],
                   y = nodes[from, "Position_y"],
                   xend = nodes[to, "Position_x"],
                   yend = nodes[to, "Position_y"]),
               color = "grey", size = edge_size, curvature = curvature)
}
network <- network +
  # Add labels at the position of the nodes
  geom_label(data = nodes,
             aes(x = Position_x,
                 y = Position_y,
                 label = ID,
                 fill = FC_thresh),
             size = label_size,
             color = "white") +
  # Add annotations
  geom_text(data = pathways,
            aes(x = Position_x,
                y = Position_y,
                label = Pathway,
                color = Color),
            size = annotation_size,
            show.legend = FALSE) +
  # # Set the x and y axis limits
  # xlim(0, 10) +
  # ylim(0, 10) +
  theme_void() +
  # Add a title and remove the x and y axis labels
  ggtitle("Example Network Plot with Colored Area Behind Curved Edges") +
  theme(plot.title = element_text(hjust = 0.5))
network
