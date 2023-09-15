#' Collection of functions for drawing.
#' 
#' 
library(openxlsx)
source("./DEgene-network.r")


# For MainFigs.
#
# The main figure was generated using function "generate.DEgeneNetwork.main()"
# and adjusted in Illustrator.
#
if (FALSE) {
    # Results of KEGG-mining.
    #
    #
    DEgene_data   <- read.xlsx("../data/DEgene_data.xlsx")
    edge_data     <- read.xlsx("../data/edge_data.xlsx")
    upstream_data <- read.xlsx("../data/upstream_data.xlsx")


    # Drawing.
    #
    #
    nw <- generate.DEgeneNetwork.main(
              DEgene_data,
              edge_data,
              upstream_data,
              DEgene_condition = "F0-3",
              edge_condition   = "F0-3",
              organ_area.alpha = 0.35,
              organ_area.radius   = c(12, 21),
              organ_area.n_points = 30L,
              organ_label.color   = "#808080",
              organ_label.radius  = 18.5,
              organ_label.size    = 17,
              DEgene.colors       = c("F0-3"="#00B0F6"),
              DEgene.radius       = c(12.2, 20.8),
              DEgene.size         = 3,
              DEgene_label.offset = 0.2,
              DEgene_label.size   = 5,
              DEgene_legend       = "DEgene condition",
              ligand.radius       = 27,
              ligand.size         = 3,
              ligand.shape        = 16,
              ligand_label.offset = 0.2,
              ligand_label.size   = 5,
              upstream_edge.width    = 1,
              upstream_edge.n_points = 50L,
              edge.colors   = c("F0-3"="#00B0F6"),
              edge.width    = 1,
              edge.n_points = 50L,
              edge_legend   = "Edge condition",
              gap.degree = 2,
              arrow      = grid::arrow(angle=15, length=unit(8,"mm"), type="closed")
    )

    ggsave("../data/Main-style.pdf", plot=nw,
           limitsize=FALSE, units="mm", height=800, width=800)
}


# For ExtDataFigs.
#
# The extended data figure was a merging of two figures generated
# using function "generate.DEgeneNetwork.ext()".
# Merging and adjustment were performed in Illustrator.
#
if (FALSE) {
    library(circlize)
    library(scales)

    # Results of KEGG-mining.
    #
    #
    DEgene_data   <- read.xlsx("../data/DEgene_data.xlsx")
    edge_data     <- read.xlsx("../data/edge_data.xlsx")
    upstream_data <- read.xlsx("../data/upstream_data.xlsx")


    # Drawing
    #
    #
    nw <- generate.DEgeneNetwork.ext(
              DEgene_data,
              edge_data,
              upstream_data,
              DEgene_condition = "F0-3",
              edge_condition   = "F0-3",
              organ_area.alpha = 0.35,
              organ_area.radius   = c(12, 21),
              organ_area.n_points = 30L,
              organ_label.color   = "#808080",
              organ_label.radius  = 18.5,
              organ_label.size    = 17,
              DEgene.radius       = c(12, 21),
              ligand.radius       = 27,
              ligand.size         = 3,
              ligand.shape        = 16,
              ligand_label.offset = 0.2,
              ligand_label.size   = 5,
              upstream_edge.width    = 1,
              upstream_edge.n_points = 50L,
              gap.degree = 2
    )

    ggsave("./data/ExtData-style_1.pdf", plot=nw$background,
           limitsize=FALSE, units="mm", height=800, width=800)


    circos.par("clock.wise"=FALSE, gap.degree=2)

    pdf("./data/ExtData-style_2.pdf")

    chordDiagram(
        nw$edge,
        directional  = TRUE,
        transparency = 0.25,
        row.col  = hue_pal()(nrow(nw$edge)),
        grid.col = hue_pal()(nrow(nw$edge)),
        scale    = TRUE,
        annotationTrack = c("name","grid"),
        annotationTrackHeight = c(mm_h(1), mm_h(4)),
        direction.type = c("diffHeight", "arrows"),
        link.arr.type  = "big.arrow",
        reduce = -1
    )

    dev.off()
    circos.clear()
}