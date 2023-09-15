#' Collection of functions for drawing.
#' 
#' 
library(tidyverse)
library(ggplot2)
library(ggtext)
library(geomtextpath)
library(ggnewscale)

source("./funcs_for_data-formatting.r")
source("./funcs_for_drawing.r")

#' Generate a network with DEgenes for MainFigs.
#' 
#' 
generate.DEgeneNetwork.main <- function (
    DEgene_data,
    edge_data,
    upstream_data,
    ...
) {
    draw_data <- generate.draw_data.main(
                     DEgene_data,
                     edge_data,
                     upstream_data,
                     ...)

    organs <- unique(draw_data$DEgene$organ)

    DEgene_data <- draw_data$DEgene[!is.na(draw_data$DEgene$id),]

    if (nrow(DEgene_data) > 0L) {
        DEgene_data <- add.new_columns.DEgene(
                           DEgene_data,
                           draw_data$edge,
                           organs,
                           ...)
    }

    if (nrow(draw_data$ligand) > 0L) {
        ligand_data <- add.new_columns.ligand.main(
                           draw_data$ligand,
                           organs,
                           ...)
    }

    # Network

    nw <- ggplot() +
          coord_fixed() +
          theme(panel.background = element_blank(),
                panel.grid       = element_blank(),
                axis.ticks       = element_blank(),
                axis.text        = element_blank(),
                axis.title       = element_blank(),
                axis.line        = element_blank())

    # Organ areas.
    nw <- draw.organ_areas(nw, organs, ...)

    # Edges from ligands to DEgenes.
    if (nrow(draw_data$upstream) > 0L) {
        nw <- draw.upstream_edges(
                  nw,
                  DEgene_data,
                  ligand_data,
                  draw_data$upstream,
                  ...)
    }

    # Edges between DEgenes.
    if (nrow(draw_data$edge) > 0L) {
        nw <- draw.DEgene_edges(
                  nw,
                  DEgene_data,
                  draw_data$edge,
                  ...)
    }

    # DEgenes.
    if (nrow(DEgene_data) > 0L) {
        nw <- draw.DEgenes(nw, DEgene_data, ...)
    }

    # Ligands
    if (nrow(ligand_data) > 0L) {
        nw <- draw.ligands.main(nw, ligand_data, ...)
    }

    # Edge legend.
    nw <- draw.edge_legend(nw, draw_data$edge, ...)

    return (nw)
}



#' Generate a network with DEgenes for extended data figures.
#' 
#' 
generate.DEgeneNetwork.ext <- function (
    DEgene_data,
    edge_data,
    upstream_data,
    ...
) {
    draw_data <- generate.draw_data.ext(
                     DEgene_data,
                     edge_data,
                     upstream_data,
                     ...)

    organs <- unique(draw_data$DEgene$organ)

    DEgene_data <- draw_data$DEgene[!is.na(draw_data$DEgene$id),]

    if (nrow(DEgene_data) > 0L) {
        DEgene_data <- add.new_columns.DEgene(
                           DEgene_data,
                           draw_data$edge,
                           organs,
                           ...)
    }

    if (nrow(draw_data$ligand) > 0L) {
        ligand_data <- add.new_columns.ligand.ext(
                           draw_data$ligand,
                           organs,
                           ...)
    }


    # Network

    # For background.
    nw <- ggplot() +
          coord_fixed() +
          theme(panel.background = element_blank(),
                panel.grid       = element_blank(),
                axis.ticks       = element_blank(),
                axis.text        = element_blank(),
                axis.title       = element_blank(),
                axis.line        = element_blank())

    # Organ areas.
    nw <- draw.organ_areas(nw, organs, ...)

    # Edges from ligands to DEgenes.
    if (nrow(draw_data$upstream) > 0L) {
        nw <- draw.upstream_edges(
                  nw,
                  DEgene_data,
                  ligand_data,
                  draw_data$upstream,
                  ...)
    }

    # Ligands
    if (nrow(ligand_data) > 0L) {
        nw <- draw.ligands.ext(nw, ligand_data, ...)
    }

    # Edge legend.
    nw <- draw.edge_legend(nw, draw_data$edge, ...)

    return (list(background=nw, edge=draw_data$matrix))
}





#' Add new columns to DEgene_data.
#' 
#' 
add.new_columns.DEgene <- function (
    DEgene_data,
    edge_data,
    organs,
    DEgene.radius = c(8.1, 9.9),
    DEgene_label.offset = 0.3,
    gap.degree    = 1,
    ...
) {
    # Find the inner nodes.
    edge_data.inter <- edge_data[edge_data$from.organ!=edge_data$to.organ,]
    DEgenes.inner     <- unique(c(edge_data.inter$from, edge_data.inter$to))

    # Add a new columns "inner".
    DEgene_data <- DEgene_data %>% mutate(inner = (DEgene_data$id %in% DEgenes.inner))


    # Calculate positions for DEgenes.
    n_organs <- length(organs)

    area_angle.width  <- 360 / n_organs
    area_angle.margin <- gap.degree / 2

    for (i in seq_len(n_organs)) {
        DEgene_data_i <- DEgene_data[DEgene_data$organ==organs[i],]

        DEgene_data_i.inner <- DEgene_data_i[DEgene_data_i$inner,]  %>% select("id")
        DEgene_data_i.outer <- DEgene_data_i[!DEgene_data_i$inner,] %>% select("id")

        n_inner <- nrow(DEgene_data_i.inner)
        n_outer <- nrow(DEgene_data_i.outer)

        # Angle range.
        area_angle.start <- area_angle.width * (i - 1L) + area_angle.margin
        area_angle.end   <- area_angle.width * i        - area_angle.margin


        if (n_inner > 0L) {
            angles <- seq(area_angle.start, area_angle.end,
                          length.out=(n_inner+2L))[2L:(n_inner+1L)]

            DEgene_data_i.inner <- DEgene_data_i.inner %>%
                               mutate(radius      = rep(DEgene.radius[1L], times=n_inner),
                                      angle       = angles,
                                      label.angle = sapply(angles, adjust.angle),
                                      label.hjust = sapply(angles, adjust.hjust, offset=DEgene_label.offset))
        } else {
            DEgene_data_i.inner <- NULL
        }

        if (n_outer > 0L) {
            angles <- seq(area_angle.start, area_angle.end, length.out=(n_outer+2L))[2L:(n_outer+1L)]

            DEgene_data_i.outer <- DEgene_data_i.outer %>%
                               mutate(radius      = rep(DEgene.radius[2L], times=n_outer),
                                      angle       = angles,
                                      label.angle = sapply(angles, adjust.angle),
                                      label.hjust = sapply(angles, adjust.hjust, offset=DEgene_label.offset))
        } else {
            DEgene_data_i.outer <- NULL
        }

        position_data_i <- rbind(DEgene_data_i.inner, DEgene_data_i.outer)


        position_data <- if (i == 1L) {
            position_data_i
        } else {
            rbind(position_data, position_data_i)
        }
    }

    # Add new columns.
    DEgene_data <- DEgene_data %>% left_join(position_data, by="id")

    return (DEgene_data)
}



#' Adjust a label angle.
#' 
#' 
adjust.angle <- function (angle) {
    if ((90 < angle) && (angle <= 270)) {
        return (angle - 180)
    } else {
        return (angle)
    }
}


#' Adjust a label hjust.
#' 
#' 
adjust.hjust <- function (angle, offset=0.3) {
    if ((90 < angle) && (angle <= 270)) {
        return (1 + offset)
    } else {
        return (0 - offset)
    }
}


#' Add new columns to ligand_data for main figures.
#' 
#' 
add.new_columns.ligand.main <- function (
    ligand_data,
    organs,
    gap.degree    = 1,
    ligand_label.offset = 0.3,
    ligand.radius = 13,
    ...
) {
    # Calculate positions for ligands.
    n_organs <- length(organs)

    area_angle.width  <- 360 / n_organs
    area_angle.margin <- gap.degree / 2

    position_data <- NULL

    for (i in seq_len(n_organs)) {
        ligand_data_i <- ligand_data[ligand_data$organ==organs[i],] %>% select("ligand")
        n_ligands     <- nrow(ligand_data_i)

        # Angle range.
        area_angle.start <- area_angle.width * (i - 1L) + area_angle.margin
        area_angle.end   <- area_angle.width * i        - area_angle.margin

        if (n_ligands > 0L) {
            angles <- seq(area_angle.start, area_angle.end, length.out=(n_ligands+2L))[2L:(n_ligands+1L)]

            position_data_i <- ligand_data_i %>%
                               mutate(radius      = rep(ligand.radius, times=n_ligands),
                                      angle       = angles,
                                      label.angle = sapply(angles, adjust.angle_2),
                                      label.vjust = sapply(angles, adjust.vjust, offset=ligand_label.offset))

            position_data <- rbind(position_data, position_data_i)
        }
    }

    # Add new columns.
    ligand_data <- ligand_data %>% left_join(position_data, by="ligand")

    return (ligand_data)
}


#' Add new columns to ligand_data for extended data figures.
#' 
#' 
add.new_columns.ligand.ext <- function (
    ligand_data,
    organs,
    gap.degree    = 1,
    ligand_label.offset = 0.3,
    ligand.radius = 13,
    ...
) {
    # Calculate positions for ligands.
    n_organs <- length(organs)

    area_angle.width  <- 360 / n_organs
    area_angle.margin <- gap.degree / 2

    position_data <- NULL

    for (i in seq_len(n_organs)) {
        ligand_data_i <- ligand_data[ligand_data$organ==organs[i],] %>% select("ligand")
        n_ligands     <- nrow(ligand_data_i)

        # Angle range.
        area_angle.start <- area_angle.width * (i - 1L) + area_angle.margin
        area_angle.end   <- area_angle.width * i        - area_angle.margin

        if (n_ligands > 0L) {
            angles <- seq(area_angle.start, area_angle.end, length.out=(n_ligands+2L))[2L:(n_ligands+1L)]

            position_data_i <- ligand_data_i %>%
                               mutate(radius      = rep(ligand.radius, times=n_ligands),
                                      angle       = angles,
                                      label.angle = sapply(angles, adjust.angle),
                                      label.hjust = sapply(angles, adjust.hjust, offset=ligand_label.offset))

            position_data <- rbind(position_data, position_data_i)
        }
    }

    # Add new columns.
    ligand_data <- ligand_data %>% left_join(position_data, by="ligand")

    return (ligand_data)
}


adjust.angle_2 <- function (angle) {
    if ((0 <= angle) && (angle <= 180)) {
        return (angle - 90)
    } else {
        return (angle + 90)
    }
}


adjust.vjust <- function (angle, offset=0.3) {
    if ((0 <= angle) && (angle <= 180)) {
        return (0 - offset)
    } else {
        return (1 + offset)
    }
}


# For main figure.
#
#
if (FALSE) {
    # Results of KEGG-mining.
    DEgene_data   <- read.xlsx("../Python/data/DEgene_data.xlsx")
    edge_data     <- read.xlsx("../Python/data/DEgene-edge_data.xlsx")
    upstream_data <- read.xlsx("../Python/data/upstream-edge_data.xlsx")

    # Drawing.

    nw <- generate.DEgeneNetwork.main(
              DEgene_data,
              edge_data,
              upstream_data,
              DEgene_condition = "F0-3 Adult",
              edge_condition   = "F0-3 Adult",
              organ_area.alpha = 0.35,
              organ_area.radius   = c(12, 21),
              organ_area.n_points = 30L,
              organ_label.color   = "#808080",
              organ_label.radius  = 18.5,
              organ_label.size    = 17,
              DEgene.colors       = c("F0-3 Adult"="#00B0F6"),
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
              edge.colors   = c("F0-3 Adult"="#00B0F6"),
              edge.width    = 1,
              edge.n_points = 50L,
              edge_legend   = "Edge condition",
              gap.degree = 2,
              arrow      = grid::arrow(angle=15, length=unit(8,"mm"), type="closed")
    )

    ggsave("./data/DEgene-network.pdf", plot=nw,
           limitsize=FALSE, units="mm", height=800, width=800)
}


# For extended data figures.
#
#
if (FALSE) {
    library(circlize)
    library(scales)

    # Results of KEGG-mining.
    DEgene_data   <- read.xlsx("../Python/data/DEgene_data.xlsx")
    edge_data     <- read.xlsx("../Python/data/DEgene-edge_data.xlsx")
    upstream_data <- read.xlsx("../Python/data/upstream-edge_data.xlsx")

    # Drawing

    nw <- generate.DEgeneNetwork.ext(
              DEgene_data,
              edge_data,
              upstream_data,
              DEgene_condition = "F0-3 Adult",
              edge_condition   = "F0-3 Adult",
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

    ggsave("./data/Background.pdf", plot=nw$background,
           limitsize=FALSE, units="mm", height=800, width=800)


    circos.par("clock.wise"=FALSE, gap.degree=2)

    pdf("./data/Edge.pdf")

    chordDiagram(
        nw$edge,
        directional  = TRUE,
        transparency = 0.25,
        row.col  = hue_pal()(nrow(result.ext$matrix)),
        grid.col = hue_pal()(nrow(result.ext$matrix)),
        scale    = TRUE,
        annotationTrack = c("name","grid"),
        annotationTrackHeight = c(mm_h(1), mm_h(4)),
        direction.type = c("diffHeight", "arrows"),
        link.arr.type  = "big.arrow",
        reduce = -1
    )

    dev.off()
}