#' Collection of functions to draw.
#' 
#' 
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(geomtextpath)
library(ggnewscale)
library(extrafont)


#' Draw organ areas.
#' 
#' 
draw.organ_areas <- function (
    nw,
    organs,
    organ_area.alpha    = 1,
    organ_area.colors   = c(),
    organ_area.n_points = 30L,
    organ_area.radius   = c(8, 10),
    organ_label.color   = "#000000",
    organ_label.radius  = 11,
    organ_label.size    = 2,
    gap.degree          = 1,
    ...
) {
    n_organs <- length(organs)

    # Angle width of each area.
    area_angle.width  <- 360 / n_organs
    area_angle.margin <- gap.degree / 2

    for (i in 1:n_organs) {
        organ <- organs[i]

        # Angle range of organs[i].
        area_angle.start <- area_angle.width * (i - 1L) + area_angle.margin
        area_angle.end   <- area_angle.width * i        - area_angle.margin

        inner_angles <- seq(area_angle.start, area_angle.end,
                            length.out = organ_area.n_points) * pi / 180
        outer_angles <- sort(inner_angles, decreasing=TRUE)

        # Data of points of organ area.
        points_i.area <- rbind(data.frame(x = organ_area.radius[1L] * cos(inner_angles),
                                          y = organ_area.radius[1L] * sin(inner_angles)),
                               data.frame(x = organ_area.radius[2L] * cos(outer_angles),
                                          y = organ_area.radius[2L] * sin(outer_angles))) %>%
                         mutate(group=i, fill=organ)

        points.area <- if (i == 1L) {
            points_i.area
        } else {
            rbind(points.area, points_i.area)
        }

        area_angle.mid <- (area_angle.start + area_angle.end) / 2

        angles <- if ((0 <= area_angle.mid) && (area_angle.mid <= 180)) {
            seq(area_angle.end, area_angle.start,
                length.out = organ_area.n_points) * pi / 180
        } else {
            seq(area_angle.start, area_angle.end,
                length.out = organ_area.n_points) * pi / 180
        }

        points_i.label <- data.frame(x = organ_label.radius * cos(angles),
                                     y = organ_label.radius * sin(angles)) %>%
                          mutate(group=i, label=organ)

        points.label <- if (i == 1L) {
            points_i.label
        } else {
            rbind(points.label, points_i.label)
        }
    }

    # Draw organ areas and labels.

    if (length(organ_area.colors) == 0L) {
        organ_area.colors <- rep(organ_label.color, times=n_organs)
        names(organ_area.colors) <- organs
    }

    nw <- nw +
          geom_polygon(
              data        = points.area,
              mapping     = aes(x=x, y=y, group=group, fill=fill),
              alpha       = organ_area.alpha,
              color       = NA,
              inherit.aes = FALSE,
              show.legend = FALSE) +
          scale_fill_manual(
              values = organ_area.colors,
              guide  = "none") +
          geom_textpath(
              data        = points.label,
              mapping     = aes(x=x, y=y, group=group, label=label),
              color       = organ_label.color,
              size        = organ_label.size,
              text_only   = TRUE,
              upright     = FALSE,
              inherit.aes = FALSE,
              show.legend = FALSE)

    return (nw)
}



#' Draw edges from ligands to DEgenes.
#' 
#' 
draw.upstream_edges <- function (
    nw,
    DEgene_data,
    ligand_data,
    upstream_data,
    upstream_edge.color    = "#808080",
    upstream_edge.width    = 0.1,
    upstream_edge.n_points = 30L,
    ligand.radius          = 13,
    organ_area.radius      = c(8, 10),
    ...
) {
    # Data of edges from ligands to DEgenes.
    upstream_edges <- inner_join(
                          select(upstream_data, c("ligand","id")),
                          select(ligand_data,   c("ligand","angle")),
                          by = "ligand") %>%
                      select("ligand.angle"="angle", everything()) %>%
                    inner_join(
                        select(DEgene_data, c("id","radius","angle","inner")),
                        by = "id") %>%
                    select("DEgene.radius"="radius", "DEgene.angle"="angle", everything())

    n_edges <- nrow(upstream_edges)

    if (n_edges > 0L) {

        for (i in 1:n_edges) {
            ligand.angle  <- upstream_edges[i,"ligand.angle"]
            DEgene.radius <- upstream_edges[i,"DEgene.radius"]
            DEgene.angle  <- upstream_edges[i,"DEgene.angle"]
            DEgene.inner  <- upstream_edges[i,"inner"]

            radius.mid <- if (DEgene.inner) {
                organ_area.radius[1L]
            } else {
                organ_area.radius[2L]
            }
            angle.limit <- acos(radius.mid / ligand.radius) * 180 / pi
            angle.diff  <- abs(ligand.angle - DEgene.angle)

            if (angle.diff < angle.limit) {
                ligand.angle <- ligand.angle * pi / 180
                DEgene.angle <- DEgene.angle * pi / 180
                xy    <- ligand.radius * c(cos(ligand.angle),
                                           sin(ligand.angle))
                xyend <- DEgene.radius * c(cos(DEgene.angle),
                                           sin(DEgene.angle))
                points_i <- data.frame(x = c(xy[1L], xyend[1L]),
                                       y = c(xy[2L], xyend[2L])) %>%
                            mutate(group=i)
            } else {
                n_points <- adjust.n_points(
                                n_max = upstream_edge.n_points,
                                angle = angle.diff,
                                ...)

                points_i <- compute.circular_points(
                                radius      = c(ligand.radius, DEgene.radius),
                                angle_range = c(ligand.angle,  DEgene.angle),
                                n_points    = n_points) %>%
                            mutate(group=i)
            }

            points <- if (i == 1L) {
                points_i
            } else {
                rbind(points, points_i)
            }
        }

        # Draw edges from ligands to DEgenes.

        nw <- nw +
              new_scale_color() +
              geom_path(
                  data        = points,
                  mapping     = aes(x=x, y=y, group=group),
                  color       = upstream_edge.color,
                  linewidth   = upstream_edge.width,
                  inherit.aes = FALSE,
                  show.legend = FALSE)
    }

    return (nw)
}



#' Adjust the number of points.
#' 
#' 
adjust.n_points <- function (n_max, angle, n_min=7L, ...) {
    angle <- abs(angle)

    if (angle > 180) {
        angle <- angle %% 180L
    }

    return (max(n_min, round(angle * n_max / 180)))
}



#' Compute points on a circle.
#' 
#' 
compute.circular_points <- function (
    radius,
    angle_range = c(0, 360),
    center      = c(0, 0),
    n_points    = 50L
) {
    # The list of radius.
    if (length(radius) > 1L) {
        radius <- seq(radius[1L], radius[2L], length.out=n_points)
    }

    angles <- seq(angle_range[1L], angle_range[2L], length.out=n_points) * pi / 180

    points <- data.frame(x = center[1L] + radius * cos(angles),
                         y = center[2L] + radius * sin(angles))

    return (points)
}


#' Draw edges between DEgenes.
#' 
#' 
draw.DEgene_edges <- function (
    nw,
    DEgene_data,
    edge_data,
    arrow = grid::arrow(angle  = 15,
                        length = unit(1, "mm"),
                        type   = "closed"),
    arrow.head_ratio = 0.75,
    DEgene.radius    = c(8.1, 9.9),
    depth_ratio      = 0.75,
    edge.colors      = c(),
    edge.width       = 0.1,
    edge.n_points    = 30L,
    ...
) {
    # Data of edges.

    edge_data.inter <- edge_data[edge_data$from.organ!=edge_data$to.organ,]
    edge_data.intra <- edge_data[edge_data$from.organ==edge_data$to.organ,]

    n_inter <- nrow(edge_data.inter)
    n_intra <- nrow(edge_data.intra)

    if ((n_inter > 0L) || (n_intra > 0L)) {
        if (length(edge.colors) == 0L) {
            edge.colors        <- rep("#202020", times=length(unique(edge_data$color)))
            names(edge.colors) <- unique(edge_data$color)
        }

        nw <- nw + new_scale_color()

        # Draw inter edges.

        if (n_inter > 0L) {
            edge_data.inter <- left_join(
                                   edge_data.inter,
                                   select(DEgene_data, c("id","radius","angle")),
                                   by = c("from"="id")) %>%
                               select("from.angle"="angle", everything()) %>%
                               left_join(
                                   select(DEgene_data, c("id","angle")),
                                   by = c("to"="id")) %>%
                               select("to.angle"="angle", everything())

            nw <- draw.inter_edges(
                      nw,
                      edge_data.inter,
                      arrow            = arrow,
                      arrow.head_ratio = arrow.head_ratio,
                      edge.width       = edge.width,
                      edge.n_points    = edge.n_points,
                      ...)
        }

        if (n_intra > 0L) {
            edge_data.intra <- left_join(
                                   edge_data.intra,
                                   select(DEgene_data, c("id","radius","angle")),
                                   by = c("from"="id")) %>%
                               select("from.radius"="radius", "from.angle"="angle", everything()) %>%
                               left_join(
                                   select(DEgene_data, c("id","radius","angle")),
                                   by = c("to"="id")) %>%
                               select("to.radius"="radius", "to.angle"="angle", everything())

            nw <- draw.intra_edges(
                nw,
                edge_data.intra,
                arrow            = arrow,
                arrow.head_ratio = arrow.head_ratio,
                depth_ratio      = depth_ratio,
                edge.width       = edge.width,
                DEgene.radius    = DEgene.radius,
                edge.n_points    = edge.n_points,
                ...)
        }

        nw <- nw + scale_color_manual(
                       values = edge.colors,
                       guide  = "none")
    }

    return (nw)
}



#' Draw inter edges.
#' 
#' 
draw.inter_edges <- function (
    nw,
    edge_data.inter,
    arrow = grid::arrow(),
    arrow.head_ratio = 0.75,
    edge.width       = 0.1,
    edge.n_points    = 30L,
    ...
) {
    n_edges <- nrow(edge_data.inter)

    for (i in 1:n_edges) {
        radius      <- edge_data.inter[i,"radius"]
        angle.start <- edge_data.inter[i,"from.angle"]
        angle.end   <- edge_data.inter[i,"to.angle"]
        edge.color  <- edge_data.inter[i,"color"]

        points_i <- compute.inter_edge_points(
                        radius,
                        angle.start,
                        angle.end,
                        arrow.head_ratio = arrow.head_ratio,
                        edge.n_points    = edge.n_points,
                        ...)

        points_i.edge <- points_i$edge %>%
                         mutate(group=i, color=edge.color)

        points_i.head <- points_i$head %>%
                         mutate(group=i, color=edge.color)

        points.edge <- if (i == 1L) {
            points_i.edge
        } else {
            rbind(points.edge, points_i.edge)
        }

        points.head <- if (i == 1L) {
            points_i.head
        } else {
            rbind(points.head, points_i.head)
        }
    }

    nw <- nw +
          geom_path(
              data        = points.edge,
              mapping     = aes(x=x, y=y, group=group, color=color),
              linewidth   = edge.width,
              inherit.aes = FALSE,
              show.legend = FALSE) +
          geom_path(
              data        = points.head,
              mapping     = aes(x=x, y=y, group=group, color=color),
              arrow       = arrow,
              linewidth   = edge.width,
              inherit.aes = FALSE,
              show.legend = FALSE)

    return (nw)
}



#' Compute points of inter edges.
#' 
#' 
compute.inter_edge_points <- function (
    radius,
    angle.start,
    angle.end,
    arrow.head_ratio = 0.75,
    edge.n_points    = 100L,
    ...
) {
    direction  <- sign(angle.end - angle.start)

    angle.diff <- abs(angle.end - angle.start)
    angle.min  <- min(angle.start, angle.end)

    if (abs(angle.diff - 180) < 1e-10) {
        angle.start <- angle.start * pi / 180
        angle.end   <- angle.end   * pi / 180

        xy    <- radius * c(cos(angle.start), sin(angle.start))
        xyend <- radius * c(cos(angle.end),   sin(angle.end))

        xymid <- if (arrow.head_ratio <= 0) {
            xy + (xyend - xy) * 0.01
        } else if (arrow.head_ratio >= 1) {
            xyend
        } else {
            xy + (xyend - xy) * arrow.head_ratio
        }

        points.edge <- data.frame(x=c(xymid[1L], xyend[1L]), y=c(xymid[2L], xyend[2L]))
        points.head <- data.frame(x=c(xy[1L],    xymid[1L]), y=c(xy[2L],    xymid[2L]))
    } else {
        edge_radius      <- radius * abs(tan(angle.diff * pi / 360))
        edge_angle.width <- abs(angle.diff - 180)

        if (angle.diff < 180) {
            if (direction > 0) {
                edge_angle.start <- angle.min - 90
                edge_angle.end   <- angle.min - 90 - edge_angle.width
            } else {
                edge_angle.start <- angle.min - 90 - edge_angle.width
                edge_angle.end   <- angle.min - 90
            }
            edge_center.angle <- (angle.min + angle.diff / 2) * pi / 180
        } else {
            if (direction > 0) {
                edge_angle.start <- angle.min + 90
                edge_angle.end   <- angle.min + 90 + edge_angle.width
            } else {
                edge_angle.start <- angle.min + 90 + edge_angle.width
                edge_angle.end   <- angle.min + 90
            }
            edge_center.angle <- (angle.min + angle.diff / 2 - 180) * pi / 180
        }

        edge_center.radius <- radius / abs(cos(angle.diff * pi / 360))
        edge_center        <- c(edge_center.radius * cos(edge_center.angle),
                                edge_center.radius * sin(edge_center.angle))

        n_points <- adjust.n_points(
                        n_max = edge.n_points,
                        angle = edge_angle.start - edge_angle.end,
                        ...)

        points.edge <- compute.circular_points(
                           radius      = edge_radius,
                           angle_range = c(edge_angle.start, edge_angle.end),
                           center      = edge_center,
                           n_points    = n_points)

        n_head <- round(n_points * arrow.head_ratio)

        points.head <- if (n_head >= n_points) {
            points.edge[(n_points-1L):n_points,]
        } else if (n_head < 1L) {
            points.edge[1L:2L,]
        } else {
            points.edge[n_head:(n_head+1L),]
        }
    }

    return (list(edge=points.edge, head=points.head))
}



#' Draw intra edges.
#' 
#' 
draw.intra_edges <- function (
    nw,
    edge_data.intra,
    arrow            = grid::arrow(),
    arrow.head_ratio = 0.75,
    depth_ratio      = 0.75,
    edge.width       = 0.1,
    edge.n_points    = 30L,
    DEgene.radius    = c(8.1, 9.9),
    ...
) {
    n_edges    <- nrow(edge_data.intra)
    radius.mid <- mean(DEgene.radius)

    for (i in seq_len(n_edges)) {
        radius      <- edge_data.intra[i,"from.radius"]
        angle.start <- edge_data.intra[i,"from.angle"]
        radius.end  <- edge_data.intra[i,"to.radius"]
        angle.end   <- edge_data.intra[i,"to.angle"]
        edge.color  <- edge_data.intra[i,"color"]

        if (radius != radius.end) {
            n_points <- adjust.n_points(
                            n_max = edge.n_points,
                            angle = angle.start - angle.end,
                            ...)

            points_i.edge <- compute.circular_points(
                                 radius      = c(radius,      radius.end),
                                 angle_range = c(angle.start, angle.end),
                                 n_points    = n_points)
        } else {
            angle.diff <- abs(angle.end - angle.start)
            angle.mid <- (angle.start + angle.end) / 2

            rect.width <- radius.mid * angle.diff * pi / 180

            if (radius == DEgene.radius[1L]) {
                direction <- sign(angle.end - angle.start)
                position  <- -1L
            } else {
                direction <- sign(angle.start - angle.end)
                position  <- 1L
            }

            arc_points <- compute.arc_points_in_rectangle(
                              radius.inner  = DEgene.radius[1L],
                              radius.outer  = DEgene.radius[2L],
                              rect.width    = rect.width,
                              depth_ratio   = depth_ratio,
                              direction     = direction,
                              edge.n_points = edge.n_points,
                              position      = position,
                              ...)

            radii  <- arc_points$y
            angles <- angle.mid * pi / 180 - arc_points$x / radius.mid

            points_i.edge <- data.frame(x = radii * cos(angles),
                                        y = radii * sin(angles))
        }

        points_i.edge <- points_i.edge %>% mutate(group=i, color=edge.color)

        n_points <- nrow(points_i.edge)
        n_head   <- round(n_points * arrow.head_ratio)

        points_i.head <- if (n_head >= n_points) {
            points_i.edge[(n_points-1L):n_points,]
        } else if (n_head < 1L) {
            points_i.edge[1L:2L,]
        } else {
            points_i.edge[n_head:(n_head+1L),]
        }

        points.edge <- if (i == 1L) {
            points_i.edge
        } else {
            rbind(points.edge, points_i.edge)
        }

        points.head <- if (i == 1L) {
            points_i.head
        } else {
            rbind(points.head, points_i.head)
        }
    }

    nw <- nw +
          geom_path(
              data        = points.edge,
              mapping     = aes(x=x, y=y, group=group, color=color),
              linewidth   = edge.width,
              inherit.aes = FALSE,
              show.legend = FALSE) +
          geom_path(
              data        = points.head,
              mapping     = aes(x=x, y=y, group=group, color=color),
              arrow       = arrow,
              linewidth   = edge.width,
              inherit.aes = FALSE,
              show.legend = FALSE)

    return (nw)
}



#' Compute points on arc in a rectangle.
#' 
#' 
compute.arc_points_in_rectangle <- function (
    radius.inner,
    radius.outer,
    rect.width,
    depth_ratio   = 0.5,
    direction     = 1L,
    edge.n_points = 30L,
    position      = -1L,
    ...
) {
    rect.height <- radius.outer - radius.inner

    if (2 * rect.height * depth_ratio >= rect.width) {
        arc_radius      <- rect.width / 2
        arc_angle.start <- -90 * direction - 90 * position
        arc_angle.end   <-  90 * direction - 90 * position

        arc_center <- if (position == -1L) {
            c(0, radius.inner)
        } else {
            c(0, radius.outer)
        }
    } else {
        depth <- rect.height * depth_ratio

        arc_radius      <- (rect.width ** 2 + 4 * (depth ** 2)) / (8 * depth)
        arc_angle       <- asin(rect.width / (2 * arc_radius)) * 180 / pi
        arc_angle.start <- -arc_angle * direction - 90 * position
        arc_angle.end   <-  arc_angle * direction - 90 * position

        arc_center <- if (position == -1L) {
            c(0, radius.inner + depth - arc_radius)
        } else {
            c(0, radius.outer - depth + arc_radius)
        }
    }

    n_points <- adjust.n_points(
                    n_max = edge.n_points,
                    angle = arc_angle.start - arc_angle.end,
                    ...)
    angles <- seq(arc_angle.start, arc_angle.end, length.out=n_points) * pi / 180

    points <- data.frame(x = arc_center[1L] + arc_radius * cos(angles),
                         y = arc_center[2L] + arc_radius * sin(angles))

    return(points)
}



#' Draw DEgenes as nodes.
#' 
#' 
draw.DEgenes <- function (
    nw,
    DEgene_data,
    DEgene.colors     = c(),
    DEgene.size       = 0.5,
    DEgene_label.size = 0.8,
    DEgene_legend     = "Condition",
    ...
) {
    angles <- DEgene_data$angle * pi / 180

    points <- DEgene_data %>%
              mutate(x = DEgene_data$radius * cos(angles),
                     y = DEgene_data$radius * sin(angles))

    if (length(DEgene.colors) == 0L) {
        DEgene.colors <- rep("#202020", times=length(unique(points$color)))
        names(DEgene.colors) <- unique(points$color)
    }

    n_colors <- length(DEgene.colors)

    nw <- nw +
        new_scale_color() +
        geom_point(
            data        = select(points, c("x","y","color")),
            mapping     = aes(x=x, y=y, color=color),
            size        = DEgene.size,
            inherit.aes = FALSE,
            show.legend = TRUE) +
        scale_color_manual(name=DEgene_legend, values=DEgene.colors) +
        guides(
            color = guide_legend(nrow  = n_colors,
                                 byrow = TRUE, 
                                 override.aes = list(shape    = rep(16,      times=n_colors),
                                                     linetype = rep("blank", times=n_colors)))) +
        geom_text(
            data        = select(points, c("x","y","label","label.angle","label.hjust")),
            mapping     = aes(x=x, y=y, label=label, angle=label.angle, hjust=label.hjust),
            size        = DEgene_label.size,
            vjust       = 0.5,
            inherit.aes = FALSE,
            show.legend = FALSE)

    return (nw)
}



#' Draw ligands as nodes for MainFigs.
#' 
#' 
draw.ligands.main <- function (
    nw,
    ligand_data,
    ligand.color      = "#202020",
    ligand.shape      = 17,
    ligand.size       = 0.5,
    ligand_label.size = 0.8,
    ...
) {
    angles <- ligand_data$angle * pi / 180

    points <- ligand_data %>%
              mutate(x = ligand_data$radius * cos(angles),
                     y = ligand_data$radius * sin(angles))

    nw <- nw +
          new_scale_color() +
          geom_point(
              data        = select(points, c("x","y")),
              mapping     = aes(x=x, y=y),
              color       = ligand.color,
              shape       = ligand.shape,
              size        = ligand.size,
              inherit.aes = FALSE,
              show.legend = FALSE) +
          geom_text(
              data        = select(points, c("x","y","label","label.angle","label.vjust")),
              mapping     = aes(x=x, y=y, label=label, angle=label.angle, vjust=label.vjust),
              size        = ligand_label.size,
              hjust       = 0.5,
              inherit.aes = FALSE,
              show.legend = FALSE)

    return (nw)
}



#' Draw ligands as nodes for ExtDataFigs.
#' 
#' 
draw.ligands.ext <- function (
    nw,
    ligand_data,
    ligand.color        = "#202020",
    ligand.shape        = 17,
    ligand.size         = 0.5,
    ligand_label.size   = 0.8,
    ...
) {
    angles <- ligand_data$angle * pi / 180

    points <- ligand_data %>%
              mutate(x = ligand_data$radius * cos(angles),
                     y = ligand_data$radius * sin(angles))

    nw <- nw +
          new_scale_color() +
          geom_point(
              data        = select(points, c("x","y")),
              mapping     = aes(x=x, y=y),
              color       = ligand.color,
              shape       = ligand.shape,
              size        = ligand.size,
              inherit.aes = FALSE,
              show.legend = FALSE) +
          geom_text(
              data        = select(points, c("x","y","label","label.angle","label.hjust")),
              mapping     = aes(x=x, y=y, label=label, angle=label.angle, hjust=label.hjust),
              size        = ligand_label.size,
              vjust       = 0.5,
              inherit.aes = FALSE,
              show.legend = FALSE)

    return (nw)
}



#' Draw a legend of edges.
#' 
#' 
draw.edge_legend <- function (
    nw,
    edge_data,
    edge.colors = c(),
    edge_legend = "Condition",
    ...
) {
    if ((nrow(edge_data) > 0L) || (length(edge.colors) > 0L)) {
        if (length(edge.colors) > 0L) {
            edge_color_data <- data.frame(color=rep(names(edge.colors), each=2L)) %>%
                               mutate(x=0, y=0)
        } else {
            edge_color_data <- data.frame(color=rep(unique(edge_data$color), each=2L)) %>%
                               mutate(x=0, y=0)

            edge.colors <- rep("#202020", times=length(unique(edge_data$color)))
            names(edge.colors) <- unique(edge_data$color)
        }

        n_colors <- length(edge.colors)

        nw <- nw +
              new_scale_color() +
              geom_path(
                  data        = edge_color_data,
                  mapping     = aes(x=x, y=y, color=color, group=color),
                  size        = 0,
                  inherit.aes = FALSE,
                  show.legend = TRUE) +
              scale_color_manual(name=edge_legend, values=edge.colors) +
              guides(color = guide_legend(nrow  = n_colors,
                                          byrow = TRUE, 
                                          override.aes = list(shape    = rep(NA,      times=n_colors),
                                                              linetype = rep("solid", times=n_colors))))
    }

    return (nw)
}