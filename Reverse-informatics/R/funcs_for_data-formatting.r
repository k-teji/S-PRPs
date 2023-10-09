#' Collection of functions to format data.
#' 
#' 
library(tidyverse)

#' Generate data for MainFigs.
#' 
#' 
generate.draw_data.main <- function (
    DEgene_data,
    edge_data,
    upstream_data,
    ...
) {
    DEgene_data.new <- generate.DEgene_data(DEgene_data, ...)

    edge_data.new <- generate.edge_data(DEgene_data.new, edge_data, ...)

    result <- generate.upstream_data.main(DEgene_data.new, upstream_data, ...)

    return (list(DEgene   = DEgene_data.new,
                 edge     = edge_data.new,
                 ligand   = result$ligand,
                 upstream = result$upstream))
}



#' Generate data for ExtDataFigs.
#' 
#' 
generate.draw_data.ext <- function (
    DEgene_data,
    edge_data,
    upstream_data,
    ...
) {
    DEgene_data.new <- generate.DEgene_data(DEgene_data, ...)

    organs   <- sort(unique(DEgene_data.new$organ))
    n_organs <- length(organs)

    edge_data.new <- generate.edge_data(DEgene_data.new, edge_data, ...)


    # Generate a matrix for edges between organs.
    #   edge_matrix[organ_i, organ_j] is
    #   the number of edges from organ_i to organ_j.

    edge_matrix <- diag(0, n_organs)
    rownames(edge_matrix) <- organs
    colnames(edge_matrix) <- organs

    for (i in 1:n_organs) {
        organ_i <- organs[i]
        for (j in 1:n_organs) {
            organ_j <- organs[j]
            n_edges <- nrow(edge_data.new[(edge_data.new$from.organ==organ_i) &
                                          (edge_data.new$to.organ==organ_j),])
            edge_matrix[organ_i, organ_j] <- n_edges
        }
    }


    edge_data.new <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0),]
    colnames(edge_data.new) <- c("from.organ","from","to.organ","to","color")

    result <- generate.upstream_data.ext(DEgene_data.new, upstream_data, ...)

    return (list(DEgene   = result$DEgene,
                 edge     = edge_data.new,
                 ligand   = result$ligand,
                 upstream = result$upstream,
                 matrix   = edge_matrix))
}



#' Generate data of DEgenes as node for drawing.
#' 
#' 
generate.DEgene_data <- function (
    DEgene_data,
    DEgene_condition = "DEgene",
    ...
) {
    data_columns <- colnames(DEgene_data)

    DEgene_data.new <- if ("condition" %in% data_columns) {
        DEgene_data %>%
        select(c("organ","DEgene","condition"))
    } else {
        DEgene_data %>%
        select(c("organ","DEgene")) %>%
        mutate(condition=DEgene_condition)
    }

    DEgene_data.new <- DEgene_data.new[order(DEgene_data.new$organ,
                                             DEgene_data.new$condition,
                                             DEgene_data.new$DEgene),] %>%
                       select("label"="DEgene", "color"="condition", everything()) %>%
                       mutate(id=1:nrow(DEgene_data.new))

    # Remove unnecessary IDs.
    empty_organs <- unique(DEgene_data.new[is.na(DEgene_data.new$label),]$organ)

    if (length(empty_organs) > 0L) {
        DEgene_data.new <- DEgene_data.new[!is.na(DEgene_data.new$label),]

        DEgene_data.add <- data.frame(organ=empty_organs) %>%
                           mutate(id=NA, label="", color="")

        DEgene_data.new <- rbind(DEgene_data.new, DEgene_data.add)

        DEgene_data.new <- DEgene_data.new[order(DEgene_data.new$organ,
                                                 DEgene_data.new$color,
                                                 DEgene_data.new$label),]
    }

    return (DEgene_data.new)
}



#' Generate data of edges between DEgenes for drawing.
#' 
#' 
generate.edge_data <- function(
    DEgene_data.new,
    edge_data,
    edge_condition = "Edge",
    ...
) {
    data_columns <- colnames(edge_data)

    edge_data.new <- if ("condition" %in% data_columns) {
        edge_data %>%
        select(c("from.organ","from.DEgene",
                 "to.organ","to.DEgene","condition")) %>%
        distinct()
    } else {
        edge_data %>%
        select(c("from.organ","from.DEgene",
                 "to.organ","to.DEgene")) %>%
        distinct() %>%
        mutate(condition=edge_condition)
    }

    edge_data.new <- inner_join(
                         edge_data.new,
                         select(DEgene_data.new, c("organ","id","label")),
                         by = c("from.organ"="organ","from.DEgene"="label")) %>%
                     select("from"="id", everything()) %>%
                     inner_join(
                         select(DEgene_data.new, c("organ","id","label")),
                         by = c("to.organ"="organ","to.DEgene"="label")) %>%
                     select("to"="id", "color"="condition", everything())

    edge_data.new <- edge_data.new %>%
                     select(c("from.organ","from","to.organ","to","color")) %>%
                     distinct()

    # Remove self-loops.
    edge_data.new <- edge_data.new[edge_data.new$from!=edge_data.new$to,]

    return (edge_data.new)
}



#' Generate data of edges and data of ligands as nodes for MainFigs.
#' 
#' 
generate.upstream_data.main <- function (
    DEgene_data.new,
    upstream_data,
    topn_ligands = 3L,
    ...
) {
    toDEgene_data <- upstream_data %>%
                     select(c("organ","DEgene")) %>%
                     distinct()

    for (i in 1:nrow(toDEgene_data)) {
        organ  <- toDEgene_data[i,"organ"]
        DEgene <- toDEgene_data[i,"DEgene"]

        ligands_i <- upstream_data[(upstream_data$organ==organ) &
                                   (upstream_data$DEgene==DEgene),]$ligand %>%
                     unique()

        if (length(ligands_i) > topn_ligands) {
            ligands_i <- c(ligands_i[1L:(topn_ligands - 1L)], "etc.")
        }

        ligands_i <- paste(ligands_i, collapse="\n")

        upstream_data.new <- if (i == 1L) {
            data.frame(ligand=ligands_i, organ=organ, DEgene=DEgene)
        } else {
            rbind(upstream_data.new,
                  data.frame(ligand=ligands_i, organ=organ, DEgene=DEgene))
        }
    }

    ligand_data.new <- upstream_data.new %>%
                       select(c("organ","ligand")) %>%
                       distinct()

    ligand_data.new <- ligand_data.new[order(ligand_data.new$organ,
                                             ligand_data.new$ligand),] %>%
                       select("label"="ligand", everything()) %>%
                       mutate(ligand=1:nrow(ligand_data.new))

    upstream_data.new <- inner_join(
                             upstream_data.new,
                             select(DEgene_data.new, c("organ","id","label")),
                             by = c("organ"="organ","DEgene"="label")) %>%
                         select(c("organ","ligand","id")) %>%
                         select("label"="ligand", everything()) %>%
                         inner_join(
                             ligand_data.new,
                             by = c("organ","label")) %>%
                         select(-c("label"))

    return (list(ligand   = ligand_data.new,
                 upstream = upstream_data.new))
}


#' Generate data of edges and data of ligands as nodes for ExtDataFigs.
#' 
#' 
generate.upstream_data.ext <- function (DEgene_data.new, upstream_data, n_ligands=10L, ...) {
    ligand_data <- upstream_data %>%
                   select(c("organ","ligand")) %>%
                   distinct() %>%
                   select("label"="ligand", everything())

    organs <- unique(ligand_data$organ)

    for (i in 1:length(organs)) {
        organ <- organs[i]

        ligand_data_i <- ligand_data[ligand_data$organ==organ,]

        if (nrow(ligand_data_i) > n_ligands) {
            ligand_data_i <- rbind(ligand_data_i[1L:(n_ligands-1L),],
                                   data.frame(organ=organ, label="etc."))
        }

        ligand_data.new <- if (i == 1L) {
            ligand_data_i
        } else {
            rbind(ligand_data.new, ligand_data_i)
        }
    }

    ligand_data.new <- ligand_data.new[order(ligand_data.new$organ,
                                             ligand_data.new$label),] %>%
                       mutate(ligand=1:nrow(ligand_data.new))

    organs_out <- setdiff(unique(DEgene_data.new$organ), organs)

    DEgene_data.new <- ligand_data.new %>%
                       select("id"="ligand", everything()) %>%
                       select(c("organ","id")) %>%
                       mutate(label="", color="none")

    if (length(organs_out) > 0L) {
        DEgene_data.add <- data.frame(organ=organs_out) %>%
                           mutate(id=NA, label="", color="none")
    }

    upstream_data.new <- ligand_data.new %>%
                         select(c("organ","ligand")) %>%
                         mutate(id=ligand_data.new$ligand)


    return (list(DEgene   = DEgene_data.new,
                 ligand   = ligand_data.new,
                 upstream = upstream_data.new))
}