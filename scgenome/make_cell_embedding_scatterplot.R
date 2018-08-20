#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gtools)
library(tidyr)


palette <- c(
    "#1d1d1d", "#ebce2b", "#702c8c", "#db6917", "#96cde6", "#ba1c30",
    "#c0bd7f", "#7f7e80", "#5fa641", "#d485b2", "#4277b6", "#df8461",
    "#463397", "#e1a11a", "#91218c", "#e8e948", "#7e1510", "#92ae31",
    "#6f340d", "#d32b1e", "#2b3514"
)

get_args <- function() {
    p <- ArgumentParser(description="Plot cell UMAP embedding")

    p$add_argument("embedding", help="UMAP output tsv file")
    p$add_argument("clusters", help="cell cluster assignments tsv file")
    p$add_argument("pdf", help="output plot pdf")

    return(p$parse_args())
}

open_file <- function(fn) {
    if (fn %in% c("-", "/dev/stdin")) {
        return(file("stdin", open="r"))
    } else if (grepl("^/dev/fd/", fn)) {
        return(fifo(fn, open="r"))
    } else {
        return(file(fn, open="r"))
    }
}

read_tsv <- function(fn, ...) {
    con <- open_file(fn)
    df <- read.delim(con, check.names=FALSE, stringsAsFactors=FALSE, ...)
    close(con)
    return(df)
}

annot_cluster_labels <- function(clusters) {
    for (cl in unique(clusters)) {
        mask <- clusters == cl

        n <- sum(mask)
        cl_label <- paste0(cl, " (", n, ")")
        clusters[mask] <- cl_label
    }

    return(clusters)
}

make_umap_embedding_plot <- function(embedding, clusters) {
    embedding <- full_join(embedding, clusters)

    noise <- subset(embedding, cluster == -1)
    embedding <- subset(embedding, cluster != -1)

    cl_labels <- annot_cluster_labels(embedding$cluster)
    embedding$cluster <- factor(cl_labels, levels=mixedsort(unique(cl_labels)))
    cluster_centres <- embedding %>%
        group_by(cluster) %>%
        summarise(centre_x=mean(x), centre_y=mean(y))

    ggplot(embedding, aes(x=x, y=y, colour=cluster)) +
        geom_point(alpha=0.5) +
        geom_point(data=noise, pch=3, colour="black", alpha=0.25) +
        geom_text_repel(
            data=cluster_centres, aes(x=centre_x, y=centre_y, label=cluster),
            box.padding=1, show.legend=FALSE
        ) +
        scale_colour_manual(values=palette) +
        guides(colour=FALSE) +
        theme_minimal() +
        theme(
            panel.grid.minor=element_blank(),
            axis.title=element_blank()
        )
}

main <- function() {
    argv <- get_args()

    embedding <- read_tsv(argv$embedding)
    clusters <- read_tsv(argv$clusters)
    make_umap_embedding_plot(embedding, clusters)

    ggsave(argv$pdf, width=6, height=5, useDingbats=FALSE)
}

main()
