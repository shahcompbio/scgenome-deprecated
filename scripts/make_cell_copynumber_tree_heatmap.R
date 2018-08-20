#!/usr/bin/env Rscript

library(ape)
library(argparse)
library(ComplexHeatmap)
library(dplyr)
library(gtools)
library(RColorBrewer)

options(expressions=50000)

cluster_palette <- c(
    "#be5f72", "#d74058", "#dc4229", "#a6552c", "#df956f", "#e47a33",
    "#d49f34", "#836e2c", "#b2ad5a", "#92b539", "#4c7d38", "#4dc041",
    "#5dba7f", "#47b8c3", "#6280ca", "#7b57db", "#ce8bd1", "#934f94",
    "#cb48cb", "#d74391"
)


get_args <- function() {
    p <- ArgumentParser(description="Plot cell copynumber heatmap with tree")

    p$add_argument("tree", help="cell newick tree file")
    p$add_argument("copynumber", help="cell copynumber tsv file")
    p$add_argument("clusters", help="cell clusters tsv file")
    p$add_argument("pdf", help="output plot pdf")

    p$add_argument("--samples", help="cell samples tsv file")

    return(p$parse_args())
}

read_tsv <- function(fn, ...) {
    df <- read.delim(fn, check.names=FALSE, stringsAsFactors=FALSE, ...)
    close(con)
    return(df)
}

make_discrete_palette <- function(pal_name, levels) {
    pal <- brewer.pal(max(length(levels), 3), pal_name)
    names(pal) <- levels
    pal <- pal[levels]
    return(pal)
}

make_cluster_palette <- function(levels) {
    pal <- cluster_palette
    names(pal) <- levels
    pal <- pal[levels]
    return(pal)
}

format_copynumber_values <- function(copynumber) {
    copynumber[copynumber > 6] <- 6
    for(col in colnames(copynumber)) {
        values <- as.character(copynumber[, col])
        values[values == "6"] <- ">=6"
        copynumber[, col] <- values
    }
    return(copynumber)
}

space_copynumber_columns <- function(copynumber, spacer_cols) {
    chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
    spacer <- as.data.frame(matrix(
        data=NA, nrow=nrow(copynumber), ncol=spacer_cols
    ))
    chrom_copynumber_dfs <- list()
    for(chrom in mixedsort(unique(chroms))) {
        chrom_copynumber <- copynumber[, chroms == chrom, drop=FALSE]
        chrom_copynumber_dfs <- c(chrom_copynumber_dfs, chrom_copynumber)
        chrom_copynumber_dfs <- c(chrom_copynumber_dfs, spacer)
    }
    chrom_copynumber_dfs[length(chrom_copynumber_dfs)] <- NULL
    copynumber <- do.call(cbind, chrom_copynumber_dfs)

    return(copynumber)
}

format_copynumber <- function(copynumber, tree, spacer_cols=40) {
    rownames(copynumber) <- paste0(
        copynumber$chr, ":", copynumber$start, ":", copynumber$end
    )
    copynumber <- subset(copynumber, select=-c(chr, start, end, width))
    copynumber <- as.data.frame(t(copynumber))
    copynumber <- copynumber[tree$tip, ]

    copynumber <- format_copynumber_values(copynumber)
    copynumber <- space_copynumber_columns(copynumber, spacer_cols)
    return(copynumber)
}

format_clusters <- function(clusters, tree) {
    cluster_counts <- clusters %>% group_by(cluster) %>% summarise(count=n())
    clusters$cluster_label
    for(i in 1:nrow(cluster_counts)) {
        cluster <- unlist(cluster_counts[i, "cluster"], use.names=FALSE)
        cluster_count <- unlist(cluster_counts[i, "count"], use.names=FALSE)
        cluster_label <- paste0(cluster, " (", cluster_count, ")")
        clusters[clusters$cluster == cluster, "cluster_label"] <- cluster_label
    }
    rownames(clusters) <- clusters$cell
    clusters <- clusters[tree$tip, ]

    return(clusters)
}

get_chrom_label_pos <- function(copynumber) {
    chrom_label_pos <- list()
    chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
    uniq_chroms <- c(as.character(1:22), "X", "Y")
    for(chrom in uniq_chroms) {
        chrom_idx <- which(chroms == chrom)
        chrom_label_pos[[chrom]] <- as.integer(round(mean(chrom_idx)))
    }
    return(chrom_label_pos)
}

get_cluster_label_pos <- function(clusters) {
    cluster_label_pos <- list()
    for(cluster in unique(clusters$cluster)) {
        cluster_idx <- which(clusters$cluster == cluster)
        cluster_label_pos[[cluster]] <- as.integer(round(mean(cluster_idx)))
    }
    return(cluster_label_pos)
}

get_row_annot <- function(tree, clusters, samples) {
    cluster_levels <- mixedsort(unique(clusters$cluster_label))
    annot_colours <- list(Clone=make_cluster_palette(cluster_levels))

    rowname_generator <- function(index) {
        cluster_label_pos <- get_cluster_label_pos(clusters)
        y_pos <- 1 - unlist(cluster_label_pos) / nrow(clusters)
        grid.text(
            names(cluster_label_pos), 0.5, y_pos,
            just=c("centre", "centre")
        )
    }

    if(!is.null(samples)) {
        rownames(samples) <- samples$cell
        samples <- samples[tree$tip, ]

        sample_levels <- mixedsort(unique(samples$sample))
        annot_colours$Sample <- make_discrete_palette("Set2", sample_levels)

        row_annot <- HeatmapAnnotation(
            col=annot_colours, annotation_width=0.5, which="row",
            show_annotation_name=TRUE,
            Clone=clusters$cluster_label, rowname=rowname_generator,
            Sample=samples$sample
        )
    } else {
        row_annot <- HeatmapAnnotation(
            col=annot_colours, annotation_width=0.5, which="row",
            show_annotation_name=TRUE,
            Clone=clusters$cluster_label, rowname=rowname_generator,
        )
    }

    return(row_annot)
}

make_cell_copynumber_tree_heatmap <- function(tree, copynumber, clusters,
                                              samples) {
    copynumber <- format_copynumber(copynumber, tree)
    clusters <- format_clusters(clusters, tree)

    cn_colours <- structure(
        c(
            "#2166ac", "#92c5de", "#bababa", "#fddbc7", "#f4a582", "#d6604d",
            "#b2182b"
        ),
        names=c("0", "1", "2", "3", "4", "5", ">=6")
    )

    bottom_annot <- HeatmapAnnotation(column_labels=function(index) {
        chrom_label_pos <- get_chrom_label_pos(copynumber)
        x_pos <- unlist(chrom_label_pos)/ncol(copynumber)
        grid.text(names(chrom_label_pos), x_pos, 1, just=c("center", "top"))
    })

    row_annot <- get_row_annot(tree, clusters, samples)
    h <- row_annot + Heatmap(
        name="Copy Number",
        copynumber,
        col=cn_colours,
        na_col="white",
        cluster_rows=as.hclust(tree),
        row_dend_width=unit(40, "mm"),
        show_row_names=FALSE,
        cluster_columns=FALSE,
        show_column_names=FALSE,
        bottom_annotation=bottom_annot,
        use_raster=TRUE,
        raster_quality=5
    )
    draw(h, row_dend_side="left")
}

main <- function() {
    argv <- get_args()

    tree <- read.tree(argv$tree)
    copynumber <- read_tsv(argv$copynumber)
    clusters <- read_tsv(argv$clusters)

    samples <- NULL
    if(!is.null(argv$samples)) {
        samples <- read_tsv(argv$samples)
    }

    pdf(argv$pdf, width=14)
    make_cell_copynumber_tree_heatmap(tree, copynumber, clusters, samples)
    dev.off()
}

main()
