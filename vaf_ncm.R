#!/bin/bash
# -*- mode:R -*-

'\' >/dev/null 2>&1 || true
# This is bash code to set up the environment

if ! which Rscript &>/dev/null; then
    module load R
fi

# Bash setup code ends here
Rscript - "$@" <<"EOF";
invisible('#')
## R code starts after this line

suppressPackageStartupMessages({
    library(conflicted)
    library(argparse)
    library(fs)
    library(readr)
    library(dplyr)
    conflict_prefer("filter", "dplyr", quiet = TRUE)
    library(magrittr)
    library(assertthat)
    library(stringr)
    library(tidyr)
})

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

## Matrix to (row, column, value)
matrix_to_long_table  <- function(mat, dnn = names(dimnames(mat)), value_name = "value", ...) {
    if (is.null(dnn)) {
        dnn <- c("row", "col")
    }
    mat %>% as_tibble(rownames = dnn[[1]]) %>%
        pivot_longer(-!!dnn[[1]], names_to = dnn[[2]], values_to = value_name, ...)
}

## (row, column, value) to matrix
long_table_to_matrix <- function(tab, dnn, value_name, ...) {
    tab[c(dnn, value_name)] %>%
        pivot_wider(id_cols = !!dnn[[1]], names_from = !!dnn[[2]], values_from = !!value_name, ...) %>%
        { set_rownames(as.matrix(select(., -!!dnn[[1]])), .[[dnn[[1]]]]) }
}

main <- function() {
    parser <- ArgumentParser(description='Compute correlations between NGSCheckMate VAF files.')
    parser$add_argument('-I','--inputDir',
                        metavar='input_dir_name', required=TRUE, dest='inputdirname', action='store',
                        help='Input directory that contains the output VAF files of ngscheckmate_fastq or ncm_fastq.py')
    parser$add_argument('-O','--outdir',
                        metavar='output_dir',required=TRUE,dest='outdir',action='store',
                        help='Output directory')
    parser$add_argument('-N','--outfilename',
                        metavar='output_filename',dest='outfilename',action='store',default="output",
                        help='Output file prefix (default : "output")')
    parser$add_argument('-f','--family_cutoff',
                        dest='family_cutoff',action='store_true',
                        help='Use strict VAF correlation cutoffs. Recommended when your data may include related individuals (parents-child, siblings)')
    parser$add_argument('-nz','--nonzero',
                        dest='nonzero_read',action='store_true',
                        help='Use the mean of non-zero depths across the SNPs as a reference depth (default: Use the mean depth across all the SNPs)')
    ## test_args <- c("-f", "-I", "/sc/arion/projects/mscic/results/Ryan/NGSCheckMate/combined/vaf", "-O", "~/temp/ncm_test/", "-N", "WGS")
    args <- parser$parse_args(commandArgs(TRUE))

    tsmsg("Finding VAF files")
    vaf_files <- dir_ls(args$inputdirname, glob="*.vaf")# %>% head
    assert_that(length(vaf_files) >= 2, msg = str_c("Need at least 2 VAF files in ", shQuote(args$inputdirname)))
    names(vaf_files) <- path_ext_remove(path_file(vaf_files))

    tsmsg("Reading VAF files")
    col_spec <- cols(
        index = col_double(),
        ref = col_double(),
        alt = col_double(),
        vaf = col_double()
    )
    vaf_tables <- lapply(vaf_files, read_tsv, col_types = col_spec)
    vaf_full <- bind_rows(vaf_tables, .id = "sample")

    tsmsg("Computing pairwise correlations")
    ## row_var <- "index"
    ## col_var <- "sample"
    mat_vars <- c("ref", "alt", "vaf")
    data_mats <- lapply(mat_vars, long_table_to_matrix, tab = vaf_full, dnn = c("index", "sample")) %>%
        set_names(mat_vars)

    cormat <- cor(data_mats$vaf, use = "pairwise.complete.obs", method = "pearson")

    tsmsg("Computing mean depths")
    ## Values copied from the Python script
    predefined_model <- read_csv("
    family, min_depth, p1V,p1S, p0V, p0S
    true, 10,  0.874546, 0.022211, 0.646256175, 0.021336239
    true, 5,   0.785249, 0.021017, 0.598277053, 0.02253561
    true, 2,   0.650573, 0.018699, 0.536020197, 0.020461932
    true, 1,   0.578386, 0.018526, 0.49497342,  0.022346597
    true, 0.5, 0.529327, 0.025785, 0.465275173, 0.028221203
    true, 0,   0.529327, 0.025785, 0.465275173, 0.028221203
    false, 10, 0.874546, 0.022211, 0.310549, 0.060058
    false, 5, 0.785249,0.021017, 0.279778, 0.054104
    false, 2, 0.650573, 0.018699,0.238972, 0.047196
    false, 1, 0.578386,0.018526, 0.222322, 0.041186
    false, 0.5, 0.529327,0.025785, 0.217839, 0.040334
    false, 0, 0.529327,0.025785, 0.217839, 0.040334") %>%
        filter(family == args$family_cutoff)

    assign_depth_group <- function(x) {
        cut(x, breaks = c(predefined_model$min_depth, Inf), right = FALSE)
    }

    predefined_model %<>%
        mutate(depth_group = assign_depth_group(min_depth))

    depth_table <- {
        vaf_full %>%
            filter(!is.na(ref), !is.na(alt)) %>%
            mutate(depth = ref + alt) %>%
            filter(ifelse(args$nonzero_read, depth > 0, TRUE)) %>%
            group_by(sample) %>%
            summarise(mean_depth = mean(depth),
                      .groups = "drop")
    }

    tsmsg("Classifying pairs")
    ## Convert to table with columns (sampleA, sampleB, cor), then add
    ## depth info, etc.
    match_table <- cormat %>% matrix_to_long_table(c("sampleA", "sampleB"), "cor") %>%
        ## Add depth for both samples
        inner_join(depth_table %>% set_colnames(str_c(colnames(.), "A")),
                   by = "sampleA") %>%
        inner_join(depth_table %>% set_colnames(str_c(colnames(.), "B")),
                   by = "sampleB") %>%
        ## Compute the smaller of the 2 depths and assign it to a
        ## parameter group
        mutate(
            smaller_mean_depth = pmin(mean_depthA, mean_depthB),
            depth_group = assign_depth_group(smaller_mean_depth)) %>%
        ## Bring in the parameters for the assigned depth group
        inner_join(predefined_model, by="depth_group") %>%
        ## Compute the classification score. I differ here from the
        ## Python script, which allows "score1" and "score0" to be
        ## negative for some reason.
        mutate(score_match = pmax(abs(p1V - cor) - p1S, 0) ,
               score_unmatch = pmax(abs(p0V - cor) - p0S, 0),
               score_ratio = score_unmatch / score_match,
               matched = score_ratio > 1)

    tsmsg("Saving output files")
    ## Produce output files matching the original script
    match_output_table <- match_table %>%
        select(sampleA, matched, sampleB, cor, smaller_mean_depth) %>%
        mutate(matched = ifelse(matched, "matched", "unmatched")) %>%
        ## Each pair is represented twice, so filter out the 2nd one.
        ## Also filter out self-pairs
        mutate(ord_sampleA = match(sampleA, names(vaf_tables)),
               ord_sampleB = match(sampleB, names(vaf_tables))) %>%
        filter(ord_sampleA < ord_sampleB) %>%
        arrange(ord_sampleA, ord_sampleB) %>%
        select(-ord_sampleA, -ord_sampleB)
    out_file <- path(args$outdir, str_c(args$outfilename, "_all.txt"))
    write_tsv(match_output_table, out_file, col_names = FALSE)

    match_mat <- match_table %>%
        select(sampleA, sampleB, cor, matched) %>%
        ## Set correlations for non-matched pairs and self-matches to
        ## 0, to be consistent with vaf_ncm.py. Unlike vaf_ncm.py, I
        ## produce a symmetrical matrix and leave the diagonal (i.e.
        ## self matches) set to 1.
        mutate(cor = ifelse(matched, cor, 0)) %>%
        long_table_to_matrix(c("sampleA", "sampleB"), "cor") %>%
        as_tibble(rownames="sample_ID")
    matrix_out_file <- path(args$outdir, str_c(args$outfilename, "_corr_matrix.txt"))
    write_tsv(match_mat, matrix_out_file)

    tsmsg("Plotting")
    ## Reproduce the plot that is generated by the python script
    ## writing an R script and ruinning it. Obviously we don't need a
    ## separate script.
    output_corr_matrix <- read.delim(matrix_out_file)
    data <- output_corr_matrix
    ## This line is new, a correction
    data[is.na(data)] <- 0
    d3 <- as.dist((1 - data[,-1]))
    clust3 <- hclust(d3, method = "average")
    pdf_file <- path(args$outdir, str_c(args$outfilename, ".pdf"))
    pdf(pdf_file, width=10, height=7)
    op = par(bg = "gray85")
    par(plt=c(0.05, 0.95, 0.5, 0.9))
    plot(clust3, lwd = 2, lty = 1,cex=0.8, xlab="Samples", sub = "",  ylab="Distance (1-Pearson correlation)",hang = -1, axes = FALSE)
    axis(side = 2, at = seq(0, 1, 0.2), labels = FALSE, lwd = 2)
    mtext(seq(0, 1, 0.2), side = 2, at = seq(0, 1, 0.2), line = 1,   las = 2)
    dev.off()

    tsmsg("Done")
}

suppressWarnings(main())

EOF <- NULL
EOF
