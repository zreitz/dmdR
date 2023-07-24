#' Filter out genes with little/no expression
#'
#' @param mat A matrix of gene expression
#' @param cutoff Minimal required expression level
#' @param fraction Fraction of samples that must meet `cutoff`, or else the gene is
#'      discarded.
#'
#' @return Filtered matrix
#' @export
filterLowExpression <- function(mat, cutoff = 10, fraction = 0.5) {

  # minimum amount of samples that need to have at least `cutoff` expression level
  nsamples <- ncol(mat)
  min_samples <- ceiling(nsamples * fraction)

  count_above_filter <- apply(mat, 1, function(x) length(which(x >= cutoff)))
  keep <- mat[which(count_above_filter >= min_samples),]
  message(paste0(nrow(mat) - nrow(keep), " low-expression items removed (",
                  floor(100*(nrow(mat) - nrow(keep))/nrow(mat)), "%)"))

  return(keep)
}

#' Normalize across samples with TMM
#'
#' Wrapper for EdgeR implementation
#'
#' @param counts A matrix or data.frame with genes in rows
#' @param replicates An integer: the number of replicates in each group
#'
#' @return A matrix with the same dimensions as `counts`
#'
#' @import edgeR
#' @export
logTMM <- function(counts, replicates) {
  group <- factor(rep(letters[1 : (ncol(counts)/replicates)], each = replicates))
  dge <- DGEList(counts = counts, group = group)
  dge <- calcNormFactors(dge, method = "TMM")
  tmm <- cpm(dge)
  ## Normalize using asinh (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9)
  log2(tmm + sqrt(tmm^2 + 1))
}


#' Bucket rows of a matrix
#'
#' @param expr A matrix
#' @param bin_nr An integer, the desired number of bins
#' @param symmetric A boolean: Should the matrix be binned in two dimensions
#'      to produce square bins along the main diagonal?
#'
#' @return A list of matrices
#' @export
binMatrix <- function(expr, bin_nr, symmetric = F) {
  indices <- split(seq_len(nrow(expr)), cut(seq_len(nrow(expr)),
                                            bin_nr, labels = FALSE))
  names(indices) <- NULL
  if (symmetric) lapply(indices, function(x) expr[x,x])
      else lapply(indices, function(x) expr[x,])
}


#' Reduce confounding artifacts by PC removal
#'
#' https://doi.org/10.1186/s13059-019-1700-9
#'
#' @param rpkm A matrix of log2-RPKM values with genes in rows
#' @param removed If non-null, manually selects the number of PCs to be removed
#'
#' @return A matrix of normalized expression with genes in rows
#'
#' @export
removePCs <- function(rpkm, removed = NULL) {
  # Determine how many PCs to be removed
  if (is.null(removed)) {
    mod <- matrix(1,nrow=dim(rpkm)[2],ncol = 1)
    colnames(mod) <- "Intercept"
    removed <- sva::num.sv(rpkm, mod, method = "be")
  }
  message(paste("Removing", removed, "principal components"))
  t(sva::sva_network(t(rpkm), removed))
}


#' Divide a correlation matrix into bins and plot the per-bin cor distribution
#'
#' @param cors A square correlation matrix pre-sorted by mean expression (in both dims)
#' @param max_sample An integer: The maximum number of correlations (per-bin) used to
#'      produce the density plots. Lower numbers saves time and memory but produces
#'      more lumpy plots
#' @param title A string, or NULL for no title
#'
#' @import ggridges
#'
#' @export
plotBinnedCorDensity <- function(cors, max_sample = 1e6, title = NULL) {
  # Appease R CMD check
  values <- ind <- x <- NULL

  # Break the correlation matrix into bins by mean expression
  cor_bins <- binMatrix(cors, 10, T)

  cor_list <- vector(mode = "list", length = 10)
  names(cor_list) <- paste0(seq(from = 0, by = 10, length.out = 10), "-",
                            seq(from = 9, by = 10, length.out = 10), "%")

  print("Medians:")
  for (i in seq_along(cor_bins)) {
    bin_cors <- cor_bins[[i]]
    ## remove duplicates, convert to vector
    bin_cors <- as.vector(bin_cors[lower.tri(bin_cors)])
    ## Sample
    if (length(bin_cors) > max_sample) {
      bin_cors <- sample(bin_cors, max_sample)
    }
    print(paste(names(cor_list)[[i]], stats::median(bin_cors), sep = ": "))
    cor_list[[i]] <- bin_cors
  }

  cor_table <- utils::stack(cor_list)

  cor_colors <- c('#67001F', '#B2182B', '#D6604D', '#F4A582',
                  '#FDDBC7', '#FFFFFF', '#D1E5F0', '#92C5DE',
                  '#4393C3', '#2166AC', '#053061')

  if (is.null(title)) title <- "Correlation density by mean gene expression"
  ggplot(cor_table, aes(x = values, y = ind, fill = stat(x))) +
     geom_segment(mapping = aes(x = -1, xend = 1, y = ind, yend = ind)) +
     geom_density_ridges_gradient(scale = 1, quantile_lines = TRUE, show.legend = F) +
     scale_fill_gradientn(colours = cor_colors, limits = c(-1,1)) +
     theme_ridges(center_axis_labels = T) +
     scale_x_continuous(limits = c(-1,1)) +
     scale_y_discrete(expand = expansion(mult = c(0.03, 0.05))) +
     geom_vline(xintercept = 0, linetype = "dashed") +
     xlab("Pearson correlation coefficient density") +
     ylab("Mean gene expression percentile") +
     ggtitle(title) +
     coord_flip()
}