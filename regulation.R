#' Parse transcription factor binding site hits
#'
#' @param hit_table A data table. Must contain named columns:
#'      Region: Hyphen-separated locus tags (ex: "SCO2781-SCO2782")
#'      Strand: + or -
#' @param min_score A numeric. Hits with scores below this value will
#'      be removed to improve speeds. Set to 0 to keep all hits.
#'
#' @return A data table. The region column is replaced by `Upstream` and `Downstream`
#'      and the table is sorted by descending score
#'
#' @export
parseTFBS <- function(hit_table, min_score = 5) {
  # Remove very low scores for speed
  hit_table <- hit_table[hit_table$Score > min_score,]
  # Ensure that the table is sorted
  hit_table <- hit_table[order(hit_table$Score, decreasing = T),]
  # Remove duplicates
  hit_table <- hit_table[!duplicated(hit_table["Region"]),]
  # Remove sites at chromosome end
  hit_table <- hit_table[grep("-", hit_table$Region),]
  # Add gene columns from $Region
  region_genes <- strsplit(hit_table$Region, "-", fixed = T)
  hit_table$Gene1 <- unlist(lapply(region_genes, `[`, 1))
  hit_table$Gene2 <- unlist(lapply(region_genes, `[`, 2))

  return(hit_table)
}

#' Filter TFBSs that aren't actually upstream of a gene and duplicates
#'
#' @param hits A data table: the output from parseTFBS()
#' @param strand_info A data table. Must contain named columns:
#'      Geneid: Locus tags
#'      Strand: + or -
#'
#' @export
filterTFBS <- function(hits, strand_info) {
  get_strand <- function(x) strand_info[strand_info$Geneid == x,"Strand"]
  # Gene1 should be negative and gene2 positive, else -> NA
  upstrands <- unlist(lapply(hits$Gene1, get_strand))
  hits$Gene1[upstrands != "-"] <- NA
  downstrands <- unlist(lapply(hits$Gene2, get_strand))
  hits$Gene2[downstrands != "+"] <- NA
  good_hits <- hits[!is.na(hits$Gene1) | !is.na(hits$Gene2), ]
  return(good_hits)
}


#' Simple operon prediction by strand and intergenic distance
#'
#' @param gene_info A data.frame of per-gene info with required columns:
#'      $Geneid: Locus tags
#'      $Start, $End: Gene boundaries (strings will be converted)
#'      $Strand: Can be given as any of boolean, strings, etc.
#' @param igd A two-length numeric vector: the min and max intergenic distances
#'      for two genes to be considered co-operonic
#' @param flip_minus A string or NULL: If non-null, gene vectors will be
#'      reversed for operons where `gene_info$Strand` equals `flip_minus`
#'
#' @return A list of gene name vectors
#'
#' @export
predictOperons <- function(gene_info, igd = c(-50, 150), flip_minus = "-") {
  min_dist <- igd[1]
  max_dist <- igd[2]

  len <- nrow(gene_info)

  # Initialize empty vector to store operon numbers
  membership <- vector(mode = "integer", length = len)

  # Get relevant columns as vectors
  starts <- suppressWarnings(as.integer(gene_info$Start)) # See next block
  ends <- suppressWarnings(as.integer(gene_info$End))
  strands <- gene_info$Strand

  # Genes with missing/bad locations get removed
  nas <- union(which(is.na(starts)), which(is.na(ends)))
  if (length(nas) > 0) {
    # Warn the user
    bads <- paste(sort(gene_info$Geneid[nas]), collapse = ", ")
    warning(paste("Dropping genes where location could not be parsed:", bads))
  }

  # Prepare for loop
  membership[1] <- 1L
  current <- 1L

  # Start loop with gene 2
  for (this in 2:len) {
    # Drop genes with bad locations
    if (this %in% nas) {
      current <- current + 1L
      next
    }

    # this gene must have the same strand as the last gene
    if (strands[this] != strands[this - 1]) {
      current <- current + 1L
      membership[this] <- current
      next
    }
    # this gene must have the right intergenic distance with the last gene
    gene_dist <- starts[this] - ends[this - 1]
    if (is.na(gene_dist) || gene_dist > max_dist || gene_dist < min_dist) {
      current <- current + 1L
      membership[this] <- current
      next
    }
    # this gene has passed both checks and is part of the same operon
    membership[this] <- current
  }

  # Convert membership vector into a list of index vectors then locus tags
  operons <- lapply(1:current, function(x) which(membership == x))
  operons <- Filter(length, operons)  # Removes NAs

  # Flip genes on negative strand if desired
  if (!is.null(flip_minus)) {
    operons <- lapply(operons, function(x) {
      if (strands[x][1] == flip_minus) rev(x)
      else x
    })
  }

  named <- lapply(operons, function(x) gene_info$Geneid[x])

  return(named)
}


#' Split operons based on coexpression data
#'
#' TODO: Consider allowing filtered genes to remain
#' TODO: vectorize for speed?
#'
#' @param operons A list of gene name vectors
#' @param cors A symmetric matrix of gene correlations with genenames as
#'      dimnames
#' @param min_cor A double: Minimal PCC between neighbors for clustering
#' @param keep_missing A boolean: Should genes with no correlation information
#'      be kept (as a one-gene operon)?
#'
#' @return A new list of gene name vectors
#' @export
splitByCoexpression <- function(operons, cors, min_cor = 0.5, keep_missing = F) {
  if (min_cor >= 1 || min_cor <= 0) stop("'min_cor' should be between 0 and 1")

  good_genes <- rownames(cors)

  in_cors <- lapply(operons, function(x) x %in% good_genes)

  # Save genes with no correlation info if it's wanted
  if (keep_missing) missing <- unlist(operons)[!unlist(in_cors)]

  # Remove operons that are all F
  operons <- operons[lapply(in_cors, sum) > 0]
  # Keep in_cors up-to-date
  in_cors <- in_cors[lapply(in_cors, sum) > 0]

  # Operons of length 1 are trivial
  bool_single <- lapply(operons, length) == 1
  singles <- operons[bool_single]

  operons <- operons[!bool_single]
  in_cors <- in_cors[!bool_single]


  new_ops <- vector(mode = "list", length = length(operons))

  # First pass: get rid of genes not in cors, and split operons around them
  for (i in seq_along(operons)) {
    op_genes <- operons[[i]]
    op_in_cors <- in_cors[[i]]

    # Only keep genes where in_cors == T, and split if there's a F in between
    new_ops[[i]] <- split(op_genes[op_in_cors],
                          data.table::rleidv(op_in_cors)[op_in_cors])
  }
  # Flatten the list of list of vectors to a list of vectors
  new_ops <- unlist(new_ops, recursive = F)

  # Remove new singles
  bool_single <- lapply(new_ops, length) == 1
  singles <- c(singles, new_ops[bool_single])

  operons <- new_ops[!bool_single]
  new_ops <- vector(mode = "list", length = length(operons))

  # These operons now are all in_cors and can be split by coexpression
  for (i in seq_along(operons)) {
    op_genes <- operons[[i]]
    club <- rep(1, length(op_genes))
    current <- 1
    for (this in 2:length(op_genes)) {
      corr <- cors[ op_genes[this], op_genes[this - 1] ]
      if (corr < min_cor) current <- current + 1
      club[this] <- current
    }
    # All genes are now assigned a 'club' and can be split accordingly
    new_ops[[i]] <- split(op_genes, club)
  }

  # Flatten the list of list of vectors to a list of vectors
  new_ops <- unlist(new_ops, recursive = F)

  # Add singles back in
  new_ops <- c(new_ops, singles)

  # If the no-correlation genes were requested, put those in
  if (keep_missing) new_ops <- c(new_ops, missing)

  # Sort the list, remove the inconsistent names, (and return)
  unname(new_ops[order(unlist(lapply(new_ops, utils::head, 1)))])
}


#' Given a vector of locus tags, cluster them based on their number
#'
#' Genes are first sorted (alphabetically).
#'
#' @param tags A character vector of locus tags
#' @param max_diff An integer: The maximum difference in locus tag numbers
#'      for two genes to be clustered
#'
#' @return A list of character vectors containing gene names
#'
#' @export
clusterTags <- function(tags, max_diff = 1) {
  tags <- tags[order(tags)]
  numbers <- extractGeneNumber(tags)
  diffs <- diff(numbers)

  # Initialize output vector
  clusters <- rep(1, length(tags))
  counter <- 1

  for (i in seq_along(diffs)) {
    if (diffs[i] > max_diff) counter <- counter + 1
    clusters[i+1] <- counter
  }

  cluster_ind <- unique(clusters)
  clusters <- lapply(cluster_ind, function(x) tags[clusters == x])
}


#' Extract gene numbers from locus tags
#'
#' Uses regex to extract the final number from the strings
#' Note: Letter suffixes will be removed
#'
#' @param genenames A character vector
#'
#' @return A numeric vector
#'
#' @export
extractGeneNumber <- function(genenames) {
  # Extract all numbers
  regex <- gregexpr("\\d+", genenames, perl = TRUE)
  # Nums is a list of lists
  nums <- regmatches(genenames, regex)

  # regex will return -1 if there are no matches - accept empty strings
  if (any(unlist(regex) == -1 & genenames != "")) {
    stop(paste("Unable to extract gene numbers from locus tag(s):",
               paste(genenames[unlist(regex) == -1], collapse = ", ")))
  }
  # Get the last item of each list (-> char vector)
  nums <- sapply(nums, utils::tail, 1)
  # Convert to numbers
  nums <- as.numeric(nums)

  return(nums)
}


#' Extract rows in a matrix by the first and last desired row
#'
#' @param mat A matrix
#' @param from,to Strings: the first and last desired rows
#' @param match_col A string or integer: the colname or index, respectively,
#'      of the column containing the names 'from' and 'to'. If NULL, then
#'      'rownames(mat)' will be used.
#' @param padding An integer: The number of extra genes to add before and after
#' @param and_col A boolean: Should the columns of a symmetric matrix be sliced
#'      as well?
#'
#' @return A matrix, equivalent to mat\[from:to,\]
#'
#' @export
getRowRange <- function(mat, from, to, match_col = NULL, padding = 0, and_col = F) {
  if (is.null(match_col)) names <- rownames(mat)
  else names <- mat[,match_col]

  ind1 <- which(names == from) - padding
  ind2 <- which(names == to) + padding

  # Check for missing rownames
  bad <- c(from, to)[c(length(ind1) == 0, length(ind2) == 0)]
  if (length(bad) != 0) {
    stop(paste0("Requested row name(s) not found in matrix: ",
                paste(bad, collapse = ", ")))
  }
  if (and_col) return(mat[ind1:ind2, ind1:ind2])
  else return(mat[ind1:ind2,])
}