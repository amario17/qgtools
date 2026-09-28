#' hybridForge: forge hybrid (F1) genotypes from parental line genotypes
#'
#' Reconstructs the SNP-level genotype of each hybrid (F1) from the genotypes
#' of its two parental inbred lines, using a strict homozygous-parent rule:
#' a hybrid genotype is only called when BOTH parents have a confirmed
#' homozygous call at that marker. Any locus where a parent is heterozygous,
#' missing, or otherwise ambiguous is returned as \code{NA} rather than
#' approximated, so the output never contains fractional dosages.
#'
#' @details
#' For each hybrid \eqn{i} and marker \eqn{j}, with parental dosages
#' \eqn{P1_{ij}} and \eqn{P2_{ij}} coded as 0 (homozygous reference),
#' 1 (heterozygous) or 2 (homozygous alternate):
#' \itemize{
#'   \item \eqn{P1 = 2} and \eqn{P2 = 2} -> hybrid = 2
#'   \item \eqn{P1 = 0} and \eqn{P2 = 0} -> hybrid = 0
#'   \item \eqn{(P1 = 2, P2 = 0)} or \eqn{(P1 = 0, P2 = 2)} -> hybrid = 1
#'   \item any other combination (heterozygous parent, missing data, etc.) -> \code{NA}
#' }
#' This mirrors the classic F1-from-inbred-parents construction used in
#' plant breeding genomic prediction pipelines (e.g. Dias, K.O.G., 2018,
#' unpublished script "Create Hybrids Genotypes"), rewritten here as a
#' vectorized, general-purpose function so it scales to genome-wide marker
#' panels (hundreds of thousands of SNPs) without a per-SNP/per-hybrid loop.
#'
#' @param geno A numeric matrix of parental line genotypes, dosage-coded as
#'   0/1/2 (missing values as \code{NA}). Rows are individuals (parental
#'   lines), columns are markers/SNPs. Row names must be the line
#'   identifiers referenced in \code{pedigree}.
#' @param pedigree A \code{data.frame} (or \code{data.table}) with one row
#'   per hybrid, containing at least the three columns named in
#'   \code{parent1_col}, \code{parent2_col} and \code{hybrid_col}.
#' @param parent1_col Character. Name of the column in \code{pedigree}
#'   holding the first parent's identifier (must match a row name of
#'   \code{geno}). Default \code{"Pai"}.
#' @param parent2_col Character. Name of the column in \code{pedigree}
#'   holding the second parent's identifier. Default \code{"Mae"}.
#' @param hybrid_col Character. Name of the column in \code{pedigree}
#'   holding the hybrid's own identifier; this becomes the row name of the
#'   returned matrix. Default \code{"ID"}.
#' @param verbose Logical. If \code{TRUE} (default), prints a short summary
#'   (how many hybrids were requested, how many had both parents genotyped,
#'   how many were dropped and why).
#'
#' @return A numeric matrix of hybrid genotypes, dosage-coded 0/1/2 with
#'   \code{NA} for undetermined markers. Rows are hybrids (named after
#'   \code{hybrid_col}), columns are the same markers as \code{geno}
#'   (same order, same column names).
#'
#' @examples
#' \dontrun{
#' geno <- readRDS("matriz_marcadores_para_G.rds")       # lines x SNPs, 0/1/2
#' pedigree <- read.csv("pedigree_para_matriz_A.csv", sep = ";")
#'
#' hybrid_geno <- hybridForge(
#'   geno      = geno,
#'   pedigree  = pedigree,
#'   parent1_col = "Pai",
#'   parent2_col = "Mae",
#'   hybrid_col  = "ID"
#' )
#'
#' dim(hybrid_geno)
#' hybrid_geno[1:5, 1:5]
#' }
#'
#' @author Adapted for general use from an original hybrid-genotype
#'   construction script by K.O.G. Dias (Jan/2018).
#' @export
hybridForge <- function(geno,
                         pedigree,
                         parent1_col = "Pai",
                         parent2_col = "Mae",
                         hybrid_col  = "ID",
                         verbose     = TRUE) {

  # ---- input validation --------------------------------------------------
  if (!is.matrix(geno)) stop("'geno' must be a matrix (individuals x markers).")
  if (is.null(rownames(geno))) stop("'geno' must have row names (line identifiers).")
  if (!all(c(parent1_col, parent2_col, hybrid_col) %in% names(pedigree))) {
    stop("'pedigree' must contain the columns: ", parent1_col, ", ", parent2_col, ", ", hybrid_col)
  }

  pedigree <- as.data.frame(pedigree)
  parent1 <- as.character(pedigree[[parent1_col]])
  parent2 <- as.character(pedigree[[parent2_col]])
  hybrid  <- as.character(pedigree[[hybrid_col]])

  n_total <- length(hybrid)

  # ---- keep only crosses with both parents present and hybrid named -----
  valido <- !is.na(parent1) & !is.na(parent2) & parent1 != "" & parent2 != "" &
            !is.na(hybrid) & hybrid != "" &
            parent1 %in% rownames(geno) & parent2 %in% rownames(geno)

  n_sem_pai <- sum(is.na(parent1) | parent1 == "" | is.na(parent2) | parent2 == "")
  n_pai_nao_genotipado <- sum(!valido) - n_sem_pai
  n_valido <- sum(valido)

  if (n_valido == 0) stop("No hybrid has both parents present in 'geno' row names.")

  parent1 <- parent1[valido]
  parent2 <- parent2[valido]
  hybrid  <- hybrid[valido]

  dup <- duplicated(hybrid)
  if (any(dup)) {
    warning(sum(dup), " duplicated hybrid identifier(s) found in '", hybrid_col,
            "' -- keeping the first occurrence of each.")
    parent1 <- parent1[!dup]
    parent2 <- parent2[!dup]
    hybrid  <- hybrid[!dup]
  }

  # ---- vectorized construction of hybrid genotypes -----------------------
  # rule: hybrid genotype is only called when BOTH parents are confirmed
  # homozygous (0 or 2); anything else (het parent, missing data) -> NA
  P1 <- geno[parent1, , drop = FALSE]
  P2 <- geno[parent2, , drop = FALSE]

  hybrid_geno <- matrix(NA_real_, nrow = nrow(P1), ncol = ncol(P1),
                         dimnames = list(hybrid, colnames(geno)))

  hybrid_geno[P1 == 2 & P2 == 2] <- 2
  hybrid_geno[P1 == 0 & P2 == 0] <- 0
  hybrid_geno[(P1 == 2 & P2 == 0) | (P1 == 0 & P2 == 2)] <- 1

  if (verbose) {
    pct_na <- round(100 * mean(is.na(hybrid_geno)), 1)
    cat("hybridForge():\n")
    cat("  Hybrids requested:                 ", n_total, "\n")
    cat("  Dropped (parent missing/unknown):  ", n_sem_pai, "\n")
    cat("  Dropped (parent not in 'geno'):    ", n_pai_nao_genotipado, "\n")
    cat("  Hybrids reconstructed:             ", nrow(hybrid_geno), "\n")
    cat("  Markers:                           ", ncol(hybrid_geno), "\n")
    cat("  Missing genotype calls (NA):       ", pct_na, "%\n")
  }

  hybrid_geno
}
