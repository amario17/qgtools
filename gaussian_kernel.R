#' Gaussian kernel from marker data, with optional bandwidth optimisation
#'
#' Builds the Gaussian kernel between individuals (lines, hybrids or
#' environments) from a numeric matrix, usually SNP markers:
#' \deqn{K_{ij} = \exp(-h \, d^2_{ij} / q)}
#' where \eqn{d^2_{ij}} is the squared Euclidean distance between rows
#' \eqn{i} and \eqn{j}, \eqn{q} is a quantile of all \eqn{d^2} values
#' (the median by default, computed over the full \eqn{n \times n} matrix,
#' diagonal included) and \eqn{h} is the bandwidth. Larger \eqn{h} makes
#' the kernel decay faster with distance (individuals look less related).
#'
#' The bandwidth can be fixed by the user or optimised from a phenotype
#' vector \code{y} (e.g. genotype BLUEs). In the second case, \eqn{h} and
#' \eqn{\phi = \sigma^2_u / \sigma^2_e} are chosen by maximising the
#' marginal likelihood of \eqn{y = \mu + u + e}, \eqn{u \sim N(0,
#' \sigma^2_u K_h)}, with a gamma prior on \eqn{h} (Pérez-Elizalde et al.
#' 2015), using the spectral decomposition of \eqn{K_h}.
#'
#' Squared distances are computed with \code{tcrossprod()} instead of
#' \code{dist()}, which is much faster for large marker panels (hundreds
#' of thousands of SNPs) and gives the same values up to rounding error.
#'
#' @param X Numeric matrix (individuals x markers), without missing values.
#'   Row names are kept in the kernel.
#' @param y Optional numeric vector with one phenotype per row of \code{X}
#'   (e.g. BLUEs). If named, it is reordered by \code{rownames(X)}. Used
#'   only when \code{h} is \code{NULL}.
#' @param h Optional fixed bandwidth. If given, it is used as is (and
#'   \code{y} is ignored). If both \code{h} and \code{y} are \code{NULL},
#'   \code{h = 1}.
#' @param prob Quantile of \eqn{d^2} used to scale the distances
#'   (default 0.5, the median). If it is zero, the 0.05 quantile is used.
#' @param prior_shape,prior_rate Shape and rate of the gamma prior on
#'   \eqn{h} (default 3 and 1.5).
#' @param start Starting values for (\eqn{h}, \eqn{\phi}) (default 1, 3).
#' @param lower,upper Bounds for (\eqn{h}, \eqn{\phi}) in
#'   \code{optim(method = "L-BFGS-B")} (default 0.05-6 and 0.05-30).
#' @param nu,Sc Hyperparameters of the residual part of the marginal
#'   likelihood (default 1e-4).
#'
#' @return A list with:
#' \describe{
#'   \item{K}{Gaussian kernel (n x n), with the row names of \code{X}.}
#'   \item{h}{Bandwidth used.}
#'   \item{h_source}{"fixed", "optimised" or "default".}
#'   \item{q}{Quantile of \eqn{d^2} used to scale the distances.}
#'   \item{phi}{Optimised \eqn{\sigma^2_u / \sigma^2_e} (NA if \eqn{h} was
#'     not optimised).}
#'   \item{optim}{Output of \code{optim()} (NULL if \eqn{h} was not
#'     optimised).}
#' }
#'
#' @note
#' The kernel and the bandwidth optimisation follow \code{GK_Kernel} by
#' Germano Costa-Neto (Costa-Neto et al. 2021), available at
#' \url{https://github.com/gcostaneto/KernelMethods}. This function is a
#' rewrite of the published method: same formula, prior, bounds and
#' starting values, with faster distances, a single matrix as input, and
#' the bandwidth, quantile and \eqn{\phi} returned with the kernel. With
#' the same inputs it reproduces \code{GK_Kernel} (maximum difference
#' around 1e-10 on 1,691 maize hybrids x 302,300 SNPs).
#'
#' @references
#' Costa-Neto, G., Fritsche-Neto, R., Crossa, J. (2021). Nonlinear kernels,
#' dominance, and envirotyping data increase the accuracy of genome-based
#' prediction in multi-environment trials. \emph{Heredity}, 126, 92-106.
#' \doi{10.1038/s41437-020-00353-1}
#'
#' Pérez-Elizalde, S., Cuevas, J., Pérez-Rodríguez, P., Crossa, J. (2015).
#' Selection of the bandwidth parameter in a Bayesian kernel regression
#' model for genomic-enabled prediction. \emph{Journal of Agricultural,
#' Biological, and Environmental Statistics}, 20, 512-532.
#' \doi{10.1007/s13253-015-0229-y}
#'
#' @author João Marcos Amario, PhD (\email{joao.m.sousa@ufv.br})
#'
#' @examples
#' \dontrun{
#' # M: hybrids x SNPs (0, 1, 2); blue: named vector of hybrid BLUEs
#' gk <- gaussian_kernel(M, y = blue)      # bandwidth optimised
#' gk$h
#' GK <- gk$K
#' Gg <- 0.99 * GK + 0.01 * diag(nrow(GK)) # invertible, for vm() in ASReml
#'
#' gk1 <- gaussian_kernel(M, h = 1)        # fixed bandwidth
#' }
#'
#' @export
gaussian_kernel <- function(X, y = NULL, h = NULL, prob = 0.5,
                            prior_shape = 3, prior_rate = 1.5,
                            start = c(1, 3), lower = c(0.05, 0.05), upper = c(6, 30),
                            nu = 1e-4, Sc = 1e-4) {

  X <- as.matrix(X)
  if (anyNA(X)) stop("X has missing values; impute them first")

  # squared Euclidean distances: d2 = diag(XX') + diag(XX')' - 2XX'
  XX <- tcrossprod(X)
  D  <- outer(diag(XX), diag(XX), "+") - 2 * XX
  rm(XX)
  D[D < 0] <- 0      # rounding error
  diag(D)  <- 0
  q <- as.numeric(quantile(D, probs = prob))
  if (q == 0) q <- as.numeric(quantile(D, probs = 0.05))

  # negative log marginal likelihood (+ gamma prior on h)
  marg <- function(theta, y, D, q) {
    Kh  <- exp(-theta[1] * D / q)
    eg  <- eigen(Kh, symmetric = TRUE)
    nr  <- sum(eg$values > 1e-10)
    S   <- eg$values[1:nr]
    d   <- crossprod(eg$vectors[, 1:nr], y - mean(y))
    ll  <- -0.5 * sum(log(1 + theta[2] * S)) - (nu + nr - 1) / 2 * log(Sc + sum(d^2 / (1 + theta[2] * S)))
    -(ll + dgamma(theta[1], prior_shape, prior_rate, log = TRUE))
  }

  sol <- NULL; phi <- NA
  if (!is.null(h)) {
    h_source <- "fixed"
  } else if (!is.null(y)) {
    if (!is.null(names(y)) && !is.null(rownames(X))) y <- y[rownames(X)]
    if (length(y) != nrow(X) || anyNA(y)) stop("y must have one non-missing value per row of X")
    sol <- optim(start, marg, y = as.numeric(y), D = D, q = q, method = "L-BFGS-B",
                 lower = lower, upper = upper)
    h <- sol$par[1]; phi <- sol$par[2]
    h_source <- "optimised"
  } else {
    h <- 1
    h_source <- "default"
  }

  K <- exp(-h * D / q)
  dimnames(K) <- list(rownames(X), rownames(X))
  list(K = K, h = h, h_source = h_source, q = q, phi = phi, optim = sol)
}
