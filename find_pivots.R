# find_pivots.R --------------------------------------------------------------
#
# Identify PIVOT locations in a Target Population of Environments (TPE):
# the subset of locations that best represents the whole trial network.
#
# Pipeline:
#   environmental covariates -> quantiles per location
#   -> PLS against a phenotypic kernel
#   -> Gaussian kernel between locations
#   -> subset selection by PEVMEAN (genetic algorithm, STPGA)
#
# NOTE: a "pivot" is the most REPRESENTATIVE location, i.e. the one whose
# information is largely shared with the others. Locations left out by the
# algorithm are the most environmentally distant (unique) ones.
#
# Dependencies: plyr, reshape2, plsdepot, EnvRtype, STPGA
# ----------------------------------------------------------------------------

find_pivots <- function(data,
                        loc,              # location column name
                        year,             # year column name
                        trait,            # phenotype column name (e.g. BLUE)
                        covs,             # environmental covariate column names
                        n_pivots = 5,
                        seed = 123) {

  set.seed(seed)

  d <- data[, c(loc, year, trait, covs)]
  names(d)[1:3] <- c("loc", "year", "y")

  # 1. Environmental matrix: location x (covariate quantiles), standardised
  S <- reshape2::melt(d, measure.vars = covs)
  S <- plyr::ddply(S, plyr::.(loc, variable), plyr::summarise,
                   q10 = quantile(value, .10, na.rm = TRUE),
                   q50 = quantile(value, .50, na.rm = TRUE),
                   q90 = quantile(value, .90, na.rm = TRUE))
  S <- reshape2::melt(S, measure.vars = c("q10", "q50", "q90"),
                      variable.name = "qt")
  S <- reshape2::acast(S, loc ~ variable + qt, mean, value.var = "value")
  S <- scale(S)
  S <- S[, colSums(is.na(S)) == 0, drop = FALSE]      # drop constant columns

  # 2. Phenotypic kernel between locations (PLS response)
  S0 <- reshape2::melt(d[, c("loc", "year", "y")], measure.vars = "y")
  S0 <- plyr::ddply(S0, plyr::.(loc, year), plyr::summarise,
                    q10 = quantile(value, .10, na.rm = TRUE),
                    q50 = quantile(value, .50, na.rm = TRUE),
                    q90 = quantile(value, .90, na.rm = TRUE))
  S0 <- reshape2::melt(S0, measure.vars = c("q10", "q50", "q90"),
                       variable.name = "qt")
  S0 <- reshape2::acast(S0, loc ~ qt, mean, value.var = "value")
  S0 <- scale(S0)

  K0 <- EnvRtype::env_kernel(env.data = S0, gaussian = TRUE)[[2]]

  # 3. PLS: environmental covariates -> phenotypic similarity
  pls <- plsdepot::plsreg2(predictors = S, responses = K0,
                           comps = min(nrow(S) - 1, ncol(S), 25),
                           crosval = FALSE)

  # 4. Kernel between locations, built from the PLS weights
  K <- EnvRtype::env_kernel(env.data = pls$std.coefs, gaussian = TRUE)[[1]]

  # 5. Pivots: subset that minimises the mean prediction error variance
  pivots <- STPGA::GenAlgForSubsetSelectionNoTest(
    P = K, ntoselect = n_pivots, InitPop = NULL,
    npop = 100, nelite = 5, mutprob = .5, mutintensity = 1,
    niterations = 200, minitbefstop = 20,
    tabu = FALSE, tabumemsize = 0, plotiters = FALSE,
    lambda = 1e-5, errorstat = "PEVMEAN", mc.cores = 1)[[1]]

  as.character(pivots)
}

# ---- Example ---------------------------------------------------------------
# data <- readRDS("my_data.rds")
# covs <- names(data)[10:110]          # environmental covariate columns
#
# find_pivots(data,
#             loc   = "Field.Location",
#             year  = "Year",
#             trait = "BLUE",
#             covs  = covs,
#             n_pivots = 5)
# ----------------------------------------------------------------------------
