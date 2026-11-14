#' Post-processing of a reduced-rank factor-analytic (FA) mixed model
#'
#' Written for the reduced-rank form \code{rr(env, k):vm(gen, G) +
#' diag(env):gen}, where the genomic relationship enters only the common
#' factors and the specific variances are estimated in the \code{diag}
#' term (ASReml keeps them at zero in the \code{rr} term). It reads them
#' from there. It also works on \code{fa()} models, with the same result
#' as \code{fa_model()}.
#'
#' From a mixed model fitted in ASReml with a factor-analytic random term,
#' \code{gen:fa(env, k)} or \code{fa(env, k):vm(gen, G)} (with a genomic
#' relationship matrix), extracts and organises the main FA results: the
#' environment loadings and genotype scores rotated to canonical form via
#' singular value decomposition (SVD), the genetic covariance matrix
#' between environments (\code{Gvcov}) and the corresponding genetic
#' correlation matrix (\code{Gcor}) implied by the FA(k) structure, plus
#' goodness-of-fit diagnostics: total and factor/environment-wise
#' explained variance, the Average Semivariance Ratio (ASVR), AIC, BIC and
#' log-likelihood.
#'
#' Genotype-by-environment BLUPs are also returned in two forms: the
#' \strong{marginal} value (\eqn{\Lambda f}), corresponding only to the
#' common part of the factor-analytic structure, without the intercept or
#' the environment fixed effect - useful as an input for prediction in
#' untested environments (e.g. GIS-FA), where the overall mean and other
#' effects are added separately; and the \strong{conditional} value
#' (intercept + environment fixed effect + genotype-by-environment random
#' effect, including the specific deviation), which is the full genotype
#' prediction within the environments actually tested.
#'
#' @param model An ASReml model object (\code{asreml}), with a
#'   factor-analytic random term for some number of factors \code{k}:
#'   \code{gen:fa(env, k)}, \code{fa(env, k):vm(gen, G)}, or the reduced-rank
#'   form \code{rr(env, k):gen + diag(env):gen} (with or without \code{vm()}).
#'   With \code{vm(gen, G)}, set the genotype factor levels in the same
#'   order as the rows of \code{G} (\code{factor(gen, levels = rownames(G))});
#'   otherwise ASReml may label the \code{diag(env):gen} effects wrongly and
#'   the conditional predictions will not match the fitted values.
#' @param name.env Name (string) of the column identifying the
#'   environment in the model (e.g. "env").
#' @param name.gen Name (string) of the column identifying the genotype
#'   in the model (e.g. "gen").
#'
#' @return A list with:
#' \describe{
#'   \item{rot.loads}{Matrix of rotated environment loadings (environment
#'     x factor), in canonical form obtained via SVD.}
#'   \item{rot.scores}{Matrix of rotated genotype scores (genotype x
#'     factor).}
#'   \item{Gvcov}{Genetic covariance matrix between environments implied
#'     by the FA(k) model (rotated loadings plus each environment's
#'     specific variance).}
#'   \item{Gcor}{Genetic correlation matrix between environments,
#'     i.e. standardised \code{Gvcov}.}
#'   \item{diagnostics}{Named vector with \code{expvar} (percentage of
#'     total genetic variance explained by the common factors),
#'     \code{ASVR} (Average Semivariance Ratio, percentage), \code{aic},
#'     \code{bic} and \code{logl} (log-likelihood) of the model.}
#'   \item{expvar_j}{Matrix (factor x environment) with the percentage of
#'     variance explained by each factor within each environment.}
#'   \item{marginal}{Data frame (env, gen, marginal) with the common part
#'     (\eqn{\Lambda f}) of the genotypic prediction, without intercept
#'     or environment effect, by environment and genotype.}
#'   \item{blups}{Data frame (env, gen, conditional, marginal) joining
#'     the conditional prediction (full, within the tested environments)
#'     and the marginal one (\eqn{\Lambda f}), by environment and
#'     genotype.}
#' }
#'
#' @references
#' Piepho, H.P. (1997). Analyzing genotype-environment data by mixed
#' models with multiplicative terms. \emph{Biometrics}, 53(2), 761-766.
#'
#' Smith, A.B., Cullis, B.R., Thompson, R. (2001). Analyzing variety by
#' environment data using multiplicative mixed models and adjustments for
#' spatial field trend. \emph{Biometrics}, 57(4), 1138-1147.
#'
#'
#' #############################################################################
#' This function was adapted and modified from the \code{fa.outs} function
#' by Prof. Saulo Chaves, ESALQ/USP.
#' #############################################################################
#'
#'
#' @author João Marcos Amario, PhD (\email{joao.m.sousa@ufv.br})
#'
#' @examples
#' \dontrun{
#' Gg <- 0.99 * GK + 0.01 * diag(nrow(GK))   # genomic kernel, invertible
#' Mod1 <- asreml(
#'   fixed   = blue ~ env,
#'   random  = ~ rr(env, 4):vm(gen, Gg) + diag(env):gen,
#'   weights = W,
#'   family  = asr_gaussian(dispersion = 1),
#'   data    = blues,
#'   maxit   = 100
#' )
#' Mod1 <- update(Mod1)
#'
#' res <- fa_model_rr(Mod1, name.env = "env", name.gen = "gen")
#' res$diagnostics
#' head(res$blups)
#' }
#'
#' @export
fa_model_rr = function(model, name.env, name.gen){
  
  require(asreml)
  
  data = as.data.frame(model$mf)
  num.env = nlevels(data[, name.env])
  num.gen = nlevels(data[, name.gen])
  env.lev = levels(data[, name.env])
  gen.lev = levels(data[, name.gen])
  
  sum = summary(model)
  load = sum$varcomp[grep('fa\\d+', rownames(sum$varcomp)), , drop = FALSE]
  var = sum$varcomp[grep('!var$', rownames(sum$varcomp)), 1]
  # rr(env, k):vm(gen, G) + diag(env):gen model -> the specific variances
  # are in the diag(env):gen term (in the rr term they are fixed at zero)
  esp = sum$varcomp[grep(paste0("^", name.env, ":", name.gen, "!", name.env, "_"),
                         rownames(sum$varcomp)), 1]
  if (length(esp) == num.env) var = var + esp
  load$fa = regmatches(rownames(load), regexpr("fa\\d+", rownames(load)))
  mat.loadings = do.call(cbind, lapply(split(load[, c('fa','component')], f = load$fa), 
                                       function(x) x[, 'component']))
  rownames(mat.loadings) = env.lev
  
  svd.lambda = svd(mat.loadings)
  D = diag(svd.lambda$d^2, nrow = length(svd.lambda$d))
  
  if(sum(svd.lambda$u[,1] < 0) / nrow(svd.lambda$u) > 0.5){
    svd.lambda$u = -1 * svd.lambda$u
    svd.lambda$v = -1 * svd.lambda$v
  }
  
  mat.loadings.star = svd.lambda$u
  colnames(mat.loadings.star) = colnames(mat.loadings)
  rownames(mat.loadings.star) = rownames(mat.loadings)
  
  coef.rand = data.frame(summary(model, coef = TRUE)$coef.random)
  scores = coef.rand[grep('Comp', rownames(coef.rand)), , drop = FALSE]
  scores$fa = regmatches(rownames(scores), regexpr("Comp\\d+", rownames(scores)))
  
  scor.vec = do.call(c, lapply(split(scores, f = scores$fa), function(x) x[, 'solution']))
  # same as kronecker(sqrt(D) V', I) %*% scores, without building the huge matrix
  scor.vec.star = as.vector(matrix(scor.vec, nrow = num.gen) %*% svd.lambda$v %*% sqrt(D))
  
  scor.mat.star = matrix(
    scor.vec.star,
    nrow = num.gen,
    ncol = length(unique(scores$fa)),
    dimnames = list(gen.lev, unique(load$fa))
  )
  
  Gvcov = mat.loadings.star %*% tcrossprod(D, mat.loadings.star) + diag(var)
  Gcor = cov2cor(Gvcov)
  
  expvar = sum(diag(mat.loadings.star %*% tcrossprod(D, mat.loadings.star))) / sum(diag(Gvcov))
  
  expvar.j = matrix(nrow = ncol(mat.loadings),
                    ncol = num.env, 
                    dimnames = list(colnames(mat.loadings), env.lev))
  
  for(i in 1:nrow(expvar.j)){
    for(j in 1:ncol(expvar.j)){
      expvar.j[i, j] = 100 * mat.loadings.star[j, i]^2 * diag(D)[i] /
        (sum(mat.loadings.star[j, ]^2 * diag(D)) + var[j])
    }
  }
  
  lambdacross = mat.loadings.star %*% tcrossprod(D, mat.loadings.star)
  
  fullsv = list()
  for(i in env.lev){
    semivar = matrix(nrow = num.env, ncol = 3,
                     dimnames = list(env.lev, c('i', 'j', 'semivar')))
    for(j in env.lev){
      semivar[j, ] = c(i, j, .5 * (lambdacross[i, i] + lambdacross[j, j]) - lambdacross[i, j])
    }
    fullsv[[i]] = semivar
  }
  
  fullsv = data.frame(do.call(rbind, fullsv))
  fullsv$semivar = as.numeric(fullsv$semivar)
  fullsv = reshape(data = fullsv, timevar = 'i', idvar = 'j', direction = 'wide')[, -1]
  colnames(fullsv) = sub('semivar.', '', colnames(fullsv))
  lambda_ASV = 2 / (nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])
  
  fullsv = list()
  for(i in env.lev){
    semivar = matrix(nrow = num.env, ncol = 3,
                     dimnames = list(env.lev, c('i', 'j', 'semivar')))
    for(j in env.lev){
      semivar[j, ] = c(i, j, 0.5 * (Gvcov[i, i] + Gvcov[j, j]) - Gvcov[i, j])
    }
    fullsv[[i]] = semivar
  }
  
  fullsv = data.frame(do.call(rbind, fullsv))
  fullsv$semivar = as.numeric(fullsv$semivar)
  fullsv = reshape(data = fullsv, timevar = 'i', idvar = 'j', direction = 'wide')[, -1]
  colnames(fullsv) = sub('semivar.', '', colnames(fullsv))
  G_ASV = 2 / (nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])
  ASVR = lambda_ASV / G_ASV
  
  marginal = data.frame(
    env = rep(env.lev, each = num.gen),
    gen = rep(gen.lev, times = num.env),
    marginal = as.vector(matrix(scor.vec.star, nrow = num.gen) %*% t(mat.loadings.star))   # = kronecker(loadings, I) %*% scores
  )
  
  cf = data.frame(summary(model, coef = TRUE)$coef.fixed)
  cf$term = rownames(cf)
  rownames(cf) = NULL
  
  cr = data.frame(summary(model, coef = TRUE)$coef.random)
  cr$term = rownames(cr)
  rownames(cr) = NULL
  
  intercept = cf$solution[cf$term == "(Intercept)"]
  if(length(intercept) == 0) intercept = 0
  
  env_fixed = data.frame(env = env.lev, env_effect = 0)
  # environment fixed effect only (not env:covariate terms)
  tmp_env = cf[grep(paste0("^", name.env, "_[^:]+$"), cf$term), , drop = FALSE]
  if(nrow(tmp_env) > 0){
    tmp_env$env = sub(paste0("^", name.env, "_"), "", tmp_env$term)
    env_fixed$env_effect[match(tmp_env$env, env_fixed$env)] = tmp_env$solution
  }
  # genotype x environment random effects, whatever the term:
  #   gen:fa(env, k), fa(env, k):vm(gen, G), rr(env, k):gen or rr(env, k):vm(gen, G)
  #   (common part) and diag(env):gen (specific part, summed to the common one)
  pe = paste0("(?:(?:fa|rr)\\(", name.env, ", \\d+\\)|", name.env, ")_")
  pg = paste0("(?:vm\\(", name.gen, ", [^)]+\\)|", name.gen, ")_")
  env_1 = paste0("^", pe, "(.+?):", pg, "(.+)$")    # environment first
  gen_1 = paste0("^", pg, "(.+?):", pe, "(.+)$")    # genotype first
  a = grepl(env_1, cr$term, perl = TRUE)
  b = grepl(gen_1, cr$term, perl = TRUE) & !a
  rand_ge = rbind(
    data.frame(env = sub(env_1, "\\1", cr$term[a], perl = TRUE),
               gen = sub(env_1, "\\2", cr$term[a], perl = TRUE),
               solution = cr$solution[a]),
    data.frame(env = sub(gen_1, "\\2", cr$term[b], perl = TRUE),
               gen = sub(gen_1, "\\1", cr$term[b], perl = TRUE),
               solution = cr$solution[b])
  )
  rand_ge = rand_ge[rand_ge$env %in% env.lev & rand_ge$gen %in% gen.lev, ]   # drops the scores (Comp)
  if(nrow(rand_ge) > 0) rand_ge = aggregate(solution ~ gen + env, data = rand_ge, FUN = base::sum)   # "sum" is the model summary inside this function
  else rand_ge = data.frame(gen = character(0), env = character(0), solution = numeric(0))
  conditional_df = expand.grid(
    env = env.lev,
    gen = gen.lev,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  conditional_df = merge(conditional_df, env_fixed, by = "env", all.x = TRUE)
  conditional_df = merge(
    conditional_df,
    rand_ge[, c("gen", "env", "solution")],
    by = c("gen", "env"),
    all.x = TRUE
  )
  
  conditional_df$env_effect[is.na(conditional_df$env_effect)] = 0
  conditional_df$solution[is.na(conditional_df$solution)] = 0
  conditional_df$conditional = intercept + conditional_df$env_effect + conditional_df$solution
  conditional_df = conditional_df[, c("env", "gen", "conditional")]
  
  blups = merge(conditional_df, marginal, by = c("env", "gen"))
  
  results = list(
    rot.loads   = mat.loadings.star,
    rot.scores  = scor.mat.star,
    Gvcov       = Gvcov,
    Gcor        = Gcor,
    diagnostics = c(
      expvar = round(expvar * 100, 3),
      ASVR   = round(ASVR * 100, 3),
      aic    = round(sum$aic, 3),
      bic    = round(sum$bic, 3),
      logl   = round(sum$loglik, 3)
    ),
    expvar_j    = expvar.j,
    marginal    = marginal,
    blups       = blups
  )
  
  return(results)
}
