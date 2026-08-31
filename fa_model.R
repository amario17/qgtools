#' Pos-processamento de um modelo misto fatorial (Factor Analytic - FA)
#'
#' A partir de um modelo misto ajustado no ASReml com termo aleatorio
#' \code{gen:fa(env, k)}, extrai e organiza os principais resultados do
#' modelo fatorial (FA): as cargas fatoriais (\emph{loadings}) e os
#' escores de genotipo (\emph{scores}) rotacionados para a forma canonica
#' via decomposicao em valores singulares (SVD), a matriz de covariancia
#' genetica entre ambientes (\code{Gvcov}) e a respectiva matriz de
#' correlacao genetica (\code{Gcor}) implicadas pela estrutura FA(k), alem
#' de medidas de qualidade de ajuste do modelo: variancia explicada total
#' e por fator/ambiente, \emph{Average Semivariance Ratio} (ASVR), AIC,
#' BIC e log-verossimilhanca.
#'
#' Tambem sao calculados os BLUPs de genotipo por ambiente em duas
#' versoes: o valor \strong{marginal} (\eqn{\Lambda f}), correspondente
#' apenas a parte comum da estrutura fatorial, sem o intercepto nem o
#' efeito fixo de ambiente - util como insumo para etapas de predicao em
#' ambientes nao testados (ex: GIS-FA), onde a media geral e demais
#' efeitos sao adicionados a parte; e o valor \strong{condicional}
#' (intercepto + efeito fixo de ambiente + desvio aleatorio especifico do
#' genotipo naquele ambiente), que representa a predicao completa do
#' genotipo dentro dos ambientes efetivamente testados. 
#'
#' @param model Objeto de modelo ajustado no ASReml (\code{asreml}), com
#'   um termo aleatorio fatorial \code{gen:fa(env, k)} para algum numero
#'   de fatores \code{k}.
#' @param name.env Nome (string) da coluna que identifica o ambiente no
#'   modelo (ex: "env").
#' @param name.gen Nome (string) da coluna que identifica o genotipo no
#'   modelo (ex: "gen").
#'
#' @return Uma lista com:
#' \describe{
#'   \item{rot.loads}{Matriz de cargas fatoriais rotacionadas (ambiente x
#'     fator), na forma canonica obtida via SVD.}
#'   \item{rot.scores}{Matriz de escores de genotipo rotacionados
#'     (genotipo x fator).}
#'   \item{Gvcov}{Matriz de covariancia genetica entre ambientes
#'     implicada pelo modelo FA(k) (cargas rotacionadas mais as
#'     variancias especificas de cada ambiente).}
#'   \item{Gcor}{Matriz de correlacao genetica entre ambientes,
#'     correspondente a \code{Gvcov} padronizada.}
#'   \item{diagnostics}{Vetor nomeado com \code{expvar} (percentual da
#'     variancia genetica total explicada pelos fatores comuns),
#'     \code{ASVR} (Average Semivariance Ratio, percentual), \code{aic},
#'     \code{bic} e \code{logl} (log-verossimilhanca) do modelo.}
#'   \item{expvar_j}{Matriz (fator x ambiente) com o percentual da
#'     variancia explicado por cada fator dentro de cada ambiente.}
#'   \item{marginal}{Data frame (env, gen, marginal) com a parte comum
#'     (\eqn{\Lambda f}) da predicao genotipica, sem intercepto nem
#'     efeito de ambiente, por ambiente e genotipo.}
#'   \item{blups}{Data frame (env, gen, conditional, marginal) unindo a
#'     predicao condicional (completa, dentro dos ambientes testados) e
#'     a marginal (\eqn{\Lambda f}), por ambiente e genotipo.}
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
#' Esta funcao foi adaptada e modificada a partir da funcao \code{fa.outs}
#' do Prof. Saulo Chaves da ESALQ/USP.
#' #############################################################################
#' 
#'
#' @author João Marcos Amario, PhD (\email{joao.m.sousa@ufv.br})
#'
#' @examples
#' \dontrun{
#' Mod1 <- asreml(
#'   fixed     = GY ~ env + block:env,
#'   random    = ~ gen:fa(env, 3),
#'   residual  = ~ dsum(~id(units) | env),
#'   data      = data_M1,
#'   na.action = na.method(x = "include", y = "include"),
#'   maxit     = 100
#' )
#' Mod1 <- update(Mod1)
#'
#' res <- fa_model(Mod1, name.env = "env", name.gen = "gen")
#' res$diagnostics
#' head(res$blups)
#' }
#'
#' @export
fa_model = function(model, name.env, name.gen){
  
  require(asreml)
  
  data = as.data.frame(model$mf)
  num.env = nlevels(data[, name.env])
  num.gen = nlevels(data[, name.gen])
  env.lev = levels(data[, name.env])
  gen.lev = levels(data[, name.gen])
  
  sum = summary(model)
  load = sum$varcomp[grep('fa\\d+', rownames(sum$varcomp)), , drop = FALSE]
  var = sum$varcomp[grep('var', rownames(sum$varcomp)), 1]
  
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
  scor.vec.star = kronecker(sqrt(D) %*% t(svd.lambda$v), diag(num.gen)) %*% scor.vec
  
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
    marginal = as.vector(kronecker(mat.loadings.star, diag(num.gen)) %*% scor.vec.star)
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
  tmp_env = cf[grep("^env_", cf$term), , drop = FALSE]
  
  if(nrow(tmp_env) > 0){
    tmp_env$env = sub("^env_", "", tmp_env$term)
    env_fixed$env_effect[match(tmp_env$env, env_fixed$env)] = tmp_env$solution
  }
  
  rand_ge = cr[grep("^gen_[^:]+:fa\\(env, [0-9]+\\)_E[0-9]+$", cr$term), , drop = FALSE]
  rand_ge$gen = sub("^gen_([^:]+):fa\\(env, [0-9]+\\)_E[0-9]+$", "\\1", rand_ge$term)
  rand_ge$env = sub("^gen_[^:]+:fa\\(env, [0-9]+\\)_(E[0-9]+)$", "\\1", rand_ge$term)
  
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