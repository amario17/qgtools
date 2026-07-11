#' Calcula BLUEs e pesos (W) por ambiente, via ASReml
#'
#' Ajusta um modelo misto independente para cada nivel de um fator ambiental
#' (ex: env, local, ano), extrai as BLUEs de genotipo e a matriz de pesos W
#' (inverso da matriz de variancia-covariancia das preditas), usada
#' tipicamente como ponderacao em analises de segundo estagio (ex: modelos
#' de estabilidade, GIS-FA, ou qualquer pipeline que utilize BLUEs
#' ponderadas por ambiente).
#'
#' Suporta dois delineamentos experimentais:
#' \itemize{
#'   \item \strong{DBC} (blocos completos casualizados): efeito aleatorio
#'     \code{~ block_var}.
#'   \item \strong{Latice} (blocos incompletos dentro de repeticao): efeito
#'     aleatorio \code{~ rep_var + rep_var:block_var} (repeticao e bloco
#'     dentro de repeticao, ambos aleatorios) - estrutura classica de
#'     analise de latice/alpha-latice por ambiente.
#' }
#'
#' Para cada nivel de \code{env_var}, o modelo ajustado eh
#' \code{response ~ gen_var} (fixo) + o efeito aleatorio correspondente ao
#' \code{design} escolhido, usando os dados restritos aquele ambiente. O
#' peso de cada genotipo naquele ambiente eh a diagonal da inversa da
#' matriz de variancia das preditas (\code{W = diag(solve(vcov))}), e
#' \code{sigma} eh a variancia residual (\code{sigma2}) do modelo daquele
#' ambiente.
#'
#' @param data Data frame com os dados fenotipicos.
#' @param response Nome (string) da variavel resposta (ex: "prod", "GY").
#' @param gen_var Nome (string) da coluna de genotipo (ex: "gen").
#' @param block_var Nome (string) da coluna de bloco. Para \code{design =
#'   "dbc"}, e o bloco completo casualizado. Para \code{design = "latice"},
#'   e o bloco incompleto dentro de repeticao.
#' @param env_var Nome (string) da coluna que define os ambientes, um
#'   modelo independente sera ajustado para cada nivel (ex: "env", "local").
#' @param design Delineamento experimental: \code{"dbc"} (blocos completos
#'   casualizados, padrao) ou \code{"latice"} (blocos incompletos dentro de
#'   repeticao).
#' @param rep_var Nome (string) da coluna de repeticao. OBRIGATORIO quando
#'   \code{design = "latice"}; ignorado quando \code{design = "dbc"}.
#' @param maxit Numero maximo de iteracoes do ASReml (padrao: 100).
#' @param workspace Memoria de workspace do ASReml (padrao: "4gb").
#' @param pworkspace Memoria de pworkspace do ASReml (padrao: "2gb").
#'
#' @return Um data frame (\code{yhat}) com uma linha por genotipo x ambiente,
#'   contendo pelo menos: \code{gen_var}, \code{env_var}, \code{predicted.value}
#'   (BLUE), \code{std.error}, \code{status}, \code{W} (peso) e \code{sigma}
#'   (variancia residual do ambiente).
#'
#' @examples
#' \dontrun{
#' # Delineamento DBC (blocos completos casualizados)
#' yhat <- get_weights_blue(
#'   data      = data,
#'   response  = "GY",
#'   gen_var   = "gen",
#'   block_var = "block",
#'   env_var   = "env",
#'   design    = "dbc"
#' )
#'
#' # Delineamento Latice (blocos incompletos dentro de repeticao)
#' yhat <- get_weights_blue(
#'   data      = data,
#'   response  = "GY",
#'   gen_var   = "gen",
#'   block_var = "bloco_incompleto",
#'   env_var   = "env",
#'   design    = "latice",
#'   rep_var   = "rep"
#' )
#' }
#'
#' @references
#' Smith, A.B., Cullis, B.R., Thompson, R. (2001). Analyzing variety by
#' environment data using multiplicative mixed models and adjustments for
#' spatial field trend. \emph{Biometrics}, 57(4), 1138-1147.
#'
#' Smith, A., Cullis, B., Gilmour, A. (2001). The analysis of crop variety
#' evaluation data in Australia. \emph{Australian & New Zealand Journal of
#' Statistics}, 43(2), 129-145.
#'
#' Möhring, J., Piepho, H.P. (2009). Comparison of Weighting in Two-Stage
#' Analysis of Plant Breeding Trials. \emph{Crop Science}, 49(6), 1977-1988.
#'
#' @author João Marcos Amario, PhD (\email{joao.m.sousa@ufv.br})
#'
#' @export
get_weights_blue <- function(data,
                             response,
                             gen_var,
                             block_var,
                             env_var,
                             design     = c("dbc", "latice"),
                             rep_var    = NULL,
                             maxit      = 100,
                             workspace  = "4gb",
                             pworkspace = "2gb") {
  
  require(asreml)
  require(dplyr)
  require(tibble)
  
  design <- match.arg(design)
  
  if (design == "latice" && is.null(rep_var)) {
    stop("Para design = 'latice', informe 'rep_var' (coluna de repeticao).")
  }
  
  asreml.options(maxit = maxit, workspace = workspace, pworkspace = pworkspace)
  
  eda <- droplevels(data)
  
  fixed_formula <- as.formula(paste(response, "~", gen_var))
  
  random_formula <- switch(
    design,
    dbc    = as.formula(paste("~", block_var)),
    latice = as.formula(paste("~", rep_var, "+", rep_var, ":", block_var))
  )
  
  pred_env  <- list()
  sigma_env <- NULL
  
  for (i in levels(eda[[env_var]])) {
    
    mod <- asreml(
      fixed     = fixed_formula,
      random    = random_formula,
      data      = eda,
      na.action = na.method(x = "exclude", y = "exclude"),
      subset    = eda[[env_var]] == i
    )
    
    sigma_env[i] <- mod$sigma2
    
    pred <- predict(mod, classify = gen_var, vcov = TRUE)
    pred$pvals$W          <- diag(solve(pred$vcov))
    pred$pvals[[env_var]] <- i
    
    pred_env[[i]] <- pred$pvals
  }
  
  yhat <- do.call(rbind, pred_env)
  
  sigma_df <- data.frame(sigma = sigma_env) %>%
    rownames_to_column(env_var)
  
  yhat <- yhat %>%
    left_join(
      eda[, c(gen_var, env_var)] %>% distinct(),
      by = c(gen_var, env_var)
    ) %>%
    left_join(sigma_df, by = env_var)
  
  yhat
}