#' Analise individual por ambiente (Single Analysis)
#'
#' Ajusta, para cada nivel de um fator ambiental, um modelo misto com efeito
#' aleatorio de genotipo e um modelo nulo (sem genotipo), via ASReml. A
#' partir dai calcula: componentes de variancia genetica e residual, teste
#' de razao de verossimilhanca (LRT) para significancia do efeito de
#' genotipo, herdabilidade no sentido de Cullis (\emph{Cullis heritability},
#' baseada no erro padrao medio da diferenca - avsed), herdabilidade no
#' sentido amplo (\emph{broad-sense heritability}) e coeficiente de
#' variacao (CV) residual, todos por ambiente.
#'
#' Ao final, os ambientes cujo efeito de genotipo NAO foi significativo
#' (LRT com p-valor >= \code{alpha}) sao removidos, e os ambientes
#' remanescentes sao renomeados sequencialmente (\code{E01}, \code{E02}, ...)
#' na ordem em que aparecem, tanto na tabela resumo quanto nos dados
#' filtrados - uteis para etapas posteriores (ex: modelo FA multiambiente)
#' que dependem de nomes de ambiente sequenciais e limpos.
#'
#' @param data Data frame com os dados fenotipicos.
#' @param response Nome (string) da variavel resposta (ex: "GY").
#' @param gen_var Nome (string) da coluna de genotipo (ex: "gen").
#' @param block_var Nome (string) da coluna de bloco, usada como efeito
#'   fixo no modelo individual (ex: "block").
#' @param env_var Nome (string) da coluna que define os ambientes, uma
#'   analise independente sera feita para cada nivel (ex: "env").
#' @param alpha Nivel de significancia do teste LRT para reter um ambiente
#'   (padrao: 0.05). Ambientes com p-valor >= alpha sao descartados.
#'
#' @return Uma lista com:
#' \describe{
#'   \item{summary_full}{Data frame com sig2g, sig2e, H2_Cullis, H2_Broad,
#'     CV, Pval e Signif, para TODOS os ambientes (antes do filtro).}
#'   \item{summary_significant}{O mesmo resumo, mas so para os ambientes
#'     significativos, com o ambiente ja renomeado (E01, E02, ...).}
#'   \item{data_significant}{O \code{data} original, filtrado para conter
#'     apenas os ambientes significativos, com \code{env_var} renomeado
#'     (E01, E02, ...) e niveis nao usados removidos (droplevels).}
#'   \item{env_map}{Vetor nomeado com o mapeamento do nome original do
#'     ambiente para o nome renomeado (E01, E02, ...).}
#' }
#'
#' @references
#' Cullis, B.R., Smith, A.B., Coombes, N.E. (2006). On the design of early
#' generation variety trials with correlated data. \emph{Journal of
#' Agricultural, Biological, and Environmental Statistics}, 11(4), 381-393.
#'
#' Piepho, H.P., Möhring, J. (2007). Computing Heritability and Selection
#' Response From Unbalanced Plant Breeding Trials. \emph{Genetics}, 177(3),
#' 1881-1888.
#'
#' @author João Marcos Amario, PhD (\email{joao.m.sousa@ufv.br})
#'
#' @examples
#' \dontrun{
#' res <- single_analysis(
#'   data      = data,
#'   response  = "GY",
#'   gen_var   = "gen",
#'   block_var = "block",
#'   env_var   = "env",
#'   alpha     = 0.05
#' )
#' str(res$data_significant)
#' res$summary_significant
#' }
#'
#' @export
single_analysis <- function(data, response, gen_var, block_var, env_var, alpha = 0.05) {
  
  require(asreml)
  require(dplyr)
  require(tidyr)
  require(tibble)
  
  fixed_formula  <- as.formula(paste(response, "~", block_var))
  random_formula <- as.formula(paste("~", gen_var))
  
  suma    <- list()
  signi   <- list()
  H2c     <- list()
  H2broad <- list()
  cv      <- list()
  
  for (i in levels(data[[env_var]])) {
    
    m  <- asreml(fixed = fixed_formula, random = random_formula,
                 data = data, subset = data[[env_var]] == i)
    m0 <- asreml(fixed = fixed_formula,
                 data = data, subset = data[[env_var]] == i)
    
    vc <- summary(m)$varcomp
    suma[[i]] <- vc
    
    sig2g <- vc[gen_var, "component"]
    sig2e <- vc["units!R", "component"]
    
    lrt_res <- lrt(m, m0)
    signi[[i]] <- lrt_res
    
    pred   <- predict(m, classify = gen_var, vcov = TRUE)
    avsed  <- pred$avsed[2]
    
    H2c[[i]] <- round(1 - (avsed^2) / (2 * sig2g), 2)
    
    nrep <- length(unique(data[[block_var]][data[[env_var]] == i]))
    H2broad[[i]] <- round(sig2g / (sig2g + (sig2e / nrep)), 2)
    
    media  <- mean(data[[response]][data[[env_var]] == i], na.rm = TRUE)
    cv[[i]] <- (sqrt(sig2e) / media) * 100
  }
  
  suma <- do.call(rbind, suma) |>
    rownames_to_column("rowname") |>
    dplyr::select(rowname, component) |>
    separate(col = rowname, into = c("Amb", "Efeito")) |>
    reshape(timevar = "Efeito", idvar = "Amb", direction = "wide")
  colnames(suma) <- c("Amb", "sig2g", "sig2e")
  
  signi <- do.call(rbind, signi)
  signi$Amb   <- rownames(signi)
  signi$Pval  <- signi$`Pr(Chisq)`
  signi$Signif <- ifelse(signi$Pval < 0.001, "***",
                         ifelse(signi$Pval < 0.01, "**",
                                ifelse(signi$Pval < 0.05, "*", "ns")))
  signi <- signi[, c("Amb", "Pval", "Signif")]
  
  H2c     <- data.frame(Amb = names(H2c),     H2_Cullis = unlist(H2c))
  H2broad <- data.frame(Amb = names(H2broad), H2_Broad  = unlist(H2broad))
  cv      <- data.frame(Amb = names(cv),      CV        = unlist(cv))
  
  resultado_final <- Reduce(function(x, y) merge(x, y, by = "Amb"),
                            list(suma, H2c, H2broad, cv, signi))
  
  # ---- Filtra ambientes significativos (p < alpha) e renomeia sequencialmente ----
  amb_sig <- resultado_final$Amb[resultado_final$Pval < alpha]
  amb_sig <- amb_sig[order(as.integer(gsub("^[^0-9]*", "", amb_sig)))]
  
  data_sig <- droplevels(subset(data, data[[env_var]] %in% amb_sig))
  map_env  <- setNames(sprintf("E%02d", seq_along(amb_sig)), amb_sig)
  data_sig[[env_var]] <- factor(map_env[as.character(data_sig[[env_var]])],
                                levels = sprintf("E%02d", seq_along(amb_sig)))
  
  resultado_final_sig <- subset(resultado_final, Amb %in% amb_sig)
  resultado_final_sig$Amb <- map_env[resultado_final_sig$Amb]
  resultado_final_sig <- resultado_final_sig[order(as.integer(gsub("^E", "", resultado_final_sig$Amb))), ]
  row.names(resultado_final_sig) <- NULL
  
  list(
    summary_full        = resultado_final,
    summary_significant = resultado_final_sig,
    data_significant    = data_sig,
    env_map             = map_env
  )
}