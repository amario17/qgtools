#' =============================================================================
#' get_blues() — BLUEs de genótipos por ambiente (análise em dois estágios)
#' =============================================================================
#'
#' Autor:  João Marcos Amario
#' E-mail: joao.m.sousa@ufv.br
#'
#' -----------------------------------------------------------------------------
#' CONTEXTO: POR QUE ESTIMAR BLUEs POR AMBIENTE?
#' -----------------------------------------------------------------------------
#' Em ensaios multiambientais (MET), a análise costuma ser feita em dois estágios:
#'
#'   1º ESTÁGIO (esta função)
#'      Cada ambiente é analisado separadamente. Ajusta-se um modelo que remove
#'      os efeitos de delineamento (blocos, repetições) e retorna uma média
#'      ajustada por genótipo naquele ambiente — o BLUE.
#'
#'   2º ESTÁGIO
#'      A matriz de BLUEs (genótipos x ambientes) alimenta o modelo GxE / genômico,
#'      onde entram a matriz de parentesco (G) e as covariáveis ambientais (W).
#'
#' A vantagem do dois estágios é computacional: em vez de ajustar um modelo único
#' e gigantesco com todos os dados de parcela, resolve-se um problema pequeno por
#' ambiente e passa-se adiante apenas o resumo. O custo é uma perda de informação
#' — as incertezas do 1º estágio não se propagam automaticamente para o 2º. Por
#' isso a função devolve também os erros-padrão, que podem ser usados como pesos
#' no 2º estágio (ver `weights` na saída).
#'
#' -----------------------------------------------------------------------------
#' BLUE vs BLUP: POR QUE GENÓTIPO É FIXO AQUI?
#' -----------------------------------------------------------------------------
#' BLUE  (Best Linear Unbiased Estimator) — genótipo como efeito FIXO.
#'       Não há encolhimento: a média de cada genótipo é estimada por si só.
#'       É o que se quer no 1º estágio, justamente para não aplicar shrinkage
#'       duas vezes (uma aqui, outra no modelo genômico do 2º estágio).
#'
#' BLUP  (Best Linear Unbiased Predictor) — genótipo como efeito ALEATÓRIO.
#'       Há encolhimento em direção à média geral, proporcional à herdabilidade.
#'       Apropriado quando se quer o valor genético em si, não um insumo para
#'       um segundo modelo.
#'
#' Como o objetivo aqui é montar a matriz de fenótipos ajustados que entrará no
#' modelo de predição genômica, o default é genótipo FIXO (BLUE). A função
#' permite `genotype_random = TRUE` para quem quiser BLUPs e a herdabilidade
#' por ambiente.
#'
#' -----------------------------------------------------------------------------
#' O MODELO AJUSTADO
#' -----------------------------------------------------------------------------
#' Por padrão, dentro de cada ambiente:
#'
#'     y_ij = mu + g_i + b_j + e_ij
#'
#'     g_i   efeito do genótipo i        (FIXO por default)
#'     b_j   efeito do bloco j           (ALEATÓRIO por default)
#'     e_ij  resíduo
#'
#' Bloco como aleatório é o default porque, com poucos blocos, tratá-lo como fixo
#' consome graus de liberdade sem ganho e impede a recuperação de informação
#' interbloco. Use `block_random = FALSE` para blocos fixos (o modelo original).
#'
#' Ambientes com um único bloco têm o termo de bloco removido automaticamente.
#'
#' -----------------------------------------------------------------------------
#' DEPENDÊNCIAS
#' -----------------------------------------------------------------------------
#'   ASReml-R (licença comercial) + dplyr
#'   install.packages("dplyr")
#'
#' =============================================================================
#'
#' -----------------------------------------------------------------------------
#' EXEMPLOS DE USO
#' -----------------------------------------------------------------------------
#'
#' library(dplyr)
#' library(asreml)
#' source("get_blues.R")
#'
#' dados <- read.csv("data_M4.txt", header = TRUE, sep = "\t")
#'
#' ## --- Uso básico: BLUEs de GY, bloco aleatório -----------------------------
#'
#' res <- get_blues(dados, trait = "GY")
#'
#' res$blues        # env, gen, blue, se
#' res$wide         # matriz genótipos x ambientes
#' res$summary      # diagnóstico por ambiente (n, convergência, CV%)
#' res$failed       # ambientes que não convergiram
#'
#' write.csv(res$blues, "blues_M4.csv", row.names = FALSE)
#'
#'
#' ## --- Bloco fixo (equivalente ao modelo original) --------------------------
#'
#' res <- get_blues(dados, trait = "GY", block_random = FALSE)
#'
#'
#' ## --- BLUPs + herdabilidade por ambiente -----------------------------------
#'
#' res <- get_blues(dados, trait = "GY", genotype_random = TRUE)
#' res$summary$h2      # herdabilidade média (Cullis) por ambiente
#'
#'
#' ## --- Delineamento com repetição dentro de bloco ---------------------------
#'
#' res <- get_blues(dados, trait = "GY",
#'                  block = "rep",
#'                  extra_random = "~ rep:block")
#'
#'
#' ## --- Loop nas cinco estruturas --------------------------------------------
#'
#' for (s in paste0("M", 1:5)) {
#'   d <- read.csv(sprintf("data_%s.txt", s), header = TRUE, sep = "\t")
#'   r <- get_blues(d, trait = "GY")
#'   write.csv(r$blues, sprintf("blues_%s.csv", s), row.names = FALSE)
#' }
#'
#' =============================================================================


#' BLUEs (ou BLUPs) de genótipos, ambiente a ambiente
#'
#' @param data data.frame. Dados de parcela, formato longo (uma linha por parcela).
#' @param trait character. Nome da coluna da variável resposta. Ex.: "GY".
#' @param gen character. Coluna do genótipo. Default "gen".
#' @param env character. Coluna do ambiente. Default "env".
#' @param block character ou NULL. Coluna do bloco/repetição. Default "block".
#'   Use NULL para ajustar apenas a média (delineamento inteiramente casualizado
#'   sem estrutura).
#' @param block_random logical. Bloco como efeito aleatório. Default TRUE.
#' @param genotype_random logical. Genótipo como efeito aleatório — devolve BLUPs
#'   e herdabilidade em vez de BLUEs. Default FALSE.
#' @param extra_fixed character ou NULL. Termos fixos adicionais, em notação de
#'   fórmula. Ex.: "+ cov1 + cov2".
#' @param extra_random character ou NULL. Termos aleatórios adicionais.
#'   Ex.: "~ rep:block".
#' @param min_gen integer. Mínimo de genótipos para o ambiente ser analisado.
#'   Default 3.
#' @param maxit integer. Máximo de iterações do ASReml. Default 100.
#' @param workspace numeric. Memória do ASReml. Default 6e8.
#' @param verbose logical. Imprimir progresso. Default TRUE.
#'
#' @return Uma lista com:
#'   \item{blues}{data.frame longo: env, gen, blue, se}
#'   \item{wide}{matriz genótipos (linhas) x ambientes (colunas)}
#'   \item{weights}{data.frame com pesos 1/se^2, para uso no 2º estágio}
#'   \item{summary}{diagnóstico por ambiente: n, n_gen, n_block, convergência,
#'     média, CV\%, e h2 se genotype_random = TRUE}
#'   \item{failed}{ambientes que falharam ou não convergiram}
#'   \item{models}{lista dos objetos asreml, um por ambiente}
#'
#' @export
get_blues <- function(data,
                      trait,
                      gen             = "gen",
                      env             = "env",
                      block           = "block",
                      block_random    = TRUE,
                      genotype_random = FALSE,
                      extra_fixed     = NULL,
                      extra_random    = NULL,
                      min_gen         = 3,
                      maxit           = 100,
                      workspace       = 6e8,
                      verbose         = TRUE) {
  
  say <- function(...) if (verbose) cat(...)
  
  # ==========================================================================
  # VALIDAÇÃO
  # ==========================================================================
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop("Pacote 'asreml' é necessário (licença comercial).")
  }
  if (!is.data.frame(data)) stop("`data` precisa ser um data.frame.")
  if (missing(trait))       stop("Informe `trait` — o nome da coluna da variável resposta.")
  
  obrigatorias <- c(trait, gen, env)
  faltando     <- setdiff(obrigatorias, names(data))
  if (length(faltando) > 0) {
    stop("Coluna(s) ausente(s) em `data`: ", paste(faltando, collapse = ", "),
         "\nDisponíveis: ", paste(names(data), collapse = ", "))
  }
  
  tem_bloco <- !is.null(block)
  if (tem_bloco && !block %in% names(data)) {
    stop("Coluna de bloco '", block, "' não encontrada. Use block = NULL se não houver.")
  }
  
  # ---- Tipagem: fatores nos efeitos, numérico na resposta ------------------
  data[[gen]] <- factor(as.character(data[[gen]]))
  data[[env]] <- factor(as.character(data[[env]]))
  if (tem_bloco) data[[block]] <- factor(as.character(data[[block]]))
  
  data[[trait]] <- suppressWarnings(as.numeric(as.character(data[[trait]])))
  
  n_na <- sum(is.na(data[[trait]]))
  if (n_na == nrow(data)) {
    stop("A coluna '", trait, "' virou toda NA ao converter para numérico. ",
         "Confira o separador decimal do arquivo.")
  }
  
  ambientes <- levels(data[[env]])
  
  say("=============================================================\n")
  say(" BLUEs por ambiente — ", trait, "\n")
  say("=============================================================\n")
  say("Ambientes  : ", length(ambientes), "\n")
  say("Genótipos  : ", nlevels(data[[gen]]), "\n")
  say("Parcelas   : ", nrow(data), "\n")
  say("Genótipo   : ", if (genotype_random) "ALEATÓRIO (BLUP)" else "FIXO (BLUE)", "\n")
  if (tem_bloco) say("Bloco      : ", if (block_random) "aleatório" else "fixo", "\n")
  if (n_na > 0)  say("NAs em ", trait, ": ", n_na, "\n")
  say("\n")
  
  resultados <- list()
  diagnostico <- list()
  modelos    <- list()
  falharam   <- character(0)
  
  # ==========================================================================
  # LOOP POR AMBIENTE
  # ==========================================================================
  for (e in ambientes) {
    
    d <- data[data[[env]] == e, , drop = FALSE]
    
    # Refatorar dentro do ambiente: níveis não presentes quebram o asreml
    d[[gen]] <- factor(as.character(d[[gen]]))
    if (tem_bloco) d[[block]] <- factor(as.character(d[[block]]))
    
    n_gen   <- nlevels(d[[gen]])
    n_block <- if (tem_bloco) nlevels(d[[block]]) else 0L
    n_obs   <- sum(!is.na(d[[trait]]))
    
    # ---- Pular ambientes inviáveis ----------------------------------------
    if (n_obs == 0) {
      say(sprintf("  %s — sem observações válidas, pulado\n", e))
      falharam <- c(falharam, e); next
    }
    if (n_gen < min_gen) {
      say(sprintf("  %s — apenas %d genótipos (< %d), pulado\n", e, n_gen, min_gen))
      falharam <- c(falharam, e); next
    }
    
    # Bloco só entra se houver mais de um nível
    usar_bloco <- tem_bloco && n_block > 1
    
    # ---- Montar as fórmulas ------------------------------------------------
    termos_fixos <- character(0)
    termos_aleat <- character(0)
    
    if (genotype_random) termos_aleat <- c(termos_aleat, gen)
    else                 termos_fixos <- c(termos_fixos, gen)
    
    if (usar_bloco) {
      if (block_random) termos_aleat <- c(termos_aleat, block)
      else              termos_fixos <- c(termos_fixos, block)
    }
    
    f_fixed <- paste(trait, "~",
                     if (length(termos_fixos) > 0) paste(termos_fixos, collapse = " + ") else "1")
    if (!is.null(extra_fixed)) f_fixed <- paste(f_fixed, extra_fixed)
    f_fixed <- stats::as.formula(f_fixed)
    
    f_random <- NULL
    if (length(termos_aleat) > 0) {
      f_random <- paste("~", paste(termos_aleat, collapse = " + "))
      if (!is.null(extra_random)) {
        f_random <- paste(f_random, "+", sub("^\\s*~\\s*", "", extra_random))
      }
      f_random <- stats::as.formula(f_random)
    } else if (!is.null(extra_random)) {
      f_random <- stats::as.formula(extra_random)
    }
    
    # ---- Ajustar ------------------------------------------------------------
    args <- list(fixed     = f_fixed,
                 data      = d,
                 maxit     = maxit,
                 workspace = workspace,
                 trace     = FALSE,
                 na.action = asreml::na.method(x = "include", y = "include"))
    if (!is.null(f_random)) args$random <- f_random
    
    mod <- try(suppressWarnings(do.call(asreml::asreml, args)), silent = TRUE)
    
    if (inherits(mod, "try-error")) {
      say(sprintf("  %s — FALHOU: %s\n", e,
                  sub("\n.*", "", conditionMessage(attr(mod, "condition")))))
      falharam <- c(falharam, e); next
    }
    
    # Uma rodada extra de update() ajuda a convergir
    if (!isTRUE(mod$converge)) {
      mod <- try(suppressWarnings(update(mod)), silent = TRUE)
      if (inherits(mod, "try-error")) {
        say(sprintf("  %s — FALHOU no update()\n", e))
        falharam <- c(falharam, e); next
      }
    }
    
    convergiu <- isTRUE(mod$converge)
    
    # ---- Extrair as médias ajustadas ---------------------------------------
    pv <- try(
      suppressWarnings(predict(mod, classify = gen, trace = FALSE)$pvals),
      silent = TRUE)
    
    if (inherits(pv, "try-error") || is.null(pv) || nrow(pv) == 0) {
      say(sprintf("  %s — FALHOU no predict()\n", e))
      falharam <- c(falharam, e); next
    }
    
    pv <- as.data.frame(pv)
    
    res_e <- data.frame(
      env  = e,
      gen  = as.character(pv[[gen]]),
      blue = pv$predicted.value,
      se   = pv$std.error,
      stringsAsFactors = FALSE)
    
    resultados[[e]] <- res_e
    modelos[[e]]    <- mod
    
    # ---- Diagnóstico --------------------------------------------------------
    media <- mean(d[[trait]], na.rm = TRUE)
    vres  <- summary(mod)$varcomp
    sigma2 <- if ("units!R" %in% rownames(vres)) vres["units!R", "component"] else
      utils::tail(vres$component, 1)
    cv <- 100 * sqrt(sigma2) / media
    
    # Herdabilidade de Cullis (só faz sentido com genótipo aleatório)
    h2 <- NA_real_
    if (genotype_random) {
      vg  <- vres[grep(paste0("^", gen), rownames(vres))[1], "component"]
      vbl <- mean(pv$std.error^2, na.rm = TRUE)   # PEV média
      h2  <- 1 - vbl / (2 * vg)
    }
    
    diagnostico[[e]] <- data.frame(
      env       = e,
      n_obs     = n_obs,
      n_gen     = n_gen,
      n_block   = n_block,
      converged = convergiu,
      mean      = media,
      cv        = cv,
      h2        = h2,
      stringsAsFactors = FALSE)
    
    say(sprintf("  %s — n=%d | gen=%d | CV=%.1f%%%s%s\n",
                e, n_obs, n_gen, cv,
                if (genotype_random) sprintf(" | h2=%.2f", h2) else "",
                if (convergiu) "" else "  [NÃO CONVERGIU]"))
  }
  
  if (length(resultados) == 0) {
    stop("Nenhum ambiente foi analisado com sucesso. Verifique os dados.")
  }
  
  # ==========================================================================
  # MONTAR A SAÍDA
  # ==========================================================================
  blues <- do.call(rbind, resultados)
  rownames(blues) <- NULL
  
  diag_df <- do.call(rbind, diagnostico)
  rownames(diag_df) <- NULL
  if (!genotype_random) diag_df$h2 <- NULL
  
  # Matriz genótipos x ambientes
  wide <- stats::reshape(blues[, c("env", "gen", "blue")],
                         idvar = "gen", timevar = "env", direction = "wide")
  names(wide) <- sub("^blue\\.", "", names(wide))
  rownames(wide) <- NULL
  
  # Pesos para o 2º estágio: inverso da variância do BLUE.
  # Ambientes/genótipos estimados com mais precisão pesam mais no modelo GxE.
  pesos <- blues
  pesos$weight <- ifelse(is.na(pesos$se) | pesos$se <= 0, NA_real_, 1 / pesos$se^2)
  pesos <- pesos[, c("env", "gen", "weight")]
  
  say("\n=============================================================\n")
  say(" CONCLUÍDO\n")
  say("=============================================================\n")
  say("Ambientes analisados : ", length(resultados), " de ", length(ambientes), "\n")
  say("BLUEs estimados      : ", nrow(blues), "\n")
  say("Matriz gen x env     : ", nrow(wide), " x ", ncol(wide) - 1, "\n")
  
  nao_conv <- diag_df$env[!diag_df$converged]
  if (length(nao_conv) > 0) {
    say("\nNão convergiram (", length(nao_conv), "): ",
        paste(nao_conv, collapse = ", "), "\n")
    say("Resultados devolvidos mesmo assim — inspecione antes de usar.\n")
  }
  if (length(falharam) > 0) {
    say("\nFalharam (", length(falharam), "): ", paste(falharam, collapse = ", "), "\n")
  }
  
  list(blues   = blues,
       wide    = wide,
       weights = pesos,
       summary = diag_df,
       failed  = falharam,
       models  = modelos)
}