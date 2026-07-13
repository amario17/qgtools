#' =============================================================================
#' rfsi_envcov() — Interpolação espacial de covariáveis ambientais via RFSI
#' =============================================================================
#'
#' Autor:  João Marcos Amario
#' E-mail: joao.m.sousa@ufv.br
#'
#' -----------------------------------------------------------------------------
#' CRÉDITOS / MÉTODO ORIGINAL
#' -----------------------------------------------------------------------------
#' Esta função é um wrapper em torno do método Random Forest Spatial Interpolation
#' (RFSI), proposto por:
#'
#'   Sekulić, A., Kilibarda, M., Heuvelink, G.B.M., Nikolić, M., Bajat, B. (2020).
#'   Random Forest Spatial Interpolation. Remote Sensing, 12(10), 1687.
#'   https://doi.org/10.3390/rs12101687
#'
#' A implementação computacional utilizada aqui é a do pacote `meteo`
#' (funções `rfsi()` e `pred.rfsi()`), desenvolvido pelos mesmos autores.
#' Todo o crédito metodológico e de implementação do RFSI pertence a eles.
#' Esta função apenas organiza o fluxo de trabalho (tuning + predição +
#' montagem da matriz final) para o contexto de covariáveis ambientais em
#' melhoramento genético de plantas.
#'
#' -----------------------------------------------------------------------------
#' O QUE ESTA FUNÇÃO FAZ
#' -----------------------------------------------------------------------------
#' O problema: temos covariáveis CLIMÁTICAS disponíveis apenas em uma grade
#' espacial grosseira (ex.: 0.5°, ~50 km, limitada pela resolução do NASA POWER),
#' mas temos covariáveis de SOLO em uma grade fina (ex.: 0.01°, ~1 km, vinda do
#' SoilGrids). Para predizer o desempenho genotípico em ambientes não avaliados
#' na resolução fina, precisamos das duas fontes na MESMA grade.
#'
#' A solução: usar o RFSI para interpolar espacialmente cada variável climática
#' da grade grosseira para a grade fina, e então juntar tudo em uma única matriz.
#'
#' Diferente de um Random Forest comum, o RFSI não usa apenas as coordenadas
#' como preditores. Ele constrói, para cada ponto, um conjunto de covariáveis
#' espaciais formadas pelos VALORES OBSERVADOS nos n vizinhos mais próximos e
#' pelas DISTÂNCIAS até esses vizinhos. Isso permite que o modelo capture a
#' autocorrelação espacial explicitamente, aproximando-se do comportamento de um
#' krigging, mas sem assumir estacionariedade ou uma forma paramétrica de
#' variograma.
#'
#' -----------------------------------------------------------------------------
#' FLUXO DE TRABALHO EM DUAS ETAPAS (executado para CADA variável climática)
#' -----------------------------------------------------------------------------
#'
#'   ETAPA 1 — TUNING via validação cruzada leave-one-out (LOO-CV)
#'   -------------------------------------------------------------
#'   Percorre uma grade de hiperparâmetros (n.obs × num.trees). Para cada
#'   combinação, deixa um ponto de treino de fora por vez, treina o RFSI nos
#'   demais, prediz o ponto excluído e guarda o par (observado, predito).
#'   Ao final, calcula a correlação preditiva (cp) entre observado e predito.
#'   A combinação com maior cp é a vencedora.
#'
#'   Isso responde: "quantos vizinhos e quantas árvores este campo espacial
#'   precisa para ser bem reconstruído?" — e a resposta pode variar entre
#'   variáveis (temperatura é mais suave espacialmente que precipitação, por
#'   exemplo, e tende a precisar de menos vizinhos).
#'
#'   ETAPA 2 — PREDIÇÃO na grade fina
#'   ---------------------------------
#'   Retreina o modelo usando TODOS os pontos de treino com os hiperparâmetros
#'   vencedores, e prediz a variável em todos os pontos da grade fina.
#'
#' -----------------------------------------------------------------------------
#' PARÂMETROS DO RFSI — O QUE CADA UM SIGNIFICA
#' -----------------------------------------------------------------------------
#'
#'   n.obs (vizinhos)
#'       Número de observações vizinhas mais próximas usadas como preditores.
#'       Para cada ponto, o RFSI cria 2 × n.obs colunas: os valores observados
#'       nos n vizinhos e as distâncias até eles. Valores baixos (3-5) capturam
#'       estrutura local; valores altos (8-10) suavizam mais. É o hiperparâmetro
#'       mais importante do método.
#'
#'   num.trees (arvores)
#'       Número de árvores da floresta. Mais árvores = predição mais estável,
#'       com custo computacional linear. Retornos decrescentes após ~500.
#'
#'   mtry
#'       Número de variáveis candidatas sorteadas em cada split. Default aqui é 5.
#'
#'   splitrule = "variance"
#'       Critério de divisão para regressão (minimiza a variância intra-nó).
#'
#'   min.node.size
#'       Tamanho mínimo do nó terminal. Valores maiores regularizam mais
#'       (árvores mais rasas, menos overfitting). Default 10.
#'
#'   sample.fraction
#'       Fração de observações amostradas para construir cada árvore. Default 0.95
#'       (quase todas — apropriado quando o número de pontos de treino é pequeno).
#'
#'   importance = "impurity"
#'       Calcula a importância das variáveis pela redução de impureza. Útil para
#'       diagnosticar quais vizinhos/distâncias mais contribuem.
#'
#'   crs
#'       Sistema de referência de coordenadas. IMPORTANTE: o RFSI calcula
#'       distâncias euclidianas, então as coordenadas devem estar em um CRS
#'       PROJETADO (em metros), não em graus. EPSG:31982 = SIRGAS 2000 / UTM 22S,
#'       adequado para boa parte do Brasil Central. Ajuste conforme sua região.
#'
#' -----------------------------------------------------------------------------
#' DADOS DE ENTRADA
#' -----------------------------------------------------------------------------
#'
#'   train (w) — grade grosseira, os "pontos de treino"
#'       Deve conter: LAT, LONG, uma coluna identificadora (env) e as variáveis
#'       climáticas a serem interpoladas. Exemplo: envcov50km_M1.csv, com ~100
#'       pontos em grade 0.5° dentro do território da estrutura de melhoramento,
#'       contendo médias climáticas de 25 anos (2000-2024) do NASA POWER,
#'       processadas com EnvRtype::processWTH().
#'
#'   target (Z) — grade fina, os "pontos alvo"
#'       Deve conter: LAT, LONG e as variáveis de solo. Exemplo:
#'       soilgrid_M1.csv, com ~252.000 pontos em grade 0.01°, contendo 9
#'       variáveis do SoilGrids (bdod, clay, nitrogen, phh, sand, silt, soc,
#'       cfvo, ocd) extraídas na camada 5-15 cm.
#'
#'   SAÍDA: um data.frame com a grade fina completa — LAT, LONG, todas as
#'   variáveis de solo (originais) E todas as variáveis climáticas (interpoladas).
#'   Pronto para uso como matriz W na predição genômica de ambientes não avaliados.
#'
#' -----------------------------------------------------------------------------
#' DEPENDÊNCIAS
#' -----------------------------------------------------------------------------
#'   install.packages(c("meteo", "ranger", "sf", "terra", "dplyr", "ggplot2"))
#'
#' =============================================================================


#' Interpolação espacial de covariáveis ambientais via RFSI
#'
#' @param train data.frame. Grade de treino (grosseira) com LAT, LONG, coluna
#'   identificadora e as variáveis a interpolar.
#' @param target data.frame. Grade alvo (fina) com LAT, LONG e, opcionalmente,
#'   variáveis próprias (ex.: solo) que serão preservadas na saída.
#' @param vars character. Nomes das variáveis a interpolar. Se NULL (default),
#'   usa todas as colunas numéricas de `train` que não estão em `target` e não
#'   são coordenadas/identificadores.
#' @param id_col character. Nome da coluna identificadora em `train`, usada no
#'   leave-one-out. Default "env".
#' @param coords character(2). Nomes das colunas de longitude e latitude.
#'   Default c("LONG", "LAT").
#' @param crs integer. Código EPSG do CRS projetado (em metros). Default 31982.
#' @param neighbors integer. Valores de n.obs a testar no tuning. Default 3:10.
#' @param trees integer. Valores de num.trees a testar. Default seq(300, 1000, 200).
#' @param mtry integer. Variáveis candidatas por split. Default 5.
#' @param min_node_size integer. Tamanho mínimo do nó terminal. Default 10.
#' @param sample_fraction numeric. Fração amostrada por árvore. Default 0.95.
#' @param seed integer. Semente para reprodutibilidade. Default 42.
#' @param cpus integer. Núcleos para o tuning. Default detectCores() - 1.
#' @param plot logical. Se TRUE, gera mapas das variáveis interpoladas. Default FALSE.
#' @param plot_vars character. Quais variáveis plotar. Se NULL, plota todas.
#' @param plot_dir character. Pasta onde salvar os mapas. Se NULL, apenas exibe.
#' @param verbose logical. Imprimir progresso. Default TRUE.
#'
#' @return Uma lista com:
#'   \item{data}{data.frame da grade fina com solo + clima interpolado}
#'   \item{tuning}{data.frame com os hiperparâmetros vencedores por variável}
#'   \item{tuning_full}{lista com a grade completa de tuning de cada variável}
#'   \item{plots}{lista de objetos ggplot (se plot = TRUE)}
#'
#' @examples
#' \dontrun{
#' w <- read.csv("envcov50km_M1.csv")
#' Z <- read.csv("soilgrids_out/soilgrid_M1.csv")
#' res <- rfsi_envcov(train = w, target = Z, plot = TRUE)
#' write.csv(res$data, "W_matrix_untested_M1.csv", row.names = FALSE)
#' }
#'
#' @export
rfsi_envcov <- function(train,
                        target,
                        vars            = NULL,
                        id_col          = "env",
                        coords          = c("LONG", "LAT"),
                        crs             = 31982,
                        neighbors       = 3:10,
                        trees           = seq(300, 1000, by = 200),
                        mtry            = 5,
                        min_node_size   = 10,
                        sample_fraction = 0.95,
                        seed            = 42,
                        cpus            = parallel::detectCores() - 1,
                        plot            = FALSE,
                        plot_vars       = NULL,
                        plot_dir        = NULL,
                        verbose         = TRUE) {

  # ---- Checagem de dependências ----
  for (pkg in c("meteo", "ranger", "sf", "terra", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Pacote '", pkg, "' é necessário. Instale com install.packages('", pkg, "')")
    }
  }
  if (plot && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Pacote 'ggplot2' é necessário para plot = TRUE.")
  }

  lon_col <- coords[1]
  lat_col <- coords[2]

  # ---- Validação das entradas ----
  if (!all(coords %in% names(train)))  stop("`train` precisa conter as colunas: ", paste(coords, collapse = ", "))
  if (!all(coords %in% names(target))) stop("`target` precisa conter as colunas: ", paste(coords, collapse = ", "))
  if (!id_col %in% names(train))       stop("`train` precisa conter a coluna identificadora '", id_col, "'.")

  # ---- Definir variáveis a interpolar ----
  if (is.null(vars)) {
    candidatas <- setdiff(names(train), c(names(target), coords, id_col))
    numericas  <- candidatas[sapply(train[candidatas], is.numeric)]
    vars <- numericas
  }
  vars <- setdiff(vars, c(coords, id_col))

  if (length(vars) == 0) stop("Nenhuma variável para interpolar. Verifique `train` e `target`.")

  # Descartar variáveis constantes ou totalmente NA (quebram o ranger)
  ok <- sapply(vars, function(v) {
    x <- train[[v]]
    is.numeric(x) && sum(!is.na(x)) > 1 && stats::var(x, na.rm = TRUE) > 0
  })
  if (any(!ok)) {
    warning("Variáveis descartadas (constantes ou sem variação): ",
            paste(vars[!ok], collapse = ", "))
    vars <- vars[ok]
  }

  if (verbose) {
    cat("=============================================================\n")
    cat(" RFSI — Interpolação de covariáveis ambientais\n")
    cat(" Método: Sekulić et al. (2020), pacote meteo\n")
    cat("=============================================================\n")
    cat("Pontos de treino :", nrow(train), "\n")
    cat("Pontos alvo      :", nrow(target), "\n")
    cat("Variáveis        :", length(vars), "->", paste(vars, collapse = ", "), "\n")
    cat("Grade de tuning  :", length(neighbors), "vizinhos x", length(trees), "árvores =",
        length(neighbors) * length(trees), "combinações\n")
    cat("CRS              : EPSG:", crs, "\n\n", sep = "")
  }

  # ---- Preparar objetos espaciais (uma vez só, fora dos loops) ----
  train_sf <- sf::st_as_sf(train, coords = coords, crs = crs, remove = FALSE)
  train_sf[[lon_col]] <- sf::st_coordinates(train_sf)[, 1]
  train_sf[[lat_col]] <- sf::st_coordinates(train_sf)[, 2]

  target_sf <- sf::st_as_sf(target, coords = coords, crs = crs, remove = FALSE)
  target_sf[[lon_col]] <- sf::st_coordinates(target_sf)[, 1]
  target_sf[[lat_col]] <- sf::st_coordinates(target_sf)[, 2]

  ids <- unique(train[[id_col]])

  # Saída: começa com a grade alvo intacta (LAT, LONG e variáveis de solo)
  out         <- target
  tuning_best <- list()
  tuning_full <- list()

  # =========================================================================
  # LOOP PRINCIPAL — uma variável de cada vez
  # =========================================================================
  for (h in seq_along(vars)) {

    v  <- vars[h]
    fm <- stats::as.formula(paste(v, "~", lon_col, "+", lat_col))

    if (verbose) cat(sprintf("[%d/%d] %s\n", h, length(vars), v))

    # -----------------------------------------------------------------------
    # ETAPA 1 — Tuning por LOO-CV
    # -----------------------------------------------------------------------
    grade <- expand.grid(n_obs = neighbors, n_trees = trees)
    grade$cp <- NA_real_

    for (g in seq_len(nrow(grade))) {

      obs_vec  <- numeric(length(ids))
      pred_vec <- numeric(length(ids))

      for (k in seq_along(ids)) {

        treino_k <- train_sf[train_sf[[id_col]] != ids[k], ]
        alvo_k   <- train_sf[train_sf[[id_col]] == ids[k], ]

        mod_k <- try(
          meteo::rfsi(formula         = fm,
                      data            = treino_k,
                      n.obs           = grade$n_obs[g],
                      s.crs           = sf::st_crs(treino_k),
                      p.crs           = sf::st_crs(treino_k),
                      cpus            = cpus,
                      progress         = FALSE,
                      importance      = "impurity",
                      seed            = seed,
                      num.trees       = grade$n_trees[g],
                      mtry            = mtry,
                      splitrule       = "variance",
                      min.node.size   = min_node_size,
                      sample.fraction = sample_fraction),
          silent = TRUE)

        if (inherits(mod_k, "try-error")) { obs_vec[k] <- NA; pred_vec[k] <- NA; next }

        pr_k <- try(
          meteo::pred.rfsi(model         = mod_k,
                           data          = treino_k,
                           obs.col       = v,
                           newdata       = alvo_k,
                           output.format = "data.frame",
                           zero.tol      = 0,
                           s.crs         = sf::st_crs(treino_k),
                           newdata.s.crs = sf::st_crs(treino_k),
                           p.crs         = sf::st_crs(treino_k),
                           cpus          = 1,
                           progress      = FALSE),
          silent = TRUE)

        if (inherits(pr_k, "try-error")) { obs_vec[k] <- NA; pred_vec[k] <- NA; next }

        obs_vec[k]  <- alvo_k[[v]][1]
        pred_vec[k] <- as.data.frame(pr_k)$pred[1]
      }

      grade$cp[g] <- suppressWarnings(
        stats::cor(obs_vec, pred_vec, use = "complete.obs")
      )
    }

    grade <- grade[order(-grade$cp), ]
    melhor <- grade[1, ]

    if (verbose) {
      cat(sprintf("      tuning: n.obs = %d | num.trees = %d | cp = %.3f\n",
                  melhor$n_obs, melhor$n_trees, melhor$cp))
    }

    tuning_best[[v]] <- data.frame(variavel  = v,
                                   n_obs     = melhor$n_obs,
                                   n_trees   = melhor$n_trees,
                                   cp        = melhor$cp)
    tuning_full[[v]] <- grade

    # -----------------------------------------------------------------------
    # ETAPA 2 — Modelo final e predição na grade fina
    # -----------------------------------------------------------------------
    mod_final <- meteo::rfsi(formula         = fm,
                             data            = train_sf,
                             n.obs           = melhor$n_obs,
                             s.crs           = sf::st_crs(train_sf),
                             p.crs           = sf::st_crs(train_sf),
                             cpus            = cpus,
                             progress        = FALSE,
                             importance      = "impurity",
                             seed            = seed,
                             num.trees       = melhor$n_trees,
                             mtry            = mtry,
                             splitrule       = "variance",
                             min.node.size   = min_node_size,
                             sample.fraction = sample_fraction)

    pred_final <- meteo::pred.rfsi(model         = mod_final,   # <- modelo FINAL, não o do LOO
                                   data          = train_sf,
                                   obs.col       = v,
                                   newdata       = target_sf,
                                   output.format = "data.frame",
                                   zero.tol      = 0,
                                   s.crs         = sf::st_crs(train_sf),
                                   newdata.s.crs = sf::st_crs(train_sf),
                                   p.crs         = sf::st_crs(train_sf),
                                   cpus          = 1,
                                   progress      = verbose)

    # A ordem das linhas de pred.rfsi acompanha `newdata`, então a atribuição
    # é posicional — evita merge por coordenadas (que falha por erro de ponto flutuante)
    out[[v]] <- as.data.frame(pred_final)$pred

    if (verbose) {
      cat(sprintf("      predito: %d pontos | range [%.2f, %.2f]\n\n",
                  sum(!is.na(out[[v]])),
                  min(out[[v]], na.rm = TRUE),
                  max(out[[v]], na.rm = TRUE)))
    }
  }

  tuning_best <- do.call(rbind, tuning_best)
  rownames(tuning_best) <- NULL

  # =========================================================================
  # PLOTS (opcional)
  # =========================================================================
  plots <- list()

  if (plot) {
    if (is.null(plot_vars)) plot_vars <- vars
    plot_vars <- intersect(plot_vars, names(out))

    if (!is.null(plot_dir)) dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

    if (verbose) cat("Gerando mapas...\n")

    for (v in plot_vars) {

      p <- ggplot2::ggplot() +
        ggplot2::geom_point(
          data = out,
          ggplot2::aes(x = .data[[lon_col]], y = .data[[lat_col]], color = .data[[v]]),
          size = 0.1) +
        ggplot2::geom_point(
          data = train,
          ggplot2::aes(x = .data[[lon_col]], y = .data[[lat_col]]),
          color = "black", shape = 4, size = 1.6, stroke = 0.7) +
        ggplot2::scale_color_viridis_c(option = "C", na.value = "grey85") +
        ggplot2::coord_fixed() +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(
          title    = v,
          subtitle = "Superfície interpolada (RFSI) — x = pontos de treino",
          x = "Longitude", y = "Latitude", color = v)

      plots[[v]] <- p

      if (!is.null(plot_dir)) {
        ggplot2::ggsave(file.path(plot_dir, paste0(v, ".png")),
                        p, width = 8, height = 8, dpi = 150)
      } else {
        print(p)
      }
    }
    if (verbose && !is.null(plot_dir)) cat("Mapas salvos em: ", plot_dir, "\n", sep = "")
  }

  if (verbose) {
    cat("\n=============================================================\n")
    cat(" CONCLUÍDO\n")
    cat("=============================================================\n")
    cat("Matriz final:", nrow(out), "pontos x", ncol(out), "colunas\n")
    cat("  - solo  :", length(setdiff(names(target), coords)), "variáveis\n")
    cat("  - clima :", length(vars), "variáveis (interpoladas)\n")
    cat("cp média do tuning:", round(mean(tuning_best$cp, na.rm = TRUE), 3), "\n")
  }

  list(data        = out,
       tuning      = tuning_best,
       tuning_full = tuning_full,
       plots       = plots)
}


# =============================================================================
# EXEMPLO DE USO
# =============================================================================
#
# library(dplyr)
# source("rfsi_envcov.R")
#
# # ---- 1) Dados de entrada -----------------------------------------------
#
# # Grade grosseira (0.5°, ~50 km): clima do NASA POWER, médias 2000-2024,
# # processado com EnvRtype::processWTH(Tbase1 = 8, Tbase2 = 45,
# #                                     Topt1 = 30, Topt2 = 35)
# w <- read.csv("envcov50km_M1.csv")
#
# # Grade fina (0.01°, ~1 km): solo do SoilGrids, camada 5-15 cm
# Z <- read.csv("soilgrids_out/soilgrid_M1.csv")
#
# # ---- 2) Interpolação ----------------------------------------------------
#
# res <- rfsi_envcov(
#   train     = w,
#   target    = Z,
#   crs       = 31982,          # SIRGAS 2000 / UTM 22S
#   neighbors = 3:10,
#   trees     = seq(300, 1000, by = 200),
#   plot      = TRUE,
#   plot_vars = c("T2M", "PRECTOT", "RH2M", "VPD"),
#   plot_dir  = "mapas_rfsi/M1"
# )
#
# # ---- 3) Resultados ------------------------------------------------------
#
# W_untested <- res$data     # solo + clima na grade fina
# res$tuning                 # hiperparâmetros vencedores por variável
# res$tuning_full$T2M        # grade completa de tuning de uma variável
# res$plots$T2M              # mapa de uma variável
#
# write.csv(W_untested, "W_matrix_untested_M1.csv", row.names = FALSE)
#
# # ---- 4) Loop para todas as estruturas de melhoramento -------------------
#
# for (s in paste0("M", 1:5)) {
#   w_s <- read.csv(sprintf("envcov50km_%s.csv", s))
#   Z_s <- read.csv(sprintf("soilgrids_out/soilgrid_%s.csv", s))
#
#   res_s <- rfsi_envcov(train = w_s, target = Z_s,
#                        plot = TRUE, plot_dir = sprintf("mapas_rfsi/%s", s))
#
#   write.csv(res_s$data, sprintf("W_matrix_untested_%s.csv", s), row.names = FALSE)
#   write.csv(res_s$tuning, sprintf("tuning_rfsi_%s.csv", s), row.names = FALSE)
# }
#
# =============================================================================