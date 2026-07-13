#' =============================================================================
#' get_envcov() — Covariáveis ambientais (clima + solo) para qualquer conjunto
#'                de pontos geográficos
#' =============================================================================
#'
#' Autor:  João Marcos Amario
#' E-mail: joao.m.sousa@ufv.br
#'
#' -----------------------------------------------------------------------------
#' O QUE ESTA FUNÇÃO FAZ
#' -----------------------------------------------------------------------------
#' Uma função única que monta a matriz de covariáveis ambientais (W) para
#' qualquer conjunto de pontos, cobrindo os dois cenários do fluxo de trabalho:
#'
#'   CENÁRIO A — Ambientes AVALIADOS (experimentos reais)
#'     Você tem os ensaios com latitude, longitude e as datas reais de plantio
#'     e colheita de cada um. Quer clima E solo, na mesma matriz.
#'     -> Passe `points` com as colunas de data. Use what = "both".
#'     -> O clima é baixado no ciclo REAL de cada experimento.
#'
#'   CENÁRIO B — Ambientes NÃO AVALIADOS (grade de predição)
#'     Você tem uma grade de pontos (0.5°, 0.01°, o que for) e quer caracterizar
#'     o ambiente. Não existem datas reais, então define-se uma janela típica
#'     (data de plantio + duração do ciclo) e uma série de anos para promediar.
#'     -> Passe `points` sem datas, mais `planting` (MM-DD), `cycle` e `years`.
#'     -> Pode pedir só clima (what = "weather") ou só solo (what = "soil"),
#'        já que cada fonte costuma usar uma grade de resolução diferente:
#'        o NASA POWER não passa de ~0.5°, enquanto o SoilGrids chega a 250 m.
#'
#' Em ambos os casos, quando what = "both", clima e solo saem juntos na mesma
#' matriz, pareados por LAT/LONG.
#'
#' -----------------------------------------------------------------------------
#' FONTES DE DADOS
#' -----------------------------------------------------------------------------
#' CLIMA  NASA POWER, via EnvRtype::get_weather(). Baixa a série diária no
#'        intervalo pedido e, com EnvRtype::processWTH(), deriva as variáveis
#'        ecofisiológicas (GDD, FRUE, PETP, VPD, etc.) a partir das temperaturas
#'        cardinais informadas. Em seguida promedia a série de cada ponto.
#'
#' SOLO   SoilGrids (ISRIC), lido de rasters .tif já baixados em disco. A função
#'        NÃO baixa os rasters — extrai os valores nos pontos com terra::extract().
#'        Para obter os rasters, use geodata::soil_world() uma única vez.
#'
#' -----------------------------------------------------------------------------
#' TEMPERATURAS CARDINAIS (parâmetros do processWTH)
#' -----------------------------------------------------------------------------
#'   Tbase1  Temperatura basal inferior. Abaixo dela não há desenvolvimento.
#'   Tbase2  Temperatura basal superior. Acima dela o desenvolvimento cessa.
#'   Topt1   Início do intervalo ótimo.
#'   Topt2   Fim do intervalo ótimo.
#'
#' Esses quatro valores definem a curva de resposta térmica usada para calcular
#' o GDD (graus-dia acumulados) e o FRUE (fator de eficiência de uso da radiação
#' limitado pela temperatura). Os defaults (8/45/30/35) são típicos para milho.
#' Ajuste para sua cultura — usar os mesmos valores nos ambientes avaliados e nos
#' não avaliados é obrigatório, senão as matrizes não são comparáveis.
#'
#' -----------------------------------------------------------------------------
#' DEPENDÊNCIAS
#' -----------------------------------------------------------------------------
#'   install.packages(c("EnvRtype", "terra", "dplyr"))
#'
#' =============================================================================
#'
#' -----------------------------------------------------------------------------
#' EXEMPLOS DE USO
#' -----------------------------------------------------------------------------
#'
#' source("get_envcov.R")
#'
#' ## --- CENÁRIO A: experimentos avaliados (clima no ciclo real + solo) -------
#'
#' exp <- read.table("data_M1.txt", header = TRUE, sep = "\t")
#'
#' pts <- exp %>%
#'   distinct(env, .keep_all = TRUE) %>%
#'   transmute(env,
#'             LAT   = as.numeric(lat),
#'             LONG  = -abs(as.numeric(long)),   # garante hemisfério oeste
#'             start = as.Date(planting_date),
#'             end   = as.Date(harvest_date))
#'
#' W_M1 <- get_envcov(
#'   points    = pts,
#'   what      = "both",
#'   soil_path = "Rasters/soil_world",
#'   Tbase1 = 8, Tbase2 = 45, Topt1 = 30, Topt2 = 35
#' )
#'
#' write.csv(W_M1, "envcov_M1.csv", row.names = FALSE)
#'
#'
#' ## --- CENÁRIO B1: grade grossa (0.5°) — só CLIMA ---------------------------
#' ##     O NASA POWER não tem resolução melhor que isso.
#'
#' grid_05 <- read.csv("grid_0.5_M1.csv")     # precisa de LAT, LONG
#'
#' clima_50km <- get_envcov(
#'   points   = grid_05,
#'   what     = "weather",
#'   planting = "10-15",        # data típica de plantio (MM-DD)
#'   cycle    = 150,            # duração do ciclo, em dias
#'   years    = 2000:2024,      # anos a promediar
#'   cache_dir = "NASA_Power/M1"
#' )
#'
#' write.csv(clima_50km, "envcov50km_M1.csv", row.names = FALSE)
#'
#'
#' ## --- CENÁRIO B2: grade fina (0.01°) — só SOLO -----------------------------
#' ##     O SoilGrids tem 250 m, então aproveita a grade fina.
#'
#' grid_001 <- read.csv("grid_0.01_M1.csv")
#'
#' solo_1km <- get_envcov(
#'   points    = grid_001,
#'   what      = "soil",
#'   soil_path = "Rasters/soil_world"
#' )
#'
#' write.csv(solo_1km, "soilgrid_M1.csv", row.names = FALSE)
#'
#'
#' ## --- Loop nas cinco estruturas de melhoramento ---------------------------
#'
#' janelas <- list(M1 = list(planting = "10-15", cycle = 150),
#'                 M2 = list(planting = "11-01", cycle = 145),
#'                 M3 = list(planting = "10-20", cycle = 160),
#'                 M4 = list(planting = "11-10", cycle = 140),
#'                 M5 = list(planting = "12-01", cycle = 135))
#'
#' for (s in paste0("M", 1:5)) {
#'   g05  <- read.csv(sprintf("grid_0.5_%s.csv", s))
#'   g001 <- read.csv(sprintf("grid_0.01_%s.csv", s))
#'
#'   clima <- get_envcov(g05, what = "weather",
#'                       planting = janelas[[s]]$planting,
#'                       cycle    = janelas[[s]]$cycle,
#'                       years    = 2000:2024,
#'                       cache_dir = file.path("NASA_Power", s))
#'
#'   solo  <- get_envcov(g001, what = "soil", soil_path = "Rasters/soil_world")
#'
#'   write.csv(clima, sprintf("envcov50km_%s.csv", s), row.names = FALSE)
#'   write.csv(solo,  sprintf("soilgrid_%s.csv", s),   row.names = FALSE)
#' }
#'
#' =============================================================================


#' Covariáveis ambientais (clima e/ou solo) para um conjunto de pontos
#'
#' @param points data.frame. Precisa conter LAT e LONG. Para clima com datas
#'   por ponto (experimentos reais), inclua também `start` e `end` (Date ou
#'   character "YYYY-MM-DD"). Colunas extras (env, id, UF...) são preservadas
#'   na saída.
#' @param what character. "both" (clima + solo), "weather" (só clima) ou
#'   "soil" (só solo). Default "both".
#' @param coords character(2). Nomes das colunas de latitude e longitude.
#'   Default c("LAT", "LONG").
#'
#' @param planting character "MM-DD". Data típica de plantio, usada quando
#'   `points` não traz datas próprias. Ex.: "10-15".
#' @param cycle integer. Duração do ciclo em dias, usada junto com `planting`.
#' @param years integer. Anos a baixar e promediar. Default 2000:2024. Só se
#'   aplica quando a janela vem de `planting`/`cycle`.
#'
#' @param Tbase1,Tbase2,Topt1,Topt2 numeric. Temperaturas cardinais da cultura,
#'   passadas a EnvRtype::processWTH(). Defaults 8/45/30/35 (milho).
#'
#' @param soil_path character. Pasta com os rasters .tif do SoilGrids.
#' @param soil_vars character. Variáveis de solo a extrair. Default: as nove
#'   padrão (bdod, clay, nitrogen, phh2o, sand, silt, soc, cfvo, ocd).
#' @param soil_depth character. Camada de profundidade, usada para montar o nome
#'   do arquivo. Default "5-15cm".
#' @param soil_res character. Sufixo de resolução no nome do arquivo.
#'   Default "30s". O padrão de nome esperado é:
#'   \code{<var>_<depth>_mean_<res>.tif}, ex.: "clay_5-15cm_mean_30s.tif".
#'
#' @param batch_size integer. Pontos por requisição ao NASA POWER. Default 200.
#' @param sleep numeric. Pausa entre requisições, em segundos. Default 10.
#' @param retries integer. Tentativas por lote em caso de falha. Default 3.
#' @param cache_dir character ou NULL. Se informado, cada lote baixado é salvo
#'   como .rds nessa pasta e reaproveitado se a função for chamada de novo.
#'   Fortemente recomendado para grades grandes.
#' @param drop_vars character. Colunas do clima a descartar ao final.
#'   Default c("DOY", "daysFromStart", "FROST_DAYS").
#' @param verbose logical. Imprimir progresso. Default TRUE.
#'
#' @return data.frame com as colunas originais de `points` mais as covariáveis
#'   pedidas. Traz o atributo "failed" com os lotes que não puderam ser baixados.
#'
#' @export
get_envcov <- function(points,
                       what       = c("both", "weather", "soil"),
                       coords     = c("LAT", "LONG"),
                       
                       # janela do clima (cenário B)
                       planting   = NULL,
                       cycle      = NULL,
                       years      = 2000:2024,
                       
                       # temperaturas cardinais
                       Tbase1 = 8, Tbase2 = 45, Topt1 = 30, Topt2 = 35,
                       
                       # solo
                       soil_path  = "Rasters/soil_world",
                       soil_vars  = c("bdod", "clay", "nitrogen", "phh2o",
                                      "sand", "silt", "soc", "cfvo", "ocd"),
                       soil_depth = "5-15cm",
                       soil_res   = "30s",
                       
                       # download
                       batch_size = 200,
                       sleep      = 10,
                       retries    = 3,
                       cache_dir  = NULL,
                       drop_vars  = c("DOY", "daysFromStart", "FROST_DAYS"),
                       verbose    = TRUE) {
  
  what    <- match.arg(what)
  lat_col <- coords[1]
  lon_col <- coords[2]
  
  say <- function(...) if (verbose) cat(...)
  
  # ==========================================================================
  # VALIDAÇÃO DAS ENTRADAS
  # ==========================================================================
  if (!is.data.frame(points))   stop("`points` precisa ser um data.frame.")
  if (nrow(points) == 0)        stop("`points` está vazio.")
  if (!all(coords %in% names(points))) {
    stop("`points` precisa das colunas ", paste(coords, collapse = " e "),
         ". Encontrado: ", paste(names(points), collapse = ", "))
  }
  
  points[[lat_col]] <- as.numeric(points[[lat_col]])
  points[[lon_col]] <- as.numeric(points[[lon_col]])
  
  if (anyNA(points[[lat_col]]) || anyNA(points[[lon_col]])) {
    stop("Há coordenadas NA em `points`. Limpe antes de prosseguir.")
  }
  if (any(abs(points[[lat_col]]) > 90) || any(abs(points[[lon_col]]) > 180)) {
    stop("Coordenadas fora do intervalo válido. Confira se LAT e LONG não estão trocadas.")
  }
  
  precisa_clima <- what %in% c("both", "weather")
  precisa_solo  <- what %in% c("both", "soil")
  
  if (precisa_clima && !requireNamespace("EnvRtype", quietly = TRUE)) {
    stop("Pacote 'EnvRtype' é necessário para baixar clima.")
  }
  if (precisa_solo && !requireNamespace("terra", quietly = TRUE)) {
    stop("Pacote 'terra' é necessário para extrair solo.")
  }
  
  out    <- points
  falhas <- NULL
  
  # ==========================================================================
  # BLOCO 1 — CLIMA (NASA POWER)
  # ==========================================================================
  if (precisa_clima) {
    
    tem_datas <- all(c("start", "end") %in% names(points))
    
    # ---- Decidir de onde vem a janela temporal --------------------------
    if (tem_datas) {
      modo <- "datas por ponto"
      if (!is.null(planting)) {
        warning("`points` já traz start/end — `planting`/`cycle` serão ignorados.")
      }
      anos_loop <- NA   # um único passe, sem loop de anos
    } else {
      modo <- "janela típica"
      if (is.null(planting) || is.null(cycle)) {
        stop("Sem as colunas start/end em `points`, é preciso informar ",
             "`planting` (\"MM-DD\") e `cycle` (dias).")
      }
      if (!grepl("^\\d{2}-\\d{2}$", planting)) {
        stop("`planting` deve estar no formato \"MM-DD\", ex.: \"10-15\".")
      }
      if (!is.numeric(cycle) || cycle <= 0) stop("`cycle` deve ser um número positivo de dias.")
      anos_loop <- years
    }
    
    say("=============================================================\n")
    say(" CLIMA — NASA POWER\n")
    say("=============================================================\n")
    say("Pontos    : ", nrow(points), "\n")
    say("Janela    : ", modo, "\n")
    if (modo == "janela típica") {
      say("            plantio ", planting, " | ciclo ", cycle, " dias | anos ",
          min(years), "-", max(years), "\n")
    }
    say("Cardinais : Tbase1=", Tbase1, " Tbase2=", Tbase2,
        " Topt1=", Topt1, " Topt2=", Topt2, "\n\n")
    
    if (!is.null(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    
    # ---- Baixa um conjunto de pontos em lotes, com retry e cache ---------
    baixar_lotes <- function(df, tag) {
      
      idx <- split(seq_len(nrow(df)), ceiling(seq_len(nrow(df)) / batch_size))
      res <- vector("list", length(idx))
      
      for (i in seq_along(idx)) {
        
        cache_file <- if (!is.null(cache_dir)) {
          file.path(cache_dir, sprintf("%s_lote%03d.rds", tag, i))
        } else NULL
        
        # Reaproveita o que já foi baixado
        if (!is.null(cache_file) && file.exists(cache_file)) {
          res[[i]] <- readRDS(cache_file)
          say(sprintf("  [%s] lote %d/%d — cache\n", tag, i, length(idx)))
          next
        }
        
        lote <- df[idx[[i]], ]
        ok   <- FALSE
        
        for (tent in seq_len(retries)) {
          
          say(sprintf("  [%s] lote %d/%d — tentativa %d\n", tag, i, length(idx), tent))
          
          w <- try(
            EnvRtype::get_weather(env.id    = lote$.env_dl,
                                  lat       = lote[[lat_col]],
                                  lon       = lote[[lon_col]],
                                  start.day = lote$.start,
                                  end.day   = lote$.end),
            silent = TRUE)
          
          if (!inherits(w, "try-error") && !is.null(w) && nrow(w) > 0) {
            res[[i]] <- w
            if (!is.null(cache_file)) saveRDS(w, cache_file)
            ok <- TRUE
            break
          }
          
          # Backoff progressivo: a API costuma bloquear por excesso de requisições
          if (tent < retries) Sys.sleep(60 * tent)
        }
        
        if (!ok) {
          warning("Falha no lote ", i, " de ", tag, " após ", retries, " tentativas.")
          falhas <<- c(falhas, sprintf("%s_lote%d", tag, i))
        }
        
        Sys.sleep(sleep)
      }
      
      res <- res[!vapply(res, is.null, logical(1))]
      if (length(res) == 0) return(NULL)
      do.call(rbind, res)
    }
    
    # ---- Executar o download ---------------------------------------------
    if (modo == "datas por ponto") {
      
      df <- points
      df$.start <- as.Date(df$start)
      df$.end   <- as.Date(df$end)
      
      if (anyNA(df$.start) || anyNA(df$.end)) stop("Há datas inválidas em start/end.")
      if (any(df$.end <= df$.start))          stop("Há linhas com end <= start.")
      
      # Identificador interno para a API (independente das colunas do usuário)
      df$.env_dl <- sprintf("P%05d", seq_len(nrow(df)))
      
      bruto <- baixar_lotes(df, "exp")
      if (is.null(bruto)) stop("Nenhum dado climático foi baixado.")
      
      proc <- EnvRtype::processWTH(env.data = bruto,
                                   Tbase1 = Tbase1, Tbase2 = Tbase2,
                                   Topt1  = Topt1,  Topt2  = Topt2)
      
      clima <- agregar_clima(proc, id_col = "env", drop_vars = drop_vars)
      names(clima)[names(clima) == "env"] <- ".env_dl"
      
      out <- merge(cbind(out, .env_dl = df$.env_dl), clima, by = ".env_dl", all.x = TRUE)
      out$.env_dl <- NULL
      
    } else {
      
      # Janela típica: baixa ano a ano e promedia a série toda
      por_ano <- vector("list", length(anos_loop))
      
      for (k in seq_along(anos_loop)) {
        
        a  <- anos_loop[k]
        say("\n--- Ano ", a, " ---\n")
        
        df <- points
        df$.start  <- as.Date(paste0(a, "-", planting))
        df$.end    <- df$.start + cycle
        df$.env_dl <- sprintf("P%05d", seq_len(nrow(df)))
        
        if (anyNA(df$.start)) {
          warning("Data de plantio inválida para o ano ", a, " — ano pulado.")
          next
        }
        
        bruto <- baixar_lotes(df, paste0("y", a))
        if (is.null(bruto)) next
        
        proc <- EnvRtype::processWTH(env.data = bruto,
                                     Tbase1 = Tbase1, Tbase2 = Tbase2,
                                     Topt1  = Topt1,  Topt2  = Topt2)
        
        por_ano[[k]] <- agregar_clima(proc, id_col = "env", drop_vars = drop_vars)
      }
      
      por_ano <- por_ano[!vapply(por_ano, is.null, logical(1))]
      if (length(por_ano) == 0) stop("Nenhum ano foi baixado com sucesso.")
      
      say("\nPromediando ", length(por_ano), " anos...\n")
      
      todos <- do.call(rbind, por_ano)
      num   <- setdiff(names(todos)[vapply(todos, is.numeric, logical(1))], "env")
      
      clima <- stats::aggregate(todos[num], by = list(env = todos$env),
                                FUN = mean, na.rm = TRUE)
      names(clima)[names(clima) == "env"] <- ".env_dl"
      
      out <- merge(cbind(out, .env_dl = sprintf("P%05d", seq_len(nrow(out)))),
                   clima, by = ".env_dl", all.x = TRUE)
      out$.env_dl <- NULL
    }
    
    n_clima <- length(setdiff(names(out), names(points)))
    say("\nClima: ", n_clima, " variáveis adicionadas\n\n")
  }
  
  # ==========================================================================
  # BLOCO 2 — SOLO (SoilGrids)
  # ==========================================================================
  if (precisa_solo) {
    
    say("=============================================================\n")
    say(" SOLO — SoilGrids\n")
    say("=============================================================\n")
    
    if (!dir.exists(soil_path)) {
      stop("Pasta de rasters não encontrada: ", soil_path,
           "\nBaixe antes com geodata::soil_world().")
    }
    
    arquivos <- file.path(soil_path,
                          sprintf("%s_%s_mean_%s.tif", soil_vars, soil_depth, soil_res))
    existe   <- file.exists(arquivos)
    
    if (!any(existe)) {
      stop("Nenhum raster encontrado em ", soil_path,
           "\nEsperado o padrão: <var>_", soil_depth, "_mean_", soil_res, ".tif",
           "\nArquivos presentes: ",
           paste(utils::head(list.files(soil_path, pattern = "\\.tif$"), 5), collapse = ", "))
    }
    if (any(!existe)) {
      warning("Rasters ausentes (ignorados): ",
              paste(soil_vars[!existe], collapse = ", "))
    }
    
    arquivos  <- arquivos[existe]
    vars_ok   <- soil_vars[existe]
    
    say("Rasters   : ", length(arquivos), " de ", length(soil_vars), "\n")
    say("Camada    : ", soil_depth, "\n")
    say("Extraindo : ", nrow(points), " pontos...\n")
    
    stk <- terra::rast(arquivos)
    # phh2o -> phh, mantendo o nome curto usado no restante do fluxo
    names(stk) <- sub("^phh2o$", "phh", vars_ok)
    
    pts  <- terra::vect(points[, c(lon_col, lat_col)],
                        geom = c(lon_col, lat_col), crs = "EPSG:4326")
    vals <- terra::extract(stk, pts)
    vals$ID <- NULL
    
    out <- cbind(out, vals)
    
    n_na <- sum(!stats::complete.cases(vals))
    if (n_na > 0) {
      say("Aviso     : ", n_na, " pontos sem dado de solo (água, área urbana ",
          "ou fora da cobertura)\n")
    }
    say("Solo      : ", ncol(vals), " variáveis adicionadas\n\n")
  }
  
  # ==========================================================================
  # RESUMO
  # ==========================================================================
  say("=============================================================\n")
  say(" CONCLUÍDO\n")
  say("=============================================================\n")
  say("Matriz final: ", nrow(out), " pontos x ", ncol(out), " colunas\n")
  say("Novas covariáveis: ", ncol(out) - ncol(points), "\n")
  
  if (!is.null(falhas)) {
    say("\nLotes que falharam (", length(falhas), "): ",
        paste(utils::head(falhas, 10), collapse = ", "), "\n")
    say("Rode a função de novo com o mesmo `cache_dir` para tentar só o que faltou.\n")
  }
  
  attr(out, "failed") <- falhas
  out
}


#' Agrega a série diária do processWTH em uma linha por ponto (auxiliar interna)
#' @keywords internal
agregar_clima <- function(proc, id_col = "env", drop_vars = NULL) {
  
  proc <- as.data.frame(proc)
  
  # LON/LAT/YYYYMMDD voltam do NASA POWER e não devem ser promediadas
  descartar <- c("LON", "LAT", "YYYYMMDD", drop_vars)
  proc      <- proc[, !names(proc) %in% descartar, drop = FALSE]
  
  # processWTH às vezes duplica colunas com sufixo ".1"
  proc <- proc[, !grepl("\\.1$", names(proc)), drop = FALSE]
  
  num <- setdiff(names(proc)[vapply(proc, is.numeric, logical(1))], id_col)
  if (length(num) == 0) stop("Nenhuma variável numérica após o processWTH.")
  
  ag <- stats::aggregate(proc[num], by = setNames(list(proc[[id_col]]), id_col),
                         FUN = mean, na.rm = TRUE)
  ag
}