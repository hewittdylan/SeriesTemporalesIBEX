# 03_asimetricos.R

# Setup

paquetes <- c(
  "quantmod", "forecast", "tseries", "lmtest", "tsoutliers",
  "FinTS", "rugarch", "xts", "zoo", "tidyverse"
)
nuevos <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]
if (length(nuevos) > 0) {
  install.packages(nuevos,
    dependencies = TRUE, repos = "https://cloud.r-project.org"
  )
}
invisible(lapply(paquetes, library, character.only = TRUE))

set.seed(42)

dir_figs <- "resultados/figuras"
dir_mods <- "resultados/modelos"
dir_tabs <- "resultados"
dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_mods, recursive = TRUE, showWarnings = FALSE)

# Configuración

activos <- list(
  list(ticker = "^IBEX",  nombre = "IBEX 35",         color_main = "darkblue",
       color_fill = "lightblue", suf = "IBEX"),
  list(ticker = "SAN.MC", nombre = "Banco Santander", color_main = "darkred",
       color_fill = "tomato",    suf = "SAN")
)

fecha_ini <- "2015-01-01"
fecha_fin <- "2025-12-31"

# preparar_datos() (idéntica a los scripts previos)

preparar_datos <- function(ticker, fecha_ini, fecha_fin) {
  pp <- quantmod::getSymbols(ticker,
    src = "yahoo", from = fecha_ini, to = fecha_fin, auto.assign = FALSE
  )
  precio_xts <- na.omit(quantmod::Ad(pp))
  colnames(precio_xts) <- "precio"
  ret_xts <- diff(log(precio_xts))[-1, ]
  ret_num <- as.numeric(ret_xts)
  ret_ts  <- ts(ret_num, frequency = 252)
  tso_ao  <- tsoutliers::tso(ret_ts,
    types = "AO", tsmethod = "arima",
    args.tsmethod = list(order = c(0, 0, 0), include.mean = TRUE),
    maxit.iloop = 4, cval = 4
  )
  ret_clean <- ret_num
  if (!is.null(tso_ao$outliers) && nrow(tso_ao$outliers) > 0) {
    idx <- tso_ao$outliers$ind
    for (i in idx) {
      if (i > 1 && i < length(ret_clean)) {
        ret_clean[i] <- (ret_clean[i - 1] + ret_clean[i + 1]) / 2
      } else if (i == 1) {
        ret_clean[i] <- ret_clean[i + 1]
      } else if (i == length(ret_clean)) {
        ret_clean[i] <- ret_clean[i - 1]
      }
    }
  }
  ret_clean_xts <- xts::xts(ret_clean, order.by = zoo::index(ret_xts))
  colnames(ret_clean_xts) <- "r"
  list(ticker = ticker, ret_clean_xts = ret_clean_xts, ret_clean = ret_clean)
}

# comparar_asimetricos()
# Ajusta los tres modelos para un activo y devuelve la tabla comparativa
# más la selección por reglas (AIC entre los significativos, o sGARCH si
# ninguno asimétrico aporta)

comparar_asimetricos <- function(activo_info, dir_figs, dir_mods) {
  ticker <- activo_info$ticker
  nombre <- activo_info$nombre
  suf    <- activo_info$suf
  color_main <- activo_info$color_main
  color_fill <- activo_info$color_fill
  message(sprintf("\n%s ", nombre))

  d <- preparar_datos(ticker, fecha_ini, fecha_fin)
  ret_clean_xts <- d$ret_clean_xts

  # Especificación de los tres modelos
  spec_sgarch <- ugarchspec(
    variance.model     = list(model = "sGARCH",   garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "std"
  )
  spec_gjr <- ugarchspec(
    variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "std"
  )
  spec_egarch <- ugarchspec(
    variance.model     = list(model = "eGARCH",   garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "std"
  )

  fit_sgarch <- ugarchfit(spec_sgarch, ret_clean_xts, solver = "hybrid")
  fit_gjr    <- ugarchfit(spec_gjr,    ret_clean_xts, solver = "hybrid")
  fit_egarch <- ugarchfit(spec_egarch, ret_clean_xts, solver = "hybrid")
  modelos <- list(sGARCH = fit_sgarch, GJR = fit_gjr, EGARCH = fit_egarch)

  # Extracción de métricas por modelo
  extraer_metricas <- function(fit, nm) {
    cf <- coef(fit)
    ic <- infocriteria(fit)
    sb <- signbias(fit)
    res_std <- as.numeric(residuals(fit, standardize = TRUE))
    lb_20    <- Box.test(res_std,   lag = 20, type = "Ljung-Box")
    lb_sq_20 <- Box.test(res_std^2, lag = 20, type = "Ljung-Box")
    # Acceso defensivo a las columnas de matcoef (varían entre versiones
    # de rugarch y a veces tienen espacios en los nombres)
    if ("gamma1" %in% names(cf)) {
      gamma_est <- as.numeric(cf["gamma1"])
      rob <- fit@fit$robust.matcoef
      mat <- if (!is.null(rob) && "gamma1" %in% rownames(rob)) rob else fit@fit$matcoef
      col_se <- grep("Std", colnames(mat), value = TRUE)[1]
      col_p  <- grep("^Pr|^p|p\\.value", colnames(mat),
                     ignore.case = TRUE, value = TRUE)[1]
      gamma_se <- as.numeric(mat["gamma1", col_se])
      gamma_p  <- as.numeric(mat["gamma1", col_p])
    } else {
      gamma_est <- NA_real_; gamma_se <- NA_real_; gamma_p <- NA_real_
    }
    persistencia <- tryCatch(persistence(fit), error = function(e) NA_real_)

    data.frame(
      Modelo = nm,
      AIC = round(ic[1, 1], 4), BIC = round(ic[2, 1], 4),
      LogLik = round(likelihood(fit), 2),
      mu     = signif(as.numeric(cf["mu"]), 3),
      omega  = signif(as.numeric(cf["omega"]), 3),
      alpha1 = round(as.numeric(cf["alpha1"]), 4),
      beta1  = round(as.numeric(cf["beta1"]), 4),
      gamma1 = round(gamma_est, 4),
      gamma1_p = signif(gamma_p, 3),
      shape_nu = round(as.numeric(cf["shape"]), 2),
      persistencia = round(persistencia, 4),
      LB_std_20_p   = round(lb_20$p.value, 4),
      LBsq_std_20_p = round(lb_sq_20$p.value, 4),
      SignBias_p      = round(sb["Sign Bias", "prob"], 4),
      SignBiasNeg_p   = round(sb["Negative Sign Bias", "prob"], 4),
      SignBiasPos_p   = round(sb["Positive Sign Bias", "prob"], 4),
      SignBiasJoint_p = round(sb["Joint Effect", "prob"], 4),
      stringsAsFactors = FALSE
    )
  }
  df_compar <- do.call(rbind, lapply(names(modelos), function(nm) {
    extraer_metricas(modelos[[nm]], nm)
  }))
  cat("\nComparativa:\n"); print(df_compar, row.names = FALSE)

  # Selección
  gjr_sig <- !is.na(df_compar$gamma1_p[df_compar$Modelo == "GJR"]) &&
             df_compar$gamma1_p[df_compar$Modelo == "GJR"] < 0.05
  egarch_sig <- !is.na(df_compar$gamma1_p[df_compar$Modelo == "EGARCH"]) &&
                df_compar$gamma1_p[df_compar$Modelo == "EGARCH"] < 0.05
  if (!gjr_sig && !egarch_sig) {
    seleccion <- "sGARCH"
    motivo <- "gamma_1 no significativo en ningún asimétrico"
  } else {
    candidatos_asim <- df_compar[df_compar$Modelo %in% c("GJR", "EGARCH"), ]
    candidatos_asim <- candidatos_asim[candidatos_asim$gamma1_p < 0.05, ]
    if (nrow(candidatos_asim) == 1) {
      seleccion <- candidatos_asim$Modelo
      motivo <- "único asimétrico con gamma_1 significativo"
    } else {
      seleccion <- candidatos_asim$Modelo[which.min(candidatos_asim$AIC)]
      motivo <- "menor AIC entre asimétricos significativos"
    }
  }
  cat(sprintf("Seleccionado: %s  (%s)\n", seleccion, motivo))

  # Volatilidad condicional comparada (de los tres modelos)
  vol_sgarch <- as.numeric(sigma(fit_sgarch)) * 100
  vol_gjr    <- as.numeric(sigma(fit_gjr))    * 100
  vol_egarch <- as.numeric(sigma(fit_egarch)) * 100
  df_vol <- data.frame(
    fecha   = rep(as.Date(zoo::index(sigma(fit_sgarch))), 3),
    vol_pct = c(vol_sgarch, vol_gjr, vol_egarch),
    modelo  = factor(rep(c("sGARCH", "GJR-GARCH", "EGARCH"),
                         each = length(vol_sgarch)),
                     levels = c("sGARCH", "GJR-GARCH", "EGARCH"))
  )
  p_vol <- ggplot(df_vol, aes(x = fecha, y = vol_pct, color = modelo)) +
    geom_line(alpha = 0.75, linewidth = 0.4) +
    scale_color_manual(values = c("sGARCH" = "gray45",
                                  "GJR-GARCH" = color_main,
                                  "EGARCH" = color_fill),
                       name = "Modelo") +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
    labs(title = sprintf("Volatilidad condicional comparada, %s", nombre),
         x = "Fecha", y = expression(sigma[t] ~ "(%)"),
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("03_Asimetricos_Volatilidad_%s.png", suf)),
    plot = p_vol, width = 11, height = 5, bg = "white"
  )

  # Curva de impacto de noticias (NIC)
  # Evalúa sigma_t como función del shock a_{t-1} sobre un rango simétrico
  # de innovaciones estandarizadas, manteniendo el resto del modelo
  # constante. Asimetría visible si la curva es más empinada por la
  # izquierda que por la derecha
  rango_a <- seq(-5, 5, by = 0.05)
  sigma_lag <- sd(d$ret_clean)
  cf_s <- coef(fit_sgarch); cf_g <- coef(fit_gjr); cf_e <- coef(fit_egarch)

  nic_sgarch <- as.numeric(cf_s["omega"]) +
    as.numeric(cf_s["alpha1"]) * (rango_a * sigma_lag)^2 +
    as.numeric(cf_s["beta1"])  * sigma_lag^2
  nic_gjr <- as.numeric(cf_g["omega"]) +
    (as.numeric(cf_g["alpha1"]) +
       as.numeric(cf_g["gamma1"]) * (rango_a < 0)) *
      (rango_a * sigma_lag)^2 +
    as.numeric(cf_g["beta1"]) * sigma_lag^2

  shape_e <- as.numeric(cf_e["shape"])
  E_abs_eps <- 2 * sqrt((shape_e - 2) / pi) *
    gamma((shape_e + 1) / 2) / ((shape_e - 1) * gamma(shape_e / 2))
  ln_sigma2_eg <- as.numeric(cf_e["omega"]) +
    as.numeric(cf_e["alpha1"]) * rango_a +
    as.numeric(cf_e["gamma1"]) * (abs(rango_a) - E_abs_eps) +
    as.numeric(cf_e["beta1"])  * log(sigma_lag^2)
  nic_egarch <- exp(ln_sigma2_eg)

  df_nic <- data.frame(
    a      = rep(rango_a * sigma_lag, 3),
    sigma2 = c(nic_sgarch, nic_gjr, nic_egarch),
    modelo = factor(rep(c("sGARCH", "GJR-GARCH", "EGARCH"),
                        each = length(rango_a)),
                    levels = c("sGARCH", "GJR-GARCH", "EGARCH"))
  )
  p_nic <- ggplot(df_nic, aes(x = a, y = sqrt(sigma2) * 100, color = modelo)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
    scale_color_manual(values = c("sGARCH" = "gray45",
                                  "GJR-GARCH" = color_main,
                                  "EGARCH" = color_fill),
                       name = "Modelo") +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
    labs(title = sprintf("Curva de impacto de noticias (NIC), %s", nombre),
         x = expression(a[t - 1]), y = expression(sigma[t] ~ "(% diario)"),
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("03_Asimetricos_NIC_%s.png", suf)),
    plot = p_nic, width = 9, height = 6, bg = "white"
  )

  saveRDS(fit_sgarch, file.path(dir_mods, sprintf("asimetrico_sGARCH_%s.rds", suf)))
  saveRDS(fit_gjr,    file.path(dir_mods, sprintf("asimetrico_GJR_%s.rds", suf)))
  saveRDS(fit_egarch, file.path(dir_mods, sprintf("asimetrico_EGARCH_%s.rds", suf)))
  write.csv(df_compar,
    file.path(dir_tabs, sprintf("tabla_asimetricos_%s.csv", suf)),
    row.names = FALSE)

  list(activo = nombre, suf = suf, df_compar = df_compar,
       seleccion = seleccion, motivo = motivo, fits = modelos)
}

# Ejecución sobre IBEX y Santander

resultados_asim <- lapply(activos, function(a) {
  tryCatch(comparar_asimetricos(a, dir_figs, dir_mods),
           error = function(e) {
             message(sprintf("ERROR en %s: %s", a$ticker, conditionMessage(e)))
             NULL
           })
})
names(resultados_asim) <- sapply(activos, function(a) a$suf)

# Tabla resumen unificada

filas_resumen <- do.call(rbind, lapply(resultados_asim, function(r) {
  if (is.null(r)) return(NULL)
  df <- r$df_compar
  df$Activo <- r$activo
  df$Optimo <- df$Modelo == r$seleccion
  df[, c("Activo", "Modelo", "AIC", "BIC", "alpha1", "beta1", "gamma1",
         "gamma1_p", "persistencia", "shape_nu", "LB_std_20_p",
         "LBsq_std_20_p", "SignBiasNeg_p", "Optimo")]
}))
cat("\nTabla resumen unificada\n")
print(filas_resumen, row.names = FALSE)
write.csv(filas_resumen,
  file.path(dir_tabs, "tabla_asimetricos_resumen.csv"), row.names = FALSE)

cat("\nConclusión por activo\n")
for (r in resultados_asim) {
  if (is.null(r)) next
  cat(sprintf("  %-18s -> %-8s  (%s)\n", r$activo, r$seleccion, r$motivo))
}

message("\nBloque asimétrico terminado")
