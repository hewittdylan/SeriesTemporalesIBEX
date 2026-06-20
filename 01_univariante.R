# 01_univariante.R

# Setup

paquetes <- c(
  "quantmod", "forecast", "tseries", "lmtest", "tsoutliers", "nortest",
  "FinTS", "rugarch", "xts", "zoo", "tidyverse", "patchwork"
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

# Configuración de los activos

activos <- list(
  list(ticker = "^IBEX",   nombre = "IBEX 35",         color_main = "darkblue",   color_fill = "lightblue",  suf = "IBEX"),
  list(ticker = "SAN.MC",  nombre = "Banco Santander", color_main = "darkred",    color_fill = "tomato",     suf = "SAN"),
  list(ticker = "BBVA.MC", nombre = "BBVA",            color_main = "navy",       color_fill = "steelblue",  suf = "BBVA"),
  list(ticker = "CABK.MC", nombre = "CaixaBank",       color_main = "darkgreen",  color_fill = "palegreen",  suf = "CABK"),
  list(ticker = "SAB.MC",  nombre = "Banco Sabadell",  color_main = "darkorange", color_fill = "navajowhite",suf = "SAB"),
  list(ticker = "BKT.MC",  nombre = "Bankinter",       color_main = "purple",     color_fill = "plum",       suf = "BKT")
)

fecha_ini <- "2015-01-01"
fecha_fin <- "2025-12-31"
h_fore    <- 12  # horizonte de predicción en días bursátiles


# Pipeline por activo
# Encapsulamos todo el análisis univariante en una función para aplicarla
# de forma uniforme a los seis activos

analizar_activo <- function(ticker, nombre, color_main, color_fill, suf,
                            fecha_ini, fecha_fin, h_fore,
                            dir_figs, dir_mods) {
  message(sprintf("\n- %s (%s)", nombre, ticker))

  # Descarga del precio
  precio_raw <- quantmod::getSymbols(ticker,
    src = "yahoo", from = fecha_ini, to = fecha_fin, auto.assign = FALSE
  )
  precio_xts <- na.omit(quantmod::Ad(precio_raw))
  colnames(precio_xts) <- "precio"
  cat(sprintf("n = %d  |  %s -> %s\n",
    nrow(precio_xts), format(start(precio_xts)), format(end(precio_xts))
  ))

  df_p <- data.frame(fecha = index(precio_xts), precio = as.numeric(precio_xts))
  p_precio <- ggplot(df_p, aes(x = fecha, y = precio)) +
    geom_line(color = color_main, alpha = 0.85, linewidth = 0.6) +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("Precio de cierre ajustado, %s", nombre),
         x = "Fecha", y = "Precio",
         caption = "Fuente: Yahoo Finance | Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("01_Precio_%s.png", suf)),
    plot = p_precio, width = 10, height = 5, bg = "white"
  )

  # Log-rendimientos: r_t = log(P_t / P_{t-1})
  ret_xts <- diff(log(precio_xts))[-1, ]
  colnames(ret_xts) <- "r"
  ret_num <- as.numeric(ret_xts)

  # Panel precio vs rendimientos (no estacionario vs estacionario)
  df_r <- data.frame(fecha = index(ret_xts), r = ret_num)
  p_ret_only <- ggplot(df_r, aes(x = fecha, y = r)) +
    geom_line(color = color_main, alpha = 0.7, linewidth = 0.4) +
    theme_minimal(base_size = 13) +
    labs(title = "Log-rendimientos diarios", x = "Fecha", y = expression(r[t]))
  p_panel <- (p_precio + labs(title = sprintf("Precio, %s", nombre),
                              subtitle = NULL, caption = NULL)) /
             (p_ret_only + labs(title = "Log-rendimientos"))
  ggsave(file.path(dir_figs, sprintf("01_Panel_Precio_Retornos_%s.png", suf)),
    plot = p_panel, width = 10, height = 8, bg = "white"
  )

  # Descriptivos básicos
  desc <- list(
    media    = mean(ret_num),
    sd       = sd(ret_num),
    skewness = mean((ret_num - mean(ret_num))^3) / sd(ret_num)^3,
    kurtosis = mean((ret_num - mean(ret_num))^4) / sd(ret_num)^4
  )
  cat(sprintf("media=%.3e  sd=%.3e  skew=%.3f  kurt=%.3f\n",
    desc$media, desc$sd, desc$skewness, desc$kurtosis))

  # Periodograma: si los rendimientos son aproximadamente ruido blanco,
  # el espectro debería ser plano (sin periodicidad latente)
  spec_ret <- spectrum(ret_num, plot = FALSE)
  p_spec <- ggplot(data.frame(freq = spec_ret$freq, spec = spec_ret$spec),
                   aes(x = freq, y = spec)) +
    geom_line(color = color_main, alpha = 0.85, linewidth = 0.5) +
    scale_y_log10() +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("Periodograma, %s", nombre),
         x = "Frecuencia", y = "Espectro",
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("01_Periodograma_%s.png", suf)),
    plot = p_spec, width = 8, height = 5, bg = "white"
  )

  # Atípicos aditivos. Fijamos tsmethod = "arima" con orden (0,0,0) por
  # eficiencia (auto.arima dentro de tso() es muy lento para 2500 obs.)
  # cval = 4 etiqueta solo rendimientos realmente extremos
  ret_ts <- ts(ret_num, frequency = 252)
  tso_ao <- tsoutliers::tso(ret_ts,
    types = "AO",
    tsmethod = "arima",
    args.tsmethod = list(order = c(0, 0, 0), include.mean = TRUE),
    maxit.iloop = 4, cval = 4
  )
  ret_clean <- ret_num
  n_ao <- if (is.null(tso_ao$outliers)) 0 else nrow(tso_ao$outliers)
  if (n_ao > 0) {
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
  cat(sprintf("AOs detectados: %d\n", n_ao))
  ret_clean_xts <- xts::xts(ret_clean, order.by = zoo::index(ret_xts))
  colnames(ret_clean_xts) <- "r"

  # Identificación ARMA
  p_acf <- ggAcf(ret_clean, lag.max = 40) + theme_minimal(base_size = 13) +
    labs(title = sprintf("FAS, %s", nombre), x = "Retardo", y = "ACF",
         caption = "Fuente: Elaboración propia.")
  p_pacf <- ggPacf(ret_clean, lag.max = 40) + theme_minimal(base_size = 13) +
    labs(title = sprintf("FAP, %s", nombre), x = "Retardo", y = "PACF",
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("01_FAS_%s.png", suf)),
    plot = p_acf, width = 8, height = 5, bg = "white"
  )
  ggsave(file.path(dir_figs, sprintf("01_FAP_%s.png", suf)),
    plot = p_pacf, width = 8, height = 5, bg = "white"
  )

  # Tres candidatos parsimoniosos cubren las lecturas razonables del
  # correlograma para series financieras: ruido blanco (sin estructura),
  # AR(1) (FAS decae, FAP corta tras lag 1) y ARMA(1,1) (ambos decaen)
  candidatos <- list(
    "ARMA(0,0)" = c(0, 0, 0),
    "ARMA(1,0)" = c(1, 0, 0),
    "ARMA(1,1)" = c(1, 0, 1)
  )
  resumen_arma <- do.call(rbind, lapply(names(candidatos), function(nom) {
    ord <- candidatos[[nom]]
    aju <- tryCatch(
      arima(ret_clean, order = ord, method = "ML",
            optim.control = list(maxit = 1000)),
      error = function(e) NULL
    )
    if (is.null(aju)) {
      data.frame(modelo = nom, AIC = NA_real_, BIC = NA_real_,
                 sigma2 = NA_real_, loglik = NA_real_)
    } else {
      data.frame(modelo = nom, AIC = AIC(aju), BIC = BIC(aju),
                 sigma2 = aju$sigma2, loglik = aju$loglik)
    }
  }))
  resumen_arma <- resumen_arma[order(resumen_arma$AIC), ]
  cat("Comparativa ARMA:\n")
  print(resumen_arma, row.names = FALSE)

  # Parsimonia: si la mejora del mejor candidato sobre el ARMA(0,0) es
  # inferior a 2 unidades de AIC, nos quedamos con el (0,0)
  mejor_aic <- resumen_arma$AIC[1]
  aic_00 <- resumen_arma$AIC[resumen_arma$modelo == "ARMA(0,0)"]
  modelo_seleccionado <- if (length(aic_00) == 1 && (aic_00 - mejor_aic) < 2) {
    "ARMA(0,0)"
  } else {
    as.character(resumen_arma$modelo[1])
  }
  ord_sel <- candidatos[[modelo_seleccionado]]
  ar_p <- ord_sel[1]; ma_q <- ord_sel[3]
  cat(sprintf("ARMA seleccionado: %s\n", modelo_seleccionado))

  # Diagnóstico del ARMA preliminar
  modelo_arma <- arima(ret_clean, order = ord_sel, method = "ML",
                       optim.control = list(maxit = 1000))
  res_arma <- as.numeric(residuals(modelo_arma))

  lb_arma_10 <- Box.test(res_arma, lag = 10, type = "Ljung-Box")
  lb_arma_20 <- Box.test(res_arma, lag = 20, type = "Ljung-Box")
  lf_arma <- lillie.test(res_arma)
  bp_arma <- bptest(res_arma ~ I(seq_along(res_arma)))
  cat(sprintf("LB(10) p=%.3f  LB(20) p=%.3f  Lill p=%.2e  BP p=%.3f\n",
    lb_arma_10$p.value, lb_arma_20$p.value, lf_arma$p.value, bp_arma$p.value))

  # Test ARCH-LM a varios retardos
  arch_pvals <- sapply(c(1, 2, 5, 20), function(k) {
    ArchTest(res_arma, lags = k, demean = TRUE)$p.value
  })
  names(arch_pvals) <- paste0("lag", c(1, 2, 5, 20))
  cat("ARCH-LM:", paste(sprintf("%s=%.2e", names(arch_pvals), arch_pvals),
                        collapse = " | "), "\n")

  # Residuos al cuadrado (visual de los clusters de volatilidad)
  p_r2 <- ggplot(data.frame(t = seq_along(res_arma), r2 = res_arma^2),
                 aes(x = t, y = r2)) +
    geom_line(color = color_main, alpha = 0.7, linewidth = 0.4) +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("Residuos² del ARMA, %s", nombre),
         x = "Índice temporal", y = expression(a[t]^2),
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("01_Residuos2_%s.png", suf)),
    plot = p_r2, width = 8, height = 5, bg = "white"
  )

  # ARMA + GARCH(1,1) con t-Student
  # Estimación conjunta por máxima verosimilitud condicional. La
  # leptocurtosis observada justifica la t-Student frente a la normal
  spec_garch <- ugarchspec(
    variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(ar_p, ma_q), include.mean = TRUE),
    distribution.model = "std"
  )
  ajuste_garch <- ugarchfit(spec = spec_garch, data = ret_clean_xts,
                            solver = "hybrid")
  coefs <- coef(ajuste_garch)
  alpha1 <- as.numeric(coefs["alpha1"])
  beta1  <- as.numeric(coefs["beta1"])
  shape_t <- as.numeric(coefs["shape"])
  persistencia <- alpha1 + beta1
  ic_garch <- infocriteria(ajuste_garch)
  cat(sprintf("alpha1=%.4f  beta1=%.4f  persist=%.4f  nu=%.2f  AIC=%.4f\n",
    alpha1, beta1, persistencia, shape_t, ic_garch[1, 1]))

  # Validación sobre residuos estandarizados
  res_std <- as.numeric(residuals(ajuste_garch, standardize = TRUE))
  lb_std_10  <- Box.test(res_std,   lag = 10, type = "Ljung-Box")
  lb_std_20  <- Box.test(res_std,   lag = 20, type = "Ljung-Box")
  lb_std2_10 <- Box.test(res_std^2, lag = 10, type = "Ljung-Box")
  lb_std2_20 <- Box.test(res_std^2, lag = 20, type = "Ljung-Box")
  arch_std_20 <- ArchTest(res_std, lags = 20, demean = TRUE)$p.value
  cat(sprintf("Residuos std: LB(20) p=%.3f  LB²(20) p=%.3f  ARCH(20) p=%.3f\n",
    lb_std_20$p.value, lb_std2_20$p.value, arch_std_20))

  # Sign Bias Test: detecta asimetría residual frente al sGARCH
  sign_b <- signbias(ajuste_garch)
  sign_neg_p <- sign_b["Negative Sign Bias", "prob"]
  cat(sprintf("Sign Bias Negativo p=%.4f\n", sign_neg_p))

  # Histograma de los residuos estandarizados con la t-Student teórica
  df_rs <- data.frame(z = res_std)
  p_hist <- ggplot(df_rs, aes(x = z)) +
    geom_histogram(aes(y = after_stat(density)), bins = 60,
                   fill = color_fill, color = "black", alpha = 0.6) +
    geom_density(color = color_main, linewidth = 1) +
    stat_function(fun = function(x) dt(x, df = shape_t),
                  color = "black", linewidth = 1, linetype = "dashed") +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("Residuos estandarizados, %s", nombre),
         subtitle = sprintf("Densidad empírica vs t-Student (nu = %.2f)", shape_t),
         x = expression(a[t] / sigma[t]), y = "Densidad",
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("01_Residuos_Std_Hist_%s.png", suf)),
    plot = p_hist, width = 8, height = 5, bg = "white"
  )

  # QQ-Plot frente a la t-Student estimada
  p_qq_std <- ggplot(df_rs, aes(sample = z)) +
    geom_qq(distribution = stats::qt, dparams = list(df = shape_t),
            color = color_main, alpha = 0.5, size = 1.2) +
    geom_qq_line(distribution = stats::qt, dparams = list(df = shape_t),
                 color = "black", linetype = "dashed", linewidth = 0.8) +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("Q-Q Plot residuos estandarizados, %s", nombre),
         subtitle = sprintf("vs t-Student (nu = %.2f)", shape_t),
         x = "Cuantiles teóricos", y = "Cuantiles muestrales",
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("01_QQPlot_Residuos_Std_%s.png", suf)),
    plot = p_qq_std, width = 7, height = 5, bg = "white"
  )

  # Volatilidad condicional histórica y bandas
  vol_hist <- as.numeric(sigma(ajuste_garch))
  df_vol <- data.frame(
    fecha = index(ret_clean_xts),
    ret = ret_clean,
    sup = 2 * vol_hist, inf = -2 * vol_hist
  )
  p_bandas <- ggplot(df_vol, aes(x = fecha)) +
    geom_line(aes(y = ret), color = "gray45", alpha = 0.6, linewidth = 0.3) +
    geom_line(aes(y = sup), color = color_main, linewidth = 0.6) +
    geom_line(aes(y = inf), color = color_main, linewidth = 0.6) +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("Rendimientos vs bandas ±2σ_t, %s", nombre),
         x = "Fecha", y = "Rendimiento / banda",
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("01_Bandas_Volatilidad_%s.png", suf)),
    plot = p_bandas, width = 10, height = 5, bg = "white"
  )

  # Predicción a h_fore jornadas
  fore <- ugarchforecast(ajuste_garch, n.ahead = h_fore)
  mu_pred    <- as.numeric(fitted(fore))
  sigma_pred <- as.numeric(sigma(fore))
  q_975 <- qdist("std", p = 0.975, shape = shape_t)

  p_vol_pred <- ggplot(data.frame(h = 1:h_fore, sigma_pct = sigma_pred * 100),
                       aes(x = h, y = sigma_pct)) +
    geom_line(color = color_main, linewidth = 1.1) +
    geom_point(color = color_main, size = 2) +
    scale_x_continuous(breaks = 1:h_fore) +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("Predicción de volatilidad, %s", nombre),
         x = "Días en adelante", y = expression(sigma[t] ~ "(%)"),
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("01_Prediccion_Volatilidad_%s.png", suf)),
    plot = p_vol_pred, width = 8, height = 5, bg = "white"
  )

  # Reconstrucción del precio: X_{T+h} = X_T * exp( Σ r_{T+j} )
  # Bandas: aplicamos el cuantil t-Student a la suma acumulada de varianzas
  # condicionales predichas
  precio_T   <- as.numeric(xts::last(precio_xts[, "precio"]))
  log_acum   <- cumsum(mu_pred)
  sigma_acum <- sqrt(cumsum(sigma_pred^2))
  precio_pred <- precio_T * exp(log_acum)
  precio_low  <- precio_T * exp(log_acum - q_975 * sigma_acum)
  precio_high <- precio_T * exp(log_acum + q_975 * sigma_acum)
  df_pred <- data.frame(h = 1:h_fore,
                        precio_pred = precio_pred,
                        precio_low  = precio_low,
                        precio_high = precio_high)

  p_pred <- ggplot(df_pred, aes(x = h)) +
    geom_ribbon(aes(ymin = precio_low, ymax = precio_high),
                fill = color_fill, alpha = 0.25) +
    geom_line(aes(y = precio_pred), color = color_main, linewidth = 1.1) +
    geom_point(aes(y = precio_pred), color = color_main, size = 2) +
    scale_x_continuous(breaks = 1:h_fore) +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("Predicción de precio, %s", nombre),
         subtitle = sprintf("ARMA(%d,%d) + GARCH(1,1), t-Student | IC 95%%",
                            ar_p, ma_q),
         x = "Días bursátiles en adelante", y = "Precio",
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("01_Prediccion_Precio_%s.png", suf)),
    plot = p_pred, width = 8, height = 5, bg = "white"
  )

  # Guardado de modelos
  saveRDS(modelo_arma,  file.path(dir_mods, sprintf("arma_%s.rds", suf)))
  saveRDS(ajuste_garch, file.path(dir_mods, sprintf("garch_%s.rds", suf)))

  # Devolución estructurada para tabla comparativa
  list(
    ticker = ticker, nombre = nombre, suf = suf, color_main = color_main,
    n_obs = nrow(precio_xts), n_ao = n_ao, desc = desc,
    arma_orden = ord_sel, coefs = coefs,
    persistencia = persistencia, shape_t = shape_t,
    aic_garch = ic_garch[1, 1], bic_garch = ic_garch[2, 1],
    lb_arma_20_p = lb_arma_20$p.value,
    bp_arma_p = bp_arma$p.value, lf_arma_p = lf_arma$p.value,
    lb_std_20_p = lb_std_20$p.value,
    lb_std2_20_p = lb_std2_20$p.value,
    arch_std_20_p = arch_std_20,
    sign_neg_p = sign_neg_p,
    vol_hist_df = data.frame(
      fecha = as.Date(index(ret_clean_xts)),
      vol_pct = vol_hist * 100
    ),
    pred_retornos = data.frame(
      h = 1:h_fore, mu = mu_pred, sigma = sigma_pred
    ),
    pred_precio = df_pred
  )
}

# Ejecución sobre los seis activos

resultados <- list()
for (a in activos) {
  resultados[[a$suf]] <- tryCatch(
    analizar_activo(a$ticker, a$nombre, a$color_main, a$color_fill, a$suf,
                    fecha_ini, fecha_fin, h_fore, dir_figs, dir_mods),
    error = function(e) {
      message(sprintf("ERROR en %s: %s", a$ticker, conditionMessage(e)))
      NULL
    }
  )
}

# Tabla comparativa univariante

filas <- lapply(resultados, function(r) {
  if (is.null(r)) return(NULL)
  data.frame(
    Activo = r$nombre, Ticker = r$ticker, N_obs = r$n_obs,
    AOs = r$n_ao,
    Skewness = round(r$desc$skewness, 3),
    Kurtosis = round(r$desc$kurtosis, 2),
    ARMA = sprintf("(%d,%d)", r$arma_orden[1], r$arma_orden[3]),
    mu     = signif(as.numeric(r$coefs["mu"]), 3),
    omega  = signif(as.numeric(r$coefs["omega"]), 3),
    alpha1 = round(as.numeric(r$coefs["alpha1"]), 4),
    beta1  = round(as.numeric(r$coefs["beta1"]), 4),
    persistencia = round(r$persistencia, 4),
    nu = round(r$shape_t, 2),
    AIC = round(r$aic_garch, 4),
    LB_std_20_p   = round(r$lb_std_20_p, 4),
    LBsq_std_20_p = round(r$lb_std2_20_p, 4),
    SignBiasNeg_p = round(r$sign_neg_p, 4),
    stringsAsFactors = FALSE
  )
})
df_resumen <- do.call(rbind, filas)
rownames(df_resumen) <- NULL

cat("\nTabla comparativa univariante\n")
print(df_resumen, row.names = FALSE)

write.csv(df_resumen,
  file.path(dir_tabs, "tabla_comparativa_univariante.csv"), row.names = FALSE
)
saveRDS(resultados, file.path(dir_mods, "resultados_univariante.rds"))

# Comparativa visual de volatilidades históricas

df_comp <- do.call(rbind, lapply(resultados, function(r) {
  if (is.null(r)) return(NULL)
  data.frame(
    fecha   = r$vol_hist_df$fecha,
    vol_pct = r$vol_hist_df$vol_pct,
    activo  = factor(r$nombre, levels = sapply(activos, function(a) a$nombre))
  )
}))
paleta <- setNames(
  sapply(activos, function(a) a$color_main),
  sapply(activos, function(a) a$nombre)
)
p_comp <- ggplot(df_comp, aes(x = fecha, y = vol_pct, color = activo)) +
  geom_line(alpha = 0.85, linewidth = 0.5) +
  facet_wrap(~activo, ncol = 2, scales = "free_y") +
  scale_color_manual(values = paleta, guide = "none") +
  theme_minimal(base_size = 12) +
  labs(title = "Volatilidad condicional GARCH(1,1) por activo",
       x = "Fecha", y = expression(sigma[t] ~ "(%)"),
       caption = "Fuente: Elaboración propia.")
ggsave(file.path(dir_figs, "01_Comparativa_Volatilidad_GARCH.png"),
  plot = p_comp, width = 12, height = 8, bg = "white"
)

message("\nPipeline univariante terminado")
