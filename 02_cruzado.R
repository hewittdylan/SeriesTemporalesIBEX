# 02_cruzado.R

# Setup

paquetes <- c(
  "quantmod", "forecast", "tseries", "lmtest", "tsoutliers", "nortest",
  "TSA", "urca", "tsDyn", "xts", "zoo", "tidyverse"
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

ticker_ibex <- "^IBEX"
nombre_ibex <- "IBEX 35"

bancos <- list(
  list(ticker = "SAN.MC",  nombre = "Banco Santander", color = "darkred",    suf = "SAN"),
  list(ticker = "BBVA.MC", nombre = "BBVA",            color = "navy",       suf = "BBVA"),
  list(ticker = "CABK.MC", nombre = "CaixaBank",       color = "darkgreen",  suf = "CABK"),
  list(ticker = "SAB.MC",  nombre = "Banco Sabadell",  color = "darkorange", suf = "SAB"),
  list(ticker = "BKT.MC",  nombre = "Bankinter",       color = "purple",     suf = "BKT")
)

fecha_ini <- "2015-01-01"
fecha_fin <- "2025-12-31"
lag_max_ccf  <- 20      # rango de retardos para la FCC
lags_granger <- c(2, 5, 10)
alpha_sig    <- 0.05
K_johansen   <- 5       # semana bursátil
K_johansen_2 <- 10      # dos semanas

# preparar_datos()
# Misma lógica que en el bloque univariante: precio -> log-rendimientos
# -> limpieza de atípicos aditivos. Mantenerla aquí garantiza que el
# análisis cruzado opera sobre exactamente las mismas series que el
# univariante

preparar_datos <- function(ticker, fecha_ini, fecha_fin) {
  pp <- quantmod::getSymbols(ticker,
    src = "yahoo", from = fecha_ini, to = fecha_fin, auto.assign = FALSE
  )
  precio_xts <- na.omit(quantmod::Ad(pp))
  colnames(precio_xts) <- "precio"
  log_precio_xts <- log(precio_xts)
  colnames(log_precio_xts) <- "logp"

  ret_xts <- diff(log_precio_xts)[-1, ]
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

  list(
    ticker = ticker,
    log_precio_xts = log_precio_xts,
    ret_clean_xts  = ret_clean_xts,
    ret_clean      = ret_clean,
    n_ao = if (is.null(tso_ao$outliers)) 0 else nrow(tso_ao$outliers)
  )
}

# Carga de los seis activos

message("Preparando datos")
datos_ibex <- preparar_datos(ticker_ibex, fecha_ini, fecha_fin)
cat(sprintf("IBEX 35: n=%d  AOs=%d\n",
  nrow(datos_ibex$ret_clean_xts), datos_ibex$n_ao))

datos_bancos <- lapply(bancos, function(b) {
  d <- preparar_datos(b$ticker, fecha_ini, fecha_fin)
  cat(sprintf("%-9s n=%d  AOs=%d\n", b$suf, nrow(d$ret_clean_xts), d$n_ao))
  d
})
names(datos_bancos) <- sapply(bancos, function(b) b$suf)

# Paneles alineados de log-precios y log-rendimientos
logp_list <- c(
  list(IBEX = datos_ibex$log_precio_xts),
  setNames(lapply(datos_bancos, function(d) d$log_precio_xts),
           sapply(bancos, function(b) b$suf))
)
logp_panel <- Reduce(function(x, y) merge(x, y, join = "inner"), logp_list)
colnames(logp_panel) <- c("IBEX", sapply(bancos, function(b) b$suf))

ret_list <- c(
  list(IBEX = datos_ibex$ret_clean_xts),
  setNames(lapply(datos_bancos, function(d) d$ret_clean_xts),
           sapply(bancos, function(b) b$suf))
)
ret_panel <- Reduce(function(x, y) merge(x, y, join = "inner"), ret_list)
colnames(ret_panel) <- c("IBEX", sapply(bancos, function(b) b$suf))

cat(sprintf("Panel alineado: log-precios n=%d  rendimientos n=%d\n",
  nrow(logp_panel), nrow(ret_panel)))

# analizar_par()
# Encapsula todo el bloque bilateral para un par banco-IBEX:
# preblanqueo, FCC, Granger, Box-Tiao (si procede), y cointegración con
# los tres contrastes habituales

analizar_par <- function(banco_info, datos_banco, datos_ibex,
                         lag_max_ccf, lags_granger, alpha_sig, K_johansen,
                         dir_figs, dir_mods) {
  suf    <- banco_info$suf
  nombre <- banco_info$nombre
  color  <- banco_info$color
  message(sprintf("\n--- %s vs IBEX 35 ---", nombre))

  # Alineación por fechas comunes
  ret_panel_par <- merge(datos_ibex$ret_clean_xts,
                         datos_banco$ret_clean_xts, join = "inner")
  colnames(ret_panel_par) <- c("ibex", "banco")
  logp_panel_par <- merge(datos_ibex$log_precio_xts,
                          datos_banco$log_precio_xts, join = "inner")
  colnames(logp_panel_par) <- c("ibex", "banco")

  ret_ibex_v   <- as.numeric(ret_panel_par[, "ibex"])
  ret_banco_v  <- as.numeric(ret_panel_par[, "banco"])
  logp_ibex_v  <- as.numeric(logp_panel_par[, "ibex"])
  logp_banco_v <- as.numeric(logp_panel_par[, "banco"])
  cat(sprintf("Observaciones alineadas: %d\n", length(ret_ibex_v)))

  # Preblanqueo del input
  # IBEX como input por su papel agregado. Como los log-rendimientos del
  # IBEX seleccionaron ARMA(0,0) en el bloque univariante, el preblanqueo
  # se reduce a centrar la serie
  arima_input  <- arima(ret_ibex_v, order = c(0, 0, 0), method = "ML",
                        include.mean = TRUE)
  input_pre    <- as.numeric(residuals(arima_input))
  filtro_out   <- forecast::Arima(ret_banco_v, model = arima_input)
  output_pre   <- as.numeric(residuals(filtro_out))

  # Función de correlación cruzada
  ccf_obj <- ccf(input_pre, output_pre, lag.max = lag_max_ccf, plot = FALSE)
  df_ccf  <- data.frame(lag = as.numeric(ccf_obj$lag),
                        acf = as.numeric(ccf_obj$acf))
  banda   <- 2 / sqrt(length(input_pre))
  i_max   <- which.max(abs(df_ccf$acf))
  ccf_max <- df_ccf$acf[i_max]
  lag_max <- df_ccf$lag[i_max]
  cat(sprintf("FCC máx |rho|=%.4f en lag=%d (banda ±%.4f)\n",
    ccf_max, lag_max, banda))

  p_ccf <- ggplot(df_ccf, aes(x = lag, y = acf)) +
    geom_segment(aes(xend = lag, yend = 0), color = color, linewidth = 0.7) +
    geom_hline(yintercept = c(-banda, banda), linetype = "dashed",
               color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "darkorange", linewidth = 1,
               linetype = "dotted") +
    theme_minimal(base_size = 13) +
    labs(title = sprintf("FCC preblanqueada: IBEX vs %s", nombre),
         subtitle = sprintf("Banda ±2/sqrt(n) = ±%.4f", banda),
         x = "Retardo (días)", y = expression(rho[XY](k)),
         caption = "Fuente: Elaboración propia.")
  ggsave(file.path(dir_figs, sprintf("02_FCC_%s.png", suf)),
    plot = p_ccf, width = 9, height = 5, bg = "white"
  )

  # Causalidad de Granger bidireccional
  resultados_granger <- list()
  for (k in lags_granger) {
    g_in_to_out <- lmtest::grangertest(ret_banco_v ~ ret_ibex_v, order = k)
    g_out_to_in <- lmtest::grangertest(ret_ibex_v ~ ret_banco_v, order = k)
    resultados_granger[[paste0("lag", k)]] <- list(
      ibex_to_banco_p = g_in_to_out$`Pr(>F)`[2],
      banco_to_ibex_p = g_out_to_in$`Pr(>F)`[2]
    )
  }
  for (k in lags_granger) {
    g <- resultados_granger[[paste0("lag", k)]]
    cat(sprintf("Granger lag=%2d  IBEX->banco p=%.4f  banco->IBEX p=%.4f\n",
      k, g$ibex_to_banco_p, g$banco_to_ibex_p))
  }

  # Clasificación a horizonte semanal
  g5 <- resultados_granger[["lag5"]]
  ibex_lidera  <- g5$ibex_to_banco_p < alpha_sig
  banco_lidera <- g5$banco_to_ibex_p < alpha_sig
  clasificacion <- if (ibex_lidera && !banco_lidera) {
    "Unidireccional IBEX -> banco"
  } else if (banco_lidera && !ibex_lidera) {
    "Unidireccional banco -> IBEX"
  } else if (ibex_lidera && banco_lidera) {
    "Feedback bidireccional"
  } else {
    "Sin causalidad significativa"
  }
  cat(sprintf("Clasificación (lag 5): %s\n", clasificacion))

  # Función de transferencia bilateral
  # Solo se ajusta si la causalidad es unidireccional clara. Bajo feedback
  # bidireccional, Box-Tiao queda mal especificado (asume input exógeno)
  modelo_arimax <- NULL
  bt_orden <- NA_character_
  bt_residuos_lb_p <- NA_real_
  bt_omega0 <- NA_real_
  bt_omega0_se <- NA_real_

  if (clasificacion %in% c("Unidireccional IBEX -> banco",
                           "Unidireccional banco -> IBEX")) {
    if (clasificacion == "Unidireccional IBEX -> banco") {
      x_in  <- ret_ibex_v; y_out <- ret_banco_v; direc <- "IBEX -> banco"
    } else {
      x_in  <- ret_banco_v; y_out <- ret_ibex_v; direc <- "banco -> IBEX"
    }
    # Identificación de b a partir del primer lag positivo significativo
    # de la FCC entre input y output preblanqueados
    ccf_id <- ccf(x_in, y_out, lag.max = lag_max_ccf, plot = FALSE)
    df_id  <- data.frame(lag = as.numeric(ccf_id$lag),
                         acf = as.numeric(ccf_id$acf))
    df_pos <- df_id[df_id$lag >= 0, ]
    sig    <- abs(df_pos$acf) > banda
    b_init <- if (any(sig)) df_pos$lag[which(sig)[1]] else 0

    n <- length(x_in)
    x_lag <- if (b_init >= 0 && b_init < n) {
      c(rep(NA, b_init), x_in[1:(n - b_init)])
    } else {
      x_in
    }
    obs_validas <- complete.cases(data.frame(x_lag = x_lag)) & !is.na(y_out)

    modelo_arimax <- tryCatch(
      arima(y_out[obs_validas], order = c(1, 0, 1),
            xreg = data.frame(x_lag = x_lag)[obs_validas, , drop = FALSE],
            method = "ML"),
      error = function(e) NULL
    )
    if (!is.null(modelo_arimax)) {
      bt_orden <- sprintf("ARMA(1,1) + Tx[b=%d, r=0, s=0]  %s",
                          b_init, direc)
      res_bt   <- as.numeric(residuals(modelo_arimax))
      lb_bt    <- Box.test(res_bt, lag = 20, type = "Ljung-Box")
      bt_residuos_lb_p <- lb_bt$p.value
      bt_omega0    <- as.numeric(coef(modelo_arimax)["x_lag"])
      bt_omega0_se <- sqrt(diag(modelo_arimax$var.coef))["x_lag"]
      cat(sprintf("Box-Tiao  omega_0=%.4f (SE=%.4f, t=%.2f)  LB(20) p=%.4f\n",
        bt_omega0, bt_omega0_se,
        bt_omega0 / bt_omega0_se, bt_residuos_lb_p))
      saveRDS(modelo_arimax,
              file.path(dir_mods, sprintf("cruzado_arimax_%s.rds", suf)))
    }
  } else if (clasificacion == "Feedback bidireccional") {
    cat("Feedback bidireccional, no se ajusta Box-Tiao.\n")
  }

  # Cointegración bilateral
  # Engle-Granger: regresión de cointegración + ADF sobre los residuos.
  modelo_eg <- lm(logp_banco_v ~ logp_ibex_v)
  resid_eg  <- as.numeric(residuals(modelo_eg))
  adf_resid <- ur.df(resid_eg, type = "none", selectlags = "AIC")
  tau_eg    <- adf_resid@teststat[1, "tau1"]
  cv_eg_5pct <- -3.34  # Davidson-MacKinnon aprox. con 1 regresor, n grande
  eg_cointegra <- tau_eg < cv_eg_5pct
  cat(sprintf("Engle-Granger  tau=%.3f (cv 5%%: %.3f)  -> %s\n",
    tau_eg, cv_eg_5pct, ifelse(eg_cointegra, "cointegra", "no cointegra")))

  # Phillips-Ouliaris (más robusto bajo autocorrelación residual)
  vec_par <- as.matrix(logp_panel_par)
  vec_par <- vec_par[complete.cases(vec_par), ]
  po_short  <- ca.po(vec_par, demean = "constant", lag = "short")
  po_stat   <- as.numeric(po_short@teststat)
  po_cv5pct <- po_short@cval[1, "5pct"]
  po_cointegra <- po_stat > po_cv5pct
  cat(sprintf("Phillips-Ouliaris  Pu=%.3f (cv 5%%: %.3f) -> %s\n",
    po_stat, po_cv5pct, ifelse(po_cointegra, "cointegra", "no cointegra")))

  # Johansen sobre el par a horizonte de una semana
  jo <- ca.jo(vec_par, type = "trace", ecdet = "const", K = K_johansen)
  jo_stat_r0 <- as.numeric(jo@teststat[2])
  jo_cv5_r0  <- jo@cval[2, "5pct"]
  jo_cointegra <- jo_stat_r0 > jo_cv5_r0
  cat(sprintf("Johansen K=%d (r=0)  stat=%.3f (cv 5%%: %.3f) -> %s\n",
    K_johansen, jo_stat_r0, jo_cv5_r0,
    ifelse(jo_cointegra, "cointegra", "no cointegra")))

  # VECM si Johansen apunta a un vector cointegrante
  modelo_vecm <- NULL
  alpha_ibex  <- NA_real_
  alpha_banco <- NA_real_
  if (jo_cointegra) {
    modelo_vecm <- tryCatch(
      tsDyn::VECM(vec_par, lag = 1, r = 1, estim = "ML",
                  include = "const", LRinclude = "none"),
      error = function(e) NULL
    )
    if (!is.null(modelo_vecm)) {
      vecm_coefs <- coefficients(modelo_vecm)
      alpha_ibex  <- vecm_coefs[1, "ECT"]
      alpha_banco <- vecm_coefs[2, "ECT"]
      cat(sprintf("VECM  alpha_IBEX=%.4f  alpha_banco=%.4f\n",
        alpha_ibex, alpha_banco))
      saveRDS(modelo_vecm,
              file.path(dir_mods, sprintf("cruzado_vecm_%s.rds", suf)))
    }
  }

  list(
    suf = suf, nombre = nombre, n_obs = length(ret_ibex_v),
    ccf_max = ccf_max, lag_ccf_max = lag_max,
    granger = resultados_granger, clasificacion = clasificacion,
    bt_orden = bt_orden, bt_residuos_lb_p = bt_residuos_lb_p,
    bt_omega0 = bt_omega0, bt_omega0_se = bt_omega0_se,
    eg_tau = tau_eg, eg_cointegra = eg_cointegra,
    po_stat = po_stat, po_cv5 = po_cv5pct, po_cointegra = po_cointegra,
    jo_stat = jo_stat_r0, jo_cv5 = jo_cv5_r0, jo_cointegra = jo_cointegra,
    alpha_ibex_vecm = alpha_ibex, alpha_banco_vecm = alpha_banco
  )
}

# Ejecución sobre los cinco pares
resultados_cruzado <- list()
for (b in bancos) {
  resultados_cruzado[[b$suf]] <- tryCatch(
    analizar_par(b, datos_bancos[[b$suf]], datos_ibex,
                 lag_max_ccf, lags_granger, alpha_sig, K_johansen,
                 dir_figs, dir_mods),
    error = function(e) {
      message(sprintf("ERROR en %s: %s", b$suf, conditionMessage(e)))
      NULL
    }
  )
}

# Tabla resumen del bloque bilateral
filas <- lapply(resultados_cruzado, function(r) {
  if (is.null(r)) return(NULL)
  g5  <- r$granger[["lag5"]]
  g10 <- r$granger[["lag10"]]
  data.frame(
    Activo = r$nombre, N_obs = r$n_obs,
    FCC_max = round(r$ccf_max, 4), Lag_FCC_max = r$lag_ccf_max,
    Granger_IBEX_lag5  = signif(g5$ibex_to_banco_p, 3),
    Granger_BANC_lag5  = signif(g5$banco_to_ibex_p, 3),
    Granger_IBEX_lag10 = signif(g10$ibex_to_banco_p, 3),
    Granger_BANC_lag10 = signif(g10$banco_to_ibex_p, 3),
    Clasificacion = r$clasificacion,
    BT_omega0   = ifelse(is.na(r$bt_omega0), NA, round(r$bt_omega0, 4)),
    BT_t_value  = ifelse(is.na(r$bt_omega0), NA,
                         round(r$bt_omega0 / r$bt_omega0_se, 2)),
    BT_LB_p     = ifelse(is.na(r$bt_residuos_lb_p), NA,
                         round(r$bt_residuos_lb_p, 4)),
    EG_tau = round(r$eg_tau, 3), EG_cointegra = r$eg_cointegra,
    PO_Pu  = round(r$po_stat, 3), PO_cointegra = r$po_cointegra,
    JO_stat = round(r$jo_stat, 3), JO_cointegra = r$jo_cointegra,
    Alpha_IBEX_VECM  = ifelse(is.na(r$alpha_ibex_vecm), NA,
                              round(r$alpha_ibex_vecm, 4)),
    Alpha_BANC_VECM  = ifelse(is.na(r$alpha_banco_vecm), NA,
                              round(r$alpha_banco_vecm, 4)),
    stringsAsFactors = FALSE
  )
})
df_resumen_cruzado <- do.call(rbind, filas)
rownames(df_resumen_cruzado) <- NULL

cat("\nTabla bilateral\n")
print(df_resumen_cruzado, row.names = FALSE)

write.csv(df_resumen_cruzado,
  file.path(dir_tabs, "tabla_cruzado.csv"), row.names = FALSE)
saveRDS(resultados_cruzado,
  file.path(dir_mods, "resultados_cruzado.rds"))

# Visualización conjunta de los Granger en escala -log10(p)
df_g_long <- do.call(rbind, lapply(resultados_cruzado, function(r) {
  if (is.null(r)) return(NULL)
  do.call(rbind, lapply(names(r$granger), function(nom) {
    g <- r$granger[[nom]]
    data.frame(
      Activo = r$nombre, Lag = as.integer(sub("lag", "", nom)),
      Direccion = c("IBEX -> banco", "banco -> IBEX"),
      p_valor = c(g$ibex_to_banco_p, g$banco_to_ibex_p)
    )
  }))
}))
p_granger <- ggplot(df_g_long,
                    aes(x = factor(Lag), y = -log10(p_valor), fill = Direccion)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = -log10(alpha_sig), linetype = "dashed",
             color = "red", linewidth = 0.5) +
  facet_wrap(~Activo, ncol = 3) +
  scale_fill_manual(values = c("IBEX -> banco" = "steelblue",
                               "banco -> IBEX" = "tomato")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(title = "Causalidad de Granger bidireccional",
       subtitle = expression("-log"[10]*"(p)  |  línea roja: significación 5%"),
       x = "Lag", y = expression("-log"[10]*"(p)"),
       fill = "Dirección",
       caption = "Fuente: Elaboración propia.")
ggsave(file.path(dir_figs, "02_Comparativa_Granger.png"),
  plot = p_granger, width = 11, height = 7, bg = "white"
)

message("\nCointegración multivariante")

vec_multi <- as.matrix(logp_panel[, c("IBEX", sapply(bancos, function(b) b$suf))])
vec_multi <- vec_multi[complete.cases(vec_multi), ]

jo_multi_5  <- ca.jo(vec_multi, type = "trace", ecdet = "const",
                     K = K_johansen)
jo_multi_10 <- ca.jo(vec_multi, type = "trace", ecdet = "const",
                     K = K_johansen_2)

cat("Johansen K=5 (seis series):\n");  print(summary(jo_multi_5))
cat("Johansen K=10 (seis series):\n"); print(summary(jo_multi_10))

# Procedimiento secuencial: la primera nula no rechazada fija el rango
# En ca.jo, teststat[p] corresponde a H0: r = 0 (más restrictiva) y
# teststat[1] a H0: r <= p-1.
determinar_rango <- function(jo) {
  stats <- as.numeric(jo@teststat)
  cvs5  <- jo@cval[, "5pct"]
  p     <- length(stats)
  rango <- 0
  for (j in 0:(p - 1)) {
    i <- p - j
    if (stats[i] > cvs5[i]) rango <- j + 1 else break
  }
  rango
}
rango_K5  <- determinar_rango(jo_multi_5)
rango_K10 <- determinar_rango(jo_multi_10)
cat(sprintf("Rango (K=5): r = %d\nRango (K=10): r = %d\n",
  rango_K5, rango_K10))

saveRDS(jo_multi_5,  file.path(dir_mods, "johansen_multi_K5.rds"))
saveRDS(jo_multi_10, file.path(dir_mods, "johansen_multi_K10.rds"))

message("\nTransferencia multi-input")

# Los hallazgos del bloque bilateral (FCC concentrada en lag 0) justifican
# fijar b_i = 0 para todos los inputs. El modelo se reduce a una regresión
# contemporánea con ruido ARMA(p,q):
# IBEX_t = beta_0 + sum_i beta_i * banco_i,t + N_t,  N_t ~ ARMA(p,q)
y_ibex   <- as.numeric(ret_panel[, "IBEX"])
X_bancos <- as.matrix(ret_panel[, sapply(bancos, function(b) b$suf)])

candidatos_tf <- list(
  "TF-Multi[ARMA(0,0)]" = c(0, 0, 0),
  "TF-Multi[ARMA(1,0)]" = c(1, 0, 0),
  "TF-Multi[ARMA(0,1)]" = c(0, 0, 1),
  "TF-Multi[ARMA(1,1)]" = c(1, 0, 1)
)
resumen_tf <- do.call(rbind, lapply(names(candidatos_tf), function(nom) {
  ord <- candidatos_tf[[nom]]
  aju <- tryCatch(
    arima(y_ibex, order = ord, xreg = X_bancos, method = "ML",
          optim.control = list(maxit = 1000)),
    error = function(e) NULL
  )
  if (is.null(aju)) {
    data.frame(modelo = nom, AIC = NA_real_, BIC = NA_real_, sigma2 = NA_real_)
  } else {
    data.frame(modelo = nom, AIC = AIC(aju), BIC = BIC(aju), sigma2 = aju$sigma2)
  }
}))
resumen_tf <- resumen_tf[order(resumen_tf$AIC), ]
cat("Comparativa de candidatos:\n"); print(resumen_tf, row.names = FALSE)

mejor_aic <- resumen_tf$AIC[1]
aic_simple <- resumen_tf$AIC[resumen_tf$modelo == "TF-Multi[ARMA(0,0)]"]
modelo_tf_sel <- if (length(aic_simple) == 1 && (aic_simple - mejor_aic) < 2) {
  "TF-Multi[ARMA(0,0)]"
} else {
  as.character(resumen_tf$modelo[1])
}
ord_tf <- candidatos_tf[[modelo_tf_sel]]
cat(sprintf("Modelo seleccionado: %s\n", modelo_tf_sel))

modelo_tf <- arima(y_ibex, order = ord_tf, xreg = X_bancos,
                   method = "ML", optim.control = list(maxit = 1000))
cat("\nCoeftest:\n"); print(coeftest(modelo_tf))

# Diagnóstico de residuos
res_tf  <- as.numeric(residuals(modelo_tf))
lb_tf_10 <- Box.test(res_tf, lag = 10, type = "Ljung-Box")
lb_tf_20 <- Box.test(res_tf, lag = 20, type = "Ljung-Box")
lf_tf    <- nortest::lillie.test(res_tf)
bp_tf    <- bptest(res_tf ~ I(seq_along(res_tf)))
cat(sprintf("LB(10) p=%.4f  LB(20) p=%.4f  Lill p=%.2e  BP p=%.4f\n",
  lb_tf_10$p.value, lb_tf_20$p.value, lf_tf$p.value, bp_tf$p.value))

# Coeficientes por banco con errores estándar asintóticos
coefs_tf  <- coef(modelo_tf)
ses_tf    <- sqrt(diag(modelo_tf$var.coef))
suf_bancos <- sapply(bancos, function(b) b$suf)
beta_b    <- coefs_tf[suf_bancos]
beta_se   <- ses_tf[suf_bancos]

df_betas <- data.frame(
  Banco   = sapply(bancos, function(b) b$nombre),
  Beta    = round(as.numeric(beta_b), 4),
  SE      = round(as.numeric(beta_se), 4),
  t_value = round(as.numeric(beta_b) / as.numeric(beta_se), 3),
  p_valor = signif(2 * pnorm(-abs(as.numeric(beta_b) / as.numeric(beta_se))), 3),
  stringsAsFactors = FALSE
)
df_betas$Suma_beta_total <- sum(df_betas$Beta)
cat("\nPesos efectivos por banco:\n"); print(df_betas, row.names = FALSE)
cat(sprintf("Suma beta_i = %.4f (elasticidad agregada)\n", sum(df_betas$Beta)))

write.csv(df_betas, file.path(dir_tabs, "tabla_transfer_multiinput.csv"),
          row.names = FALSE)
saveRDS(modelo_tf, file.path(dir_mods, "transfer_multiinput_ibex.rds"))

# Gráfico de pesos con IC al 95%
df_plot <- df_betas
df_plot$Banco_factor <- factor(df_plot$Banco,
                               levels = sapply(bancos, function(b) b$nombre))
df_plot$IC_low  <- df_plot$Beta - 1.96 * df_plot$SE
df_plot$IC_high <- df_plot$Beta + 1.96 * df_plot$SE
p_betas <- ggplot(df_plot, aes(x = Banco_factor, y = Beta)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = IC_low, ymax = IC_high),
                width = 0.2, color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(title = "Pesos efectivos de los bancos en el IBEX 35",
       subtitle = sprintf("Coeficientes beta_i con IC 95%%  |  suma = %.3f",
                          sum(df_plot$Beta)),
       x = NULL, y = expression(beta[i]),
       caption = "Fuente: Elaboración propia.") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(dir_figs, "02_Pesos_Bancos_IBEX.png"),
  plot = p_betas, width = 9, height = 6, bg = "white"
)

message("\nAnálisis cruzado terminado")
