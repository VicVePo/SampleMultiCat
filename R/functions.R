#' Safe ceiling that handles floating-point imprecision
#'
#' Rounds the input to 8 decimal places before applying ceiling(),
#' preventing off-by-one errors caused by IEEE 754 floating-point
#' arithmetic (e.g. ceiling(60 / 0.09999999...) returning 601 instead of 600).
#'
#' @param x Numeric value.
#' @return Integer (ceiling of x after rounding to 8 decimals).
#' @keywords internal
safe_ceiling <- function(x) {
  ceiling(round(x, 8))
}


#' Distribute N into integer counts preserving proportions
#'
#' Internal function to distribute a total N into integer counts
#' for each category, ensuring sum(counts) == N exactly.
#'
#' @param N Total sample size.
#' @param probs Vector of probabilities (must sum to 1).
#'
#' @return Integer vector of counts per category.
#' @keywords internal
distribute_counts <- function(N, probs) {
  raw_counts <- N * probs
  n_floor    <- floor(raw_counts)
  remainder  <- N - sum(n_floor)

  if (remainder > 0) {
    frac      <- raw_counts - n_floor
    idx_order <- order(frac, decreasing = TRUE)
    n_cat     <- n_floor
    n_cat[idx_order[seq_len(remainder)]] <- n_cat[idx_order[seq_len(remainder)]] + 1L
  } else {
    n_cat <- n_floor
  }

  as.integer(n_cat)
}


#' Sample size for proportional odds ordinal logistic regression
#'
#' Calculates sample size for proportional odds ordinal logistic regression
#' with J ordered categories using the worst cumulative cut EPP approach.
#' Optionally includes power-based calculation via Hmisc::posamsize.
#'
#' @param k Integer. Total number of regression parameters in the model
#'   (degrees of freedom consumed by all predictors). For continuous variables,
#'   count 1 per variable. For categorical variables with m levels, count m-1.
#'   For example: 3 continuous + 1 categorical with 4 levels = 3 + 3 = 6.
#' @param probs Numeric vector. Expected probabilities for each ordered category
#'   (length J >= 3). Must sum to 1 (or will be rescaled).
#' @param EPP Numeric. Target events per parameter (default 20).
#' @param OR_min Numeric or NULL. Minimum clinically relevant OR to detect.
#' @param power Numeric. Desired power for OR_min component (default 0.80).
#' @param alpha Numeric. Significance level (default 0.05).
#' @param fraction Numeric. Fraction of subjects in group 1 for Hmisc::posamsize (default 0.5).
#' @param scenarios Logical. If TRUE, shows scenarios for EPP = 10, 20, 30, 40, 50.
#' @param language Character. 'es' (Spanish) or 'en' (English).
#'
#' @return Object of class 'SampleOrdinal' (list).
#'
#' @details
#' The parameter \code{k} must reflect the total degrees of freedom consumed
#' by all predictors, not simply the number of variables. This distinction is
#' critical: a model with 2 continuous variables and 1 categorical variable
#' with 5 levels requires k = 2 + 4 = 6, not k = 3.
#'
#' Sample size calculations use \code{safe_ceiling()} instead of raw
#' \code{ceiling()} to avoid off-by-one errors caused by IEEE 754
#' floating-point imprecision (e.g. \code{1 - 0.9} not being exactly 0.1).
#'
#' @references
#' Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR (1996).
#' A simulation study of the number of events per variable in logistic
#' regression analysis. J Clin Epidemiol 49:1373-1379.
#'
#' Harrell FE (2015). Regression Modeling Strategies. Springer.
#'
#' @examples
#' # Depression scale: none, mild, moderate, severe
#' # Model with 5 continuous + 1 categorical (3 levels) = 5 + 2 = 7 parameters
#' probs_dep <- c(0.50, 0.20, 0.20, 0.10)
#' SampleOrdinal(k = 7, probs = probs_dep, EPP = 20)
#'
#' # With power component
#' SampleOrdinal(k = 7, probs = probs_dep, EPP = 20, OR_min = 1.5)
#'
#' # Multiple scenarios
#' SampleOrdinal(k = 7, probs = probs_dep, scenarios = TRUE)
#'
#' @export
SampleOrdinal <- function(k,
                          probs,
                          EPP = 20,
                          OR_min = NULL,
                          power = 0.8,
                          alpha = 0.05,
                          fraction = 0.5,
                          scenarios = FALSE,
                          language = 'es') {

  if (!language %in% c('es', 'en')) {
    stop('language must be "es" or "en"')
  }

  if (missing(k) || k <= 0 || k != round(k)) {
    stop(ifelse(language == 'es',
                'k debe ser un entero positivo (parametros de regresion, no variables)',
                'k must be a positive integer (regression parameters, not variables)'))
  }

  if (missing(probs)) {
    stop(ifelse(language == 'es',
                'Debe especificar probs (vector de probabilidades)',
                'You must specify probs (vector of probabilities)'))
  }

  if (!is.numeric(probs) || any(probs <= 0)) {
    stop(ifelse(language == 'es',
                'probs debe ser numerico y > 0 en todas las categorias',
                'probs must be numeric and > 0 in all categories'))
  }

  J <- length(probs)
  if (J < 3) {
    stop(ifelse(language == 'es',
                'Se requieren al menos 3 categorias para un modelo ordinal',
                'At least 3 categories are required for an ordinal model'))
  }

  s <- sum(probs)
  if (abs(s - 1) > 1e-6) {
    probs <- probs / s
    warning(ifelse(language == 'es',
                   'probs no sumaba 1. Se ha reescalado para que sume 1.',
                   'probs did not sum to 1. It has been rescaled to sum to 1.'))
  }

  if (EPP <= 0) {
    stop(ifelse(language == 'es',
                'EPP debe ser positivo',
                'EPP must be positive'))
  }

  cum_probs      <- cumsum(probs)[1:(J - 1)]
  event_props    <- cum_probs
  nonevent_props <- 1 - cum_probs
  q_j            <- pmin(event_props, nonevent_props)
  p_min          <- min(q_j)

  if (p_min <= 0) {
    stop(ifelse(language == 'es',
                'Algun corte acumulativo tiene probabilidad cero. Revise probs.',
                'Some cumulative cut has zero probability. Check probs.'))
  }

  N_power_single <- NA_integer_
  power_object   <- NULL

  if (!is.null(OR_min)) {
    if (OR_min <= 0) {
      stop(ifelse(language == 'es',
                  'OR_min debe ser > 0',
                  'OR_min must be > 0'))
    }

    if (!requireNamespace('Hmisc', quietly = TRUE)) {
      stop(ifelse(language == 'es',
                  'Debe instalar el paquete Hmisc para usar OR_min/poder',
                  'You must install the Hmisc package to use OR_min/power'))
    }

    power_object <- Hmisc::posamsize(
      p          = probs,
      odds.ratio = OR_min,
      fraction   = fraction,
      alpha      = alpha,
      power      = power
    )

    N_power_single <- ceiling(power_object$n)
  }

  if (scenarios) {
    EPP_vals  <- c(10, 20, 30, 40, 50)
    N_EPP_vec <- sapply(EPP_vals, function(e) safe_ceiling(e * k / p_min))

    if (!is.null(OR_min)) {
      N_power_vec <- rep(N_power_single, length(EPP_vals))
      N_total_vec <- pmax(N_EPP_vec, N_power_vec)
    } else {
      N_power_vec <- rep(NA_integer_, length(EPP_vals))
      N_total_vec <- N_EPP_vec
    }

    counts <- t(sapply(N_total_vec, function(N) distribute_counts(N, probs)))
    colnames(counts) <- paste0('cat', seq_len(J))

    res <- data.frame(
      EPP     = EPP_vals,
      N_EPP   = N_EPP_vec,
      N_power = N_power_vec,
      N_total = N_total_vec,
      p_min   = rep(p_min, length(EPP_vals)),
      counts,
      row.names = NULL
    )

    if (language == 'es') {
      cat('\nESCENARIOS DE TAMANO MUESTRAL - MODELO ORDINAL (ODDS PROPORCIONALES)\n')
      cat('Categorias (J):', J, '\n')
      cat('Parametros (k):', k, '\n')
      cat('p_min (peor corte acumulativo):', round(p_min, 4), '\n')
      if (!is.null(OR_min)) {
        cat('OR_min para componente de poder:', OR_min, '\n')
        cat('Power objetivo:', power, ' | Alpha:', alpha, '\n')
      }
      cat('\n')
      print(res, row.names = FALSE)
      cat('\nNota: N_EPP controla sobreajuste; N_power controla el poder para OR_min.\n')
      cat('      N_total es el maximo de ambos criterios.\n')
      cat('      Se asume odds proporcionales; si se viola, considerar modelos\n')
      cat('      de odds proporcionales parciales o modelos multinomiales.\n')
      cat('      k = parametros de regresion (grados de libertad), no variables.\n\n')
    } else {
      cat('\nSAMPLE SIZE SCENARIOS - ORDINAL MODEL (PROPORTIONAL ODDS)\n')
      cat('Categories (J):', J, '\n')
      cat('Parameters (k):', k, '\n')
      cat('p_min (worst cumulative cut):', round(p_min, 4), '\n')
      if (!is.null(OR_min)) {
        cat('OR_min for power component:', OR_min, '\n')
        cat('Target power:', power, ' | Alpha:', alpha, '\n')
      }
      cat('\n')
      print(res, row.names = FALSE)
      cat('\nNote: N_EPP controls overfitting; N_power controls power for OR_min.\n')
      cat('      N_total is the maximum of both criteria.\n')
      cat('      Proportional odds are assumed; if violated, consider partial\n')
      cat('      proportional odds or multinomial models.\n')
      cat('      k = regression parameters (degrees of freedom), not variables.\n\n')
    }

    return(invisible(res))
  }

  N_EPP <- safe_ceiling(EPP * k / p_min)
  if (!is.null(OR_min)) {
    N_total <- max(N_EPP, N_power_single)
  } else {
    N_total <- N_EPP
  }

  n_cat        <- distribute_counts(N_total, probs)
  EPP_achieved <- (p_min * N_total) / k

  resultados <- list(
    design        = ifelse(language == 'es',
                           'Modelo con desenlace ordinal (odds proporcionales)',
                           'Model with ordinal outcome (proportional odds)'),
    model         = ifelse(language == 'es',
                           'Regresion logistica ordinal (odds proporcionales)',
                           'Proportional odds ordinal logistic regression'),
    J_categories   = J,
    k_parameters = k,
    probs          = probs,
    EPP_target     = EPP,
    EPP_achieved   = EPP_achieved,
    p_min          = p_min,
    OR_min         = OR_min,
    power          = if (!is.null(OR_min)) power else NA_real_,
    alpha          = if (!is.null(OR_min)) alpha else NA_real_,
    fraction       = if (!is.null(OR_min)) fraction else NA_real_,
    N_EPP          = N_EPP,
    N_power        = N_power_single,
    N_total        = N_total,
    n_by_category  = n_cat,
    power_object   = power_object,
    language       = language
  )
  class(resultados) <- c('SampleOrdinal', 'list')

  if (language == 'es') {
    cat('\n=== TAMANO MUESTRAL - MODELO ORDINAL (ODDS PROPORCIONALES) ===\n')
    cat('Categorias (J):', J, '\n')
    cat('Parametros (k):', k, '\n')
    cat('EPP objetivo:', EPP, '\n')
    cat('p_min (peor corte acumulativo):', round(p_min, 4), '\n')
    cat('Tamano por EPP (N_EPP):', N_EPP, '\n')
    if (!is.null(OR_min)) {
      cat('OR_min:', OR_min, ' | Power:', power, ' | Alpha:', alpha, '\n')
      cat('Tamano por poder (N_power):', N_power_single, '\n')
    }
    cat('----------------------------------------------------------\n')
    cat('>>> TAMANO MUESTRAL RECOMENDADO (N_total):', N_total, '<<<\n')
    cat('EPP efectivo con N_total:', round(EPP_achieved, 2), '\n')
    cat('Frecuencia esperada por categoria:\n')
    for (j in seq_len(J)) {
      cat('  Cat', j, ':', n_cat[j], '\n')
    }
    if (EPP_achieved < 10) {
      cat('\n[X] EPP efectivo < 10: alto riesgo de sobreajuste.\n')
    } else if (EPP_achieved < 20) {
      cat('\n[!] EPP efectivo entre 10 y 20: interpretar con cautela.\n')
    } else if (EPP_achieved < 30) {
      cat('\n[OK] EPP efectivo >= 20: rango aceptable.\n')
    } else {
      cat('\n[OK+] EPP efectivo >= 30: muy robusto.\n')
    }
    if (N_total > 10000) {
      cat('\n[AVISO] N_total muy grande (>', N_total, '). Considere:\n')
      cat('        - Fusionar categorias raras\n')
      cat('        - Reducir el numero de parametros\n')
      cat('        - Aceptar un EPP menor (con cautela)\n')
    }
    cat('\nNota: k = parametros de regresion (grados de libertad), no variables.\n')
    cat('      Continuas: 1 parametro por variable.\n')
    cat('      Categoricas: (numero de niveles - 1) parametros.\n')
    cat('      EPP calculado sobre el peor corte acumulativo.\n')
    cat('      La evidencia para EPP en modelos ordinales es limitada;\n')
    cat('      tomar estos valores como guias pragmaticas.\n\n')
  } else {
    cat('\n=== SAMPLE SIZE - ORDINAL MODEL (PROPORTIONAL ODDS) ===\n')
    cat('Categories (J):', J, '\n')
    cat('Parameters (k):', k, '\n')
    cat('Target EPP:', EPP, '\n')
    cat('p_min (worst cumulative cut):', round(p_min, 4), '\n')
    cat('EPP-based size (N_EPP):', N_EPP, '\n')
    if (!is.null(OR_min)) {
      cat('OR_min:', OR_min, ' | Power:', power, ' | Alpha:', alpha, '\n')
      cat('Power-based size (N_power):', N_power_single, '\n')
    }
    cat('----------------------------------------------------------\n')
    cat('>>> RECOMMENDED TOTAL SAMPLE SIZE (N_total):', N_total, '<<<\n')
    cat('Achieved EPP with N_total:', round(EPP_achieved, 2), '\n')
    cat('Expected counts per category:\n')
    for (j in seq_len(J)) {
      cat('  Cat', j, ':', n_cat[j], '\n')
    }
    if (EPP_achieved < 10) {
      cat('\n[X] Effective EPP < 10: high risk of overfitting.\n')
    } else if (EPP_achieved < 20) {
      cat('\n[!] Effective EPP between 10 and 20: interpret cautiously.\n')
    } else if (EPP_achieved < 30) {
      cat('\n[OK] Effective EPP >= 20: acceptable.\n')
    } else {
      cat('\n[OK+] Effective EPP >= 30: very robust.\n')
    }
    if (N_total > 10000) {
      cat('\n[WARNING] Very large N_total (>', N_total, '). Consider:\n')
      cat('          - Merging rare categories\n')
      cat('          - Reducing the number of parameters\n')
      cat('          - Accepting a lower EPP (with caution)\n')
    }
    cat('\nNote: k = regression parameters (degrees of freedom), not variables.\n')
    cat('      Continuous: 1 parameter per variable.\n')
    cat('      Categorical: (number of levels - 1) parameters.\n')
    cat('      EPP calculated on the worst cumulative cut.\n')
    cat('      Evidence for EPP in ordinal models is limited;\n')
    cat('      treat thresholds as pragmatic guides.\n\n')
  }

  invisible(resultados)
}

#' @export
print.SampleOrdinal <- function(x, ...) {
  if (x$language == 'es') {
    cat('\nRESUMEN - MODELO ORDINAL (ODDS PROPORCIONALES)\n')
    cat('Tamano muestral recomendado (N_total):', x$N_total, '\n')
    cat('  N_EPP   :', x$N_EPP, '\n')
    if (!is.na(x$N_power)) {
      cat('  N_power :', x$N_power, '\n')
    }
    cat('EPP efectivo (peor corte):', round(x$EPP_achieved, 2), '\n')
    cat('Categorias (J):', x$J_categories, '\n')
    cat('Parametros (k):', x$k_parameters, '\n')
  } else {
    cat('\nSUMMARY - ORDINAL MODEL (PROPORTIONAL ODDS)\n')
    cat('Recommended sample size (N_total):', x$N_total, '\n')
    cat('  N_EPP   :', x$N_EPP, '\n')
    if (!is.na(x$N_power)) {
      cat('  N_power :', x$N_power, '\n')
    }
    cat('Effective EPP (worst cut):', round(x$EPP_achieved, 2), '\n')
    cat('Categories (J):', x$J_categories, '\n')
    cat('Parameters (k):', x$k_parameters, '\n')
  }
  invisible(x)
}


#' Sample size for multinomial logistic regression (EPPm approach)
#'
#' Calculates sample size for multinomial logistic regression with K nominal
#' categories using the EPPm approach from de Jong et al.
#'
#' @param k Integer. Total number of regression parameters in the model
#'   (degrees of freedom consumed by all predictors). For continuous variables,
#'   count 1 per variable. For categorical variables with m levels, count m-1.
#'   For example: 3 continuous + 1 categorical with 4 levels = 3 + 3 = 6.
#' @param probs Numeric vector. Expected probabilities for each category (K >= 3).
#' @param EPPm Numeric. Target EPPm (recommended 20-30).
#' @param scenarios Logical. If TRUE, shows scenarios for EPPm = 10, 20, 30, 40, 50.
#' @param language Character. 'es' or 'en'.
#'
#' @return Object of class 'SampleMultinomial' (list).
#'
#' @details
#' The parameter \code{k} must reflect the total degrees of freedom consumed
#' by all predictors, not simply the number of variables. The total number of
#' model parameters is P = (K-1) * k, where K is the number of outcome
#' categories.
#'
#' Sample size calculations use \code{safe_ceiling()} instead of raw
#' \code{ceiling()} to avoid off-by-one errors caused by IEEE 754
#' floating-point imprecision.
#'
#' @references
#' de Jong VMT, Eijkemans MJC, Reitsma JB, et al. (2019).
#' Sample size considerations and predictive performance of multinomial
#' logistic prediction models. Stat Med 38:1601-1619.
#'
#' @examples
#' # Tumor subtypes: 5 nominal categories
#' # Model: 8 continuous + 2 categorical (3 and 4 levels) = 8 + 2 + 3 = 13 parameters
#' probs_tumor <- c(0.40, 0.25, 0.20, 0.10, 0.05)
#' SampleMultinomial(k = 13, probs = probs_tumor, EPPm = 20)
#'
#' # Multiple scenarios
#' SampleMultinomial(k = 13, probs = probs_tumor, scenarios = TRUE)
#'
#' @export
SampleMultinomial <- function(k,
                              probs,
                              EPPm = 20,
                              scenarios = FALSE,
                              language = 'es') {

  if (!language %in% c('es', 'en')) {
    stop('language must be "es" or "en"')
  }

  if (missing(k) || k <= 0 || k != round(k)) {
    stop(ifelse(language == 'es',
                'k debe ser un entero positivo (parametros de regresion, no variables)',
                'k must be a positive integer (regression parameters, not variables)'))
  }

  if (missing(probs)) {
    stop(ifelse(language == 'es',
                'Debe especificar probs (vector de probabilidades)',
                'You must specify probs (vector of probabilities)'))
  }

  if (!is.numeric(probs) || any(probs <= 0)) {
    stop(ifelse(language == 'es',
                'probs debe ser numerico y > 0 en todas las categorias',
                'probs must be numeric and > 0 in all categories'))
  }

  K <- length(probs)
  if (K < 3) {
    stop(ifelse(language == 'es',
                'Se requieren al menos 3 categorias para un modelo multinomial',
                'At least 3 categories are required for a multinomial model'))
  }

  s <- sum(probs)
  if (abs(s - 1) > 1e-6) {
    probs <- probs / s
    warning(ifelse(language == 'es',
                   'probs no sumaba 1. Se ha reescalado para que sume 1.',
                   'probs did not sum to 1. It has been rescaled to sum to 1.'))
  }

  if (EPPm <= 0) {
    stop(ifelse(language == 'es',
                'EPPm debe ser positivo',
                'EPPm must be positive'))
  }

  P     <- (K - 1) * k
  p_min <- min(probs)

  if (p_min <= 0) {
    stop(ifelse(language == 'es',
                'Alguna categoria tiene probabilidad cero. Revise probs.',
                'Some category has zero probability. Check probs.'))
  }

  if (scenarios) {
    EPPm_vals   <- c(10, 20, 30, 40, 50)
    N_total_vec <- sapply(EPPm_vals, function(e) safe_ceiling(e * P / p_min))

    counts <- t(sapply(N_total_vec, function(N) distribute_counts(N, probs)))
    colnames(counts) <- paste0('cat', seq_len(K))

    res <- data.frame(
      EPPm           = EPPm_vals,
      K_categories   = rep(K, length(EPPm_vals)),
      k_parameters = rep(k, length(EPPm_vals)),
      parameters     = rep(P, length(EPPm_vals)),
      p_min          = rep(p_min, length(EPPm_vals)),
      N_total        = N_total_vec,
      counts,
      row.names      = NULL
    )

    if (language == 'es') {
      cat('\nESCENARIOS DE TAMANO MUESTRAL - LOGISTICA MULTINOMIAL (EPPm)\n')
      cat('Categorias (J):', K, '\n')
      cat('Parametros (k):', k, '\n')
      cat('Parametros totales del modelo (P):', P, '\n')
      cat('Proporcion minima esperada (p_min):', round(p_min, 4), '\n\n')
      print(res, row.names = FALSE)
      cat('\nNota: EPPm controla sobreajuste y estabilidad del modelo;\n')
      cat('      no garantiza un poder especifico.\n')
      cat('      k = parametros de regresion (grados de libertad), no variables.\n\n')
    } else {
      cat('\nSAMPLE SIZE SCENARIOS - MULTINOMIAL LOGISTIC (EPPm)\n')
      cat('Categories (J):', K, '\n')
      cat('Parameters (k):', k, '\n')
      cat('Total model parameters (P):', P, '\n')
      cat('Smallest expected proportion (p_min):', round(p_min, 4), '\n\n')
      print(res, row.names = FALSE)
      cat('\nNote: EPPm controls overfitting and model stability;\n')
      cat('      it does not guarantee a specific power.\n')
      cat('      k = regression parameters (degrees of freedom), not variables.\n\n')
    }

    return(invisible(res))
  }

  N_total       <- safe_ceiling(EPPm * P / p_min)
  n_cat         <- distribute_counts(N_total, probs)
  EPPm_achieved <- min(n_cat) / P

  resultados <- list(
    design         = ifelse(language == 'es',
                            'Desenlace multinomial (nominal)',
                            'Multinomial outcome (nominal)'),
    model          = ifelse(language == 'es',
                            'Regresion logistica multinomial',
                            'Multinomial logistic regression'),
    K_categories   = K,
    k_parameters = k,
    parameters     = P,
    probs          = probs,
    EPPm_target    = EPPm,
    EPPm_achieved  = EPPm_achieved,
    p_min          = p_min,
    N_total        = N_total,
    n_by_category  = n_cat,
    language       = language
  )
  class(resultados) <- c('SampleMultinomial', 'list')

  if (language == 'es') {
    cat('\n=== TAMANO MUESTRAL - LOGISTICA MULTINOMIAL (EPPm) ===\n')
    cat('Categorias (J):', K, '\n')
    cat('Parametros (k):', k, '\n')
    cat('Parametros totales del modelo (P):', P, '\n')
    cat('EPPm objetivo:', EPPm, '\n')
    cat('Proporcion minima esperada (p_min):', round(p_min, 4), '\n')
    cat('----------------------------------------------------------\n')
    cat('>>> TAMANO MUESTRAL RECOMENDADO (N_total):', N_total, '<<<\n')
    cat('EPPm efectivo con N_total:', round(EPPm_achieved, 2), '\n')
    cat('Frecuencia esperada por categoria:\n')
    for (j in seq_len(K)) {
      cat('  Cat', j, ':', n_cat[j], '\n')
    }
    if (EPPm_achieved < 10) {
      cat('\n[X] EPPm efectivo < 10: alto riesgo de sobreajuste.\n\n')
    } else if (EPPm_achieved < 20) {
      cat('\n[!] EPPm efectivo entre 10 y 20: interpretar con cautela.\n\n')
    } else if (EPPm_achieved < 30) {
      cat('\n[OK] EPPm efectivo >= 20: rango aceptable.\n\n')
    } else {
      cat('\n[OK+] EPPm efectivo >= 30: muy robusto.\n\n')
    }
    if (N_total > 10000) {
      cat('[AVISO] N_total muy grande (>', N_total, '). Considere:\n')
      cat('        - Fusionar categorias raras\n')
      cat('        - Usar modelo ordinal si tiene sentido clinico\n')
      cat('        - Reducir el numero de parametros\n')
      cat('        - Aceptar un EPPm menor (con cautela)\n\n')
    }
    cat('Nota: k = parametros de regresion (grados de libertad), no variables.\n')
    cat('      Continuas: 1 parametro por variable.\n')
    cat('      Categoricas: (numero de niveles - 1) parametros.\n\n')
  } else {
    cat('\n=== SAMPLE SIZE - MULTINOMIAL LOGISTIC (EPPm) ===\n')
    cat('Categories (J):', K, '\n')
    cat('Parameters (k):', k, '\n')
    cat('Total model parameters (P):', P, '\n')
    cat('Target EPPm:', EPPm, '\n')
    cat('Smallest expected proportion (p_min):', round(p_min, 4), '\n')
    cat('----------------------------------------------------------\n')
    cat('>>> RECOMMENDED TOTAL SAMPLE SIZE (N_total):', N_total, '<<<\n')
    cat('Effective EPPm with N_total:', round(EPPm_achieved, 2), '\n')
    cat('Expected counts per category:\n')
    for (j in seq_len(K)) {
      cat('  Cat', j, ':', n_cat[j], '\n')
    }
    if (EPPm_achieved < 10) {
      cat('\n[X] Effective EPPm < 10: high risk of overfitting.\n\n')
    } else if (EPPm_achieved < 20) {
      cat('\n[!] Effective EPPm between 10 and 20: interpret cautiously.\n\n')
    } else if (EPPm_achieved < 30) {
      cat('\n[OK] Effective EPPm >= 20: acceptable.\n\n')
    } else {
      cat('\n[OK+] Effective EPPm >= 30: very robust.\n\n')
    }
    if (N_total > 10000) {
      cat('[WARNING] Very large N_total (>', N_total, '). Consider:\n')
      cat('          - Merging rare categories\n')
      cat('          - Using an ordinal model if clinically reasonable\n')
      cat('          - Reducing the number of parameters\n')
      cat('          - Accepting a lower EPPm (with caution)\n\n')
    }
    cat('Note: k = regression parameters (degrees of freedom), not variables.\n')
    cat('      Continuous: 1 parameter per variable.\n')
    cat('      Categorical: (number of levels - 1) parameters.\n\n')
  }

  invisible(resultados)
}

#' @export
print.SampleMultinomial <- function(x, ...) {
  if (x$language == 'es') {
    cat('\nRESUMEN - LOGISTICA MULTINOMIAL (EPPm)\n')
    cat('Tamano muestral recomendado (N_total):', x$N_total, '\n')
    cat('EPPm objetivo   :', x$EPPm_target, '\n')
    cat('EPPm efectivo   :', round(x$EPPm_achieved, 2), '\n')
    cat('Categorias (J)  :', x$K_categories, '\n')
    cat('Parametros (k) :', x$k_parameters, '\n')
    cat('Parametros (P)  :', x$parameters, '\n')
  } else {
    cat('\nSUMMARY - MULTINOMIAL LOGISTIC (EPPm)\n')
    cat('Recommended sample size (N_total):', x$N_total, '\n')
    cat('Target EPPm      :', x$EPPm_target, '\n')
    cat('Effective EPPm   :', round(x$EPPm_achieved, 2), '\n')
    cat('Categories (J)   :', x$K_categories, '\n')
    cat('Parameters (k)  :', x$k_parameters, '\n')
    cat('Parameters (P)   :', x$parameters, '\n')
  }
  invisible(x)
}


#' Post-hoc validation of EPP / EPPm
#'
#' Evaluates whether a study with ordinal or multinomial outcome meets
#' reasonable EPP (or EPPm) criteria with the observed sample.
#'
#' @param observed_counts Numeric vector. Observed counts per category.
#' @param k Integer. Total number of regression parameters (degrees of
#'   freedom consumed by all predictors). For continuous variables, count 1
#'   per variable. For categorical variables with m levels, count m-1.
#' @param model Character. 'ordinal' or 'multinomial'.
#' @param language Character. 'es' or 'en'.
#'
#' @return List with metrics and risk level.
#'
#' @examples
#' # Validate an ordinal study (8 regression parameters)
#' ValidateSample(observed_counts = c(500, 200, 150, 50), k = 8, model = 'ordinal')
#'
#' # Validate a multinomial study (10 regression parameters)
#' ValidateSample(observed_counts = c(400, 250, 200, 100, 50), k = 10, model = 'multinomial')
#'
#' @export
ValidateSample <- function(observed_counts,
                           k,
                           model = c('ordinal', 'multinomial'),
                           language = 'es') {

  if (!language %in% c('es', 'en')) {
    stop('language must be "es" or "en"')
  }

  model <- match.arg(model)

  if (missing(k) || k <= 0 || k != round(k)) {
    stop(ifelse(language == 'es',
                'k debe ser un entero positivo (parametros de regresion, no variables)',
                'k must be a positive integer (regression parameters, not variables)'))
  }

  if (missing(observed_counts)) {
    stop(ifelse(language == 'es',
                'Debe especificar observed_counts',
                'You must specify observed_counts'))
  }

  counts <- as.numeric(observed_counts)
  if (any(counts < 0)) {
    stop(ifelse(language == 'es',
                'observed_counts no puede contener valores negativos',
                'observed_counts cannot contain negative values'))
  }

  N_total <- sum(counts)
  if (N_total <= 0) {
    stop(ifelse(language == 'es',
                'La suma de observed_counts debe ser > 0',
                'The sum of observed_counts must be > 0'))
  }

  J <- length(counts)

  risk_level <- function(x) {
    if (x < 10) return('critical')
    if (x < 20) return('caution')
    if (x < 30) return('acceptable')
    return('excellent')
  }

  if (model == 'ordinal') {
    if (J < 3) {
      stop(ifelse(language == 'es',
                  'Para model = "ordinal" se requieren al menos 3 categorias',
                  'For model = "ordinal" at least 3 categories are required'))
    }

    p         <- counts / N_total
    cum_probs <- cumsum(p)[1:(J - 1)]
    q_j       <- pmin(cum_probs, 1 - cum_probs)
    p_min     <- min(q_j)
    EPP_eff   <- (p_min * N_total) / k
    lvl       <- risk_level(EPP_eff)

    res <- list(
      model          = 'ordinal',
      N_total        = N_total,
      J_categories   = J,
      k_parameters = k,
      p_min          = p_min,
      EPP_effective  = EPP_eff,
      risk_level     = lvl,
      counts         = counts,
      language       = language
    )

    if (language == 'es') {
      cat('\n=== VALIDACION POST-HOC - MODELO ORDINAL (EPP) ===\n')
      cat('Muestra total:', N_total, '\n')
      cat('Categorias (J):', J, '\n')
      cat('Parametros (k):', k, '\n')
      cat('EPP efectivo (peor corte acumulativo):', round(EPP_eff, 2), '\n')
      cat('Nivel de riesgo:', lvl, '\n')
      cat('\nInterpretacion:\n')
      cat('  critical   : EPP < 10 (resultados muy inestables)\n')
      cat('  caution    : 10 <= EPP < 20 (interpretar con cautela)\n')
      cat('  acceptable : 20 <= EPP < 30\n')
      cat('  excellent  : EPP >= 30\n')
      cat('\nNota: k = parametros de regresion (grados de libertad), no variables.\n\n')
    } else {
      cat('\n=== POST-HOC VALIDATION - ORDINAL MODEL (EPP) ===\n')
      cat('Total sample:', N_total, '\n')
      cat('Categories (J):', J, '\n')
      cat('Parameters (k):', k, '\n')
      cat('Effective EPP (worst cumulative cut):', round(EPP_eff, 2), '\n')
      cat('Risk level:', lvl, '\n')
      cat('\nInterpretation:\n')
      cat('  critical   : EPP < 10 (very unstable results)\n')
      cat('  caution    : 10 <= EPP < 20 (interpret cautiously)\n')
      cat('  acceptable : 20 <= EPP < 30\n')
      cat('  excellent  : EPP >= 30\n')
      cat('\nNote: k = regression parameters (degrees of freedom), not variables.\n\n')
    }

    return(invisible(res))
  }

  P        <- (J - 1) * k
  n_min    <- min(counts)
  p_min    <- n_min / N_total
  EPPm_eff <- n_min / P
  lvl      <- risk_level(EPPm_eff)

  res <- list(
    model          = 'multinomial',
    N_total        = N_total,
    J_categories   = J,
    k_parameters = k,
    parameters     = P,
    n_min          = n_min,
    p_min          = p_min,
    EPPm_effective = EPPm_eff,
    risk_level     = lvl,
    counts         = counts,
    language       = language
  )

  if (language == 'es') {
    cat('\n=== VALIDACION POST-HOC - MODELO MULTINOMIAL (EPPm) ===\n')
    cat('Muestra total:', N_total, '\n')
    cat('Categorias (J):', J, '\n')
    cat('Parametros (k):', k, '\n')
    cat('Parametros totales del modelo (P):', P, '\n')
    cat('n minimo de categoria:', n_min, '\n')
    cat('EPPm efectivo:', round(EPPm_eff, 2), '\n')
    cat('Nivel de riesgo:', lvl, '\n')
    cat('\nInterpretacion:\n')
    cat('  critical   : EPPm < 10 (resultados muy inestables)\n')
    cat('  caution    : 10 <= EPPm < 20 (interpretar con cautela)\n')
    cat('  acceptable : 20 <= EPPm < 30\n')
    cat('  excellent  : EPPm >= 30\n')
    cat('\nNota: k = parametros de regresion (grados de libertad), no variables.\n\n')
  } else {
    cat('\n=== POST-HOC VALIDATION - MULTINOMIAL MODEL (EPPm) ===\n')
    cat('Total sample:', N_total, '\n')
    cat('Categories (J):', J, '\n')
    cat('Parameters (k):', k, '\n')
    cat('Total model parameters (P):', P, '\n')
    cat('Minimum category count:', n_min, '\n')
    cat('Effective EPPm:', round(EPPm_eff, 2), '\n')
    cat('Risk level:', lvl, '\n')
    cat('\nInterpretation:\n')
    cat('  critical   : EPPm < 10 (very unstable results)\n')
    cat('  caution    : 10 <= EPPm < 20 (interpret cautiously)\n')
    cat('  acceptable : 20 <= EPPm < 30\n')
    cat('  excellent  : EPPm >= 30\n')
    cat('\nNote: k = regression parameters (degrees of freedom), not variables.\n\n')
  }

  invisible(res)
}
