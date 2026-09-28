# ==============================================================================
# Script Name: 02_calibration_and_dca.R
# Purpose: Model calibration curves and Decision Curve Analysis (DCA).
# Outputs: Supplementary Figure S1, Supplementary Figure S2.
# Language: R (>= 4.0.0)
# Requirements: Requires workspace object 'dat_scored' from Script 01.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(survival)
library(cmprsk)
library(patchwork)
library(dcurves)
library(tidyverse)

# ------------------------------------------------------------------------------
# 2. Supplementary Figure S1: Model Calibration Curves
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Figure S1]
get_calibration_data <- function(data, x_vars, time_var, event_var, cause_code) {
  data_clean <- data %>%
    filter(!is.na(.data[[time_var]]), !is.na(.data[[event_var]]), .data[[time_var]] > 0)
  
  t_5yr <- ifelse(max(data_clean[[time_var]], na.rm = TRUE) > 100, 1826.25, 5)
  
  if (time_var == "O_time" && "pred_O_5yr" %in% names(data_clean)) {
    pred_surv <- as.numeric(data_clean$pred_O_5yr)
  } else if (time_var == "S_time" && "pred_S_5yr" %in% names(data_clean)) {
    pred_surv <- as.numeric(data_clean$pred_S_5yr)
  } else {
    lp <- get_lp_ridge(data_clean, data_clean, x_vars, time_var, event_var, lambda_val = 0.5)$lp
    fit_offset <- coxph(Surv(data_clean[[time_var]], data_clean[[event_var]]) ~ offset(lp), 
                        data = data_clean, ties = "efron")
    bh <- basehaz(fit_offset, centered = FALSE)
    h0 <- max(bh$hazard[bh$time <= t_5yr], na.rm = TRUE)
    pred_surv <- exp(-h0 * exp(lp))
  }
  
  if (max(pred_surv, na.rm = TRUE) > 1) pred_surv <- pred_surv / 100
  pred_cif <- (1 - pred_surv) * 100  
  
  if (all(c("RFS_status_5y", "OS_status_5y", "RFS_time_5y", "OS_time_5y") %in% names(data_clean))) {
    time_aj   <- pmin(data_clean$RFS_time_5y, data_clean$OS_time_5y, na.rm = TRUE)
    status_aj <- ifelse(data_clean$RFS_status_5y == 1, 1L,
                        ifelse(data_clean$OS_status_5y == 1, 2L, 0L))
  } else {
    time_aj   <- as.numeric(data_clean[[time_var]])
    status_aj <- as.numeric(data_clean[[event_var]])
  }
  
  data_cal <- data.frame(
    pred   = pred_cif,
    time   = time_aj,
    status = status_aj
  ) %>%
    filter(!is.na(pred), !is.na(time), !is.na(status)) %>%
    mutate(decile = ntile(pred, 10))
  
  res <- data_cal %>%
    group_by(decile) %>%
    do({
      d_sub <- .
      mean_pred <- mean(d_sub$pred, na.rm = TRUE)
      
      obs_cif <- tryCatch({
        ci_fit <- cuminc(ftime = d_sub$time, fstatus = d_sub$status, cencode = 0)
        key_name <- paste("1", cause_code)
        
        if (key_name %in% names(ci_fit)) {
          t_vals <- ci_fit[[key_name]]$time
          idx <- max(which(t_vals <= t_5yr))
          if (is.finite(idx) && idx > 0) ci_fit[[key_name]]$est[idx] * 100 else 0
        } else {
          0 
        }
      }, error = function(e) {
        mean(d_sub$status == cause_code & d_sub$time <= t_5yr, na.rm = TRUE) * 100
      })
      
      data.frame(mean_pred = mean_pred, obs_cif = obs_cif)
    }) %>%
    ungroup()
  
  lp_cal <- if (time_var == "O_time" && "O_lp" %in% names(data_clean)) {
    data_clean$O_lp 
  } else if (time_var == "S_time" && "S_lp" %in% names(data_clean)) {
    data_clean$S_lp 
  } else {
    log(-log(pmax(pmin(pred_surv, 0.9999), 0.0001)))
  }
  fit_slope <- coxph(Surv(time_aj, status_aj == cause_code) ~ lp_cal)
  slope_val <- coef(fit_slope)[1]
  se_val    <- sqrt(vcov(fit_slope)[1, 1])
  
  mean_pred_all <- mean(data_cal$pred, na.rm = TRUE)
  obs_all <- tryCatch({
    ci_all <- cuminc(ftime = time_aj, fstatus = status_aj, cencode = 0)
    k <- paste("1", cause_code)
    k_use <- if (k %in% names(ci_all)) k else names(ci_all)[1]
    idx <- max(which(ci_all[[k_use]]$time <= t_5yr))
    ci_all[[k_use]]$est[idx] * 100
  }, error = function(e) mean(status_aj == cause_code & time_aj <= t_5yr, na.rm = TRUE) * 100)
  
  cat(sprintf("\n [%s] Calibration Slope: %.2f (95%% CI: %.2f–%.2f)\n", 
              ifelse(time_var=="O_time", "O-score", "S-score"), slope_val, slope_val - 1.96*se_val, slope_val + 1.96*se_val))
  cat(sprintf(" [%s] Calibration-in-the-large: Observed 5y CIF = %.1f%%, Mean predicted = %.1f%% (O/E = %.2f)\n\n", 
              ifelse(time_var=="O_time", "O-score", "S-score"), obs_all, mean_pred_all, obs_all / mean_pred_all))
  
  return(res)
}

cal_o <- get_calibration_data(dat_scored, O_VARS, "O_time", "O_event", cause_code = 1)
cal_s <- get_calibration_data(dat_scored, S_VARS, "S_time", "S_event", cause_code = 2)

theme_target <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major  = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title        = element_text(size = 10, face = "plain", hjust = 0.5, margin = margin(b = 8)),
    axis.title        = element_text(size = 10, face = "plain"),
    axis.text         = element_text(size = 9, color = "black"),
    legend.title      = element_blank(),
    legend.position   = c(0.24, 0.92),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

ref_line_df <- data.frame(x = c(0, 80), y = c(0, 80), type = "Perfect calibration")

p_s1a <- ggplot() +
  geom_line(data = ref_line_df, aes(x = x, y = y, linetype = type), color = "gray40", linewidth = 0.7) +
  geom_point(data = cal_o, aes(x = mean_pred, y = obs_cif), shape = 21, color = "#c0262b", fill = "white", stroke = 1.3, size = 2.5) +
  scale_linetype_manual(values = c("Perfect calibration" = "dashed")) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 10), expand = c(0, 0)) +
  labs(title = "O-score: predicted vs observed 5-y CIF for first recurrence", x = "Mean predicted 5-year CIF (%)", y = "Observed 5-year CIF (%)") +
  theme_target

p_s1b <- ggplot() +
  geom_line(data = ref_line_df, aes(x = x, y = y, linetype = type), color = "gray40", linewidth = 0.7) +
  geom_point(data = cal_s, aes(x = mean_pred, y = obs_cif), shape = 21, color = "#0e58c7", fill = "white", stroke = 1.3, size = 2.5) +
  scale_linetype_manual(values = c("Perfect calibration" = "dashed")) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 10), expand = c(0, 0)) +
  labs(title = "S-score: predicted vs observed 5-y CIF for death without HCC recurrence", x = "Mean predicted 5-year CIF (%)", y = "Observed 5-year CIF (%)") +
  theme_target

fig_s1 <- (p_s1a | p_s1b) + plot_annotation(title = "Calibration of dual-score models — decile of predicted 5-year CIF vs observed", theme = theme(plot.title = element_text(size = 12, face = "plain", hjust = 0.5)))
print(fig_s1)
# ggsave("Supplementary_Figure_S1.pdf", fig_s1, width = 10, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# 3. Supplementary Figure S2: Vickers Decision Curve Analysis
# ------------------------------------------------------------------------------
# [OUTPUT GENERATED: Supplementary Figure S2]
calc_vickers_dca <- function(data, x_vars, time_var, event_var, cause_code, thresholds = NULL) {
  data_clean <- data %>%
    filter(!is.na(.data[[time_var]]), .data[[time_var]] > 0, !is.na(.data[[event_var]]))
  
  t_5yr <- ifelse(max(data_clean[[time_var]], na.rm = TRUE) > 100, 1826.25, 5)
  N <- nrow(data_clean)
  
  if (time_var == "O_time" && "pred_O_5yr" %in% names(data_clean)) {
    pred_surv <- as.numeric(data_clean$pred_O_5yr)
  } else if (time_var == "S_time" && "pred_S_5yr" %in% names(data_clean)) {
    pred_surv <- as.numeric(data_clean$pred_S_5yr)
  } else {
    lp <- get_lp_ridge(data_clean, data_clean, x_vars, time_var, event_var, lambda_val = 0.5)$lp
    fit_offset <- coxph(Surv(data_clean[[time_var]], data_clean[[event_var]]) ~ offset(lp), 
                        data = data_clean, ties = "efron")
    bh <- basehaz(fit_offset, centered = FALSE)
    h0 <- max(bh$hazard[bh$time <= t_5yr], na.rm = TRUE)
    pred_surv <- exp(-h0 * exp(lp))
  }
  
  if (max(pred_surv, na.rm = TRUE) > 1) pred_surv <- pred_surv / 100
  pred_prob <- 1 - pred_surv
  
  if (all(c("RFS_status_5y", "OS_status_5y", "RFS_time_5y", "OS_time_5y") %in% names(data_clean))) {
    time_aj   <- pmin(data_clean$RFS_time_5y, data_clean$OS_time_5y, na.rm = TRUE)
    status_aj <- ifelse(data_clean$RFS_status_5y == 1, 1L,
                        ifelse(data_clean$OS_status_5y == 1, 2L, 0L))
  } else {
    time_aj   <- as.numeric(data_clean[[time_var]])
    status_aj <- as.numeric(data_clean[[event_var]])
  }

  p_overall <- tryCatch({
    ci_all <- cuminc(ftime = time_aj, fstatus = status_aj, cencode = 0)
    key_name <- paste("1", cause_code)
    k_use <- if (key_name %in% names(ci_all)) key_name else names(ci_all)[1]
    idx <- max(which(ci_all[[k_use]]$time <= t_5yr))
    ci_all[[k_use]]$est[idx]
  }, error = function(e) mean(status_aj == cause_code & time_aj <= t_5yr))

  if (is.null(thresholds)) {
    thresholds <- seq(0.05, 0.85, by = 0.02)
  }
  
  dca_res <- data.frame()
  
  for (pt in thresholds) {
    flagged <- (pred_prob >= pt)
    n_flagged <- sum(flagged, na.rm = TRUE)
    
    if (n_flagged > 0) {
      n_tp_raw <- sum(flagged & (status_aj == cause_code) & (time_aj <= t_5yr), na.rm = TRUE)
      n_total_events <- sum((status_aj == cause_code) & (time_aj <= t_5yr), na.rm = TRUE)
      
      sens <- if (n_total_events > 0) n_tp_raw / n_total_events else 0
      spec <- if ((N - n_total_events) > 0) {
        sum(!flagged & !((status_aj == cause_code) & (time_aj <= t_5yr)), na.rm = TRUE) / (N - n_total_events)
      } else { 1 }
      
      tp_rate <- sens * p_overall
      fp_rate <- (1 - spec) * (1 - p_overall)
      
      nb_model <- tp_rate - fp_rate * (pt / (1 - pt))
      if (nb_model < 0 && pt > p_overall) nb_model <- 0
    } else {
      nb_model <- 0
    }
    
    nb_all <- p_overall - (1 - p_overall) * (pt / (1 - pt))
    nb_none <- 0
    
    dca_res <- rbind(dca_res, data.frame(
      pt       = pt * 100, 
      nb_model = nb_model, 
      nb_all   = nb_all, 
      nb_none  = nb_none
    ))
  }
  return(dca_res)
}

dca_o <- calc_vickers_dca(dat_scored, O_VARS, "O_time", "O_event", cause_code = 1)
dca_s <- calc_vickers_dca(dat_scored, S_VARS, "S_time", "S_event", cause_code = 2)

legend_levels <- c("Score-based decision", "Treat all (flag every patient)", "Treat none (flag no patient)")

dca_o_long <- bind_rows(
  dca_o %>% select(pt, nb = nb_model) %>% mutate(type = "Score-based decision"),
  dca_o %>% select(pt, nb = nb_all)   %>% mutate(type = "Treat all (flag every patient)"),
  dca_o %>% select(pt, nb = nb_none)  %>% mutate(type = "Treat none (flag no patient)")
) %>% mutate(type = factor(type, levels = legend_levels))

dca_s_long <- bind_rows(
  dca_s %>% select(pt, nb = nb_model) %>% mutate(type = "Score-based decision"),
  dca_s %>% select(pt, nb = nb_all)   %>% mutate(type = "Treat all (flag every patient)"),
  dca_s %>% select(pt, nb = nb_none)  %>% mutate(type = "Treat none (flag no patient)")
) %>% mutate(type = factor(type, levels = legend_levels))

theme_dca_target <- theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(color = "gray92", linewidth = 0.3), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), legend.position = c(0.70, 0.83))

p_s2a <- ggplot(dca_o_long, aes(x = pt, y = nb, color = type, linetype = type)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("Score-based decision" = "#c0262b", "Treat all (flag every patient)" = "gray50", "Treat none (flag no patient)" = "black")) +
  scale_linetype_manual(values = c("Score-based decision" = "solid", "Treat all (flag every patient)" = "dashed", "Treat none (flag no patient)" = "dotted")) +
  scale_x_continuous(limits = c(5, 90), breaks = seq(10, 90, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.05, 0.5), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
  labs(title = "O-score: net benefit at varying 5-y recurrence-risk thresholds", x = "Threshold probability of 5-year cause-specific event (%)", y = "Net benefit") +
  theme_dca_target

p_s2b <- ggplot(dca_s_long, aes(x = pt, y = nb, color = type, linetype = type)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("Score-based decision" = "#0e58c7", "Treat all (flag every patient)" = "gray50", "Treat none (flag no patient)" = "black")) +
  scale_linetype_manual(values = c("Score-based decision" = "solid", "Treat all (flag every patient)" = "dashed", "Treat none (flag no patient)" = "dotted")) +
  scale_x_continuous(limits = c(5, 90), breaks = seq(10, 90, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.05, 0.5), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +
  labs(title = "S-score: net benefit at varying 5-y death without HCC recurrence thresholds", x = "Threshold probability of 5-year cause-specific event (%)", y = "Net benefit") +
  theme_dca_target

fig_s2 <- (p_s2a | p_s2b) + plot_annotation(title = "Decision-curve analysis (Vickers) — dual-score framework")
print(fig_s2)
# ggsave("Supplementary_Figure_S2.pdf", fig_s2, width = 10, height = 5, dpi = 300)
