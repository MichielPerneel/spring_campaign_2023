#=========================== carbonate_gamm_analysis.R =========================
# GAMM analysis for O2 and DIC at stations 51 & 130
# Instead of O2 concentration, we model O2' = O2 saturation anomaly (µmol/kg)
# This is to remove physical solubility effects, enriching the biological signal,
# and to make it more easily interpretable in terms of photosynthesis vs respiration
# (positive O2' = supersaturation = photosynthesis,
#  negative O2' = undersaturation = respiration)
#
# Abiotic predictors: Salinity and wind speed
#
# Check https://ecogambler.netlify.app/blog/interpreting-gams/ for GAM interpretation
#
# Outputs:
#  - Model selection tables (AIC, CV-RMSE) per station/response
#  - Model choice per station (min mean CV-RMSE across O2 & DIC)
#  - Stoichiometric slope from residuals (ΔDIC/ΔO2), residual correlation
#  - 3-panel plots (observed vs predicted; abiotic explanatory vars, residuals & d/dt)
#===============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(purrr)
library(data.table)
library(seacarb)
library(mgcv)
library(gratia)
library(suncalc)
library(cowplot)
library(patchwork)
library(scales)
library(oce)
library(Smisc)
library(ggpubr)

Sys.setenv(TZ = "UTC")
theme_set(theme_minimal())

# ------------------ Some helper functions -------------------------------------
safe_time <- function(x) as.POSIXct(x, tz = "UTC")

make_time_blocks <- function(n, k = 5) {
  # Create k contiguous time blocks for CV
  # Returns a list of index vectors defining contiguous folds
  stopifnot(n >= k)
  fold_size <- floor(n/k)
  starts <- seq(1, by = fold_size, length.out = k)
  ends   <- c(starts[-1] - 1, n)
  Map(seq.int, starts, ends)
}

match_oxygen <- function(data, oxygen_data) {
  # Match nearest oxygen data to data by timestamp (rolling join)
  oxygen_dt <- as.data.table(oxygen_data)
  data_dt   <- as.data.table(data)
  oxygen_dt[, Date := safe_time(Date)]
  data_dt[, Date.Time := safe_time(Date.Time)]
  setnames(data_dt, "Date.Time", "Date")
  setkey(oxygen_dt, Date)
  matched <- oxygen_dt[data_dt, roll = "nearest", on = "Date", nomatch = 0]
  data_dt$O2uM <- matched$O2uM
  as_tibble(data_dt)
}

interp_to_continuous <- function(df_cont, df_sparse, time_col = "Date", value_col = "value") {
  # linear interpolation of sparse data to continuous timestamps
  # Can be used to align spot samples to continuous time series
  approx(x = as.numeric(df_sparse[[time_col]]),
         y = df_sparse[[value_col]],
         xout = as.numeric(df_cont[[time_col]]),
         rule = 2)$y
}

# Sunrise/sunset shading
get_station_sun_times <- function(lat, lon, dates) {
  # Based on time and location, get sun times for shading
  getSunlightTimes(date = dates, lat = lat, lon = lon,
                   keep = c("sunrise","sunset","dawn","dusk","nauticalDawn","nauticalDusk","night","nightEnd"),
                   tz = "UTC") |>
    mutate(date = as.Date(date)) |>
    dplyr::select(date, night, nightEnd, nauticalDawn, nauticalDusk, dawn, dusk, sunrise, sunset)
}

add_day_moment <- function(df, sun_times) {
  # Add day_moment column to df based on sun_times
  df %>%
    mutate(date = as.Date(Date)) %>%
    left_join(sun_times, by = "date") %>%
    rowwise() %>%
    mutate(day_moment = case_when(
      is.na(night) ~ NA_character_,
      Date < nightEnd ~ "Night",
      Date < nauticalDawn ~ "Astronomical twilight",
      Date < dawn ~ "Nautical twilight",
      Date < sunrise ~ "Civil twilight",
      Date < sunset ~ "Day",
      Date < dusk ~ "Civil twilight",
      Date < nauticalDusk ~ "Nautical twilight",
      Date < night ~ "Astronomical twilight",
      TRUE ~ "Night"
    )) %>%
    ungroup() %>% dplyr::select(-date,-night,-nightEnd,-nauticalDawn,-nauticalDusk,-dawn,-dusk,-sunrise,-sunset)
}

make_model_df <- function(df) {
  df %>%
    arrange(Date) %>%
    mutate(
      Time_numeric    = as.numeric(Date - min(Date)),
      Salinity_scaled = scale(Salinity)[,1],
      Temp_scaled     = scale(Temp)[,1],
      Wind_scaled     = scale(Wind.Speed)[,1],
      gap_start       = c(TRUE, diff(Date) > 1800),
      block           = cumsum(gap_start)
    )
}

fit_gamm <- function(response, rhs, dat) {
  form <- reformulate(termlabels = as.character(rhs)[-1], response = response)
  mgcv::gamm(
    formula = form,
    data = dat,
    method = "REML",
    correlation = nlme::corCAR1(form = ~ Time_numeric | block)
  )
}

cv_rmse <- function(response, rhs, dat, folds) {
  # Compute CV-RMSE for a GAMM with given formula
  rmses <- vapply(folds, function(idx) {
    tr <- dat[-idx, ]
    te <- dat[idx, ]
    m  <- fit_gamm(response, rhs, tr)
    p  <- predict(m$gam, newdata = te)
    sqrt(mean((te[[response]] - p)^2, na.rm = TRUE))
  }, numeric(1))
  mean(rmses)
}

run_grid <- function(response, dat, folds, grid) {
  # Run a module grid (different formula) for one response var
  # Returns a list with model, AIC, CV-RMSE per formula
  lapply(names(grid), function(nm) {
    rhs <- grid[[nm]]
    m   <- fit_gamm(response, rhs, dat)
    list(
      name = nm,
      model = m,
      AIC = AIC(m$lme),
      cv_rmse = cv_rmse(response, rhs, dat, folds)
    )
  }) |> setNames(names(grid))
}

# Extract residuals aligned by date
extract_resids <- function(mod, dat, resp) {
  pr <- predict(mod$gam, newdata = dat)
  tibble(Date = dat$Date, resid = dat[[resp]] - pr)
}

var_parts <- function(model, data, response) {
  # Variance partitioning: total variance, abiotic (model) variance, residual variance
  preds <- predict(model$gam, newdata = data)
  residuals <- data[[response]] - preds
  tibble(
    total_var        = var(data[[response]], na.rm = TRUE),
    abiotic_var      = var(preds, na.rm = TRUE),
    resid_var        = var(residuals, na.rm = TRUE),
    abiotic_fraction = abiotic_var / total_var,
    resid_fraction   = resid_var   / total_var
  )
}

# Plot
plot_key_panels <- function(df_model, response, pred_col, resid_col,
                            station_name, start_time, end_time, base_color,
                            sal_col = "Salinity",
                            wind_col = "Wind.Speed",
                            ylim_raw = NULL, ylim_resid = NULL) {
  light_colors <- c("Night"="#d9d9d9","Astronomical twilight"="#ffb347",
                    "Nautical twilight"="#ffc870","Civil twilight"="#ffe0a3","Day"="#ffffb3")

  # ----------------- PANEL A: Observed vs Predicted (original units) -----------------
  p_top <- ggplot(df_model, aes(Date)) +
    geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date),
                  ymin = -Inf, ymax = Inf, fill = day_moment),
              alpha = 0.3, color = NA, na.rm = TRUE) +
    geom_line(aes(y = .data[[response]], color = "Observed"), linewidth = 1.1) +
    geom_line(aes(y = .data[[pred_col]], color = "Predicted\n(abiotic model)"), linewidth = 0.95) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Observed" = base_color,
                                  "Predicted\n(abiotic model)" = "grey40")) +
    scale_fill_manual(values = light_colors) +
    coord_cartesian(xlim = c(safe_time(start_time), safe_time(end_time)),
                    ylim = ylim_raw) +
    labs(y = paste0(response), x = NULL, color = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "right",
      axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    guides(fill = "none")

  # ---------------- PANEL B: Salinity ----------------
  p_sal <- ggplot(df_model, aes(Date, .data[[sal_col]])) +
    geom_rect(aes(xmin = Date, xmax = lead(Date), ymin = -Inf, ymax = Inf,
                  fill = day_moment), alpha = 0.3, color = NA) +
    geom_line(color = "#fdae61", linewidth = 0.9) +
    scale_fill_manual(values = light_colors) +
    coord_cartesian(xlim = c(safe_time(start_time), safe_time(end_time)),
                    ylim = c(29.5, 33)) +
    labs(y = "Salinity\n(psu)", x = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # ---------------- PANEL C: Wind speed ----------------
  p_wind <- ggplot(df_model, aes(Date, .data[[wind_col]])) +
    geom_rect(aes(xmin = Date, xmax = lead(Date), ymin = -Inf, ymax = Inf,
                  fill = day_moment), alpha = 0.3, color = NA) +
    geom_line(color = "#74add1", linewidth = 0.9) +
    scale_fill_manual(values = light_colors) +
    coord_cartesian(xlim = c(safe_time(start_time), safe_time(end_time)),
                    ylim = c(0,14)) +
    labs(y = "Wind speed\n(m/s)", x = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # ----------------- PANEL C: Residuals (µmol/kg) -----------------
  df_bot <- df_model %>%
    dplyr::select(Date, day_moment, resid = .data[[resid_col]]) %>%
    dplyr::arrange(Date)

  p_bot <- ggplot(df_bot, aes(Date, resid)) +
    geom_rect(aes(xmin = Date, xmax = lead(Date), ymin = -Inf, ymax = Inf, fill = day_moment),
              alpha = 0.25, color = NA, na.rm = TRUE) +
    geom_line(color = base_color, linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = light_colors) +
    coord_cartesian(xlim = c(safe_time(start_time), safe_time(end_time)),
                    ylim = ylim_resid) +
    labs(y = paste0(response, " Residuals\n(µmol/kg)"), x = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  cowplot::plot_grid(
    cowplot::ggdraw() + draw_label(
      paste0("Station ", station_name, " - ", response, " model"),
      fontface = "bold", size = 12
    ),
    cowplot::plot_grid(
      p_top, p_sal, p_wind, p_bot,
      ncol = 1,
      align = "v",
      rel_heights = c(1.0, 0.7, 0.7, 1.3),
      greedy = TRUE
    ),
    ncol = 1, rel_heights = c(0.12, 1)
  )
}

# ---- Read data ----------------------------------------------------------------
SiSt_202304 <- read.csv("data/raw/11SS20230401_DO.csv")
SiSt_202304$Date.Time <- as.POSIXct(strptime(SiSt_202304$Date.Time,"%Y-%m-%dT%H:%M:%S", tz="UTC"))

TA_sist2 <- -43.158 * SiSt_202304$P_sal..psu. + 3827.1

# Spot samples for validation
Spot_samples <- read.csv("data/raw/Spot_carb_jn2023.csv")
Spot_samples$Date.Time <- as.POSIXct(strptime(Spot_samples$Date.Time,"%Y-%m-%d %H:%M:%S", tz="UTC"))

# ICOS oxygen
oxygen_files <- list.files("data/raw/ICOS_Oxygen_processed-HT", pattern = "*.csv", full.names = TRUE)
oxygen_data  <- bind_rows(lapply(oxygen_files, read_csv, show_col_types = FALSE))
oxygen_data$Date <- safe_time(oxygen_data$Date)

# Sample metadata
sample_metadata <- read_csv("data/samples_env.csv", show_col_types = FALSE)

Sil <- mean(na.omit(sample_metadata$Si)) * 1e-6
PO4 <- mean(na.omit(sample_metadata$PO4)) * 1e-6

# Extract station sampling windows
data_51  <- SiSt_202304 %>% filter(Date.Time >= safe_time("2023-04-18 11:00:00"),
                                   Date.Time <= safe_time("2023-04-19 07:30:00"))
data_130 <- SiSt_202304 %>% filter(Date.Time >= safe_time("2023-04-20 08:00:00"),
                                   Date.Time <= safe_time("2023-04-21 08:30:00"))

# Build clean frames per station
mk_station_df <- function(d) {
  out <- tibble(
    Date.Time  = d$Date.Time,
    Salinity   = d$P_sal..psu.,
    Longitude  = d$Longitude,
    Latitude   = d$Latitude,
    Temp       = d$Temp..degC.,
    pCO2       = d$pCO2..uatm.,
    Wind.Speed = d$AWSWindSpeed,
    TA         = -43.158 * d$P_sal..psu. + 3827.1,

  ) %>% drop_na()
  out
}

data_51_clean  <- mk_station_df(data_51)
data_130_clean <- mk_station_df(data_130)

# Compute DIC (µmol/kg) from pCO2 and TA
calc_dic <- function(df) {
  D <- seacarb::carb(flag = 24, var1 = df$pCO2, var2 = df$TA * 1e-6,
                     S = df$Salinity, T = df$Temp, P = 0, Patm = 1.0,
                     Pt = PO4, Sit = Sil, pHscale = "T", kf = "pf", k1k2 = "l", ks = "d", b = "u74", warn = FALSE)
  D$DIC * 1e6
}

data_51_clean$DIC  <- calc_dic(data_51_clean)
data_130_clean$DIC <- calc_dic(data_130_clean)

# Match O2
data_51_clean  <- match_oxygen(data_51_clean,  oxygen_data)
data_130_clean <- match_oxygen(data_130_clean, oxygen_data)

# Aggregate sensor data into 5-min bins, to remove high frequency noise
aggregate_5min <- function(df, datetime_col = "Date", vars = NULL) {
  df <- df %>% dplyr::arrange(.data[[datetime_col]])
  df <- df %>% dplyr::mutate(Time_5min = lubridate::floor_date(.data[[datetime_col]], "5 minutes"))

  # Identify numeric columns to average if vars not specified
  if (is.null(vars)) {
    vars <- df %>% dplyr::select(where(is.numeric)) %>% colnames()
  }

  df_agg <- df %>%
    dplyr::group_by(Time_5min) %>%
    dplyr::summarise(across(all_of(vars), mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(Date = Time_5min)

  return(df_agg)
}

data_51_clean  <- aggregate_5min(data_51_clean)
data_130_clean <- aggregate_5min(data_130_clean)

# --------------- DIC comparison between stations ------------------------------
# Compute mean, SD, n per station
dic_stats <- tibble(
  Station = c("51", "130"),
  Mean_DIC = c(mean(data_51_clean$DIC, na.rm = TRUE),
               mean(data_130_clean$DIC, na.rm = TRUE)),
  SD_DIC = c(sd(data_51_clean$DIC, na.rm = TRUE),
             sd(data_130_clean$DIC, na.rm = TRUE)),
  N = c(sum(!is.na(data_51_clean$DIC)),
        sum(!is.na(data_130_clean$DIC)))
)
print(dic_stats)

# Perform t-test
dic_ttest <- t.test(data_51_clean$DIC, data_130_clean$DIC)
print(dic_ttest)

# Report difference
mean_diff <- mean(data_130_clean$DIC, na.rm = TRUE) - mean(data_51_clean$DIC, na.rm = TRUE)
cat(sprintf("\nMean difference (Station130 - Station51) = %.2f µmol/kg\n", mean_diff))

# --------------- Plot DIC and the Spot samples for validation -----------------
# Extract DIC
dic_cont <- bind_rows(
  data_51_clean  %>% dplyr::select(Date, DIC) %>% mutate(Station = "51"),
  data_130_clean %>% dplyr::select(Date, DIC) %>% mutate(Station = "130")
) %>%
  mutate(Date = as.POSIXct(Date, tz = "UTC")) %>%
  arrange(Station, Date)

# spot sample DIC
spots <- Spot_samples %>%
  mutate(
    Date = as.POSIXct(Date.Time, tz = "UTC"),
    Station = case_when(
      grepl("^JN51-", Sampling.Station) ~ "51",
      grepl("^130-",  Sampling.Station) ~ "130",
      TRUE ~ NA_character_
    ),
    DIC_spot = as.numeric(DIC),
    DIC_spot_sd = as.numeric(stdev.DIC)
  ) %>%
  filter(Station %in% c("51", "130")) %>%
  filter(!is.na(Station), !is.na(Date), is.finite(DIC_spot))

# Interpolate continuous DIC to exact spot times (per station)
interp_dic_at <- function(cont_df, t_query) {
  # rule = 1: linear interpolation; returns NA outside range
  approx(
    x = as.numeric(cont_df$Date),
    y = cont_df$DIC,
    xout = as.numeric(t_query),
    rule = 1
  )$y
}

# check nearest time offset to quantify matching success
nearest_dt_min <- function(cont_times, t_query) {
  # returns minimum absolute time difference (minutes) between query and cont_times
  vapply(t_query, function(tt) {
    min(abs(difftime(cont_times, tt, units = "mins")), na.rm = TRUE)
  }, numeric(1))
}

matched_spots <- spots %>%
  group_by(Station) %>%
  group_modify(~{
    st <- .y$Station[[1]]
    cont <- dic_cont %>% filter(Station == st)

    .x %>%
      mutate(
        DIC_calc = interp_dic_at(cont, Date),
        dt_min   = nearest_dt_min(cont$Date, Date)
      )
  }) %>%
  ungroup() %>%
  filter(is.finite(DIC_calc))

# Check time tolerance
max_dt_min <- 10  # minutes
matched_spots_f <- matched_spots %>%
  filter(dt_min <= max_dt_min) # Only 1 spot sample is not within 5 min

# Diagnostics
diag_tbl <- matched_spots_f %>%
  group_by(Station) %>%
  summarise(
    n = n(),
    r = cor(DIC_calc, DIC_spot, use = "complete.obs"),
    bias = mean(DIC_calc - DIC_spot, na.rm = TRUE),
    rmse = sqrt(mean((DIC_calc - DIC_spot)^2, na.rm = TRUE)),
    mae  = mean(abs(DIC_calc - DIC_spot), na.rm = TRUE),
    median_dt_min = median(dt_min, na.rm = TRUE),
    max_dt_min = max(dt_min, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Spot samples vs. continuous DIC diagnostics ---\n")
print(diag_tbl)

# helper to add station-wise annotation text in facets
ann <- diag_tbl %>%
  mutate(
    label = sprintf("n=%d\nr=%.2f\nbias=%.1f\nRMSE=%.1f\nmedian |Δt|=%.1f min",
                    n, r, bias, rmse, median_dt_min)
  )

# Scatter with 1:1 line + error bars
p_scatter <- ggplot(matched_spots_f, aes(x = DIC_calc, y = DIC_spot)) +
  geom_errorbar(aes(ymin = DIC_spot - DIC_spot_sd, ymax = DIC_spot + DIC_spot_sd),
                width = 0, alpha = 0.6, na.rm = TRUE) +
  geom_point(alpha = 0.85, size = 2.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ Station, scales = "free") +
  geom_text(
    data = ann,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.05,
    size = 3.2
  ) +
  labs(
    x = expression(paste("DIC from pCO"[2]*" sensor + TA(S) (", mu, "mol kg"^{-1}, ")")),
    y = expression(paste("Spot-sample DIC (", mu, "mol kg"^{-1}, ")"))
  ) +
  theme_minimal(base_size = 11)

# Difference vs mean
ba <- matched_spots_f %>%
  mutate(
    mean_DIC = (DIC_calc + DIC_spot) / 2,
    diff_DIC = DIC_calc - DIC_spot
  )

p_ba <- ggplot(ba, aes(x = mean_DIC, y = diff_DIC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.8, size = 2.3) +
  facet_wrap(~ Station, scales = "free_x") +
  # Add mean and limits of agreement per station
  stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.9) +
  labs(
    x = expression(paste("Mean(DIC"[calc]*", DIC"[spot]*") (", mu, "mol kg"^{-1}, ")")),
    y = expression(paste("DIC"[calc]*" - DIC"[spot]*" (", mu, "mol kg"^{-1}, ")"))
  ) +
  theme_minimal(base_size = 11)

fig_validation <- (p_scatter / p_ba) +
  plot_annotation(tag_levels = "A")

fig_validation

ggsave("figures/environmental/DIC_validation_spot_vs_continuous.png", fig_validation,
       width = 18, height = 16, units = "cm", dpi = 800)

# ----- O2 comparison between stations -----------------------------------------
# Compute mean, SD, n per station
o2_stats <- tibble(
  Station = c("51", "130"),
  Mean_O2 = c(mean(data_51_clean$O2uM, na.rm = TRUE),
                   mean(data_130_clean$O2uM, na.rm = TRUE)),
  SD_O2 = c(sd(data_51_clean$O2uM, na.rm = TRUE),
                 sd(data_130_clean$O2uM, na.rm = TRUE)),
  N = c(sum(!is.na(data_51_clean$O2uM)),
        sum(!is.na(data_130_clean$O2uM)))
)
print(o2_stats)

# t-test
o2_ttest <- t.test(data_51_clean$O2uM, data_130_clean$O2uM)
print(o2_ttest)

# Report difference
mean_diff_o2 <- mean(data_130_clean$O2uM, na.rm = TRUE) - mean(data_51_clean$O2uM, na.rm = TRUE)
cat(sprintf("\nMean difference (Station130 - Station51) = %.2f µmol/kg\n", mean_diff_o2))

#---------- Calculate O2', the O2 saturation anomaly (µmol/kg) -----------------
## As in Garcia & Gordon, 1992 (https://doi.org/10.4319/lo.1992.37.6.1307)
## Basically, we compute O2 solubility at in situ T/S/P, then subtract that from observed O2
## to get the anomaly (positive = supersaturated, negative = undersaturated)
## The supersaturation is assumed to be driven by photosynthesis, the undersaturation by respiration
add_O2sat_anomaly <- function(df,
                              O2_col = "O2uM",
                              S_col  = "Salinity",
                              T_col  = "Temp",
                              p = 0,
                              verbose = TRUE,
                              digits = 6) {
  # Compute oxygen solubility (µmol/kg) using Garcia & Gordon (1992)
  # First, extract the relevant columns
  O2 <- df[[O2_col]]
  S  <- df[[S_col]]
  T  <- df[[T_col]]

  # Define coefficients used in Garcia & Gordon 1992
  A0 <-  5.80818;  A1 <-  3.20684;  A2 <-  4.11890;  A3 <-  4.93845
  A4 <-  1.01567;  A5 <-  1.41575
  B0 <- -0.00701211; B1 <- -0.00725958; B2 <- -0.00793334; B3 <- -0.00554491
  C0 <- -1.32412e-7

  # Temperature transform used in the formulation
  Ts <- log((298.15 - T) / (273.15 + T))

  # ln(solubility)
  lnC <- A0 + A1*Ts + A2*Ts^2 + A3*Ts^3 + A4*Ts^4 + A5*Ts^5 +
    S * (B0 + B1*Ts + B2*Ts^2 + B3*Ts^3) + C0 * S^2

  O2sat_umolkg <- exp(lnC)

  # Convert µmol/L -> µmol/kg using in-situ density
  rho_kg_m3 <- oce::swRho(salinity = S, temperature = T, pressure = p)
  rho_kg_L  <- rho_kg_m3 / 1000

  # Calculate O2' anomaly (µmol/kg)
  O2_obs_umolkg <- O2 / rho_kg_L
  O2prime_umolkg <- O2_obs_umolkg - O2sat_umolkg

  # Attach results
  df$O2sat_umolkg <- O2sat_umolkg
  df$O2prime      <- O2prime_umolkg
  df$O2_obs_umolkg  <- O2_obs_umolkg

  # Debug info for checking intermediate values and correlations
  dbg <- list(
    inputs = list(S = S, T = T, O2 = O2, p = p),
    intermediates = list(
      Ts = Ts,
      lnC = lnC,
      rho_kg_m3 = rho_kg_m3,
      rho_kg_L  = rho_kg_L,
      O2_obs_umolkg = O2_obs_umolkg
    ),
    outputs = list(
      O2sat_umolkg = O2sat_umolkg,
      O2prime_umolkg = O2prime_umolkg
    )
  )
  attr(df, "debug_O2prime") <- dbg

  if (verbose) {
    op <- options(digits = digits)
    on.exit(options(op), add = TRUE)

    cat("\n=== Some debug info to verify the calculation ===\n")
    cat("Columns:", O2_col, "(O2),", S_col, "(S),", T_col, "(T)\n")
    cat("N rows:", nrow(df), "\n")

    summ <- function(x) {
      c(
        n = sum(is.finite(x)),
        min = min(x, na.rm = TRUE),
        q25 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
        med = median(x, na.rm = TRUE),
        mean = mean(x, na.rm = TRUE),
        q75 = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
        max = max(x, na.rm = TRUE)
      )
    }

    cat("\n-- Inputs --\n")
    print(rbind(
      Salinity = summ(S),
      Temp     = summ(T),
      O2uM     = summ(O2)
    ))

    cat("\n-- Intermediate --\n")
    print(rbind(
      Ts        = summ(Ts),
      lnC       = summ(lnC),
      rho_kg_L  = summ(rho_kg_L),
      O2_umolkg = summ(O2_obs_umolkg)
    ))

    cat("\n-- Outputs --\n")
    print(rbind(
      O2sat_umolkg = summ(O2sat_umolkg),
      O2prime      = summ(O2prime_umolkg)
    ))

    # Does O2sat change with salinity at all?
    ok <- is.finite(S) & is.finite(O2sat_umolkg)
    if (sum(ok) > 5) {
      r <- suppressWarnings(cor(S[ok], O2sat_umolkg[ok]))
      cat("\nCor(S, O2sat_umolkg) =", r, "\n")
      cat("Range(O2sat_umolkg) =", range(O2sat_umolkg[ok]), "\n")
    }

    # Is O2' much different from O2uM?
    ok <- is.finite(O2prime_umolkg) & is.finite(O2_obs_umolkg)
    if (sum(ok) > 5) {
      r <- suppressWarnings(cor(O2prime_umolkg[ok], O2_obs_umolkg[ok]))
      cat("\nCor(O2prime, O2_obs_umolkg) =", r, "\n")
    }

    # “Delta if set S=35” quick check (recompute O2sat only, same T)
    if (sum(is.finite(S) & is.finite(T)) > 5) {
      S35 <- rep(35, length(S))
      lnC_35 <- A0 + A1*Ts + A2*Ts^2 + A3*Ts^3 + A4*Ts^4 + A5*Ts^5 +
        S35 * (B0 + B1*Ts + B2*Ts^2 + B3*Ts^3) + C0 * S35^2
      O2sat_35 <- exp(lnC_35)
      d_sat <- O2sat_35 - O2sat_umolkg
      cat("\nSummary of ΔO2sat if salinity would be a constant 35 PSU:\n")
      print(summ(d_sat))
    }

    cat("=== end debug ===\n\n")
  }

  df
}

data_130_clean <- add_O2sat_anomaly(data_130_clean, verbose = TRUE)
data_51_clean  <- add_O2sat_anomaly(data_51_clean, verbose = TRUE)

# Check the standard deviation of O2', O2sat, and O2obs to see how much variability is in each componentr
sd(attr(data_130_clean, "debug_O2prime")$intermediates$O2_obs_umolkg, na.rm=TRUE)
sd(data_130_clean$O2sat_umolkg, na.rm=TRUE)
sd(data_130_clean$O2prime, na.rm=TRUE)

sd(attr(data_51_clean, "debug_O2prime")$intermediates$O2_obs_umolkg, na.rm=TRUE)
sd(data_51_clean$O2sat_umolkg, na.rm=TRUE)
sd(data_51_clean$O2prime, na.rm=TRUE)

# Check correlations to see if O2' is more related to T/S or O2sat
with(data_130_clean, cor(O2prime, Temp,     use="complete.obs"))
with(data_130_clean, cor(O2prime, Salinity, use="complete.obs"))
with(data_130_clean, cor(O2prime, O2sat_umolkg, use="complete.obs"))

with(data_51_clean, cor(O2prime, Temp,     use="complete.obs"))
with(data_51_clean, cor(O2prime, Salinity, use="complete.obs"))
with(data_51_clean, cor(O2prime, O2sat_umolkg, use="complete.obs"))

# QC: Plot O2obs vs O2sat
df <- data_51_clean %>%
  arrange(Date) %>%
  mutate(
    # divergence between observed and O2sat
    div = O2_obs_umolkg - O2sat_umolkg,
    ymin = pmin(O2_obs_umolkg, O2sat_umolkg),
    ymax = pmax(O2_obs_umolkg, O2sat_umolkg)
  )

# Define shared y limits based on the range of both stations
y_min <- min(c(data_51_clean$O2_obs_umolkg, data_130_clean$O2_obs_umolkg, data_51_clean$O2sat_umolkg, data_130_clean$O2sat_umolkg), na.rm = TRUE)
y_max <- max(c(data_51_clean$O2_obs_umolkg, data_130_clean$O2_obs_umolkg, data_51_clean$O2sat_umolkg, data_130_clean$O2sat_umolkg), na.rm = TRUE)
ylim_shared <- c(floor(y_min / 10) * 10, ceiling(y_max / 10) * 10)

p1 <- ggplot(df, aes(x = Date)) +
  # shaded band between the two series
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.20) +
  # two time series
  geom_line(aes(y = O2_obs_umolkg, colour = "Observed O2"), linewidth = 0.9) +
  geom_line(aes(y = O2sat_umolkg,      colour = "Saturation O2"), linewidth = 0.9) +
  scale_colour_manual(values = c("Observed O2" = "black",
                                 "Saturation O2" = "darkred")) +
  labs(
    title = "Station 51",
    x = NULL,
    y = NULL,
    colour = NULL
  ) +
  ylim(ylim_shared) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# Repeat for Station 130
df <- data_130_clean %>%
  arrange(Date) %>%
  mutate(
    # divergence between observed and O2sat
    div = O2_obs_umolkg - O2sat_umolkg,
    ymin = pmin(O2_obs_umolkg, O2sat_umolkg),
    ymax = pmax(O2_obs_umolkg, O2sat_umolkg)
  )

p2 <- ggplot(df, aes(x = Date)) +
  # shaded band between the two series
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.20) +
  # two time series
  geom_line(aes(y = O2_obs_umolkg, colour = "Observed O2"), linewidth = 0.9) +
  geom_line(aes(y = O2sat_umolkg,      colour = "Saturation O2"), linewidth = 0.9) +
  scale_colour_manual(values = c("Observed O2" = "black",
                                 "Saturation O2" = "darkred")) +
  labs(
    title = "Station 130",
    x = NULL,
    y = expression(paste("O"[2], " (", mu, "mol kg"^{-1}, ")")),
    colour = NULL
  ) +
  ylim(ylim_shared) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Save composite plot with station names
plot <- p1 + p2
plot

ggsave("figures/environmental/O2obs_O2sat_comparison.pdf", plot,
       width = 20, height = 12, units = "cm", dpi = 800)
ggsave("figures/environmental/O2obs_O2sat_comparison.png", plot,
       width = 20, height = 12, units = "cm", dpi = 800)

# ---- O2' comparison between stations -----------------------------------------
cat("\n--- Mean O2' comparison between Station 51 and 130 ---\n")
# Compute mean, SD, n per station
o2_stats <- tibble(
  Station = c("51", "130"),
  Mean_O2prime = c(mean(data_51_clean$O2prime, na.rm = TRUE),
                   mean(data_130_clean$O2prime, na.rm = TRUE)),
  SD_O2prime = c(sd(data_51_clean$O2prime, na.rm = TRUE),
                 sd(data_130_clean$O2prime, na.rm = TRUE)),
  N = c(sum(!is.na(data_51_clean$O2prime)),
        sum(!is.na(data_130_clean$O2prime)))
)
print(o2_stats)

# t-test
o2_ttest <- t.test(data_51_clean$O2prime, data_130_clean$O2prime)
print(o2_ttest)

# Report difference
mean_diff_o2 <- mean(data_130_clean$O2prime, na.rm = TRUE) - mean(data_51_clean$O2prime, na.rm = TRUE)
cat(sprintf("\nMean difference (Station130 - Station51) = %.2f µmol/kg\n", mean_diff_o2))

# ---- Sun shading fields -------------------------------------------------------
lat_130 <- median(data_130_clean$Latitude, na.rm = TRUE); lon_130 <- median(data_130_clean$Longitude, na.rm = TRUE)
lat_51  <- median(data_51_clean$Latitude,  na.rm = TRUE); lon_51  <- median(data_51_clean$Longitude,  na.rm = TRUE)

sun_130 <- get_station_sun_times(lat_130, lon_130, unique(as.Date(data_130_clean$Date)))
sun_51  <- get_station_sun_times(lat_51,  lon_51,  unique(as.Date(data_51_clean$Date)))

data_51_clean  <- add_day_moment(data_51_clean,  sun_51)
data_130_clean <- add_day_moment(data_130_clean, sun_130)

# ---- Build model data (scaled + AR.start) ------------------------------------
data_51_model      <- make_model_df(data_51_clean)
data_130_model     <- make_model_df(data_130_clean)

# Add day_moment for plotting
data_51_model$day_moment      <- data_51_clean$day_moment
data_130_model$day_moment     <- data_130_clean$day_moment

# ---- Model grids -----------------------------------------------------
model_grid <- list(
  minimal = ~ s(Salinity_scaled, k = 20, bs = "cs"),
  full    = ~ s(Salinity_scaled, k = 20, bs = "cs") +
              s(Wind_scaled, k = 20, bs = "cs")
)

# ---- Time-block CV folds per response ----------------------------------
folds_51   <- make_time_blocks(nrow(data_51_model), k = 5)
folds_130  <- make_time_blocks(nrow(data_130_model), k = 5)

# ---- Run grids ---------------------------------------------------------------
grid_51_O2   <- run_grid("O2prime", data_51_model, folds_51, model_grid)
grid_51_DIC  <- run_grid("DIC", data_51_model, folds_51, model_grid)
grid_130_O2  <- run_grid("O2prime", data_130_model, folds_130, model_grid)
grid_130_DIC <- run_grid("DIC", data_130_model, folds_130, model_grid)

tab_from_grid <- function(G) {
  tibble(
    model   = names(G),
    AIC     = sapply(G, \(x) x$AIC),
    CV_RMSE = sapply(G, \(x) x$cv_rmse)
  ) %>% arrange(match(model, names(model_grid)))
}

cat("\n--- Station 51 O2' ---\n");  print(tab_from_grid(grid_51_O2))
cat("\n--- Station 51 DIC ---\n"); print(tab_from_grid(grid_51_DIC))
cat("\n--- Station 130 O2' ---\n"); print(tab_from_grid(grid_130_O2))
cat("\n--- Station 130 DIC ---\n"); print(tab_from_grid(grid_130_DIC))

# Choose the best model per response
best_51_O2   <- names(grid_51_O2)[which.min(sapply(grid_51_O2, \(x) x$cv_rmse))]
best_51_DIC  <- names(grid_51_DIC)[which.min(sapply(grid_51_DIC, \(x) x$cv_rmse))]
best_130_O2  <- names(grid_130_O2)[which.min(sapply(grid_130_O2, \(x) x$cv_rmse))]
best_130_DIC <- names(grid_130_DIC)[which.min(sapply(grid_130_DIC, \(x) x$cv_rmse))]

# Fit models
final_51_O2   <- grid_51_O2 [[best_51_O2]]$model
final_51_DIC  <- grid_51_DIC[[best_51_DIC]]$model
final_130_O2  <- grid_130_O2[[best_130_O2]]$model
final_130_DIC <- grid_130_DIC[[best_130_DIC]]$model

# ---- Predictions, residuals, derivatives -------------------------------------
augment_block <- function(dat, model, resp, label) {
  # Predict abiotic component
  pred <- as.numeric(predict(model$gam, newdata = dat))
  # Compute residuals
  res  <- dat[[resp]] - pred
  # Return additional columns in original dataframe
  dat %>%
    mutate(
      !!paste0(label, "_pred")  := pred,
      !!paste0(label, "_resid") := res
    )
}

data_51_model <- data_51_model %>%
  augment_block(final_51_O2,  "O2prime", "O2") %>%
  augment_block(final_51_DIC, "DIC",     "DIC")

data_130_model <- data_130_model %>%
  augment_block(final_130_O2,  "O2prime", "O2") %>%
  augment_block(final_130_DIC, "DIC",     "DIC")

# ----------------- Variance partitioning summary ------------------------------
vp_51_O2   <- var_parts(final_51_O2,   data_51_model, "O2prime")
vp_51_DIC  <- var_parts(final_51_DIC,  data_51_model, "DIC")
vp_130_O2  <- var_parts(final_130_O2,  data_130_model, "O2prime")
vp_130_DIC <- var_parts(final_130_DIC, data_130_model, "DIC")

cat("\n--- Variance partition (abiotic vs residual) ---\n")
cat("Station 51 O2:\n");   print(vp_51_O2)
cat("Station 51 DIC:\n");  print(vp_51_DIC)
cat("Station 130 O2:\n");  print(vp_130_O2)
cat("Station 130 DIC:\n"); print(vp_130_DIC)

# ---- Interval stats for O2' super/undersaturation ----------------------------
# Signed trapezoid AUC over an interval (units: µmol/kg * hours)
trapz_auc <- function(time_posix, y, units = "hours") {
  ok <- !is.na(time_posix) & is.finite(y)
  t <- time_posix[ok]; y <- y[ok]
  o <- order(t); t <- t[o]; y <- y[o]
  tnum <- as.numeric(difftime(t, t[1], units = units))
  sum(diff(tnum) * (y[-1] + y[-length(y)]) / 2)
}

interval_stats_O2prime <- function(df, station_lab,
                                   time_col = "Date",
                                   value_col = "O2prime",
                                   units = "hours",
                                   min_duration_min = 0) {

  d <- df %>%
    dplyr::select(time = all_of(time_col), value = all_of(value_col)) %>%
    dplyr::filter(!is.na(time), is.finite(value)) %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(state = dplyr::case_when(
      value > 0 ~ "super",
      value < 0 ~ "under",
      TRUE      ~ "zero"
    ))

  if (nrow(d) < 2) return(tibble::tibble())

  # contiguous runs of same state
  d <- d %>%
    dplyr::mutate(run_id = cumsum(state != dplyr::lag(state, default = first(state))))

  out <- d %>%
    dplyr::group_by(run_id, state) %>%
    dplyr::summarise(
      Station = station_lab,
      start = first(time),
      end   = last(time),
      duration_h = as.numeric(difftime(end, start, units = units)),
      n = dplyr::n(),

      mean_O2prime = mean(value, na.rm = TRUE),
      min_O2prime  = min(value, na.rm = TRUE),
      max_O2prime  = max(value, na.rm = TRUE),

      t_min = time[which.min(value)][1],
      t_max = time[which.max(value)][1],

      auc_signed = trapz_auc(time, value, units = units),

      .groups = "drop"
    ) %>%
    dplyr::mutate(duration_min = duration_h * 60) %>%
    dplyr::filter(duration_min >= min_duration_min)

  out
}

# Summary per station/state
summarize_intervals <- function(interval_tbl) {
  interval_tbl %>%
    dplyr::group_by(Station, state) %>%
    dplyr::summarise(
      n_intervals = dplyr::n(),
      total_duration_h = sum(duration_h, na.rm = TRUE),
      mean_interval_h  = mean(duration_h, na.rm = TRUE),
      median_interval_h= median(duration_h, na.rm = TRUE),
      # total signed AUC across all intervals of that state
      total_auc_signed = sum(auc_signed, na.rm = TRUE),
      .groups = "drop"
    )
}

int_51  <- interval_stats_O2prime(data_51_clean,  station_lab = "51",  min_duration_min = 30)
int_130 <- interval_stats_O2prime(data_130_clean, station_lab = "130", min_duration_min = 30)

# Combine and view
intervals_all <- dplyr::bind_rows(int_51, int_130)
print(intervals_all)

# Summary by station and saturation state
intervals_summary <- summarize_intervals(intervals_all)
print(intervals_summary)

# ---- Plotting --------------------------------------
# Time windows
start_51  <- "2023-04-18 11:00:00"; end_51   <- "2023-04-19 08:00:00"
start_130 <- "2023-04-20 08:00:00"; end_130  <- "2023-04-21 08:00:00"

# O2 panels
p51_O2prime  <- plot_key_panels(
  df_model   = data_51_model,
  response   = "O2prime",
  pred_col   = "O2_pred",
  resid_col  = "O2_resid",
  station_name = "51",
  start_time   = start_51,
  end_time     = end_51,
  base_color   = "darkgreen",
  ylim_raw   = c(-25, 30),
  ylim_resid = c(-12, 12)
)

p130_O2prime <- plot_key_panels(
  df_model   = data_130_model,
  response   = "O2prime",
  pred_col   = "O2_pred",
  resid_col  = "O2_resid",
  station_name = "130",
  start_time   = start_130,
  end_time     = end_130,
  base_color   = "darkgreen",
  ylim_raw   = c(-25, 30),
  ylim_resid = c(-12, 12)
)

# DIC panels
p51_DIC <- plot_key_panels(
  df_model   = data_51_model,
  response   = "DIC",
  pred_col   = "DIC_pred",
  resid_col  = "DIC_resid",
  station_name = "51",
  start_time   = start_51,
  end_time     = end_51,
  base_color   = "#0d5f9e",
  ylim_raw   = c(2100, 2400)
)

p130_DIC <- plot_key_panels(
  df_model   = data_130_model,
  response   = "DIC",
  pred_col   = "DIC_pred",
  resid_col  = "DIC_resid",
  station_name = "130",
  start_time   = start_130,
  end_time     = end_130,
  base_color   = "#0d5f9e",
  ylim_raw   = c(2100, 2400)
)

# Plot together
(p51_O2prime | p130_O2prime) / (p51_DIC | p130_DIC) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Save plots
ggsave("figures/environmental/Station51_O2.pdf",  p51_O2prime,
       width = 12, height = 14, units = "cm")
ggsave("figures/environmental/Station130_O2.pdf", p130_O2prime,
       width = 12, height = 14, units = "cm")
ggsave("figures/environmental/Station51_DIC.pdf", p51_DIC,
       width = 12, height = 14, units = "cm")
ggsave("figures/environmental/Station130_DIC.pdf", p130_DIC,
       width = 12, height = 14, units = "cm")

# Save composite figure
ggsave("figures/environmental/Stations51_130_O2_DIC.pdf",
       (p51_O2prime | p130_O2prime) / (p51_DIC | p130_DIC) +
         plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 18, height = 16, units = "cm")

# ---- 2-panel O2 / O2' + residuals figure -------------------------------------
# Light phases
light_colors <- c("Night" = "#d9d9d9",
                  "Astronomical twilight" = "#ffb347",
                  "Nautical twilight" = "#ffc870",
                  "Civil twilight" = "#ffe0a3",
                  "Day" = "#ffffb3")

# Build background rectangles from Date + day_moment
make_bg_df <- function(df) {
  df %>%
    dplyr::select(Date, day_moment) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Date) %>%
    dplyr::mutate(Date_start = Date,
                  Date_end   = dplyr::lead(Date)) %>%
    dplyr::filter(!is.na(Date_end))
}

# Update time windows as POSIXct
start_51  <- as.POSIXct("2023-04-18 11:00:00", tz = "UTC")
end_51    <- as.POSIXct("2023-04-19 08:00:00", tz = "UTC")
start_130 <- as.POSIXct("2023-04-20 08:00:00", tz = "UTC")
end_130   <- as.POSIXct("2023-04-21 08:00:00", tz = "UTC")

# Add station names to plot data
df51  <- data_51_model  %>%
  dplyr::mutate(Station = "51") %>%
  dplyr::filter(Date >= start_51, Date <= end_51)

df130 <- data_130_model %>%
  dplyr::mutate(Station = "130") %>%
  dplyr::filter(Date >= start_130, Date <= end_130)

# Combine into a long df
O2_long <- dplyr::bind_rows(df51, df130) %>%
  dplyr::mutate(
    Station = factor(Station, levels = c("51","130"),
                     labels = c("Station 51", "Station 130")),
    day_moment = factor(day_moment,
                        levels = c("Night","Astronomical twilight","Nautical twilight",
                                   "Civil twilight","Day"))
  )

# Add z-scored O2', O2, pCO2 and DIC for comparison
O2_long <- O2_long %>%
  dplyr::mutate(
    O2prime_z = as.numeric(scale(O2prime)),
    O2_z      = as.numeric(scale(O2_obs_umolkg)),
    DIC_z  = as.numeric(scale(DIC)),
    pCO2_z = as.numeric(scale(pCO2))
  )

# Background rectangles per station
bg_51  <- make_bg_df(df51)  %>% mutate(Station = "Station 51")
bg_130 <- make_bg_df(df130) %>% mutate(Station = "Station 130")
bg_all <- bind_rows(bg_51, bg_130) %>%
  mutate(Station = factor(Station, levels = c("Station 51","Station 130")))

# ---------- O2 over time  ----------
p_O2 <- ggplot(O2_long, aes(x = Date)) +
  geom_rect(
    data = bg_all,
    aes(xmin = Date_start, xmax = Date_end, ymin = -Inf, ymax = Inf, fill = day_moment),
    inherit.aes = FALSE, alpha = 0.25, color = NA
  ) +
  geom_line(aes(y = O2uM), color = "black", linewidth = 0.8) +
  geom_point(aes(y = O2uM), color = "black", size = 0.35, alpha = 0.6) +
  facet_grid(. ~ Station, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = light_colors, name = "Diel phase") +
  scale_x_datetime(
    date_breaks = "2 hours",
    date_labels = "%H:%M",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_minimal(base_size = 11) +
  labs(x = NULL, y = expression(paste("O"[2], " (", mu, "mol L"^{-1}, ")"))) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

# ---------- Z-scored O2' and O2 to compare dynamics ----------
p_O2prime_z_O2_z <- ggplot(O2_long, aes(x = Date)) +
  geom_rect(
    data = bg_all,
    aes(xmin = Date_start, xmax = Date_end, ymin = -Inf, ymax = Inf, fill = day_moment),
    inherit.aes = FALSE, alpha = 0.25, color = NA
  ) +
  geom_line(aes(y = O2_z, color = "O2 (z)"), linewidth = 0.9, alpha = 0.8) +
  geom_line(aes(y = O2prime_z, color = "O2' (z)"), linewidth = 0.9, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40", linewidth = 0.6) +
  facet_grid(. ~ Station, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = light_colors, name = "Diel phase") +
  scale_color_manual(values = c("O2 (z)" = "black", "O2\' (z)" = "darkgreen"), name = NULL) +
  scale_x_datetime(
    date_breaks = "2 hours",
    date_labels = "%H:%M",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_minimal(base_size = 11) +
  labs(x = "Time of day", y = "Standardized units (z-score)") +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

# ---------- O2' + residual bars + abiotic smooth ----------
p_O2prime <- ggplot(O2_long, aes(x = Date)) +
  geom_rect(
    data = bg_all,
    aes(xmin = Date_start, xmax = Date_end, ymin = -Inf, ymax = Inf, fill = day_moment),
    inherit.aes = FALSE, alpha = 0.25, color = NA
  ) +
  geom_bar(aes(y = O2_resid), stat = "identity", fill = "#CEDCA0", alpha = 0.9) +
  geom_line(aes(y = O2prime), color = "darkgreen", linewidth = 0.9, alpha = 0.8) +
  geom_line(aes(y = O2_pred), color = "grey40", linewidth = 0.7, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40", linewidth = 0.6) +
  facet_grid(. ~ Station, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = light_colors, name = "Diel phase") +
  scale_x_datetime(
    date_breaks = "2 hours",
    date_labels = "%H:%M",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_minimal(base_size = 11) +
  labs(x = "Time of day", y = expression(paste("O"[2], "'" , " (", mu, "mol kg"^{-1}, ")"))) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

fig_O2_suite <- p_O2 / p_O2prime_z_O2_z + plot_annotation(tag_levels = "A")
fig_O2_suite

# Save
ggsave("figures/environmental/O2_O2prime_zscore.pdf",
       fig_O2_suite, width = 18, height = 16, units = "cm")
ggsave("figures/environmental/O2_O2prime_zscore.png",
       fig_O2_suite, width = 18, height = 16, units = "cm", dpi = 800)

# ---------- pCO2 over time  ----------
p_pCO2 <- ggplot(O2_long, aes(x = Date)) +
  geom_rect(
    data = bg_all,
    aes(xmin = Date_start, xmax = Date_end, ymin = -Inf, ymax = Inf, fill = day_moment),
    inherit.aes = FALSE, alpha = 0.25, color = NA
  ) +
  geom_line(aes(y = pCO2), color = "black", linewidth = 0.8) +
  facet_grid(. ~ Station, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = light_colors, name = "Diel phase") +
  scale_x_datetime(
    date_breaks = "2 hours",
    date_labels = "%H:%M",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_minimal(base_size = 11) +
  labs(x = NULL, y = expression(paste("pCO"[2], " (", mu, "atm)"))) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

p_pCO2

p_pCO2_DIC <- ggplot(O2_long, aes(x = Date)) +
  geom_rect(
    data = bg_all,
    aes(xmin = Date_start, xmax = Date_end, ymin = -Inf, ymax = Inf, fill = day_moment),
    inherit.aes = FALSE, alpha = 0.25, color = NA
  ) +
  geom_line(aes(y = pCO2_z, color = "pCO2 (z)"), linewidth = 0.8) +
  geom_line(aes(y = DIC_z, color = "DIC (z)"), linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40", linewidth = 0.6) +
  facet_grid(. ~ Station, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = light_colors, name = "Diel phase") +
  scale_color_manual(values = c("pCO2 (z)" = "black", "DIC (z)" = "#0d5f9e"), name = NULL) +
  scale_x_datetime(
    date_breaks = "2 hours",
    date_labels = "%H:%M",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_minimal(base_size = 11) +
  labs(x = "Time of day", y = "Standardized units (z-score)") +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

p_pCO2_DIC

plot <- p_pCO2 / p_pCO2_DIC + plot_annotation(tag_levels = "A")
plot

# Save
ggsave("figures/environmental/pCO2_DIC_comparison.pdf", plot,
       width = 18, height = 16, units = "cm")
ggsave("figures/environmental/pCO2_DIC_comparison.png", plot,
       width = 18, height = 16, units = "cm")

# ---- Statistical summary table ------------------------------------------------
# Sensitivity table function: How does the relationship between DIC and O2 residuals
# vary depending on model choice? Also includes CV-RMSE and AIC for each model
sens_table <- function(G_O2, G_DIC, dat_O2, dat_DIC, station_lab) {
  rows <- lapply(names(model_grid), function(nm) {
    mO2  <- G_O2[[nm]]$model
    mDIC <- G_DIC[[nm]]$model
    rO2  <- extract_resids(mO2,  dat_O2,  "O2prime") %>% dplyr::rename(O2_resid = resid)
    rDIC <- extract_resids(mDIC, dat_DIC, "DIC")  %>% dplyr::rename(DIC_resid = resid)
    R    <- inner_join(rO2, rDIC, by = "Date") %>% arrange(Date)
    lm1  <- lm(DIC_resid ~ O2_resid, data = R)
    tibble(
      station = station_lab,
      model   = nm,
      slope   = unname(coef(lm1)["O2_resid"]),
      slope_se= summary(lm1)$coefficients["O2_resid","Std. Error"],
      slope_p = summary(lm1)$coefficients["O2_resid","Pr(>|t|)"],
      resid_r = cor(R$O2_resid, R$DIC_resid, use = "complete.obs"),
      O2_CV   = G_O2[[nm]]$cv_rmse,
      DIC_CV  = G_DIC[[nm]]$cv_rmse,
      O2_AIC  = G_O2[[nm]]$AIC,
      DIC_AIC = G_DIC[[nm]]$AIC,
      n       = nrow(R)
    )
  })
  bind_rows(rows) %>% arrange(match(model, names(model_grid)))
}

cat("\nSUMMARY\n")
sens_51  <- sens_table(grid_51_O2,  grid_51_DIC,  data_51_model,  data_51_model,  "51")
print(sens_51)
sens_130 <- sens_table(grid_130_O2, grid_130_DIC, data_130_model, data_130_model, "130")
print(sens_130)

phi_51_O2  <- coef(grid_51_O2[[sym_51]]$model$lme$modelStruct$corStruct, unconstrained = FALSE)
phi_130_O2 <- coef(grid_130_O2[[sym_130]]$model$lme$modelStruct$corStruct, unconstrained = FALSE)
cat("Estimated Φ (AR1) values: Station 51 =", phi_51_O2, ", Station 130 =", phi_130_O2, "\n")

# Print out the Φ AR1 for all models to check if they are similar across models
cat("\nAR1 Φ values for all models:\n")
print(data.frame(
  Station = c(rep("51", length(grid_51_O2)), rep("130", length(grid_130_O2))),
  Model   = c(names(grid_51_O2), names(grid_130_O2)),
  Phi_O2  = c(sapply(grid_51_O2, function(x) coef(x$model$lme$modelStruct$corStruct, unconstrained = FALSE)),
             sapply(grid_130_O2, function(x) coef(x$model$lme$modelStruct$corStruct, unconstrained = FALSE)))
))

summ_line <- function(st, var, rhs, AIC, CV) sprintf(
  "Station %s %s: RHS = %s | AIC = %.1f | CV-RMSE = %.3f",
  st, var, rhs, AIC, CV
)
cat(summ_line("51","O2",  sym_51,
              grid_51_O2[[sym_51]]$AIC,  grid_51_O2[[sym_51]]$cv_rmse), "\n")
cat(summ_line("51","DIC", sym_51,
              grid_51_DIC[[sym_51]]$AIC, grid_51_DIC[[sym_51]]$cv_rmse), "\n")
cat(summ_line("130","O2",  sym_130,
              grid_130_O2[[sym_130]]$AIC,  grid_130_O2[[sym_130]]$cv_rmse), "\n")
cat(summ_line("130","DIC", sym_130,
              grid_130_DIC[[sym_130]]$AIC, grid_130_DIC[[sym_130]]$cv_rmse), "\n")

fmt_vp <- function(vp) sprintf("abiotic_fraction=%.3f, resid_fraction=%.3f",
                               vp$abiotic_fraction, vp$resid_fraction)
cat("Station 51 O2 var parts:  ", fmt_vp(vp_51_O2),  "\n")
cat("Station 51 DIC var parts: ", fmt_vp(vp_51_DIC), "\n")
cat("Station 130 O2 var parts: ", fmt_vp(vp_130_O2), "\n")
cat("Station 130 DIC var parts:", fmt_vp(vp_130_DIC), "\n")

cat('Station 51 GAM check')
gam.check(final_51_O2$gam)
gam.check(final_51_DIC$gam)
cat('Station 130 GAM check')
gam.check(final_130_O2$gam)
gam.check(final_130_DIC$gam)

cat(sprintf("Station 51 residual stoichiometry: slope = %.3f (SE=%.3f, p=%s), r=%.3f, n=%d\n",
            stoich_51$slope, stoich_51$slope_se, signif(stoich_51$p_value,3), stoich_51$cor_resid, stoich_51$n))
cat(sprintf("Station 130 residual stoichiometry: slope = %.3f (SE=%.3f, p=%s), r=%.3f, n=%d\n",
            stoich_130$slope, stoich_130$slope_se, signif(stoich_130$p_value,3), stoich_130$cor_resid, stoich_130$n))

# Combine both stations and export the O2' anomaly data and residuals for downstream analysis
data_export <- bind_rows(
  data_51_model %>%
    dplyr::select(Date, Salinity, Temp, pCO2, Wind.Speed, TA, DIC, O2uM, O2sat_umolkg, O2prime,
           O2_pred, O2_resid) %>%
    mutate(Station = "51"),
  data_130_model %>%
    dplyr::select(Date, Salinity, Temp, pCO2, Wind.Speed, TA, DIC, O2uM, O2sat_umolkg, O2prime,
           O2_pred, O2_resid) %>%
    mutate(Station = "130")
) %>%
  relocate(Station, .before = Date)
write_csv(data_export, "data/analysis/O2prime_resids.csv")

# Plot correlation between O2' and DIC, and O2' residuals and DIC residuals
data_export$Station <- factor(data_export$Station, levels = c("51", "130"))

correlation_plot <- ggplot(data_export, aes(x = O2prime, y = DIC)) +
  geom_point(alpha = 0.6, color = "grey") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           method = "pearson", label.x = min(data_export$O2prime), label.y = max(data_export$DIC), size = 4) +
  facet_wrap(~ Station, scales = "free", labeller = labeller(Station = c("51" = "Station 51", "130" = "Station 130"))) +
  labs(x = expression(paste(O[2]*"'"~" ("*mu*"mol kg"^{-1}*")")),
       y = expression(paste("DIC ("*mu*"mol kg"^{-1}*")"))) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 11, face = "bold"))

correlation_plot

# Save the figure as an pdf and png file
output_dir <- "figures/environmental"
output_path <- file.path(output_dir, "O2prime_DIC_correlation.pdf")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "pdf")
output_path <- file.path(output_dir, "O2prime_DIC_correlation.png")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "png")

# Fit the linear model: DIC as a function of O2
model <- lm(DIC ~ O2prime, data = data_130_model)

# Extract coefficients
slope <- coef(model)["O2prime"]
intercept <- coef(model)["(Intercept)"]

# Extract R-squared and p-value
summary_model <- summary(model)
r_value <- sqrt(summary_model$r.squared)  # R value (square root of R²)
p_value <- summary_model$coefficients[2, 4]  # P-value for O2' coefficient

# Print results
cat(sprintf("Equation: DIC = %.3f * O2 + %.3f\n", slope, intercept))
cat(sprintf("R-value: %.3f\n", r_value))
cat(sprintf("P-value: %.3e\n", p_value))
