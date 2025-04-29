# load packages -----------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(xts)
library(zoo)
library(data.table)
library(plyr)
library(dplyr)
library(grid)
library(broom)
library(scales)
library(stats)
library(DBI)
library(RSQLite)
library(TTR)
library(tibble)
library(lubridate)
library(seacarb)
library(readxl)
library(gridExtra)
library(cowplot)
library(svglite)
library(lme4)
library(lmerTest)
library(gamm4)
library(readr)
library(ggpubr)
library(suncalc)
library(caret)
library(mgcv)
library(gratia)
library(purrr)

Sys.setenv(TZ='UTC')


# read data ---------------------------------------------------------------

SiSt_202304 <- read.csv("data/raw/11SS20230401_DO.csv",sep = ",", dec=".", header=TRUE )

SiSt_202304$Date.Time <- as.POSIXct(strptime(SiSt_202304$Date.Time, format="%Y-%m-%dT%H:%M:%S", tz="UTC"), format="%Y-%m-%d %H:%M:%S",tz="UTC")

#salinity - alkalinity relationship y = -43.158x+3827.1 all data

TA_sist2 <- -43.158*SiSt_202304$P_sal..psu.+3827.1

Spot_samples <- read.csv("data/raw/Spot_carb_jn2023.csv",sep = ",", dec=".", header=TRUE )

Spot_samples$Date.Time <- as.POSIXct(strptime(Spot_samples$Date.Time, format="%Y-%m-%d %H:%M:%S", tz="UTC"), format="%Y-%m-%d %H:%M:%S",tz="UTC")

# Extract data for Station 51: 2023-04-18 13:00 to 2023-05-19 09:00 (local time, UTC+2)
data_51 <- SiSt_202304 %>%
  filter(Date.Time >= as.POSIXct("2023-04-18 11:00:00", tz = "UTC") &
           Date.Time <= as.POSIXct("2023-04-19 07:30:00", tz = "UTC"))

# Extract data for Station 130: 2023-04-20 10:00 to 2023-04-21 10:00 (local time, UTC+2)
data_130 <- SiSt_202304 %>%
  filter(Date.Time >= as.POSIXct("2023-04-20 08:00:00", tz = "UTC") &
           Date.Time <= as.POSIXct("2023-04-21 08:30:00", tz = "UTC"))

# Read ICOS oxygen data files
oxygen_files <- list.files("data/raw/ICOS_Oxygen_processed-HT", pattern = "*.csv", full.names = TRUE)
oxygen_data <- bind_rows(lapply(oxygen_files, read_csv))
oxygen_data$Date <- as.POSIXct(oxygen_data$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Visualize ---------------------------------------------------------------
plot_go1<-ggplot()   +
  geom_point(aes(x = SiSt_202304$Date.Time, y =TA_sist2, colour = "TA_UW2"), size = 0.2, alpha = 0.75) +
  geom_point(aes(x = Spot_samples$Date.Time, y= Spot_samples$Total.Alkalinity..Vindta., colour="Spot"), size = 2, alpha = 0.75)+
  ylab("TA umol/kg") + xlab("")+ylim(2300,2610)+xlim(as.POSIXct("2023-04-17 06:00:00"), as.POSIXct("2023-04-21 23:59:59"))

plot_go1

plot_go2<-ggplot()   +
  geom_point(aes(x = SiSt_202304$Date.Time, y = SiSt_202304$P_sal..psu.), size = 0.2, alpha = 0.75) +
  geom_point(aes(x = Spot_samples$Date.Time, y= Spot_samples$Saliniteit, colour="Spot"), size = 2, alpha = 0.75)+
  ylab("Salinity") +
  xlab("")+
  ylim(25,35)+
  xlim(as.POSIXct("2023-04-17 06:00:00"), as.POSIXct("2023-04-21 23:59:59"))

plot_go2

plot_go3 <-ggplot()   +
  geom_point(aes(x = SiSt_202304$Date.Time, y = SiSt_202304$Temp..degC.), size = 0.2, alpha = 0.75) +
  geom_point(aes(x = Spot_samples$Date.Time, y= Spot_samples$Temperature), size = 2, alpha = 0.75)+
  ylab("Temperature") + xlab("")+ylim(0,15)+xlim(as.POSIXct("2023-04-17 06:00:00"), as.POSIXct("2023-04-21 23:59:59"))

plot_go3

# reconstruct carbonate system parameters using seacarb-------------------------------------------------------------

#calculate DIC UW based on pCO2 and reconstructed TA_UW2

# average of PO4 and Si based on spot samples form JN 2023 cruise
sample_metadata <- read_csv("data/samples_env.csv")
Sil <- mean(na.omit(sample_metadata$Si))*10^(-6)
PO4 <- mean(na.omit(sample_metadata$PO4))*10^(-6)

#calculate DIC from pCO2 record and reconstructed TA from sal (from TA/Sal of all data!)
SiSt_clean <- data.frame(SiSt_202304$Date.Time, SiSt_202304$P_sal..psu., SiSt_202304$Temp..degC., SiSt_202304$pCO2..uatm., TA_sist2)

SiSt_clean <- na.omit(SiSt_clean)


DIC_UW <- carb(flag=24, var1=SiSt_clean$SiSt_202304.pCO2..uatm.,
               var2=SiSt_clean$TA_sist2*10^(-6),
               S=SiSt_clean$SiSt_202304.P_sal..psu.,
               T=SiSt_clean$SiSt_202304.Temp..degC.,
               P=0, Patm=1.0, Pt=PO4, Sit=Sil, pHscale="T",
               kf="pf", k1k2="l", ks="d", b="u74", warn=FALSE)

plot_go4 <-ggplot()   +
  geom_point(aes(x = SiSt_clean$SiSt_202304.Date.Time, y =DIC_UW$DIC*10^6, colour = "DIC_UW"), size = 0.2, alpha = 0.75) +
  ylab("DICumol/kg") +
  xlab("") +
  ylim(2000,2500) +
  xlim(as.POSIXct("2023-04-17 06:00:00"), as.POSIXct("2023-04-21 23:59:59"))

plot_go4

# Per station -------------------------------------------------------------
# Recalculate DIC for each station based on pCO2 and reconstructed TA
data_51_clean <- data.frame(Date.Time = data_51$Date.Time,
                            Salinity = data_51$P_sal..psu.,
                            Longitude = data_51$Longitude,
                            Latitude = data_51$Latitude,
                            Temp = data_51$Temp..degC.,
                            pCO2 = data_51$pCO2..uatm.,
                            Wind.Speed = data_51$AWSWindSpeed,
                            TA = -43.158 * data_51$P_sal..psu. + 3827.1)

data_130_clean <- data.frame(Date.Time = data_130$Date.Time,
                             Salinity = data_130$P_sal..psu.,
                             Longitude = data_130$Longitude,
                             Latitude = data_130$Latitude,
                             Temp = data_130$Temp..degC.,
                             pCO2 = data_130$pCO2..uatm.,
                             Wind.Speed = data_130$AWSWindSpeed,
                             TA = -43.158 * data_130$P_sal..psu. + 3827.1)

# Omit NAs
data_51_clean <- na.omit(data_51_clean)
data_130_clean <- na.omit(data_130_clean)

# Calculate DIC using seacarb
DIC_51 <- carb(flag = 24, var1 = data_51_clean$pCO2, var2 = data_51_clean$TA * 1e-6,
               S = data_51_clean$Salinity, T = data_51_clean$Temp, P = 0, Patm = 1.0,
               Pt = PO4, Sit = Sil, pHscale = "T", kf = "pf", k1k2 = "l", ks = "d", b = "u74",
               warn = FALSE)

DIC_130 <- carb(flag = 24, var1 = data_130_clean$pCO2, var2 = data_130_clean$TA * 1e-6,
                S = data_130_clean$Salinity, T = data_130_clean$Temp, P = 0, Patm = 1.0,
                Pt = PO4, Sit = Sil, pHscale = "T", kf = "pf", k1k2 = "l", ks = "d", b = "u74",
                warn = FALSE)

# Add DIC to data frames
data_51_clean$DIC <- DIC_51$DIC * 1e6
data_130_clean$DIC <- DIC_130$DIC * 1e6

# Function to match oxygen data based on Date-Time -----------------------------
match_oxygen <- function(data, oxygen_data) {
  # Convert data to data.table format for efficient joining
  oxygen_dt <- as.data.table(oxygen_data)
  data_dt <- as.data.table(data)

  # Ensure time columns are of the same class
  oxygen_dt[, Date := as.POSIXct(Date, tz = "UTC")]
  data_dt[, Date.Time := as.POSIXct(Date.Time, tz = "UTC")]
  # Rename column for joining
  setnames(data_dt, "Date.Time", "Date")

  # Perform a rolling join to get the closest match
  setkey(oxygen_dt, Date)
  matched_data <- oxygen_dt[data_dt, roll = "nearest", on = "Date", nomatch = 0]

  # Ensure original dataset structure is preserved
  data_dt$O2uM <- matched_data$O2uM

  # Convert back to tibble
  data_clean <- as_tibble(data_dt)
  # Return the cleaned data
  return(data_clean)
}

# Match oxygen data for Station 51
data_51_clean <- match_oxygen(data_51_clean, oxygen_data)

# Match oxygen data for Station 130
data_130_clean <- match_oxygen(data_130_clean, oxygen_data)

# Create plot for pCO2
pco2_plot <- ggplot() +
  geom_point(data = data_51_clean, aes(x = Date, y = pCO2, color = "Station 51"), size = 1, alpha = 0.75) +
  geom_point(data = data_130_clean, aes(x = Date, y = pCO2, color = "Station 130"), size = 1, alpha = 0.75) +
  labs(y = expression(paste("pCO"[2], " (", mu, "atm)")), x = "Date-Time", color = "Station") +
  scale_color_manual(values = c("Station 51" = "blue", "Station 130" = "red")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pco2_plot

# Create plot for DIC
dic_plot <- ggplot() +
  geom_point(data = data_51_clean, aes(x = Date, y = DIC, color = "Station 51"), size = 1, alpha = 0.75) +
  geom_point(data = data_130_clean, aes(x = Date, y = DIC, color = "Station 130"), size = 1, alpha = 0.75) +
  labs(y = expression(paste("DIC (", mu, "mol/kg)")), x = "Date-Time", color = "Station") +
  scale_color_manual(values = c("Station 51" = "blue", "Station 130" = "red")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dic_plot

# Create plot for Oxygen
o2_plot <- ggplot() +
  geom_point(data = data_51_clean, aes(x = Date, y = O2uM, color = "Station 51"), size = 1, alpha = 0.75) +
  geom_point(data = data_130_clean, aes(x = Date, y = O2uM, color = "Station 130"), size = 1, alpha = 0.75) +
  labs(y = expression(paste("Oxygen (", mu, "mol/L)")), x = "Date-Time", color = "Station") +
  scale_color_manual(values = c("Station 51" = "blue", "Station 130" = "red")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

o2_plot

# Function to create a dual y-axis plot with smoother
dual_axis_plot <- function(data, station_name, start_time, end_time, y_limits_pco2, y_limits_dic) {
  # Create the primary plot for pCO2
  p1 <- ggplot(data, aes(x = Date, y = pCO2)) +
    geom_point(size = 1, alpha = 0.75) +
    geom_smooth(method = "loess", color = "grey", se = TRUE) +
    labs(y = expression(paste("pCO"[2], " (", mu, "atm)")), x = NULL) +
    scale_x_datetime(limits = c(as.POSIXct(start_time, tz = "UTC"),
                                as.POSIXct(end_time, tz = "UTC")),
                     date_breaks = "6 hours", date_labels = "%d-%b %H:%M") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  # Create the secondary plot for DIC
  p2 <- ggplot(data, aes(x = Date, y = DIC)) +
    geom_point(size = 1, alpha = 0.75) +
    geom_smooth(method = "loess", color = "grey", se = TRUE) +
    labs(y = expression(paste("DIC (", mu, "mol/kg)")), x = "Date-Time") +
    scale_x_datetime(limits = c(as.POSIXct(start_time, tz = "UTC"),
                                as.POSIXct(end_time, tz = "UTC")),
                     date_breaks = "6 hours", date_labels = "%d-%b %H:%M") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Combine the two plots using gridExtra
  grid.arrange(p1, p2, ncol = 1, heights = c(3, 3),
               top = paste("Station", station_name))
}

# Set start and end times for the plots (time we were anchored)
start_time_51 <- "2023-04-18 11:00:00"
end_time_51 <- "2023-04-19 08:15:00"

start_time_130 <- "2023-04-20 08:00:00"
end_time_130 <- "2023-04-21 08:15:00"

# Plot for Station 51
dual_axis_plot(data_51_clean, station_name = "51", start_time = start_time_51, end_time = end_time_51)

# Plot for Station 130
dual_axis_plot(data_130_clean, station_name = "130", start_time = start_time_130, end_time = end_time_130)

# Function to create a triple-panel plot for pCO2, DIC, and Oxygen -------------
triple_axis_plot <- function(data, station_name, start_time, end_time) {
  # Plot for pCO2
  p1 <- ggplot(data, aes(x = Date, y = pCO2)) +
    geom_point(size = 1, alpha = 0.75) +
    geom_smooth(method = "loess", color = "grey", se = TRUE) +
    labs(y = expression(paste("pCO"[2], " (", mu, "atm)")), x = NULL) +
    scale_x_datetime(limits = c(as.POSIXct(start_time, tz = "UTC"),
                                as.POSIXct(end_time, tz = "UTC")),
                     date_breaks = "4 hours", date_labels = "%d-%b %H:%M") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  # Plot for DIC
  p2 <- ggplot(data, aes(x = Date, y = DIC)) +
    geom_point(size = 1, alpha = 0.75) +
    geom_smooth(method = "loess", color = "grey", se = TRUE) +
    labs(y = expression(paste("DIC (", mu, "mol/kg)")), x = NULL) +
    scale_x_datetime(limits = c(as.POSIXct(start_time, tz = "UTC"),
                                as.POSIXct(end_time, tz = "UTC")),
                     date_breaks = "4 hours", date_labels = "%d-%b %H:%M") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  # Plot for Oxygen (O2uM)
  p3 <- ggplot(data, aes(x = Date, y = O2uM)) +
    geom_point(size = 1, alpha = 0.75) +
    geom_smooth(method = "loess", color = "grey", se = TRUE) +
    labs(y = expression(paste("Oxygen (", mu, "mol/L)")), x = "Date-Time") +
    scale_x_datetime(limits = c(as.POSIXct(start_time, tz = "UTC"),
                                as.POSIXct(end_time, tz = "UTC")),
                     date_breaks = "4 hours", date_labels = "%d-%b %H:%M") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Combine the three plots using gridExtra
  grid.arrange(p1, p2, p3, ncol = 1, heights = c(3, 3, 3),
               top = paste("Station", station_name))
}

# Plot for Station 51
triple_axis_plot(data_51_clean, station_name = "51", start_time = start_time_51, end_time = end_time_51)

# Plot for Station 130
triple_axis_plot(data_130_clean, station_name = "130", start_time = start_time_130, end_time = end_time_130)

# Generalized Additive Model for DIC and O2 -----------------------------------
## We will use the GAM function to fit a linear mixed effects model to the data for each station.
## We are interested in the relationship between O2/Dic and abiotic factors (Salinity, Temp, wind speed) while accounting for the autoregressive structure of the data.
## The goal is to extract the residuals from the model to analyze points in time
## where the DIC/O2 values deviate from the expected values based on the abiotic factors, possibly indicating biological activity.
data_130_model <- data_130_clean %>%
  filter(!is.na(O2uM), !is.na(Temp), !is.na(Salinity), !is.na(Wind.Speed)) %>%
  arrange(Date) %>%
  mutate(
    Time_numeric = as.numeric(Date - min(Date)),
    Salinity_scaled = scale(Salinity)[, 1],
    Temp_scaled = scale(Temp)[, 1],
    Wind_scaled = scale(Wind.Speed)[, 1],
    AR.start = c(TRUE, diff(Date) > 1800)  # mark start of new time block for AR structure (e.g. > 30min gap)
  )

# Fit base model with autocorrelation
# Salinity-only model
mod_salin <- mgcv::bam(O2uM ~ s(Salinity_scaled),
                 data = data_130_model,
                 method = "fREML",
                 rho = 0.9, AR.start = data_130_model$AR.start)

# Salinity + Wind model
mod_salin_wind <- mgcv::bam(O2uM ~ s(Salinity_scaled) + s(Wind_scaled),
                      data = data_130_model,
                      method = "fREML",
                      rho = 0.9, AR.start = data_130_model$AR.start)

# Full model
mod_full <- mgcv::bam(O2uM ~ s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled),
                data = data_130_model,
                method = "fREML",
                rho = 0.9, AR.start = data_130_model$AR.start)

# Compare AIC
AIC(mod_salin, mod_salin_wind, mod_full)

# Cross-validation: time-based blocks
## Split data into 5 time-based blocks
n <- nrow(data_130_model)
fold_size <- floor(n / 5)
breaks <- seq(1, n + 1, by = fold_size)
folds <- map2(breaks[-6], breaks[-1] - 1, ~ .x:.y)

cv_models <- list(mod_salin = "s(Salinity_scaled)",
                  mod_salin_wind = "s(Salinity_scaled) + s(Wind_scaled)",
                  mod_full = "s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled)")

cv_results <- map(cv_models, function(formula_str) {
  rmses <- map_dbl(folds, function(test_idx) {
    train_data <- data_130_model[-test_idx, ]
    test_data <- data_130_model[test_idx, ]

    model <- mgcv::bam(as.formula(paste("O2uM ~", formula_str)),
                 data = train_data,
                 method = "fREML",
                 rho = 0.9,
                 AR.start = train_data$AR.start)

    preds <- predict(model, newdata = test_data)
    sqrt(mean((test_data$O2uM - preds)^2, na.rm = TRUE))
  })
  mean(rmses)
})

cv_results  # RMSE per model

# 5. Residual analysis from best model
data_130_model <- data_130_model %>%
  mutate(
    O2_pred = mgcv::predict.bam(mod_full),
    resid = O2uM - O2_pred
  )

# Smooth residuals over time
gam_resid <- gam(resid ~ s(Time_numeric, k = 10), data = data_130_model)
derivs <- derivatives(gam_resid, term = "s(Time_numeric)")

# Interpolate to match time steps
data_130_model$dO2_resid_dt <- approx(
  x = derivs$`Time_numeric`,
  y = derivs$.derivative,
  xout = data_130_model$Time_numeric,
  rule = 2
)$y

# Plot residuals + derivative (z-scored)
data_130_model <- data_130_model %>%
  mutate(
    resid_z = scale(resid)[, 1],
    dO2_resid_dt_z = scale(dO2_resid_dt)[, 1]
  )

ggplot(data_130_model, aes(x = Date)) +
  geom_line(aes(y = resid_z), color = "forestgreen", size = 1) +
  geom_line(aes(y = dO2_resid_dt_z), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Standardized residual O2 and first derivative",
    subtitle = "green = residuals, purple = d[Residual]/dt",
    y = "Z-score",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_130_model, aes(x = Date)) +
  geom_line(aes(y = resid), color = "forestgreen", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residual O2",
    y = "Residual O2 (μmol/L)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_130_model, aes(x = Date)) +
  geom_line(aes(y = dO2_resid_dt), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "First derivative of residual O2",
    y = "d[Residual O2]/dt (μmol/L/hr)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

# Export output, this can then be used for the GLSAR
data_130_model %>%
  mutate(O2_pred = as.numeric(O2_pred),
         resid = as.numeric(resid),
         dO2_resid_dt = as.numeric(dO2_resid_dt)) %>%
  select(Date, O2uM, Salinity, Temp, Wind.Speed,
         Salinity_scaled, Temp_scaled, Wind_scaled,
         O2_pred, resid, dO2_resid_dt) %>%
  write_csv("data/analysis/O2_model_station130_autocorr.csv")

# DIC
data_130_model_DIC <- data_130_clean %>%
  filter(!is.na(DIC), !is.na(Temp), !is.na(Salinity), !is.na(Wind.Speed)) %>%
  arrange(Date) %>%
  mutate(
    Time_numeric = as.numeric(Date - min(Date)),
    Salinity_scaled = scale(Salinity)[, 1],
    Temp_scaled = scale(Temp)[, 1],
    Wind_scaled = scale(Wind.Speed)[, 1],
    AR.start = c(TRUE, diff(Date) > 1800)  # mark start of new time block for AR structure (e.g. > 30min gap)
  )

# Fit base model with autocorrelation
# Salinity-only model
mod_salin <- mgcv::bam(DIC ~ s(Salinity_scaled),
                       data = data_130_model_DIC,
                       method = "fREML",
                       rho = 0.9, AR.start = data_130_model_DIC$AR.start)

# Salinity + Wind model
mod_salin_wind <- mgcv::bam(DIC ~ s(Salinity_scaled) + s(Wind_scaled),
                            data = data_130_model_DIC,
                            method = "fREML",
                            rho = 0.9, AR.start = data_130_model_DIC$AR.start)

# Full model
mod_full <- mgcv::bam(DIC ~ s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled),
                      data = data_130_model_DIC,
                      method = "fREML",
                      rho = 0.9, AR.start = data_130_model_DIC$AR.start)

# Compare AIC
AIC(mod_salin, mod_salin_wind, mod_full)

# Cross-validation: time-based blocks
## Split data into 5 time-based blocks
n <- nrow(data_130_model_DIC)
fold_size <- floor(n / 5)
breaks <- seq(1, n + 1, by = fold_size)
folds <- map2(breaks[-6], breaks[-1] - 1, ~ .x:.y)

cv_models <- list(mod_salin = "s(Salinity_scaled)",
                  mod_salin_wind = "s(Salinity_scaled) + s(Wind_scaled)",
                  mod_full = "s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled)")

cv_results <- map(cv_models, function(formula_str) {
  rmses <- map_dbl(folds, function(test_idx) {
    train_data <- data_130_model_DIC[-test_idx, ]
    test_data <- data_130_model_DIC[test_idx, ]

    model <- mgcv::bam(as.formula(paste("DIC ~", formula_str)),
                       data = train_data,
                       method = "fREML",
                       rho = 0.9,
                       AR.start = train_data$AR.start)

    preds <- predict(model, newdata = test_data)
    sqrt(mean((test_data$DIC - preds)^2, na.rm = TRUE))
  })
  mean(rmses)
})

cv_results  # RMSE per model

# 5. Residual analysis from best model
data_130_model_DIC <- data_130_model_DIC %>%
  mutate(
    DIC_pred = predict(mod_salin),
    resid = DIC - DIC_pred
  )

# Smooth residuals over time
gam_resid <- gam(resid ~ s(Time_numeric, k = 10), data = data_130_model_DIC)
derivs <- derivatives(gam_resid, term = "s(Time_numeric)")

# Interpolate to match time steps
data_130_model_DIC$dDIC_resid_dt <- approx(
  x = derivs$`Time_numeric`,
  y = derivs$.derivative,
  xout = data_130_model_DIC$Time_numeric,
  rule = 2
)$y

# Plot residuals + derivative (z-scored)
data_130_model_DIC <- data_130_model_DIC %>%
  mutate(
    resid_z = scale(resid)[, 1],
    dDIC_resid_dt_z = scale(dDIC_resid_dt)[, 1]
  )

ggplot(data_130_model_DIC, aes(x = Date)) +
  geom_line(aes(y = resid_z), color = "#3275a8", size = 1) +
  geom_line(aes(y = dDIC_resid_dt_z), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Standardized residual DIC and first derivative",
    subtitle = "green = residuals, purple = d[Residual]/dt",
    y = "Z-score",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_130_model_DIC, aes(x = Date)) +
  geom_line(aes(y = resid), color = "#3275a8", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residual DIC",
    y = "Residual DIC (μmol/L)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_130_model_DIC, aes(x = Date)) +
  geom_line(aes(y = dDIC_resid_dt), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "First derivative of residual DIC",
    y = "d[Residual DIC]/dt (μmol/L/hr)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

# Check Pearson correlation bettween residual DIC and O2
cor.test(data_130_model_DIC$resid, data_130_model$resid, method = "pearson")

# Quantify diel amplitude of residuals
print(range(data_130_model_DIC$resid))
print(range(data_130_model$resid))

# Repeat for Station 51 ---------------------------------------------------
data_51_model <- data_51_clean %>%
  filter(!is.na(O2uM), !is.na(Temp), !is.na(Salinity), !is.na(Wind.Speed)) %>%
  arrange(Date) %>%
  mutate(
    Time_numeric = as.numeric(Date - min(Date)),
    Salinity_scaled = scale(Salinity)[, 1],
    Temp_scaled = scale(Temp)[, 1],
    Wind_scaled = scale(Wind.Speed)[, 1],
    AR.start = c(TRUE, diff(Date) > 1800)  # mark start of new time block for AR structure (e.g. > 30min gap)
  )

# Fit base model with autocorrelation
# Salinity-only model
mod_salin <- mgcv::bam(O2uM ~ s(Salinity_scaled),
                       data = data_51_model,
                       method = "fREML",
                       rho = 0.9, AR.start = data_51_model$AR.start)

# Salinity + Wind model
mod_salin_wind <- mgcv::bam(O2uM ~ s(Salinity_scaled) + s(Wind_scaled),
                            data = data_51_model,
                            method = "fREML",
                            rho = 0.9, AR.start = data_51_model$AR.start)

# Full model
mod_full <- mgcv::bam(O2uM ~ s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled),
                      data = data_51_model,
                      method = "fREML",
                      rho = 0.9, AR.start = data_51_model$AR.start)

# Compare AIC
AIC(mod_salin, mod_salin_wind, mod_full)

# Cross-validation: time-based blocks
## Split data into 5 time-based blocks
n <- nrow(data_51_model)
fold_size <- floor(n / 5)
breaks <- seq(1, n + 1, by = fold_size)
folds <- map2(breaks[-6], breaks[-1] - 1, ~ .x:.y)

cv_models <- list(mod_salin = "s(Salinity_scaled)",
                  mod_salin_wind = "s(Salinity_scaled) + s(Wind_scaled)",
                  mod_full = "s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled)")

cv_results <- map(cv_models, function(formula_str) {
  rmses <- map_dbl(folds, function(test_idx) {
    train_data <- data_51_model[-test_idx, ]
    test_data <- data_51_model[test_idx, ]

    model <- mgcv::bam(as.formula(paste("O2uM ~", formula_str)),
                       data = train_data,
                       method = "fREML",
                       rho = 0.9,
                       AR.start = train_data$AR.start)

    preds <- predict(model, newdata = test_data)
    sqrt(mean((test_data$O2uM - preds)^2, na.rm = TRUE))
  })
  mean(rmses)
})

cv_results  # RMSE per model

# 5. Residual analysis from best model
data_51_model <- data_51_model %>%
  mutate(
    O2_pred = predict(mod_salin_wind),
    resid = O2uM - O2_pred
  )

# Smooth residuals over time
gam_resid <- gam(resid ~ s(Time_numeric, k = 10), data = data_51_model)
derivs <- derivatives(gam_resid, term = "s(Time_numeric)")

# Interpolate to match time steps
data_51_model$dO2_resid_dt <- approx(
  x = derivs$`Time_numeric`,
  y = derivs$.derivative,
  xout = data_51_model$Time_numeric,
  rule = 2
)$y

# Plot residuals + derivative (z-scored)
data_51_model <- data_51_model %>%
  mutate(
    resid_z = scale(resid)[, 1],
    dO2_resid_dt_z = scale(dO2_resid_dt)[, 1]
  )

ggplot(data_51_model, aes(x = Date)) +
  geom_line(aes(y = resid_z), color = "forestgreen", size = 1) +
  geom_line(aes(y = dO2_resid_dt_z), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Standardized residual O2 and first derivative",
    subtitle = "green = residuals, purple = d[Residual]/dt",
    y = "Z-score",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_51_model, aes(x = Date)) +
  geom_line(aes(y = resid), color = "forestgreen", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residual O2",
    y = "Residual O2 (μmol/L)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_51_model, aes(x = Date)) +
  geom_line(aes(y = dO2_resid_dt), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "First derivative of residual O2",
    y = "d[Residual O2]/dt (μmol/L/hr)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

data_51_model_DIC <- data_51_clean %>%
  filter(!is.na(DIC), !is.na(Temp), !is.na(Salinity), !is.na(Wind.Speed)) %>%
  arrange(Date) %>%
  mutate(
    Time_numeric = as.numeric(Date - min(Date)),
    Salinity_scaled = scale(Salinity)[, 1],
    Temp_scaled = scale(Temp)[, 1],
    Wind_scaled = scale(Wind.Speed)[, 1],
    AR.start = c(TRUE, diff(Date) > 1800)  # mark start of new time block for AR structure (e.g. > 30min gap)
  )

# Fit base model with autocorrelation
# Salinity-only model
mod_salin <- mgcv::bam(DIC ~ s(Salinity_scaled),
                       data = data_51_model_DIC,
                       method = "fREML",
                       rho = 0.9, AR.start = data_51_model_DIC$AR.start)

# Salinity + Wind model
mod_salin_wind <- mgcv::bam(DIC ~ s(Salinity_scaled) + s(Wind_scaled),
                            data = data_51_model_DIC,
                            method = "fREML",
                            rho = 0.9, AR.start = data_51_model_DIC$AR.start)

# Full model
mod_full <- mgcv::bam(DIC ~ s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled),
                      data = data_51_model_DIC,
                      method = "fREML",
                      rho = 0.9, AR.start = data_51_model_DIC$AR.start)

# Compare AIC
AIC(mod_salin, mod_salin_wind, mod_full)

# Cross-validation: time-based blocks
## Split data into 5 time-based blocks
n <- nrow(data_51_model_DIC)
fold_size <- floor(n / 5)
breaks <- seq(1, n + 1, by = fold_size)
folds <- map2(breaks[-6], breaks[-1] - 1, ~ .x:.y)

cv_models <- list(mod_salin = "s(Salinity_scaled)",
                  mod_salin_wind = "s(Salinity_scaled) + s(Wind_scaled)",
                  mod_full = "s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled)")

cv_results <- map(cv_models, function(formula_str) {
  rmses <- map_dbl(folds, function(test_idx) {
    train_data <- data_51_model_DIC[-test_idx, ]
    test_data <- data_51_model_DIC[test_idx, ]

    model <- mgcv::bam(as.formula(paste("DIC ~", formula_str)),
                       data = train_data,
                       method = "fREML",
                       rho = 0.9,
                       AR.start = train_data$AR.start)

    preds <- predict(model, newdata = test_data)
    sqrt(mean((test_data$DIC - preds)^2, na.rm = TRUE))
  })
  mean(rmses)
})

cv_results  # RMSE per model

# 5. Residual analysis from best model
data_51_model_DIC <- data_51_model_DIC %>%
  mutate(
    DIC_pred = predict(mod_salin),
    resid = DIC - DIC_pred
  )

# Smooth residuals over time
gam_resid <- gam(resid ~ s(Time_numeric, k = 10), data = data_51_model_DIC)
derivs <- derivatives(gam_resid, term = "s(Time_numeric)")

# Interpolate to match time steps
data_51_model_DIC$dDIC_resid_dt <- approx(
  x = derivs$`Time_numeric`,
  y = derivs$.derivative,
  xout = data_51_model_DIC$Time_numeric,
  rule = 2
)$y

# Plot residuals + derivative (z-scored)
data_51_model_DIC <- data_51_model_DIC %>%
  mutate(
    resid_z = scale(resid)[, 1],
    dDIC_resid_dt_z = scale(dDIC_resid_dt)[, 1]
  )

ggplot(data_51_model_DIC, aes(x = Date)) +
  geom_line(aes(y = resid_z), color = "#3275a8", size = 1) +
  geom_line(aes(y = dDIC_resid_dt_z), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Standardized residual DIC and first derivative",
    subtitle = "green = residuals, purple = d[Residual]/dt",
    y = "Z-score",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_51_model_DIC, aes(x = Date)) +
  geom_line(aes(y = resid), color = "#3275a8", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residual DIC",
    y = "Residual DIC (μmol/L)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_51_model_DIC, aes(x = Date)) +
  geom_line(aes(y = dDIC_resid_dt), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "First derivative of residual DIC",
    y = "d[Residual DIC]/dt (μmol/L/hr)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

# DIC
data_51_model_DIC <- data_51_clean %>%
  filter(!is.na(DIC), !is.na(Temp), !is.na(Salinity), !is.na(Wind.Speed)) %>%
  arrange(Date) %>%
  mutate(
    Time_numeric = as.numeric(Date - min(Date)),
    Salinity_scaled = scale(Salinity)[, 1],
    Temp_scaled = scale(Temp)[, 1],
    Wind_scaled = scale(Wind.Speed)[, 1],
    AR.start = c(TRUE, diff(Date) > 1800)  # mark start of new time block for AR structure (e.g. > 30min gap)
  )

# Fit base model with autocorrelation
# Salinity-only model
mod_salin <- mgcv::bam(DIC ~ s(Salinity_scaled),
                       data = data_51_model_DIC,
                       method = "fREML",
                       rho = 0.9, AR.start = data_51_model_DIC$AR.start)

# Salinity + Wind model
mod_salin_wind <- mgcv::bam(DIC ~ s(Salinity_scaled) + s(Wind_scaled),
                            data = data_51_model_DIC,
                            method = "fREML",
                            rho = 0.9, AR.start = data_51_model_DIC$AR.start)

# Full model
mod_full <- mgcv::bam(DIC ~ s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled),
                      data = data_51_model_DIC,
                      method = "fREML",
                      rho = 0.9, AR.start = data_51_model_DIC$AR.start)

# Compare AIC
AIC(mod_salin, mod_salin_wind, mod_full)

# Cross-validation: time-based blocks
## Split data into 5 time-based blocks
n <- nrow(data_51_model_DIC)
fold_size <- floor(n / 5)
breaks <- seq(1, n + 1, by = fold_size)
folds <- map2(breaks[-6], breaks[-1] - 1, ~ .x:.y)

cv_models <- list(mod_salin = "s(Salinity_scaled)",
                  mod_salin_wind = "s(Salinity_scaled) + s(Wind_scaled)",
                  mod_full = "s(Salinity_scaled) + s(Temp_scaled) + s(Wind_scaled)")

cv_results <- map(cv_models, function(formula_str) {
  rmses <- map_dbl(folds, function(test_idx) {
    train_data <- data_51_model_DIC[-test_idx, ]
    test_data <- data_51_model_DIC[test_idx, ]

    model <- mgcv::bam(as.formula(paste("DIC ~", formula_str)),
                       data = train_data,
                       method = "fREML",
                       rho = 0.9,
                       AR.start = train_data$AR.start)

    preds <- predict(model, newdata = test_data)
    sqrt(mean((test_data$DIC - preds)^2, na.rm = TRUE))
  })
  mean(rmses)
})

cv_results  # RMSE per model

# 5. Residual analysis from best model
data_51_model_DIC <- data_51_model_DIC %>%
  mutate(
    DIC_pred = predict(mod_full),
    resid = DIC - DIC_pred
  )

# Smooth residuals over time
gam_resid <- gam(resid ~ s(Time_numeric, k = 10), data = data_51_model_DIC)
derivs <- derivatives(gam_resid, term = "s(Time_numeric)")

# Interpolate to match time steps
data_51_model_DIC$dDIC_resid_dt <- approx(
  x = derivs$`Time_numeric`,
  y = derivs$.derivative,
  xout = data_51_model_DIC$Time_numeric,
  rule = 2
)$y

# Plot residuals + derivative (z-scored)
data_51_model_DIC <- data_51_model_DIC %>%
  mutate(
    resid_z = scale(resid)[, 1],
    dDIC_resid_dt_z = scale(dDIC_resid_dt)[, 1]
  )

ggplot(data_51_model_DIC, aes(x = Date)) +
  geom_line(aes(y = resid_z), color = "#3275a8", size = 1) +
  geom_line(aes(y = dDIC_resid_dt_z), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Standardized residual DIC and first derivative",
    subtitle = "green = residuals, purple = d[Residual]/dt",
    y = "Z-score",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_51_model_DIC, aes(x = Date)) +
  geom_line(aes(y = resid), color = "#3275a8", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Residual DIC",
    y = "Residual DIC (μmol/L)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

ggplot(data_51_model_DIC, aes(x = Date)) +
  geom_line(aes(y = dDIC_resid_dt), color = "purple", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "First derivative of residual DIC",
    y = "d[Residual DIC]/dt (μmol/L/hr)",
    x = "Time"
  ) +
  theme_minimal(base_size = 13)

# Check Pearson correlation bettween residual DIC and O2
cor.test(data_51_model_DIC$resid, data_51_model$resid, method = "pearson")

# Quantify diel amplitude of residuals
print(range(data_51_model_DIC$resid))
print(range(data_51_model$resid))

# Function to plot z-scaled O2/DIC and flipped Salinity ---------------------------
## Function to add day moment column based on sunrise/sunset times
get_station_sun_times <- function(lat, lon, dates) {
  sun_times <- getSunlightTimes(
    date = dates,
    lat = lat,
    lon = lon,
    keep = c("sunrise", "sunset",
             "dawn", "dusk",
             "nauticalDawn", "nauticalDusk",
             "night", "nightEnd"),
    tz = "UTC"
  ) %>%
    mutate(date = as.Date(date)) %>%
    select(date, night, nightEnd, nauticalDawn, nauticalDusk,
           dawn, dusk, sunrise, sunset)

  return(sun_times)
}

lat_130 <- median(data_130_clean$Latitude, na.rm = TRUE)
lon_130 <- median(data_130_clean$Longitude, na.rm = TRUE)

lat_51 <- median(data_51_clean$Latitude, na.rm = TRUE)
lon_51 <- median(data_51_clean$Longitude, na.rm = TRUE)

sun_times_130 <- get_station_sun_times(
  lat = lat_130,
  lon = lon_130,
  dates = unique(as.Date(data_130_clean$Date))
)

sun_times_51 <- get_station_sun_times(
  lat = lat_51,
  lon = lon_51,
  dates = unique(as.Date(data_51_clean$Date))
)

data_51_clean <- data_51_clean %>%
  mutate(date = as.Date(Date)) %>%
  left_join(sun_times_51, by = "date")

data_51_clean <- data_51_clean %>%
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
  ungroup()

data_130_clean <- data_130_clean %>%
  mutate(date = as.Date(Date)) %>%
  left_join(sun_times_130, by = "date")

data_130_clean <- data_130_clean %>%
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
  ungroup()

# Cleanup
data_51_clean <- data_51_clean %>%
  select(-date, -night, -nightEnd, -nauticalDawn, -nauticalDusk,
         -dawn, -dusk, -sunrise, -sunset)
data_130_clean <- data_130_clean %>%
  select(-date, -night, -nightEnd, -nauticalDawn, -nauticalDusk,
         -dawn, -dusk, -sunrise, -sunset)

# Plot z-scaled O2 vs Salinity and shade by day moment
plot_gamm_scaled_O2_sal <- function(data, station_name, start_time, end_time, smoother_span = 0.5) {
  # Define light phase colors
  light_colors <- c(
    "Night" = "#d9d9d9",
    "Astronomical twilight" = "#ffb347",
    "Nautical twilight" = "#ffc870",
    "Civil twilight" = "#ffe0a3",
    "Day" = "#ffffb3"
  )

  # ---- Top plot: Residuals and d[Residual]/dt ----
  res_deriv_plot <- ggplot(data, aes(x = Date)) +
    geom_rect(aes(xmin = Date, xmax = lead(Date), ymin = -Inf, ymax = Inf, fill = day_moment),
              alpha = 0.4, color = NA) +
    geom_line(aes(y = resid_z), color = "forestgreen", size = 0.9) +
    geom_line(aes(y = dO2_resid_dt_z), color = "purple", size = 0.9, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = light_colors) +
    labs(
      y = "Z-score", x = NULL
    ) +
    scale_x_datetime(
      limits = c(as.POSIXct(start_time, tz = "UTC"),
                 as.POSIXct(end_time, tz = "UTC")),
      date_breaks = "4 hours", date_labels = "%d-%b %H:%M"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")

  # ---- Bottom plot: Salinity and O₂ (Z-scored) for context ----
  time_data_long <- data %>%
    dplyr::select(Date, O2uM, Salinity) %>%
    mutate(
      O2_z = scale(O2uM)[, 1],
      Salinity_z = scale(Salinity)[, 1]
    ) %>%
    tidyr::pivot_longer(cols = c(O2_z, Salinity_z), names_to = "Variable", values_to = "Z_value") %>%
    mutate(
      Variable = recode(Variable,
                        "O2_z" = "O₂",
                        "Salinity_z" = "Salinity")
    ) %>%
    left_join(data %>% select(Date, day_moment), by = "Date")

  sal_o2_plot <- ggplot(time_data_long, aes(x = Date, y = Z_value)) +
    geom_rect(aes(xmin = Date, xmax = lead(Date), ymin = -Inf, ymax = Inf, fill = day_moment),
              alpha = 0.4, color = NA) +
    geom_point(aes(color = Variable), size = 1, alpha = 0.75) +
    geom_smooth(aes(color = Variable), method = "loess", se = TRUE, span = smoother_span) +
    scale_color_manual(values = c("O₂" = "darkgreen", "Salinity" = "#ff7f0e")) +
    scale_fill_manual(values = light_colors) +
    labs(y = "Z-score\n(mean = 0, SD = 1)", color = "", fill = NULL) +
    scale_x_datetime(
      limits = c(as.POSIXct(start_time, tz = "UTC"),
                 as.POSIXct(end_time, tz = "UTC")),
      date_breaks = "4 hours", date_labels = "%d-%b %H:%M"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Combine and add title
  combined <- cowplot::plot_grid(
    res_deriv_plot,
    sal_o2_plot,
    ncol = 1,
    rel_heights = c(1.2, 2)
  )

  final_plot <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(paste("Station", station_name), fontface = "bold", size = 14),
    combined,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )

  return(final_plot)
}

# Plot z-scaled DIC vs -Salinity and shade by day moment
plot_gamm_scaled_DIC_sal <- function(data, station_name, start_time, end_time, smoother_span = 0.5) {
  # Define light phase colors
  light_colors <- c(
    "Night" = "#d9d9d9",
    "Astronomical twilight" = "#ffb347",
    "Nautical twilight" = "#ffc870",
    "Civil twilight" = "#ffe0a3",
    "Day" = "#ffffb3"
  )

  # ---- Top plot: Residuals and d[Residual]/dt ----
  res_deriv_plot <- ggplot(data, aes(x = Date)) +
    geom_rect(aes(xmin = Date, xmax = lead(Date), ymin = -Inf, ymax = Inf, fill = day_moment),
              alpha = 0.4, color = NA) +
    geom_line(aes(y = resid_z), color = "#3275a8", size = 0.9) +
    geom_line(aes(y = dDIC_resid_dt_z), color = "purple", size = 0.9, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = light_colors) +
    labs(
      y = "Z-score", x = NULL
    ) +
    scale_x_datetime(
      limits = c(as.POSIXct(start_time, tz = "UTC"),
                 as.POSIXct(end_time, tz = "UTC")),
      date_breaks = "4 hours", date_labels = "%d-%b %H:%M"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")

  # ---- Bottom plot: Salinity and O₂ (Z-scored) for context ----
  time_data_long <- data %>%
    dplyr::select(Date, DIC, Salinity) %>%
    mutate(
      DIC_z = scale(DIC)[, 1],
      Salinity_z = scale(Salinity)[, 1] * -1
    ) %>%
    tidyr::pivot_longer(cols = c(DIC_z, Salinity_z), names_to = "Variable", values_to = "Z_value") %>%
    mutate(
      Variable = recode(Variable,
                        "DIC_z" = "DIC",
                        "Salinity_z" = "-Salinity")
    ) %>%
    left_join(data %>% select(Date, day_moment), by = "Date")

  sal_DIC_plot <- ggplot(time_data_long, aes(x = Date, y = Z_value)) +
    geom_rect(aes(xmin = Date, xmax = lead(Date), ymin = -Inf, ymax = Inf, fill = day_moment),
              alpha = 0.4, color = NA) +
    geom_point(aes(color = Variable), size = 1, alpha = 0.75) +
    geom_smooth(aes(color = Variable), method = "loess", se = TRUE, span = smoother_span) +
    scale_color_manual(values = c("DIC" = "#0d5f9e", "-Salinity" = "#ff7f0e")) +
    scale_fill_manual(values = light_colors) +
    labs(y = "Z-score\n(mean = 0, SD = 1)", color = "", fill = NULL) +
    scale_x_datetime(
      limits = c(as.POSIXct(start_time, tz = "UTC"),
                 as.POSIXct(end_time, tz = "UTC")),
      date_breaks = "4 hours", date_labels = "%d-%b %H:%M"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Combine and add title
  combined <- cowplot::plot_grid(
    res_deriv_plot,
    sal_DIC_plot,
    ncol = 1,
    rel_heights = c(1.2, 2)
  )

  final_plot <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(paste("Station", station_name), fontface = "bold", size = 14),
    combined,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )

  return(final_plot)
}

# Define time window for the figure
start_time <- "2023-04-20 08:00:00"
end_time   <- "2023-04-21 08:00:00"

data_130_model <- data_130_model %>%
  left_join(data_130_clean %>% select(Date, day_moment), by = "Date")

# Generate and print the O2 plot for station 130
plot_130 <- plot_gamm_scaled_O2_sal(
  data = data_130_model,
  station_name = "130",
  start_time = start_time,
  end_time = end_time,
  smoother_span = 0.4  # adjust for visual smoothness if needed
)

ggsave("figures/environmental/Zscaled_O2_vs_Salinity_Station130.svg", plot_130, width = 16, height = 12, units = "cm")
plot_130

data_130_model_DIC <- data_130_model_DIC %>%
  left_join(data_130_clean %>% select(Date, day_moment), by = "Date")

# Generate and print the DIC plot for station 130
plot_130 <- plot_gamm_scaled_DIC_sal(
  data = data_130_model_DIC,
  station_name = "130",
  start_time = start_time,
  end_time = end_time,
  smoother_span = 0.4  # adjust for visual smoothness if needed
)

ggsave("figures/environmental/Zscaled_DIC_vs_Salinity_Station130.svg", plot_130, width = 16, height = 12, units = "cm")
plot_130

# For Station 51
start_time_51 <- "2023-04-18 11:00:00"
end_time_51   <- "2023-04-19 08:00:00"

data_51_model <- data_51_model %>%
  left_join(data_51_clean %>% select(Date, day_moment), by = "Date")

# Generate and print the O2 plot for station 51
plot_51 <- plot_gamm_scaled_O2_sal(
  data = data_51_model,
  station_name = "51",
  start_time = start_time_51,
  end_time = end_time_51,
  smoother_span = 0.4  # adjust for visual smoothness if needed
)

ggsave("figures/environmental/Zscaled_O2_vs_Salinity_Station51.svg", plot_51, width = 16, height = 12, units = "cm")
plot_51

data_51_model_DIC <- data_51_model_DIC %>%
  left_join(data_51_clean %>% select(Date, day_moment), by = "Date")

# Generate and print the DIC plot for station 51
plot_51 <- plot_gamm_scaled_DIC_sal(
  data = data_51_model_DIC,
  station_name = "51",
  start_time = start_time_51,
  end_time = end_time_51,
  smoother_span = 0.4  # adjust for visual smoothness if needed
)

ggsave("figures/environmental/Zscaled_DIC_vs_Salinity_Station51.svg", plot_51, width = 16, height = 12, units = "cm")
plot_51

# Function to create a single plot for a given parameter ------------------------
create_plot <- function(data, y_var, y_label, color, y_limits, start_time, end_time, show_x_axis = TRUE, smoother_span = 0.5) {
  p <- ggplot(data, aes(x = Date, y = !!sym(y_var))) +
    geom_point(size = 1, alpha = 0.75, color = "black") +
    geom_smooth(method = "loess", color = color, se = TRUE, span = smoother_span) +
    labs(y = y_label, x = NULL) +
    scale_x_datetime(limits = c(as.POSIXct(start_time, tz = "UTC"),
                                as.POSIXct(end_time, tz = "UTC")),
                     date_breaks = "4 hours", date_labels = "%d-%b %H:%M",
                     date_minor_breaks = "2 hours") +
    scale_y_continuous(limits = y_limits) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = if (show_x_axis) element_text(angle = 45, hjust = 1) else element_blank(),
          axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
          axis.title.x = element_blank())

  return(p)
}

# Dynamic y-axis limits based on min/max values ---------------------------------
y_limits_pco2 <- c(
  min(c(data_51_clean$pCO2, data_130_clean$pCO2), na.rm = TRUE),
  max(c(data_51_clean$pCO2, data_130_clean$pCO2), na.rm = TRUE)
)

y_limits_dic <- c(
  min(c(data_51_clean$DIC, data_130_clean$DIC), na.rm = TRUE),
  max(c(data_51_clean$DIC, data_130_clean$DIC), na.rm = TRUE)
)

y_limits_o2 <- c(
  min(c(data_51_clean$O2uM, data_130_clean$O2uM), na.rm = TRUE),
  max(c(data_51_clean$O2uM, data_130_clean$O2uM), na.rm = TRUE)
)

y_limits_ta <- c(
  min(c(data_51_clean$TA, data_130_clean$TA), na.rm = TRUE),
  max(c(data_51_clean$TA, data_130_clean$TA), na.rm = TRUE)
)

# Create plots for Station 51 ---------------------------------------------------
pco2_51 <- create_plot(data_51_clean, "pCO2", expression(paste("pCO"[2], " (", mu, "atm)")), "grey", y_limits_pco2, start_time_51, end_time_51, show_x_axis = FALSE)
dic_51 <- create_plot(data_51_clean, "DIC", expression(paste("DIC (", mu, "mol/kg)")), "grey", y_limits_dic, start_time_51, end_time_51,  show_x_axis = FALSE)
o2_51 <- create_plot(data_51_clean, "O2uM", expression(paste("Oxygen (", mu, "mol/L)")), "grey", y_limits_o2, start_time_51, end_time_51, show_x_axis = TRUE)
ta_51 <- create_plot(data_51_clean, "TA", expression(paste("TA (", mu, "mol/kg)")), "grey", y_limits_ta, start_time_51, end_time_51,  show_x_axis = TRUE)

# Create plots for Station 130 --------------------------------------------------
pco2_130 <- create_plot(data_130_clean, "pCO2", expression(paste("pCO"[2], " (", mu, "atm)")), "grey", y_limits_pco2, start_time_130, end_time_130, show_x_axis = FALSE)
dic_130 <- create_plot(data_130_clean, "DIC", expression(paste("DIC (", mu, "mol/kg)")), "grey", y_limits_dic, start_time_130, end_time_130, show_x_axis = FALSE)
o2_130 <- create_plot(data_130_clean, "O2uM", expression(paste("Oxygen (", mu, "mol/L)")), "grey", y_limits_o2, start_time_130, end_time_130, show_x_axis = TRUE)
ta_130 <- create_plot(data_130_clean, "TA", expression(paste("TA (", mu, "mol/kg)")), "grey", y_limits_ta, start_time_130, end_time_130,  show_x_axis = TRUE)

## Extract the smoothers for station 130 for subsequent analysis
extract_smoother_values <- function(plot) {
  plot_data <- ggplot_build(plot)
  smoothed_values <- plot_data$data[[2]]  # Layer 2 is geom_smooth
  return(smoothed_values)
}

# Process and clean smoother data
process_smoother_data <- function(smoother_data) {
  smoother_data_clean <- smoother_data %>%
    dplyr::select(x, y, ymin, ymax, se) %>%   # Only keep relevant columns
    dplyr::mutate(Date = as.POSIXct(x, origin = "1970-01-01", tz = "UTC")) %>%  # Convert epoch to Date.Time
    dplyr::select(Date, y, ymin, ymax, se)
  return(smoother_data_clean)
}

# Extract smoother values
dic_smoother_130 <- extract_smoother_values(dic_130)
o2_smoother_130 <- extract_smoother_values(o2_130)

dic_smoother_51 <- extract_smoother_values(dic_51)
o2_smoother_51 <- extract_smoother_values(o2_51)

# Process smoother data
dic_smoother_130_clean <- process_smoother_data(dic_smoother_130)
o2_smoother_130_clean <- process_smoother_data(o2_smoother_130)

dic_smoother_51_clean <- process_smoother_data(dic_smoother_51)
o2_smoother_51_clean <- process_smoother_data(o2_smoother_51)

# Save smoother values to CSV
write.csv(dic_smoother_130_clean, "data/analysis/DIC_smoother_station_130.csv", row.names = FALSE)
write.csv(o2_smoother_130_clean, "data/analysis/O2_smoother_station_130.csv", row.names = FALSE)

write.csv(dic_smoother_51_clean, "data/analysis/DIC_smoother_station_51.csv", row.names = FALSE)
write.csv(o2_smoother_51_clean, "data/analysis/O2_smoother_station_51.csv", row.names = FALSE)

# Combine plots into two columns with equal heights -----------------------------
left_column <- plot_grid(pco2_51, dic_51, o2_51, ncol = 1, align = "v", rel_heights = c(1, 1, 1))
right_column <- plot_grid(pco2_130, dic_130, o2_130, ncol = 1, align = "v", rel_heights = c(1, 1, 1))

final_plot <- plot_grid(
  left_column,
  right_column,
  ncol = 2,
  labels = c("Station 51", "Station 130"),
  label_size = 11
)

# Save the final figure as an SVG file -----------------------------------------
output_dir <- "figures/environmental/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_path <- file.path(output_dir, "carbonate_chemistry.svg")
ggsave(output_path, plot = final_plot, width = 18, height = 18, units = "cm", device = "svg")

# Display the final plot
print(final_plot)

# Combine TA plots into two columns with equal heights -----------------------------
final_plot <- plot_grid(
  ta_51,
  ta_130,
  ncol = 2,
  labels = c("Station 51", "Station 130"),
  label_size = 11
)

# Save the final figure as an SVG file -----------------------------------------
output_path <- file.path(output_dir, "TA.svg")
ggsave(output_path, plot = final_plot, width = 18, height = 6, units = "cm", device = "svg")

print(final_plot)

# Linear model to describe the relation ---------------------------------
# Combine data from both stations
combined_data <- rbind(
  data.frame(Station = "51", Time = data_51_clean$Date, DIC = data_51_clean$DIC),
  data.frame(Station = "130", Time = data_130_clean$Date, DIC = data_130_clean$DIC)
)

# Fit the model for DIC
model <- lm(DIC ~ Station * Time, data = combined_data)

# Summary of the model
summary(model)

# ANOVA for the model
anova(model)

# Find DIC maximum and minimum for Station 51
max_DIC_51 <- data_51_clean[which.max(data_51_clean$DIC), ]
min_DIC_51 <- data_51_clean[which.min(data_51_clean$DIC), ]

cat(sprintf("For station 51, the maximum DIC value is %.2f at %s\n", max_DIC_51$DIC, max_DIC_51$Date))
cat(sprintf("For station 51, the minimum DIC value is %.2f at %s\n", min_DIC_51$DIC, min_DIC_51$Date))

# Find DIC maximum and minimum for Station 130
max_DIC_130 <- data_130_clean[which.max(data_130_clean$DIC), ]
min_DIC_130 <- data_130_clean[which.min(data_130_clean$DIC), ]

cat(sprintf("For station 130, the maximum DIC value is %.2f at %s\n", max_DIC_130$DIC, max_DIC_130$Date))
cat(sprintf("For station 130, the minimum DIC value is %.2f at %s\n", min_DIC_130$DIC, min_DIC_130$Date))

# Fit the model for O2
combined_data <- rbind(
  data.frame(Station = "51", Time = data_51_clean$Date, O2 = data_51_clean$O2uM),
  data.frame(Station = "130", Time = data_130_clean$Date, O2 = data_130_clean$O2uM)
)

model <- lm(O2 ~ Station * Time, data = combined_data)

# Summary of the model
summary(model)

# ANOVA for the model
anova(model)

# Find O2 maximum and minimum for Station 51
max_O2_51 <- data_51_clean[which.max(data_51_clean$O2uM), ]
min_O2_51 <- data_51_clean[which.min(data_51_clean$O2uM), ]

cat(sprintf("For station 51, the maximum O2 value is %.2f μmol/L at %s\n", max_O2_51$O2uM, max_O2_51$Date))
cat(sprintf("For station 51, the minimum O2 value is %.2f μmol/L at %s\n", min_O2_51$O2uM, min_O2_51$Date))

# Find O2 maximum and minimum for Station 130
max_O2_130 <- data_130_clean[which.max(data_130_clean$O2uM), ]
min_O2_130 <- data_130_clean[which.min(data_130_clean$O2uM), ]

cat(sprintf("For station 130, the maximum O2 value is %.2f μmol/L at %s\n", max_O2_130$O2uM, max_O2_130$Date))
cat(sprintf("For station 130, the minimum O2 value is %.2f μmol/L at %s\n", min_O2_130$O2uM, min_O2_130$Date))

# Correlation between O2 and DIC?
# Combine data from both stations with O2 and DIC
combined_corr_data <- rbind(
  data.frame(Station = "51", O2 = data_51_clean$O2uM, DIC = data_51_clean$DIC),
  data.frame(Station = "130", O2 = data_130_clean$O2uM, DIC = data_130_clean$DIC)
)

# Create scatter plot with regression lines and correlation coefficients
correlation_plot <- ggplot(combined_corr_data, aes(x = O2, y = DIC)) +
  geom_point(alpha = 0.6, color = "grey") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           method = "pearson", label.x = min(combined_corr_data$O2), label.y = max(combined_corr_data$DIC), size = 4) +
  facet_wrap(~ Station, scales = "free", labeller = labeller(Station = c("51" = "Station 51", "130" = "Station 130"))) +
  labs(x = expression(paste("Oxygen (", mu, "mol/L)")),
       y = expression(paste("DIC (", mu, "mol/kg)"))) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 11, face = "bold"))

correlation_plot

# Fit the linear model: DIC as a function of O2
model <- lm(DIC ~ O2uM, data = data_130_clean)

# Extract coefficients
slope <- coef(model)["O2uM"]
intercept <- coef(model)["(Intercept)"]

# Extract R-squared and p-value
summary_model <- summary(model)
r_value <- sqrt(summary_model$r.squared)  # R value (square root of R²)
p_value <- summary_model$coefficients[2, 4]  # P-value for O2 coefficient

# Print results
cat(sprintf("Equation: DIC = %.3f * O2 + %.3f\n", slope, intercept))
cat(sprintf("R-value: %.3f\n", r_value))
cat(sprintf("P-value: %.3e\n", p_value))

# Save the figure as an SVG file
output_path <- file.path(output_dir, "O2_DIC_correlation.svg")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "svg")
output_path <- file.path(output_dir, "O2_DIC_correlation.png")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "png")

# Correlation between TA and DIC?
# Combine data from both stations with O2 and DIC
combined_corr_data <- rbind(
  data.frame(Station = "51", TA = data_51_clean$TA, DIC = data_51_clean$DIC),
  data.frame(Station = "130", TA = data_130_clean$TA, DIC = data_130_clean$DIC)
)

correlation_plot <- ggplot(combined_corr_data, aes(x = TA, y = DIC)) +
  geom_point(alpha = 0.6, color = "grey") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           method = "pearson", label.x = min(combined_corr_data$TA), label.y = max(combined_corr_data$DIC), size = 4) +
  facet_wrap(~ Station, scales = "free", labeller = labeller(Station = c("51" = "Station 51", "130" = "Station 130"))) +
  labs(x = expression(paste("TA (", mu, "mol/L)")),
       y = expression(paste("DIC (", mu, "mol/kg)"))) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 11, face = "bold"))

correlation_plot

# Save the figure as an SVG file
output_path <- file.path(output_dir, "TA_DIC_correlation.svg")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "svg")
output_path <- file.path(output_dir, "TA_DIC_correlation.png")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "png")

# Correlation between TA and Salinity?
# Combine data from both stations with O2 and DIC
combined_corr_data <- rbind(
  data.frame(Station = "51", TA = data_51_clean$TA, Salinity = data_51_clean$Salinity),
  data.frame(Station = "130", TA = data_130_clean$TA, Salinity = data_130_clean$Salinity)
)

correlation_plot <- ggplot(combined_corr_data, aes(x = TA, y = Salinity)) +
  geom_point(alpha = 0.6, color = "grey") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           method = "pearson", label.x = min(combined_corr_data$TA), label.y = max(combined_corr_data$Salinity), size = 4) +
  facet_wrap(~ Station, scales = "free", labeller = labeller(Station = c("51" = "Station 51", "130" = "Station 130"))) +
  labs(x = expression(paste("TA (", mu, "mol/L)")),
       y = expression("Salinity (PSU)")) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 11, face = "bold"))

correlation_plot

# Save the figure as an SVG file
output_path <- file.path(output_dir, "TA_Salinity_correlation.svg")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "svg")
output_path <- file.path(output_dir, "TA_Salinity_correlation.png")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "png")
