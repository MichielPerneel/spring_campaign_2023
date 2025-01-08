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
library(readr)
library(ggpubr)

Sys.setenv(TZ='UTC')


# read data ---------------------------------------------------------------

SiSt_202304 <- read.csv("data/raw/11SS20230401.csv",sep = ",", dec=".", header=TRUE )

SiSt_202304$Date.Time <- as.POSIXct(strptime(SiSt_202304$Date.Time, format="%Y-%m-%dT%H:%M:%S", tz="UTC"), format="%Y-%m-%d %H:%M:%S",tz="UTC")

#salinity - alkalinity relationship y = -43.158x+3827.1 all data

TA_sist2 <- -43.158*SiSt_202304$P_sal..psu.+3827.1

Spot_samples <- read.csv("data/raw/Spot_carb_jn2023.csv",sep = ",", dec=".", header=TRUE )

Spot_samples$Date.Time <- as.POSIXct(strptime(Spot_samples$Date.Time, format="%Y-%m-%d %H:%M:%S", tz="UTC"), format="%Y-%m-%d %H:%M:%S",tz="UTC")

# Extract data for Station 51: 2023-04-18 13:00 to 2023-05-19 09:00
data_51 <- SiSt_202304 %>%
  filter(Date.Time >= as.POSIXct("2023-04-18 13:00:00", tz = "UTC") &
           Date.Time <= as.POSIXct("2023-04-19 09:30:00", tz = "UTC"))

# Extract data for Station 130: 2023-04-20 10:00 to 2023-04-21 10:00
data_130 <- SiSt_202304 %>%
  filter(Date.Time >= as.POSIXct("2023-04-20 10:00:00", tz = "UTC") &
           Date.Time <= as.POSIXct("2023-04-21 10:30:00", tz = "UTC"))

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

# reconstruct carbonate system parameters using seacarb-------------------------------------------------------------

#calculate DIC UW based on pCO2 and reconstructed TA_UW2

# average of PO4 and Si based on spot samples form JN 2023 cruise
sample_metadata <- read_csv("data/samples_env.csv")
Sil <- mean(na.omit(sample_metadata$Si))*10^(-6)
PO4 <- mean(na.omit(sample_metadata$PO4))*10^(-6)

#calculate DIC from pCO2 record and reconstructed TA from sal (from TA/Sal of all data!)

SiSt_clean <- data.frame(SiSt_202304$Date.Time, SiSt_202304$P_sal..psu., SiSt_202304$Temp..degC., SiSt_202304$pCO2..uatm., TA_sist2)

SiSt_clean <- na.omit(SiSt_clean)


DIC_UW <- carb(flag=24, var1=SiSt_clean$SiSt_202304.pCO2..uatm., var2=SiSt_clean$TA_sist2*10^(-6), S=SiSt_clean$SiSt_202304.P_sal..psu., T=SiSt_clean$SiSt_202304.Temp..degC., P=0, Patm=1.0, Pt=PO4, Sit=Sil, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74")


plot_go2<-ggplot()   +
  geom_point(aes(x = SiSt_clean$SiSt_202304.Date.Time, y =DIC_UW$DIC*10^6, colour = "DIC_UW"), size = 0.2, alpha = 0.75) +
  ylab("DICumol/kg") + xlab("")+ylim(2000,2500)+xlim(as.POSIXct("2023-04-17 06:00:00"), as.POSIXct("2023-04-21 23:59:59"))

plot_go2

# Per station -------------------------------------------------------------
# Recalculate DIC for each station based on pCO2 and reconstructed TA
data_51_clean <- data.frame(Date.Time = data_51$Date.Time,
                            Salinity = data_51$P_sal..psu.,
                            Temp = data_51$Temp..degC.,
                            pCO2 = data_51$pCO2..uatm.,
                            TA = -43.158 * data_51$P_sal..psu. + 3827.1)

data_130_clean <- data.frame(Date.Time = data_130$Date.Time,
                             Salinity = data_130$P_sal..psu.,
                             Temp = data_130$Temp..degC.,
                             pCO2 = data_130$pCO2..uatm.,
                             TA = -43.158 * data_130$P_sal..psu. + 3827.1)

# Omit NAs
data_51_clean <- na.omit(data_51_clean)
data_130_clean <- na.omit(data_130_clean)

# Calculate DIC using seacarb
DIC_51 <- carb(flag = 24, var1 = data_51_clean$pCO2, var2 = data_51_clean$TA * 1e-6,
               S = data_51_clean$Salinity, T = data_51_clean$Temp, P = 0, Patm = 1.0,
               Pt = PO4, Sit = Sil, pHscale = "T", kf = "pf", k1k2 = "l", ks = "d", b = "u74",
               )

DIC_130 <- carb(flag = 24, var1 = data_130_clean$pCO2, var2 = data_130_clean$TA * 1e-6,
                S = data_130_clean$Salinity, T = data_130_clean$Temp, P = 0, Patm = 1.0,
                Pt = PO4, Sit = Sil, pHscale = "T", kf = "pf", k1k2 = "l", ks = "d", b = "u74")

# Add DIC to data frames
data_51_clean$DIC <- DIC_51$DIC * 1e6
data_130_clean$DIC <- DIC_130$DIC * 1e6

# Function to match oxygen data based on Date-Time -----------------------------
match_oxygen <- function(data, oxygen_data) {
  data %>%
    rowwise() %>%
    mutate(
      O2uM = oxygen_data$O2uM[which.min(abs(difftime(oxygen_data$Date, Date.Time, units = "secs")))]
    )
}

# Match oxygen data for Station 51
data_51_clean <- match_oxygen(data_51_clean, oxygen_data)

# Match oxygen data for Station 130
data_130_clean <- match_oxygen(data_130_clean, oxygen_data)

# Create plot for pCO2
pco2_plot <- ggplot() +
  geom_point(data = data_51_clean, aes(x = Date.Time, y = pCO2, color = "Station 51"), size = 1, alpha = 0.75) +
  geom_point(data = data_130_clean, aes(x = Date.Time, y = pCO2, color = "Station 130"), size = 1, alpha = 0.75) +
  labs(y = expression(paste("pCO"[2], " (", mu, "atm)")), x = "Date-Time", color = "Station") +
  scale_color_manual(values = c("Station 51" = "blue", "Station 130" = "red")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pco2_plot

# Create plot for DIC
dic_plot <- ggplot() +
  geom_point(data = data_51_clean, aes(x = Date.Time, y = DIC, color = "Station 51"), size = 1, alpha = 0.75) +
  geom_point(data = data_130_clean, aes(x = Date.Time, y = DIC, color = "Station 130"), size = 1, alpha = 0.75) +
  labs(y = expression(paste("DIC (", mu, "mol/kg)")), x = "Date-Time", color = "Station") +
  scale_color_manual(values = c("Station 51" = "blue", "Station 130" = "red")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dic_plot

# Function to create a dual y-axis plot with smoother
dual_axis_plot <- function(data, station_name, start_time, end_time, y_limits_pco2, y_limits_dic) {
  # Create the primary plot for pCO2
  p1 <- ggplot(data, aes(x = Date.Time, y = pCO2)) +
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
  p2 <- ggplot(data, aes(x = Date.Time, y = DIC)) +
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

# Set start and end times for the plots (10:00 to 10:00 the following day)
start_time_51 <- "2023-04-18 10:00:00"
end_time_51 <- "2023-04-19 10:00:00"

start_time_130 <- "2023-04-20 10:00:00"
end_time_130 <- "2023-04-21 10:00:00"

# Plot for Station 51
dual_axis_plot(data_51_clean, station_name = "51", start_time = start_time_51, end_time = end_time_51)

# Plot for Station 130
dual_axis_plot(data_130_clean, station_name = "130", start_time = start_time_130, end_time = end_time_130)

# Function to create a triple-panel plot for pCO2, DIC, and Oxygen -------------
triple_axis_plot <- function(data, station_name, start_time, end_time) {
  # Plot for pCO2
  p1 <- ggplot(data, aes(x = Date.Time, y = pCO2)) +
    geom_point(size = 1, alpha = 0.75) +
    geom_smooth(method = "loess", color = "grey", se = TRUE) +
    labs(y = expression(paste("pCO"[2], " (", mu, "atm)")), x = NULL) +
    scale_x_datetime(limits = c(as.POSIXct(start_time, tz = "UTC"),
                                as.POSIXct(end_time, tz = "UTC")),
                     date_breaks = "6 hours", date_labels = "%d-%b %H:%M") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  # Plot for DIC
  p2 <- ggplot(data, aes(x = Date.Time, y = DIC)) +
    geom_point(size = 1, alpha = 0.75) +
    geom_smooth(method = "loess", color = "grey", se = TRUE) +
    labs(y = expression(paste("DIC (", mu, "mol/kg)")), x = NULL) +
    scale_x_datetime(limits = c(as.POSIXct(start_time, tz = "UTC"),
                                as.POSIXct(end_time, tz = "UTC")),
                     date_breaks = "6 hours", date_labels = "%d-%b %H:%M") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  # Plot for Oxygen (O2uM)
  p3 <- ggplot(data, aes(x = Date.Time, y = O2uM)) +
    geom_point(size = 1, alpha = 0.75) +
    geom_smooth(method = "loess", color = "grey", se = TRUE) +
    labs(y = expression(paste("Oxygen (", mu, "mol/L)")), x = "Date-Time") +
    scale_x_datetime(limits = c(as.POSIXct(start_time, tz = "UTC"),
                                as.POSIXct(end_time, tz = "UTC")),
                     date_breaks = "6 hours", date_labels = "%d-%b %H:%M") +
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

# Function to create a single plot for a given parameter ------------------------
create_plot <- function(data, y_var, y_label, color, y_limits, show_x_axis = TRUE, smoother_span = 0.5) {
  p <- ggplot(data, aes(x = Date.Time, y = !!sym(y_var))) +
    geom_point(size = 1, alpha = 0.75, color = "black") +
    geom_smooth(method = "loess", color = color, se = TRUE, span = smoother_span) +
    labs(y = y_label, x = NULL) +
    scale_x_datetime(date_breaks = "6 hours", date_labels = "%d-%b %H:%M") +
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

# Create plots for Station 51 ---------------------------------------------------
pco2_51 <- create_plot(data_51_clean, "pCO2", expression(paste("pCO"[2], " (", mu, "atm)")), "grey", y_limits_pco2, show_x_axis = FALSE)
dic_51 <- create_plot(data_51_clean, "DIC", expression(paste("DIC (", mu, "mol/kg)")), "grey", y_limits_dic, show_x_axis = FALSE)
o2_51 <- create_plot(data_51_clean, "O2uM", expression(paste("Oxygen (", mu, "mol/L)")), "grey", y_limits_o2, show_x_axis = TRUE)

# Create plots for Station 130 --------------------------------------------------
pco2_130 <- create_plot(data_130_clean, "pCO2", expression(paste("pCO"[2], " (", mu, "atm)")), "grey", y_limits_pco2, show_x_axis = FALSE)
dic_130 <- create_plot(data_130_clean, "DIC", expression(paste("DIC (", mu, "mol/kg)")), "grey", y_limits_dic, show_x_axis = FALSE)
o2_130 <- create_plot(data_130_clean, "O2uM", expression(paste("Oxygen (", mu, "mol/L)")), "grey", y_limits_o2, show_x_axis = TRUE)

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
    dplyr::mutate(Date.Time = as.POSIXct(x, origin = "1970-01-01", tz = "UTC")) %>%  # Convert epoch to Date.Time
    dplyr::select(Date.Time, y, ymin, ymax, se)
  return(smoother_data_clean)
}

# Extract smoother values
dic_smoother_130 <- extract_smoother_values(dic_130)
o2_smoother_130 <- extract_smoother_values(o2_130)

# Process smoother data for Station 130
dic_smoother_130_clean <- process_smoother_data(dic_smoother_130)
o2_smoother_130_clean <- process_smoother_data(o2_smoother_130)

# Save smoother values to CSV
write.csv(dic_smoother_130_clean, "data/analysis/DIC_smoother_station_130.csv", row.names = FALSE)
write.csv(o2_smoother_130_clean, "data/analysis/O2_smoother_station_130.csv", row.names = FALSE)

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

# Linear model to describe the relation ---------------------------------
# Combine data from both stations
combined_data <- rbind(
  data.frame(Station = "51", Time = data_51_clean$Date.Time, DIC = data_51_clean$DIC),
  data.frame(Station = "130", Time = data_130_clean$Date.Time, DIC = data_130_clean$DIC)
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

cat(sprintf("For station 51, the maximum DIC value is %.2f at %s\n", max_DIC_51$DIC, max_DIC_51$Date.Time))
cat(sprintf("For station 51, the minimum DIC value is %.2f at %s\n", min_DIC_51$DIC, min_DIC_51$Date.Time))

# Find DIC maximum and minimum for Station 130
max_DIC_130 <- data_130_clean[which.max(data_130_clean$DIC), ]
min_DIC_130 <- data_130_clean[which.min(data_130_clean$DIC), ]

cat(sprintf("For station 130, the maximum DIC value is %.2f at %s\n", max_DIC_130$DIC, max_DIC_130$Date.Time))
cat(sprintf("For station 130, the minimum DIC value is %.2f at %s\n", min_DIC_130$DIC, min_DIC_130$Date.Time))

# Fit the model for O2
model <- lm(DIC ~ Station * Time, data = combined_data)

# Summary of the model
summary(model)

# ANOVA for the model
anova(model)

# Find O2 maximum and minimum for Station 51
max_O2_51 <- data_51_clean[which.max(data_51_clean$O2uM), ]
min_O2_51 <- data_51_clean[which.min(data_51_clean$O2uM), ]

cat(sprintf("For station 51, the maximum O2 value is %.2f μmol/L at %s\n", max_O2_51$O2uM, max_O2_51$Date.Time))
cat(sprintf("For station 51, the minimum O2 value is %.2f μmol/L at %s\n", min_O2_51$O2uM, min_O2_51$Date.Time))

# Find O2 maximum and minimum for Station 130
max_O2_130 <- data_130_clean[which.max(data_130_clean$O2uM), ]
min_O2_130 <- data_130_clean[which.min(data_130_clean$O2uM), ]

cat(sprintf("For station 130, the maximum O2 value is %.2f μmol/L at %s\n", max_O2_130$O2uM, max_O2_130$Date.Time))
cat(sprintf("For station 130, the minimum O2 value is %.2f μmol/L at %s\n", min_O2_130$O2uM, min_O2_130$Date.Time))

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

# Save the figure as an SVG file
output_path <- file.path(output_dir, "O2_DIC_correlation.svg")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "svg")
output_path <- file.path(output_dir, "O2_DIC_correlation.png")
ggsave(output_path, plot = correlation_plot, width = 18, height = 10, units = "cm", device = "png")
