############################################################
# Diel expression analysis (GAM) + WGCNA
# Phaeocystis genome-resolved data
#
# Input: long-format CSV with columns:
#   gene_id | sample | Station | Date | time_of_day_hours | TPM | day_moment
#
# Output:
#   - GAM stats & significant diel gene list (per station)
#   - Summary: fraction of cyclic genes per station
#   - Optional day/night test results
#   - WGCNA modules for rhythmic genes + eigengene GAM tests
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(mgcv); library(pbapply); library(WGCNA)
  library(ggplot2); library(lubridate); library(flashClust);
  library(ggforce); library(patchwork)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
enableWGCNAThreads()

#------------------------------------------------------------
# 1. GAM for single gene
#------------------------------------------------------------
fit_gam_gene <- function(dfg, k = 8) {
  dfg <- as.data.frame(dfg)
  dfg$TPM <- suppressWarnings(as.numeric(dfg$TPM))
  dfg$td_norm <- suppressWarnings(as.numeric(dfg$time_of_day_hours) / 24)

  if (any(!is.finite(dfg$TPM)) || any(!is.finite(dfg$td_norm))) {
    return(list(p_value=NA_real_, td_max=NA_real_, amplitude=NA_real_,
                aic0=NA_real_, aic1=NA_real_, deltaAIC=NA_real_, devexp=NA_real_))
  }
  if (nrow(dfg) < 4 || length(unique(dfg$TPM)) < 2 || length(unique(dfg$td_norm)) < 3) {
    return(list(p_value=NA_real_, td_max=NA_real_, amplitude=NA_real_,
                aic0=NA_real_, aic1=NA_real_, deltaAIC=NA_real_, devexp=NA_real_))
  }

  res <- tryCatch({
    m0 <- gam(TPM ~ 1, data = dfg, method = "REML")
    m1 <- gam(TPM ~ s(td_norm, bs = "cc", k = k), data = dfg, method = "REML")
    p  <- anova(m0, m1, test = "F")$`Pr(>F)`[2]

    aic0 <- AIC(m0); aic1 <- AIC(m1); dAIC <- aic0 - aic1
    dev0 <- deviance(m0); dev1 <- deviance(m1)
    devexp <- if (is.finite(dev0) && dev0 > 0) max(0, (dev0 - dev1)/dev0) else NA_real_

    grid <- data.frame(td_norm = seq(0, 1, length.out = 200))
    pred <- predict(m1, newdata = grid, type = "response")
    td_max <- grid$td_norm[which.max(pred)] * 24
    amplitude <- max(pred, na.rm = TRUE) - min(pred, na.rm = TRUE)

    list(p_value=p, td_max=td_max, amplitude=amplitude,
         aic0=aic0, aic1=aic1, deltaAIC=dAIC, devexp=devexp)
  }, error = function(e) {
    list(p_value=NA_real_, td_max=NA_real_, amplitude=NA_real_,
         aic0=NA_real_, aic1=NA_real_, deltaAIC=NA_real_, devexp=NA_real_)
  })
  res
}

#------------------------------------------------------------
# 2. Run GAM per station
#------------------------------------------------------------
run_station_gam <- function(expr_station, station_label, outdir="figures/metatranscriptomics/genome_resolved") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
  message("=== Running GAM at Station ", station_label, " ===")

  expr_filt <- expr_station %>%
    group_by(gene_id) %>%
    filter(max(TPM, na.rm=TRUE) > 1,
           var(TPM, na.rm=TRUE) > 0.5) %>%
    ungroup()
  message("Kept ", n_distinct(expr_filt$gene_id), " genes after filtering")

  by_gene <- split(expr_filt, expr_filt$gene_id)
  gam_list <- pbapply::pblapply(by_gene, fit_gam_gene, k=8)

  gam_df <- bind_rows(lapply(names(gam_list), function(g) {
    as.data.frame(cbind(gene_id=g, t(unlist(gam_list[[g]]))))
  }))
  num_cols <- c("p_value","td_max","amplitude","aic0","aic1","deltaAIC","devexp")
  gam_df[num_cols] <- lapply(gam_df[num_cols], as.numeric)
  gam_df$padj <- p.adjust(gam_df$p_value, method="BH")

  # Significant genes
  sig <- gam_df %>%
    filter(!is.na(padj),
           padj < 0.05,
           deltaAIC > 2,
           devexp >= 0.10,
           amplitude > 1)

  write_csv(gam_df, file.path(outdir, paste0("station", station_label, "_gam_all.csv")))
  write_csv(sig, file.path(outdir, paste0("station", station_label, "_gam_sig.csv")))

  list(gam=gam_df, sig=sig)
}

#Summarize cyclicity per station
summarize_cyclicity <- function(gam_df, sig_df, station_label) {
  total <- nrow(gam_df)
  sig   <- nrow(sig_df)
  prop  <- round(100 * sig / total, 1)
  message("Station ", station_label, ": ", sig, "/", total, " cyclic (", prop, "%)")
  data.frame(station=station_label, total_tested=total,
             cyclic_genes=sig, proportion_cyclic=prop)
}

#------------------------------------------------------------
# 3. Binary day/night test
#------------------------------------------------------------
fit_daynight_gene <- function(dfg) {
  if (!"day_moment" %in% colnames(dfg)) return(NULL)
  # Combine "Civil twilight", "Nautical twilight", and "Astronomical twilight" into "Twilight"
  dfg$day_moment <- ifelse(dfg$day_moment %in% c("Civil twilight", "Nautical twilight", "Astronomical twilight"),
                           "Twilight", dfg$day_moment)
  dfg <- dfg %>% mutate(day_moment=factor(day_moment, levels=c("Night", "Twilight", "Day")))
  tryCatch({
    m <- lm(TPM ~ day_moment, data=dfg)
    p <- anova(m)$`Pr(>F)`[1]
    logFC <- diff(tapply(dfg$TPM, dfg$day_moment, mean))
    list(p_value=p, logFC=logFC)
  }, error=function(e) NULL)
}

#------------------------------------------------------------
# 4. WGCNA on cyclic genes
#------------------------------------------------------------
run_wgcna_on_sig <- function(expr_station, sig_table,
                             outdir="figures/metatranscriptomics/genome_resolved/WGCNA/",
                             datadir="data/analysis/WGCNA_130/") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
  if (!dir.exists(datadir)) dir.create(datadir, recursive=TRUE)

  keep_genes <- unique(sig_table$gene_id)
  message("Running WGCNA on ", length(keep_genes), " rhythmic genes")

  mat <- expr_station %>%
    filter(gene_id %in% keep_genes) %>%
    select(gene_id, sample, TPM) %>%
    pivot_wider(names_from=gene_id, values_from=TPM, values_fill=0)

  sample_ids <- mat$sample
  X <- as.matrix(mat[,-1, drop=FALSE])
  datExpr <- log2(X+1); rownames(datExpr) <- sample_ids

  gsg <- goodSamplesGenes(datExpr, verbose=3)
  if (!gsg$allOK) datExpr <- datExpr[, gsg$goodGenes, drop=FALSE]

  powers <- c(1:20, seq(20,60,5))
  disableWGCNAThreads()
  sft <- pickSoftThreshold(datExpr, powerVector=powers, networkType="signed", blockSize=500, verbose=2)
  enableWGCNAThreads()
  softPower <- min(ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate), 20)
  message("Chosen soft threshold: ", softPower)

  adjacencyM <- adjacency(datExpr, power=softPower, type="signed")
  TOM <- TOMsimilarity(adjacencyM, TOMType="signed")
  dissTOM <- 1 - TOM

  geneTree <- flashClust(as.dist(dissTOM), method="average")
  dynamicMods <- cutreeDynamic(dendro=geneTree, distM=dissTOM,
                               deepSplit=2, pamRespectsDendro=FALSE,
                               minClusterSize=70)
  dynamicColors <- labels2colors(dynamicMods)

  svg(file.path(outdir,"gene_dendrogram_and_module_colors.svg"))
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
  dev.off()

  merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight=0.5, verbose=3)
  mergedColors <- merge$colors
  mergedMEs <- merge$newMEs

  # Map each gene to its module
  gene_module <- data.frame(
    gene_id = colnames(datExpr),
    module  = mergedColors
  )

  # Correlate each gene with each module eigengene
  gene_module_corr <- cor(datExpr, mergedMEs, use="p")

  # Get p-values for the correlations
  gene_module_p <- corPvalueStudent(gene_module_corr, nSamples = nrow(datExpr))

  # Combine correlations and p-values into a long table
  ## First, map each gene to its assigned module color
  gene_module <- tibble::tibble(
    gene_id = colnames(datExpr),
    assigned_color = mergedColors
  )

  # Get gene–ME correlations and p-values
  gene_module_corr <- cor(datExpr, mergedMEs, use = "p")
  gene_module_p    <- corPvalueStudent(gene_module_corr, nSamples = nrow(datExpr))

  # Reshape to long format
  corr_long <- as.data.frame(gene_module_corr) %>%
    tibble::rownames_to_column("gene_id") %>%
    tidyr::pivot_longer(-gene_id, names_to = "ME", values_to = "correlation") %>%
    dplyr::mutate(module = gsub("^ME", "", ME), .keep = "unused")

  p_long <- as.data.frame(gene_module_p) %>%
    tibble::rownames_to_column("gene_id") %>%
    tidyr::pivot_longer(-gene_id, names_to = "ME", values_to = "p_value") %>%
    dplyr::mutate(module = gsub("^ME", "", ME), .keep = "unused")

  # Merge r and p
  gene_me_long <- corr_long %>%
    dplyr::left_join(p_long, by = c("gene_id", "module"))

  # Per-gene best module by correlation
  best_by_gene <- gene_me_long %>%
    dplyr::group_by(gene_id) %>%
    dplyr::slice_max(order_by = abs(correlation), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(gene_id,
                     best_module = module,
                     best_corr   = correlation)

  # Final table: all r/p + assignment + best call + consistency flag
  gene_module_corr_df <- gene_me_long %>%
    dplyr::left_join(best_by_gene, by = "gene_id") %>%
    dplyr::left_join(gene_module,   by = "gene_id") %>%
    dplyr::mutate(consistent = (assigned_color == best_module))

  # Save
  write.csv(mergedMEs, file = file.path(datadir, "module_eigengenes.csv"), row.names = TRUE)
  write.csv(gene_module_corr_df, file = file.path(datadir, "gene_module_correlations.csv"), row.names = FALSE)

  # Write 1 file per module with the “consistent” members ranked by |kME|
  dir.create(file.path(datadir, "modules"), showWarnings = FALSE)
  gene_module_corr_df %>%
    dplyr::filter(module == assigned_color) %>%
    dplyr::group_by(module) %>%
    dplyr::arrange(dplyr::desc(abs(correlation)), .by_group = TRUE) %>%
    dplyr::group_walk(~ readr::write_csv(.x, file.path(datadir, "modules", paste0(.y$module, "_content.csv"))))

  write.csv(mergedMEs, file=file.path(datadir, "module_eigengenes.csv"), row.names=TRUE)
  write.csv(gene_module_corr_df, file=file.path(datadir, "gene_module_correlations.csv"), row.names=FALSE)
  list(colors=mergedColors, MEs=mergedMEs, gene_module_corr=gene_module_corr_df)
}

# ===========================
# =========== RUN ===========
# ===========================
expr <- read_csv("data/phaglo1_mapping/gene_expression/phaeo_long_tpm_tpl.csv", show_col_types=FALSE) %>%
  mutate(Date=ymd_hms(Date),
         time_of_day_hours=as.numeric(time_of_day_hours),
         TPM=as.numeric(TPM),
         TPL=as.numeric(TPL))

# Station 130 (bloom)
expr_130 <- expr %>% filter(Station=="130")
gam130 <- run_station_gam(expr_130, "130")
summary130 <- summarize_cyclicity(gam130$gam, gam130$sig, "130")

# Station 51 (low abundance)
expr_51 <- expr %>% filter(Station=="51")
gam51 <- run_station_gam(expr_51, "51")
summary51 <- summarize_cyclicity(gam51$gam, gam51$sig, "51")

# Compare
cyclicity_summary <- bind_rows(summary130, summary51)
print(cyclicity_summary)

# Day/night test for station 130
daynight_130 <- expr_130 %>%
  dplyr::filter(gene_id %in% unique(gam130$gam$gene_id)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::group_modify(~ {
    out <- fit_daynight_gene(.x)
    if (is.null(out)) tibble::tibble(p_value = NA_real_, logFC = NA_real_)
    else tibble::as_tibble(out)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(padj = p.adjust(p_value, method = "BH"))

# WGCNA on bloom station rhythmic genes
wgcna130 <- run_wgcna_on_sig(expr_130, gam130$sig,
                             outdir="figures/metatranscriptomics/genome_resolved/WGCNA/130",
                             datadir="data/analysis/WGCNA_130/")

# Analyze WGCNA modules
# Sizes
module_sizes <- tibble::tibble(module = unique(wgcna130$colors)) %>%
  dplyr::left_join(
    tibble::tibble(module = wgcna130$colors) %>%
      dplyr::count(module, name = "n_genes"),
    by = "module"
  )

# Read environmental data
env <- read.csv("data/samples_env.csv", stringsAsFactors = FALSE)

# Make sure the sample IDs match the expression data
env <- env %>%
  filter(Station %in% rownames(wgcna130$MEs)) %>%
  arrange(match(Station, rownames(wgcna130$MEs))) %>%
  mutate(Date = ymd_hms(Date))

# Build long table and label
MElong <- wgcna130$MEs %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::left_join(expr_130 %>% dplyr::select(sample, Date, time_of_day_hours) %>% dplyr::distinct(), by = "sample") %>%
  tidyr::pivot_longer(-c(sample, Date, time_of_day_hours), names_to = "ME", values_to = "eigengene") %>%
  dplyr::mutate(module = gsub("^ME", "", ME)) %>%
  dplyr::left_join(module_sizes, by = "module") %>%
  dplyr::mutate(facet_lab = paste0(module, " (n=", n_genes, ")"),
                module = factor(module, levels = unique(wgcna130$colors))) %>%
  dplyr::left_join(env)

module_colors <- setNames(as.character(unique(wgcna130$colors)),
                          unique(wgcna130$colors))
ggplot(MElong, aes(x = Date, y = eigengene, color = module)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = 8), se = FALSE) +
  facet_wrap(~ facet_lab, scales = "free_y") +
  scale_color_manual(values = module_colors) +
  theme_minimal() +
  labs(x = "Time", y = "Module eigengene")

# Plot module expression with diel background shading
# First define light phases
light_colors <- c("Night" = "#d9d9d9",
                  "Astronomical twilight" = "#ffb347",
                  "Nautical twilight" = "#ffc870",
                  "Civil twilight" = "#ffe0a3",
                  "Day" = "#ffffb3")

# Get unique Date + day_moment entries sorted by time
bg_df <- MElong %>%
  dplyr::select(Date, day_moment) %>%
  dplyr::distinct() %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    Date_start = Date,
    Date_end = lead(Date)
  ) %>%
  dplyr::filter(!is.na(Date_end))

start_time <- as.POSIXct("2023-04-20 07:00:00", tz = "UTC")
end_time   <- as.POSIXct("2023-04-21 09:00:00", tz = "UTC")

p_modules <- ggplot(MElong, aes(x = Date, y = eigengene, color = module)) +
  geom_rect(data = bg_df,
            aes(xmin = Date_start, xmax = Date_end, ymin = -Inf, ymax = Inf, fill = day_moment),
            inherit.aes = FALSE,
            alpha = 0.3,
            color = NA) +
  geom_point(size = 1) +
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cc", k = 8),
              se = FALSE,
              linewidth = 0.9) +
  facet_wrap(~ facet_lab, scales = "free_y", ncol = 2) +
  scale_color_manual(values = module_colors, guide = "none") +
  scale_fill_manual(values = light_colors, name = "Diel phase") +
  scale_x_datetime(
    limits = c(start_time, end_time),
    date_breaks = "2 hours",
    date_labels = "%H:%M",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_minimal(base_size = 11) +
  labs(x = "Time of day", y = "Module eigengene expression") +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  )

p_modules

ggsave("figures/metatranscriptomics/genome_resolved/WGCNA/130/module_eigengene_expression.pdf", p_modules, width = 18, height = 14, units = "cm")

# Keep only numeric parameters
env_params <- env %>%
  # Station to rownames
  tibble::column_to_rownames("Station") %>%
  select(-Date)  %>%
  # Select relevant parameters
  select("Temperature","Salinity","Oxygen","Fluorescence","NH4","NO2","NO3",
         "NOX","PO4","Si","TEP","sea_surface_height_above_sea_level") %>%
  # Z-score normalization
  mutate(across(everything(), ~ (.-mean(., na.rm=TRUE))/sd(., na.rm=TRUE)))

moduleTraitCor <- cor(wgcna130$MEs, env_params[rownames(wgcna130$MEs), ], use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(wgcna130$MEs))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(env_params),
               yLabels = colnames(wgcna130$MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = paste(signif(moduleTraitCor,2), "\n(", signif(moduleTraitPvalue,1), ")", sep=""))

# Plot all module eigengenes (except the grey bin) with O2' dynamics
O2prime <- read_csv("data/analysis/O2prime_resids.csv")
head(O2prime)

O2prime_station130 <- O2prime %>%
  filter(Station == "130") %>%
  arrange(Date) %>%
  mutate(time_of_day_hours = as.numeric(format(Date, "%H")) +
           as.numeric(format(Date, "%M")) / 60 +
           as.numeric(format(Date, "%S")) / 3600)

# Convert to data.table
melong_dt <- as.data.table(MElong %>% filter(module != "grey"))
o2_dt <- as.data.table(O2prime_station130 %>% select(Date, O2prime, O2_pred, O2_resid))

# Set keys
setkey(melong_dt, Date)
setkey(o2_dt, Date)

# Perform rolling join to match nearest O2prime value
MElong_O2 <- o2_dt[melong_dt, roll = "nearest"]

# Check
head(MElong_O2)

# Plot the O2' dynamics over time
p1 <- ggplot(O2prime_station130, aes(x = Date, y = O2prime)) +
  geom_rect(data = bg_df,
            aes(xmin = Date_start, xmax = Date_end, ymin = -Inf, ymax = Inf, fill = day_moment),
            inherit.aes = FALSE, alpha = 0.25, color = NA) +
  # Add bar plot with O2_resid values over time
  geom_bar(aes(y = O2_resid), stat = "identity", fill = "#CEDCA0", alpha = 0.9) +
  geom_point(color = "darkgreen", size = 0.5) +
  # Add smoother line based on O2_pred
  geom_smooth(aes(y = O2_pred), method = "gam",
              formula = y ~ s(x, bs = "cc", k = 8),
              color = "darkgreen", size = 0.9, se = FALSE) +
  # Add horizontal line at y=0
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40", size = 0.6) +
  scale_fill_manual(values = light_colors) +
  scale_x_datetime(
    limits = c(start_time, end_time),
    date_breaks = "2 hours",
    date_labels = "%H:%M",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_minimal(base_size = 11) +
  labs(x = "Time of day", y = "O2'") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10)
  )
p1

ggsave("figures/environmental/O2prime_resid_dynamics.pdf", p1, width = 9, height = 7, units = "cm")

# Plot O2' residuals vs. module eigengenes
ggplot(MElong_O2, aes(x = O2_resid, y = eigengene)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ module, scales = "free_y") +
  labs(x = "O₂′", y = "Module eigengene expression") +
  theme_minimal()

# WGCNA on station 51 rhythmic genes
wgcna51 <- run_wgcna_on_sig(expr_51, gam51$sig,
                            outdir = "figures/metatranscriptomics/genome_resolved/WGCNA/51",
                            datadir = "data/analysis/WGCNA_51/")
# Only 1 module detected (turquoise), so no further analysis performed
