library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(patchwork)
library(scales)
library(grid)

root <- normalizePath(Sys.getenv("FANLAB_RESEARCH2_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig2_milk_bioactive_remodeling")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

milk_path <- file.path(root, "direction1_publication_figures", "source_data", "all_milk_long_values.csv")
dyn_path <- file.path(root, "direction3_multiomics_blueprint", "tables", "direction3_top_milk_dynamic_features.csv")
hit_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_hmo_ltf_to_scfa_lagged_model_highlights.csv")

milk <- read_csv(milk_path, show_col_types = FALSE) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("D05", "D14", "D30", "D90")),
    Assay = factor(Assay, levels = c("LTF", "HMO", "LCFA")),
    FeatureClean = as.character(Feature),
    FeatureClean = if_else(FeatureClean == "Lactoferrin", "Lactoferrin", FeatureClean)
  ) %>%
  filter(
    !tolower(as.character(Include)) %in% c("false", "0", "no"),
    !is.na(Value),
    Timepoint %in% c("D05", "D14", "D30")
  )

dyn <- read_csv(dyn_path, show_col_types = FALSE) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("D05", "D14", "D30")),
    Assay = factor(Assay, levels = c("LTF", "HMO", "LCFA")),
    FeatureClean = as.character(FeatureClean)
  )

feature_lookup <- milk %>%
  distinct(Assay, FeatureClean)

join_diagnostic <- feature_lookup %>%
  full_join(dyn %>% distinct(Assay, FeatureClean, FDR), by = c("Assay", "FeatureClean"), suffix = c("_milk", "_dyn")) %>%
  mutate(join_status = case_when(
    !is.na(FDR) ~ "matched_dynamic_fdr",
    TRUE ~ "milk_without_dynamic_fdr"
  ))
write_csv(join_diagnostic, file.path(out_dir, "Fig2_join_diagnostic.csv"))

hits <- read_csv(hit_path, show_col_types = FALSE) %>%
  mutate(
    ExposureAssay = str_extract(ExposureLabel, "^[^:]+"),
    ExposureFeature = str_remove(ExposureLabel, "^[^:]+:"),
    ExposureFeature = str_replace_all(ExposureFeature, "LNFP-$", "LNFP-Ⅲ"),
    ExposureFeature = str_replace_all(ExposureFeature, "Lactoferrin", "Lactoferrin")
  )

candidate_counts <- hits %>%
  count(ExposureAssay, ExposureFeature, name = "LaggedSCFAHits") %>%
  rename(Assay = ExposureAssay, FeatureClean = ExposureFeature) %>%
  right_join(feature_lookup, by = c("Assay", "FeatureClean")) %>%
  filter(!is.na(LaggedSCFAHits))

feature_stats <- milk %>%
  group_by(Assay, FeatureClean, Timepoint) %>%
  summarise(
    n = n_distinct(DyadID),
    median = median(Value, na.rm = TRUE),
    q25 = quantile(Value, 0.25, na.rm = TRUE),
    q75 = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(dyn %>% distinct(Assay, FeatureClean, FDR), by = c("Assay", "FeatureClean")) %>%
  group_by(Assay, FeatureClean) %>%
  mutate(
    z_median = as.numeric(scale(log10(median + 1))),
    d05_median = median[Timepoint == "D05"][1],
    d30_median = median[Timepoint == "D30"][1],
    log2_d30_d05 = log2((d30_median + 1) / (d05_median + 1)),
    direction = case_when(
      log2_d30_d05 > 0.1 ~ "Higher at D30",
      log2_d30_d05 < -0.1 ~ "Lower at D30",
      TRUE ~ "Stable"
    )
  ) %>%
  ungroup() %>%
  left_join(candidate_counts, by = c("Assay", "FeatureClean")) %>%
  mutate(
    LaggedSCFAHits = replace_na(LaggedSCFAHits, 0L),
    Candidate = LaggedSCFAHits > 0,
    neglog10_fdr = -log10(FDR),
    sig_class = case_when(
      is.na(FDR) ~ "not ranked",
      FDR < 0.001 ~ "FDR < 0.001",
      FDR < 0.05 ~ "FDR < 0.05",
      TRUE ~ "FDR >= 0.05"
    )
  )

write_csv(feature_stats, file.path(out_dir, "Fig2_debug_feature_stats.csv"))

rank_df <- feature_stats %>%
  distinct(Assay, FeatureClean, FDR, log2_d30_d05, direction, LaggedSCFAHits, Candidate, neglog10_fdr) %>%
  filter(!is.na(FDR)) %>%
  arrange(FDR, desc(abs(log2_d30_d05))) %>%
  mutate(
    FeatureRank = factor(FeatureClean, levels = rev(unique(FeatureClean))),
    AssayLabel = recode(as.character(Assay), LTF = "Immune protein", HMO = "HMO", LCFA = "LCFA")
  )

atlas_order <- rank_df %>%
  arrange(Assay, desc(neglog10_fdr), desc(abs(log2_d30_d05))) %>%
  pull(FeatureClean) %>%
  unique()

atlas_df <- feature_stats %>%
  filter(FeatureClean %in% atlas_order) %>%
  mutate(
    FeatureAtlas = factor(FeatureClean, levels = rev(atlas_order)),
    AssayLabel = recode(as.character(Assay), LTF = "Immune protein", HMO = "HMO", LCFA = "LCFA"),
    TimeIndex = as.numeric(Timepoint)
  )

traj_features <- c("Lactoferrin", "6-SL", "3-SL", "LNnT", "3-FL", "FA(18:3)", "FA(18:1)")
traj_df <- milk %>%
  filter(FeatureClean %in% traj_features) %>%
  mutate(FeaturePlot = factor(FeatureClean, levels = traj_features)) %>%
  left_join(dyn %>% distinct(Assay, FeatureClean, FDR), by = c("Assay", "FeatureClean")) %>%
  group_by(FeaturePlot) %>%
  mutate(
    UnitScaled = as.numeric(scale(log10(Value + 1))),
    FDRLabel = paste0("FDR=", scientific(first(FDR), digits = 2))
  ) %>%
  ungroup()

module_score <- milk %>%
  group_by(Assay, FeatureClean) %>%
  mutate(z_value = as.numeric(scale(log10(Value + 1)))) %>%
  ungroup() %>%
  group_by(DyadID, Assay, Timepoint) %>%
  summarise(ModuleScore = mean(z_value, na.rm = TRUE), n_features = n_distinct(FeatureClean), .groups = "drop") %>%
  mutate(AssayLabel = recode(as.character(Assay), LTF = "Immune protein", HMO = "HMO", LCFA = "LCFA"))

module_summary <- module_score %>%
  group_by(AssayLabel, Timepoint) %>%
  summarise(
    n = n(),
    median = median(ModuleScore, na.rm = TRUE),
    q25 = quantile(ModuleScore, 0.25, na.rm = TRUE),
    q75 = quantile(ModuleScore, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

candidate_df <- rank_df %>%
  filter(Candidate) %>%
  mutate(
    FeatureCandidate = fct_reorder(FeatureClean, LaggedSCFAHits),
    CandidateLabel = paste0(FeatureClean, " (", LaggedSCFAHits, ")")
  )

sample_counts <- milk %>%
  group_by(Assay, Timepoint) %>%
  summarise(n = n_distinct(DyadID), .groups = "drop") %>%
  mutate(
    AssayLabel = recode(as.character(Assay), LTF = "Immune protein", HMO = "HMO", LCFA = "LCFA"),
    label = paste0("n=", n)
  )

pal_assay <- c("Immune protein" = "#D68422", "HMO" = "#168E87", "LCFA" = "#356AA0")
pal_dir <- c("Lower at D30" = "#1F77B4", "Stable" = "#8A8F98", "Higher at D30" = "#C9493D")
pal_heat <- c("#245C9A", "#F7F8FA", "#C9463E")

base_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#222222"),
    axis.ticks = element_line(linewidth = 0.3, colour = "#222222"),
    axis.text = element_text(colour = "#222222"),
    strip.background = element_rect(fill = "#F2F4F7", colour = NA),
    strip.text = element_text(face = "bold", colour = "#222222"),
    plot.title = element_text(face = "bold", size = 7.6, hjust = 0),
    plot.subtitle = element_blank(),
    legend.title = element_text(size = 6.3),
    legend.text = element_text(size = 5.8),
    legend.key.size = unit(2.6, "mm")
  )

pA <- ggplot(atlas_df, aes(TimeIndex, FeatureAtlas)) +
  geom_tile(aes(fill = z_median), width = 0.78, height = 0.66, colour = "white", linewidth = 0.28) +
  geom_point(
    data = atlas_df %>% distinct(FeatureAtlas, AssayLabel),
    aes(x = 0.55, y = FeatureAtlas, colour = AssayLabel),
    inherit.aes = FALSE,
    size = 1.6
  ) +
  geom_point(
    data = atlas_df %>% distinct(FeatureAtlas, TimeIndex, Candidate, LaggedSCFAHits) %>% filter(Candidate),
    aes(size = LaggedSCFAHits),
    shape = 21, fill = "white", colour = "#222222", stroke = 0.32
  ) +
  scale_x_continuous(limits = c(0.15, 3.45), breaks = c(1, 2, 3), labels = c("D05", "D14", "D30"), expand = c(0, 0)) +
  scale_fill_gradient2(low = pal_heat[1], mid = pal_heat[2], high = pal_heat[3], midpoint = 0, name = "Row z\nmedian") +
  scale_colour_manual(values = pal_assay, name = "Module") +
  scale_size_continuous(range = c(1.1, 2.6), breaks = c(1, 2, 4), name = "Lagged\nSCFA hits") +
  labs(
    title = "a  Dynamic milk bioactive atlas",
    x = NULL,
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 5.8),
    panel.grid = element_blank()
  ) +
  guides(
    colour = "none",
    fill = guide_colourbar(order = 2, barheight = unit(16, "mm"), barwidth = unit(2.5, "mm")),
    size = guide_legend(order = 3)
  )

pB <- ggplot(rank_df, aes(x = log2_d30_d05, y = FeatureRank)) +
  geom_vline(xintercept = 0, colour = "#C7CCD4", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2_d30_d05, yend = FeatureRank, colour = direction), linewidth = 0.65, alpha = 0.88) +
  geom_point(aes(size = neglog10_fdr, fill = AssayLabel), shape = 21, colour = "white", stroke = 0.25) +
  geom_point(data = rank_df %>% filter(Candidate), shape = 23, size = 2.2, fill = "white", colour = "#111111", stroke = 0.35) +
  scale_colour_manual(values = pal_dir, guide = "none") +
  scale_fill_manual(values = pal_assay, name = "Module") +
  scale_size_continuous(range = c(1.4, 4), name = "-log10 FDR") +
  labs(
    title = "b  Ranked temporal effect",
    x = "log2 median change (D30 / D05)",
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 5.7)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 2.2)),
    size = guide_legend(order = 2)
  )

pC <- ggplot(traj_df, aes(Timepoint, UnitScaled, group = DyadID)) +
  geom_line(colour = "#C7CCD4", linewidth = 0.22, alpha = 0.55) +
  geom_point(colour = "#C7CCD4", size = 0.45, alpha = 0.55) +
  stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 0.8, colour = "#102A43") +
  stat_summary(aes(group = 1), fun = median, geom = "point", size = 1.45, colour = "#102A43", fill = "white", shape = 21, stroke = 0.35) +
  facet_wrap(~ FeaturePlot, nrow = 1, scales = "free_y") +
  labs(
    title = "c  Representative trajectories",
    x = NULL,
    y = "Within-feature z score"
  ) +
  base_theme +
  theme(
    strip.text = element_text(size = 6.1),
    axis.text.x = element_text(angle = 0),
    panel.spacing.x = unit(1.4, "mm")
  )

pD <- ggplot(module_score, aes(Timepoint, ModuleScore, colour = AssayLabel, fill = AssayLabel)) +
  geom_hline(yintercept = 0, colour = "#D7DBE2", linewidth = 0.25) +
  geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.18, linewidth = 0.35, colour = "#4A5568") +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 0.7, alpha = 0.35, stroke = 0) +
  geom_line(
    data = module_summary,
    aes(Timepoint, median, group = AssayLabel, colour = AssayLabel),
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = module_summary,
    aes(Timepoint, median, colour = AssayLabel),
    size = 1.7,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ AssayLabel, nrow = 1) +
  scale_colour_manual(values = pal_assay, guide = "none") +
  scale_fill_manual(values = pal_assay, guide = "none") +
  labs(
    title = "d  Module-level temporal drift",
    x = NULL,
    y = "Module score"
  ) +
  base_theme +
  theme(panel.spacing.x = unit(2, "mm"))

pE <- ggplot(candidate_df, aes(x = LaggedSCFAHits, y = FeatureCandidate, fill = AssayLabel)) +
  geom_col(width = 0.56, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = paste0("n=", LaggedSCFAHits)), hjust = -0.12, size = 2, colour = "#263238") +
  scale_fill_manual(values = pal_assay, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)), breaks = pretty_breaks(4)) +
  labs(
    title = "e  Prioritized exposure candidates",
    x = "Number of lagged SCFA hits",
    y = NULL
  ) +
  base_theme +
  theme(axis.text.y = element_text(size = 6.1))

title_plot <- ggplot() +
  annotate("text", x = 0, y = 0.78, label = "Fig. 2 | Milk bioactive remodeling", hjust = 0, vjust = 1, size = 3.25, fontface = "bold", family = "Arial", colour = "#102A43") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void()

count_text <- sample_counts %>%
  arrange(AssayLabel, Timepoint) %>%
  mutate(part = paste0(AssayLabel, " ", Timepoint, " ", label)) %>%
  group_by(AssayLabel) %>%
  summarise(line = paste(part, collapse = "  |  "), .groups = "drop") %>%
  pull(line) %>%
  paste(collapse = "\n")

note_plot <- ggplot() +
  annotate("text", x = 0, y = 0.72, label = count_text, hjust = 0, vjust = 1, size = 1.45, family = "Arial", colour = "#4A5568", lineheight = 1.08) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void()

fig <- title_plot / ((pA | pB) / pC / (pD | pE) / note_plot) +
  plot_layout(heights = c(0.045, 1), widths = c(1)) &
  theme(plot.margin = margin(3, 4, 3, 4))

base_file <- file.path(out_dir, "Fig2_milk_bioactive_remodeling_redrawn")
ggsave(paste0(base_file, ".svg"), fig, width = 183, height = 210, units = "mm", device = svglite::svglite)
ggsave(paste0(base_file, ".pdf"), fig, width = 183, height = 210, units = "mm", device = cairo_pdf)
ragg::agg_png(paste0(base_file, ".png"), width = 183, height = 210, units = "mm", res = 450, background = "white")
print(fig)
dev.off()
ragg::agg_tiff(paste0(base_file, ".tiff"), width = 183, height = 210, units = "mm", res = 600, compression = "lzw", background = "white")
print(fig)
dev.off()

source_data <- list(
  atlas = atlas_df,
  ranked_features = rank_df,
  trajectories = traj_df,
  module_scores = module_score,
  candidates = candidate_df,
  sample_counts = sample_counts
)
write_csv(atlas_df, file.path(out_dir, "Fig2_source_atlas.csv"))
write_csv(rank_df, file.path(out_dir, "Fig2_source_ranked_dynamic_features.csv"))
write_csv(module_score, file.path(out_dir, "Fig2_source_module_scores.csv"))
write_csv(candidate_df, file.path(out_dir, "Fig2_source_candidate_exposures.csv"))
writeLines(c(
  "# Fig. 2 redraw notes",
  "",
  "Conclusion: early lactation shows coordinated and module-specific remodeling of milk bioactive features.",
  "Evidence chain: atlas of dynamic features; ranked D30/D05 changes; representative dyad trajectories; module-level drift; candidate exposure markers linking this figure to lagged SCFA analyses.",
  "Boundary: candidate exposure panel is prioritization for downstream association models, not causal evidence.",
  paste0("Rendered: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Source data: ", milk_path)
), file.path(out_dir, "Fig2_redraw_notes.md"))

message("Saved Fig. 2 redraw to: ", out_dir)
