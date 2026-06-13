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
out_dir <- file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig4_scfa_maturation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scfa_path <- file.path(root, "clinical_scfa_metagenome_report", "inputs", "bf_scfa_clean_ascii.csv")
score_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_scfa_maturation_score.csv")
trend_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_scfa_time_trend_for_maturation_score.csv")
model_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_scfa_maturation_score_model_summary.csv")

time_levels <- c("D05", "D14", "D30", "D90")
pal_time <- c(D05 = "#5778A4", D14 = "#79A88D", D30 = "#D89C3D", D90 = "#B56576")
pal_effect <- c("FDR < 0.05" = "#1F6FAE", "P < 0.05" = "#6A9BC8", "NS" = "#9B9B9B")
pal_scfa <- c(
  Succinate = "#245F9E",
  Isobutyrate = "#2F9C95",
  Butyrate = "#D99019",
  Methylbutyrate2 = "#C03D4D",
  Valerate = "#8058C7",
  Isovalerate = "#7A7A7A",
  Acetate = "#8BAF3D",
  Propionate = "#B97952",
  Hexanoate = "#5A6F82"
)
scfa_order <- c("Succinate", "Isobutyrate", "Butyrate", "Methylbutyrate2", "Valerate", "Isovalerate", "Acetate", "Hexanoate", "Propionate")
scfa_labels <- c(
  Succinate = "Succinate",
  Isobutyrate = "Isobutyrate",
  Butyrate = "Butyrate",
  Methylbutyrate2 = "2-Methylbutyrate",
  Valerate = "Valerate",
  Isovalerate = "Isovalerate",
  Acetate = "Acetate",
  Hexanoate = "Hexanoate",
  Propionate = "Propionate"
)

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

scfa_raw <- read_csv(scfa_path, show_col_types = FALSE) %>%
  mutate(
    DyadID = as.character(DyadID),
    Timepoint = factor(TimepointCode, levels = time_levels),
    Day = recode(as.character(TimepointCode), D05 = 5, D14 = 14, D30 = 30, D90 = 90) %>% as.numeric()
  )

score <- read_csv(score_path, show_col_types = FALSE) %>%
  mutate(
    DyadID = as.character(DyadID),
    Timepoint = factor(TimepointCode, levels = time_levels),
    Day = as.numeric(Day)
  )

trend <- read_csv(trend_path, show_col_types = FALSE) %>%
  mutate(
    SCFA = factor(SCFA, levels = scfa_order),
    SCFALabel = recode(as.character(SCFA), !!!scfa_labels),
    Evidence = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      PValue < 0.05 ~ "P < 0.05",
      TRUE ~ "NS"
    ),
    Evidence = factor(Evidence, levels = c("FDR < 0.05", "P < 0.05", "NS"))
  )

model <- read_csv(model_path, show_col_types = FALSE)
model_label <- paste0(
  "N=", model$N[1],
  ", dyads=", model$Dyads[1],
  ", slope/day=", sprintf("%.4f", model$SlopePerDay[1]),
  ", P=", fmt_p(model$PValue[1]),
  ", R2=", sprintf("%.2f", model$R2[1])
)

scfa_long <- scfa_raw %>%
  select(Sample, DyadID, Timepoint, Day, all_of(scfa_order)) %>%
  pivot_longer(cols = all_of(scfa_order), names_to = "SCFA", values_to = "Value") %>%
  mutate(
    SCFA = factor(SCFA, levels = scfa_order),
    SCFALabel = recode(as.character(SCFA), !!!scfa_labels),
    LogValue = log10(Value + 1)
  ) %>%
  left_join(trend %>% select(SCFA, Slope, PValue, FDR, Evidence), by = "SCFA")

sample_counts <- score %>%
  distinct(Sample, DyadID, Timepoint) %>%
  count(Timepoint, name = "n_samples") %>%
  mutate(label = paste0(as.character(Timepoint), "\nn=", n_samples))

sample_order <- score %>%
  group_by(Timepoint) %>%
  arrange(SCFA_MaturationScore, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(SamplePlot = factor(Sample, levels = Sample))

heat_df <- scfa_long %>%
  left_join(sample_order %>% select(Sample, SamplePlot, SCFA_MaturationScore), by = "Sample") %>%
  group_by(SCFA) %>%
  mutate(Z = as.numeric(scale(LogValue))) %>%
  ungroup() %>%
  mutate(
    SCFALabel = factor(SCFALabel, levels = rev(scfa_labels[scfa_order])),
    Timepoint = factor(Timepoint, levels = time_levels)
  )

time_strip_df <- sample_order %>%
  mutate(y = 1)

score_summary <- score %>%
  group_by(Timepoint) %>%
  summarise(
    n = n_distinct(Sample),
    median = median(SCFA_MaturationScore, na.rm = TRUE),
    q25 = quantile(SCFA_MaturationScore, 0.25, na.rm = TRUE),
    q75 = quantile(SCFA_MaturationScore, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(TimeIndex = as.numeric(Timepoint))

traj_keep <- c("Succinate", "Isobutyrate", "Butyrate", "Methylbutyrate2", "Valerate", "Acetate")
traj_df <- scfa_long %>%
  filter(SCFA %in% traj_keep) %>%
  group_by(Timepoint, SCFA, SCFALabel) %>%
  summarise(
    n = n_distinct(Sample),
    median = median(LogValue, na.rm = TRUE),
    q25 = quantile(LogValue, 0.25, na.rm = TRUE),
    q75 = quantile(LogValue, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    TimeIndex = as.numeric(Timepoint),
    SCFALabel = factor(SCFALabel, levels = scfa_labels[traj_keep])
  )

source_data <- list(
  sample_counts = sample_counts,
  scfa_sample_landscape = heat_df,
  score_values = score,
  score_summary = score_summary,
  score_model = model,
  scfa_time_trends = trend,
  scfa_trajectories = traj_df
)

for (nm in names(source_data)) {
  write_csv(source_data[[nm]], file.path(out_dir, paste0("Fig4_source_", nm, ".csv")))
}

theme_set(
  theme_classic(base_size = 6.2, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#1F1F1F"),
      axis.ticks = element_line(linewidth = 0.28, colour = "#1F1F1F"),
      axis.text = element_text(colour = "#1F1F1F"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "#1F1F1F"),
      plot.title = element_text(face = "bold", size = 7.2, hjust = 0, colour = "#111111"),
      legend.title = element_text(size = 6.0),
      legend.text = element_text(size = 5.5),
      legend.key.size = unit(2.5, "mm"),
      plot.margin = margin(3, 4, 3, 4)
    )
)

pA_strip <- ggplot(time_strip_df, aes(SamplePlot, y, fill = Timepoint)) +
  geom_tile(width = 0.95, height = 1) +
  facet_grid(. ~ Timepoint, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_time, guide = "none") +
  labs(title = "a  SCFA sample landscape", x = NULL, y = NULL) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 7.2, hjust = 0, margin = margin(0, 0, 2, 0)),
    strip.text = element_blank(),
    panel.spacing.x = unit(0.8, "mm"),
    plot.margin = margin(1, 4, 0, 4)
  )

pA_heat <- ggplot(heat_df, aes(SamplePlot, SCFALabel, fill = Z)) +
  geom_tile(width = 0.95, height = 0.90) +
  facet_grid(. ~ Timepoint, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#356C9E",
    mid = "#F7F7F7",
    high = "#B44A42",
    midpoint = 0,
    limits = c(-2.4, 2.4),
    oob = squish,
    name = "Row z\nlog10"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 6.0) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 5.6, colour = "#222222"),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.8, "mm"),
    strip.text = element_text(face = "bold", size = 5.8),
    legend.position = "right",
    legend.key.height = unit(13, "mm"),
    plot.margin = margin(0, 4, 2, 4)
  )

pA <- pA_strip / pA_heat + plot_layout(heights = c(0.10, 1))

pB <- ggplot(score, aes(as.numeric(Timepoint), SCFA_MaturationScore)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#9A9A9A") +
  geom_line(aes(group = DyadID), colour = "#CFCFCF", alpha = 0.34, linewidth = 0.28) +
  geom_jitter(aes(colour = Timepoint), width = 0.06, height = 0, size = 1.35, alpha = 0.78, show.legend = FALSE) +
  geom_ribbon(data = score_summary, aes(x = TimeIndex, ymin = q25, ymax = q75), inherit.aes = FALSE, fill = "#2A6FBB", alpha = 0.12) +
  geom_line(data = score_summary, aes(TimeIndex, median), inherit.aes = FALSE, linewidth = 0.7, colour = "#245F9E") +
  geom_point(data = score_summary, aes(TimeIndex, median), inherit.aes = FALSE, size = 2.0, colour = "#245F9E") +
  annotate("text", x = 1.05, y = Inf, label = model_label, hjust = 0, vjust = 1.25, size = 1.75, family = "Arial", colour = "#4A4A4A") +
  geom_text(
    data = sample_counts %>% mutate(TimeIndex = as.numeric(Timepoint), y = min(score$SCFA_MaturationScore, na.rm = TRUE) - 0.13),
    aes(TimeIndex, y, label = paste0("n=", n_samples)),
    inherit.aes = FALSE,
    size = 1.65,
    family = "Arial",
    colour = "#4A4A4A"
  ) +
  scale_colour_manual(values = pal_time) +
  scale_x_continuous(breaks = 1:4, labels = time_levels, expand = expansion(mult = c(0.04, 0.05))) +
  labs(
    title = "b  Composite SCFA maturation score",
    x = NULL,
    y = "Maturation score"
  ) +
  theme(panel.grid.major.y = element_line(linewidth = 0.15, colour = "#E9E9E9"))

pC <- ggplot(trend, aes(Slope, fct_reorder(SCFALabel, Slope))) +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "#777777") +
  geom_segment(aes(x = 0, xend = Slope, y = fct_reorder(SCFALabel, Slope), yend = fct_reorder(SCFALabel, Slope), colour = Evidence), linewidth = 0.55) +
  geom_point(aes(colour = Evidence, shape = Evidence), size = 2.1, stroke = 0.5) +
  scale_colour_manual(values = pal_effect, name = NULL) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "P < 0.05" = 17, "NS" = 1), name = NULL) +
  labs(
    title = "c  Feature-level age effects",
    x = "Linear age slope",
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 5.5),
    panel.grid.major.y = element_line(linewidth = 0.12, colour = "#EFEFEF"),
    legend.position = "right"
  )

pD <- ggplot(traj_df, aes(TimeIndex, median, colour = SCFA, group = SCFA)) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = SCFA), colour = NA, alpha = 0.13, show.legend = FALSE) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.6) +
  facet_wrap(~ SCFALabel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal_scfa, guide = "none") +
  scale_fill_manual(values = pal_scfa, guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = time_levels, expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "d  Representative SCFA trajectories",
    x = NULL,
    y = "Median log10(value + 1)"
  ) +
  theme(
    axis.text.x = element_text(size = 5.1),
    strip.text = element_text(size = 5.3),
    panel.grid.major.y = element_line(linewidth = 0.15, colour = "#E9E9E9")
  )

pE <- ggplot(sample_counts, aes(Timepoint, n_samples, fill = Timepoint)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = n_samples), vjust = -0.35, size = 1.9, family = "Arial", colour = "#333333") +
  scale_fill_manual(values = pal_time, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "e  SCFA sample coverage",
    x = NULL,
    y = "Samples"
  ) +
  theme(panel.grid.major.y = element_line(linewidth = 0.15, colour = "#E9E9E9"))

pF <- trend %>%
  mutate(Weight = if_else(Direction > 0, "Positive", "Negative")) %>%
  ggplot(aes(fct_reorder(SCFALabel, Slope), Slope, fill = Evidence)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.22) +
  coord_flip() +
  scale_fill_manual(values = pal_effect, guide = "none") +
  labs(
    title = "f  Score-aligned SCFA set",
    x = NULL,
    y = "Slope used for alignment"
  ) +
  theme(
    axis.text.y = element_text(size = 5.2),
    panel.grid.major.x = element_line(linewidth = 0.15, colour = "#E9E9E9")
  )

fig <- wrap_plots(
  A = pA,
  B = pB,
  C = pC,
  D = pD,
  E = pE,
  F = pF,
  design = "
AAAA
BBCC
DDDD
EEFF
"
) +
  plot_layout(heights = c(1.05, 0.98, 0.74, 0.78)) &
  theme(plot.margin = margin(3, 4, 3, 4))

base <- file.path(out_dir, "Fig4_scfa_maturation_redrawn")
w <- 183 / 25.4
h <- 255 / 25.4

ggsave(paste0(base, ".png"), fig, width = w, height = h, dpi = 600, bg = "white")
ggsave(paste0(base, ".pdf"), fig, width = w, height = h, device = cairo_pdf, bg = "white")
svglite::svglite(paste0(base, ".svg"), width = w, height = h, bg = "white")
print(fig)
dev.off()
ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = 600, background = "white")
print(fig)
dev.off()

message("Exported Fig. 4 to: ", out_dir)
