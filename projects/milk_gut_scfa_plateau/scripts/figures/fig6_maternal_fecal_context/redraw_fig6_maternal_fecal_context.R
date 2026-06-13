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
out_dir <- file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig6_maternal_fecal_context")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_mf_bf_all_pairwise_similarity.csv")
summary_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_mf_bf_own_vs_other_similarity_summary.csv")

time_levels <- c("D05", "D14", "D30", "D90")
mf_levels <- c("D30", "MF_PRE")
mf_labels <- c(D30 = "MF post-1m", MF_PRE = "MF pre")
pal_time <- c(D05 = "#5778A4", D14 = "#79A88D", D30 = "#D89C3D", D90 = "#B56576")
pal_mf <- c("MF post-1m" = "#5F7FA3", "MF pre" = "#9B8B76")
pal_enrich <- c("Empirical P <= 0.10" = "#2C6E91", "P > 0.10" = "#B8B8B8")

theme_pub <- function(base_size = 6.3) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks = element_line(linewidth = 0.28, colour = "black"),
      axis.text = element_text(colour = "#222222"),
      plot.title = element_text(size = base_size + 1.1, face = "bold", hjust = 0),
      legend.title = element_text(size = base_size - 0.1),
      legend.text = element_text(size = base_size - 0.4),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(3, 4, 3, 4)
    )
}

all_sim <- read_csv(all_path, show_col_types = FALSE) %>%
  filter(MF_Window %in% mf_levels, BF_Timepoint %in% time_levels) %>%
  mutate(
    BF_Timepoint = factor(BF_Timepoint, levels = time_levels),
    MF_Window = factor(MF_Window, levels = mf_levels),
    MF_Label = factor(recode(as.character(MF_Window), !!!mf_labels), levels = mf_labels[mf_levels]),
    PairType = factor(PairType, levels = c("Other dyads", "Own dyad"))
  )

own <- read_csv(summary_path, show_col_types = FALSE) %>%
  filter(MF_Window %in% mf_levels, BF_Timepoint %in% time_levels) %>%
  mutate(
    BF_Timepoint = factor(BF_Timepoint, levels = time_levels),
    MF_Window = factor(MF_Window, levels = mf_levels),
    MF_Label = factor(recode(as.character(MF_Window), !!!mf_labels), levels = mf_labels[mf_levels]),
    Enriched = if_else(EmpiricalP_BrayGreater <= 0.10, "Empirical P <= 0.10", "P > 0.10"),
    Enriched = factor(Enriched, levels = c("Empirical P <= 0.10", "P > 0.10")),
    DyadLabel = paste0("D", DyadID_BF)
  )

window_summary <- own %>%
  group_by(MF_Label, BF_Timepoint) %>%
  summarise(
    OwnN = n(),
    MedianOwn = median(OwnBray, na.rm = TRUE),
    MedianOther = median(OtherMedianBray, na.rm = TRUE),
    MedianZ = median(OwnVsOtherZ, na.rm = TRUE),
    EnrichedN = sum(EmpiricalP_BrayGreater <= 0.10, na.rm = TRUE),
    MinEmpiricalP = min(EmpiricalP_BrayGreater, na.rm = TRUE),
    MaxZ = max(OwnVsOtherZ, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    BF_Timepoint = factor(BF_Timepoint, levels = time_levels),
    MF_Label = factor(MF_Label, levels = mf_labels[mf_levels])
  )

z_summary <- own %>%
  group_by(MF_Label, BF_Timepoint) %>%
  summarise(
    OwnN = n(),
    MedianZ = median(OwnVsOtherZ, na.rm = TRUE),
    Q25Z = quantile(OwnVsOtherZ, 0.25, na.rm = TRUE),
    Q75Z = quantile(OwnVsOtherZ, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    BF_Timepoint = factor(BF_Timepoint, levels = time_levels),
    MF_Label = factor(MF_Label, levels = mf_labels[mf_levels])
  )

other_density <- all_sim %>%
  filter(PairType == "Other dyads") %>%
  group_by(MF_Label, BF_Timepoint) %>%
  summarise(
    OtherN = n(),
    OtherMedian = median(BraySimilarity, na.rm = TRUE),
    OtherQ25 = quantile(BraySimilarity, 0.25, na.rm = TRUE),
    OtherQ75 = quantile(BraySimilarity, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

top_dyads <- own %>%
  arrange(EmpiricalP_BrayGreater, desc(OwnVsOtherZ)) %>%
  slice_head(n = 10) %>%
  mutate(
    RankLabel = paste0(DyadLabel, " | ", as.character(BF_Timepoint), " / ", as.character(MF_Label)),
    RankLabel = factor(RankLabel, levels = rev(RankLabel))
  )

coverage <- all_sim %>%
  group_by(MF_Label, BF_Timepoint, PairType) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(names_from = PairType, values_from = N, values_fill = 0) %>%
  left_join(window_summary %>% select(MF_Label, BF_Timepoint, EnrichedN, MinEmpiricalP), by = c("MF_Label", "BF_Timepoint")) %>%
  mutate(
    Label = paste0("own ", `Own dyad`, "\nother ", `Other dyads`),
    PointSize = pmax(`Own dyad`, 1)
  )

write_csv(all_sim, file.path(out_dir, "Fig6_source_all_pairwise_similarity.csv"))
write_csv(own, file.path(out_dir, "Fig6_source_own_vs_other_summary.csv"))
write_csv(window_summary, file.path(out_dir, "Fig6_source_window_summary.csv"))
write_csv(z_summary, file.path(out_dir, "Fig6_source_z_iqr_summary.csv"))
write_csv(top_dyads, file.path(out_dir, "Fig6_source_top_dyads.csv"))
write_csv(coverage, file.path(out_dir, "Fig6_source_coverage.csv"))

pA <- ggplot() +
  geom_boxplot(
    data = filter(all_sim, PairType == "Other dyads"),
    aes(BF_Timepoint, BraySimilarity, fill = BF_Timepoint),
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.22,
    colour = "#6F6F6F",
    linewidth = 0.32
  ) +
  geom_point(
    data = filter(all_sim, PairType == "Other dyads"),
    aes(BF_Timepoint, BraySimilarity),
    position = position_jitter(width = 0.13, height = 0),
    size = 0.32,
    alpha = 0.14,
    colour = "#5F5F5F"
  ) +
  geom_point(
    data = own,
    aes(BF_Timepoint, OwnBray, colour = BF_Timepoint, shape = Enriched),
    position = position_jitter(width = 0.045, height = 0),
    size = 1.35,
    stroke = 0.42,
    alpha = 0.95
  ) +
  facet_wrap(~MF_Label, nrow = 1) +
  scale_fill_manual(values = pal_time, guide = "none") +
  scale_colour_manual(values = pal_time, name = NULL) +
  scale_shape_manual(values = c("Empirical P <= 0.10" = 17, "P > 0.10" = 16), name = NULL) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0.01, 0.03))) +
  coord_cartesian(ylim = c(0, 0.88)) +
  labs(
    title = "a  Pairwise MF-BF similarity",
    x = NULL,
    y = "Bray-Curtis similarity"
  ) +
  theme_pub(5.9) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(2.8, "mm"),
    panel.grid.major.y = element_line(linewidth = 0.14, colour = "#E8E8E8")
  )

pB <- ggplot(own, aes(BF_Timepoint, OwnVsOtherZ)) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "#9A9A9A") +
  geom_linerange(
    data = z_summary,
    aes(BF_Timepoint, ymin = Q25Z, ymax = Q75Z),
    inherit.aes = FALSE,
    linewidth = 0.55,
    colour = "#555555",
    alpha = 0.82
  ) +
  geom_point(
    data = z_summary,
    aes(BF_Timepoint, MedianZ),
    inherit.aes = FALSE,
    shape = 95,
    size = 7,
    colour = "#333333"
  ) +
  geom_point(
    aes(fill = Enriched),
    shape = 21,
    colour = "#303030",
    stroke = 0.22,
    position = position_jitter(width = 0.065, height = 0),
    size = 1.35,
    alpha = 0.94
  ) +
  facet_wrap(~MF_Label, nrow = 1) +
  scale_fill_manual(values = pal_enrich, name = NULL) +
  scale_y_continuous(breaks = c(-1, 0, 1, 3, 5), expand = expansion(mult = c(0.03, 0.05))) +
  coord_cartesian(ylim = c(-1.25, 6.45)) +
  labs(
    title = "b  Own-vs-other Z",
    x = NULL,
    y = "Z score"
  ) +
  theme_pub(5.9) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(2.8, "mm"),
    panel.grid.major.y = element_line(linewidth = 0.14, colour = "#E8E8E8")
  )

pC <- ggplot(top_dyads, aes(OwnVsOtherZ, RankLabel)) +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "#A0A0A0") +
  geom_segment(aes(x = 0, xend = OwnVsOtherZ, yend = RankLabel, colour = MF_Label), linewidth = 0.55) +
  geom_point(aes(size = -log10(EmpiricalP_BrayGreater), colour = MF_Label, shape = BF_Timepoint), alpha = 0.95) +
  scale_colour_manual(values = pal_mf, name = NULL) +
  scale_shape_manual(values = c(D05 = 16, D14 = 17, D30 = 15, D90 = 18), guide = "none") +
  scale_size_continuous(
    range = c(1.7, 3.2),
    breaks = -log10(c(0.10, 0.05)),
    labels = c("0.10", "0.05"),
    name = "Empirical P"
  ) +
  scale_x_continuous(breaks = c(0, 2, 4, 6), expand = expansion(mult = c(0.01, 0.04))) +
  coord_cartesian(xlim = c(-0.25, 6.5)) +
  labs(
    title = "c  Dyad-specific outliers",
    x = "Own-vs-other Z",
    y = NULL
  ) +
  theme_pub(5.7) +
  theme(
    axis.text.y = element_text(size = 4.8),
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(2.8, "mm"),
    panel.grid.major.x = element_line(linewidth = 0.14, colour = "#E8E8E8")
  )

pD <- ggplot(window_summary, aes(BF_Timepoint, MedianZ, fill = MF_Label)) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "#9A9A9A") +
  geom_col(position = position_dodge(width = 0.7), width = 0.58, colour = "white", linewidth = 0.18, alpha = 0.92) +
  geom_text(
    aes(label = paste0("n=", OwnN, "\nP<=.10 ", EnrichedN)),
    position = position_dodge(width = 0.7),
    vjust = if_else(window_summary$MedianZ >= 0, -0.25, 1.15),
    size = 1.55,
    family = "Arial",
    colour = "#4F4F4F"
  ) +
  scale_fill_manual(values = pal_mf, name = NULL) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5), expand = expansion(mult = c(0.08, 0.10))) +
  coord_cartesian(ylim = c(-0.8, 0.8), clip = "off") +
  labs(
    title = "d  Window summary",
    x = NULL,
    y = "Median Z"
  ) +
  theme_pub(5.7) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(2.8, "mm"),
    panel.grid.major.y = element_line(linewidth = 0.14, colour = "#E8E8E8")
  )

fig <- (pA / (pB | (pC / pD))) +
  plot_layout(heights = c(0.95, 1.15), widths = c(1.02, 0.98)) &
  theme(plot.margin = margin(3, 4, 3, 4))

base_file <- file.path(out_dir, "Fig6_maternal_fecal_context_redrawn")
ggsave(paste0(base_file, ".png"), fig, width = 183, height = 165, units = "mm", dpi = 600, bg = "white")
ggsave(paste0(base_file, ".pdf"), fig, width = 183, height = 165, units = "mm", device = cairo_pdf, bg = "white")
svglite::svglite(paste0(base_file, ".svg"), width = 183 / 25.4, height = 165 / 25.4, bg = "white")
print(fig)
dev.off()
ragg::agg_tiff(paste0(base_file, ".tiff"), width = 183 / 25.4, height = 165 / 25.4, units = "in", res = 600, background = "white")
print(fig)
dev.off()

message("Wrote Fig6 redraw outputs to: ", out_dir)
