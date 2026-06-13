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
out_dir <- file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig5_milk_scfa_lagged_candidates")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

screen_path <- file.path(root, "direction3_multiomics_blueprint", "tables", "direction3_lagged_milk_infant_association_screen.csv")
model_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_hmo_ltf_to_scfa_lagged_model_highlights.csv")

time_windows <- c("D05 -> D14", "D05 -> D30", "D05 -> D90", "D14 -> D30", "D14 -> D90")
scfa_order <- c("SCFA:AA", "SCFA:PA", "SCFA:BA", "SCFA:IBA", "SCFA:IVA", "SCFA:2MBA_2-", "SCFA:VA", "SCFA:HA", "SCFA:SA")
scfa_labels <- c(
  "SCFA:AA" = "Acetate",
  "SCFA:PA" = "Propionate",
  "SCFA:BA" = "Butyrate",
  "SCFA:IBA" = "Isobutyrate",
  "SCFA:IVA" = "Isovalerate",
  "SCFA:2MBA_2-" = "2-Methylbutyrate",
  "SCFA:VA" = "Valerate",
  "SCFA:HA" = "Hexanoate",
  "SCFA:SA" = "Succinate"
)
scfa_abbrev <- c(
  "SCFA:AA" = "AA",
  "SCFA:PA" = "PA",
  "SCFA:BA" = "BA",
  "SCFA:IBA" = "IBA",
  "SCFA:IVA" = "IVA",
  "SCFA:2MBA_2-" = "2-MBA",
  "SCFA:VA" = "VA",
  "SCFA:HA" = "HA",
  "SCFA:SA" = "SA"
)
exposure_order <- c("LTF:Lactoferrin", "HMO:3-FL", "HMO:LNT", "HMO:LNFP-", "HMO:LDFT", "HMO:3-SL", "HMO:6-SL", "HMO:LNnT")
exposure_labels <- c(
  "LTF:Lactoferrin" = "Lactoferrin",
  "HMO:3-FL" = "3-FL",
  "HMO:LNT" = "LNT",
  "HMO:LNFP-" = "LNFP-",
  "HMO:LDFT" = "LDFT",
  "HMO:3-SL" = "3-SL",
  "HMO:6-SL" = "6-SL",
  "HMO:LNnT" = "LNnT"
)

pal_assay <- c(HMO = "#35618F", LTF = "#9A6A33")
pal_model <- c("Basic model" = "#305F8F", "Clinical-adjusted" = "#7A6A9E")
pal_count <- c(Screen = "#8BA6C4", Formal = "#D09A45", Adjusted = "#7A6A9E")

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

theme_pub <- function(base_size = 6.3) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks = element_line(linewidth = 0.28, colour = "black"),
      axis.text = element_text(colour = "#222222"),
      plot.title = element_text(size = base_size + 1.1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, colour = "#555555", hjust = 0),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.4),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(3, 4, 3, 4)
    )
}

screen <- read_csv(screen_path, show_col_types = FALSE) %>%
  filter(OutcomeType == "SCFA", ExposureAssay %in% c("HMO", "LTF")) %>%
  mutate(
    Contrast = factor(Contrast, levels = time_windows),
    ExposureLabel = factor(ExposureLabel, levels = exposure_order),
    OutcomeLabel = factor(OutcomeLabel, levels = scfa_order),
    ExposureClean = recode(as.character(ExposureLabel), !!!exposure_labels),
    OutcomeClean = recode(as.character(OutcomeLabel), !!!scfa_labels),
    ScreenEvidence = case_when(
      FDR < 0.10 ~ "q < 0.10",
      PValue < 0.01 ~ "P < 0.01",
      PValue < 0.05 ~ "P < 0.05",
      TRUE ~ "NS"
    ),
    ScreenEvidence = factor(ScreenEvidence, levels = c("q < 0.10", "P < 0.01", "P < 0.05", "NS")),
    NegLogP = pmin(pmax(-log10(PValue), -log10(0.05)), -log10(0.002))
  )

models <- read_csv(model_path, show_col_types = FALSE) %>%
  filter(OutcomeType == "SCFA", ModelOk) %>%
  mutate(
    Contrast = factor(Contrast, levels = time_windows),
    ExposureAssay = if_else(str_detect(ExposureLabel, "^HMO:"), "HMO", "LTF"),
    ExposureClean = recode(ExposureLabel, !!!exposure_labels),
    OutcomeClean = recode(OutcomeLabel, !!!scfa_labels),
    LinkClean = paste0(ExposureClean, " -> ", OutcomeClean),
    ModelClean = recode(Model, lagged_lm_basic = "Basic model", lagged_lm_clinical_adjusted = "Clinical-adjusted"),
    ModelClean = factor(ModelClean, levels = c("Basic model", "Clinical-adjusted")),
    ModelEvidence = case_when(
      PValueModel < 0.01 ~ "P < 0.01",
      PValueModel < 0.05 ~ "P < 0.05",
      PValueModel < 0.10 ~ "P < 0.10",
      TRUE ~ "NS"
    ),
    ModelEvidence = factor(ModelEvidence, levels = c("P < 0.01", "P < 0.05", "P < 0.10", "NS")),
    CI_low = Beta - 1.96 * SE,
    CI_high = Beta + 1.96 * SE,
    LinkWindow = paste0(LinkClean, "\n", as.character(Contrast)),
    ModelRank = row_number()
  )

candidate_ids <- models %>%
  filter(PValueModel < 0.10) %>%
  arrange(PValueModel) %>%
  slice_head(n = 16) %>%
  pull(ModelRank)

forest_df <- models %>%
  filter(ModelRank %in% candidate_ids) %>%
  arrange(Beta) %>%
  mutate(
    LinkWindow = factor(LinkWindow, levels = LinkWindow),
    PointSize = rescale(NModel, to = c(1.6, 3.4)),
    Label = paste0("N=", NModel, "; screen q=", sprintf("%.2f", FDRScreen))
  )

model_overlay <- models %>%
  filter(ModelRank %in% candidate_ids) %>%
  mutate(
    ExposureLabel = factor(ExposureLabel, levels = exposure_order),
    OutcomeLabel = factor(OutcomeLabel, levels = scfa_order)
  )

count_df <- screen %>%
  group_by(Contrast) %>%
  summarise(`Screen P<0.05` = sum(PValue < 0.05, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    models %>%
      group_by(Contrast) %>%
      summarise(
        `Formal model P<0.05` = sum(PValueModel < 0.05, na.rm = TRUE),
        `Adjusted model P<0.05` = sum(PValueModel < 0.05 & ModelClean == "Clinical-adjusted", na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Contrast"
  ) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
  pivot_longer(-Contrast, names_to = "Layer", values_to = "Count") %>%
  mutate(
    Layer = factor(Layer, levels = c("Screen P<0.05", "Formal model P<0.05", "Adjusted model P<0.05")),
    LayerShort = recode(as.character(Layer), "Screen P<0.05" = "Screen", "Formal model P<0.05" = "Formal", "Adjusted model P<0.05" = "Adjusted"),
    LayerShort = factor(LayerShort, levels = c("Screen", "Formal", "Adjusted")),
    Contrast = factor(Contrast, levels = time_windows)
  )

coverage_df <- models %>%
  mutate(Window = as.character(Contrast)) %>%
  group_by(Window) %>%
  summarise(
    MinN = min(NModel, na.rm = TRUE),
    MedianN = median(NModel, na.rm = TRUE),
    MaxN = max(NModel, na.rm = TRUE),
    Candidates = n(),
    MinScreenFDR = min(FDRScreen, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Window = factor(Window, levels = time_windows))

write_csv(screen, file.path(out_dir, "Fig5_source_lagged_screen_hmo_ltf_scfa.csv"))
write_csv(models, file.path(out_dir, "Fig5_source_formal_model_candidates.csv"))
write_csv(forest_df, file.path(out_dir, "Fig5_source_forest_candidates.csv"))
write_csv(count_df, file.path(out_dir, "Fig5_source_evidence_counts.csv"))
write_csv(coverage_df, file.path(out_dir, "Fig5_source_candidate_coverage.csv"))

pA <- ggplot(screen, aes(OutcomeLabel, ExposureLabel)) +
  geom_tile(aes(fill = Rho), colour = "white", linewidth = 0.18) +
  geom_point(
    data = filter(screen, PValue < 0.05),
    aes(size = NegLogP),
    shape = 21,
    colour = "black",
    fill = "white",
    stroke = 0.18,
    alpha = 0.88,
    show.legend = TRUE
  ) +
  geom_point(
    data = model_overlay,
    colour = "#202020",
    size = 1.15,
    stroke = 0.42,
    inherit.aes = FALSE,
    aes(OutcomeLabel, ExposureLabel, shape = "Formal model")
  ) +
  facet_wrap(~Contrast, nrow = 1) +
  scale_fill_gradient2(low = "#9B4F5A", mid = "#F5F5F2", high = "#2C6E91", midpoint = 0, limits = c(-0.8, 0.8), name = "Spearman rho") +
  scale_size_continuous(
    range = c(0.55, 1.85),
    limits = c(-log10(0.05), -log10(0.002)),
    breaks = c(-log10(0.05), -log10(0.01), -log10(0.002)),
    labels = c("0.05", "0.01", "0.002"),
    name = "Screen P"
  ) +
  scale_shape_manual(values = c("Formal model" = 4), name = NULL) +
  scale_x_discrete(labels = scfa_abbrev) +
  scale_y_discrete(labels = exposure_labels) +
  labs(
    title = "a  Lagged HMO/LTF-SCFA screen",
    x = NULL,
    y = NULL
  ) +
  theme_pub(5.6) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 5.0),
    axis.text.y = element_text(size = 5.0),
    panel.spacing.x = unit(2.2, "mm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(7, "mm")
  ) +
  guides(
    fill = guide_colourbar(order = 1, barwidth = 4.5, barheight = 0.35, title.position = "top"),
    size = guide_legend(order = 2, title.position = "top", override.aes = list(shape = 21, fill = "white")),
    shape = guide_legend(order = 3, override.aes = list(size = 2.0, colour = "#202020"))
  )

pB <- ggplot(forest_df, aes(Beta, LinkWindow)) +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "#A9A9A9") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, colour = ExposureAssay), height = 0, linewidth = 0.36, alpha = 0.82) +
  geom_point(aes(colour = ExposureAssay, shape = ModelClean, size = NModel), alpha = 0.96, stroke = 0.4) +
  geom_text(aes(label = if_else(PValueModel < 0.05, "*", "."), x = if_else(Beta >= 0, pmin(CI_high + 0.06, 1.06), pmax(CI_low - 0.06, -1.06))), size = 2.25, family = "Arial", colour = "#202020") +
  scale_colour_manual(values = pal_assay, name = NULL) +
  scale_shape_manual(values = c("Basic model" = 16, "Clinical-adjusted" = 17), name = NULL) +
  scale_size_continuous(range = c(1.45, 3.0), breaks = c(15, 21, 26), name = "N") +
  scale_x_continuous(limits = c(-1.12, 1.22), breaks = seq(-1, 1, 0.5)) +
  labs(
    title = "b  Candidate beta estimates",
    x = "Standardized beta",
    y = NULL
  ) +
  theme_pub(5.9) +
  theme(
    axis.text.y = element_text(size = 4.9, lineheight = 0.88),
    legend.position = "top",
    legend.justification = "left",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 1, 0),
    legend.key.size = unit(2.8, "mm"),
    legend.text = element_text(size = 5.1),
    panel.grid.major.y = element_line(linewidth = 0.14, colour = "#E9E9E9")
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 2.0)),
    shape = guide_legend(order = 2, override.aes = list(size = 2.1)),
    size = guide_legend(order = 3, override.aes = list(shape = 16, colour = "#555555"))
  )

pC <- ggplot(count_df, aes(Contrast, Count, fill = LayerShort)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "white", linewidth = 0.18) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.72), vjust = -0.35, size = 1.8, family = "Arial") +
  scale_fill_manual(values = pal_count, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16)), breaks = pretty_breaks(n = 4)) +
  labs(
    title = "c  Signal attrition",
    x = NULL,
    y = "Links"
  ) +
  theme_pub(5.8) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 5.0),
    panel.grid.major.y = element_line(linewidth = 0.14, colour = "#E6E6E6")
  )

pD <- ggplot(coverage_df, aes(Window, MedianN)) +
  geom_linerange(aes(ymin = MinN, ymax = MaxN), linewidth = 0.52, colour = "#595959") +
  geom_point(aes(size = Candidates, colour = MinScreenFDR), alpha = 0.95) +
  geom_text(aes(label = paste0("n=", Candidates, "\nq=", sprintf("%.2f", MinScreenFDR))), nudge_y = 2.2, size = 1.85, family = "Arial", colour = "#4F4F4F") +
  scale_colour_gradient(low = "#4E79A7", high = "#B7B7B7", limits = c(0.10, 1.0), oob = squish, guide = "none") +
  scale_size_continuous(range = c(2.0, 4.5), guide = "none") +
  scale_y_continuous(limits = c(12, 31), breaks = seq(12, 28, 4)) +
  labs(
    title = "d  Coverage",
    x = NULL,
    y = "Model N"
  ) +
  theme_pub(5.8) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.major.y = element_line(linewidth = 0.14, colour = "#E6E6E6"),
    legend.position = "none"
  )

fig <- (pA / (pB | (pC / pD))) +
  plot_layout(heights = c(1.04, 1.16), widths = c(1.18, 0.82)) &
  theme(plot.margin = margin(3, 4, 3, 4))

base_file <- file.path(out_dir, "Fig5_milk_scfa_lagged_candidates_redrawn")
ggsave(paste0(base_file, ".png"), fig, width = 183, height = 165, units = "mm", dpi = 600, bg = "white")
ggsave(paste0(base_file, ".pdf"), fig, width = 183, height = 165, units = "mm", device = cairo_pdf, bg = "white")
svglite::svglite(paste0(base_file, ".svg"), width = 183 / 25.4, height = 165 / 25.4, bg = "white")
print(fig)
dev.off()
ragg::agg_tiff(paste0(base_file, ".tiff"), width = 183 / 25.4, height = 165 / 25.4, units = "in", res = 600, background = "white")
print(fig)
dev.off()

message("Wrote Fig5 redraw outputs to: ", out_dir)
