library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(patchwork)
library(svglite)
library(ragg)

root <- normalizePath(Sys.getenv("FANLAB_RESEARCH2_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig1_cohort_design_overview")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

time_levels <- c("D05", "D14", "D30", "D90")

add_record <- function(dyad, row_key, module_family, row_label, timepoint = NA_character_) {
  tibble(
    DyadID = as.character(dyad),
    RowKey = row_key,
    ModuleFamily = module_family,
    RowLabel = row_label,
    Timepoint = timepoint,
    Present = TRUE
  )
}

milk <- read_csv(file.path(root, "p2_bridge_analyses", "tables", "p2_milk_module_clinical_matrix.csv"), show_col_types = FALSE) %>%
  mutate(DyadID = as.character(DyadID), Timepoint = factor(Timepoint, levels = time_levels))

scfa <- read_csv(file.path(root, "p0_multiomics_extension", "tables", "p0_scfa_maturation_score.csv"), show_col_types = FALSE) %>%
  mutate(DyadID = as.character(DyadID), TimepointCode = factor(TimepointCode, levels = time_levels))

clinical <- read_csv(file.path(root, "clinical_scfa_metagenome_report", "tables", "clinical_source_complete_utf8.csv"), show_col_types = FALSE) %>%
  mutate(DyadID = as.character(DyadID))

metag <- read_csv(file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig3_infant_microbiome_maturation", "Fig3_source_pcoa.csv"), show_col_types = FALSE) %>%
  mutate(
    DyadID = str_match(SampleID, "^BF_([0-9]+)_")[, 2],
    Timepoint = factor(Timepoint, levels = time_levels)
  ) %>%
  filter(!is.na(DyadID))

mf <- read_csv(file.path(root, "p0_multiomics_extension", "tables", "p0_mf_bf_own_vs_other_similarity_summary.csv"), show_col_types = FALSE) %>%
  transmute(
    DyadID = as.character(DyadID_BF),
    BF_Timepoint = factor(BF_Timepoint, levels = time_levels),
    MF_Window = as.character(MF_Window)
  )

baseline_dyads <- clinical %>%
  filter(!is.na(DyadID), DyadID != "") %>%
  distinct(DyadID)

records <- bind_rows(
  baseline_dyads %>%
    rowwise() %>%
    do(add_record(.$DyadID, "00_clinical_baseline", "Clinical", "Clinical baseline")) %>%
    ungroup(),
  milk %>%
    filter(Module == "LTF", !is.na(Value), Timepoint %in% time_levels) %>%
    distinct(DyadID, Timepoint) %>%
    transmute(DyadID, RowKey = paste0("10_milk_ltf_", Timepoint), ModuleFamily = "Milk", RowLabel = paste0("Milk LTF ", Timepoint), Timepoint = as.character(Timepoint), Present = TRUE),
  milk %>%
    filter(Module == "HMO_total", !is.na(Value), Timepoint %in% c("D05", "D14", "D30")) %>%
    distinct(DyadID, Timepoint) %>%
    transmute(DyadID, RowKey = paste0("11_milk_hmo_", Timepoint), ModuleFamily = "Milk", RowLabel = paste0("Milk HMO ", Timepoint), Timepoint = as.character(Timepoint), Present = TRUE),
  milk %>%
    filter(Module == "LCFA_total", !is.na(Value), Timepoint %in% time_levels) %>%
    distinct(DyadID, Timepoint) %>%
    transmute(DyadID, RowKey = paste0("12_milk_lcfa_", Timepoint), ModuleFamily = "Milk", RowLabel = paste0("Milk LCFA ", Timepoint), Timepoint = as.character(Timepoint), Present = TRUE),
  scfa %>%
    filter(!is.na(SCFA_MaturationScore), TimepointCode %in% time_levels) %>%
    distinct(DyadID, TimepointCode) %>%
    transmute(DyadID, RowKey = paste0("20_infant_scfa_", TimepointCode), ModuleFamily = "Infant feces", RowLabel = paste0("Infant SCFA ", TimepointCode), Timepoint = as.character(TimepointCode), Present = TRUE),
  metag %>%
    distinct(DyadID, Timepoint) %>%
    transmute(DyadID, RowKey = paste0("21_infant_metag_", Timepoint), ModuleFamily = "Infant feces", RowLabel = paste0("Infant metagenome ", Timepoint), Timepoint = as.character(Timepoint), Present = TRUE),
  mf %>%
    filter(MF_Window == "MF_PRE") %>%
    distinct(DyadID) %>%
    transmute(DyadID, RowKey = "30_maternal_metag_pre", ModuleFamily = "Maternal feces", RowLabel = "Maternal metagenome pre", Timepoint = "MF_PRE", Present = TRUE),
  mf %>%
    filter(MF_Window == "D30") %>%
    distinct(DyadID) %>%
    transmute(DyadID, RowKey = "31_maternal_metag_post1m", ModuleFamily = "Maternal feces", RowLabel = "Maternal metagenome post-1m", Timepoint = "D30", Present = TRUE)
) %>%
  filter(!is.na(DyadID), DyadID != "")

row_order <- records %>%
  distinct(RowKey, RowLabel, ModuleFamily) %>%
  arrange(RowKey) %>%
  mutate(RowLabel = factor(RowLabel, levels = rev(RowLabel)))

dyad_order <- records %>%
  distinct(DyadID, RowKey) %>%
  count(DyadID, name = "CoverageCells") %>%
  arrange(desc(CoverageCells), suppressWarnings(as.integer(DyadID)), DyadID) %>%
  mutate(DyadPlot = factor(DyadID, levels = DyadID))

plot_data <- tidyr::expand_grid(DyadID = dyad_order$DyadID, RowKey = row_order$RowKey) %>%
  left_join(row_order, by = "RowKey") %>%
  left_join(records %>% distinct(DyadID, RowKey, Present), by = c("DyadID", "RowKey")) %>%
  left_join(dyad_order, by = "DyadID") %>%
  mutate(
    Present = !is.na(Present),
    DyadPlot = factor(DyadID, levels = dyad_order$DyadID),
    CellFamily = if_else(Present, ModuleFamily, "Missing")
  )

row_counts <- plot_data %>%
  filter(Present) %>%
  count(RowKey, RowLabel, ModuleFamily, name = "n") %>%
  mutate(RowLabel = factor(RowLabel, levels = levels(row_order$RowLabel)))

dyad_counts <- plot_data %>%
  filter(Present) %>%
  count(DyadID, DyadPlot, name = "n_present")

module_counts <- records %>%
  distinct(ModuleFamily, RowKey, RowLabel, DyadID) %>%
  count(ModuleFamily, RowKey, RowLabel, name = "n_dyads") %>%
  arrange(RowKey)

write_csv(records %>% arrange(RowKey, DyadID), file.path(out_dir, "Fig1_source_coverage_long.csv"))
write_csv(plot_data %>% arrange(RowKey, DyadID), file.path(out_dir, "Fig1_source_coverage_matrix.csv"))
write_csv(row_counts %>% arrange(RowKey), file.path(out_dir, "Fig1_source_coverage_row_counts.csv"))
write_csv(dyad_counts %>% arrange(desc(n_present), DyadID), file.path(out_dir, "Fig1_source_coverage_dyad_counts.csv"))
write_csv(module_counts, file.path(out_dir, "Fig1_source_coverage_module_counts.csv"))

pal <- c(
  "Clinical" = "#6D7A86",
  "Milk" = "#D9822B",
  "Infant feces" = "#168C88",
  "Maternal feces" = "#5E6FA3",
  "Missing" = "#F3F5F7"
)

p_top <- ggplot(dyad_counts, aes(x = DyadPlot, y = n_present)) +
  geom_col(width = 0.82, fill = "#174A67") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Cells per dyad") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(size = 6.5),
    axis.text.y = element_text(size = 5.8),
    plot.margin = margin(2, 2, 0, 2)
  )

p_heat <- ggplot(plot_data, aes(x = DyadPlot, y = RowLabel, fill = CellFamily)) +
  geom_tile(width = 0.92, height = 0.86, colour = "white", linewidth = 0.12) +
  scale_fill_manual(values = pal, breaks = c("Clinical", "Milk", "Infant feces", "Maternal feces"), name = NULL) +
  labs(x = "Dyads ordered by multi-omics completeness", y = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6.3, colour = "black"),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(10, "pt"),
    legend.key.height = unit(6, "pt"),
    legend.text = element_text(size = 6.2),
    plot.margin = margin(0, 2, 2, 2)
  )

p_right <- ggplot(row_counts, aes(x = n, y = RowLabel, fill = ModuleFamily)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n), hjust = -0.12, size = 1.9, family = "Arial") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "n dyads", y = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(size = 6.5),
    axis.text.x = element_text(size = 5.8),
    plot.margin = margin(0, 2, 2, 0)
  )

coverage_panel <- (p_top + plot_spacer()) / (p_heat + p_right) +
  plot_layout(widths = c(1, 0.18), heights = c(0.28, 1)) +
  plot_annotation(
    title = "B  Modular longitudinal coverage supports analysis by data layer",
    subtitle = "Each column is a mother-infant dyad; filled cells indicate available module-window data.",
    theme = theme(
      plot.title = element_text(size = 9, face = "bold", family = "Arial"),
      plot.subtitle = element_text(size = 6.5, family = "Arial", colour = "#4C5560"),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

base_name <- file.path(out_dir, "Fig1B_coverage_heatmap")
ggsave(paste0(base_name, ".png"), coverage_panel, width = 183, height = 95, units = "mm", dpi = 600, bg = "white")
ggsave(paste0(base_name, ".pdf"), coverage_panel, width = 183, height = 95, units = "mm", device = cairo_pdf, bg = "white")
svglite(paste0(base_name, ".svg"), width = 183 / 25.4, height = 95 / 25.4, bg = "white")
print(coverage_panel)
dev.off()
ragg::agg_tiff(paste0(base_name, ".tiff"), width = 183 / 25.4, height = 95 / 25.4, units = "in", res = 600, background = "white")
print(coverage_panel)
dev.off()

message("Wrote Fig1B coverage heatmap to: ", out_dir)
