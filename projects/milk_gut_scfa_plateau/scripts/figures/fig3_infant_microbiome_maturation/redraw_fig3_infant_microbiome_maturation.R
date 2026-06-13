library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(patchwork)
library(scales)
library(vegan)
library(grid)

root <- normalizePath(Sys.getenv("FANLAB_RESEARCH2_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig3_infant_microbiome_maturation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genus_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_metagenome_genus_abundance_long.csv")
species_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_metagenome_species_abundance_long.csv")
species_summary_path <- file.path(root, "p0_multiomics_extension", "tables", "p0_bifidobacterium_species_trajectory_summary.csv")

time_levels <- c("D05", "D14", "D30", "D90")
pal_time <- c(D05 = "#5778A4", D14 = "#79A88D", D30 = "#D89C3D", D90 = "#B56576")
pal_key <- c(
  Bifidobacterium = "#2A6FBB",
  Escherichia = "#C84D4D",
  Klebsiella = "#8E66AA",
  Streptococcus = "#D8A13B",
  Enterococcus = "#7A7A7A",
  Bacteroides = "#6E8B3D"
)
pal_species <- c(
  "Bifidobacterium longum" = "#245F9E",
  "Bifidobacterium bifidum" = "#2F9C95",
  "Bifidobacterium breve" = "#D99019",
  "Bifidobacterium catenulatum" = "#C03D4D",
  "Bifidobacterium pseudocatenulatum" = "#8058C7",
  "Bifidobacterium animalis" = "#666666"
)

safe_kw <- function(df) {
  if (length(unique(df$Timepoint)) < 2 || sd(df$RelAbundance, na.rm = TRUE) == 0) return(NA_real_)
  tryCatch(kruskal.test(RelAbundance ~ Timepoint, data = df)$p.value, error = function(e) NA_real_)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

genus_raw <- read_csv(genus_path, show_col_types = FALSE) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = time_levels),
    DyadID = as.character(DyadID),
    Abundance = as.numeric(Abundance)
  ) %>%
  filter(Source == "BF", !is.na(Timepoint), !is.na(Abundance))

species_raw <- read_csv(species_path, show_col_types = FALSE) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = time_levels),
    DyadID = as.character(DyadID),
    Abundance = as.numeric(Abundance)
  ) %>%
  filter(Source == "BF", !is.na(Timepoint), !is.na(Abundance), Genus == "Bifidobacterium")

species_summary <- read_csv(species_summary_path, show_col_types = FALSE) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = time_levels),
    PrevalencePercent = 100 * Prevalence
  )

sample_meta <- genus_raw %>%
  distinct(SampleID, DyadID, Timepoint)

sample_counts <- sample_meta %>%
  count(Timepoint, name = "n_samples") %>%
  mutate(label = paste0(as.character(Timepoint), "\nn=", n_samples))

all_genera <- sort(unique(genus_raw$Genus))

genus_summed <- genus_raw %>%
  group_by(SampleID, DyadID, Timepoint, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

genus_complete <- sample_meta %>%
  crossing(Genus = all_genera) %>%
  left_join(genus_summed, by = c("SampleID", "DyadID", "Timepoint", "Genus")) %>%
  mutate(Abundance = replace_na(Abundance, 0)) %>%
  group_by(SampleID) %>%
  mutate(
    SampleTotal = sum(Abundance, na.rm = TRUE),
    RelAbundance = if_else(SampleTotal > 0, 100 * Abundance / SampleTotal, 0)
  ) %>%
  ungroup()

genus_stats <- genus_complete %>%
  group_by(Genus) %>%
  summarise(
    mean_abundance = mean(RelAbundance, na.rm = TRUE),
    median_abundance = median(RelAbundance, na.rm = TRUE),
    prevalence = mean(RelAbundance > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance))

top_heat_genera <- genus_stats %>%
  slice_head(n = 14) %>%
  pull(Genus)

bifido_order <- genus_complete %>%
  filter(Genus == "Bifidobacterium") %>%
  select(SampleID, BifidoAbundance = RelAbundance)

sample_order <- sample_meta %>%
  left_join(bifido_order, by = "SampleID") %>%
  arrange(Timepoint, desc(BifidoAbundance), SampleID) %>%
  mutate(SamplePlot = factor(SampleID, levels = SampleID))

heat_df <- genus_complete %>%
  mutate(GenusPlot = if_else(Genus %in% top_heat_genera, Genus, "Other")) %>%
  group_by(SampleID, Timepoint, GenusPlot) %>%
  summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(sample_order %>% select(SampleID, SamplePlot), by = "SampleID") %>%
  mutate(
    GenusPlot = factor(GenusPlot, levels = rev(c(top_heat_genera, "Other"))),
    Timepoint = factor(Timepoint, levels = time_levels)
  )

time_strip_df <- sample_order %>%
  mutate(y = 1)

pcoa_matrix <- genus_complete %>%
  select(SampleID, Timepoint, Genus, RelAbundance) %>%
  pivot_wider(names_from = Genus, values_from = RelAbundance, values_fill = 0) %>%
  arrange(match(SampleID, sample_order$SampleID))

pcoa_meta <- pcoa_matrix %>% select(SampleID, Timepoint)
pcoa_abund <- pcoa_matrix %>% select(-SampleID, -Timepoint) %>% as.data.frame()
rownames(pcoa_abund) <- pcoa_matrix$SampleID
bray <- vegdist(pcoa_abund, method = "bray")
pcoa_fit <- cmdscale(bray, k = 2, eig = TRUE)
var_exp <- round(100 * pcoa_fit$eig[1:2] / sum(abs(pcoa_fit$eig)), 1)
set.seed(20260606)
adon <- adonis2(bray ~ Timepoint, data = pcoa_meta, permutations = 999)
adon_label <- paste0("PERMANOVA R2=", sprintf("%.2f", adon$R2[1]), ", P=", fmt_p(adon$`Pr(>F)`[1]))

pcoa_df <- pcoa_meta %>%
  bind_cols(as_tibble(pcoa_fit$points, .name_repair = ~ c("PCoA1", "PCoA2"))) %>%
  left_join(sample_counts, by = "Timepoint")

centroid_df <- pcoa_df %>%
  group_by(Timepoint) %>%
  summarise(
    PCoA1 = mean(PCoA1),
    PCoA2 = mean(PCoA2),
    .groups = "drop"
  ) %>%
  mutate(TimeIndex = as.numeric(Timepoint))

genus_effect <- bind_rows(lapply(all_genera, function(gn) {
  df <- genus_complete %>% filter(Genus == gn)
  med <- df %>%
    group_by(Timepoint) %>%
    summarise(median = median(RelAbundance, na.rm = TRUE), .groups = "drop")
  tibble(
    Genus = gn,
    mean_abundance = mean(df$RelAbundance, na.rm = TRUE),
    prevalence = mean(df$RelAbundance > 0, na.rm = TRUE),
    kw_p = safe_kw(df),
    d05 = med$median[med$Timepoint == "D05"][1],
    d30 = med$median[med$Timepoint == "D30"][1],
    d90 = med$median[med$Timepoint == "D90"][1]
  )
})) %>%
  mutate(
    d05 = replace_na(d05, 0),
    d30 = replace_na(d30, 0),
    d90 = replace_na(d90, 0),
    delta_d30_d05 = d30 - d05,
    delta_d90_d05 = d90 - d05,
    kw_fdr = p.adjust(kw_p, method = "BH"),
    Direction = case_when(delta_d30_d05 > 0 ~ "Higher at D30", delta_d30_d05 < 0 ~ "Lower at D30", TRUE ~ "No median shift")
  )

effect_df <- genus_effect %>%
  filter(mean_abundance >= 0.45 | Genus %in% names(pal_key)) %>%
  arrange(desc(abs(delta_d30_d05))) %>%
  slice_head(n = 18) %>%
  mutate(
    Genus = factor(Genus, levels = rev(Genus)),
    Direction = factor(Direction, levels = c("Higher at D30", "Lower at D30", "No median shift")),
    FDRClass = case_when(kw_fdr < 0.05 ~ "FDR < 0.05", kw_p < 0.05 ~ "P < 0.05", TRUE ~ "NS")
  )

trajectory_genera <- c("Bifidobacterium", "Escherichia", "Enterococcus", "Klebsiella", "Streptococcus", "Bacteroides")
genus_traj <- genus_complete %>%
  filter(Genus %in% trajectory_genera) %>%
  group_by(Timepoint, Genus) %>%
  summarise(
    n = n_distinct(SampleID),
    median = median(RelAbundance, na.rm = TRUE),
    q25 = quantile(RelAbundance, 0.25, na.rm = TRUE),
    q75 = quantile(RelAbundance, 0.75, na.rm = TRUE),
    prevalence = mean(RelAbundance > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Genus = factor(Genus, levels = trajectory_genera),
    TimeIndex = as.numeric(Timepoint)
  )

bifido_samples <- genus_complete %>%
  filter(Genus == "Bifidobacterium") %>%
  mutate(TimeIndex = as.numeric(Timepoint))

bifido_summary <- genus_traj %>%
  filter(Genus == "Bifidobacterium")

species_keep <- c(
  "Bifidobacterium longum",
  "Bifidobacterium bifidum",
  "Bifidobacterium breve",
  "Bifidobacterium catenulatum",
  "Bifidobacterium pseudocatenulatum",
  "Bifidobacterium animalis"
)

species_panel <- species_summary %>%
  filter(Species %in% species_keep) %>%
  mutate(
    Species = factor(Species, levels = species_keep),
    TimeIndex = as.numeric(Timepoint),
    SpeciesShort = str_remove(as.character(Species), "^Bifidobacterium ")
  )

species_labels <- species_panel %>%
  filter(Timepoint == "D90") %>%
  mutate(LabelY = MeanPercent) %>%
  arrange(LabelY) %>%
  mutate(LabelY = if_else(row_number() <= 3 & LabelY < 4, LabelY + (row_number() - 1) * 1.1, LabelY))

source_data <- list(
  sample_counts = sample_counts,
  genus_sample_landscape = heat_df,
  pcoa = pcoa_df,
  pcoa_centroids = centroid_df,
  permanova = tibble(R2 = adon$R2[1], P = adon$`Pr(>F)`[1], F = adon$F[1]),
  genus_temporal_effect = genus_effect,
  genus_trajectories = genus_traj,
  bifidobacterium_samples = bifido_samples,
  bifidobacterium_species_summary = species_panel
)

for (nm in names(source_data)) {
  write_csv(source_data[[nm]], file.path(out_dir, paste0("Fig3_source_", nm, ".csv")))
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
  labs(title = "a  Genus-level sample landscape", x = NULL, y = NULL) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 7.2, hjust = 0, margin = margin(0, 0, 2, 0)),
    strip.text = element_blank(),
    panel.spacing.x = unit(0.8, "mm"),
    plot.margin = margin(1, 4, 0, 4)
  )

pA_heat <- ggplot(heat_df, aes(SamplePlot, GenusPlot, fill = RelAbundance)) +
  geom_tile(width = 0.95, height = 0.92) +
  facet_grid(. ~ Timepoint, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#D6E6F5", "#7AA6D8", "#245F9E", "#08306B"),
    values = rescale(c(0, 5, 20, 55, 100)),
    limits = c(0, 100),
    oob = squish,
    name = "Relative\nabundance"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 6.0) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 5.5, colour = "#222222"),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.8, "mm"),
    strip.text = element_text(face = "bold", size = 5.8),
    legend.position = "right",
    legend.key.height = unit(12, "mm"),
    plot.margin = margin(0, 4, 2, 4)
  )

pA <- pA_strip / pA_heat + plot_layout(heights = c(0.09, 1))

pB <- ggplot(pcoa_df, aes(PCoA1, PCoA2, colour = Timepoint)) +
  stat_ellipse(aes(group = Timepoint), type = "norm", linewidth = 0.34, alpha = 0.6, show.legend = FALSE) +
  geom_point(size = 1.55, alpha = 0.76) +
  geom_path(
    data = centroid_df %>% arrange(TimeIndex),
    aes(PCoA1, PCoA2),
    inherit.aes = FALSE,
    linewidth = 0.42,
    colour = "#222222",
    arrow = arrow(length = unit(1.4, "mm"), type = "closed")
  ) +
  geom_point(data = centroid_df, aes(PCoA1, PCoA2), inherit.aes = FALSE, size = 2.2, shape = 21, fill = "white", colour = "#222222", stroke = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = adon_label, hjust = -0.02, vjust = 1.25, size = 1.85, family = "Arial", colour = "#4A4A4A") +
  scale_colour_manual(values = pal_time, name = NULL) +
  labs(
    title = "b  Bray-Curtis community shift",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.0))) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.box.margin = margin(-6, 0, -4, 0),
    panel.grid.major = element_line(linewidth = 0.15, colour = "#E9E9E9")
  )

pC <- ggplot(effect_df, aes(delta_d30_d05, Genus)) +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "#777777") +
  geom_segment(aes(x = 0, xend = delta_d30_d05, y = Genus, yend = Genus, colour = Direction), linewidth = 0.5) +
  geom_point(aes(size = prevalence * 100, shape = FDRClass, colour = Direction), stroke = 0.42) +
  scale_colour_manual(values = c("Higher at D30" = "#2A6FBB", "Lower at D30" = "#C84D4D", "No median shift" = "#9A9A9A"), name = NULL) +
  scale_shape_manual(values = c("FDR < 0.05" = 21, "P < 0.05" = 24, "NS" = 21), name = NULL) +
  scale_size_continuous(range = c(1.1, 3.1), breaks = c(25, 50, 75, 100), name = "Prevalence (%)") +
  labs(
    title = "c  Genus temporal effect",
    x = "Median change (D30 - D05)",
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 5.4),
    panel.grid.major.y = element_line(linewidth = 0.12, colour = "#EFEFEF"),
    legend.position = "right"
  )

pD <- ggplot(genus_traj, aes(TimeIndex, median, colour = Genus, group = Genus)) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = Genus), colour = NA, alpha = 0.13, show.legend = FALSE) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.65) +
  facet_wrap(~ Genus, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal_key, guide = "none") +
  scale_fill_manual(values = pal_key, guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = time_levels, expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "d  Key genus trajectories",
    x = NULL,
    y = "Median relative abundance"
  ) +
  theme(
    axis.text.x = element_text(size = 5.2),
    strip.text = element_text(size = 5.4),
    panel.grid.major.y = element_line(linewidth = 0.15, colour = "#E9E9E9")
  )

pE <- ggplot(bifido_samples, aes(TimeIndex, RelAbundance)) +
  geom_jitter(aes(colour = Timepoint), width = 0.08, height = 0, size = 1.15, alpha = 0.50, show.legend = FALSE) +
  geom_ribbon(data = bifido_summary, aes(x = TimeIndex, ymin = q25, ymax = q75), inherit.aes = FALSE, fill = "#2A6FBB", alpha = 0.13) +
  geom_line(data = bifido_summary, aes(TimeIndex, median), inherit.aes = FALSE, linewidth = 0.68, colour = "#245F9E") +
  geom_point(data = bifido_summary, aes(TimeIndex, median), inherit.aes = FALSE, size = 2.0, colour = "#245F9E") +
  geom_text(
    data = sample_counts %>% mutate(TimeIndex = as.numeric(Timepoint), y = 99),
    aes(TimeIndex, y, label = paste0("n=", n_samples)),
    inherit.aes = FALSE,
    size = 1.65,
    family = "Arial",
    colour = "#4A4A4A"
  ) +
  scale_colour_manual(values = pal_time) +
  scale_x_continuous(breaks = 1:4, labels = time_levels, expand = expansion(mult = c(0.04, 0.05))) +
  scale_y_continuous(limits = c(0, 104), labels = label_percent(scale = 1), expand = c(0, 0)) +
  labs(
    title = "e  Bifidobacterium dominance",
    x = NULL,
    y = "Relative abundance"
  ) +
  theme(panel.grid.major.y = element_line(linewidth = 0.15, colour = "#E9E9E9"))

pF <- ggplot(species_panel, aes(TimeIndex, MeanPercent, colour = Species, group = Species)) +
  geom_line(linewidth = 0.52) +
  geom_point(aes(size = PrevalencePercent), alpha = 0.92) +
  geom_text(
    data = species_labels,
    aes(x = 4.18, y = LabelY, label = SpeciesShort, colour = Species),
    inherit.aes = FALSE,
    hjust = 0,
    size = 1.65,
    family = "Arial",
    show.legend = FALSE
  ) +
  scale_colour_manual(values = pal_species, name = NULL) +
  scale_size_continuous(range = c(1.0, 3.4), guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = time_levels, limits = c(1, 5.25), expand = expansion(mult = c(0.04, 0.02))) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0.02, 0.07))) +
  labs(
    title = "f  Species-level refinement",
    x = NULL,
    y = "Mean relative abundance"
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(linewidth = 0.15, colour = "#E9E9E9")
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
  plot_layout(heights = c(1.05, 0.85, 0.72, 0.95)) &
  theme(plot.margin = margin(3, 4, 3, 4))

base <- file.path(out_dir, "Fig3_infant_microbiome_maturation_redrawn")
w <- 183 / 25.4
h <- 260 / 25.4

ggsave(paste0(base, ".png"), fig, width = w, height = h, dpi = 600, bg = "white")
ggsave(paste0(base, ".pdf"), fig, width = w, height = h, device = cairo_pdf, bg = "white")
svglite::svglite(paste0(base, ".svg"), width = w, height = h, bg = "white")
print(fig)
dev.off()
ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = 600, background = "white")
print(fig)
dev.off()

message("Exported Fig. 3 to: ", out_dir)
