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
out_dir <- file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig7_clinical_modifier_layer")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(root, "p1_clinical_modifiers", "tables", "p1_scfa_maturation_clinical_analysis_matrix.csv")
summary_path <- file.path(root, "p1_clinical_modifiers", "tables", "p1_scfa_maturation_by_clinical_group_summary.csv")
model_path <- file.path(root, "p1_clinical_modifiers", "tables", "p1_scfa_maturation_clinical_lme_results.csv")
coverage_path <- file.path(root, "p1_clinical_modifiers", "tables", "p1_clinical_coverage_by_time.csv")

time_levels <- c("D05", "D14", "D30", "D90")
time_labels <- c(D05 = "D05", D14 = "D14", D30 = "D30", D90 = "D90")
pal_bmi <- c(Lower = "#6F8FB4", Higher = "#C18470")
pal_feeding <- c("Exclusive BF" = "#6F8FB4", "Non-exclusive BF" = "#C18470")
pal_delivery <- c(Vaginal = "#7BAA8F", Cesarean = "#B28E6B")
pal_bmi_pred <- c("BMI z P10" = "#6F8FB4", "BMI z median" = "#777777", "BMI z P90" = "#C18470")
pal_model <- c("P < 0.05" = "#B84A56", "P >= 0.05" = "#5E5E5E")

theme_pub <- function(base_size = 6.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks = element_line(linewidth = 0.28, colour = "black"),
      axis.text = element_text(colour = "#222222"),
      plot.title = element_text(size = base_size + 1.1, face = "bold", hjust = 0),
      legend.title = element_text(size = base_size - 0.1),
      legend.text = element_text(size = base_size - 0.35),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(3, 4, 3, 4)
    )
}

dat <- read_csv(matrix_path, show_col_types = FALSE) %>%
  mutate(
    TimepointCode = factor(TimepointCode, levels = time_levels),
    MaternalBMIGroup = factor(MaternalBMIGroup, levels = c("Lower", "Higher")),
    FeedingMode = factor(FeedingMode, levels = c("Exclusive breastfeeding", "Non-exclusive breastfeeding")),
    DeliveryMode = factor(DeliveryMode, levels = c("Vaginal", "Cesarean"))
  )

group_summary <- read_csv(summary_path, show_col_types = FALSE) %>%
  mutate(TimepointCode = factor(TimepointCode, levels = time_levels))

bmi_summary <- group_summary %>%
  filter(ClinicalVariable == "MaternalBMIGroup") %>%
  mutate(Group = factor(Group, levels = c("Lower", "Higher")))

clinical_summary <- group_summary %>%
  filter(ClinicalVariable %in% c("FeedingMode", "DeliveryMode")) %>%
  mutate(
    Panel = recode(ClinicalVariable, FeedingMode = "Feeding", DeliveryMode = "Delivery"),
    Group = case_when(
      ClinicalVariable == "FeedingMode" ~ recode(Group, "Exclusive breastfeeding" = "Exclusive BF", "Non-exclusive breastfeeding" = "Non-exclusive BF"),
      ClinicalVariable == "DeliveryMode" ~ factor(Group, levels = c("Vaginal", "Cesarean")) %>% as.character(),
      TRUE ~ Group
    ),
    Group = factor(Group, levels = c("Exclusive BF", "Non-exclusive BF", "Vaginal", "Cesarean")),
    Panel = factor(Panel, levels = c("Feeding", "Delivery"))
  )

model_raw <- read_csv(model_path, show_col_types = FALSE)

time_lookup <- dat %>%
  group_by(TimepointCode) %>%
  summarise(Day = median(Day, na.rm = TRUE), DayScaled = median(DayScaled, na.rm = TRUE), .groups = "drop") %>%
  mutate(TimepointCode = factor(TimepointCode, levels = time_levels))

m7_coef <- model_raw %>%
  filter(ModelID == "M7_bmi_interaction") %>%
  select(term, estimate)
m7_beta <- setNames(m7_coef$estimate, m7_coef$term)
bmi_quantiles <- quantile(dat$MaternalBMI_z, probs = c(0.10, 0.50, 0.90), na.rm = TRUE, names = FALSE)

bmi_prediction <- tidyr::expand_grid(
  time_lookup,
  BMI_Level = factor(names(pal_bmi_pred), levels = names(pal_bmi_pred))
) %>%
  mutate(
    MaternalBMI_z = recode(
      as.character(BMI_Level),
      `BMI z P10` = bmi_quantiles[1],
      `BMI z median` = bmi_quantiles[2],
      `BMI z P90` = bmi_quantiles[3]
    ),
    PredictedScore =
      m7_beta["(Intercept)"] +
      m7_beta["DayScaled"] * DayScaled +
      m7_beta["MaternalBMI_z"] * MaternalBMI_z +
      m7_beta["DayScaled:MaternalBMI_z"] * DayScaled * MaternalBMI_z
  )

make_contrast <- function(data, variable, group1, group2, label) {
  rows <- lapply(time_levels, function(tp) {
    d <- data %>%
      filter(TimepointCode == tp, !is.na(.data[[variable]]), .data[[variable]] %in% c(group1, group2))
    x <- d %>% filter(.data[[variable]] == group1) %>% pull(SCFA_MaturationScore)
    y <- d %>% filter(.data[[variable]] == group2) %>% pull(SCFA_MaturationScore)
    n1 <- length(x)
    n2 <- length(y)
    med1 <- if (n1 > 0) median(x, na.rm = TRUE) else NA_real_
    med2 <- if (n2 > 0) median(y, na.rm = TRUE) else NA_real_
    p <- if (n1 > 0 && n2 > 0) suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value) else NA_real_
    tibble(
      Contrast = label,
      TimepointCode = tp,
      Group1 = group1,
      Group2 = group2,
      N1 = n1,
      N2 = n2,
      Median1 = med1,
      Median2 = med2,
      Difference = med2 - med1,
      PValue = p
    )
  })
  bind_rows(rows)
}

contrast_summary <- bind_rows(
  make_contrast(dat, "MaternalBMIGroup", "Lower", "Higher", "BMI: higher - lower"),
  make_contrast(dat, "FeedingMode", "Exclusive breastfeeding", "Non-exclusive breastfeeding", "Feeding: non-exclusive - exclusive"),
  make_contrast(dat, "DeliveryMode", "Vaginal", "Cesarean", "Delivery: cesarean - vaginal"),
  make_contrast(dat, "BirthWeightGroup", "Lower", "Higher", "Birth weight: higher - lower")
) %>%
  mutate(
    TimepointCode = factor(TimepointCode, levels = time_levels),
    Contrast = factor(Contrast, levels = rev(c(
      "BMI: higher - lower",
      "Feeding: non-exclusive - exclusive",
      "Delivery: cesarean - vaginal",
      "Birth weight: higher - lower"
    ))),
    TileLabel = if_else(
      is.na(PValue),
      "NA",
      paste0(sprintf("%.2f", Difference), "\nP=", sprintf("%.2f", PValue), "\n", N1, "/", N2)
    )
  )

forest_x_min <- -0.95
forest_x_max <- 1.25

forest_terms <- model_raw %>%
  filter(
    (ModelID == "M1_age_only" & term == "DayScaled") |
      (ModelID == "M7_bmi_interaction" & term %in% c("MaternalBMI_z", "DayScaled:MaternalBMI_z")) |
      (ModelID == "M3_feeding_interaction" & term %in% c("FeedingModeNon-exclusive breastfeeding", "DayScaled:FeedingModeNon-exclusive breastfeeding")) |
      (ModelID == "M5_delivery_interaction" & term %in% c("DeliveryModeCesarean", "DayScaled:DeliveryModeCesarean")) |
      (ModelID == "M9_birthweight_interaction" & term %in% c("BirthWeightKg_z", "DayScaled:BirthWeightKg_z")) |
      (ModelID == "M10_minimal_multivariable" & term %in% c(
        "FeedingModeNon-exclusive breastfeeding",
        "DeliveryModeCesarean",
        "MaternalBMI_z",
        "BirthWeightKg_z"
      ))
  ) %>%
  mutate(
    TermLabel = case_when(
      ModelID == "M1_age_only" & term == "DayScaled" ~ "Age",
      ModelID == "M7_bmi_interaction" & term == "DayScaled:MaternalBMI_z" ~ "Age x maternal BMI",
      ModelID == "M7_bmi_interaction" & term == "MaternalBMI_z" ~ "Maternal BMI",
      ModelID == "M3_feeding_interaction" & term == "DayScaled:FeedingModeNon-exclusive breastfeeding" ~ "Age x feeding",
      ModelID == "M3_feeding_interaction" & term == "FeedingModeNon-exclusive breastfeeding" ~ "Non-exclusive feeding",
      ModelID == "M5_delivery_interaction" & term == "DayScaled:DeliveryModeCesarean" ~ "Age x cesarean",
      ModelID == "M5_delivery_interaction" & term == "DeliveryModeCesarean" ~ "Cesarean delivery",
      ModelID == "M9_birthweight_interaction" & term == "DayScaled:BirthWeightKg_z" ~ "Age x birth weight",
      ModelID == "M9_birthweight_interaction" & term == "BirthWeightKg_z" ~ "Birth weight",
      ModelID == "M10_minimal_multivariable" & term == "MaternalBMI_z" ~ "Maternal BMI, minimal",
      ModelID == "M10_minimal_multivariable" & term == "FeedingModeNon-exclusive breastfeeding" ~ "Feeding, minimal",
      ModelID == "M10_minimal_multivariable" & term == "DeliveryModeCesarean" ~ "Delivery, minimal",
      ModelID == "M10_minimal_multivariable" & term == "BirthWeightKg_z" ~ "Birth weight, minimal",
      TRUE ~ term
    ),
    Layer = case_when(
      ModelID == "M1_age_only" ~ "Age baseline",
      ModelID == "M7_bmi_interaction" ~ "BMI model",
      ModelID == "M3_feeding_interaction" ~ "Feeding model",
      ModelID == "M5_delivery_interaction" ~ "Delivery model",
      ModelID == "M9_birthweight_interaction" ~ "Birth weight model",
      ModelID == "M10_minimal_multivariable" ~ "Minimal model",
      TRUE ~ InterpretationLayer
    ),
    Sig = if_else(PValue < 0.05, "P < 0.05", "P >= 0.05"),
    TermLabel = factor(TermLabel, levels = rev(c(
      "Age",
      "Age x maternal BMI", "Maternal BMI",
      "Age x feeding", "Non-exclusive feeding",
      "Age x cesarean", "Cesarean delivery",
      "Age x birth weight", "Birth weight",
      "Maternal BMI, minimal", "Feeding, minimal", "Delivery, minimal", "Birth weight, minimal"
    ))),
    Layer = factor(Layer, levels = c("Age baseline", "BMI model", "Feeding model", "Delivery model", "Birth weight model", "Minimal model")),
    Sig = factor(Sig, levels = c("P < 0.05", "P >= 0.05")),
    PlotLow = pmax(conf.low, forest_x_min),
    PlotHigh = pmin(conf.high, forest_x_max),
    WideCI = conf.low < forest_x_min | conf.high > forest_x_max,
    Text = paste0(
      "P=", sprintf("%.3g", PValue), " | n=", N, ", d=", N_Dyads,
      if_else(WideCI, " | wide CI", "")
    )
  )

coverage <- read_csv(coverage_path, show_col_types = FALSE) %>%
  mutate(TimepointCode = factor(TimepointCode, levels = time_levels)) %>%
  rename(
    samples = N_samples,
    dyads = N_dyads,
    delivery = N_delivery,
    feeding = N_feeding,
    bmi = N_bmi,
    birthweight = N_birthweight
  )

coverage_long <- coverage %>%
  select(TimepointCode, samples, delivery, feeding, bmi, birthweight) %>%
  pivot_longer(-TimepointCode, names_to = "Layer", values_to = "N") %>%
  mutate(
    Layer = factor(Layer, levels = c("samples", "bmi", "feeding", "delivery", "birthweight")),
    LayerLabel = recode(Layer, samples = "SCFA", bmi = "BMI", feeding = "Feeding", delivery = "Delivery", birthweight = "Birth weight"),
    LayerLabel = factor(LayerLabel, levels = c("SCFA", "BMI", "Feeding", "Delivery", "Birth weight"))
  )

write_csv(dat, file.path(out_dir, "Fig7_source_clinical_analysis_matrix.csv"))
write_csv(bmi_summary, file.path(out_dir, "Fig7_source_bmi_group_summary.csv"))
write_csv(bmi_prediction, file.path(out_dir, "Fig7_source_bmi_model_prediction.csv"))
write_csv(contrast_summary, file.path(out_dir, "Fig7_source_timepoint_contrasts.csv"))
write_csv(clinical_summary, file.path(out_dir, "Fig7_source_feeding_delivery_summary.csv"))
write_csv(forest_terms, file.path(out_dir, "Fig7_source_model_forest_terms.csv"))
write_csv(coverage_long, file.path(out_dir, "Fig7_source_coverage_by_time.csv"))

pA <- ggplot() +
  geom_point(
    data = filter(dat, !is.na(MaternalBMI_z)),
    aes(Day, SCFA_MaturationScore, fill = MaternalBMI_z),
    shape = 21,
    size = 1.35,
    stroke = 0.25,
    colour = "#333333",
    alpha = 0.74
  ) +
  geom_line(
    data = bmi_prediction,
    aes(Day, PredictedScore, group = BMI_Level, colour = BMI_Level),
    linewidth = 0.68
  ) +
  geom_point(
    data = bmi_prediction,
    aes(Day, PredictedScore, colour = BMI_Level),
    size = 1.35
  ) +
  scale_fill_gradient2(
    low = "#6F8FB4",
    mid = "#F4F4F4",
    high = "#C18470",
    midpoint = 0,
    breaks = c(-0.5, 0, 1, 2, 3),
    name = "BMI z"
  ) +
  scale_colour_manual(values = pal_bmi_pred, name = NULL) +
  scale_x_continuous(breaks = c(5, 14, 30, 90), labels = time_labels) +
  scale_y_continuous(breaks = c(-1, 0, 1), expand = expansion(mult = c(0.04, 0.08))) +
  coord_cartesian(ylim = c(-1.65, 1.25)) +
  labs(
    title = "a  BMI x age model",
    x = "Infant age",
    y = "SCFA maturation score"
  ) +
  theme_pub(5.9) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(2.8, "mm"),
    panel.grid.major = element_line(linewidth = 0.14, colour = "#E8E8E8")
  )

pB <- ggplot(contrast_summary, aes(TimepointCode, Contrast, fill = Difference)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = TileLabel), size = 1.65, family = "Arial", lineheight = 0.88, colour = "#222222") +
  scale_fill_gradient2(low = "#6F8FB4", mid = "#F5F5F5", high = "#C18470", midpoint = 0, limits = c(-1.1, 1.1), name = "Median diff.") +
  labs(
    title = "b  Timepoint contrasts",
    x = NULL,
    y = NULL
  ) +
  theme_pub(5.8) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.text.y = element_text(size = 5.1),
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(2.7, "mm"),
    panel.grid = element_blank()
  )

pC <- ggplot(forest_terms, aes(estimate, TermLabel)) +
  geom_vline(xintercept = 0, linewidth = 0.30, linetype = "dashed", colour = "#8A8A8A") +
  geom_errorbarh(aes(xmin = PlotLow, xmax = PlotHigh, colour = Sig), height = 0, linewidth = 0.46) +
  geom_point(
    data = filter(forest_terms, conf.low < forest_x_min),
    aes(x = forest_x_min, y = TermLabel, colour = Sig),
    inherit.aes = FALSE,
    shape = 60,
    size = 1.65
  ) +
  geom_point(
    data = filter(forest_terms, conf.high > forest_x_max),
    aes(x = forest_x_max, y = TermLabel, colour = Sig),
    inherit.aes = FALSE,
    shape = 62,
    size = 1.65
  ) +
  geom_point(aes(fill = Sig), shape = 21, size = 1.85, stroke = 0.25, colour = "#333333") +
  geom_text(aes(x = 1.36, label = Text), hjust = 0, size = 1.40, family = "Arial", colour = "#4F4F4F") +
  facet_grid(Layer ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = pal_model, name = NULL) +
  scale_fill_manual(values = pal_model, name = NULL) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1.0), expand = expansion(mult = c(0.03, 0.25))) +
  coord_cartesian(xlim = c(forest_x_min, 1.85), clip = "off") +
  labs(
    title = "c  Mixed-model terms",
    x = "Coefficient",
    y = NULL
  ) +
  theme_pub(5.7) +
  theme(
    axis.text.y = element_text(size = 5.1),
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(2.8, "mm"),
    panel.grid.major.x = element_line(linewidth = 0.14, colour = "#E8E8E8"),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 5.4, face = "bold")
  )

pD <- ggplot(coverage_long, aes(TimepointCode, LayerLabel, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = N), size = 1.95, family = "Arial", colour = "#222222") +
  scale_fill_gradient(low = "#F4F4F4", high = "#6F8FB4", name = "n", breaks = c(10, 20)) +
  labs(
    title = "d  Clinical coverage",
    x = NULL,
    y = NULL
  ) +
  theme_pub(5.7) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "top",
    legend.justification = "left",
    legend.key.size = unit(2.8, "mm"),
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5.3),
    panel.grid = element_blank()
  )

fig <- (pA | pC) / (pB | pD) +
  plot_layout(widths = c(1.0, 1.08), heights = c(1.05, 0.82)) &
  theme(plot.margin = margin(3, 4, 3, 4))

base_file <- file.path(out_dir, "Fig7_clinical_modifier_layer_redrawn")
ggsave(paste0(base_file, ".png"), fig, width = 183, height = 145, units = "mm", dpi = 600, bg = "white")
ggsave(paste0(base_file, ".pdf"), fig, width = 183, height = 145, units = "mm", device = cairo_pdf, bg = "white")
svglite::svglite(paste0(base_file, ".svg"), width = 183 / 25.4, height = 145 / 25.4, bg = "white")
print(fig)
dev.off()
ragg::agg_tiff(paste0(base_file, ".tiff"), width = 183 / 25.4, height = 145 / 25.4, units = "in", res = 600, background = "white")
print(fig)
dev.off()

message("Wrote Fig7 redraw outputs to: ", out_dir)
