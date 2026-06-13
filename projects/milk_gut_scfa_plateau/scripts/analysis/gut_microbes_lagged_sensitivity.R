suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

root <- normalizePath(Sys.getenv("FANLAB_RESEARCH2_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "final_manuscript_planning_20260605", "gut_microbes_revision")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lag_matrix_path <- file.path(root, "direction3_multiomics_blueprint", "tables", "direction3_analysis_matrix_lagged_long.csv")
formal_path <- file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig5_milk_scfa_lagged_candidates", "Fig5_source_formal_model_candidates.csv")

lag_matrix <- read_csv(lag_matrix_path, show_col_types = FALSE)
formal <- read_csv(formal_path, show_col_types = FALSE)

scale_numeric <- function(x) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 3 || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

fit_basic <- function(df) {
  df <- df %>%
    mutate(
      ExposureLog = log10(ExposureValue + 1),
      OutcomeLog = log10(OutcomeValue + 1),
      ExposureZ = scale_numeric(ExposureLog),
      OutcomeZ = scale_numeric(OutcomeLog)
    ) %>%
    filter(is.finite(ExposureZ), is.finite(OutcomeZ))

  if (nrow(df) < 8 || length(unique(df$DyadID)) < 8) {
    return(tibble(
      N = nrow(df), Dyads = length(unique(df$DyadID)), Beta = NA_real_,
      SE = NA_real_, Statistic = NA_real_, PValue = NA_real_, Error = "Insufficient N"
    ))
  }

  tryCatch({
    fit <- lm(OutcomeZ ~ ExposureZ, data = df)
    s <- summary(fit)$coefficients
    tibble(
      N = nrow(df),
      Dyads = length(unique(df$DyadID)),
      Beta = unname(s["ExposureZ", "Estimate"]),
      SE = unname(s["ExposureZ", "Std. Error"]),
      Statistic = unname(s["ExposureZ", "t value"]),
      PValue = unname(s["ExposureZ", "Pr(>|t|)"]),
      Error = NA_character_
    )
  }, error = function(e) {
    tibble(
      N = nrow(df), Dyads = length(unique(df$DyadID)), Beta = NA_real_,
      SE = NA_real_, Statistic = NA_real_, PValue = NA_real_, Error = conditionMessage(e)
    )
  })
}

candidate_hits <- formal %>%
  filter(
    Screen == "lagged",
    OutcomeType == "SCFA",
    Model == "lagged_lm_basic",
    ModelOk,
    !is.na(PValueModel),
    PValueModel < 0.05
  ) %>%
  arrange(PValueModel) %>%
  select(HitID, Contrast, ExposureLabel, OutcomeLabel, ExposureClean, OutcomeClean, NModel, Dyads, Beta, PValueModel, FDRScreen) %>%
  slice_head(n = 14)

loo_details <- bind_rows(lapply(seq_len(nrow(candidate_hits)), function(i) {
  hit <- candidate_hits[i, ]
  df <- lag_matrix %>% filter(HitID == hit$HitID)
  full_fit <- fit_basic(df) %>% mutate(DroppedDyadID = NA_integer_, FitType = "full")
  loo <- bind_rows(lapply(sort(unique(df$DyadID)), function(d) {
    fit_basic(df %>% filter(DyadID != d)) %>%
      mutate(DroppedDyadID = d, FitType = "leave_one_dyad_out")
  }))
  bind_rows(full_fit, loo) %>%
    mutate(
      HitID = hit$HitID,
      Contrast = hit$Contrast,
      ExposureLabel = hit$ExposureLabel,
      OutcomeLabel = hit$OutcomeLabel,
      ExposureClean = hit$ExposureClean,
      OutcomeClean = hit$OutcomeClean,
      FullBetaOriginal = hit$Beta,
      FullPOriginal = hit$PValueModel,
      FDRScreen = hit$FDRScreen
    )
}))

loo_summary <- loo_details %>%
  filter(FitType == "leave_one_dyad_out") %>%
  group_by(HitID, Contrast, ExposureLabel, OutcomeLabel, ExposureClean, OutcomeClean, FullBetaOriginal, FullPOriginal, FDRScreen) %>%
  summarise(
    LeaveOneOutFits = n(),
    DirectionStableN = sum(sign(Beta) == sign(FullBetaOriginal), na.rm = TRUE),
    DirectionStableFraction = DirectionStableN / LeaveOneOutFits,
    NominalP05N = sum(PValue < 0.05, na.rm = TRUE),
    NominalP10N = sum(PValue < 0.10, na.rm = TRUE),
    MinBeta = min(Beta, na.rm = TRUE),
    MaxBeta = max(Beta, na.rm = TRUE),
    MedianBeta = median(Beta, na.rm = TRUE),
    MaxPValue = max(PValue, na.rm = TRUE),
    MedianPValue = median(PValue, na.rm = TRUE),
    RobustnessClass = case_when(
      DirectionStableFraction == 1 & NominalP05N / LeaveOneOutFits >= 0.75 ~ "direction_and_nominal_signal_stable",
      DirectionStableFraction == 1 & NominalP10N / LeaveOneOutFits >= 0.75 ~ "direction_stable_nominal_signal_sensitive",
      DirectionStableFraction >= 0.90 ~ "direction_mostly_stable_p_sensitive",
      TRUE ~ "sensitive_to_individual_dyads"
    ),
    .groups = "drop"
  ) %>%
  arrange(FullPOriginal)

write_csv(loo_details, file.path(out_dir, "gut_microbes_lagged_leave_one_dyad_out_details.csv"))
write_csv(loo_summary, file.path(out_dir, "gut_microbes_lagged_leave_one_dyad_out_summary.csv"))

writeLines(
  c(
    "# Gut Microbes lagged LTF/HMO-SCFA sensitivity analysis",
    "",
    "Sensitivity test: for each nominally significant lagged basic SCFA model, the model `OutcomeZ ~ ExposureZ` was refitted after dropping one mother-infant dyad at a time.",
    "Exposure and outcome variables were transformed as `log10(x + 1)` and then standardized within each fitted subset, matching the original formal-model workflow.",
    "This analysis tests whether the candidate direction is dominated by a single dyad. It does not convert exploratory associations into confirmatory causal evidence.",
    "",
    paste0("Models tested: ", nrow(candidate_hits)),
    paste0("Output summary: ", file.path(out_dir, "gut_microbes_lagged_leave_one_dyad_out_summary.csv")),
    paste0("Output details: ", file.path(out_dir, "gut_microbes_lagged_leave_one_dyad_out_details.csv"))
  ),
  con = file.path(out_dir, "gut_microbes_lagged_sensitivity_notes.md")
)

message("Wrote Gut Microbes sensitivity outputs to: ", out_dir)
