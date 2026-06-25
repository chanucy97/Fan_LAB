suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pROC)
})

option_list <- list(
  make_option("--predictions", type = "character"),
  make_option("--out", type = "character"),
  make_option("--label-col", type = "character", default = "label"),
  make_option("--score-col", type = "character", default = "score"),
  make_option("--model-col", type = "character", default = "model")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$predictions) || is.null(opt$out)) stop("--predictions and --out are required")

df <- read_csv(opt$predictions, show_col_types = FALSE)
if (!(opt$model_col %in% names(df))) df[[opt$model_col]] <- "model"

roc_df <- bind_rows(lapply(split(df, df[[opt$model_col]]), function(part) {
  roc_obj <- pROC::roc(part[[opt$label_col]], part[[opt$score_col]], quiet = TRUE, direction = "<", levels = c(0, 1))
  data.frame(
    model = unique(part[[opt$model_col]])[1],
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    auc = as.numeric(pROC::auc(roc_obj))
  )
}))
labels <- roc_df %>% group_by(model) %>% summarize(auc = first(auc), .groups = "drop") %>% mutate(label = sprintf("%s (AUC %.3f)", model, auc))
roc_df <- roc_df %>% left_join(labels, by = "model")

p <- ggplot(roc_df, aes(x = fpr, y = tpr, color = label)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.8) +
  coord_equal() +
  theme_classic(base_size = 8) +
  labs(x = "False-positive rate", y = "True-positive rate", color = NULL)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(opt$out, ".png"), p, width = 90, height = 75, units = "mm", dpi = 600, bg = "white")
ggsave(paste0(opt$out, ".pdf"), p, width = 90, height = 75, units = "mm")
