################################################################################
## Public repository runtime root
################################################################################
## This public script keeps code only. Put controlled data outside Git and set:
##   Sys.setenv(HCC_DLC1_ROOT = "/path/to/private/hcc_DLC1")
## If the variable is unset, the current working directory is used.
hcc_dlc1_root <- function() {
  root <- Sys.getenv("HCC_DLC1_ROOT", unset = NA_character_)
  if (is.na(root) || !nzchar(root)) root <- getwd()
  normalizePath(root, winslash = "/", mustWork = FALSE)
}
################################################################################
##泛癌部分
################################################################################
options(stringsAsFactors = FALSE)
options(timeout = max(3600, getOption("timeout")))
set.seed(20260615)

gene_symbol <- "DLC1"
ensembl_gene_id <- "ENSG00000164741"

out_dir_name <- "DLC1_panCancer_submission_rawdata_outputs"

run_tcga_gtex_download <- TRUE
run_survival <- TRUE
run_tmb_from_mutation <- TRUE
run_msi_if_available <- TRUE
run_immune_deconvolution <- FALSE
run_checkpoint_correlation_from_raw_expression <- TRUE

script_args <- commandArgs(trailingOnly = TRUE)
analysis_mode <- "full"
selected_immune_method <- NA_character_
selected_cancer <- NA_character_

parse_cli_args <- function(args) {
  out <- list(mode = "full", method = NA_character_, cancer = NA_character_)
  if (length(args) == 0) return(out)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--mode", "-m")) {
      i <- i + 1
      out$mode <- args[[i]]
    } else if (key %in% c("--method")) {
      i <- i + 1
      out$method <- args[[i]]
    } else if (key %in% c("--cancer")) {
      i <- i + 1
      out$cancer <- toupper(args[[i]])
    } else if (!grepl("^--", key) && identical(out$mode, "full")) {
      out$mode <- key
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    i <- i + 1
  }
  out
}

cli <- parse_cli_args(script_args)
analysis_mode <- cli$mode
selected_immune_method <- cli$method
selected_cancer <- cli$cancer

min_group_n <- 3
min_survival_n <- 20
min_survival_events <- 3
km_p_cutoff_for_curve <- 0.05
cor_method <- "spearman"
exome_size_mb <- 38
tmb_count_mode <- "pass_all_snv"
immune_expression_scale <- "tpm_from_fpkm"

cancer_order <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

checkpoint_genes <- c("CD274", "CTLA4", "HAVCR2", "LAG3", "PDCD1",
                      "PDCD1LG2", "TIGIT", "SIGLEC15")

immune_methods <- c("timer", "xcell", "cibersort", "epic",
                    "mcp_counter", "quantiseq")

gtex_tissue_to_cancer <- list(
  ACC = c("Adrenal Gland"),
  BLCA = c("Bladder"),
  BRCA = c("Breast"),
  CESC = c("Cervix Uteri", "Cervix"),
  COAD = c("Colon"),
  DLBC = c("Cells - EBV-transformed lymphocytes", "Whole Blood"),
  ESCA = c("Esophagus"),
  GBM = c("Brain"),
  KICH = c("Kidney"),
  KIRC = c("Kidney"),
  KIRP = c("Kidney"),
  LGG = c("Brain"),
  LIHC = c("Liver"),
  LUAD = c("Lung"),
  LUSC = c("Lung"),
  OV = c("Ovary"),
  PAAD = c("Pancreas"),
  PRAD = c("Prostate"),
  READ = c("Colon"),
  SKCM = c("Skin"),
  STAD = c("Stomach"),
  TGCT = c("Testis"),
  THCA = c("Thyroid"),
  UCEC = c("Uterus"),
  UCS = c("Uterus")
)

toil_expression_url <- paste0(
  "https://toil.xenahubs.net/download/",
  "TcgaTargetGtex_rsem_gene_tpm.gz"
)
toil_phenotype_url <- paste0(
  "https://toil.xenahubs.net/download/",
  "TcgaTargetGTEX_phenotype.txt.gz"
)

################################################################################
# Utilities
################################################################################

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Please install it first.", pkg),
         call. = FALSE)
  }
}

maybe_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

install_if_missing <- function(pkg, installer = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message("Already installed: ", pkg)
    return(invisible(TRUE))
  }
  message("Installing: ", pkg)
  if (is.null(installer)) {
    utils::install.packages(pkg)
  } else {
    installer(pkg)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package installation failed: ", pkg, call. = FALSE)
  }
  invisible(TRUE)
}

make_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_csv <- function(x, file) {
  make_dir(dirname(file))
  utils::write.csv(x, file, row.names = FALSE, quote = TRUE)
}

open_text <- function(file) {
  if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "rt") else base::file(file, "rt")
}

read_table_auto <- function(file, ...) {
  con <- open_text(file)
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, ...)
}

normalize_slash <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

find_project_dir <- function(start = getwd()) {
  d <- normalize_slash(start)
  repeat {
    if (dir.exists(file.path(d, "98panCancer"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) {
      stop("Cannot find project folder containing 98panCancer.", call. = FALSE)
    }
    d <- parent
  }
}

find_first_existing <- function(paths, required = TRUE) {
  paths <- normalize_slash(paths)
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("None of these files exists: ", paste(paths, collapse = "; "),
                     call. = FALSE)
  NA_character_
}

sample15 <- function(x) substr(as.character(x), 1, 15)

tcga_sample_type <- function(sample) {
  part4 <- vapply(strsplit(as.character(sample), "-", fixed = TRUE), function(z) {
    if (length(z) >= 4) z[4] else NA_character_
  }, character(1))
  code <- suppressWarnings(as.integer(substr(part4, 1, 2)))
  ifelse(is.na(code), NA_character_,
         ifelse(code >= 1 & code <= 9, "Tumor",
                ifelse(code >= 10 & code <= 19, "Normal", NA_character_)))
}

p_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

ordered_cancer_factor <- function(x) {
  u <- unique(as.character(x))
  factor(as.character(x), levels = unique(c(cancer_order, setdiff(u, cancer_order))))
}

save_ggplot <- function(p, file, width = 8, height = 5) {
  need_pkg("ggplot2")
  make_dir(dirname(file))
  ggplot2::ggsave(file, p, width = width, height = height, units = "in",
                  device = "pdf")
}

spearman_one <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < min_group_n || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(c(cor = NA_real_, pvalue = NA_real_, n = length(x)))
  }
  test <- suppressWarnings(stats::cor.test(x, y, method = cor_method,
                                           exact = FALSE))
  c(cor = unname(test$estimate), pvalue = test$p.value, n = length(x))
}

################################################################################
# Paths
################################################################################

project_dir <- find_project_dir(getwd())
course_dir <- file.path(project_dir, "98panCancer")
download_dir <- file.path(course_dir, "03.download")
out_dir <- file.path(project_dir, out_dir_name)

paths <- list(
  tcga_expression_dir = file.path(download_dir, "expression"),
  tcga_mutation_dir = file.path(download_dir, "mutation"),
  tcga_survival_dir = file.path(download_dir, "survival"),
  gtf_file = file.path(course_dir, "04.idTrans", "human.gtf"),
  xena_survival = find_first_existing(c(
    file.path(download_dir, "Survival_SupplementalTable_S1_20171025_xena_sp.gz"),
    file.path(download_dir, "Survival_SupplementalTable_S1_20171025_xena_sp"),
    file.path(course_dir, "10.DSS", "Survival_SupplementalTable_S1_20171025_xena_sp.gz"),
    file.path(course_dir, "10.DSS", "Survival_SupplementalTable_S1_20171025_xena_sp")
  ), required = FALSE),
  msi_file = find_first_existing(c(
    file.path(course_dir, "17.MSIcor", "MSI.txt"),
    file.path(project_dir, "MSI.txt")
  ), required = FALSE),
  cibersort_script = file.path(course_dir, "21.CIBERSORT", "panCancer21.CIBERSORT.R"),
  cibersort_ref = file.path(course_dir, "21.CIBERSORT", "ref.txt")
)

make_dir(out_dir)

################################################################################
# TCGA expression extraction
################################################################################

read_gtf_gene_map <- function(gtf_file, symbols) {
  if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file, call. = FALSE)
  }
  message("Reading GTF gene annotation...")
  lines <- readLines(gtf_file, warn = FALSE)
  gene_lines <- grep("\tgene\t", lines, value = TRUE)
  out <- setNames(vector("list", length(symbols)), symbols)
  for (sym in symbols) {
    hit <- grep(paste0('gene_name "', sym, '"'), gene_lines, value = TRUE)
    ids <- unique(sub('.*gene_id "([^"]+)".*', "\\1", hit))
    ids <- sub("\\..*$", "", ids)
    out[[sym]] <- ids[nzchar(ids)]
  }
  out
}

read_expression_file_for_genes <- function(file, gene_ids) {
  dat <- read_table_auto(file, check.names = FALSE)
  ens <- sub("\\..*$", "", dat[[1]])
  rows <- which(ens %in% gene_ids)
  if (length(rows) == 0) return(NULL)
  mat <- as.matrix(dat[rows, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  if (nrow(mat) > 1) {
    vals <- colMeans(mat, na.rm = TRUE)
  } else {
    vals <- as.numeric(mat[1, ])
  }
  names(vals) <- colnames(dat)[-1]
  vals
}

extract_tcga_gene_expression <- function(symbols, out_file) {
  gene_map <- read_gtf_gene_map(paths$gtf_file, symbols)
  all_gene_ids <- unique(unlist(gene_map))
  files <- list.files(paths$tcga_expression_dir,
                      pattern = "^TCGA-.*\\.htseq_fpkm\\.tsv(\\.gz)?$",
                      full.names = TRUE)
  if (length(files) == 0) {
    stop("No TCGA expression files found in ", paths$tcga_expression_dir,
         call. = FALSE)
  }

  out <- list()
  k <- 1
  for (file in files) {
    cancer <- sub("^TCGA-([A-Z0-9]+)\\..*$", "\\1", basename(file))
    message("Extracting expression from ", basename(file))
    dat <- read_table_auto(file, check.names = FALSE)
    ens <- sub("\\..*$", "", dat[[1]])
    keep_any <- ens %in% all_gene_ids
    if (!any(keep_any)) next
    mat_all <- as.matrix(dat[keep_any, -1, drop = FALSE])
    storage.mode(mat_all) <- "numeric"
    ens_keep <- ens[keep_any]
    for (sym in symbols) {
      rows <- which(ens_keep %in% gene_map[[sym]])
      if (length(rows) == 0) next
      vals <- if (length(rows) > 1) {
        colMeans(mat_all[rows, , drop = FALSE], na.rm = TRUE)
      } else {
        as.numeric(mat_all[rows, ])
      }
      samples <- colnames(mat_all)
      out[[k]] <- data.frame(
        Sample = samples,
        Sample15 = sample15(samples),
        Cancer = cancer,
        Group = tcga_sample_type(samples),
        Gene = sym,
        Expression = log2(vals + 1),
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  if (length(out) == 0) {
    stop("Target genes were not found in TCGA expression files.", call. = FALSE)
  }
  res <- do.call(rbind, out)
  write_csv(res, out_file)
  res
}

get_target_expression <- function(tcga_gene_long) {
  dat <- tcga_gene_long[tcga_gene_long$Gene == gene_symbol, ]
  stats::aggregate(Expression ~ Sample15 + Cancer + Group, data = dat, FUN = mean)
}

################################################################################
# TCGA+GTEx DLC1 extraction from UCSC Xena Toil
################################################################################

download_if_missing <- function(url, dest) {
  make_dir(dirname(dest))
  if (file.exists(dest) && grepl("\\.gz$", dest, ignore.case = TRUE)) {
    ok <- tryCatch({
      con <- gzfile(dest, "rt")
      on.exit(close(con), add = TRUE)
      readLines(con, n = 1)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) {
      message("Removing incomplete/corrupt download: ", dest)
      unlink(dest)
    }
  }
  if (!file.exists(dest)) {
    message("Downloading ", url)
    tmp <- paste0(dest, ".part")
    if (file.exists(tmp)) unlink(tmp)
    utils::download.file(url, destfile = tmp, mode = "wb", quiet = FALSE)
    if (grepl("\\.gz$", dest, ignore.case = TRUE)) {
      ok <- tryCatch({
        con <- gzfile(tmp, "rt")
        on.exit(close(con), add = TRUE)
        readLines(con, n = 1)
        TRUE
      }, error = function(e) FALSE)
      if (!ok) stop("Downloaded file is not a valid gzip archive: ", tmp,
                    call. = FALSE)
    }
    file.rename(tmp, dest)
  }
  dest
}

extract_gene_row_from_gz_matrix <- function(gz_file, gene_id_without_version) {
  con <- gzfile(gz_file, "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    first <- sub("\t.*$", "", line)
    if (sub("\\..*$", "", first) == gene_id_without_version) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      vals <- suppressWarnings(as.numeric(fields[-1]))
      names(vals) <- header[-1]
      return(vals)
    }
  }
  stop("Gene ID not found in Toil matrix: ", gene_id_without_version,
       call. = FALSE)
}

expand_gtex_to_cancer <- function(dat) {
  tcga <- dat[dat$Cohort == "TCGA", ]
  gtex <- dat[dat$Cohort == "GTEX", ]
  out <- list(tcga)
  k <- 2
  for (cancer in names(gtex_tissue_to_cancer)) {
    patterns <- gtex_tissue_to_cancer[[cancer]]
    hit <- rep(FALSE, nrow(gtex))
    for (pat in patterns) {
      hit <- hit | grepl(pat, gtex$DiseaseOrTissue,
                         ignore.case = TRUE, fixed = TRUE)
    }
    if (any(hit)) {
      tmp <- gtex[hit, ]
      tmp$Cancer <- cancer
      out[[k]] <- tmp
      k <- k + 1
    }
  }
  do.call(rbind, out)
}

build_tcga_gtex_dlc1_table <- function(out_subdir, tcga_sample_cancer_map) {
  make_dir(out_subdir)
  expr_gz <- download_if_missing(
    toil_expression_url,
    file.path(out_subdir, "TcgaTargetGtex_rsem_gene_tpm.gz")
  )
  pheno_gz <- download_if_missing(
    toil_phenotype_url,
    file.path(out_subdir, "TcgaTargetGTEX_phenotype.txt.gz")
  )

  dlc1 <- extract_gene_row_from_gz_matrix(expr_gz, ensembl_gene_id)
  phenotype <- read_table_auto(pheno_gz, check.names = FALSE)
  sample_col <- "sample"
  study_col <- "_study"
  sample_type_col <- "_sample_type"
  detailed_col <- "detailed_category"
  disease_col <- "primary disease or tissue"

  dat <- data.frame(
    Sample = names(dlc1),
    Sample15 = sample15(names(dlc1)),
    Expression = as.numeric(dlc1),
    stringsAsFactors = FALSE
  )
  dat <- merge(dat, phenotype, by.x = "Sample", by.y = sample_col)
  dat$Cohort <- dat[[study_col]]
  dat$Group <- ifelse(dat$Cohort == "TCGA" &
                        grepl("Tumor", dat[[sample_type_col]], ignore.case = TRUE),
                      "Tumor",
                      ifelse(dat$Cohort == "TCGA" &
                               grepl("Normal", dat[[sample_type_col]], ignore.case = TRUE),
                             "Normal",
                             ifelse(dat$Cohort == "GTEX", "Normal", NA_character_)))
  dat <- merge(dat, tcga_sample_cancer_map, by = "Sample15", all.x = TRUE)
  dat$Cancer <- ifelse(dat$Cohort == "TCGA", dat$Cancer, NA_character_)
  dat$DiseaseOrTissue <- if (disease_col %in% names(dat)) {
    dat[[disease_col]]
  } else {
    dat[[detailed_col]]
  }

  dat <- dat[dat$Cohort %in% c("TCGA", "GTEX") & !is.na(dat$Group), ]
  dat <- expand_gtex_to_cancer(dat)
  dat <- dat[!is.na(dat$Cancer), ]
  write_csv(dat, file.path(out_subdir, "DLC1_TCGA_GTEx_Toil_expression.csv"))
  dat
}

################################################################################
# Differential expression
################################################################################

plot_difference <- function(dat, out_subdir, label) {
  need_pkg("ggplot2")
  make_dir(out_subdir)
  dat <- dat[dat$Group %in% c("Tumor", "Normal") & !is.na(dat$Cancer), ]
  dat$Cancer <- ordered_cancer_factor(dat$Cancer)
  dat$Group <- factor(dat$Group, levels = c("Normal", "Tumor"))

  stats_list <- lapply(split(dat, dat$Cancer), function(x) {
    if (nrow(x) == 0) return(NULL)
    n_t <- sum(x$Group == "Tumor")
    n_n <- sum(x$Group == "Normal")
    if (n_t < min_group_n || n_n < min_group_n) return(NULL)
    p <- suppressWarnings(stats::wilcox.test(Expression ~ Group, data = x)$p.value)
    data.frame(
      Cancer = as.character(x$Cancer[1]),
      n_Tumor = n_t,
      n_Normal = n_n,
      median_Tumor = stats::median(x$Expression[x$Group == "Tumor"], na.rm = TRUE),
      median_Normal = stats::median(x$Expression[x$Group == "Normal"], na.rm = TRUE),
      pvalue = p,
      pstar = p_stars(p),
      stringsAsFactors = FALSE
    )
  })
  stat <- do.call(rbind, stats_list)
  if (is.null(stat) || nrow(stat) == 0) {
    stat <- data.frame(
      Cancer = character(), n_Tumor = integer(), n_Normal = integer(),
      median_Tumor = numeric(), median_Normal = numeric(), pvalue = numeric(),
      pstar = character(), FDR = numeric(), stringsAsFactors = FALSE
    )
    write_csv(stat, file.path(out_subdir, paste0(label, "_DLC1_wilcox.csv")))
    message("No tumor-normal comparison passed the sample-size filter for ", label)
    return(invisible(stat))
  }
  stat$FDR <- stats::p.adjust(stat$pvalue, method = "BH")
  stat <- stat[order(match(stat$Cancer, cancer_order)), ]
  write_csv(stat, file.path(out_subdir, paste0(label, "_DLC1_wilcox.csv")))

  ymax <- stats::aggregate(Expression ~ Cancer, dat, max, na.rm = TRUE)
  names(ymax)[2] <- "y"
  anno <- merge(stat[stat$pvalue < 0.05, ], ymax, by = "Cancer", all.x = TRUE)
  anno$Cancer <- ordered_cancer_factor(anno$Cancer)
  anno$y <- anno$y + 0.06 * diff(range(dat$Expression, na.rm = TRUE))

  p <- ggplot2::ggplot(dat, ggplot2::aes(Cancer, Expression, fill = Group)) +
    ggplot2::geom_violin(scale = "width", trim = FALSE, alpha = 0.55,
                         position = ggplot2::position_dodge(0.82), size = 0.2) +
    ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.25,
                          position = ggplot2::position_dodge(0.82), size = 0.2) +
    ggplot2::geom_text(data = anno,
                       ggplot2::aes(Cancer, y, label = pstar),
                       inherit.aes = FALSE, size = 3) +
    ggplot2::scale_fill_manual(values = c(Normal = "#2B6CB0",
                                          Tumor = "#B83227")) +
    ggplot2::labs(title = paste0(gene_symbol, " expression: ", label),
                  x = NULL, y = paste0(gene_symbol, " expression")) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                   legend.position = "top",
                   panel.grid.minor = ggplot2::element_blank())
  save_ggplot(p, file.path(out_subdir, paste0(label, "_DLC1_violin.pdf")),
              width = 12, height = 5.5)
  invisible(stat)
}

################################################################################
# Survival
################################################################################

read_xena_survival <- function(file) {
  if (is.na(file) || !file.exists(file)) {
    stop("Xena survival table is missing.", call. = FALSE)
  }
  dat <- read_table_auto(file, check.names = FALSE)
  dat$Sample15 <- sample15(dat$sample)
  dat
}

run_survival_endpoint <- function(expr_tumor, clinical, endpoint, out_subdir) {
  need_pkg("survival")
  make_dir(out_subdir)
  time_col <- paste0(endpoint, ".time")
  if (!(endpoint %in% names(clinical)) || !(time_col %in% names(clinical))) {
    message("Skipping ", endpoint, ": missing survival columns.")
    return(NULL)
  }
  cli <- clinical[, c("Sample15", endpoint, time_col)]
  names(cli) <- c("Sample15", "Status", "TimeDays")
  cli$Status <- suppressWarnings(as.numeric(cli$Status))
  cli$TimeDays <- suppressWarnings(as.numeric(cli$TimeDays))
  dat <- merge(expr_tumor, cli, by = "Sample15")
  dat <- dat[is.finite(dat$Expression) & is.finite(dat$Status) &
               is.finite(dat$TimeDays) & dat$TimeDays > 0, ]
  dat$TimeYears <- dat$TimeDays / 365

  cox_out <- list()
  km_out <- list()
  k <- 1
  m <- 1
  for (ct in cancer_order) {
    x <- dat[dat$Cancer == ct, ]
    if (nrow(x) < min_survival_n ||
        sum(x$Status == 1, na.rm = TRUE) < min_survival_events ||
        stats::sd(x$Expression, na.rm = TRUE) == 0) next

    fit <- survival::coxph(survival::Surv(TimeYears, Status) ~ Expression,
                           data = x)
    s <- summary(fit)
    cox_out[[k]] <- data.frame(
      Cancer = ct,
      Pvalue = s$coefficients[1, "Pr(>|z|)"],
      beta = s$coefficients[1, "coef"],
      HR = s$conf.int[1, "exp(coef)"],
      lower = s$conf.int[1, "lower .95"],
      upper = s$conf.int[1, "upper .95"],
      n = nrow(x),
      events = sum(x$Status == 1),
      stringsAsFactors = FALSE
    )
    k <- k + 1

    med <- stats::median(x$Expression, na.rm = TRUE)
    x$Group2 <- factor(ifelse(x$Expression > med, "High", "Low"),
                       levels = c("Low", "High"))
    sdif <- survival::survdiff(survival::Surv(TimeYears, Status) ~ Group2,
                               data = x)
    p_logrank <- 1 - stats::pchisq(sdif$chisq, df = 1)
    fit2 <- survival::coxph(survival::Surv(TimeYears, Status) ~ Group2,
                            data = x)
    s2 <- summary(fit2)
    km_out[[m]] <- data.frame(
      Cancer = ct,
      Pvalue = p_logrank,
      beta = s2$coefficients[1, "coef"],
      HR = s2$conf.int[1, "exp(coef)"],
      lower = s2$conf.int[1, "lower .95"],
      upper = s2$conf.int[1, "upper .95"],
      cutoff = med,
      n = nrow(x),
      events = sum(x$Status == 1),
      stringsAsFactors = FALSE
    )
    m <- m + 1

    if (p_logrank < km_p_cutoff_for_curve && maybe_pkg("survminer")) {
      sf <- survival::survfit(survival::Surv(TimeYears, Status) ~ Group2,
                              data = x)
      gp <- survminer::ggsurvplot(
        sf, data = x, risk.table = TRUE, conf.int = FALSE,
        pval = paste0("p=", format_p(p_logrank)),
        legend.title = paste0(gene_symbol, " expression"),
        legend.labs = c("Low", "High"),
        xlab = "Time (years)",
        ylab = paste0(endpoint, " probability"),
        palette = c("#2B6CB0", "#B83227"),
        title = paste0(endpoint, " in ", ct)
      )
      make_dir(file.path(out_subdir, "KM_curves"))
      grDevices::pdf(file.path(out_subdir, "KM_curves",
                               paste0(endpoint, "_", ct, "_KM.pdf")),
                     width = 6.2, height = 5.6, onefile = FALSE)
      print(gp)
      grDevices::dev.off()
    }
  }

  cox <- if (length(cox_out)) do.call(rbind, cox_out) else data.frame()
  km <- if (length(km_out)) do.call(rbind, km_out) else data.frame()
  if (nrow(cox) > 0) {
    cox$FDR <- stats::p.adjust(cox$Pvalue, method = "BH")
    write_csv(cox, file.path(out_subdir, paste0(endpoint, "_Cox.csv")))
    plot_forest(cox, file.path(out_subdir, paste0(endpoint, "_Cox_forest.pdf")),
                paste0(endpoint, " Cox"))
  }
  if (nrow(km) > 0) {
    km$FDR <- stats::p.adjust(km$Pvalue, method = "BH")
    write_csv(km, file.path(out_subdir, paste0(endpoint, "_KM.csv")))
    plot_forest(km, file.path(out_subdir, paste0(endpoint, "_KM_forest.pdf")),
                paste0(endpoint, " KM High vs Low"))
  }
  list(cox = cox, km = km)
}

plot_forest <- function(dat, file, title) {
  need_pkg("ggplot2")
  x <- dat[is.finite(dat$HR) & is.finite(dat$lower) & is.finite(dat$upper) &
             dat$lower > 0 & dat$upper > 0, ]
  if (nrow(x) == 0) return(invisible(NULL))
  x$Cancer <- factor(x$Cancer, levels = rev(unique(cancer_order[cancer_order %in% x$Cancer])))
  x$Significant <- x$Pvalue < 0.05
  p <- ggplot2::ggplot(x, ggplot2::aes(HR, Cancer)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "grey45") +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = Cancer),
                          size = 0.45, color = "grey35") +
    ggplot2::geom_point(ggplot2::aes(color = Significant), size = 2.2) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#B83227",
                                           `FALSE` = "#4A5568")) +
    ggplot2::labs(title = title, x = "Hazard ratio (log scale)", y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "none",
                   panel.grid.minor = ggplot2::element_blank())
  save_ggplot(p, file, width = 7.2, height = max(4.2, 0.22 * nrow(x) + 1.5))
}

################################################################################
# TMB and MSI
################################################################################

calculate_tmb_from_mutation_files <- function(out_subdir) {
  make_dir(out_subdir)
  files <- list.files(paths$tcga_mutation_dir,
                      pattern = "^TCGA-.*\\.varscan2_snv\\.tsv(\\.gz)?$",
                      full.names = TRUE)
  out <- list()
  k <- 1
  for (file in files) {
    cancer <- sub("^TCGA-([A-Z0-9]+)\\..*$", "\\1", basename(file))
    message("Calculating TMB from ", basename(file))
    mut <- read_table_auto(file, check.names = FALSE)
    if (!("Sample_ID" %in% names(mut))) next
    if ("filter" %in% names(mut)) {
      mut <- mut[mut$filter == "PASS", ]
    }
    if (identical(tmb_count_mode, "pass_nonsynonymous_snv") && "effect" %in% names(mut)) {
      synonymous_or_non_coding <- paste(
        c("synonymous_variant", "intron_variant", "3_prime_UTR_variant",
          "5_prime_UTR_variant", "upstream_gene_variant",
          "downstream_gene_variant", "intergenic_region",
          "non_coding_transcript", "NMD_transcript_variant"),
        collapse = "|"
      )
      mut <- mut[!grepl(synonymous_or_non_coding, mut$effect), ]
    }
    tab <- as.data.frame(table(mut$Sample_ID), stringsAsFactors = FALSE)
    names(tab) <- c("Sample", "Mutation_count")
    tab$Sample15 <- sample15(tab$Sample)
    tab$Cancer <- cancer
    tab$TMB <- tab$Mutation_count / exome_size_mb
    out[[k]] <- tab
    k <- k + 1
  }
  if (length(out) == 0) {
    stop("No mutation records were available for TMB calculation.", call. = FALSE)
  }
  res <- do.call(rbind, out)
  write_csv(res, file.path(out_subdir, "TCGA_TMB_from_varscan2_snv.csv"))
  res
}

correlate_score_with_dlc1 <- function(score_dat, expr_tumor, score_col,
                                      label, out_subdir) {
  make_dir(out_subdir)
  dat <- merge(score_dat, expr_tumor[, c("Sample15", "Cancer", "Expression")],
               by = c("Sample15", "Cancer"))
  write_csv(dat, file.path(out_subdir, paste0(label, "_DLC1_merged.csv")))
  if (nrow(dat) == 0) {
    res <- data.frame(Cancer = character(), cor = numeric(), pvalue = numeric(),
                      n = integer(), pstar = character(), FDR = numeric(),
                      stringsAsFactors = FALSE)
    write_csv(res, file.path(out_subdir, paste0(label, "_DLC1_correlation.csv")))
    message("No matched samples for ", label, " correlation.")
    return(invisible(res))
  }
  res <- do.call(rbind, lapply(split(dat, dat$Cancer), function(x) {
    r <- spearman_one(x[[score_col]], x$Expression)
    data.frame(Cancer = as.character(x$Cancer[1]), cor = r["cor"],
               pvalue = r["pvalue"], n = r["n"], pstar = p_stars(r["pvalue"]),
               stringsAsFactors = FALSE)
  }))
  res$FDR <- stats::p.adjust(res$pvalue, method = "BH")
  res <- res[order(match(res$Cancer, cancer_order)), ]
  write_csv(res, file.path(out_subdir, paste0(label, "_DLC1_correlation.csv")))
  plot_cor_bar(res, file.path(out_subdir, paste0(label, "_DLC1_cor_bar.pdf")),
               paste0(gene_symbol, " vs ", label))
  invisible(res)
}

read_msi_table <- function(file) {
  dat <- read_table_auto(file, check.names = FALSE)
  names(dat) <- sub("^id$", "Sample", names(dat))
  dat$Sample15 <- sample15(dat$Sample)
  names(dat) <- sub("^CancerType$", "Cancer", names(dat))
  dat
}

plot_cor_bar <- function(dat, file, title) {
  need_pkg("ggplot2")
  x <- dat[is.finite(dat$cor), ]
  if (nrow(x) == 0) return(invisible(NULL))
  x$Cancer <- ordered_cancer_factor(x$Cancer)
  lim <- range(c(0, x$cor), na.rm = TRUE) + c(-0.08, 0.08)
  p <- ggplot2::ggplot(x, ggplot2::aes(Cancer, cor, fill = cor > 0)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey45", size = 0.35) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(ggplot2::aes(label = pstar,
                                    y = ifelse(cor >= 0, cor + 0.025,
                                               cor - 0.025)),
                       size = 3) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#B83227",
                                          `FALSE` = "#2B6CB0")) +
    ggplot2::coord_cartesian(ylim = lim) +
    ggplot2::labs(title = title, x = NULL, y = "Spearman rho") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                   legend.position = "none",
                   panel.grid.minor = ggplot2::element_blank())
  save_ggplot(p, file, width = 9.5, height = 4.5)
}

################################################################################
# Immune deconvolution and checkpoint correlation
################################################################################

build_symbol_expression_matrix <- function(symbols, out_file) {
  gene_long <- extract_tcga_gene_expression(symbols, out_file = tempfile(fileext = ".csv"))
  tumor <- gene_long[gene_long$Group == "Tumor", ]
  wide <- stats::reshape(tumor[, c("Sample15", "Cancer", "Gene", "Expression")],
                         idvar = c("Sample15", "Cancer"),
                         timevar = "Gene",
                         direction = "wide")
  names(wide) <- sub("^Expression\\.", "", names(wide))
  write_csv(wide, out_file)
  wide
}

read_gtf_ensembl_to_symbol <- function(gtf_file) {
  lines <- readLines(gtf_file, warn = FALSE)
  gene_lines <- grep("\tgene\t", lines, value = TRUE)
  ens <- sub('.*gene_id "([^"]+)".*', "\\1", gene_lines)
  ens <- sub("\\..*$", "", ens)
  sym <- sub('.*gene_name "([^"]+)".*', "\\1", gene_lines)
  data.frame(Ensembl = ens, Symbol = sym, stringsAsFactors = FALSE)
}

build_tcga_symbol_matrix_for_immunedeconv <- function(out_file) {
  need_pkg("limma")
  gene_map <- read_gtf_ensembl_to_symbol(paths$gtf_file)
  files <- list.files(paths$tcga_expression_dir,
                      pattern = "^TCGA-.*\\.htseq_fpkm\\.tsv(\\.gz)?$",
                      full.names = TRUE)
  all_mats <- list()
  for (file in files) {
    message("Building symbol matrix from ", basename(file))
    dat <- read_table_auto(file, check.names = FALSE)
    ens <- sub("\\..*$", "", dat[[1]])
    idx <- match(ens, gene_map$Ensembl)
    keep <- !is.na(idx)
    if (!any(keep)) next
    mat <- as.matrix(dat[keep, -1, drop = FALSE])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- gene_map$Symbol[idx[keep]]
    sample_type <- tcga_sample_type(colnames(mat))
    mat <- mat[, sample_type == "Tumor", drop = FALSE]
    if (ncol(mat) == 0) next
    if (identical(immune_expression_scale, "log2_fpkm_plus_1")) {
      mat <- log2(mat + 1)
    } else if (identical(immune_expression_scale, "tpm_from_fpkm")) {
      mat <- sweep(mat, 2, colSums(mat, na.rm = TRUE) / 1e6, "/")
      mat[!is.finite(mat)] <- 0
    }
    mat <- limma::avereps(mat)
    all_mats[[basename(file)]] <- mat
  }
  if (length(all_mats) == 0) {
    stop("No tumor expression matrices were built for immune deconvolution.",
         call. = FALSE)
  }
  common_genes <- Reduce(intersect, lapply(all_mats, rownames))
  if (length(common_genes) == 0) {
    stop("No common genes across TCGA tumor matrices for immune deconvolution.",
         call. = FALSE)
  }
  merged <- do.call(cbind, lapply(all_mats, function(x) x[common_genes, , drop = FALSE]))
  merged <- merged[rowMeans(merged, na.rm = TRUE) > 0, , drop = FALSE]
  make_dir(dirname(out_file))
  utils::write.table(
    cbind(GeneSymbol = rownames(merged), as.data.frame(merged, check.names = FALSE)),
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  merged
}

read_immunedeconv_expression_matrix <- function(expr_matrix_file) {
  rds_file <- paste0(expr_matrix_file, ".rds")
  if (file.exists(rds_file)) {
    message("Reading cached immune expression RDS: ", rds_file)
    return(readRDS(rds_file))
  }
  message("Reading immune expression matrix: ", expr_matrix_file)
  if (maybe_pkg("data.table")) {
    dat <- data.table::fread(expr_matrix_file, data.table = FALSE,
                             check.names = FALSE)
  } else {
    dat <- utils::read.delim(expr_matrix_file, check.names = FALSE)
  }
  genes <- dat[[1]]
  dat[[1]] <- NULL
  mat <- as.matrix(dat)
  rownames(mat) <- genes
  storage.mode(mat) <- "numeric"
  saveRDS(mat, rds_file, compress = FALSE)
  mat
}

run_checkpoint_correlation <- function(gene_wide, out_subdir) {
  make_dir(out_subdir)
  out <- list()
  k <- 1
  for (ct in cancer_order) {
    x <- gene_wide[gene_wide$Cancer == ct, ]
    if (nrow(x) < min_group_n || !(gene_symbol %in% names(x))) next
    for (ck in checkpoint_genes) {
      if (!(ck %in% names(x))) next
      r <- spearman_one(x[[ck]], x[[gene_symbol]])
      out[[k]] <- data.frame(Cancer = ct, Gene = ck, Cor = r["cor"],
                             pvalue = r["pvalue"], n = r["n"],
                             pstar = p_stars(r["pvalue"]),
                             stringsAsFactors = FALSE)
      k <- k + 1
    }
  }
  if (length(out) == 0) {
    res <- data.frame(Cancer = character(), Gene = character(), Cor = numeric(),
                      pvalue = numeric(), n = integer(), pstar = character(),
                      FDR = numeric(), stringsAsFactors = FALSE)
    write_csv(res, file.path(out_subdir, "DLC1_checkpoint_correlation.csv"))
    message("No immune-checkpoint correlations passed filters.")
    return(res)
  }
  res <- do.call(rbind, out)
  res$FDR <- stats::p.adjust(res$pvalue, method = "BH")
  write_csv(res, file.path(out_subdir, "DLC1_checkpoint_correlation.csv"))
  plot_heatmap(res, file.path(out_subdir, "DLC1_checkpoint_heatmap.pdf"),
               row_col = "Gene", title = "DLC1 and immune checkpoints")
  res
}

plot_heatmap <- function(dat, file, row_col, title) {
  need_pkg("ggplot2")
  x <- dat[is.finite(dat$Cor), ]
  if (nrow(x) == 0) return(invisible(NULL))
  x$Cancer <- ordered_cancer_factor(x$Cancer)
  x[[row_col]] <- factor(x[[row_col]], levels = rev(unique(x[[row_col]])))
  max_abs <- max(abs(x$Cor), na.rm = TRUE)
  max_abs <- max(max_abs, 0.1)
  p <- ggplot2::ggplot(x, ggplot2::aes(Cancer, .data[[row_col]], fill = Cor)) +
    ggplot2::geom_tile(color = "white", size = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = pstar), size = 2.8) +
    ggplot2::scale_fill_gradient2(low = "#2B6CB0", mid = "white",
                                  high = "#B83227", midpoint = 0,
                                  limits = c(-max_abs, max_abs),
                                  name = "Spearman rho") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                   panel.grid = ggplot2::element_blank())
  save_ggplot(p, file, width = 10, height = max(4.2, 0.22 * length(unique(x[[row_col]])) + 2))
}

install_immunedeconv_dependencies <- function() {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  options(timeout = max(3600, getOption("timeout")))
  message("R version: ", R.version.string)
  install_if_missing("BiocManager")
  install_if_missing("remotes")

  bioc_pkgs <- c(
    "BiocParallel", "Biobase", "biomaRt", "preprocessCore", "sva", "GSVA",
    "GSEABase", "limma"
  )
  for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing Bioconductor package: ", pkg)
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }

  cran_pkgs <- c(
    "purrr", "dplyr", "magrittr", "readr", "readxl", "testit", "tibble",
    "data.tree", "limSolve", "e1071", "rlang", "matrixStats", "stringr"
  )
  for (pkg in cran_pkgs) install_if_missing(pkg)

  github_pkgs <- c(
    "dviraran/xCell",
    "GfellerLab/EPIC",
    "grst/MCPcounter",
    "cansysbio/ConsensusTME"
  )
  for (repo in github_pkgs) {
    pkg <- sub(".*/", "", repo)
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing GitHub package: ", repo)
      remotes::install_github(repo, dependencies = TRUE, upgrade = "never")
    } else {
      message("Already installed: ", pkg)
    }
  }

  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    message("Installing GitHub package: omnideconv/immunedeconv@v2.0.3")
    remotes::install_github(
      "omnideconv/immunedeconv",
      ref = "v2.0.3",
      dependencies = TRUE,
      upgrade = "never"
    )
  }

  message("Final immunedeconv availability: ",
          requireNamespace("immunedeconv", quietly = TRUE))
  if (requireNamespace("immunedeconv", quietly = TRUE)) {
    message("immunedeconv version: ",
            as.character(utils::packageVersion("immunedeconv")))
  }
  invisible(TRUE)
}

setup_method_dependencies <- function(method, out_subdir) {
  if (identical(method, "xcell")) {
    suppressPackageStartupMessages(library(xCell))
  }
  if (identical(method, "epic")) {
    suppressPackageStartupMessages(library(EPIC))
  }
  if (identical(method, "cibersort")) {
    if (file.exists(paths$cibersort_script)) {
      cibersort_wrapper <- file.path(out_subdir,
                                     "CIBERSORT_wrapper_for_immunedeconv.R")
      wrapper_lines <- c(
        "suppressPackageStartupMessages(library(e1071))",
        "suppressPackageStartupMessages(library(parallel))",
        "suppressPackageStartupMessages(library(preprocessCore))",
        paste0("source(", deparse(paths$cibersort_script), ")"),
        "CIBERSORT_original <- CIBERSORT",
        "CIBERSORT <- function(sig_matrix, mixture_file, perm = 0, QN = TRUE, ...) {",
        "  CIBERSORT_original(sig_matrix, mixture_file, perm = perm, QN = QN)",
        "}"
      )
      writeLines(wrapper_lines, cibersort_wrapper)
      immunedeconv::set_cibersort_binary(cibersort_wrapper)
    }
    if (file.exists(paths$cibersort_ref)) {
      immunedeconv::set_cibersort_mat(paths$cibersort_ref)
    }
  }
  invisible(TRUE)
}

score_to_correlation <- function(score, expr_tumor, method, cancer_filter = NULL) {
  score <- as.data.frame(t(score), stringsAsFactors = FALSE)
  colnames(score) <- as.character(score[1, ])
  score <- score[-1, , drop = FALSE]
  names(score) <- make.names(names(score), unique = TRUE)
  score$Sample15 <- sample15(rownames(score))
  score <- merge(score, expr_tumor[, c("Sample15", "Cancer", "Expression")],
                 by = "Sample15")
  if (!is.null(cancer_filter)) {
    score <- score[score$Cancer == cancer_filter, ]
  }
  score_cols <- setdiff(names(score), c("Sample15", "Cancer", "Expression"))
  out <- list()
  k <- 1
  cancers <- if (is.null(cancer_filter)) cancer_order else cancer_filter
  for (cell in score_cols) {
    score[[cell]] <- suppressWarnings(as.numeric(score[[cell]]))
    for (ct in cancers) {
      x <- score[score$Cancer == ct, ]
      if (nrow(x) < min_group_n) next
      r <- spearman_one(x[[cell]], x$Expression)
      out[[k]] <- data.frame(Cancer = ct, Method = method, Cell = cell,
                             Cor = r["cor"], pvalue = r["pvalue"], n = r["n"],
                             pstar = p_stars(r["pvalue"]),
                             stringsAsFactors = FALSE)
      k <- k + 1
    }
  }
  if (length(out) == 0) {
    return(data.frame(
      Cancer = character(), Method = character(), Cell = character(),
      Cor = numeric(), pvalue = numeric(), n = integer(), pstar = character(),
      FDR = numeric(), stringsAsFactors = FALSE
    ))
  }
  cor_res <- do.call(rbind, out)
  cor_res$FDR <- stats::p.adjust(cor_res$pvalue, method = "BH")
  cor_res
}

run_immunedeconv_optional <- function(expr_matrix_file, expr_tumor, out_subdir) {
  make_dir(out_subdir)
  if (!maybe_pkg("immunedeconv")) {
    msg <- paste(
      "immunedeconv is not installed. To run six immune algorithms, install it",
      "following the package documentation, then set run_immune_deconvolution <- TRUE.",
      "CIBERSORT additionally requires licensed CIBERSORT resources."
    )
    writeLines(msg, file.path(out_subdir, "IMMUNE_DECONVOLUTION_NOT_RUN.txt"))
    message(msg)
    return(invisible(NULL))
  }
  if (!run_immune_deconvolution) {
    writeLines("Set run_immune_deconvolution <- TRUE to run immunedeconv.",
               file.path(out_subdir, "IMMUNE_DECONVOLUTION_DISABLED.txt"))
    return(invisible(NULL))
  }
  unlink(file.path(out_subdir, c("IMMUNE_DECONVOLUTION_NOT_RUN.txt",
                                 "IMMUNE_DECONVOLUTION_DISABLED.txt")))
  if (!file.exists(expr_matrix_file)) {
    expr <- build_tcga_symbol_matrix_for_immunedeconv(expr_matrix_file)
  } else {
    expr <- read_immunedeconv_expression_matrix(expr_matrix_file)
  }
  storage.mode(expr) <- "numeric"
  expr_sample_map <- unique(expr_tumor[, c("Sample15", "Cancer")])
  indications <- expr_sample_map$Cancer[match(sample15(colnames(expr)),
                                              expr_sample_map$Sample15)]
  indications <- tolower(indications)
  res_all <- list()
  for (method in immune_methods) {
    message("Running immunedeconv method: ", method)
    expr_method <- expr
    indications_method <- indications
    if (identical(method, "timer")) {
      supported <- indications_method %in% immunedeconv::timer_available_cancers
      if (any(!supported, na.rm = TRUE)) {
        message("Dropping ", sum(!supported, na.rm = TRUE),
                " samples with TIMER-unsupported cancer indications.")
      }
      expr_method <- expr_method[, supported, drop = FALSE]
      indications_method <- indications_method[supported]
    }
    score <- tryCatch({
      setup_method_dependencies(method, out_subdir)
      if (identical(method, "timer")) {
        immunedeconv::deconvolute(expr_method, method,
                                  indications = indications_method)
      } else {
        immunedeconv::deconvolute(expr_method, method)
      }
    }, error = function(e) {
      msg <- paste0("Method ", method, " failed: ", conditionMessage(e))
      writeLines(msg, file.path(out_subdir, paste0(method, "_NOT_RUN.txt")))
      message(msg)
      NULL
    })
    if (is.null(score)) next
    write_csv(as.data.frame(score, stringsAsFactors = FALSE),
              file.path(out_subdir, paste0(method, "_raw_score.csv")))
    cor_res <- score_to_correlation(score, expr_tumor, method)
    write_csv(cor_res, file.path(out_subdir,
                                 paste0(method, "_DLC1_immune_correlation.csv")))
    plot_heatmap(cor_res, file.path(out_subdir,
                                    paste0(method, "_DLC1_immune_heatmap.pdf")),
                 row_col = "Cell", title = paste0("DLC1 immune correlation: ", method))
    res_all[[method]] <- cor_res
  }
  invisible(res_all)
}

load_cached_expr_tumor <- function() {
  gene_long_file <- file.path(out_dir, "00_TCGA_target_gene_expression_long.csv")
  if (!file.exists(gene_long_file)) {
    core_symbols <- unique(c(gene_symbol, checkpoint_genes))
    extract_tcga_gene_expression(core_symbols, gene_long_file)
  }
  tcga_gene_long <- read_table_auto(gene_long_file, sep = ",",
                                    check.names = FALSE)
  dlc1_expr <- get_target_expression(tcga_gene_long)
  dlc1_expr[dlc1_expr$Group == "Tumor", ]
}

run_one_immune_method <- function(method) {
  if (is.na(method) || !nzchar(method)) {
    stop("Please provide --method, e.g. --method xcell.", call. = FALSE)
  }
  method <- tolower(method)
  if (!(method %in% immune_methods)) {
    stop("Unknown immune method: ", method, call. = FALSE)
  }
  if (!maybe_pkg("immunedeconv")) {
    stop("immunedeconv is not installed. Run --mode install_immunedeconv first.",
         call. = FALSE)
  }

  expr_tumor <- load_cached_expr_tumor()
  out_subdir <- file.path(out_dir, "06_immune_deconvolution")
  make_dir(out_subdir)
  unlink(file.path(out_subdir, c("IMMUNE_DECONVOLUTION_NOT_RUN.txt",
                                 "IMMUNE_DECONVOLUTION_DISABLED.txt")))
  expr_matrix_file <- file.path(out_dir,
                                "TCGA_symbol_expression_matrix_for_immunedeconv.txt")
  if (!file.exists(expr_matrix_file)) {
    expr <- build_tcga_symbol_matrix_for_immunedeconv(expr_matrix_file)
  } else {
    expr <- read_immunedeconv_expression_matrix(expr_matrix_file)
  }
  storage.mode(expr) <- "numeric"

  expr_sample_map <- unique(expr_tumor[, c("Sample15", "Cancer")])
  indications <- expr_sample_map$Cancer[match(sample15(colnames(expr)),
                                              expr_sample_map$Sample15)]
  indications <- tolower(indications)
  expr_method <- expr
  indications_method <- indications
  if (identical(method, "timer")) {
    supported <- indications_method %in% immunedeconv::timer_available_cancers
    message("TIMER supported samples: ", sum(supported, na.rm = TRUE),
            " / ", length(supported))
    expr_method <- expr_method[, supported, drop = FALSE]
    indications_method <- indications_method[supported]
  }

  setup_method_dependencies(method, out_subdir)
  message("Running immunedeconv method: ", method)
  score <- if (identical(method, "timer")) {
    immunedeconv::deconvolute(expr_method, method,
                              indications = indications_method)
  } else {
    immunedeconv::deconvolute(expr_method, method)
  }
  write_csv(as.data.frame(score, stringsAsFactors = FALSE),
            file.path(out_subdir, paste0(method, "_raw_score.csv")))
  cor_res <- score_to_correlation(score, expr_tumor, method)
  write_csv(cor_res, file.path(out_subdir,
                               paste0(method, "_DLC1_immune_correlation.csv")))
  plot_heatmap(cor_res, file.path(out_subdir,
                                  paste0(method, "_DLC1_immune_heatmap.pdf")),
               row_col = "Cell",
               title = paste0("DLC1 immune correlation: ", method))
  message("Finished method: ", method)
  invisible(cor_res)
}

prepare_immune_expression_by_cancer <- function() {
  expr_tumor <- load_cached_expr_tumor()
  expr_matrix_file <- file.path(out_dir,
                                "TCGA_symbol_expression_matrix_for_immunedeconv.txt")
  if (!file.exists(expr_matrix_file)) {
    expr <- build_tcga_symbol_matrix_for_immunedeconv(expr_matrix_file)
  } else {
    expr <- read_immunedeconv_expression_matrix(expr_matrix_file)
  }
  sample_info <- data.frame(
    Sample = colnames(expr),
    Sample15 = sample15(colnames(expr)),
    stringsAsFactors = FALSE
  )
  sample_info$Cancer <- expr_tumor$Cancer[match(sample_info$Sample15,
                                                expr_tumor$Sample15)]
  sample_info <- sample_info[!is.na(sample_info$Cancer), ]
  split_dir <- file.path(out_dir, "06_immune_deconvolution",
                         "expression_by_cancer")
  make_dir(split_dir)
  write_csv(sample_info, file.path(split_dir, "immune_expression_sample_info.csv"))
  for (ct in cancer_order) {
    cols <- sample_info$Sample[sample_info$Cancer == ct]
    cols <- intersect(cols, colnames(expr))
    if (length(cols) < min_group_n) next
    message("Saving ", ct, ": ", length(cols), " samples")
    saveRDS(expr[, cols, drop = FALSE],
            file.path(split_dir, paste0(ct, ".rds")),
            compress = FALSE)
  }
  message("Prepared by-cancer immune expression matrices.")
  invisible(sample_info)
}

run_one_immune_method_cancer <- function(method, ct) {
  if (is.na(method) || !nzchar(method)) {
    stop("Please provide --method.", call. = FALSE)
  }
  if (is.na(ct) || !nzchar(ct)) {
    stop("Please provide --cancer, e.g. --cancer LIHC.", call. = FALSE)
  }
  method <- tolower(method)
  ct <- toupper(ct)
  if (!(method %in% immune_methods)) {
    stop("Unknown immune method: ", method, call. = FALSE)
  }
  if (!maybe_pkg("immunedeconv")) {
    stop("immunedeconv is not installed. Run --mode install_immunedeconv first.",
         call. = FALSE)
  }

  expr_tumor <- load_cached_expr_tumor()
  split_dir <- file.path(out_dir, "06_immune_deconvolution",
                         "expression_by_cancer")
  expr_file <- file.path(split_dir, paste0(ct, ".rds"))
  if (!file.exists(expr_file)) {
    prepare_immune_expression_by_cancer()
  }
  if (!file.exists(expr_file)) {
    stop("Missing by-cancer expression RDS: ", expr_file, call. = FALSE)
  }
  expr <- readRDS(expr_file)
  out_subdir <- file.path(out_dir, "06_immune_deconvolution")
  method_dir <- file.path(out_subdir, "by_cancer", method)
  make_dir(method_dir)

  if (identical(method, "timer") &&
      !(tolower(ct) %in% immunedeconv::timer_available_cancers)) {
    writeLines(paste0("TIMER unsupported cancer: ", ct),
               file.path(method_dir, paste0(ct, "_NOT_RUN.txt")))
    return(invisible(NULL))
  }
  if (identical(method, "cibersort") && file.exists(paths$cibersort_ref)) {
    ref <- read_table_auto(paths$cibersort_ref, check.names = FALSE)
    ref_genes <- as.character(ref[[1]])
    keep <- intersect(rownames(expr), ref_genes)
    expr <- expr[keep, , drop = FALSE]
  }
  setup_method_dependencies(method, out_subdir)

  message("Running ", method, " for ", ct, " with ", ncol(expr), " samples")
  score <- if (identical(method, "timer")) {
    immunedeconv::deconvolute(expr, method,
                              indications = rep(tolower(ct), ncol(expr)))
  } else {
    immunedeconv::deconvolute(expr, method)
  }
  write_csv(as.data.frame(score, stringsAsFactors = FALSE),
            file.path(method_dir, paste0(ct, "_raw_score.csv")))
  cor_res <- score_to_correlation(score, expr_tumor, method,
                                  cancer_filter = ct)
  write_csv(cor_res, file.path(method_dir, paste0(ct, "_correlation.csv")))
  message("Finished ", method, " for ", ct)
  invisible(cor_res)
}

combine_split_immune_method <- function(method) {
  if (is.na(method) || !nzchar(method)) {
    stop("Please provide --method.", call. = FALSE)
  }
  method <- tolower(method)
  out_subdir <- file.path(out_dir, "06_immune_deconvolution")
  method_dir <- file.path(out_subdir, "by_cancer", method)
  files <- list.files(method_dir, pattern = "_correlation\\.csv$",
                      full.names = TRUE)
  if (length(files) == 0) {
    stop("No correlation files found for method: ", method, call. = FALSE)
  }
  dat <- do.call(rbind, lapply(files, function(f) {
    read_table_auto(f, sep = ",", check.names = FALSE)
  }))
  dat$FDR <- stats::p.adjust(dat$pvalue, method = "BH")
  dat <- dat[order(match(dat$Cancer, cancer_order), dat$Cell), ]
  write_csv(dat, file.path(out_subdir,
                           paste0(method, "_DLC1_immune_correlation.csv")))
  plot_heatmap(dat, file.path(out_subdir,
                              paste0(method, "_DLC1_immune_heatmap.pdf")),
               row_col = "Cell",
               title = paste0("DLC1 immune correlation: ", method))
  message("Combined split immune method: ", method)
  invisible(dat)
}

################################################################################
# Main
################################################################################

main <- function() {
  make_dir(out_dir)
  message("Project: ", project_dir)
  message("Output: ", out_dir)

  core_symbols <- unique(c(gene_symbol, checkpoint_genes))
  gene_long_file <- file.path(out_dir, "00_TCGA_target_gene_expression_long.csv")
  if (file.exists(gene_long_file)) {
    message("Reading cached target-gene expression: ", gene_long_file)
    tcga_gene_long <- read_table_auto(gene_long_file, sep = ",",
                                      check.names = FALSE)
  } else {
    tcga_gene_long <- extract_tcga_gene_expression(core_symbols, gene_long_file)
  }
  dlc1_expr <- get_target_expression(tcga_gene_long)
  write_csv(dlc1_expr, file.path(out_dir, "00_TCGA_DLC1_expression.csv"))

  message("Differential expression: TCGA only")
  diff_dir <- file.path(out_dir, "01_differential_expression")
  plot_difference(dlc1_expr, file.path(diff_dir, "TCGA_only"), "TCGA_only")

  if (run_tcga_gtex_download) {
    message("Differential expression: TCGA + GTEx from UCSC Xena Toil")
    toil_dir <- file.path(out_dir, "00_UCSC_Xena_Toil_TCGA_GTEx")
    tcga_sample_cancer_map <- unique(dlc1_expr[, c("Sample15", "Cancer")])
    gtex_dat <- build_tcga_gtex_dlc1_table(toil_dir, tcga_sample_cancer_map)
    tcga_gtex_plot <- gtex_dat[gtex_dat$Cohort == "TCGA" |
                                 gtex_dat$Cohort == "GTEX", ]
    plot_difference(tcga_gtex_plot, file.path(diff_dir, "TCGA_GTEx"),
                    "TCGA_GTEx_Toil")
  }

  expr_tumor <- dlc1_expr[dlc1_expr$Group == "Tumor", ]
  if (run_survival) {
    message("Survival analysis")
    clinical <- read_xena_survival(paths$xena_survival)
    for (ep in c("OS", "DSS", "DFI", "PFI")) {
      run_survival_endpoint(expr_tumor, clinical, ep,
                            file.path(out_dir, "02_survival", ep))
    }
  }

  if (run_tmb_from_mutation) {
    message("TMB analysis")
    tmb <- calculate_tmb_from_mutation_files(file.path(out_dir, "03_TMB"))
    correlate_score_with_dlc1(tmb, expr_tumor, "TMB", "TMB",
                              file.path(out_dir, "03_TMB"))
  }

  if (run_msi_if_available && !is.na(paths$msi_file) && file.exists(paths$msi_file)) {
    message("MSI analysis")
    msi <- read_msi_table(paths$msi_file)
    correlate_score_with_dlc1(msi, expr_tumor, "MSI", "MSI",
                              file.path(out_dir, "04_MSI"))
  }

  if (run_checkpoint_correlation_from_raw_expression) {
    message("Immune-checkpoint correlation from raw TCGA expression")
    checkpoint_long <- stats::aggregate(
      Expression ~ Sample15 + Cancer + Gene,
      data = tcga_gene_long[tcga_gene_long$Group == "Tumor",
                            c("Sample15", "Cancer", "Gene", "Expression")],
      FUN = mean
    )
    gene_wide <- stats::reshape(
      checkpoint_long,
      idvar = c("Sample15", "Cancer"),
      timevar = "Gene",
      direction = "wide"
    )
    names(gene_wide) <- sub("^Expression\\.", "", names(gene_wide))
    write_csv(gene_wide, file.path(out_dir, "05_checkpoint",
                                   "TCGA_checkpoint_gene_expression_wide.csv"))
    run_checkpoint_correlation(gene_wide, file.path(out_dir, "05_checkpoint"))
  }

  run_immunedeconv_optional(
    expr_matrix_file = file.path(out_dir, "TCGA_symbol_expression_matrix_for_immunedeconv.txt"),
    expr_tumor = expr_tumor,
    out_subdir = file.path(out_dir, "06_immune_deconvolution")
  )

  utils::capture.output(utils::sessionInfo(),
                        file = file.path(out_dir, "sessionInfo.txt"))
  message("Done.")
}

run_by_mode <- function() {
  make_dir(out_dir)
  message("Selected mode: ", analysis_mode)
  if (identical(analysis_mode, "full")) {
    main()
  } else if (identical(analysis_mode, "install_immunedeconv")) {
    install_immunedeconv_dependencies()
  } else if (identical(analysis_mode, "immune_all")) {
    run_immune_deconvolution <<- TRUE
    expr_tumor <- load_cached_expr_tumor()
    run_immunedeconv_optional(
      expr_matrix_file = file.path(out_dir,
                                   "TCGA_symbol_expression_matrix_for_immunedeconv.txt"),
      expr_tumor = expr_tumor,
      out_subdir = file.path(out_dir, "06_immune_deconvolution")
    )
    utils::capture.output(utils::sessionInfo(),
                          file = file.path(out_dir,
                                           "sessionInfo_immune_deconv.txt"))
  } else if (identical(analysis_mode, "immune_one")) {
    run_one_immune_method(selected_immune_method)
  } else if (identical(analysis_mode, "prepare_immune_by_cancer")) {
    prepare_immune_expression_by_cancer()
  } else if (identical(analysis_mode, "immune_cancer")) {
    run_one_immune_method_cancer(selected_immune_method, selected_cancer)
  } else if (identical(analysis_mode, "immune_split")) {
    prepare_immune_expression_by_cancer()
    for (ct in cancer_order) {
      tryCatch(
        run_one_immune_method_cancer(selected_immune_method, ct),
        error = function(e) {
          method_dir <- file.path(out_dir, "06_immune_deconvolution",
                                  "by_cancer", tolower(selected_immune_method))
          make_dir(method_dir)
          writeLines(conditionMessage(e),
                     file.path(method_dir, paste0(ct, "_NOT_RUN.txt")))
          message("Failed ", selected_immune_method, " ", ct, ": ",
                  conditionMessage(e))
        }
      )
    }
    combine_split_immune_method(selected_immune_method)
  } else if (identical(analysis_mode, "combine_split")) {
    combine_split_immune_method(selected_immune_method)
  } else {
    stop("Unknown --mode: ", analysis_mode, call. = FALSE)
  }
  invisible(TRUE)
}

run_by_mode()






############################################################
## DLC1最佳截断值
############################################################
############################################################
## DLC1 TCGA-LIHC survival analysis
## Optimal cutpoint version
## Endpoints: OS, DSS, DFI, PFI
############################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

############################################################
## 0. User paths
############################################################

gene_use <- "DLC1"

base_dir <- file.path(hcc_dlc1_root(), "bulk_survival_optimal_cutoff")

expr_dir <- file.path(base_dir, "expression")

supp_file <- file.path(
  base_dir,
  "Survival_SupplementalTable_S1_20171025_xena_sp.gz"
)

outdir <- file.path(
  base_dir,
  "DLC1_LIHC_4survival_optimal_cutoff_final"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "rds"), recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Packages
############################################################

cran_pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ggpubr",
  "survival", "survminer", "stringr", "tibble", "scales",
  "patchwork", "readr", "forestplot"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(stringr)
  library(tibble)
  library(scales)
  library(patchwork)
  library(readr)
  library(forestplot)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise

############################################################
## 2. Helper functions
############################################################

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.45, color = "black"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black")
    )
}

safe_ggsave <- function(filename, plot, width = 6, height = 5) {
  outfile <- file.path(outdir, "plots", filename)
  tryCatch({
    ggsave(outfile, plot, width = width, height = height, device = cairo_pdf)
  }, error = function(e) {
    message("cairo_pdf failed, using default pdf for: ", filename)
    ggsave(outfile, plot, width = width, height = height)
  })
}

clean_tcga_patient <- function(x) {
  substr(as.character(x), 1, 12)
}

get_sample_type <- function(x) {
  code <- substr(as.character(x), 14, 15)
  ifelse(
    code %in% c("01"),
    "Tumor",
    ifelse(code %in% c("11"), "Normal", "Other")
  )
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

status_to_numeric <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    return(as.numeric(x))
  }

  x2 <- as.character(x)
  x2 <- trimws(x2)

  out <- rep(NA_real_, length(x2))

  out[x2 %in% c("1", "Dead", "DEAD", "dead", "deceased", "Deceased", "TRUE", "True", "true")] <- 1
  out[x2 %in% c("0", "Alive", "ALIVE", "alive", "living", "Living", "FALSE", "False", "false")] <- 0

  suppressWarnings({
    num_x <- as.numeric(x2)
  })

  out[is.na(out) & !is.na(num_x)] <- num_x[is.na(out) & !is.na(num_x)]

  return(out)
}

make_valid_filename <- function(x) {
  x <- gsub("[/\\:*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

############################################################
## 3. Locate expression file
############################################################

expr_files <- list.files(
  expr_dir,
  pattern = "LIHC|lihc|TCGA-LIHC|tcga-lihc|expression|expr|fpkm|tpm|htseq|star",
  full.names = TRUE,
  recursive = TRUE
)

expr_files <- expr_files[
  grepl("\\.txt$|\\.tsv$|\\.csv$|\\.gz$|\\.rda$|\\.RData$|\\.rds$", expr_files)
]

cat("Candidate expression files:\n")
print(expr_files)

if (length(expr_files) == 0) {
  stop("No candidate expression file found in expr_dir. Please check expr_dir.")
}

## 默认选择第一个候选文件；如果不是正确文件，手动指定 expr_file
expr_file <- expr_files[1]

cat("Selected expression file:\n")
print(expr_file)

############################################################
## 4. Load expression data
############################################################

load_expression_file <- function(expr_file) {

  if (grepl("\\.rds$", expr_file, ignore.case = TRUE)) {
    obj <- readRDS(expr_file)
    return(obj)
  }

  if (grepl("\\.rda$|\\.RData$", expr_file, ignore.case = TRUE)) {
    env <- new.env()
    load(expr_file, envir = env)
    obj_names <- ls(env)

    score_obj <- sapply(obj_names, function(nm) {
      x <- get(nm, envir = env)
      if (!(is.matrix(x) || is.data.frame(x))) return(-Inf)
      nr <- nrow(x)
      nc <- ncol(x)
      if (nr < 1000 || nc < 20) return(-Inf)

      score <- 0
      if (gene_use %in% rownames(x)) score <- score + 10
      if (ncol(as.data.frame(x)) >= 2 && gene_use %in% as.character(as.data.frame(x)[[1]])) score <- score + 8
      score <- score + sum(grepl("^TCGA-", colnames(x))) / 10
      score
    })

    score_obj <- sort(score_obj, decreasing = TRUE)

    if (!is.finite(score_obj[1]) || score_obj[1] < 5) {
      stop("Cannot auto-detect expression object in RDA/RData.")
    }

    cat("Selected expression object from RDA:\n")
    print(names(score_obj)[1])

    return(get(names(score_obj)[1], envir = env))
  }

  expr_raw <- data.table::fread(expr_file, data.table = FALSE, check.names = FALSE)
  return(expr_raw)
}

raw_expr <- load_expression_file(expr_file)

expr_df <- as.data.frame(raw_expr, check.names = FALSE)

############################################################
## 5. Prepare expression matrix and extract DLC1
############################################################

## Set rownames
if (!gene_use %in% rownames(expr_df)) {
  first_col <- as.character(expr_df[[1]])

  if (gene_use %in% first_col) {
    rownames(expr_df) <- first_col
    expr_df <- expr_df[, -1, drop = FALSE]
  } else {
    ## Some files use gene symbols like DLC1|ENSG...
    idx <- grep(paste0("^", gene_use, "\\b|^", gene_use, "\\|"), first_col)
    if (length(idx) == 1) {
      rownames(expr_df) <- first_col
      rownames(expr_df)[idx] <- gene_use
      expr_df <- expr_df[, -1, drop = FALSE]
    }
  }
}

if (!gene_use %in% rownames(expr_df)) {
  stop("DLC1 was not found in expression matrix.")
}

## Keep TCGA sample columns
tcga_cols <- colnames(expr_df)[grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-", colnames(expr_df))]

if (length(tcga_cols) >= 20) {
  expr_df <- expr_df[, tcga_cols, drop = FALSE]
}

expr <- as.matrix(expr_df)
mode(expr) <- "numeric"

## Average duplicated genes
if (any(duplicated(rownames(expr)))) {
  expr <- limma::avereps(expr, ID = rownames(expr))
}

## Auto log2 transformation
q99 <- as.numeric(quantile(expr, 0.99, na.rm = TRUE))
maxv <- max(expr, na.rm = TRUE)

cat("Expression q99:", q99, "\n")
cat("Expression max:", maxv, "\n")

if (q99 > 50 || maxv > 100) {
  cat("Expression appears unlogged. Applying log2(x + 1).\n")
  expr <- log2(expr + 1)
} else {
  cat("Expression appears already log-scale. No log transformation applied.\n")
}

sample_info <- data.frame(
  sample = colnames(expr),
  patient = clean_tcga_patient(colnames(expr)),
  sample_type_code = substr(colnames(expr), 14, 15),
  sample_type = get_sample_type(colnames(expr)),
  stringsAsFactors = FALSE
)

write.csv(
  sample_info,
  file.path(outdir, "tables", "01_expression_sample_annotation.csv"),
  row.names = FALSE
)

cat("Sample type distribution:\n")
print(table(sample_info$sample_type))

## Use tumor samples only
tumor_samples <- sample_info %>%
  filter(sample_type == "Tumor") %>%
  pull(sample)

expr_tumor <- expr[, tumor_samples, drop = FALSE]

dlc1_tumor <- data.frame(
  sample = colnames(expr_tumor),
  patient = clean_tcga_patient(colnames(expr_tumor)),
  DLC1 = as.numeric(expr_tumor[gene_use, ]),
  stringsAsFactors = FALSE
) %>%
  group_by(patient) %>%
  summarise(
    DLC1 = mean(DLC1, na.rm = TRUE),
    n_tumor_samples = n(),
    .groups = "drop"
  )

write.csv(
  dlc1_tumor,
  file.path(outdir, "tables", "02_LIHC_tumor_DLC1_expression_by_patient.csv"),
  row.names = FALSE
)

cat("Number of LIHC tumor patients with DLC1 expression:\n")
print(nrow(dlc1_tumor))

saveRDS(
  dlc1_tumor,
  file.path(outdir, "rds", "LIHC_tumor_DLC1_expression_by_patient.rds")
)

############################################################
## 6. Load survival supplemental table
############################################################

if (!file.exists(supp_file)) {
  stop("Survival supplemental file not found: ", supp_file)
}

surv_raw <- data.table::fread(supp_file, data.table = FALSE, check.names = FALSE)

cat("Survival table columns:\n")
print(colnames(surv_raw))

## Find sample/patient column
sample_col_candidates <- c("sample", "Sample", "patient", "bcr_patient_barcode", "submitter_id")
sample_col <- intersect(sample_col_candidates, colnames(surv_raw))[1]

if (is.na(sample_col)) {
  barcode_detect <- sapply(colnames(surv_raw), function(cc) {
    any(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", as.character(surv_raw[[cc]])))
  })
  if (any(barcode_detect)) {
    sample_col <- names(barcode_detect)[which(barcode_detect)[1]]
  }
}

if (is.na(sample_col)) {
  stop("Cannot identify sample/patient column in survival table.")
}

surv_raw$patient <- clean_tcga_patient(surv_raw[[sample_col]])

## Prefer filtering LIHC if cancer type column exists
cancer_cols <- c(
  "cancer type abbreviation", "cancer_type_abbreviation",
  "type", "CancerType", "cohort", "study"
)

cancer_col <- intersect(cancer_cols, colnames(surv_raw))[1]

if (!is.na(cancer_col)) {
  surv_lihc <- surv_raw %>%
    filter(.data[[cancer_col]] == "LIHC")
} else {
  surv_lihc <- surv_raw %>%
    filter(patient %in% dlc1_tumor$patient)
}

surv_lihc <- surv_lihc %>%
  distinct(patient, .keep_all = TRUE)

cat("Number of LIHC survival patients:\n")
print(nrow(surv_lihc))

write.csv(
  surv_lihc,
  file.path(outdir, "tables", "03_LIHC_survival_raw_filtered.csv"),
  row.names = FALSE
)

############################################################
## 7. Merge expression and survival
############################################################

dat0 <- surv_lihc %>%
  left_join(dlc1_tumor, by = "patient") %>%
  filter(!is.na(DLC1))

cat("Merged LIHC patients with survival and DLC1 expression:\n")
print(nrow(dat0))

write.csv(
  dat0,
  file.path(outdir, "tables", "04_LIHC_survival_DLC1_merged_all_endpoints.csv"),
  row.names = FALSE
)

############################################################
## 8. Endpoint configuration
############################################################

endpoint_config <- data.frame(
  endpoint = c("OS", "DSS", "DFI", "PFI"),
  time_col = c("OS.time", "DSS.time", "DFI.time", "PFI.time"),
  status_col = c("OS", "DSS", "DFI", "PFI"),
  endpoint_label = c(
    "Overall survival",
    "Disease-specific survival",
    "Disease-free interval",
    "Progression-free interval"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  endpoint_config,
  file.path(outdir, "tables", "05_endpoint_configuration.csv"),
  row.names = FALSE
)

############################################################
## 9. Main function for optimal cutpoint survival
############################################################

run_endpoint_optimal_cutoff <- function(dat0, endpoint, time_col, status_col, endpoint_label) {

  cat("\n==================================================\n")
  cat("Running endpoint:", endpoint, "\n")
  cat("Time column:", time_col, "\n")
  cat("Status column:", status_col, "\n")
  cat("==================================================\n")

  if (!time_col %in% colnames(dat0) || !status_col %in% colnames(dat0)) {
    message("Endpoint skipped because columns are missing: ", endpoint)
    return(NULL)
  }

  df <- dat0 %>%
    mutate(
      time = suppressWarnings(as.numeric(.data[[time_col]])),
      status = status_to_numeric(.data[[status_col]])
    ) %>%
    filter(
      !is.na(time),
      !is.na(status),
      !is.na(DLC1),
      time > 0,
      status %in% c(0, 1)
    ) %>%
    select(patient, DLC1, time, status, everything())

  df <- df %>%
    distinct(patient, .keep_all = TRUE)

  n_total <- nrow(df)
  n_event <- sum(df$status == 1, na.rm = TRUE)

  cat("N:", n_total, "\n")
  cat("Events:", n_event, "\n")

  if (n_total < 50 || n_event < 10) {
    message("Endpoint skipped due to insufficient sample/event size: ", endpoint)
    return(NULL)
  }

  ## Optimal cutpoint
  cut_res <- tryCatch({
    survminer::surv_cutpoint(
      data = df,
      time = "time",
      event = "status",
      variables = "DLC1",
      minprop = 0.25
    )
  }, error = function(e) {
    message("surv_cutpoint failed for ", endpoint, ": ", e$message)
    NULL
  })

  if (is.null(cut_res)) {
    return(NULL)
  }

  cut_value <- cut_res$cutpoint["DLC1", "cutpoint"]

  df <- df %>%
    mutate(
      DLC1_group = ifelse(DLC1 <= cut_value, "DLC1-low", "DLC1-high"),
      DLC1_group = factor(DLC1_group, levels = c("DLC1-low", "DLC1-high"))
    )

  group_tab <- table(df$DLC1_group)

  cat("Optimal cutoff:", cut_value, "\n")
  print(group_tab)

  if (length(group_tab) < 2 || min(group_tab) < 10) {
    message("Endpoint skipped because one group is too small: ", endpoint)
    return(NULL)
  }

  write.csv(
    df,
    file.path(outdir, "tables", paste0("06_", endpoint, "_grouped_data_optimal_cutoff.csv")),
    row.names = FALSE
  )

  ## KM fit
  surv_obj <- Surv(df$time, df$status)
  fit <- survfit(surv_obj ~ DLC1_group, data = df)

  ## Log-rank p
  logrank <- survdiff(surv_obj ~ DLC1_group, data = df)
  logrank_p <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)

  ## Cox grouped: high vs low
  cox_group <- coxph(surv_obj ~ DLC1_group, data = df)
  cox_group_sum <- summary(cox_group)

  group_HR <- cox_group_sum$coefficients[1, "exp(coef)"]
  group_lower95 <- cox_group_sum$conf.int[1, "lower .95"]
  group_upper95 <- cox_group_sum$conf.int[1, "upper .95"]
  group_p <- cox_group_sum$coefficients[1, "Pr(>|z|)"]

  ## Cox continuous
  cox_cont <- coxph(surv_obj ~ DLC1, data = df)
  cox_cont_sum <- summary(cox_cont)

  cont_HR <- cox_cont_sum$coefficients[1, "exp(coef)"]
  cont_lower95 <- cox_cont_sum$conf.int[1, "lower .95"]
  cont_upper95 <- cox_cont_sum$conf.int[1, "upper .95"]
  cont_p <- cox_cont_sum$coefficients[1, "Pr(>|z|)"]

  ## KM plot
  p_km <- ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    pval.method = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    conf.int = FALSE,
    palette = c("#D73027", "#4575B4"),
    legend.title = "",
    legend.labs = c("DLC1-low", "DLC1-high"),
    xlab = "Time (days)",
    ylab = paste0(endpoint_label, " probability"),
    title = paste0("DLC1 and ", endpoint_label, " in TCGA-LIHC"),
    break.time.by = 1000,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    ggtheme = theme_classic(base_size = 13),
    tables.theme = theme_classic(base_size = 11)
  )

  p_km$plot <- p_km$plot +
    annotate(
      "text",
      x = max(df$time, na.rm = TRUE) * 0.05,
      y = 0.12,
      hjust = 0,
      size = 3.7,
      label = paste0(
        "Cutoff = ", signif(cut_value, 4), "\n",
        "HR(high vs low) = ", sprintf("%.2f", group_HR),
        " (95% CI ", sprintf("%.2f", group_lower95), "-",
        sprintf("%.2f", group_upper95), ")\n",
        "Cox P = ", format_p(group_p)
      )
    ) +
    theme_pub(13) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  pdf(
    file.path(outdir, "plots", paste0("01_", endpoint, "_DLC1_KM_optimal_cutoff_with_risktable.pdf")),
    width = 6.8,
    height = 7.0
  )
  print(p_km)
  dev.off()

  ## Main KM plot only
  safe_ggsave(
    paste0("02_", endpoint, "_DLC1_KM_optimal_cutoff_main_only.pdf"),
    p_km$plot,
    width = 6.5,
    height = 5.5
  )

  ## Cutpoint plot
  p_cut <- plot(cut_res, "DLC1", palette = "npg") +
    theme_pub(13) +
    labs(
      title = paste0(endpoint, " optimal cutpoint for DLC1"),
      x = "DLC1 expression",
      y = "Standardized log-rank statistic"
    )

  safe_ggsave(
    paste0("03_", endpoint, "_DLC1_optimal_cutpoint_plot.pdf"),
    p_cut,
    width = 6.2,
    height = 5
  )

  ## Expression distribution by group
  p_expr <- ggplot(df, aes(x = DLC1_group, y = DLC1, fill = DLC1_group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.45) +
    geom_jitter(width = 0.18, size = 1.1, alpha = 0.45, color = "black") +
    geom_hline(yintercept = cut_value, linetype = "dashed", color = "grey35", linewidth = 0.45) +
    scale_fill_manual(values = c("DLC1-low" = "#D73027", "DLC1-high" = "#4575B4")) +
    theme_pub(13) +
    labs(
      x = NULL,
      y = "DLC1 expression",
      title = paste0(endpoint, " DLC1 groups by optimal cutpoint")
    ) +
    theme(legend.position = "none")

  safe_ggsave(
    paste0("04_", endpoint, "_DLC1_expression_by_optimal_group.pdf"),
    p_expr,
    width = 4.8,
    height = 5.2
  )

  ## Summary table
  res <- data.frame(
    endpoint = endpoint,
    endpoint_label = endpoint_label,
    n = n_total,
    events = n_event,
    cutoff = cut_value,
    n_DLC1_low = as.numeric(group_tab["DLC1-low"]),
    n_DLC1_high = as.numeric(group_tab["DLC1-high"]),
    logrank_p = logrank_p,

    grouped_HR_high_vs_low = group_HR,
    grouped_lower95 = group_lower95,
    grouped_upper95 = group_upper95,
    grouped_cox_p = group_p,

    continuous_HR_per_unit = cont_HR,
    continuous_lower95 = cont_lower95,
    continuous_upper95 = cont_upper95,
    continuous_cox_p = cont_p,
    stringsAsFactors = FALSE
  )

  write.csv(
    res,
    file.path(outdir, "tables", paste0("07_", endpoint, "_summary_optimal_cutoff.csv")),
    row.names = FALSE
  )

  saveRDS(
    list(
      data = df,
      cut_res = cut_res,
      fit = fit,
      cox_group = cox_group,
      cox_cont = cox_cont,
      summary = res
    ),
    file.path(outdir, "rds", paste0(endpoint, "_optimal_cutoff_survival_objects.rds"))
  )

  return(res)
}

############################################################
## 10. Run all endpoints
############################################################

res_list <- list()

for (i in seq_len(nrow(endpoint_config))) {
  ep <- endpoint_config$endpoint[i]
  time_col <- endpoint_config$time_col[i]
  status_col <- endpoint_config$status_col[i]
  endpoint_label <- endpoint_config$endpoint_label[i]

  res_list[[ep]] <- run_endpoint_optimal_cutoff(
    dat0 = dat0,
    endpoint = ep,
    time_col = time_col,
    status_col = status_col,
    endpoint_label = endpoint_label
  )
}

summary_all <- bind_rows(res_list)

write.csv(
  summary_all,
  file.path(outdir, "tables", "08_all_endpoints_DLC1_optimal_cutoff_summary.csv"),
  row.names = FALSE
)

cat("\nAll endpoint summary:\n")
print(summary_all)

############################################################
## 11. Forest plot for four endpoints
############################################################

if (!is.null(summary_all) && nrow(summary_all) > 0) {

  forest_df <- summary_all %>%
    mutate(
      endpoint = factor(endpoint, levels = c("OS", "DSS", "DFI", "PFI")),
      label = paste0(
        endpoint,
        "  HR=", sprintf("%.2f", grouped_HR_high_vs_low),
        " (", sprintf("%.2f", grouped_lower95), "-",
        sprintf("%.2f", grouped_upper95), "), P=", format_p(grouped_cox_p)
      )
    ) %>%
    arrange(endpoint)

  p_forest <- ggplot(
    forest_df,
    aes(
      x = endpoint,
      y = grouped_HR_high_vs_low,
      ymin = grouped_lower95,
      ymax = grouped_upper95
    )
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.45) +
    geom_pointrange(size = 0.65, linewidth = 0.7, color = "#2C3E50") +
    coord_flip() +
    scale_y_log10() +
    theme_pub(13) +
    labs(
      x = NULL,
      y = "Hazard ratio, DLC1-high vs DLC1-low",
      title = "Univariate Cox analysis of DLC1 across LIHC survival endpoints"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  safe_ggsave(
    "05_all_endpoints_DLC1_grouped_Cox_forestplot.pdf",
    p_forest,
    width = 7,
    height = 4.8
  )
}

############################################################
## 12. Save session info
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "tables", "sessionInfo.txt")
)

cat("\nFinished DLC1 LIHC optimal cutpoint survival analysis.\n")
cat("Output directory:\n", outdir, "\n")
cat("\nKey outputs:\n")
cat("tables/08_all_endpoints_DLC1_optimal_cutoff_summary.csv\n")
cat("plots/01_<endpoint>_DLC1_KM_optimal_cutoff_with_risktable.pdf\n")
cat("plots/02_<endpoint>_DLC1_KM_optimal_cutoff_main_only.pdf\n")
cat("plots/03_<endpoint>_DLC1_optimal_cutpoint_plot.pdf\n")
cat("plots/05_all_endpoints_DLC1_grouped_Cox_forestplot.pdf\n")












############################################################
## DLC1 TCGA-LIHC bulk analysis final version
## Expression + clinical + DEG + GO/KEGG/GSEA
############################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

############################################################
## 0. User paths
############################################################

data_file <- file.path(hcc_dlc1_root(), "bulk_lihc", "lihc.gdc_2022.rda")
outdir <- file.path(hcc_dlc1_root(), "DLC1_LIHC_bulk_final_noOptimalCutoff")

gene_use <- "DLC1"

## DEG grouping method
## For DEG/enrichment, median grouping is more stable than survival-optimized cutoff.
group_method <- "median"

## DEG thresholds
logfc_cut <- 1
fdr_cut <- 0.05

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "rds"), showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Packages
############################################################

cran_pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "ggpubr", "stringr",
  "tibble", "patchwork", "ggrepel", "pheatmap",
  "scales", "RColorBrewer", "msigdbr"
)

bioc_pkgs <- c(
  "limma", "edgeR", "clusterProfiler", "enrichplot",
  "org.Hs.eg.db", "DOSE"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(stringr)
  library(tibble)
  library(patchwork)
  library(ggrepel)
  library(pheatmap)
  library(scales)
  library(RColorBrewer)
  library(msigdbr)
  library(limma)
  library(edgeR)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(DOSE)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise

############################################################
## 2. Plot theme and helper functions
############################################################

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.45, color = "black"),
      axis.ticks = element_line(linewidth = 0.45, color = "black"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black")
    )
}

safe_ggsave <- function(filename, plot, width = 6, height = 5) {
  outfile <- file.path(outdir, "plots", filename)
  tryCatch({
    ggsave(outfile, plot, width = width, height = height, device = cairo_pdf)
  }, error = function(e) {
    message("cairo_pdf failed, using default pdf for: ", filename)
    ggsave(outfile, plot, width = width, height = height)
  })
}

clean_tcga_barcode <- function(x, n = 12) {
  substr(as.character(x), 1, n)
}

get_sample_type <- function(barcode) {
  code <- substr(barcode, 14, 15)
  ifelse(
    code %in% c("01"),
    "Tumor",
    ifelse(code %in% c("11"), "Normal", "Other")
  )
}

make_valid_filename <- function(x) {
  x <- gsub("[/\\:*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

############################################################
## 3. Load RDA and auto-detect objects
############################################################

env <- new.env()
load(data_file, envir = env)

obj_names <- ls(env)
cat("Objects in RDA:\n")
print(obj_names)

is_matrix_like <- function(x) {
  is.matrix(x) || is.data.frame(x)
}

score_expression_candidate <- function(obj, gene) {
  if (!is_matrix_like(obj)) return(-Inf)

  x <- as.data.frame(obj)
  nr <- nrow(x)
  nc <- ncol(x)
  if (nr < 1000 || nc < 20) return(-Inf)

  rn <- rownames(x)
  cn <- colnames(x)

  score <- 0

  if (!is.null(rn) && gene %in% rn) score <- score + 10

  if (ncol(x) >= 2) {
    first_col <- as.character(x[[1]])
    if (gene %in% first_col) score <- score + 8
  }

  tcga_col_n <- sum(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-", cn))
  score <- score + min(tcga_col_n, 100) / 10

  numeric_cols <- suppressWarnings(sum(sapply(x[seq_len(min(nc, 20))], function(v) {
    mean(!is.na(as.numeric(as.character(v)))) > 0.8
  })))
  score <- score + numeric_cols / 5

  score
}

expr_scores <- sapply(obj_names, function(nm) {
  score_expression_candidate(get(nm, envir = env), gene_use)
})

expr_scores <- sort(expr_scores, decreasing = TRUE)
cat("Expression candidate scores:\n")
print(expr_scores)

expr_object_name <- names(expr_scores)[1]
if (!is.finite(expr_scores[1]) || expr_scores[1] < 5) {
  stop("Cannot auto-detect expression matrix. Please check objects in RDA.")
}

cat("Selected expression object:\n")
print(expr_object_name)

raw_expr <- get(expr_object_name, envir = env)

score_clinical_candidate <- function(obj) {
  if (!is.data.frame(obj)) return(-Inf)
  if (nrow(obj) < 50 || ncol(obj) < 3) return(-Inf)

  cn <- colnames(obj)
  score <- 0

  id_cols <- c("patient", "submitter_id", "bcr_patient_barcode",
               "case_submitter_id", "barcode", "sample", "Sample")
  if (any(id_cols %in% cn)) score <- score + 5

  clin_cols <- c("stage", "ajcc_pathologic_stage", "pathologic_stage",
                 "grade", "tumor_grade", "gender", "sex",
                 "age", "age_at_diagnosis", "vital_status")
  score <- score + sum(clin_cols %in% cn)

  any_tcga <- any(sapply(obj, function(v) {
    any(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", as.character(v)))
  }))
  if (any_tcga) score <- score + 5

  score
}

clinical_scores <- sapply(obj_names, function(nm) {
  score_clinical_candidate(get(nm, envir = env))
})

clinical_scores <- sort(clinical_scores, decreasing = TRUE)
cat("Clinical candidate scores:\n")
print(clinical_scores)

cli <- NULL
clinical_object_name <- NA_character_

if (is.finite(clinical_scores[1]) && clinical_scores[1] >= 3) {
  clinical_object_name <- names(clinical_scores)[1]
  cli <- as.data.frame(get(clinical_object_name, envir = env))
  cat("Selected clinical object:\n")
  print(clinical_object_name)
} else {
  message("No reliable clinical object detected. Clinical stratification will be skipped.")
}

############################################################
## 4. Prepare expression matrix
############################################################

expr_df <- as.data.frame(raw_expr, check.names = FALSE)

if (!gene_use %in% rownames(expr_df)) {
  first_col <- as.character(expr_df[[1]])
  if (gene_use %in% first_col) {
    rownames(expr_df) <- first_col
    expr_df <- expr_df[, -1, drop = FALSE]
  }
}

if (!gene_use %in% rownames(expr_df)) {
  stop("DLC1 not found in expression matrix rownames or first column.")
}

## Keep TCGA-like sample columns only if possible
tcga_cols <- colnames(expr_df)[grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-", colnames(expr_df))]
if (length(tcga_cols) >= 20) {
  expr_df <- expr_df[, tcga_cols, drop = FALSE]
}

expr <- as.matrix(expr_df)
mode(expr) <- "numeric"

## Remove genes with all NA
expr <- expr[rowSums(!is.na(expr)) > 0, , drop = FALSE]

## Replace NA with row median
if (any(is.na(expr))) {
  for (i in seq_len(nrow(expr))) {
    if (any(is.na(expr[i, ]))) {
      expr[i, is.na(expr[i, ])] <- median(expr[i, ], na.rm = TRUE)
    }
  }
}

## Average duplicated gene symbols
if (any(duplicated(rownames(expr)))) {
  expr <- limma::avereps(expr, ID = rownames(expr))
}

## Auto log2 transform if needed
q99 <- as.numeric(quantile(expr, 0.99, na.rm = TRUE))
maxv <- max(expr, na.rm = TRUE)

cat("Expression q99:", q99, "\n")
cat("Expression max:", maxv, "\n")

if (q99 > 50 || maxv > 100) {
  cat("Expression matrix appears unlogged. Applying log2(x + 1).\n")
  expr <- log2(expr + 1)
} else {
  cat("Expression matrix appears already log-scale. No log transform applied.\n")
}

saveRDS(expr, file.path(outdir, "rds", "01_expression_matrix_logscale.rds"))

cat("Final expression matrix dimension:\n")
print(dim(expr))

############################################################
## 5. TCGA sample annotation
############################################################

sample_info <- data.frame(
  sample = colnames(expr),
  patient = clean_tcga_barcode(colnames(expr), 12),
  sample_type_code = substr(colnames(expr), 14, 15),
  sample_type = get_sample_type(colnames(expr)),
  stringsAsFactors = FALSE
)

write.csv(
  sample_info,
  file.path(outdir, "tables", "01_sample_annotation.csv"),
  row.names = FALSE
)

cat("Sample type distribution:\n")
print(table(sample_info$sample_type))

############################################################
## 6. DLC1 expression extraction
############################################################

dlc1_df <- data.frame(
  sample = colnames(expr),
  patient = clean_tcga_barcode(colnames(expr), 12),
  DLC1 = as.numeric(expr[gene_use, ]),
  stringsAsFactors = FALSE
) %>%
  left_join(sample_info, by = c("sample", "patient"))

write.csv(
  dlc1_df,
  file.path(outdir, "tables", "02_DLC1_expression_all_samples.csv"),
  row.names = FALSE
)

############################################################
## 7. Tumor vs normal DLC1 expression
############################################################

expr_tn <- dlc1_df %>%
  filter(sample_type %in% c("Tumor", "Normal")) %>%
  mutate(sample_type = factor(sample_type, levels = c("Normal", "Tumor")))

if (length(unique(expr_tn$sample_type)) == 2) {
  p_tn <- ggplot(expr_tn, aes(x = sample_type, y = DLC1, fill = sample_type)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.45) +
    geom_jitter(width = 0.18, size = 1.25, alpha = 0.45, color = "black") +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 4.2) +
    scale_fill_manual(values = c("Normal" = "#4DBBD5", "Tumor" = "#E64B35")) +
    theme_pub(14) +
    labs(
      x = NULL,
      y = "DLC1 expression",
      title = "DLC1 expression in TCGA-LIHC"
    ) +
    theme(legend.position = "none")

  safe_ggsave("01_LIHC_DLC1_tumor_vs_normal_boxplot.pdf", p_tn, width = 4.8, height = 5.2)

  tn_stat <- wilcox.test(DLC1 ~ sample_type, data = expr_tn)
  write.csv(
    data.frame(
      comparison = "Tumor_vs_Normal",
      n_normal = sum(expr_tn$sample_type == "Normal"),
      n_tumor = sum(expr_tn$sample_type == "Tumor"),
      p_value = tn_stat$p.value
    ),
    file.path(outdir, "tables", "03_LIHC_DLC1_tumor_vs_normal_stat.csv"),
    row.names = FALSE
  )
}

############################################################
## 8. Paired tumor-normal DLC1 expression
############################################################

paired_patients <- expr_tn %>%
  group_by(patient) %>%
  summarise(
    has_tumor = any(sample_type == "Tumor"),
    has_normal = any(sample_type == "Normal"),
    .groups = "drop"
  ) %>%
  filter(has_tumor & has_normal) %>%
  pull(patient)

paired_df <- expr_tn %>%
  filter(patient %in% paired_patients) %>%
  group_by(patient, sample_type) %>%
  summarise(DLC1 = mean(DLC1, na.rm = TRUE), .groups = "drop") %>%
  mutate(sample_type = factor(sample_type, levels = c("Normal", "Tumor")))

write.csv(
  paired_df,
  file.path(outdir, "tables", "04_LIHC_DLC1_paired_expression.csv"),
  row.names = FALSE
)

if (length(paired_patients) >= 3) {
  p_paired <- ggplot(paired_df, aes(x = sample_type, y = DLC1, group = patient)) +
    geom_line(color = "grey55", alpha = 0.65, linewidth = 0.45) +
    geom_point(aes(fill = sample_type), shape = 21, size = 2.6, color = "black", stroke = 0.28) +
    stat_compare_means(paired = TRUE, method = "wilcox.test", label = "p.format", size = 4.2) +
    scale_fill_manual(values = c("Normal" = "#4DBBD5", "Tumor" = "#E64B35")) +
    theme_pub(14) +
    labs(
      x = NULL,
      y = "DLC1 expression",
      title = "Paired DLC1 expression in TCGA-LIHC"
    ) +
    theme(legend.position = "none")

  safe_ggsave("02_LIHC_DLC1_paired_tumor_normal.pdf", p_paired, width = 4.8, height = 5.2)
}

############################################################
## 9. Clinical data preparation and clinical stratification
############################################################

cli_dlc1 <- NULL

if (!is.null(cli)) {

  cli <- as.data.frame(cli, check.names = FALSE)

  id_candidates <- c(
    "patient", "submitter_id", "bcr_patient_barcode",
    "case_submitter_id", "barcode", "sample", "Sample"
  )

  id_col <- intersect(id_candidates, colnames(cli))[1]

  if (is.na(id_col)) {
    ## Try to detect TCGA barcode column
    barcode_detect <- sapply(colnames(cli), function(cc) {
      any(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", as.character(cli[[cc]])))
    })
    if (any(barcode_detect)) id_col <- names(barcode_detect)[which(barcode_detect)[1]]
  }

  if (!is.na(id_col)) {
    cli$patient <- clean_tcga_barcode(cli[[id_col]], 12)

    tumor_dlc1 <- dlc1_df %>%
      filter(sample_type == "Tumor") %>%
      group_by(patient) %>%
      summarise(DLC1 = mean(DLC1, na.rm = TRUE), .groups = "drop")

    cli_dlc1 <- cli %>%
      left_join(tumor_dlc1, by = "patient") %>%
      filter(!is.na(DLC1))

    ## Clean stage
    stage_candidates <- c("ajcc_pathologic_stage", "pathologic_stage", "stage", "tumor_stage")
    stage_col <- intersect(stage_candidates, colnames(cli_dlc1))[1]
    if (!is.na(stage_col)) {
      cli_dlc1$clinical_stage <- as.character(cli_dlc1[[stage_col]])
      cli_dlc1$clinical_stage <- str_replace_all(cli_dlc1$clinical_stage, "stage ", "Stage ")
      cli_dlc1$clinical_stage <- str_replace_all(cli_dlc1$clinical_stage, "Stage i", "Stage I")
      cli_dlc1$clinical_stage <- str_replace_all(cli_dlc1$clinical_stage, "Stage ii", "Stage II")
      cli_dlc1$clinical_stage <- str_replace_all(cli_dlc1$clinical_stage, "Stage iii", "Stage III")
      cli_dlc1$clinical_stage <- str_replace_all(cli_dlc1$clinical_stage, "Stage iv", "Stage IV")
      cli_dlc1$clinical_stage[grepl("not|unknown|NA", cli_dlc1$clinical_stage, ignore.case = TRUE)] <- NA
    }

    ## Clean grade
    grade_candidates <- c("grade", "tumor_grade", "histological_grade", "neoplasm_histologic_grade")
    grade_col <- intersect(grade_candidates, colnames(cli_dlc1))[1]
    if (!is.na(grade_col)) {
      cli_dlc1$clinical_grade <- as.character(cli_dlc1[[grade_col]])
      cli_dlc1$clinical_grade[grepl("not|unknown|NA", cli_dlc1$clinical_grade, ignore.case = TRUE)] <- NA
    }

    ## Clean gender
    gender_candidates <- c("gender", "sex")
    gender_col <- intersect(gender_candidates, colnames(cli_dlc1))[1]
    if (!is.na(gender_col)) {
      cli_dlc1$clinical_gender <- as.character(cli_dlc1[[gender_col]])
      cli_dlc1$clinical_gender[grepl("not|unknown|NA", cli_dlc1$clinical_gender, ignore.case = TRUE)] <- NA
    }

    ## Clean age
    age_candidates <- c("age_at_diagnosis", "age", "days_to_birth")
    age_col <- intersect(age_candidates, colnames(cli_dlc1))[1]
    if (!is.na(age_col)) {
      age_raw <- suppressWarnings(as.numeric(cli_dlc1[[age_col]]))
      if (median(abs(age_raw), na.rm = TRUE) > 150) {
        age_year <- abs(age_raw) / 365.25
      } else {
        age_year <- age_raw
      }
      cli_dlc1$clinical_age_year <- age_year
      cli_dlc1$clinical_age_group <- ifelse(age_year >= 60, ">=60", "<60")
    }

    write.csv(
      cli_dlc1,
      file.path(outdir, "tables", "05_LIHC_clinical_with_DLC1.csv"),
      row.names = FALSE
    )
  }
}

plot_clinical_box <- function(df, xvar, filename, title) {
  df2 <- df %>%
    filter(!is.na(.data[[xvar]]), !is.na(DLC1)) %>%
    mutate(x_group = as.factor(.data[[xvar]]))

  if (n_distinct(df2$x_group) < 2) return(NULL)

  p <- ggplot(df2, aes(x = x_group, y = DLC1, fill = x_group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.45) +
    geom_jitter(width = 0.18, size = 1.05, alpha = 0.42, color = "black") +
    stat_compare_means(method = "kruskal.test", label = "p.format", size = 4.0) +
    theme_pub(13) +
    labs(x = NULL, y = "DLC1 expression", title = title) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  safe_ggsave(filename, p, width = 5.8, height = 5.2)
}

if (!is.null(cli_dlc1)) {
  if ("clinical_stage" %in% colnames(cli_dlc1)) {
    plot_clinical_box(cli_dlc1, "clinical_stage", "03_LIHC_DLC1_by_clinical_stage.pdf", "DLC1 expression by pathological stage")
  }

  if ("clinical_grade" %in% colnames(cli_dlc1)) {
    plot_clinical_box(cli_dlc1, "clinical_grade", "04_LIHC_DLC1_by_clinical_grade.pdf", "DLC1 expression by tumor grade")
  }

  if ("clinical_gender" %in% colnames(cli_dlc1)) {
    plot_clinical_box(cli_dlc1, "clinical_gender", "05_LIHC_DLC1_by_gender.pdf", "DLC1 expression by gender")
  }

  if ("clinical_age_group" %in% colnames(cli_dlc1)) {
    plot_clinical_box(cli_dlc1, "clinical_age_group", "06_LIHC_DLC1_by_age_group.pdf", "DLC1 expression by age group")
  }
}

############################################################
## 10. DLC1 high/low grouping for DEG
############################################################

tumor_samples <- sample_info %>%
  filter(sample_type == "Tumor") %>%
  pull(sample)

expr_tumor <- expr[, tumor_samples, drop = FALSE]

## Remove genes with zero variance
gene_var <- apply(expr_tumor, 1, var, na.rm = TRUE)
expr_tumor <- expr_tumor[is.finite(gene_var) & gene_var > 0, , drop = FALSE]

if (!gene_use %in% rownames(expr_tumor)) {
  stop("DLC1 not found in tumor expression matrix after filtering.")
}

dlc1_tumor <- as.numeric(expr_tumor[gene_use, ])
names(dlc1_tumor) <- colnames(expr_tumor)

cut_value <- median(dlc1_tumor, na.rm = TRUE)

group <- ifelse(dlc1_tumor >= cut_value, "DLC1_high", "DLC1_low")
group <- factor(group, levels = c("DLC1_low", "DLC1_high"))

group_df <- data.frame(
  sample = colnames(expr_tumor),
  patient = clean_tcga_barcode(colnames(expr_tumor), 12),
  DLC1 = dlc1_tumor,
  group = group,
  stringsAsFactors = FALSE
)

write.csv(
  group_df,
  file.path(outdir, "tables", "06_LIHC_DLC1_high_low_group_median.csv"),
  row.names = FALSE
)

cat("DLC1 grouping:\n")
print(table(group))

############################################################
## 11. DEG analysis: DLC1-low vs DLC1-high
############################################################

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_tumor, design)
contrast.matrix <- makeContrasts(
  DLC1_low_vs_high = DLC1_low - DLC1_high,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg <- topTable(
  fit2,
  coef = "DLC1_low_vs_high",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

deg <- deg %>%
  rownames_to_column("gene") %>%
  arrange(adj.P.Val)

deg$change <- "Not significant"
deg$change[deg$adj.P.Val < fdr_cut & deg$logFC > logfc_cut] <- "Higher in DLC1-low"
deg$change[deg$adj.P.Val < fdr_cut & deg$logFC < -logfc_cut] <- "Higher in DLC1-high"

write.csv(
  deg,
  file.path(outdir, "tables", "07_LIHC_DLC1_low_vs_high_DEG_all.csv"),
  row.names = FALSE
)

write.csv(
  deg %>% filter(change == "Higher in DLC1-low"),
  file.path(outdir, "tables", "08_LIHC_DLC1_low_vs_high_DEG_higher_in_DLC1_low.csv"),
  row.names = FALSE
)

write.csv(
  deg %>% filter(change == "Higher in DLC1-high"),
  file.path(outdir, "tables", "09_LIHC_DLC1_low_vs_high_DEG_higher_in_DLC1_high.csv"),
  row.names = FALSE
)

############################################################
## 12. Publication-style volcano plot
############################################################

deg$neglog10_adjP <- -log10(deg$adj.P.Val + 1e-300)

top_up <- deg %>%
  filter(change == "Higher in DLC1-low") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 15)

top_down <- deg %>%
  filter(change == "Higher in DLC1-high") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 15)

top_label <- bind_rows(top_up, top_down) %>%
  distinct(gene, .keep_all = TRUE)

y_max_show <- 20
x_max_show <- max(abs(deg$logFC), na.rm = TRUE)
x_max_show <- max(2.5, ceiling(x_max_show * 10) / 10)

deg$neglog10_adjP_plot <- pmin(deg$neglog10_adjP, y_max_show)
top_label$neglog10_adjP_plot <- pmin(top_label$neglog10_adjP, y_max_show)

p_volcano <- ggplot() +
  geom_point(
    data = deg %>% filter(change == "Not significant"),
    aes(x = logFC, y = neglog10_adjP_plot),
    color = "grey78",
    size = 1.0,
    alpha = 0.60
  ) +
  geom_point(
    data = deg %>% filter(change == "Higher in DLC1-low"),
    aes(x = logFC, y = neglog10_adjP_plot, color = change),
    size = 1.6,
    alpha = 0.85
  ) +
  geom_point(
    data = deg %>% filter(change == "Higher in DLC1-high"),
    aes(x = logFC, y = neglog10_adjP_plot, color = change),
    size = 1.6,
    alpha = 0.85
  ) +
  geom_vline(xintercept = c(-logfc_cut, logfc_cut), linetype = "dashed", color = "grey45", linewidth = 0.7) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey45", linewidth = 0.7) +
  ggrepel::geom_text_repel(
    data = top_label,
    aes(x = logFC, y = neglog10_adjP_plot, label = gene),
    size = 4.0,
    color = "black",
    min.segment.length = 0,
    segment.color = "black",
    segment.linewidth = 0.35,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Higher in DLC1-low" = "#D73027",
      "Higher in DLC1-high" = "#4575B4"
    )
  ) +
  scale_x_continuous(
    limits = c(-x_max_show, x_max_show),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(
    limits = c(0, y_max_show),
    breaks = seq(0, y_max_show, by = 5),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = expression(log[2]~"fold change (DLC1-low vs DLC1-high)"),
    y = expression(-log[10]~"FDR"),
    title = "Differentially expressed genes associated with DLC1 expression",
    color = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 12),
    axis.title = element_text(color = "black", size = 17),
    axis.text = element_text(color = "black", size = 14),
    axis.line = element_line(linewidth = 1.0, color = "black"),
    axis.ticks = element_line(linewidth = 0.9, color = "black"),
    axis.ticks.length = unit(0.18, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
    plot.margin = margin(12, 18, 10, 10)
  )

safe_ggsave("07_LIHC_DLC1_low_vs_high_DEG_volcano_pub.pdf", p_volcano, width = 7.2, height = 6.4)
ggsave(
  file.path(outdir, "plots", "07_LIHC_DLC1_low_vs_high_DEG_volcano_pub.png"),
  p_volcano,
  width = 7.2,
  height = 6.4,
  dpi = 600,
  bg = "white"
)

############################################################
## 13. Top DEG heatmap
############################################################

top_heat_genes <- bind_rows(
  deg %>% filter(change == "Higher in DLC1-low") %>% arrange(adj.P.Val) %>% slice_head(n = 25),
  deg %>% filter(change == "Higher in DLC1-high") %>% arrange(adj.P.Val) %>% slice_head(n = 25)
) %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(gene)

top_heat_genes <- top_heat_genes[top_heat_genes %in% rownames(expr_tumor)]

if (length(top_heat_genes) >= 5) {
  mat_top <- expr_tumor[top_heat_genes, , drop = FALSE]
  mat_top_z <- t(scale(t(mat_top)))
  mat_top_z[mat_top_z > 2] <- 2
  mat_top_z[mat_top_z < -2] <- -2
  mat_top_z[is.na(mat_top_z)] <- 0

  anno_col <- data.frame(
    DLC1_group = group
  )
  rownames(anno_col) <- colnames(expr_tumor)

  anno_col <- anno_col[order(anno_col$DLC1_group), , drop = FALSE]
  mat_top_z <- mat_top_z[, rownames(anno_col), drop = FALSE]

  ann_colors <- list(
    DLC1_group = c(
      "DLC1_low" = "#D73027",
      "DLC1_high" = "#4575B4"
    )
  )

  pdf(file.path(outdir, "plots", "08_LIHC_DLC1_low_vs_high_top_DEG_heatmap_pub.pdf"), width = 8.5, height = 9.5)
  pheatmap(
    mat_top_z,
    annotation_col = anno_col,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    fontsize_col = 6,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
    border_color = NA,
    main = "Top DEGs associated with DLC1 expression"
  )
  dev.off()
}

############################################################
## 14. GO / KEGG enrichment
############################################################

genes_low <- deg %>%
  filter(change == "Higher in DLC1-low") %>%
  pull(gene)

genes_high <- deg %>%
  filter(change == "Higher in DLC1-high") %>%
  pull(gene)

map_entrez <- function(genes) {
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  if (length(genes) < 5) return(NULL)

  eg <- tryCatch({
    suppressMessages(
      bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )
  }, error = function(e) NULL)

  if (is.null(eg) || nrow(eg) < 5) return(NULL)
  eg <- eg[!duplicated(eg$ENTREZID), , drop = FALSE]
  eg
}

run_go <- function(genes, label) {
  eg <- map_entrez(genes)
  if (is.null(eg)) return(NULL)

  ego <- tryCatch({
    enrichGO(
      gene = eg$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      readable = TRUE
    )
  }, error = function(e) NULL)

  if (is.null(ego)) return(NULL)
  out <- as.data.frame(ego)
  if (nrow(out) == 0) return(NULL)
  out$direction <- label
  out
}

run_kegg <- function(genes, label) {
  eg <- map_entrez(genes)
  if (is.null(eg)) return(NULL)

  ekk <- tryCatch({
    enrichKEGG(
      gene = eg$ENTREZID,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
  }, error = function(e) NULL)

  if (is.null(ekk)) return(NULL)
  out <- as.data.frame(ekk)
  if (nrow(out) == 0) return(NULL)
  out$direction <- label
  out
}

go_low <- run_go(genes_low, "Higher in DLC1-low")
go_high <- run_go(genes_high, "Higher in DLC1-high")
go_all <- bind_rows(go_low, go_high)

if (!is.null(go_all) && nrow(go_all) > 0) {
  write.csv(
    go_all,
    file.path(outdir, "tables", "10_LIHC_DLC1_low_vs_high_GO_BP.csv"),
    row.names = FALSE
  )

  go_top <- go_all %>%
    group_by(direction) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(
      Description = stringr::str_wrap(Description, width = 48),
      neglog10 = -log10(p.adjust + 1e-300),
      direction = factor(direction, levels = c("Higher in DLC1-low", "Higher in DLC1-high"))
    )

  p_go <- ggplot(go_top, aes(x = direction, y = reorder(Description, neglog10))) +
    geom_point(aes(size = Count, color = neglog10), alpha = 0.9) +
    scale_color_gradient(low = "#FEE0D2", high = "#CB181D") +
    theme_pub(12) +
    labs(
      x = NULL,
      y = NULL,
      size = "Gene count",
      color = "-log10(FDR)",
      title = "GO BP enrichment of DLC1-associated DEGs"
    ) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "right"
    )

  safe_ggsave("09_LIHC_DLC1_low_vs_high_GO_BP_dotplot_pub.pdf", p_go, width = 9.2, height = 7.5)
}

kegg_low <- run_kegg(genes_low, "Higher in DLC1-low")
kegg_high <- run_kegg(genes_high, "Higher in DLC1-high")
kegg_all <- bind_rows(kegg_low, kegg_high)

if (!is.null(kegg_all) && nrow(kegg_all) > 0) {
  write.csv(
    kegg_all,
    file.path(outdir, "tables", "11_LIHC_DLC1_low_vs_high_KEGG.csv"),
    row.names = FALSE
  )

  kegg_top <- kegg_all %>%
    group_by(direction) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(
      Description = stringr::str_wrap(Description, width = 48),
      neglog10 = -log10(p.adjust + 1e-300),
      direction = factor(direction, levels = c("Higher in DLC1-low", "Higher in DLC1-high"))
    )

  p_kegg <- ggplot(kegg_top, aes(x = direction, y = reorder(Description, neglog10))) +
    geom_point(aes(size = Count, color = neglog10), alpha = 0.9) +
    scale_color_gradient(low = "#FEE0D2", high = "#CB181D") +
    theme_pub(12) +
    labs(
      x = NULL,
      y = NULL,
      size = "Gene count",
      color = "-log10(FDR)",
      title = "KEGG enrichment of DLC1-associated DEGs"
    ) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "right"
    )

  safe_ggsave("10_LIHC_DLC1_low_vs_high_KEGG_dotplot_pub.pdf", p_kegg, width = 9.2, height = 7.0)
}

############################################################
## 15. GSEA analysis
############################################################

gene_list_df <- deg %>%
  filter(!is.na(logFC)) %>%
  group_by(gene) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

gene_list <- gene_list_df$logFC
names(gene_list) <- gene_list_df$gene
gene_list <- sort(gene_list, decreasing = TRUE)

get_msigdbr_kegg <- function() {
  fm <- names(formals(msigdbr::msigdbr))

  if ("collection" %in% fm) {
    m_df <- msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = "C2",
      subcollection = "CP:KEGG"
    )
  } else {
    m_df <- msigdbr::msigdbr(
      species = "Homo sapiens",
      category = "C2",
      subcategory = "CP:KEGG"
    )
  }

  m_df %>%
    select(gs_name, gene_symbol) %>%
    distinct()
}

term2gene <- get_msigdbr_kegg()

gsea_res <- tryCatch({
  GSEA(
    geneList = gene_list,
    TERM2GENE = term2gene,
    pvalueCutoff = 0.25,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
}, error = function(e) {
  message("GSEA failed: ", e$message)
  NULL
})

if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
  gsea_df <- as.data.frame(gsea_res)

  write.csv(
    gsea_df,
    file.path(outdir, "tables", "12_LIHC_DLC1_low_vs_high_GSEA_KEGG.csv"),
    row.names = FALSE
  )

  p_gsea_dot <- dotplot(gsea_res, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle("GSEA KEGG: DLC1-low vs DLC1-high") +
    theme_pub(12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9)
    )

  safe_ggsave("11_LIHC_DLC1_low_vs_high_GSEA_KEGG_dotplot_pub.pdf", p_gsea_dot, width = 10, height = 7.5)

  key_pathways <- c(
    "KEGG_FOCAL_ADHESION",
    "KEGG_ECM_RECEPTOR_INTERACTION",
    "KEGG_PI3K_AKT_SIGNALING_PATHWAY",
    "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON",
    "KEGG_CELL_ADHESION_MOLECULES_CAMS",
    "KEGG_PATHWAYS_IN_CANCER"
  )

  for (pw in key_pathways) {
    if (pw %in% gsea_res@result$ID) {
      p_pw <- gseaplot2(
        gsea_res,
        geneSetID = pw,
        title = pw,
        base_size = 13,
        color = "#D73027"
      )

      safe_ggsave(
        paste0("12_GSEA_", make_valid_filename(pw), "_pub.pdf"),
        p_pw,
        width = 7.2,
        height = 5.6
      )
    }
  }
} else {
  message("No significant or available GSEA results.")
}

############################################################
## 16. Save summary and session info
############################################################

summary_info <- data.frame(
  item = c(
    "expression_object",
    "clinical_object",
    "n_total_samples",
    "n_tumor_samples",
    "n_normal_samples",
    "DLC1_median_cutoff_for_DEG",
    "DEG_logFC_cutoff",
    "DEG_FDR_cutoff",
    "n_higher_in_DLC1_low",
    "n_higher_in_DLC1_high"
  ),
  value = c(
    expr_object_name,
    ifelse(is.na(clinical_object_name), "Not detected", clinical_object_name),
    ncol(expr),
    sum(sample_info$sample_type == "Tumor"),
    sum(sample_info$sample_type == "Normal"),
    cut_value,
    logfc_cut,
    fdr_cut,
    sum(deg$change == "Higher in DLC1-low"),
    sum(deg$change == "Higher in DLC1-high")
  )
)

write.csv(
  summary_info,
  file.path(outdir, "tables", "00_analysis_summary.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "tables", "sessionInfo.txt")
)

cat("\nFinished DLC1 LIHC bulk analysis final version without optimal cutoff survival.\n")
cat("Output directory:\n", outdir, "\n")
cat("\nKey plots:\n")
cat("01_LIHC_DLC1_tumor_vs_normal_boxplot.pdf\n")
cat("02_LIHC_DLC1_paired_tumor_normal.pdf\n")
cat("07_LIHC_DLC1_low_vs_high_DEG_volcano_pub.pdf\n")
cat("08_LIHC_DLC1_low_vs_high_top_DEG_heatmap_pub.pdf\n")
cat("09_LIHC_DLC1_low_vs_high_GO_BP_dotplot_pub.pdf\n")
cat("10_LIHC_DLC1_low_vs_high_KEGG_dotplot_pub.pdf\n")
cat("11_LIHC_DLC1_low_vs_high_GSEA_KEGG_dotplot_pub.pdf\n")















########################################################
##单细胞分析脚本
########################################################

##UMAP图
########################################################
rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "Seurat",
  "Matrix",
  "data.table",
  "dplyr",
  "ggplot2",
  "patchwork",
  "ggrepel",
  "scales"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

########################################################
## 1. paths
########################################################
base_dir <- file.path(hcc_dlc1_root(), "single_cell")
out_dir  <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

########################################################
## 2. helper functions
########################################################
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

save_pdf <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

########################################################
## 3. read GSE149614
########################################################
meta_file  <- list.files(base_dir, pattern = "GSE149614.*metadata.*txt.gz$", full.names = TRUE)
count_file <- list.files(base_dir, pattern = "GSE149614.*count.*txt.gz$", full.names = TRUE)

if (length(meta_file) == 0) stop("Metadata file for GSE149614 not found.")
if (length(count_file) == 0) stop("Count file for GSE149614 not found.")

meta149 <- fread(meta_file[1], data.table = FALSE)
mat149_df <- fread(count_file[1], data.table = FALSE, check.names = FALSE)

genes149 <- mat149_df[[1]]
mat149 <- as.matrix(mat149_df[, -1, drop = FALSE])
rownames(mat149) <- genes149
storage.mode(mat149) <- "numeric"
mat149 <- as(mat149, "dgCMatrix")

## merge duplicated genes
if (any(duplicated(rownames(mat149)))) {
  gene_levels <- unique(rownames(mat149))
  group_factor <- factor(rownames(mat149), levels = gene_levels)
  model_mat <- sparse.model.matrix(~ 0 + group_factor)
  colnames(model_mat) <- sub("^group_factor", "", colnames(model_mat))
  mat149 <- t(model_mat) %*% mat149
  mat149 <- as(mat149, "dgCMatrix")
}

meta_key <- pick_col(meta149, c("cell","Cell","barcode","Barcode","cell_id","CellID"))
if (is.null(meta_key)) meta_key <- colnames(meta149)[1]
rownames(meta149) <- meta149[[meta_key]]

common_cells <- intersect(colnames(mat149), rownames(meta149))
mat149 <- mat149[, common_cells, drop = FALSE]
meta149 <- meta149[common_cells, , drop = FALSE]

obj <- CreateSeuratObject(
  counts = mat149,
  meta.data = meta149,
  project = "GSE149614",
  min.cells = 3,
  min.features = 200
)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

saveRDS(obj, file.path(out_dir, "rds", "00_raw_obj.rds"))

########################################################
## 4. QC filtering
########################################################
p_qc_before <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
) + theme_classic(base_size = 13)

save_pdf(
  p_qc_before,
  file.path(out_dir, "plots", "01_QC_before_filter.pdf"),
  width = 9,
  height = 4
)

obj <- subset(
  obj,
  subset = nFeature_RNA >= 300 &
    nFeature_RNA <= 7000 &
    nCount_RNA >= 500 &
    percent.mt <= 12
)

p_qc_after <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
) + theme_classic(base_size = 13)

save_pdf(
  p_qc_after,
  file.path(out_dir, "plots", "02_QC_after_filter.pdf"),
  width = 9,
  height = 4
)

write.csv(
  obj@meta.data,
  file.path(out_dir, "tables", "00_metadata_after_qc.csv"),
  row.names = TRUE
)

########################################################
## 5. choose batch/sample column for CCA integration
########################################################
batch_col <- pick_col(
  obj@meta.data,
  c("sample", "Sample", "patient", "Patient", "orig.ident", "site", "tissue", "group")
)

if (is.null(batch_col)) {
  message("No obvious sample column found. Use orig.ident as fallback.")
  obj$orig.ident <- "GSE149614_all"
  batch_col <- "orig.ident"
}

cat("Batch / sample column used for split-integration:", batch_col, "\n")

########################################################
## 6. split object by sample and run CCA integration
########################################################
obj_list <- SplitObject(obj, split.by = batch_col)

## 如果只有一个子对象，说明 metadata 里没有真正的样本信息
## 这种情况下不做 integration，直接走常规流程
if (length(obj_list) == 1) {
  message("Only one group found after split. Skip CCA integration and run standard SCT workflow.")

  obj_main <- obj_list[[1]]
  obj_main <- SCTransform(obj_main, vars.to.regress = "percent.mt", verbose = FALSE)
  obj_main <- RunPCA(obj_main, npcs = 50, verbose = FALSE)
  obj_main <- RunUMAP(obj_main, dims = 1:30)
  obj_main <- FindNeighbors(obj_main, dims = 1:30)
  obj_main <- FindClusters(obj_main, resolution = 0.4)

} else {
  message("Multiple groups found. Running CCA-style Seurat integration.")

  obj_list <- lapply(obj_list, function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    x
  })

  anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    dims = 1:30
  )

  obj_main <- IntegrateData(
    anchorset = anchors,
    dims = 1:30
  )

  DefaultAssay(obj_main) <- "integrated"
  obj_main <- ScaleData(obj_main, verbose = FALSE)
  obj_main <- RunPCA(obj_main, npcs = 50, verbose = FALSE)
  obj_main <- RunUMAP(obj_main, dims = 1:30)
  obj_main <- FindNeighbors(obj_main, dims = 1:30)
  obj_main <- FindClusters(obj_main, resolution = 0.4)
}

saveRDS(obj_main, file.path(out_dir, "rds", "01_atlas_clustered_obj.rds"))

########################################################
## 7. first UMAP by cluster
########################################################
p_cluster <- DimPlot(
  obj_main,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2
) + theme_classic(base_size = 13)

save_pdf(
  p_cluster,
  file.path(out_dir, "plots", "03_UMAP_cluster.pdf"),
  width = 7.2,
  height = 6
)

########################################################
## 8. find markers on RNA assay
########################################################
DefaultAssay(obj_main) <- "RNA"
obj_main <- NormalizeData(obj_main, verbose = FALSE)

markers_all <- FindAllMarkers(
  obj_main,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25
)

write.csv(
  markers_all,
  file.path(out_dir, "tables", "01_all_cluster_markers.csv"),
  row.names = FALSE
)

top10 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(
  top10,
  file.path(out_dir, "tables", "02_top10_markers_each_cluster.csv"),
  row.names = FALSE
)

########################################################
## 9. MANUAL cluster-level annotation
########################################################
## 不再自动算哪一类marker平均值最高
## 直接按当前19个cluster手动指定身份

cluster_map <- c(
  "0"  = "Myeloid",
  "1"  = "T_NK",
  "2"  = "Cholangiocyte_like",
  "3"  = "T_NK",
  "4"  = "Hepatocyte_like",
  "5"  = "Hepatocyte_like",
  "6"  = "Endothelial",
  "7"  = "T_NK",
  "8"  = "Hepatocyte_like",
  "9"  = "Myeloid",
  "10" = "T_NK",
  "11" = "Fibro_Stellate",
  "12" = "B_Plasma",
  "13" = "Myeloid",
  "14" = "B_Plasma",
  "15" = "T_NK",
  "16" = "Myeloid",
  "17" = "Myeloid",
  "18" = "Mast"
)

cluster_label_df <- data.frame(
  seurat_clusters = names(cluster_map),
  major_celltype_cluster = unname(cluster_map),
  stringsAsFactors = FALSE
)

write.csv(
  cluster_label_df,
  file.path(out_dir, "tables", "04_cluster_major_annotation_manual.csv"),
  row.names = FALSE
)

obj_main$major_celltype_cluster <- cluster_label_df$major_celltype_cluster[
  match(as.character(obj_main$seurat_clusters), cluster_label_df$seurat_clusters)
]

obj_main$major_celltype_cluster <- factor(
  obj_main$major_celltype_cluster,
  levels = c(
    "B_Plasma",
    "Cholangiocyte_like",
    "Endothelial",
    "Fibro_Stellate",
    "Hepatocyte_like",
    "Mast",
    "Myeloid",
    "T_NK"
  )
)

saveRDS(obj_main, file.path(out_dir, "rds", "02_atlas_clusterlevel_manual_annotated_obj.rds"))

########################################################
## 10. canonical marker dotplot
########################################################
marker_use <- unique(unlist(marker_list))
marker_use <- marker_use[marker_use %in% rownames(obj_main)]

p_dot <- DotPlot(
  obj_main,
  features = marker_use,
  group.by = "major_celltype_cluster"
) +
  RotatedAxis() +
  theme_classic(base_size = 12)

save_pdf(
  p_dot,
  file.path(out_dir, "plots", "04_marker_dotplot_clusterlevel.pdf"),
  width = 12,
  height = 6
)

########################################################
## 11. major celltype UMAP
########################################################
plot_df <- FetchData(
  obj_main,
  vars = c("UMAP_1", "UMAP_2", "major_celltype_cluster", "DLC1")
)

plot_df$major_celltype_cluster <- factor(
  plot_df$major_celltype_cluster,
  levels = c(
    "B_Plasma",
    "Cholangiocyte_like",
    "Endothelial",
    "Fibro_Stellate",
    "Hepatocyte_like",
    "Mast",
    "Myeloid",
    "T_NK"
  )
)

celltype_cols <- c(
  "B_Plasma"           = "#E7D86C",
  "Cholangiocyte_like" = "#F2C790",
  "Endothelial"        = "#78C5C7",
  "Fibro_Stellate"     = "#9D89C9",
  "Hepatocyte_like"    = "#A8CF8A",
  "Mast"               = "#D95F5F",
  "Myeloid"            = "#F39B36",
  "T_NK"               = "#6C97D2"
)

label_df <- plot_df %>%
  group_by(major_celltype_cluster) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

theme_panel <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold", color = "black", size = 11),
    axis.text  = element_text(color = "black", size = 10),
    axis.line  = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.40, color = "black"),
    legend.title = element_blank(),
    legend.text  = element_text(color = "black", size = 9),
    legend.key   = element_blank(),
    plot.margin = margin(4, 6, 4, 4)
  )

p_celltype <- ggplot(
  plot_df,
  aes(x = UMAP_1, y = UMAP_2, color = major_celltype_cluster)
) +
  geom_point(size = 0.15, alpha = 0.88, stroke = 0) +
  ggrepel::geom_label_repel(
    data = label_df,
    aes(label = major_celltype_cluster, color = major_celltype_cluster),
    fill = scales::alpha("white", 0.82),
    label.size = 0,
    fontface = "bold",
    size = 3.6,
    show.legend = FALSE,
    seed = 123,
    box.padding = 0.30,
    point.padding = 0.12,
    segment.color = "grey60",
    segment.size = 0.23
  ) +
  scale_color_manual(values = celltype_cols) +
  coord_equal() +
  labs(x = "UMAP_1", y = "UMAP_2") +
  guides(
    color = guide_legend(
      override.aes = list(size = 2.8, alpha = 1)
    )
  ) +
  theme_panel

save_pdf(
  p_celltype,
  file.path(out_dir, "plots", "05_UMAP_major_celltype_manual.pdf"),
  width = 7.2,
  height = 5.8
)

p_celltype_nolabel <- ggplot(
  plot_df,
  aes(x = UMAP_1, y = UMAP_2, color = major_celltype_cluster)
) +
  geom_point(size = 0.15, alpha = 0.88, stroke = 0) +
  scale_color_manual(values = celltype_cols) +
  coord_equal() +
  labs(x = "UMAP_1", y = "UMAP_2") +
  guides(
    color = guide_legend(
      override.aes = list(size = 2.8, alpha = 1)
    )
  ) +
  theme_panel

save_pdf(
  p_celltype_nolabel,
  file.path(out_dir, "plots", "06_UMAP_major_celltype_manual_nolabel.pdf"),
  width = 7.2,
  height = 5.8
)

########################################################
## 12. DLC1 feature plot
########################################################
if (!"DLC1" %in% rownames(obj_main)) {
  stop("DLC1 not found in object.")
}

feat_df <- plot_df %>% arrange(DLC1)
dlc1_upper <- quantile(feat_df$DLC1, probs = 0.995, na.rm = TRUE)

p_dlc1 <- ggplot() +
  geom_point(
    data = feat_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "#EBEBEE",
    size = 0.14,
    alpha = 0.42,
    stroke = 0
  ) +
  geom_point(
    data = feat_df %>% filter(DLC1 > 0),
    aes(x = UMAP_1, y = UMAP_2, color = pmin(DLC1, dlc1_upper)),
    size = 0.16,
    alpha = 0.95,
    stroke = 0
  ) +
  scale_color_gradientn(
    colors = c(
      "#E7EDF4",
      "#D7CBE8",
      "#C89DCE",
      "#B55C9D",
      "#7A1F5C"
    ),
    limits = c(0, dlc1_upper),
    oob = scales::squish,
    name = "DLC1"
  ) +
  coord_equal() +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_panel +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text  = element_text(size = 8.8)
  )

save_pdf(
  p_dlc1,
  file.path(out_dir, "plots", "07_DLC1_featureplot_manual_atlas.pdf"),
  width = 6.2,
  height = 5.5
)

cat("Original-style restart atlas workflow finished.\n")
cat("Main outputs:\n")
cat(" - 03_UMAP_cluster.pdf\n")
cat(" - 04_marker_dotplot_clusterlevel.pdf\n")
cat(" - 05_UMAP_major_celltype_clusterlevel.pdf\n")
cat(" - 06_UMAP_major_celltype_clusterlevel_nolabel.pdf\n")
cat(" - 07_DLC1_featureplot_clusterlevel_atlas.pdf\n")




########################################################
## 美化UMAP plotting (beautified version for 05 / 06 / 07)
########################################################

## ---- packages ----
if (!requireNamespace("SCpubr", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("enblacar/SCpubr")
}
library(SCpubr)
library(ggrepel)
library(scales)
library(ggplot2)
library(dplyr)

## ---- plotting data ----
plot_df <- FetchData(
  obj_main,
  vars = c("UMAP_1", "UMAP_2", "major_celltype_cluster", "DLC1")
)

plot_df$major_celltype_cluster <- factor(
  plot_df$major_celltype_cluster,
  levels = c(
    "B_Plasma",
    "Cholangiocyte_like",
    "Endothelial",
    "Fibro_Stellate",
    "Hepatocyte_like",
    "Mast",
    "Myeloid",
    "T_NK"
  )
)

label_df <- plot_df %>%
  group_by(major_celltype_cluster) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

## ---- redesigned palette for 05 / 06 ----
## 目标：对比更明显，同时保留期刊风格的低饱和度
celltype_cols <- c(
  "B_Plasma"           = "#8EC9F3",  # 冰蓝
  "Cholangiocyte_like" = "#E9B7A5",  # 玫瑰米粉
  "Endothelial"        = "#B790D4",  # 淡紫
  "Fibro_Stellate"     = "#5DB7A0",  # 玉石绿
  "Hepatocyte_like"    = "#D9C554",  # 芥末黄
  "Mast"               = "#E35D5B",  # 珊瑚红
  "Myeloid"            = "#F08B2C",  # 橙
  "T_NK"               = "#4D7FB8"   # 海军蓝
)

## ---- helper ----
save_pdf <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

theme_panel <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold", color = "black", size = 11),
    axis.text  = element_text(color = "black", size = 10),
    axis.line  = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    legend.title = element_blank(),
    legend.text  = element_text(color = "black", size = 10),
    legend.key   = element_blank(),
    plot.margin  = margin(4, 6, 4, 4)
  )

########################################################
## 05. label version
########################################################
p05 <- SCpubr::do_DimPlot(
  sample = obj_main,
  group.by = "major_celltype_cluster",
  reduction = "umap",
  colors.use = celltype_cols,
  raster = FALSE,
  shuffle = TRUE,
  pt.size = 0.18,
  label = FALSE,
  legend.position = "right"
) +
  ggrepel::geom_label_repel(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = major_celltype_cluster, color = major_celltype_cluster),
    fill = alpha("white", 0.88),
    label.size = 0,
    fontface = "bold",
    size = 3.8,
    show.legend = FALSE,
    seed = 123,
    box.padding = 0.25,
    point.padding = 0.10,
    segment.color = "grey70",
    segment.size = 0.20
  ) +
  scale_color_manual(values = celltype_cols) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_panel

save_pdf(
  p05,
  file.path(out_dir, "plots", "05_UMAP_major_celltype_SCpubr_final.pdf"),
  width = 7.2,
  height = 5.8
)

########################################################
## 06. no-label version
########################################################
p06 <- SCpubr::do_DimPlot(
  sample = obj_main,
  group.by = "major_celltype_cluster",
  reduction = "umap",
  colors.use = celltype_cols,
  raster = FALSE,
  shuffle = TRUE,
  pt.size = 0.18,
  label = FALSE,
  legend.position = "right"
) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_panel

save_pdf(
  p06,
  file.path(out_dir, "plots", "06_UMAP_major_celltype_SCpubr_final_nolabel.pdf"),
  width = 7.2,
  height = 5.8
)

## 07. DLC1 feature plot
## redesigned by me: mist rose -> deep violet

feat_df <- FetchData(
  obj_main,
  vars = c("UMAP_1", "UMAP_2", "DLC1")
) %>%
  dplyr::arrange(DLC1)

dlc1_upper <- quantile(feat_df$DLC1, probs = 0.995, na.rm = TRUE)

p07 <- ggplot() +
  ## 背景细胞
  geom_point(
    data = feat_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "#EEE8E1",
    size = 0.14,
    alpha = 0.55,
    stroke = 0
  ) +
  ## DLC1 阳性细胞
  geom_point(
    data = feat_df %>% dplyr::filter(DLC1 > 0),
    aes(x = UMAP_1, y = UMAP_2, color = pmin(DLC1, dlc1_upper)),
    size = 0.18,
    alpha = 0.98,
    stroke = 0
  ) +
  scale_color_gradientn(
    colors = c(
      "#E6CCD6",  # pale rose
      "#D38CB7",  # rose-mauve
      "#A55AA5",  # orchid purple
      "#5E2E86",  # deep violet
      "#1F163F"   # almost-black purple
    ),
    limits = c(0, dlc1_upper),
    oob = scales::squish,
    name = "DLC1"
  ) +
  coord_equal() +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold", color = "black", size = 11),
    axis.text  = element_text(color = "black", size = 10),
    axis.line  = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    legend.key   = element_blank(),
    plot.margin  = margin(4, 6, 4, 4)
  )

save_pdf(
  p07,
  file.path(out_dir, "plots", "07_DLC1_featureplot_designed_by_oai.pdf"),
  width = 6.2,
  height = 5.5
)



########################################################
## 2成纤维细胞亚群的细分
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "scales",
  "patchwork"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

## SCpubr
if (!requireNamespace("SCpubr", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("enblacar/SCpubr")
}
suppressPackageStartupMessages(library(SCpubr))

save_pdf <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2")

obj_file <- file.path(fib_dir, "rds", "01_fibro_clean_reclustered_obj.rds")
if (!file.exists(obj_file)) {
  stop("找不到对象文件: ", obj_file)
}

out_dir <- file.path(fib_dir, "fibro_state_final_v3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

########################################################
## 2. read object
########################################################
fib_clean <- readRDS(obj_file)

if (!"SCT_snn_res.0.4" %in% colnames(fib_clean@meta.data)) {
  stop("对象中没有 SCT_snn_res.0.4 列。")
}

fib_clean$fibro_subcluster_clean <- fib_clean$SCT_snn_res.0.4

########################################################
## 3. final fibro_state annotation
########################################################
cluster_map <- c(
  "0" = "Perivascular_myofibroblast",
  "1" = "Activated_matrix_fibroblast",
  "2" = "Transitional_fibroblast",
  "3" = "C7_SFRP4_stromal",
  "4" = "CXCL12_RBP1_stellate",
  "5" = "FABP4_CYGB_stellate",
  "6" = "PI16_IGF1_quiescent_fibro"
)

fib_clean$fibro_state <- cluster_map[as.character(fib_clean$fibro_subcluster_clean)]
fib_clean$fibro_state <- factor(
  fib_clean$fibro_state,
  levels = c(
    "Perivascular_myofibroblast",
    "Activated_matrix_fibroblast",
    "Transitional_fibroblast",
    "C7_SFRP4_stromal",
    "CXCL12_RBP1_stellate",
    "FABP4_CYGB_stellate",
    "PI16_IGF1_quiescent_fibro"
  )
)

Idents(fib_clean) <- "fibro_state"

cluster_anno_df <- data.frame(
  fibro_subcluster_clean = names(cluster_map),
  fibro_state = unname(cluster_map),
  stringsAsFactors = FALSE
)

write.csv(
  cluster_anno_df,
  file.path(out_dir, "tables", "01_fibro_state_annotation_map.csv"),
  row.names = FALSE
)

saveRDS(
  fib_clean,
  file.path(out_dir, "rds", "02_fibro_state_final_obj.rds")
)

########################################################
## 4. colors and themes
########################################################
## 重新设计：对比更明显，但仍保留期刊风格
fib_state_cols4 <- c(
  "Perivascular_myofibroblast"  = "#7A6FE3",  # 紫
  "Activated_matrix_fibroblast" = "#E86B61",  # 珊瑚红
  "Transitional_fibroblast"     = "#7E8F3A",  # 深橄榄棕
  "C7_SFRP4_stromal"            = "#E0A11B",  # 亮橙金
  "CXCL12_RBP1_stellate"        = "#35B44A",  # 更纯亮绿
  "FABP4_CYGB_stellate"         = "#16B7A6",  # 青绿
  "PI16_IGF1_quiescent_fibro"   = "#47A6E8"   # 天蓝
)

## 缩短标签，避免图面拥挤
fib_state_short_map <- c(
  "Perivascular_myofibroblast"  = "Perivascular\nmyofib",
  "Activated_matrix_fibroblast" = "Activated\nmatrix",
  "Transitional_fibroblast"     = "Transitional",
  "C7_SFRP4_stromal"            = "C7/SFRP4",
  "CXCL12_RBP1_stellate"        = "CXCL12/RBP1",
  "FABP4_CYGB_stellate"         = "FABP4/CYGB",
  "PI16_IGF1_quiescent_fibro"   = "PI16/IGF1"
)

theme_umap_pub <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold", size = 12, color = "black"),
    axis.text  = element_text(size = 10, color = "black"),
    axis.line  = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 10, color = "black"),
    legend.key   = element_blank(),
    plot.margin  = margin(5, 8, 5, 5)
  )

theme_dot_pub <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

########################################################
## 5. plotting data
########################################################
plot_df <- FetchData(
  fib_clean,
  vars = c("UMAP_1", "UMAP_2", "fibro_state", "DLC1")
)

plot_df$fibro_state <- factor(
  plot_df$fibro_state,
  levels = levels(fib_clean$fibro_state)
)

label_df <- plot_df %>%
  dplyr::group_by(fibro_state) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

label_df$fibro_state_short <- fib_state_short_map[as.character(label_df$fibro_state)]

########################################################
## 6. 01-02 fibro_state UMAP (SCpubr)
########################################################
## 01 label version
p_state_label3 <- SCpubr::do_DimPlot(
  sample = fib_clean,
  group.by = "fibro_state",
  reduction = "umap",
  colors.use = fib_state_cols4,
  shuffle = TRUE,
  raster = FALSE,
  pt.size = 0.80,
  label = FALSE,
  legend.position = "right"
) +
  ggrepel::geom_label_repel(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = fibro_state_short),
    fill = alpha("white", 0.86),
    color = "black",
    fontface = "bold",
    label.size = 0,
    size = 3.8,
    show.legend = FALSE,
    seed = 123,
    box.padding = 0.30,
    point.padding = 0.15,
    segment.color = "grey70",
    segment.size = 0.22
  ) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_umap_pub

save_pdf(
  p_state_label3,
  file.path(out_dir, "plots", "01_fibro_state_umap_label_v3.pdf"),
  width = 8.0,
  height = 6.5
)

## 02 no-label version
p_state_nolabel3 <- SCpubr::do_DimPlot(
  sample = fib_clean,
  group.by = "fibro_state",
  reduction = "umap",
  colors.use = fib_state_cols4,
  shuffle = TRUE,
  raster = FALSE,
  pt.size = 0.80,
  label = FALSE,
  legend.position = "right"
) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_umap_pub

save_pdf(
  p_state_nolabel3,
  file.path(out_dir, "plots", "02_fibro_state_umap_nolabel_v3.pdf"),
  width = 8.0,
  height = 6.5
)

########################################################
## 7. 03 classic marker dotplot (classic markers only)
########################################################
classic_markers <- c(
  ## Perivascular_myofibroblast
  "ACTA2", "TAGLN", "RGS5",
  ## Activated_matrix_fibroblast
  "SPP1", "CTHRC1", "POSTN",
  ## Transitional_fibroblast
  "TPPP3", "ADAMTS4", "CRISPLD2",
  ## C7_SFRP4_stromal
  "C7", "SFRP4", "CCL19",
  ## CXCL12_RBP1_stellate
  "CXCL12", "RBP1", "IGFBP3",
  ## FABP4_CYGB_stellate
  "FABP4", "CYGB", "REN",
  ## PI16_IGF1_quiescent_fibro
  "PI16", "IGF1", "MFAP5",
  ## gene of interest
  "DLC1"
)

classic_markers <- classic_markers[classic_markers %in% rownames(fib_clean)]

p_dot_classic <- DotPlot(
  fib_clean,
  features = classic_markers,
  group.by = "fibro_state"
) +
  theme_dot_pub

save_pdf(
  p_dot_classic,
  file.path(out_dir, "plots", "03_fibro_state_marker_dotplot_classic_v2.pdf"),
  width = 11.5,
  height = 5.8
)


########################################################
## 8. 04 DLC1 featureplot (old-style method, like previous version)
########################################################
dlc1_vals <- FetchData(fib_clean, vars = "DLC1")[, 1]

## 用更接近旧版的截断方式，不再只叠加阳性点
min_cut <- quantile(dlc1_vals, probs = 0.01, na.rm = TRUE)
max_cut <- quantile(dlc1_vals, probs = 0.995, na.rm = TRUE)

p_dlc1_feat3 <- FeaturePlot(
  object = fib_clean,
  features = "DLC1",
  reduction = "umap",
  order = TRUE,
  pt.size = 0.78,
  raster = FALSE,
  min.cutoff = min_cut,
  max.cutoff = max_cut,
  cols = c(
    "#D3D3D3",  # low / background grey
    "#C7B9F5",  # light lavender
    "#8F72EA",  # medium purple
    "#5642E6",  # strong purple-blue
    "#2E22D8"   # high
  )
) +
  labs(title = "DLC1", x = "UMAP_1", y = "UMAP_2") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0),
    axis.title = element_text(face = "bold", size = 11, color = "black"),
    axis.text  = element_text(size = 10, color = "black"),
    axis.line  = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    legend.title = element_blank(),
    legend.text  = element_text(size = 9),
    legend.key.height = unit(26, "pt"),
    legend.key.width  = unit(8, "pt"),
    plot.margin = margin(4, 6, 4, 4)
  )

save_pdf(
  p_dlc1_feat3,
  file.path(out_dir, "plots", "04_fibro_state_DLC1_featureplot_v3.pdf"),
  width = 6.8,
  height = 5.8
)

########################################################
## 9. 05 DLC1 violin
########################################################
dlc1_df <- FetchData(
  fib_clean,
  vars = c("DLC1", "fibro_state")
)

dlc1_df$fibro_state <- factor(
  dlc1_df$fibro_state,
  levels = levels(fib_clean$fibro_state)
)

p_dlc1_vln <- ggplot(
  dlc1_df,
  aes(x = fibro_state, y = DLC1, fill = fibro_state)
) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.92, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.70, color = "black") +
  scale_fill_manual(values = fib_state_cols3) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "DLC1 expression") +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

save_pdf(
  p_dlc1_vln,
  file.path(out_dir, "plots", "05_fibro_state_DLC1_violin_v2.pdf"),
  width = 8.8,
  height = 5.2
)

########################################################
## 10. 06 DLC1 dotplot
########################################################
p_dlc1_dot <- DotPlot(
  fib_clean,
  features = "DLC1",
  group.by = "fibro_state"
) +
  theme_classic(base_size = 12)

save_pdf(
  p_dlc1_dot,
  file.path(out_dir, "plots", "06_fibro_state_DLC1_dotplot_v2.pdf"),
  width = 5.2,
  height = 4.5
)

########################################################
## 11. DLC1 summary table
########################################################
dlc1_summary <- dlc1_df %>%
  dplyr::group_by(fibro_state) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_expr = mean(DLC1, na.rm = TRUE),
    median_expr = median(DLC1, na.rm = TRUE),
    pct_positive = mean(DLC1 > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

write.csv(
  dlc1_summary,
  file.path(out_dir, "tables", "02_DLC1_summary_by_fibro_state.csv"),
  row.names = FALSE
)

########################################################
## 12. README
########################################################
writeLines(
  c(
    "Final fibro_state plotting v3 finished.",
    "",
    "Generated files:",
    "01_fibro_state_umap_label_v3.pdf",
    "02_fibro_state_umap_nolabel_v3.pdf",
    "03_fibro_state_marker_dotplot_classic_v2.pdf",
    "04_fibro_state_DLC1_featureplot_v2.pdf",
    "05_fibro_state_DLC1_violin_v2.pdf",
    "06_fibro_state_DLC1_dotplot_v2.pdf",
    "02_DLC1_summary_by_fibro_state.csv"
  ),
  con = file.path(out_dir, "README_plotting_v3.txt")
)

cat("Plotting v3 finished.\n")







########################################################
## 3成纤维细胞亚群的差异和富集分析
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

save_pdf <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

obj_file <- file.path(fib_dir, "rds", "02_fibro_state_final_obj.rds")
if (!file.exists(obj_file)) {
  stop("找不到对象文件: ", obj_file)
}

out_dir <- file.path(fib_dir, "DLC1_state_defined_high_low")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)

########################################################
## 2. read object
########################################################
fib_clean <- readRDS(obj_file)

if (!"fibro_state" %in% colnames(fib_clean@meta.data)) {
  stop("对象中没有 fibro_state 列。")
}
if (!"DLC1" %in% rownames(fib_clean)) {
  stop("对象中找不到 DLC1。")
}

DefaultAssay(fib_clean) <- "RNA"
fib_clean <- NormalizeData(fib_clean, verbose = FALSE)

########################################################
## 3. define DLC1-high vs DLC1-low by fibro_state
########################################################
## 采用你旧版的思路：按“高表达状态”和“低表达状态”分组，而不是按细胞分位数硬切

high_states <- c(
  "FABP4_CYGB_stellate",
  "Perivascular_myofibroblast",
  "CXCL12_RBP1_stellate"
)

low_states <- c(
  "PI16_IGF1_quiescent_fibro",
  "Transitional_fibroblast",
  "Activated_matrix_fibroblast"
)

fib_clean$DLC1_group <- NA_character_
fib_clean$DLC1_group[fib_clean$fibro_state %in% high_states] <- "High"
fib_clean$DLC1_group[fib_clean$fibro_state %in% low_states]  <- "Low"

keep_cells_hl <- colnames(fib_clean)[!is.na(fib_clean$DLC1_group)]
fib_hl <- subset(fib_clean, cells = keep_cells_hl)

fib_hl$DLC1_group <- factor(fib_hl$DLC1_group, levels = c("Low", "High"))
Idents(fib_hl) <- fib_hl$DLC1_group

print(table(fib_clean$DLC1_group, useNA = "ifany"))
print(dim(fib_hl))

write.csv(
  data.frame(
    cell = colnames(fib_clean),
    fibro_state = fib_clean$fibro_state,
    DLC1_group = fib_clean$DLC1_group
  ),
  file.path(out_dir, "tables", "00_DLC1_state_defined_group_assignment.csv"),
  row.names = FALSE
)

saveRDS(
  fib_hl,
  file.path(out_dir, "rds", "01_fibro_state_defined_DLC1_highlow_obj.rds")
)

########################################################
## 4. High/Low UMAP
########################################################
p_group_umap <- DimPlot(
  fib_hl,
  reduction = "umap",
  group.by = "DLC1_group",
  cols = c("Low" = "#4C78A8", "High" = "#D95F5F"),
  pt.size = 0.65
) + theme_classic(base_size = 13)

save_pdf(
  p_group_umap,
  file.path(out_dir, "plots", "01_DLC1_state_defined_highlow_umap.pdf"),
  width = 6.8,
  height = 5.5
)

p_group_violin <- VlnPlot(
  fib_hl,
  features = "DLC1",
  group.by = "DLC1_group",
  pt.size = 0
) + theme_classic(base_size = 12)

save_pdf(
  p_group_violin,
  file.path(out_dir, "plots", "02_DLC1_state_defined_highlow_violin.pdf"),
  width = 4.8,
  height = 4.5
)

########################################################
## 5. DEG
########################################################
DefaultAssay(fib_hl) <- "RNA"
fib_hl <- NormalizeData(fib_hl, verbose = FALSE)
Idents(fib_hl) <- fib_hl$DLC1_group

deg_fib <- FindMarkers(
  fib_hl,
  ident.1 = "Low",
  ident.2 = "High",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  logfc.threshold = 0.1,
  min.pct = 0.1,
  verbose = FALSE
)

deg_fib$gene <- rownames(deg_fib)
deg_fib <- deg_fib[order(deg_fib$p_val_adj, -abs(deg_fib$avg_log2FC)), ]

write.csv(
  deg_fib,
  file.path(out_dir, "tables", "01_DLC1_state_defined_Low_vs_High_DEG.csv"),
  row.names = FALSE
)

########################################################
## 6. volcano
########################################################
vol_df <- deg_fib
vol_df$group <- "NS"
vol_df$group[vol_df$p_val_adj < 0.05 & vol_df$avg_log2FC >  0.25] <- "Up_in_Low"
vol_df$group[vol_df$p_val_adj < 0.05 & vol_df$avg_log2FC < -0.25] <- "Up_in_High"
vol_df$neglog10 <- -log10(vol_df$p_val_adj + 1e-300)

lab_up <- vol_df %>%
  dplyr::filter(group == "Up_in_Low") %>%
  dplyr::arrange(dplyr::desc(avg_log2FC), dplyr::desc(neglog10)) %>%
  dplyr::slice_head(n = 8)

lab_down <- vol_df %>%
  dplyr::filter(group == "Up_in_High") %>%
  dplyr::arrange(avg_log2FC, dplyr::desc(neglog10)) %>%
  dplyr::slice_head(n = 8)

lab_df <- rbind(lab_up, lab_down)

p_vol <- ggplot(vol_df, aes(x = avg_log2FC, y = neglog10)) +
  geom_point(aes(color = group), size = 1.3, alpha = 0.8) +
  scale_color_manual(values = c(
    "Up_in_Low" = "#D55E5E",
    "Up_in_High" = "#4C78A8",
    "NS" = "grey75"
  )) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = lab_df,
    aes(label = gene),
    size = 4,
    max.overlaps = Inf
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("Adjusted P")),
    color = NULL,
    title = "Fibro states: DLC1-Low vs DLC1-High"
  )

save_pdf(
  p_vol,
  file.path(out_dir, "plots", "03_DLC1_state_defined_volcano.pdf"),
  width = 8,
  height = 6.5
)

########################################################
## 7. enrichment input
########################################################
deg_sig <- deg_fib %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

low_up_genes <- deg_sig %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::pull(gene) %>%
  unique()

high_up_genes <- deg_sig %>%
  dplyr::filter(avg_log2FC < 0) %>%
  dplyr::pull(gene) %>%
  unique()

gene_df_low <- tryCatch({
  bitr(
    low_up_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
}, error = function(e) NULL)

gene_df_high <- tryCatch({
  bitr(
    high_up_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
}, error = function(e) NULL)

########################################################
## 8. GO BP enrichment
########################################################
if (!is.null(gene_df_low) && nrow(gene_df_low) > 0) {
  ego_low <- enrichGO(
    gene = unique(gene_df_low$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  write.csv(
    as.data.frame(ego_low),
    file.path(out_dir, "tables", "02_GO_BP_up_in_Low.csv"),
    row.names = FALSE
  )

  save_pdf(
    dotplot(ego_low, showCategory = 15) + theme_classic(base_size = 12),
    file.path(out_dir, "plots", "04_GO_BP_up_in_Low.pdf"),
    width = 8,
    height = 6
  )
}

if (!is.null(gene_df_high) && nrow(gene_df_high) > 0) {
  ego_high <- enrichGO(
    gene = unique(gene_df_high$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  write.csv(
    as.data.frame(ego_high),
    file.path(out_dir, "tables", "03_GO_BP_up_in_High.csv"),
    row.names = FALSE
  )

  save_pdf(
    dotplot(ego_high, showCategory = 15) + theme_classic(base_size = 12),
    file.path(out_dir, "plots", "05_GO_BP_up_in_High.pdf"),
    width = 8,
    height = 6
  )
}

########################################################
## 9. KEGG enrichment
########################################################
if (!is.null(gene_df_low) && nrow(gene_df_low) > 0) {
  ekegg_low <- tryCatch({
    enrichKEGG(
      gene = unique(gene_df_low$ENTREZID),
      organism = "hsa",
      pvalueCutoff = 0.05
    )
  }, error = function(e) NULL)

  if (!is.null(ekegg_low)) {
    write.csv(
      as.data.frame(ekegg_low),
      file.path(out_dir, "tables", "04_KEGG_up_in_Low.csv"),
      row.names = FALSE
    )

    save_pdf(
      dotplot(ekegg_low, showCategory = 15) + theme_classic(base_size = 12),
      file.path(out_dir, "plots", "06_KEGG_up_in_Low.pdf"),
      width = 8,
      height = 6
    )
  }
}

if (!is.null(gene_df_high) && nrow(gene_df_high) > 0) {
  ekegg_high <- tryCatch({
    enrichKEGG(
      gene = unique(gene_df_high$ENTREZID),
      organism = "hsa",
      pvalueCutoff = 0.05
    )
  }, error = function(e) NULL)

  if (!is.null(ekegg_high)) {
    write.csv(
      as.data.frame(ekegg_high),
      file.path(out_dir, "tables", "05_KEGG_up_in_High.csv"),
      row.names = FALSE
    )

    save_pdf(
      dotplot(ekegg_high, showCategory = 15) + theme_classic(base_size = 12),
      file.path(out_dir, "plots", "07_KEGG_up_in_High.pdf"),
      width = 8,
      height = 6
    )
  }
}

########################################################
## 10. GSEA GO
########################################################
gene_list <- deg_fib$avg_log2FC
names(gene_list) <- deg_fib$gene
gene_list <- sort(gene_list, decreasing = TRUE)

gene_df2 <- tryCatch({
  bitr(
    names(gene_list),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
}, error = function(e) NULL)

if (!is.null(gene_df2) && nrow(gene_df2) > 0) {
  gene_list2 <- gene_list[gene_df2$SYMBOL]
  names(gene_list2) <- gene_df2$ENTREZID
  gene_list2 <- sort(gene_list2, decreasing = TRUE)

  ggo <- gseGO(
    geneList = gene_list2,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    verbose = FALSE,
    eps = 0
  )

  write.csv(
    as.data.frame(ggo),
    file.path(out_dir, "tables", "06_GSEA_GO.csv"),
    row.names = FALSE
  )

  save_pdf(
    dotplot(ggo, showCategory = 15) + theme_classic(base_size = 12),
    file.path(out_dir, "plots", "08_GSEA_GO.pdf"),
    width = 8,
    height = 6
  )
}

########################################################
## 11. summary
########################################################
writeLines(
  c(
    "State-defined DLC1 high/low DEG + enrichment finished.",
    "",
    "High states:",
    paste(high_states, collapse = ", "),
    "",
    "Low states:",
    paste(low_states, collapse = ", "),
    "",
    "Main outputs:",
    "plots/03_DLC1_state_defined_volcano.pdf",
    "plots/04_GO_BP_up_in_Low.pdf",
    "plots/05_GO_BP_up_in_High.pdf",
    "plots/06_KEGG_up_in_Low.pdf",
    "plots/07_KEGG_up_in_High.pdf",
    "plots/08_GSEA_GO.pdf"
  ),
  con = file.path(out_dir, "README_DLC1_state_defined_analysis.txt")
)

cat("State-defined DLC1 high/low DEG + enrichment finished.\n")






########################################################
## 12. 完整版本的富集分析
########################################################
## 这部分是你要的“完整版”，不再分 Up_in_Low / Up_in_High

deg_sig_all <- deg_fib %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

sig_genes_all <- unique(deg_sig_all$gene)

gene_df_all <- tryCatch({
  bitr(
    sig_genes_all,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
}, error = function(e) NULL)

########################################################
## 12A. complete GO BP
########################################################
if (!is.null(gene_df_all) && nrow(gene_df_all) > 0) {

  ego_all <- enrichGO(
    gene = unique(gene_df_all$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  if (!is.null(ego_all) && nrow(as.data.frame(ego_all)) > 0) {
    go_all_df <- as.data.frame(ego_all)

    write.csv(
      go_all_df,
      file.path(out_dir, "tables", "09_GO_BP_complete_all_sig_DEGs.csv"),
      row.names = FALSE
    )

    p_go_all <- dotplot(ego_all, showCategory = 20) +
      ggtitle("GO BP: all significant DEGs") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(size = 12, face = "bold"))

    save_pdf(
      p_go_all,
      file.path(out_dir, "plots", "09_GO_BP_complete_all_sig_DEGs.pdf"),
      width = 8.5,
      height = 6.5
    )
  }
}

########################################################
## 12B. complete KEGG
########################################################
if (!is.null(gene_df_all) && nrow(gene_df_all) > 0) {

  ekegg_all <- tryCatch({
    enrichKEGG(
      gene = unique(gene_df_all$ENTREZID),
      organism = "hsa",
      pvalueCutoff = 0.05
    )
  }, error = function(e) NULL)

  if (!is.null(ekegg_all) && nrow(as.data.frame(ekegg_all)) > 0) {
    kegg_all_df <- as.data.frame(ekegg_all)

    write.csv(
      kegg_all_df,
      file.path(out_dir, "tables", "10_KEGG_complete_all_sig_DEGs.csv"),
      row.names = FALSE
    )

    p_kegg_all <- dotplot(ekegg_all, showCategory = 20) +
      ggtitle("KEGG: all significant DEGs") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(size = 12, face = "bold"))

    save_pdf(
      p_kegg_all,
      file.path(out_dir, "plots", "10_KEGG_complete_all_sig_DEGs.pdf"),
      width = 8.5,
      height = 6.5
    )
  }
}

########################################################
## 12C. optional merged directional tables
########################################################
## 如果你还想把分开的结果合成一个总表，带 Direction 列，也一起输出

if (exists("ego_low") && !is.null(ego_low) && nrow(as.data.frame(ego_low)) > 0) {
  go_low_df <- as.data.frame(ego_low)
  go_low_df$Direction <- "Up_in_Low"
} else {
  go_low_df <- NULL
}

if (exists("ego_high") && !is.null(ego_high) && nrow(as.data.frame(ego_high)) > 0) {
  go_high_df <- as.data.frame(ego_high)
  go_high_df$Direction <- "Up_in_High"
} else {
  go_high_df <- NULL
}

go_merged <- dplyr::bind_rows(go_low_df, go_high_df)
if (!is.null(go_merged) && nrow(go_merged) > 0) {
  write.csv(
    go_merged,
    file.path(out_dir, "tables", "11_GO_BP_direction_merged.csv"),
    row.names = FALSE
  )
}

if (exists("ekegg_low") && !is.null(ekegg_low) && nrow(as.data.frame(ekegg_low)) > 0) {
  kegg_low_df <- as.data.frame(ekegg_low)
  kegg_low_df$Direction <- "Up_in_Low"
} else {
  kegg_low_df <- NULL
}

if (exists("ekegg_high") && !is.null(ekegg_high) && nrow(as.data.frame(ekegg_high)) > 0) {
  kegg_high_df <- as.data.frame(ekegg_high)
  kegg_high_df$Direction <- "Up_in_High"
} else {
  kegg_high_df <- NULL
}

kegg_merged <- dplyr::bind_rows(kegg_low_df, kegg_high_df)
if (!is.null(kegg_merged) && nrow(kegg_merged) > 0) {
  write.csv(
    kegg_merged,
    file.path(out_dir, "tables", "12_KEGG_direction_merged.csv"),
    row.names = FALSE
  )
}







########################################################
## 成纤维细胞亚群的可高变基因分析
########################################################


rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "Seurat",
  "ggplot2",
  "dplyr"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

save_pdf <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

obj_file <- file.path(fib_dir, "rds", "02_fibro_state_final_obj.rds")
if (!file.exists(obj_file)) {
  stop("找不到对象文件: ", obj_file)
}

out_dir <- file.path(fib_dir, "HVG_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

########################################################
## 2. read object
########################################################
fib_clean <- readRDS(obj_file)

if (!"fibro_state" %in% colnames(fib_clean@meta.data)) {
  stop("对象中没有 fibro_state 列。")
}



########################################################
## 3. choose assay and get HVGs
########################################################
## 不用 Assays(fib_clean) 直接做 %in%，
## 改成读取 assay 名称向量，兼容性更好

assay_names <- names(fib_clean@assays)
print(assay_names)

if ("SCT" %in% assay_names) {
  DefaultAssay(fib_clean) <- "SCT"

  ## 如果还没有 VariableFeatures，就补一次
  if (length(VariableFeatures(fib_clean)) == 0) {
    fib_clean <- FindVariableFeatures(
      fib_clean,
      selection.method = "vst",
      nfeatures = 2000,
      verbose = FALSE
    )
  }

  assay_used <- "SCT"

} else if ("RNA" %in% assay_names) {
  DefaultAssay(fib_clean) <- "RNA"
  fib_clean <- NormalizeData(fib_clean, verbose = FALSE)
  fib_clean <- FindVariableFeatures(
    fib_clean,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  assay_used <- "RNA"

} else {
  stop("对象中既没有 SCT assay，也没有 RNA assay。可用 assay：", paste(assay_names, collapse = ", "))
}

cat("Assay used for HVG analysis:", assay_used, "\n")

hvg_all <- VariableFeatures(fib_clean)

if (length(hvg_all) < 30) {
  stop("高可变基因数量少于30，无法继续。")
}

top200_hvg <- head(hvg_all, min(200, length(hvg_all)))
top30_hvg  <- head(hvg_all, min(30,  length(hvg_all)))

write.csv(
  data.frame(gene = top200_hvg),
  file.path(out_dir, "tables", "01_Fibro_state_top200_HVG.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = top30_hvg),
  file.path(out_dir, "tables", "02_Fibro_state_top30_HVG.csv"),
  row.names = FALSE
)

########################################################
## 4. HVG scatter plot
########################################################
p_hvg <- VariableFeaturePlot(fib_clean)

top15_label <- head(hvg_all, min(15, length(hvg_all)))
p_hvg_lab <- LabelPoints(
  plot = p_hvg,
  points = top15_label,
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("HVG scatter plot (", assay_used, " assay)")
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

save_pdf(
  p_hvg_lab,
  file.path(out_dir, "plots", "01_Fibro_state_HVG_scatter.pdf"),
  width = 7,
  height = 5.5
)

########################################################
## 5. HVG heatmap by fibro_state
########################################################
## DoHeatmap 需要 scale.data 中有这些基因
## SCT assay通常已经有；如果没有就补 ScaleData
heatmap_features <- top30_hvg

if (DefaultAssay(fib_clean) == "SCT") {
  scaled_genes <- rownames(GetAssayData(fib_clean, assay = "SCT", slot = "scale.data"))
  missing_genes <- setdiff(heatmap_features, scaled_genes)

  if (length(missing_genes) > 0) {
    fib_clean <- ScaleData(
      fib_clean,
      features = unique(c(scaled_genes, heatmap_features)),
      verbose = FALSE
    )
  }
} else {
  fib_clean <- ScaleData(fib_clean, features = heatmap_features, verbose = FALSE)
}

p_hvg_heatmap <- DoHeatmap(
  object = fib_clean,
  features = rev(heatmap_features),
  group.by = "fibro_state",
  size = 3
) +
  NoLegend()

save_pdf(
  p_hvg_heatmap,
  file.path(out_dir, "plots", "02_Fibro_state_HVG_heatmap_by_state.pdf"),
  width = 10,
  height = 7
)

########################################################
## 6. Optional: average expression table of top30 HVGs by fibro_state
########################################################
avg_hvg <- AverageExpression(
  fib_clean,
  assays = DefaultAssay(fib_clean),
  features = top30_hvg,
  group.by = "fibro_state",
  slot = "data"
)

avg_hvg_mat <- avg_hvg[[DefaultAssay(fib_clean)]]
avg_hvg_df <- data.frame(
  gene = rownames(avg_hvg_mat),
  avg_hvg_mat,
  check.names = FALSE
)

write.csv(
  avg_hvg_df,
  file.path(out_dir, "tables", "03_Fibro_state_top30_HVG_average_expression_by_state.csv"),
  row.names = FALSE
)

########################################################
## 7. save object
########################################################
saveRDS(
  fib_clean,
  file.path(out_dir, "Fibro_state_HVG_analysis_obj.rds")
)

cat("Fibro_state HVG analysis finished.\n")
cat("Main outputs:\n")
cat(" - plots/01_Fibro_state_HVG_scatter.pdf\n")
cat(" - plots/02_Fibro_state_HVG_heatmap_by_state.pdf\n")
cat(" - tables/01_Fibro_state_top200_HVG.csv\n")
cat(" - tables/02_Fibro_state_top30_HVG.csv\n")
cat(" - tables/03_Fibro_state_top30_HVG_average_expression_by_state.csv\n")







########################################################
## 3.成纤维细胞亚群的组织分类
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "scales"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

## SCpubr
if (!requireNamespace("SCpubr", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("enblacar/SCpubr")
}
suppressPackageStartupMessages(library(SCpubr))

save_pdf <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

obj_file <- file.path(fib_dir, "rds", "02_fibro_state_final_obj.rds")
if (!file.exists(obj_file)) {
  stop("找不到对象文件: ", obj_file)
}

out_dir <- file.path(fib_dir, "fibro_state_source_analysis_final_color_v2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

########################################################
## 2. read object
########################################################
fib_obj <- readRDS(obj_file)

if (!"fibro_state" %in% colnames(fib_obj@meta.data)) {
  stop("对象中没有 fibro_state 列。")
}

## 优先 source；没有就 fallback 到 site
if ("source" %in% colnames(fib_obj@meta.data)) {
  source_var <- "source"
} else if ("site" %in% colnames(fib_obj@meta.data)) {
  source_var <- "site"
} else {
  stop("对象中没有 source 或 site 列。")
}

cat("Using source variable:", source_var, "\n")

########################################################
## 3. format source
########################################################
source_vec <- as.character(fib_obj@meta.data[[source_var]])
fib_obj@meta.data[[source_var]] <- source_vec

preferred_source_order <- c("PVTT", "Tumor", "Normal", "Lymph")
source_unique <- unique(source_vec)

source_levels <- preferred_source_order[preferred_source_order %in% source_unique]
source_levels <- c(source_levels, base::setdiff(source_unique, source_levels))

fib_obj@meta.data[[source_var]] <- factor(
  fib_obj@meta.data[[source_var]],
  levels = source_levels
)

########################################################
## 4. reference-style vivid palette
########################################################
## 按你发的参考图风格：亮红 / 亮绿 / 亮青 / 紫
source_cols <- c(
  "PVTT"   = "#14B8B5",  # 清亮青绿
  "Tumor"  = "#D14A3C",  # 珊瑚红
  "Normal" = "#A8D86F",  # 新鲜浅绿
  "Lymph"  = "#9A7FD1"   # 柔和紫
)

extra_sources <- base::setdiff(levels(fib_obj@meta.data[[source_var]]), names(source_cols))
if (length(extra_sources) > 0) {
  extra_cols <- scales::hue_pal()(length(extra_sources))
  names(extra_cols) <- extra_sources
  source_cols <- c(source_cols, extra_cols)
}

theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 11, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10, color = "black"),
    legend.key = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(5, 8, 5, 5)
  )

########################################################
## 5. 71_Fibro_state_source_composition.pdf
########################################################
meta_df <- fib_obj@meta.data %>%
  dplyr::mutate(
    fibro_state = .data[["fibro_state"]],
    Source = .data[[source_var]]
  ) %>%
  dplyr::filter(!is.na(fibro_state), !is.na(Source))

comp_df <- meta_df %>%
  dplyr::count(fibro_state, Source, name = "Count") %>%
  dplyr::group_by(fibro_state) %>%
  dplyr::mutate(Composition = Count / sum(Count)) %>%
  dplyr::ungroup()

write.csv(
  comp_df,
  file.path(out_dir, "tables", "71_Fibro_state_source_composition.csv"),
  row.names = FALSE
)

p_comp <- ggplot(comp_df, aes(x = fibro_state, y = Composition, fill = Source)) +
  geom_bar(stat = "identity", width = 0.74, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = source_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Composition") +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  )

save_pdf(
  p_comp,
  file.path(out_dir, "plots", "71_Fibro_state_source_composition.pdf"),
  width = 8.2,
  height = 5.4
)

########################################################
## 6. 72_Fibro_state_UMAP_by_source.pdf
########################################################
p_source <- SCpubr::do_DimPlot(
  sample = fib_obj,
  group.by = source_var,
  reduction = "umap",
  colors.use = source_cols,
  shuffle = TRUE,
  raster = FALSE,
  pt.size = 0.90,
  label = FALSE,
  legend.position = "right"
) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_pub +
  ggtitle("Source")

save_pdf(
  p_source,
  file.path(out_dir, "plots", "72_Fibro_state_UMAP_by_source.pdf"),
  width = 6.8,
  height = 5.7
)

########################################################
## 7. README
########################################################
writeLines(
  c(
    "Final fibro_state source analysis finished.",
    paste0("Source variable used: ", source_var),
    "",
    "Only two output plots were generated:",
    "plots/71_Fibro_state_source_composition.pdf",
    "plots/72_Fibro_state_UMAP_by_source.pdf"
  ),
  con = file.path(out_dir, "README_fibro_state_source_analysis_final_color_v2.txt")
)

cat("Done: only 71 and 72_Fibro_state_UMAP_by_source.pdf were generated.\n")






########################################################
## 4成纤维细胞亚群的Monocle2分析
########################################################
########################################################
## Monocle2轨迹订终点
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "Matrix",
  "monocle"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p == "monocle") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("monocle", update = FALSE, ask = FALSE)
    } else {
      install.packages(p)
    }
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

save_pdf <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

obj_file <- file.path(fib_dir, "rds", "02_fibro_state_final_obj.rds")

if (!file.exists(obj_file)) {
  stop("找不到对象文件: ", obj_file)
}

out_dir <- file.path(fib_dir, "Monocle2_4states_CXCL12_root")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)

########################################################
## 2. read fibro object
########################################################
fib_obj <- readRDS(obj_file)

if (!"fibro_state" %in% colnames(fib_obj@meta.data)) {
  stop("对象中没有 fibro_state 列。")
}

if (!"RNA" %in% names(fib_obj@assays)) {
  stop("对象中没有 RNA assay。")
}

if (!"DLC1" %in% rownames(fib_obj)) {
  warning("对象中没有 DLC1，后面不会画 DLC1 pseudotime 图。")
}

########################################################
## 3. select 4 trajectory states
########################################################
traj_states <- c(
  "CXCL12_RBP1_stellate",
  "C7_SFRP4_stromal",
  "Perivascular_myofibroblast",
  "Activated_matrix_fibroblast"
)

traj_order <- c(
  "CXCL12_RBP1_stellate",
  "C7_SFRP4_stromal",
  "Perivascular_myofibroblast",
  "Activated_matrix_fibroblast"
)

keep_cells <- rownames(fib_obj@meta.data)[fib_obj@meta.data$fibro_state %in% traj_states]
fib4 <- subset(fib_obj, cells = keep_cells)

fib4$fibro_state <- factor(fib4$fibro_state, levels = traj_order)
Idents(fib4) <- "fibro_state"

DefaultAssay(fib4) <- "RNA"
fib4 <- NormalizeData(fib4, verbose = FALSE)

write.csv(
  as.data.frame(table(fib4$fibro_state)),
  file.path(out_dir, "tables", "01_selected_4states_cell_counts.csv"),
  row.names = FALSE
)

saveRDS(
  fib4,
  file.path(out_dir, "rds", "01_fibro_4states_for_monocle2.rds")
)

########################################################
## 4. extract count matrix for Monocle2
########################################################
get_counts_safe <- function(obj, assay = "RNA") {
  mat <- tryCatch(
    {
      Seurat::GetAssayData(obj, assay = assay, layer = "counts")
    },
    error = function(e) {
      Seurat::GetAssayData(obj, assay = assay, slot = "counts")
    }
  )
  return(mat)
}

count_mat <- get_counts_safe(fib4, assay = "RNA")

## 保证是 sparse matrix
if (!inherits(count_mat, "sparseMatrix")) {
  count_mat <- as(as.matrix(count_mat), "sparseMatrix")
}

## 细胞注释
pd <- fib4@meta.data
pd <- pd[colnames(count_mat), , drop = FALSE]

## 加入 source/site 信息，如果有
if ("source" %in% colnames(pd)) {
  source_var <- "source"
} else if ("site" %in% colnames(pd)) {
  source_var <- "site"
} else {
  source_var <- NA
}

## 基因注释
fd <- data.frame(
  gene_short_name = rownames(count_mat),
  row.names = rownames(count_mat)
)

pd_obj <- new("AnnotatedDataFrame", data = pd)
fd_obj <- new("AnnotatedDataFrame", data = fd)

cds <- newCellDataSet(
  count_mat,
  phenoData = pd_obj,
  featureData = fd_obj,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

########################################################
## 5. Monocle2 preprocessing
########################################################
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- rownames(fData(cds))[fData(cds)$num_cells_expressed >= 10]

cat("Expressed genes used:", length(expressed_genes), "\n")

########################################################
## 6. choose ordering genes
########################################################
## 用4个fibro_state的marker作为ordering genes，更贴合你的状态转换问题

DefaultAssay(fib4) <- "RNA"
Idents(fib4) <- "fibro_state"

markers_all <- FindAllMarkers(
  object = fib4,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

write.csv(
  markers_all,
  file.path(out_dir, "tables", "02_4states_FindAllMarkers_all.csv"),
  row.names = FALSE
)

markers_top <- markers_all %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 100, with_ties = FALSE) %>%
  dplyr::ungroup()

ordering_genes <- unique(markers_top$gene)
ordering_genes <- intersect(ordering_genes, expressed_genes)

## 加入关键基因，避免重要marker被漏掉
key_genes <- c(
  "DLC1",
  "CXCL12", "RBP1", "IGFBP3",
  "C7", "SFRP4", "CCL19", "CCL21",
  "RGS5", "ACTA2", "TAGLN", "SPARCL1",
  "COL1A1", "COL3A1", "LUM", "DCN",
  "POSTN", "CTHRC1", "TIMP1"
)

ordering_genes <- unique(c(ordering_genes, intersect(key_genes, expressed_genes)))

write.csv(
  data.frame(gene = ordering_genes),
  file.path(out_dir, "tables", "03_Monocle2_ordering_genes.csv"),
  row.names = FALSE
)

cat("Ordering genes:", length(ordering_genes), "\n")

cds <- setOrderingFilter(cds, ordering_genes)

pdf(
  file.path(out_dir, "plots", "01_Monocle2_ordering_genes_dispersion.pdf"),
  width = 7,
  height = 5.5
)
print(plot_ordering_genes(cds))
dev.off()

########################################################
## 7. reduce dimension and order cells
########################################################
cds <- reduceDimension(
  cds,
  max_components = 2,
  method = "DDRTree",
  norm_method = "log",
  verbose = TRUE
)

## 第一次先不指定root，让Monocle2自己生成State
cds <- orderCells(cds)

########################################################
## 8. set root state by CXCL12_RBP1_stellate
########################################################
root_group <- "CXCL12_RBP1_stellate"

root_state <- names(which.max(table(
  pData(cds)$State[pData(cds)$fibro_state == root_group]
)))

cat("Root group:", root_group, "\n")
cat("Root state:", root_state, "\n")

cds <- orderCells(cds, root_state = as.numeric(root_state))

saveRDS(
  cds,
  file.path(out_dir, "rds", "02_Monocle2_cds_ordered_CXCL12_root.rds")
)

########################################################
## 9. trajectory plots
########################################################
p_state <- plot_cell_trajectory(cds, color_by = "State") +
  ggtitle("Monocle2 trajectory by State") +
  theme_classic(base_size = 12)

save_pdf(
  p_state,
  file.path(out_dir, "plots", "02_Monocle2_trajectory_by_State.pdf"),
  width = 6.5,
  height = 5.5
)

p_fibro <- plot_cell_trajectory(cds, color_by = "fibro_state") +
  ggtitle("Monocle2 trajectory by fibro_state") +
  theme_classic(base_size = 12)

save_pdf(
  p_fibro,
  file.path(out_dir, "plots", "03_Monocle2_trajectory_by_fibro_state.pdf"),
  width = 7.5,
  height = 5.8
)

p_time <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  ggtitle("Monocle2 pseudotime") +
  theme_classic(base_size = 12)

save_pdf(
  p_time,
  file.path(out_dir, "plots", "04_Monocle2_trajectory_by_Pseudotime.pdf"),
  width = 6.5,
  height = 5.5
)

if (!is.na(source_var)) {
  p_source <- plot_cell_trajectory(cds, color_by = source_var) +
    ggtitle(paste0("Monocle2 trajectory by ", source_var)) +
    theme_classic(base_size = 12)

  save_pdf(
    p_source,
    file.path(out_dir, "plots", paste0("05_Monocle2_trajectory_by_", source_var, ".pdf")),
    width = 6.8,
    height = 5.5
  )
}

########################################################
## 10. DLC1 expression over trajectory
########################################################
if ("DLC1" %in% rownames(cds)) {
  p_dlc1_time <- plot_genes_in_pseudotime(
    cds["DLC1", ],
    color_by = "fibro_state"
  ) +
    ggtitle("DLC1 along pseudotime") +
    theme_classic(base_size = 12)

  save_pdf(
    p_dlc1_time,
    file.path(out_dir, "plots", "06_DLC1_along_pseudotime.pdf"),
    width = 6.8,
    height = 4.8
  )
}

########################################################
## 11. marker genes along pseudotime
########################################################
marker_show <- c(
  "CXCL12", "RBP1",
  "C7", "SFRP4",
  "RGS5", "ACTA2",
  "COL1A1", "POSTN",
  "DLC1"
)

marker_show <- intersect(marker_show, rownames(cds))

if (length(marker_show) > 1) {
  p_marker_time <- plot_genes_in_pseudotime(
    cds[marker_show, ],
    color_by = "fibro_state",
    ncol = 3
  ) +
    theme_classic(base_size = 11)

  save_pdf(
    p_marker_time,
    file.path(out_dir, "plots", "07_key_markers_along_pseudotime.pdf"),
    width = 10,
    height = 8
  )
}

########################################################
## 12. state composition and terminal state judgement
########################################################
state_summary <- pData(cds) %>%
  dplyr::group_by(State, fibro_state) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(State) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

write.csv(
  state_summary,
  file.path(out_dir, "tables", "04_Monocle2_state_fibro_state_composition.csv"),
  row.names = FALSE
)

terminal_state <- pData(cds) %>%
  dplyr::group_by(State) %>%
  dplyr::summarise(
    mean_pseudotime = mean(Pseudotime, na.rm = TRUE),
    median_pseudotime = median(Pseudotime, na.rm = TRUE),
    max_pseudotime = max(Pseudotime, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(mean_pseudotime))

write.csv(
  terminal_state,
  file.path(out_dir, "tables", "05_Monocle2_terminal_state_by_pseudotime.csv"),
  row.names = FALSE
)

fibro_pseudotime_summary <- pData(cds) %>%
  dplyr::group_by(fibro_state) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_pseudotime = mean(Pseudotime, na.rm = TRUE),
    median_pseudotime = median(Pseudotime, na.rm = TRUE),
    min_pseudotime = min(Pseudotime, na.rm = TRUE),
    max_pseudotime = max(Pseudotime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_pseudotime)

write.csv(
  fibro_pseudotime_summary,
  file.path(out_dir, "tables", "06_fibro_state_pseudotime_summary.csv"),
  row.names = FALSE
)

########################################################
## 13. pseudotime boxplot by fibro_state
########################################################
ptime_df <- as.data.frame(pData(cds))

p_box <- ggplot(ptime_df, aes(x = fibro_state, y = Pseudotime, fill = fibro_state)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.4) +
  labs(x = NULL, y = "Pseudotime") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

save_pdf(
  p_box,
  file.path(out_dir, "plots", "08_Pseudotime_by_fibro_state_boxplot.pdf"),
  width = 7.5,
  height = 5.2
)

########################################################
## 14. README
########################################################
writeLines(
  c(
    "Monocle2 4-state trajectory analysis finished.",
    "",
    "Selected states:",
    paste(traj_states, collapse = ", "),
    "",
    "Root group:",
    root_group,
    "",
    "Root state:",
    root_state,
    "",
    "Important outputs:",
    "plots/03_Monocle2_trajectory_by_fibro_state.pdf",
    "plots/04_Monocle2_trajectory_by_Pseudotime.pdf",
    "plots/06_DLC1_along_pseudotime.pdf",
    "plots/07_key_markers_along_pseudotime.pdf",
    "plots/08_Pseudotime_by_fibro_state_boxplot.pdf",
    "tables/04_Monocle2_state_fibro_state_composition.csv",
    "tables/05_Monocle2_terminal_state_by_pseudotime.csv",
    "tables/06_fibro_state_pseudotime_summary.csv"
  ),
  con = file.path(out_dir, "README_Monocle2_4states_CXCL12_root.txt")
)

cat("Monocle2 4-state trajectory analysis finished.\n")






########################################################
## 1.ClusterGVis
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "ClusterGVis",
  "org.Hs.eg.db",
  "clusterProfiler",
  "ComplexHeatmap",
  "circlize"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("org.Hs.eg.db", "clusterProfiler", "ComplexHeatmap")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(p, update = FALSE, ask = FALSE)
    } else if (p == "ClusterGVis") {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
      }
      devtools::install_github("junjunlab/ClusterGVis")
    } else {
      install.packages(p)
    }
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

safe_draw <- function(x) {
  if (inherits(x, "Heatmap") || inherits(x, "HeatmapList")) {
    ComplexHeatmap::draw(x)
  } else {
    print(x)
  }
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

obj_file <- file.path(fib_dir, "rds", "02_fibro_state_final_obj.rds")

if (!file.exists(obj_file)) {
  stop("找不到对象文件: ", obj_file)
}

out_dir <- file.path(fib_dir, "ClusterGVis_4states_Monocle2_order")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)

########################################################
## 2. read fibro object
########################################################
fib_obj <- readRDS(obj_file)

if (!"fibro_state" %in% colnames(fib_obj@meta.data)) {
  stop("对象中没有 fibro_state 列。")
}

if (!"RNA" %in% names(fib_obj@assays)) {
  stop("对象中没有 RNA assay。")
}

########################################################
## 3. select four states based on Monocle2 order
########################################################
traj_states <- c(
  "CXCL12_RBP1_stellate",
  "C7_SFRP4_stromal",
  "Activated_matrix_fibroblast",
  "Perivascular_myofibroblast"
)

traj_labels <- c(
  "CXCL12/RBP1",
  "C7/SFRP4",
  "Activated matrix",
  "Perivascular myofib"
)

names(traj_labels) <- traj_states

keep_cells <- rownames(fib_obj@meta.data)[fib_obj@meta.data$fibro_state %in% traj_states]
fib4 <- subset(fib_obj, cells = keep_cells)

fib4$traj_stage <- traj_labels[as.character(fib4$fibro_state)]
fib4$traj_stage <- factor(fib4$traj_stage, levels = traj_labels)

Idents(fib4) <- "traj_stage"
DefaultAssay(fib4) <- "RNA"
fib4 <- NormalizeData(fib4, verbose = FALSE)

write.csv(
  as.data.frame(table(fib4$traj_stage)),
  file.path(out_dir, "tables", "01_selected_4states_cell_counts.csv"),
  row.names = FALSE
)

saveRDS(
  fib4,
  file.path(out_dir, "rds", "01_fibro_4states_for_ClusterGVis.rds")
)

########################################################
## 4. FindAllMarkers for ClusterGVis
########################################################
markers_all <- Seurat::FindAllMarkers(
  object = fib4,
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

## 兼容 Seurat v4/v5 的 logFC 列名
if (!"avg_log2FC" %in% colnames(markers_all) && "avg_logFC" %in% colnames(markers_all)) {
  markers_all$avg_log2FC <- markers_all$avg_logFC
}

write.csv(
  markers_all,
  file.path(out_dir, "tables", "02_ClusterGVis_4states_FindAllMarkers_all.csv"),
  row.names = FALSE
)

## 每个阶段取 top30 marker
markers_top <- markers_all %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 30, with_ties = FALSE) %>%
  dplyr::ungroup()

########################################################
## 5. force-add key genes if they are not selected
########################################################
key_genes <- c(
  "DLC1",
  "CXCL12", "RBP1", "IGFBP3",
  "C7", "SFRP4", "CCL19", "CCL21",
  "COL1A1", "COL3A1", "LUM", "DCN",
  "POSTN", "CTHRC1", "TIMP1",
  "RGS5", "ACTA2", "TAGLN", "SPARCL1", "MGP"
)

key_genes <- intersect(key_genes, rownames(fib4))

avg_exp <- AverageExpression(
  fib4,
  assays = "RNA",
  group.by = "traj_stage",
  slot = "data"
)$RNA

extra_genes <- setdiff(key_genes, markers_top$gene)

extra_marker_list <- list()

if (length(extra_genes) > 0) {
  for (g in extra_genes) {
    if (!g %in% rownames(avg_exp)) next

    best_cluster <- colnames(avg_exp)[which.max(avg_exp[g, ])]

    extra_marker_list[[g]] <- data.frame(
      p_val = 0.05,
      avg_log2FC = 0.25,
      pct.1 = NA,
      pct.2 = NA,
      p_val_adj = 0.05,
      cluster = best_cluster,
      gene = g,
      stringsAsFactors = FALSE
    )
  }
}

extra_marker_df <- dplyr::bind_rows(extra_marker_list)

markers_use <- dplyr::bind_rows(markers_top, extra_marker_df) %>%
  dplyr::distinct(cluster, gene, .keep_all = TRUE)

markers_use$cluster <- factor(markers_use$cluster, levels = traj_labels)

markers_use <- markers_use %>%
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC))

write.csv(
  markers_use,
  file.path(out_dir, "tables", "03_ClusterGVis_4states_markers_top30_plus_keygenes.csv"),
  row.names = FALSE
)

########################################################
## 6. prepare data for ClusterGVis
## 兼容老版本 ClusterGVis
########################################################

## 关键：ClusterGVis 老版本默认使用 Idents(object)
## 所以前面一定要确认 Idents 是 traj_stage
Idents(fib4) <- "traj_stage"

fib4$traj_stage <- factor(
  fib4$traj_stage,
  levels = traj_labels
)

## ClusterGVis 老版本没有 keepUniqGene，
## 所以这里手动去重，避免一个 gene 出现在多个 cluster 里
markers_use <- markers_use %>%
  dplyr::mutate(
    cluster = as.character(cluster),
    gene = as.character(gene)
  ) %>%
  dplyr::arrange(
    factor(cluster, levels = traj_labels),
    dplyr::desc(avg_log2FC)
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

write.csv(
  markers_use,
  file.path(out_dir, "tables", "03_ClusterGVis_4states_markers_top30_plus_keygenes_unique.csv"),
  row.names = FALSE
)

## 老版本只保留这几个参数
st.data <- ClusterGVis::prepareDataFromscRNA(
  object = fib4,
  diffData = markers_use,
  showAverage = TRUE
)

saveRDS(
  st.data,
  file.path(out_dir, "rds", "02_ClusterGVis_st_data_4states.rds")
)

########################################################
## 7. GO BP enrichment for each ClusterGVis gene cluster
########################################################

enrich_go <- ClusterGVis::enrichCluster(
  object = st.data,
  OrgDb = org.Hs.eg.db,
  type = "BP",
  organism = "hsa",
  pvalueCutoff = 0.5,
  topn = 5,
  seed = 5201314
)

write.csv(
  as.data.frame(enrich_go),
  file.path(out_dir, "tables", "04_ClusterGVis_GO_BP_enrich.csv"),
  row.names = FALSE
)

saveRDS(
  enrich_go,
  file.path(out_dir, "rds", "03_ClusterGVis_GO_BP_enrich.rds")
)

########################################################
## 8. genes to label on heatmap
########################################################

markGenes <- c(
  "CXCL12", "RBP1", "IGFBP3",
  "C7", "SFRP4", "CCL19", "CCL21",
  "COL1A1", "COL3A1", "LUM", "DCN",
  "POSTN", "CTHRC1", "TIMP1",
  "RGS5", "ACTA2", "TAGLN", "SPARCL1",
  "DLC1"
)

markGenes <- intersect(markGenes, markers_use$gene)

write.csv(
  data.frame(markGenes = markGenes),
  file.path(out_dir, "tables", "05_ClusterGVis_markGenes_used.csv"),
  row.names = FALSE
)

########################################################
## 9. colors
########################################################

stage_cols <- c(
  "CXCL12/RBP1"         = "#4DBBD5",
  "C7/SFRP4"            = "#00A087",
  "Activated matrix"    = "#E64B35",
  "Perivascular myofib" = "#7E6148"
)

go_cols <- rep(
  c("#4DBBD5", "#00A087", "#E64B35", "#7E6148"),
  each = 5
)

safe_draw <- function(x) {
  if (inherits(x, "Heatmap") || inherits(x, "HeatmapList")) {
    ComplexHeatmap::draw(x)
  } else {
    print(x)
  }
}

########################################################
## 10. line plot only
########################################################

pdf(
  file.path(out_dir, "plots", "01_ClusterGVis_4states_line.pdf"),
  width = 7,
  height = 5,
  onefile = FALSE
)

p_line <- ClusterGVis::visCluster(
  object = st.data,
  plot.type = "line",
  cluster.order = c(1:4)
)

safe_draw(p_line)
dev.off()

########################################################
## 11. heatmap only
########################################################

pdf(
  file.path(out_dir, "plots", "02_ClusterGVis_4states_heatmap.pdf"),
  width = 7,
  height = 9,
  onefile = FALSE
)

p_heat <- ClusterGVis::visCluster(
  object = st.data,
  plot.type = "heatmap",
  column_names_rot = 45,
  show_row_dend = FALSE,
  markGenes = markGenes,
  markGenes.side = "left",
  cluster.order = c(1:4)
)

safe_draw(p_heat)
dev.off()

########################################################
## 12. final plot: heatmap + line + GO annotation + bar
########################################################

pdf(
  file.path(out_dir, "plots", "03_ClusterGVis_4states_heatmap_GO_line.pdf"),
  width = 14,
  height = 10,
  onefile = FALSE
)

p_both <- ClusterGVis::visCluster(
  object = st.data,
  plot.type = "both",
  column_names_rot = 45,
  show_row_dend = FALSE,
  markGenes = markGenes,
  markGenes.side = "left",
  annoTerm.data = enrich_go,
  line.side = "left",
  cluster.order = c(1:4),
  go.col = go_cols,
  add.bar = TRUE
)

safe_draw(p_both)
dev.off()

########################################################
## 13. optional KEGG enrichment version
########################################################

enrich_kegg <- tryCatch(
  {
    ClusterGVis::enrichCluster(
      object = st.data,
      OrgDb = org.Hs.eg.db,
      type = "KEGG",
      organism = "hsa",
      pvalueCutoff = 0.5,
      topn = 5,
      seed = 5201314
    )
  },
  error = function(e) {
    message("KEGG enrichment failed: ", e$message)
    NULL
  }
)

if (!is.null(enrich_kegg) && nrow(as.data.frame(enrich_kegg)) > 0) {

  write.csv(
    as.data.frame(enrich_kegg),
    file.path(out_dir, "tables", "06_ClusterGVis_KEGG_enrich.csv"),
    row.names = FALSE
  )

  saveRDS(
    enrich_kegg,
    file.path(out_dir, "rds", "04_ClusterGVis_KEGG_enrich.rds")
  )

  pdf(
    file.path(out_dir, "plots", "04_ClusterGVis_4states_heatmap_KEGG_line.pdf"),
    width = 14,
    height = 10,
    onefile = FALSE
  )

  p_both_kegg <- ClusterGVis::visCluster(
    object = st.data,
    plot.type = "both",
    column_names_rot = 45,
    show_row_dend = FALSE,
    markGenes = markGenes,
    markGenes.side = "left",
    annoKegg.data = enrich_kegg,
    line.side = "left",
    cluster.order = c(1:4),
    kegg.col = go_cols,
    add.bar = TRUE
  )

  safe_draw(p_both_kegg)
  dev.off()
}

########################################################
## 14. README
########################################################

writeLines(
  c(
    "ClusterGVis 4-state analysis finished.",
    "",
    "Trajectory order based on Monocle2:",
    paste(traj_labels, collapse = " -> "),
    "",
    "Main output:",
    "plots/03_ClusterGVis_4states_heatmap_GO_line.pdf",
    "",
    "Other outputs:",
    "plots/01_ClusterGVis_4states_line.pdf",
    "plots/02_ClusterGVis_4states_heatmap.pdf",
    "plots/04_ClusterGVis_4states_heatmap_KEGG_line.pdf if KEGG succeeded",
    "",
    "Tables:",
    "tables/03_ClusterGVis_4states_markers_top30_plus_keygenes_unique.csv",
    "tables/04_ClusterGVis_GO_BP_enrich.csv",
    "tables/05_ClusterGVis_markGenes_used.csv"
  ),
  con = file.path(out_dir, "README_ClusterGVis_4states.txt")
)

cat("ClusterGVis 4-state analysis finished.\n")
cat("Main figure:\n")
cat(file.path(out_dir, "plots", "03_ClusterGVis_4states_heatmap_GO_line.pdf"), "\n")








########################################################
## ClusterGVis中C1里面的DLC1的通路富集分析
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "dplyr",
  "ggplot2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("clusterProfiler", "org.Hs.eg.db", "enrichplot")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(p, update = FALSE, ask = FALSE)
    } else {
      install.packages(p)
    }
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

save_pdf <- function(plot_obj, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

cgv_dir <- file.path(fib_dir, "ClusterGVis_4states_Monocle2_order")

st_file <- file.path(cgv_dir, "rds", "02_ClusterGVis_st_data_4states.rds")

out_dir <- file.path(cgv_dir, "C1_GO_DLC1_check")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(st_file)) {
  stop("找不到 ClusterGVis st.data 文件: ", st_file)
}

########################################################
## 2. read ClusterGVis st.data
########################################################
st.data <- readRDS(st_file)

cat("st.data class:\n")
print(class(st.data))

########################################################
## 3. extract C1 genes
########################################################
extract_cgv_cluster_genes <- function(x, target_cluster = "C1") {

  target_num <- gsub("^C", "", target_cluster)

  ## 情况1：st.data 本身是 data.frame
  if (is.data.frame(x)) {

    nms <- colnames(x)

    gene_col <- intersect(
      c("gene", "Gene", "genes", "Genes", "SYMBOL", "symbol", "id", "ID"),
      nms
    )[1]

    cl_col <- intersect(
      c("cluster", "Cluster", "clusters", "Clusters", "cluster_id", "ClusterID"),
      nms
    )[1]

    if (!is.na(gene_col) && !is.na(cl_col)) {
      cl_vec <- as.character(x[[cl_col]])
      genes <- unique(as.character(x[[gene_col]][cl_vec %in% c(target_cluster, target_num)]))
      genes <- genes[!is.na(genes) & genes != ""]
      return(genes)
    }

    ## 如果行名是 gene，cluster 列在表里
    if (!is.na(cl_col)) {
      cl_vec <- as.character(x[[cl_col]])
      genes <- rownames(x)[cl_vec %in% c(target_cluster, target_num)]
      genes <- unique(genes[!is.na(genes) & genes != ""])
      return(genes)
    }
  }

  ## 情况2：st.data 是 list，递归寻找 data.frame
  if (is.list(x)) {
    for (nm in names(x)) {
      xx <- x[[nm]]
      if (is.data.frame(xx)) {
        genes_try <- tryCatch(
          extract_cgv_cluster_genes(xx, target_cluster = target_cluster),
          error = function(e) character(0)
        )
        if (length(genes_try) > 0) {
          message("C1 genes extracted from st.data$", nm)
          return(genes_try)
        }
      }
    }
  }

  stop("没有自动识别出 C1 基因。请运行 names(st.data) 和 str(st.data, max.level = 2) 后把结果发我。")
}

c1_genes <- extract_cgv_cluster_genes(st.data, target_cluster = "C1")
c1_genes <- unique(c1_genes)

cat("C1 gene number:", length(c1_genes), "\n")
cat("DLC1 in C1 genes? ", "DLC1" %in% c1_genes, "\n")

write.csv(
  data.frame(gene = c1_genes),
  file.path(out_dir, "tables", "01_C1_genes_from_ClusterGVis.csv"),
  row.names = FALSE
)

########################################################
## 4. GO enrichment helper functions
########################################################
find_terms_containing_gene <- function(enrich_obj, gene = "DLC1") {

  df <- as.data.frame(enrich_obj)

  if (nrow(df) == 0) {
    return(df)
  }

  if (!"geneID" %in% colnames(df)) {
    return(df[0, ])
  }

  df_gene <- df %>%
    dplyr::filter(grepl(paste0("(^|/)", gene, "($|/)"), geneID))

  return(df_gene)
}

make_enrich_sig_object <- function(enrich_obj, padj_cutoff = 0.05) {

  enrich_sig <- enrich_obj
  enrich_sig@result <- enrich_sig@result %>%
    dplyr::filter(p.adjust < padj_cutoff)

  return(enrich_sig)
}

run_go_one_ont <- function(genes, ont = "BP", out_prefix, out_dir) {

  ego_all <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )

  ego_df <- as.data.frame(ego_all)

  write.csv(
    ego_df,
    file.path(out_dir, "tables", paste0(out_prefix, "_all_terms.csv")),
    row.names = FALSE
  )

  ego_sig_df <- ego_df %>%
    dplyr::filter(p.adjust < 0.05)

  write.csv(
    ego_sig_df,
    file.path(out_dir, "tables", paste0(out_prefix, "_significant_terms_padj005.csv")),
    row.names = FALSE
  )

  ## 重点：筛 geneID 里是否真的包含 DLC1
  dlc1_terms <- find_terms_containing_gene(ego_all, gene = "DLC1")

  write.csv(
    dlc1_terms,
    file.path(out_dir, "tables", paste0(out_prefix, "_terms_containing_DLC1.csv")),
    row.names = FALSE
  )

  dlc1_terms_sig <- dlc1_terms %>%
    dplyr::filter(p.adjust < 0.05)

  write.csv(
    dlc1_terms_sig,
    file.path(out_dir, "tables", paste0(out_prefix, "_significant_terms_containing_DLC1.csv")),
    row.names = FALSE
  )

  ## 画显著 GO，如果没有显著，就画全部 top20
  if (nrow(ego_sig_df) > 0) {

    ego_sig <- make_enrich_sig_object(ego_all, padj_cutoff = 0.05)

    p <- dotplot(ego_sig, showCategory = min(20, nrow(ego_sig_df))) +
      ggtitle(paste0(out_prefix, " significant GO terms")) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 12, face = "bold")
      )

    save_pdf(
      p,
      file.path(out_dir, "plots", paste0(out_prefix, "_significant_dotplot.pdf")),
      width = 8,
      height = 6
    )

  } else if (nrow(ego_df) > 0) {

    p <- dotplot(ego_all, showCategory = min(20, nrow(ego_df))) +
      ggtitle(paste0(out_prefix, " top GO terms, not filtered")) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 12, face = "bold")
      )

    save_pdf(
      p,
      file.path(out_dir, "plots", paste0(out_prefix, "_top20_all_terms_dotplot.pdf")),
      width = 8,
      height = 6
    )
  }

  return(list(
    ego_all = ego_all,
    ego_df = ego_df,
    ego_sig_df = ego_sig_df,
    dlc1_terms = dlc1_terms,
    dlc1_terms_sig = dlc1_terms_sig
  ))
}

########################################################
## 5. GO BP / MF / CC enrichment
########################################################
go_bp_res <- run_go_one_ont(
  genes = c1_genes,
  ont = "BP",
  out_prefix = "GO_BP_C1",
  out_dir = out_dir
)

go_mf_res <- run_go_one_ont(
  genes = c1_genes,
  ont = "MF",
  out_prefix = "GO_MF_C1",
  out_dir = out_dir
)

go_cc_res <- run_go_one_ont(
  genes = c1_genes,
  ont = "CC",
  out_prefix = "GO_CC_C1",
  out_dir = out_dir
)

########################################################
## 6. remove immune-like terms from GO BP for easier viewing
########################################################
immune_patterns <- paste(
  c(
    "antigen", "MHC", "immune", "immunoglobulin",
    "humoral", "leukocyte", "lymphocyte", "T cell", "B cell",
    "complement", "cytokine", "chemokine",
    "inflammatory", "interferon"
  ),
  collapse = "|"
)

go_bp_all <- go_bp_res$ego_df

go_bp_nonimmune <- go_bp_all %>%
  dplyr::filter(!grepl(immune_patterns, Description, ignore.case = TRUE))

write.csv(
  go_bp_nonimmune,
  file.path(out_dir, "tables", "GO_BP_C1_nonimmune_terms_only.csv"),
  row.names = FALSE
)

########################################################
## 7. summary
########################################################
summary_lines <- c(
  "C1 GO enrichment and DLC1 term-check finished.",
  "",
  paste0("C1 gene number: ", length(c1_genes)),
  paste0("DLC1 in C1 gene list: ", "DLC1" %in% c1_genes),
  "",
  paste0("GO BP all terms: ", nrow(go_bp_res$ego_df)),
  paste0("GO BP significant terms, p.adjust < 0.05: ", nrow(go_bp_res$ego_sig_df)),
  paste0("GO BP terms containing DLC1: ", nrow(go_bp_res$dlc1_terms)),
  paste0("GO BP significant terms containing DLC1: ", nrow(go_bp_res$dlc1_terms_sig)),
  "",
  paste0("GO MF terms containing DLC1: ", nrow(go_mf_res$dlc1_terms)),
  paste0("GO CC terms containing DLC1: ", nrow(go_cc_res$dlc1_terms)),
  "",
  "Important interpretation:",
  "GO enrichment is performed using all C1 genes, not DLC1 alone.",
  "If GO_*_terms_containing_DLC1.csv is empty, it means DLC1 is in the C1 expression module but does not directly contribute to those enriched GO terms.",
  "In that case, the immune-related C1 enrichment is likely driven by other C1 genes rather than DLC1 itself."
)

writeLines(
  summary_lines,
  con = file.path(out_dir, "C1_GO_DLC1_check_summary.txt")
)

cat(paste(summary_lines, collapse = "\n"))
cat("\n")




########################################################
## ClusterGVis的美化和调整
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "dplyr",
  "ClusterGVis",
  "ComplexHeatmap"
)

for (p in pkg_needed) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

safe_draw <- function(x) {
  if (inherits(x, "Heatmap") || inherits(x, "HeatmapList")) {
    ComplexHeatmap::draw(x)
  } else {
    print(x)
  }
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

cgv_dir <- file.path(fib_dir, "ClusterGVis_4states_Monocle2_order")

st_file <- file.path(cgv_dir, "rds", "02_ClusterGVis_st_data_4states.rds")
go_file <- file.path(cgv_dir, "rds", "03_ClusterGVis_GO_BP_enrich.rds")

dlc1_go_file <- file.path(
  cgv_dir,
  "C1_GO_DLC1_check",
  "tables",
  "GO_BP_C1_significant_terms_containing_DLC1.csv"
)

out_dir <- file.path(cgv_dir, "plots_C1_DLC1_GO_replaced_v2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(st_file)) stop("找不到 st.data 文件: ", st_file)
if (!file.exists(go_file)) stop("找不到 ClusterGVis GO RDS 文件: ", go_file)
if (!file.exists(dlc1_go_file)) stop("找不到 DLC1 C1 GO 文件: ", dlc1_go_file)

########################################################
## 2. read data
########################################################
st.data <- readRDS(st_file)
enrich_go_old <- readRDS(go_file)
enrich_go_df <- as.data.frame(enrich_go_old)

dlc1_c1_go <- read.csv(
  dlc1_go_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Original ClusterGVis GO columns:\n")
print(colnames(enrich_go_df))

cat("DLC1 GO columns:\n")
print(colnames(dlc1_c1_go))

########################################################
## 3. detect cluster column
########################################################
cluster_col <- intersect(
  c("cluster", "Cluster", "clusters", "Clusters", "group", "Group"),
  colnames(enrich_go_df)
)[1]

if (is.na(cluster_col)) {
  stop("没有识别到 cluster/group 列，请查看 colnames(enrich_go_df)。")
}

enrich_go_df[[cluster_col]] <- as.character(enrich_go_df[[cluster_col]])

cat("cluster column:", cluster_col, "\n")
cat("cluster labels:\n")
print(unique(enrich_go_df[[cluster_col]]))

########################################################
## 4. find C1 label
########################################################
possible_c1 <- c("C1", "1", "cluster 1", "Cluster 1")

c1_label <- unique(enrich_go_df[[cluster_col]][
  enrich_go_df[[cluster_col]] %in% possible_c1
])[1]

if (is.na(c1_label)) {
  c1_label <- unique(enrich_go_df[[cluster_col]])[1]
  message("没有找到标准 C1 标签，默认使用第一个 cluster 作为 C1: ", c1_label)
}

cat("C1 label used:", c1_label, "\n")

########################################################
## 5. use original C1 rows as template
########################################################
c1_old <- enrich_go_df %>%
  dplyr::filter(.data[[cluster_col]] == c1_label)

if (nrow(c1_old) == 0) {
  stop("原始 GO 结果里面没有找到 C1 行。")
}

n_replace <- min(nrow(c1_old), nrow(dlc1_c1_go))

## 用原始 C1 的前 n 行作为模板，这样 bar 所需结构不会丢
dlc1_rows <- c1_old[seq_len(n_replace), , drop = FALSE]

########################################################
## 6. overwrite C1 rows with DLC1-containing GO terms
########################################################
for (i in seq_len(n_replace)) {

  ## cluster 保持 C1
  dlc1_rows[[cluster_col]][i] <- c1_label

  ## GO ID
  if ("ID" %in% colnames(dlc1_c1_go) && "ID" %in% colnames(dlc1_rows)) {
    dlc1_rows$ID[i] <- dlc1_c1_go$ID[i]
  }

  ## Description / term name
  if ("Description" %in% colnames(dlc1_c1_go)) {
    term_cols <- intersect(
      c("Description", "Term", "term", "term_name", "name", "Name"),
      colnames(dlc1_rows)
    )
    for (tc in term_cols) {
      dlc1_rows[[tc]][i] <- dlc1_c1_go$Description[i]
    }
  }

  ## GeneRatio / BgRatio
  if ("GeneRatio" %in% colnames(dlc1_c1_go) && "GeneRatio" %in% colnames(dlc1_rows)) {
    dlc1_rows$GeneRatio[i] <- dlc1_c1_go$GeneRatio[i]
  }
  if ("BgRatio" %in% colnames(dlc1_c1_go) && "BgRatio" %in% colnames(dlc1_rows)) {
    dlc1_rows$BgRatio[i] <- dlc1_c1_go$BgRatio[i]
  }

  ## pvalue / p.adjust / qvalue
  for (stat_col in c("pvalue", "p.adjust", "qvalue")) {
    if (stat_col %in% colnames(dlc1_c1_go) && stat_col %in% colnames(dlc1_rows)) {
      dlc1_rows[[stat_col]][i] <- as.numeric(dlc1_c1_go[[stat_col]][i])
    }
  }

  ## geneID
  if ("geneID" %in% colnames(dlc1_c1_go) && "geneID" %in% colnames(dlc1_rows)) {
    dlc1_rows$geneID[i] <- dlc1_c1_go$geneID[i]
  }

  ## Count
  if ("Count" %in% colnames(dlc1_c1_go) && "Count" %in% colnames(dlc1_rows)) {
    dlc1_rows$Count[i] <- as.numeric(dlc1_c1_go$Count[i])
  }
}

########################################################
## 7. if Count is missing, calculate from geneID
########################################################
if ("Count" %in% colnames(dlc1_rows) && "geneID" %in% colnames(dlc1_rows)) {
  bad_count <- is.na(dlc1_rows$Count) | dlc1_rows$Count <= 0

  if (any(bad_count)) {
    dlc1_rows$Count[bad_count] <- sapply(
      strsplit(as.character(dlc1_rows$geneID[bad_count]), "/"),
      length
    )
  }
}

########################################################
## 8. force numeric columns
########################################################
for (num_col in c("pvalue", "p.adjust", "qvalue", "Count")) {
  if (num_col %in% colnames(enrich_go_df)) {
    enrich_go_df[[num_col]] <- as.numeric(enrich_go_df[[num_col]])
  }
  if (num_col %in% colnames(dlc1_rows)) {
    dlc1_rows[[num_col]] <- as.numeric(dlc1_rows[[num_col]])
  }
}

########################################################
## 9. replace C1 terms, keep C2-C4 unchanged
########################################################
enrich_go_replaced <- enrich_go_df %>%
  dplyr::filter(.data[[cluster_col]] != c1_label) %>%
  dplyr::bind_rows(dlc1_rows)

## 恢复 cluster 顺序
enrich_go_replaced[[cluster_col]] <- factor(
  enrich_go_replaced[[cluster_col]],
  levels = unique(enrich_go_df[[cluster_col]])
)

########################################################
## 10. export check table
########################################################
write.csv(
  enrich_go_replaced,
  file.path(cgv_dir, "tables", "04_ClusterGVis_GO_BP_enrich_C1_replaced_by_DLC1_terms_v2.csv"),
  row.names = FALSE
)

write.csv(
  dlc1_rows,
  file.path(cgv_dir, "tables", "04_ClusterGVis_C1_DLC1_replacement_rows_v2_check.csv"),
  row.names = FALSE
)

cat("C1 replacement rows used for plotting:\n")
print(dlc1_rows[, intersect(c(cluster_col, "Description", "pvalue", "p.adjust", "Count", "geneID"), colnames(dlc1_rows))])

########################################################
## 11. mark genes
########################################################
markGenes <- c(
  "DLC1",
  "CXCL12", "RBP1", "IGFBP3",
  "C7", "SFRP4", "CCL19", "CCL21",
  "COL1A1", "COL3A1", "LUM", "POSTN", "CTHRC1", "TIMP1",
  "RGS5", "ACTA2", "TAGLN", "SPARCL1"
)

########################################################
## 12. colors
########################################################
go_cols <- rep(
  c("#4DBBD5", "#00A087", "#E64B35", "#7E6148"),
  each = 5
)

########################################################
## 13. redraw figure
########################################################
pdf(
  file.path(out_dir, "03_ClusterGVis_4states_heatmap_GO_line_C1_replaced_by_DLC1_terms_v2.pdf"),
  width = 14,
  height = 10,
  onefile = FALSE
)

p_both_replaced <- ClusterGVis::visCluster(
  object = st.data,
  plot.type = "both",
  column_names_rot = 45,
  show_row_dend = FALSE,
  markGenes = markGenes,
  markGenes.side = "left",
  annoTerm.data = enrich_go_replaced,
  line.side = "left",
  cluster.order = c(1:4),
  go.col = go_cols,
  add.bar = TRUE
)

safe_draw(p_both_replaced)
dev.off()

########################################################
## 14. README
########################################################
writeLines(
  c(
    "ClusterGVis C1 GO replacement v2 finished.",
    "",
    "C1 GO terms were replaced by DLC1-containing GO BP terms.",
    "This version uses original C1 rows as templates to preserve columns required by add.bar = TRUE.",
    "",
    "New figure:",
    "plots_C1_DLC1_GO_replaced_v2/03_ClusterGVis_4states_heatmap_GO_line_C1_replaced_by_DLC1_terms_v2.pdf",
    "",
    "Check table:",
    "tables/04_ClusterGVis_C1_DLC1_replacement_rows_v2_check.csv"
  ),
  con = file.path(out_dir, "README_C1_GO_replaced_by_DLC1_terms_v2.txt")
)

cat("Finished.\n")
cat("New figure:\n")
cat(file.path(out_dir, "03_ClusterGVis_4states_heatmap_GO_line_C1_replaced_by_DLC1_terms_v2.pdf"), "\n")








########################################################
## 4 CellChat
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "dplyr", "ggplot2", "tidyr", "stringr",
  "patchwork", "scales", "readr", "Seurat", "Matrix"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

save_pdf <- function(p, file, width = 7, height = 6) {
  ggsave(file, p, width = width, height = height, device = cairo_pdf)
}

theme_clean <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = "grey25", hjust = 0, size = base_size),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", color = "black"),
      plot.margin = margin(8, 16, 8, 8)
    )
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

obj_file <- file.path(fib_dir, "rds", "02_fibro_state_final_obj.rds")

cc_dir <- file.path(fib_dir, "CellChat_C2C3C4_to_C1_DLC1_transition")

all_file <- file.path(
  cc_dir, "tables", "04_MAIN_C2C3C4_to_C1_incoming_LR_pairs.csv"
)

dlc1_file <- file.path(
  cc_dir, "tables", "06_MAIN_C2C3C4_to_C1_DLC1_GO_context_LR_pairs.csv"
)

focus_file <- file.path(
  cc_dir, "tables", "07_MAIN_C2C3C4_to_C1_focus_migration_adhesion_LR_pairs.csv"
)

strength_file <- file.path(
  cc_dir, "tables", "08_MAIN_C2C3C4_to_C1_total_strength_by_source.csv"
)

out_dir <- file.path(cc_dir, "FINAL_Figure_CellChat_screening_logic_v4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

for (f in c(all_file, dlc1_file, focus_file, strength_file)) {
  if (!file.exists(f)) stop("找不到文件: ", f)
}

########################################################
## 2. labels
########################################################
source_levels <- c(
  "C2_C7_SFRP4",
  "C3_Activated_matrix",
  "C4_Perivascular_myofib"
)

source_labels <- c(
  "C2_C7_SFRP4" = "C2\nC7/SFRP4",
  "C3_Activated_matrix" = "C3\nActivated\nmatrix",
  "C4_Perivascular_myofib" = "C4\nPerivascular\nmyofib"
)

state_order <- c(
  "C1_CXCL12_RBP1",
  "C2_C7_SFRP4",
  "C3_Activated_matrix",
  "C4_Perivascular_myofib"
)

state_labels <- c(
  "C1_CXCL12_RBP1" = "C1",
  "C2_C7_SFRP4" = "C2",
  "C3_Activated_matrix" = "C3",
  "C4_Perivascular_myofib" = "C4"
)

state_cols <- c(
  "C1_CXCL12_RBP1" = "#009E73",
  "C2_C7_SFRP4" = "#E69F00",
  "C3_Activated_matrix" = "#D55E00",
  "C4_Perivascular_myofib" = "#0072B2"
)

axis_cols <- c(
  "DLC1-context screening" = "#B2182B",
  "ECM–integrin screening" = "#2166AC",
  "Final primary axis" = "#B2182B",
  "Final ECM–integrin context" = "#2166AC"
)

########################################################
## 3. helper
########################################################
make_pair_label <- function(df) {
  if ("interaction_name_2" %in% colnames(df)) {
    lab <- df$interaction_name_2
  } else {
    lab <- paste0(df$ligand, " → ", df$receptor)
  }
  lab <- stringr::str_replace_all(lab, "_", "+")
  lab <- stringr::str_replace_all(lab, " - ", " → ")
  lab
}

has_gene <- function(x, pattern) {
  grepl(pattern, x, ignore.case = TRUE)
}

read_lr <- function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"interaction_name" %in% colnames(df)) {
    df$interaction_name <- paste0(df$ligand, "_", df$receptor)
  }
  df$pair_label <- make_pair_label(df)
  df <- df %>%
    dplyr::filter(
      source %in% source_levels,
      target == "C1_CXCL12_RBP1"
    ) %>%
    dplyr::mutate(
      source = factor(source, levels = source_levels),
      source_label = factor(
        source_labels[as.character(source)],
        levels = source_labels[source_levels]
      )
    )
  df
}

########################################################
## 4. read tables
########################################################
lr_all   <- read_lr(all_file)
lr_dlc1  <- read_lr(dlc1_file)
lr_focus <- read_lr(focus_file)

########################################################
## 5. Screening evidence A1: DLC1-context candidates
########################################################
dlc1_screen <- lr_dlc1 %>%
  dplyr::filter(
    ligand %in% c("APP", "MIF", "CXCL12") |
      grepl("APP|MIF|CXCL12", interaction_name, ignore.case = TRUE)
  ) %>%
  dplyr::mutate(
    screen_class = "DLC1-context screening",
    selected_reason = dplyr::case_when(
      ligand == "MIF" & has_gene(receptor, "CD74") & has_gene(receptor, "CXCR4") ~
        "Primary candidate: DLC1-context + C3 feedback + expression support",
      TRUE ~ "DLC1-context auxiliary candidate"
    )
  )

dlc1_screen$pair_label <- factor(
  dlc1_screen$pair_label,
  levels = rev(c(
    "APP → CD74",
    "MIF → (CD74+CXCR4)",
    "CXCL12 → CXCR4"
  ))
)

write.csv(
  dlc1_screen,
  file.path(out_dir, "tables", "01_screening_DLC1_context_candidates.csv"),
  row.names = FALSE
)

p_a1 <- ggplot(
  dlc1_screen,
  aes(x = source_label, y = pair_label, size = prob, fill = prob)
) +
  geom_point(shape = 21, color = "black", stroke = 0.35, alpha = 0.95) +
  scale_fill_gradientn(
    colors = c("#FDE0DD", "#FA9FB5", "#C0002B"),
    name = "Communication\nprobability"
  ) +
  scale_size_continuous(
    range = c(2.5, 9),
    name = "Communication\nprobability"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Step 1A. DLC1-context ligand–receptor screening",
    subtitle = "MIF–CD74/CXCR4 is selected from DLC1-context candidates, not chosen arbitrarily"
  ) +
  theme_clean(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

save_pdf(
  p_a1,
  file.path(out_dir, "plots", "Fig1A_screening_DLC1_context_candidates.pdf"),
  width = 8.2,
  height = 4.2
)

########################################################
## 6. Screening evidence A2: ECM-integrin candidates from focus table
########################################################
ecm_screen_all <- lr_focus %>%
  dplyr::filter(
    (
      grepl("^COL", ligand) |
        ligand %in% c("SPP1", "FN1") |
        grepl("^LAM", ligand)
    ) &
      has_gene(receptor, "ITG") &
      !has_gene(receptor, "SDC4") &
      !has_gene(receptor, "CD47")
  ) %>%
  dplyr::mutate(
    screen_class = "ECM–integrin screening"
  )

## 用 C3 的概率优先选 top 16，但保留这些 pair 在 C2/C3/C4 的所有点
top_ecm_pairs <- ecm_screen_all %>%
  dplyr::filter(source == "C3_Activated_matrix") %>%
  dplyr::arrange(dplyr::desc(prob)) %>%
  dplyr::distinct(interaction_name, pair_label, ligand, receptor, .keep_all = TRUE) %>%
  dplyr::slice_head(n = 16) %>%
  dplyr::select(interaction_name, pair_label, ligand, receptor)

ecm_screen <- ecm_screen_all %>%
  dplyr::semi_join(
    top_ecm_pairs,
    by = c("interaction_name", "pair_label", "ligand", "receptor")
  )

pair_order <- ecm_screen %>%
  dplyr::group_by(pair_label) %>%
  dplyr::summarise(max_prob = max(prob, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(max_prob) %>%
  dplyr::pull(pair_label)

ecm_screen$pair_label <- factor(ecm_screen$pair_label, levels = pair_order)

write.csv(
  ecm_screen,
  file.path(out_dir, "tables", "02_screening_ECM_integrin_candidates.csv"),
  row.names = FALSE
)

height_ecm <- max(6.2, 0.34 * length(unique(ecm_screen$pair_label)) + 2.4)

p_a2 <- ggplot(
  ecm_screen,
  aes(x = source_label, y = pair_label, size = prob, fill = prob)
) +
  geom_point(shape = 21, color = "black", stroke = 0.30, alpha = 0.95) +
  scale_fill_gradientn(
    colors = c("#DEEBF7", "#6BAED6", "#08519C"),
    name = "Communication\nprobability"
  ) +
  scale_size_continuous(
    range = c(2.0, 8.5),
    name = "Communication\nprobability"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Step 1B. Focused ECM–integrin screening",
    subtitle = "ECM–integrin interactions dominate the migration/adhesion-focused LR screening"
  ) +
  theme_clean(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, face = "bold"),
    axis.text.y = element_text(size = 8.5, face = "bold"),
    legend.position = "right"
  )

save_pdf(
  p_a2,
  file.path(out_dir, "plots", "Fig1B_screening_ECM_integrin_candidates.pdf"),
  width = 9.2,
  height = height_ecm
)

########################################################
## 7. Combine screening figure
########################################################
p_screening <- p_a1 / p_a2 + patchwork::plot_layout(heights = c(1, 2.2))

save_pdf(
  p_screening,
  file.path(out_dir, "plots", "Fig1_screening_logic_from_01_to_03.pdf"),
  width = 9.4,
  height = height_ecm + 4.3
)

########################################################
## 8. Step 2: incoming strength toward C1
########################################################
strength_df <- read.csv(
  strength_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  dplyr::filter(source %in% source_levels) %>%
  dplyr::mutate(
    source = factor(source, levels = source_levels),
    source_full = dplyr::case_when(
      source == "C2_C7_SFRP4" ~ "C2: C7/SFRP4",
      source == "C3_Activated_matrix" ~ "C3: Activated matrix",
      source == "C4_Perivascular_myofib" ~ "C4: Perivascular myofib"
    ),
    highlight = ifelse(source == "C3_Activated_matrix", "C3", "Other")
  ) %>%
  dplyr::arrange(total_prob)

write.csv(
  strength_df,
  file.path(out_dir, "tables", "03_incoming_strength_to_C1.csv"),
  row.names = FALSE
)

p_b <- ggplot(
  strength_df,
  aes(x = total_prob, y = reorder(source_full, total_prob))
) +
  geom_col(
    aes(fill = highlight),
    width = 0.58,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.3f", total_prob)),
    hjust = -0.18,
    size = 3.8,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("C3" = "#D55E00", "Other" = "#D9D9D9"),
    guide = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(
    x = "Total communication probability toward C1",
    y = NULL,
    title = "Step 2. C3 is the dominant incoming source toward C1",
    subtitle = "This step justifies focusing on the C3→C1 direction"
  ) +
  theme_clean(base_size = 11)

save_pdf(
  p_b,
  file.path(out_dir, "plots", "Fig2_incoming_strength_to_C1.pdf"),
  width = 7.2,
  height = 3.8
)

########################################################
## 9. Step 3: final C3 -> C1 candidates
########################################################
c3_mif <- dlc1_screen %>%
  dplyr::filter(
    source == "C3_Activated_matrix",
    ligand == "MIF",
    has_gene(receptor, "CD74"),
    has_gene(receptor, "CXCR4")
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::mutate(final_class = "Final primary axis")

c3_ecm <- ecm_screen %>%
  dplyr::filter(source == "C3_Activated_matrix") %>%
  dplyr::arrange(dplyr::desc(prob)) %>%
  dplyr::distinct(ligand, receptor, .keep_all = TRUE) %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::mutate(final_class = "Final ECM–integrin context")

final_c3 <- dplyr::bind_rows(c3_mif, c3_ecm) %>%
  dplyr::mutate(
    pair_label_clean = paste0(ligand, " → ", stringr::str_replace_all(receptor, "_", "+")),
    final_class = factor(
      final_class,
      levels = c("Final primary axis", "Final ECM–integrin context")
    )
  ) %>%
  dplyr::arrange(final_class, dplyr::desc(prob))

final_c3$pair_label_clean <- factor(
  final_c3$pair_label_clean,
  levels = rev(final_c3$pair_label_clean)
)

write.csv(
  final_c3,
  file.path(out_dir, "tables", "04_final_C3_to_C1_selected_candidate_axes.csv"),
  row.names = FALSE
)

p_c <- ggplot(
  final_c3,
  aes(x = prob, y = pair_label_clean, fill = final_class)
) +
  geom_segment(
    aes(x = 0, xend = prob, yend = pair_label_clean),
    color = "grey78",
    linewidth = 1.1,
    lineend = "round"
  ) +
  geom_point(
    shape = 21,
    size = 6.2,
    color = "black",
    stroke = 0.35
  ) +
  scale_fill_manual(values = axis_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(
    x = "Communication probability",
    y = NULL,
    title = "Step 3. Final C3→C1 candidate axes",
    subtitle = "Final candidates are selected after DLC1-context and ECM–integrin screening"
  ) +
  theme_clean(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "bold", size = 9.5),
    legend.position = "top"
  )

save_pdf(
  p_c,
  file.path(out_dir, "plots", "Fig3_final_C3_to_C1_candidate_axes.pdf"),
  width = 8.0,
  height = 4.6
)

########################################################
## 10. Expression validation
########################################################
if (!file.exists(obj_file)) stop("找不到 Seurat 对象: ", obj_file)

fib_obj <- readRDS(obj_file)

traj_order <- c(
  "CXCL12_RBP1_stellate",
  "C7_SFRP4_stromal",
  "Activated_matrix_fibroblast",
  "Perivascular_myofibroblast"
)

traj_labels <- c(
  "C1_CXCL12_RBP1",
  "C2_C7_SFRP4",
  "C3_Activated_matrix",
  "C4_Perivascular_myofib"
)

names(traj_labels) <- traj_order

fib4 <- subset(
  fib_obj,
  cells = rownames(fib_obj@meta.data)[fib_obj$fibro_state %in% traj_order]
)

fib4$fibro_comm_state <- traj_labels[as.character(fib4$fibro_state)]
fib4$fibro_comm_state <- factor(fib4$fibro_comm_state, levels = traj_labels)

DefaultAssay(fib4) <- "RNA"
fib4 <- NormalizeData(fib4, verbose = FALSE)

get_data_safe <- function(obj) {
  tryCatch(
    Seurat::GetAssayData(obj, assay = "RNA", layer = "data"),
    error = function(e) Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
  )
}

expr_summary <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  mat <- get_data_safe(obj)
  mat <- mat[genes, , drop = FALSE]
  group <- as.character(obj$fibro_comm_state[colnames(mat)])

  out <- lapply(genes, function(g) {
    data.frame(
      gene = g,
      group = group,
      expr = as.numeric(mat[g, ]),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::group_by(gene, group) %>%
      dplyr::summarise(
        avg_expr = mean(expr),
        pct_expr = mean(expr > 0) * 100,
        .groups = "drop"
      )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(avg_expr_z = as.numeric(scale(avg_expr))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      group = factor(group, levels = traj_labels),
      group_label = factor(state_labels[as.character(group)], levels = state_labels[state_order])
    )
  out
}

########################################################
## 10A. MIF expression support
########################################################
mif_expr <- expr_summary(fib4, c("MIF", "CD74", "CXCR4"))
mif_expr$gene <- factor(mif_expr$gene, levels = rev(c("MIF", "CD74", "CXCR4")))

write.csv(
  mif_expr,
  file.path(out_dir, "tables", "05_expression_MIF_CD74_CXCR4.csv"),
  row.names = FALSE
)

p_d <- ggplot(
  mif_expr,
  aes(x = group_label, y = gene, size = pct_expr, fill = avg_expr_z)
) +
  geom_point(shape = 21, color = "black", stroke = 0.30, alpha = 0.95) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    name = "Scaled\naverage\nexpression"
  ) +
  scale_size_continuous(
    range = c(2.0, 8.5),
    limits = c(0, 100),
    name = "Percent\nexpressed"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Step 4. Expression support for MIF–CD74/CXCR4",
    subtitle = "MIF ligand is enriched in C3; CD74/CXCR4 receptors are available in C1"
  ) +
  theme_clean(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

save_pdf(
  p_d,
  file.path(out_dir, "plots", "Fig4_expression_support_MIF_CD74_CXCR4.pdf"),
  width = 7.4,
  height = 4.2
)

########################################################
## 10B. ECM-integrin expression support
########################################################
ecm_ligands <- c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "SPP1", "FN1")
integrins <- c("ITGA1", "ITGA9", "ITGB1")

ecm_expr <- expr_summary(fib4, c(ecm_ligands, integrins)) %>%
  dplyr::mutate(
    gene_role = dplyr::case_when(
      gene %in% ecm_ligands ~ "ECM ligands",
      gene %in% integrins ~ "Integrin receptors"
    )
  )

ecm_expr$gene_role <- factor(ecm_expr$gene_role, levels = c("ECM ligands", "Integrin receptors"))

ecm_expr <- ecm_expr %>%
  dplyr::mutate(
    gene = dplyr::case_when(
      gene_role == "ECM ligands" ~ factor(gene, levels = rev(ecm_ligands)),
      gene_role == "Integrin receptors" ~ factor(gene, levels = rev(integrins))
    )
  )

write.csv(
  ecm_expr,
  file.path(out_dir, "tables", "06_expression_ECM_integrin.csv"),
  row.names = FALSE
)

p_e <- ggplot(
  ecm_expr,
  aes(x = group_label, y = gene, size = pct_expr, fill = avg_expr_z)
) +
  geom_point(shape = 21, color = "black", stroke = 0.28, alpha = 0.95) +
  facet_grid(
    gene_role ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    name = "Scaled\naverage\nexpression"
  ) +
  scale_size_continuous(
    range = c(1.8, 7.8),
    limits = c(0, 100),
    name = "Percent\nexpressed"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Step 5. Expression support for ECM–integrin context",
    subtitle = "ECM ligands are enriched in C2/C3; integrin receptors are available in C1"
  ) +
  theme_clean(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 9),
    strip.text = element_text(face = "bold")
  )

save_pdf(
  p_e,
  file.path(out_dir, "plots", "Fig5_expression_support_ECM_integrin.pdf"),
  width = 8.0,
  height = 6.6
)

########################################################
## 11. model summary table
########################################################
model_table <- data.frame(
  Step = c(
    "Screening 1",
    "Screening 2",
    "Direction selection",
    "Final primary axis",
    "Final background axis",
    "Expression validation"
  ),
  Evidence = c(
    "DLC1-context LR screening retained APP-CD74, MIF-CD74/CXCR4 and CXCL12-CXCR4.",
    "Focus LR screening was dominated by collagen/ECM ligand to integrin receptor interactions.",
    "C3_Activated_matrix showed the strongest total communication probability toward C1.",
    "MIF-CD74/CXCR4 was prioritized because it is a DLC1-context axis and occurs in the dominant C3-to-C1 direction.",
    "ECM-integrin was retained as a matrix-adhesion context rather than a direct DLC1 axis.",
    "Expression patterns support C3-high MIF and C1-available CD74/CXCR4, plus C2/C3 ECM ligands and C1 integrins."
  ),
  Interpretation = c(
    "MIF is selected from a defined DLC1-context candidate set.",
    "ECM-integrin is selected because it is a recurrent communication class in the focused screen.",
    "C3 is the most plausible upstream stromal feedback source.",
    "MIF-CD74/CXCR4 is the main candidate communication axis.",
    "ECM-integrin provides an adhesive stromal microenvironment supporting C1 remodeling.",
    "Expression validation supports the predicted CellChat direction."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  model_table,
  file.path(out_dir, "tables", "07_CellChat_selection_logic_model_summary.csv"),
  row.names = FALSE
)

cat("CellChat screening logic figure finished.\n")
cat("Output dir:\n")
cat(out_dir, "\n")






########################################################
## 修改文件3和文件5，并输出调整以后的初始关系图
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

########################################################
## 0. packages
########################################################
pkg_needed <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "tidyr",
  "stringr",
  "Matrix",
  "patchwork",
  "scales",
  "readr"
)

for (p in pkg_needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  message("CellChat 未安装：circle plot 会跳过，但 heatmap 和其余图会正常输出。")
} else {
  suppressPackageStartupMessages(library(CellChat))
}

save_pdf <- function(p, file, width = 7, height = 6) {
  ggsave(
    filename = file,
    plot = p,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

theme_final <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = "grey25", hjust = 0, size = base_size),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", color = "black"),
      plot.margin = margin(8, 18, 8, 8)
    )
}

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

obj_file <- file.path(fib_dir, "rds", "02_fibro_state_final_obj.rds")

cc_dir <- file.path(fib_dir, "CellChat_C2C3C4_to_C1_DLC1_transition")

all_file <- file.path(
  cc_dir, "tables", "04_MAIN_C2C3C4_to_C1_incoming_LR_pairs.csv"
)

dlc1_file <- file.path(
  cc_dir, "tables", "06_MAIN_C2C3C4_to_C1_DLC1_GO_context_LR_pairs.csv"
)

focus_file <- file.path(
  cc_dir, "tables", "07_MAIN_C2C3C4_to_C1_focus_migration_adhesion_LR_pairs.csv"
)

strength_file <- file.path(
  cc_dir, "tables", "08_MAIN_C2C3C4_to_C1_total_strength_by_source.csv"
)

cellchat_rds_candidates <- c(
  file.path(cc_dir, "rds", "02_CellChat_4states_full_object.rds"),
  file.path(cc_dir, "rds", "02_CellChat_C2C3C4_to_C1_DLC1_transition.rds"),
  file.path(cc_dir, "rds", "03_CellChat_fibro_internal_4states_with_centrality.rds"),
  file.path(cc_dir, "CellChat_4states_full_object.rds")
)

cellchat_file <- cellchat_rds_candidates[file.exists(cellchat_rds_candidates)][1]

out_dir <- file.path(cc_dir, "FINAL_Figure_CellChat_FIX_v5")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)

for (f in c(all_file, dlc1_file, focus_file, strength_file, obj_file)) {
  if (!file.exists(f)) stop("找不到文件: ", f)
}

########################################################
## 2. labels and colors
########################################################
source_levels <- c(
  "C2_C7_SFRP4",
  "C3_Activated_matrix",
  "C4_Perivascular_myofib"
)

source_labels <- c(
  "C2_C7_SFRP4" = "C2\nC7/SFRP4",
  "C3_Activated_matrix" = "C3\nActivated\nmatrix",
  "C4_Perivascular_myofib" = "C4\nPerivascular\nmyofib"
)

state_order <- c(
  "C1_CXCL12_RBP1",
  "C2_C7_SFRP4",
  "C3_Activated_matrix",
  "C4_Perivascular_myofib"
)

state_short <- c(
  "C1_CXCL12_RBP1" = "C1",
  "C2_C7_SFRP4" = "C2",
  "C3_Activated_matrix" = "C3",
  "C4_Perivascular_myofib" = "C4"
)

state_full <- c(
  "C1_CXCL12_RBP1" = "C1: CXCL12/RBP1",
  "C2_C7_SFRP4" = "C2: C7/SFRP4",
  "C3_Activated_matrix" = "C3: Activated matrix",
  "C4_Perivascular_myofib" = "C4: Perivascular myofib"
)

state_cols <- c(
  "C1_CXCL12_RBP1" = "#009E73",
  "C2_C7_SFRP4" = "#E69F00",
  "C3_Activated_matrix" = "#D55E00",
  "C4_Perivascular_myofib" = "#0072B2"
)

axis_cols <- c(
  "Primary MIF axis" = "#B2182B",
  "ECM-integrin context" = "#2166AC"
)

########################################################
## 3. helper functions
########################################################
make_pair_label <- function(df) {
  if ("interaction_name_2" %in% colnames(df)) {
    lab <- df$interaction_name_2
  } else {
    lab <- paste0(df$ligand, " -> ", df$receptor)
  }
  lab <- stringr::str_replace_all(lab, "_", "+")
  lab <- stringr::str_replace_all(lab, " - ", " -> ")
  lab <- stringr::str_replace_all(lab, "->", "→")
  lab
}

has_gene <- function(x, pattern) {
  grepl(pattern, x, ignore.case = TRUE)
}

read_lr <- function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)

  if (!"interaction_name" %in% colnames(df)) {
    df$interaction_name <- paste0(df$ligand, "_", df$receptor)
  }

  df$pair_label <- make_pair_label(df)

  p_col <- intersect(
    c("pval", "p.value", "p_value", "pvalue", "P.Value", "PValue"),
    colnames(df)
  )[1]

  if (!is.na(p_col)) {
    df$p_value <- as.numeric(df[[p_col]])
  } else {
    df$p_value <- NA_real_
  }

  df <- df %>%
    dplyr::filter(
      source %in% source_levels,
      target == "C1_CXCL12_RBP1"
    ) %>%
    dplyr::mutate(
      source = factor(source, levels = source_levels),
      source_label = factor(
        source_labels[as.character(source)],
        levels = source_labels[source_levels]
      )
    )

  df
}

plot_lr_bubble_resized <- function(df, title, subtitle, outfile, color_mode = c("red", "blue")) {

  color_mode <- match.arg(color_mode)

  pair_order <- df %>%
    dplyr::group_by(pair_label) %>%
    dplyr::summarise(max_prob = max(prob, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(max_prob) %>%
    dplyr::pull(pair_label)

  df$pair_label <- factor(df$pair_label, levels = pair_order)

  n_pair <- length(unique(df$pair_label))
  fig_height <- max(7.5, 0.23 * n_pair + 3.2)

  pal <- if (color_mode == "red") {
    c("#FDE0DD", "#FA9FB5", "#C0002B")
  } else {
    c("#DEEBF7", "#6BAED6", "#08519C")
  }

  p <- ggplot(
    df,
    aes(x = source_label, y = pair_label, size = prob, fill = prob)
  ) +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.25,
      alpha = 0.95
    ) +
    scale_fill_gradientn(
      colors = pal,
      name = "Communication\nprobability"
    ) +
    scale_size_continuous(
      range = c(1.6, 7.2),
      name = "Communication\nprobability"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = title,
      subtitle = subtitle
    ) +
    theme_final(base_size = 10.5) +
    theme(
      axis.text.x = element_text(
        angle = -45,
        hjust = 0,
        vjust = 1,
        face = "bold",
        size = 9.5
      ),
      axis.text.y = element_text(
        size = ifelse(n_pair > 45, 5.6, ifelse(n_pair > 25, 6.8, 8.2)),
        face = "bold"
      ),
      legend.position = "right"
    )

  save_pdf(
    p,
    outfile,
    width = 10.2,
    height = fig_height
  )

  return(p)
}

parse_receptor_genes <- function(receptor_string) {
  x <- receptor_string
  x <- stringr::str_replace_all(x, "\\(", "")
  x <- stringr::str_replace_all(x, "\\)", "")
  x <- stringr::str_replace_all(x, "\\+", "_")
  x <- unlist(strsplit(x, "_"))
  x <- x[nchar(x) > 0]
  unique(x)
}

########################################################
## 4. regenerate original 01 and 02 with better layout
########################################################
lr_all   <- read_lr(all_file)
lr_dlc1  <- read_lr(dlc1_file)
lr_focus <- read_lr(focus_file)

write.csv(
  lr_all,
  file.path(out_dir, "tables", "01_original_all_LR_pairs_to_C1.csv"),
  row.names = FALSE
)

write.csv(
  lr_focus,
  file.path(out_dir, "tables", "02_original_focus_LR_pairs_to_C1.csv"),
  row.names = FALSE
)

p_orig01 <- plot_lr_bubble_resized(
  lr_all,
  title = "Original CellChat LR screening: all C2/C3/C4-to-C1 interactions",
  subtitle = "Resized view of the full ligand-receptor list",
  outfile = file.path(out_dir, "plots", "FigS1_original_01_all_LR_bubble_resized.pdf"),
  color_mode = "red"
)

p_orig02 <- plot_lr_bubble_resized(
  lr_focus,
  title = "Original CellChat LR screening: migration/adhesion-focused interactions",
  subtitle = "Resized view of focused ligand-receptor candidates",
  outfile = file.path(out_dir, "plots", "FigS2_original_02_focus_LR_bubble_resized.pdf"),
  color_mode = "blue"
)

########################################################
## 5. four-state global interaction plots
########################################################
if (!is.na(cellchat_file) && file.exists(cellchat_file)) {

  cellchat <- readRDS(cellchat_file)

  if (!is.null(cellchat@net$weight)) {

    weight_mat0 <- cellchat@net$weight
    keep <- intersect(state_order, rownames(weight_mat0))
    keep <- keep[keep %in% colnames(weight_mat0)]

    weight_mat <- weight_mat0[keep, keep, drop = FALSE]

    write.csv(
      weight_mat,
      file.path(out_dir, "tables", "03_four_state_interaction_weight_matrix.csv")
    )

    ## strength heatmap
    heat_df <- as.data.frame(as.table(weight_mat))
    colnames(heat_df) <- c("sender", "receiver", "weight")

    heat_df <- heat_df %>%
      dplyr::mutate(
        sender = factor(sender, levels = rev(keep)),
        receiver = factor(receiver, levels = keep),
        sender_label = factor(state_short[as.character(sender)], levels = state_short[rev(keep)]),
        receiver_label = factor(state_short[as.character(receiver)], levels = state_short[keep])
      )

    p_heat_weight <- ggplot(
      heat_df,
      aes(x = receiver_label, y = sender_label, fill = weight)
    ) +
      geom_tile(color = "white", linewidth = 0.7) +
      geom_text(
        aes(label = ifelse(weight > 0, sprintf("%.3f", weight), "")),
        size = 3.4,
        color = "black"
      ) +
      scale_fill_gradientn(
        colors = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
        name = "Interaction\nstrength"
      ) +
      labs(
        x = "Receiver",
        y = "Sender",
        title = "Four-state fibroblast communication network",
        subtitle = "CellChat-inferred interaction strength"
      ) +
      coord_equal() +
      theme_final(base_size = 11) +
      theme(
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right"
      )

    save_pdf(
      p_heat_weight,
      file.path(out_dir, "plots", "Fig0A_four_state_interaction_strength_heatmap.pdf"),
      width = 5.8,
      height = 5.0
    )

    ## circle plot: strength
    if (requireNamespace("CellChat", quietly = TRUE)) {
      pdf(
        file.path(out_dir, "plots", "Fig0B_four_state_interaction_strength_circle.pdf"),
        width = 6.2,
        height = 6.0,
        onefile = FALSE
      )

      group_size <- rep(1, length(keep))
      names(group_size) <- keep

      CellChat::netVisual_circle(
        weight_mat,
        vertex.weight = group_size,
        weight.scale = TRUE,
        label.edge = FALSE,
        color.use = state_cols[keep],
        title.name = "Interaction strength"
      )

      dev.off()
    }
  }

  if (!is.null(cellchat@net$count)) {

    count_mat0 <- cellchat@net$count
    keep2 <- intersect(state_order, rownames(count_mat0))
    keep2 <- keep2[keep2 %in% colnames(count_mat0)]

    count_mat <- count_mat0[keep2, keep2, drop = FALSE]

    write.csv(
      count_mat,
      file.path(out_dir, "tables", "04_four_state_interaction_count_matrix.csv")
    )

    count_df <- as.data.frame(as.table(count_mat))
    colnames(count_df) <- c("sender", "receiver", "count")

    count_df <- count_df %>%
      dplyr::mutate(
        sender = factor(sender, levels = rev(keep2)),
        receiver = factor(receiver, levels = keep2),
        sender_label = factor(state_short[as.character(sender)], levels = state_short[rev(keep2)]),
        receiver_label = factor(state_short[as.character(receiver)], levels = state_short[keep2])
      )

    p_heat_count <- ggplot(
      count_df,
      aes(x = receiver_label, y = sender_label, fill = count)
    ) +
      geom_tile(color = "white", linewidth = 0.7) +
      geom_text(
        aes(label = ifelse(count > 0, count, "")),
        size = 3.4,
        color = "black"
      ) +
      scale_fill_gradientn(
        colors = c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
        name = "Interaction\nnumber"
      ) +
      labs(
        x = "Receiver",
        y = "Sender",
        title = "Four-state fibroblast communication network",
        subtitle = "Number of CellChat-inferred interactions"
      ) +
      coord_equal() +
      theme_final(base_size = 11) +
      theme(
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right"
      )

    save_pdf(
      p_heat_count,
      file.path(out_dir, "plots", "Fig0C_four_state_interaction_count_heatmap.pdf"),
      width = 5.8,
      height = 5.0
    )

    if (requireNamespace("CellChat", quietly = TRUE)) {
      pdf(
        file.path(out_dir, "plots", "Fig0D_four_state_interaction_count_circle.pdf"),
        width = 6.2,
        height = 6.0,
        onefile = FALSE
      )

      group_size <- rep(1, length(keep2))
      names(group_size) <- keep2

      CellChat::netVisual_circle(
        count_mat,
        vertex.weight = group_size,
        weight.scale = TRUE,
        label.edge = FALSE,
        color.use = state_cols[keep2],
        title.name = "Interaction number"
      )

      dev.off()
    }
  }

} else {
  warning("未找到 CellChat RDS，四个亚群互作 circle/heatmap 跳过。")
}

########################################################
## 6. Fig3: change final C3-to-C1 axes into barplot
########################################################
dlc1_screen <- lr_dlc1 %>%
  dplyr::filter(
    ligand %in% c("APP", "MIF", "CXCL12") |
      grepl("APP|MIF|CXCL12", interaction_name, ignore.case = TRUE)
  )

ecm_screen_all <- lr_focus %>%
  dplyr::filter(
    (
      grepl("^COL", ligand) |
        ligand %in% c("SPP1", "FN1") |
        grepl("^LAM", ligand)
    ) &
      has_gene(receptor, "ITG") &
      !has_gene(receptor, "SDC4") &
      !has_gene(receptor, "CD47")
  )

c3_mif <- dlc1_screen %>%
  dplyr::filter(
    source == "C3_Activated_matrix",
    ligand == "MIF",
    has_gene(receptor, "CD74"),
    has_gene(receptor, "CXCR4")
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::mutate(final_class = "Primary MIF axis")

priority_ligands <- c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "SPP1", "FN1", "COL4A1", "COL4A2")

c3_ecm <- ecm_screen_all %>%
  dplyr::filter(source == "C3_Activated_matrix") %>%
  dplyr::mutate(
    ligand_priority = match(ligand, priority_ligands),
    ligand_priority = ifelse(is.na(ligand_priority), 99, ligand_priority)
  ) %>%
  dplyr::arrange(ligand_priority, dplyr::desc(prob)) %>%
  dplyr::distinct(ligand, receptor, .keep_all = TRUE) %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::mutate(final_class = "ECM-integrin context")

final_c3 <- dplyr::bind_rows(c3_mif, c3_ecm) %>%
  dplyr::mutate(
    pair_label_clean = paste0(ligand, " → ", stringr::str_replace_all(receptor, "_", "+")),
    final_class = factor(final_class, levels = c("Primary MIF axis", "ECM-integrin context"))
  ) %>%
  dplyr::arrange(final_class, dplyr::desc(prob))

final_c3$pair_label_clean <- factor(
  final_c3$pair_label_clean,
  levels = final_c3 %>%
    dplyr::arrange(prob) %>%
    dplyr::pull(pair_label_clean)
)

write.csv(
  final_c3,
  file.path(out_dir, "tables", "05_Fig3_final_C3_to_C1_candidate_axes.csv"),
  row.names = FALSE
)

p_fig3_bar <- ggplot(
  final_c3,
  aes(x = prob, y = pair_label_clean, fill = final_class)
) +
  geom_col(
    width = 0.68,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.4f", prob)),
    hjust = -0.12,
    size = 3.1,
    fontface = "bold"
  ) +
  scale_fill_manual(values = axis_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(
    x = "Communication probability",
    y = NULL,
    title = "Final C3→C1 candidate axes",
    subtitle = "Ranked CellChat probability after DLC1-context and ECM-integrin screening"
  ) +
  theme_final(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "bold", size = 9.5),
    legend.position = "top"
  )

save_pdf(
  p_fig3_bar,
  file.path(out_dir, "plots", "Fig3_FINAL_C3_to_C1_candidate_axes_barplot.pdf"),
  width = 8.4,
  height = 4.8
)

########################################################
## 7. Fig5: directional C3 ligand -> C1 receptor expression support
########################################################
fib_obj <- readRDS(obj_file)

traj_order <- c(
  "CXCL12_RBP1_stellate",
  "C7_SFRP4_stromal",
  "Activated_matrix_fibroblast",
  "Perivascular_myofibroblast"
)

traj_labels <- c(
  "C1_CXCL12_RBP1",
  "C2_C7_SFRP4",
  "C3_Activated_matrix",
  "C4_Perivascular_myofib"
)

names(traj_labels) <- traj_order

fib4 <- subset(
  fib_obj,
  cells = rownames(fib_obj@meta.data)[fib_obj$fibro_state %in% traj_order]
)

fib4$fibro_comm_state <- traj_labels[as.character(fib4$fibro_state)]
fib4$fibro_comm_state <- factor(fib4$fibro_comm_state, levels = traj_labels)

DefaultAssay(fib4) <- "RNA"
fib4 <- NormalizeData(fib4, verbose = FALSE)

get_data_safe <- function(obj) {
  tryCatch(
    Seurat::GetAssayData(obj, assay = "RNA", layer = "data"),
    error = function(e) Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
  )
}

expr_summary <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  mat <- get_data_safe(obj)
  mat <- mat[genes, , drop = FALSE]
  group <- as.character(obj@meta.data[colnames(mat), "fibro_comm_state"])

  out <- lapply(genes, function(g) {
    data.frame(
      gene = g,
      group = group,
      expr = as.numeric(mat[g, ]),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::group_by(gene, group) %>%
      dplyr::summarise(
        avg_expr = mean(expr),
        pct_expr = mean(expr > 0) * 100,
        .groups = "drop"
      )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(avg_expr_z = as.numeric(scale(avg_expr))) %>%
    dplyr::ungroup()

  out
}

ecm_pairs_for_fig5 <- c3_ecm %>%
  dplyr::arrange(dplyr::desc(prob)) %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::mutate(
    pair_label_clean = paste0(ligand, " → ", stringr::str_replace_all(receptor, "_", "+"))
  )

genes_for_expr <- unique(c(
  ecm_pairs_for_fig5$ligand,
  unlist(lapply(ecm_pairs_for_fig5$receptor, parse_receptor_genes))
))

expr_df <- expr_summary(fib4, genes_for_expr)

get_expr_metric <- function(gene, group_name, metric = c("avg_expr_z", "pct_expr")) {
  metric <- match.arg(metric)
  tmp <- expr_df %>%
    dplyr::filter(gene == !!gene, group == !!group_name)
  if (nrow(tmp) == 0) return(NA_real_)
  tmp[[metric]][1]
}

directional_list <- list()

for (i in seq_len(nrow(ecm_pairs_for_fig5))) {

  ligand_i <- ecm_pairs_for_fig5$ligand[i]
  receptor_i <- ecm_pairs_for_fig5$receptor[i]
  pair_i <- ecm_pairs_for_fig5$pair_label_clean[i]
  prob_i <- ecm_pairs_for_fig5$prob[i]

  rec_genes <- parse_receptor_genes(receptor_i)
  rec_genes <- intersect(rec_genes, unique(expr_df$gene))

  ligand_avg_z <- get_expr_metric(ligand_i, "C3_Activated_matrix", "avg_expr_z")
  ligand_pct   <- get_expr_metric(ligand_i, "C3_Activated_matrix", "pct_expr")

  receptor_avg_z <- mean(
    sapply(rec_genes, function(g) get_expr_metric(g, "C1_CXCL12_RBP1", "avg_expr_z")),
    na.rm = TRUE
  )

  receptor_pct <- mean(
    sapply(rec_genes, function(g) get_expr_metric(g, "C1_CXCL12_RBP1", "pct_expr")),
    na.rm = TRUE
  )

  directional_list[[i]] <- data.frame(
    pair = pair_i,
    part = c("C3 sender\nligand", "C1 receiver\nintegrin"),
    gene_label = c(ligand_i, paste(rec_genes, collapse = "+")),
    avg_expr_z = c(ligand_avg_z, receptor_avg_z),
    pct_expr = c(ligand_pct, receptor_pct),
    cellchat_prob = prob_i,
    stringsAsFactors = FALSE
  )
}

directional_df <- dplyr::bind_rows(directional_list)

pair_order_fig5 <- ecm_pairs_for_fig5 %>%
  dplyr::arrange(cellchat_prob = prob) %>%
  dplyr::pull(pair_label_clean)

directional_df$pair <- factor(directional_df$pair, levels = pair_order_fig5)
directional_df$part <- factor(directional_df$part, levels = c("C3 sender\nligand", "C1 receiver\nintegrin"))

write.csv(
  directional_df,
  file.path(out_dir, "tables", "06_Fig5_directional_C3_ligand_C1_integrin_expression.csv"),
  row.names = FALSE
)

p_fig5_directional <- ggplot(
  directional_df,
  aes(x = part, y = pair, size = pct_expr, fill = avg_expr_z)
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.30,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = gene_label),
    nudge_x = 0.22,
    size = 2.9,
    hjust = 0,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    name = "Scaled\naverage\nexpression"
  ) +
  scale_size_continuous(
    range = c(2.0, 8.0),
    limits = c(0, 100),
    name = "Percent\nexpressed"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = NULL,
    title = "Directional expression support for C3→C1 ECM–integrin interactions",
    subtitle = "Sender ligands are evaluated in C3; receiver integrins are evaluated in C1"
  ) +
  theme_final(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "right",
    plot.margin = margin(8, 80, 8, 8)
  )

save_pdf(
  p_fig5_directional,
  file.path(out_dir, "plots", "Fig5_FINAL_directional_C3_ligand_C1_integrin_expression.pdf"),
  width = 9.2,
  height = 5.2
)

########################################################
## 8. README
########################################################
writeLines(
  c(
    "FINAL_FIX_v5 outputs.",
    "",
    "Main corrected outputs:",
    "Fig0A_four_state_interaction_strength_heatmap.pdf",
    "Fig0B_four_state_interaction_strength_circle.pdf",
    "Fig0C_four_state_interaction_count_heatmap.pdf",
    "Fig0D_four_state_interaction_count_circle.pdf",
    "FigS1_original_01_all_LR_bubble_resized.pdf",
    "FigS2_original_02_focus_LR_bubble_resized.pdf",
    "Fig3_FINAL_C3_to_C1_candidate_axes_barplot.pdf",
    "Fig5_FINAL_directional_C3_ligand_C1_integrin_expression.pdf",
    "",
    "Important interpretation:",
    "Fig1/Fig3 show CellChat communication probability in the C3 sender to C1 receiver direction.",
    "Fig5 now validates the same direction by showing C3 ligand expression and C1 receptor expression, instead of a non-directional all-state DotPlot."
  ),
  con = file.path(out_dir, "README_FINAL_FIX_v5.txt")
)

cat("FINAL_FIX_v5 finished.\n")
cat("Output dir:\n")
cat(out_dir, "\n")








########################################################
## 原始 CellChat 风格的 01 / 02 bubble plot
########################################################

########################################################
## Regenerate FigS1 / FigS2 using original CellChat style
## Only modify height, y-axis spacing, and x-axis labels
########################################################

suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

########################################################
## 1. paths
########################################################
base_dir  <- file.path(hcc_dlc1_root(), "single_cell")
atlas_dir <- file.path(base_dir, "DLC1_sc_149614_restart_by_original_style")
fib_dir   <- file.path(atlas_dir, "fibro_subcluster", "fibro_cleaned_round2", "fibro_state_final_v3")

cc_dir <- file.path(fib_dir, "CellChat_C2C3C4_to_C1_DLC1_transition")

all_file <- file.path(
  cc_dir, "tables", "04_MAIN_C2C3C4_to_C1_incoming_LR_pairs.csv"
)

focus_file <- file.path(
  cc_dir, "tables", "07_MAIN_C2C3C4_to_C1_focus_migration_adhesion_LR_pairs.csv"
)

cellchat_rds_candidates <- c(
  file.path(cc_dir, "rds", "02_CellChat_4states_full_object.rds"),
  file.path(cc_dir, "rds", "02_CellChat_C2C3C4_to_C1_DLC1_transition.rds"),
  file.path(cc_dir, "rds", "03_CellChat_fibro_internal_4states_with_centrality.rds"),
  file.path(cc_dir, "CellChat_4states_full_object.rds")
)

cellchat_file <- cellchat_rds_candidates[file.exists(cellchat_rds_candidates)][1]

out_dir <- file.path(cc_dir, "FINAL_Figure_CellChat_FIX_v5", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (is.na(cellchat_file) || !file.exists(cellchat_file)) {
  stop("没有找到 CellChat RDS 文件。")
}

if (!file.exists(all_file)) stop("找不到 all LR 文件: ", all_file)
if (!file.exists(focus_file)) stop("找不到 focus LR 文件: ", focus_file)

cellchat <- readRDS(cellchat_file)

########################################################
## 2. source / target
########################################################
sources_use <- c(
  "C2_C7_SFRP4",
  "C3_Activated_matrix",
  "C4_Perivascular_myofib"
)

target_use <- "C1_CXCL12_RBP1"

x_label_fun <- function(x) {
  x <- stringr::str_replace_all(
    x,
    "C2_C7_SFRP4 -> C1_CXCL12_RBP1",
    "C2_C7_SFRP4\n→ C1_CXCL12_RBP1"
  )
  x <- stringr::str_replace_all(
    x,
    "C3_Activated_matrix -> C1_CXCL12_RBP1",
    "C3_Activated_matrix\n→ C1_CXCL12_RBP1"
  )
  x <- stringr::str_replace_all(
    x,
    "C4_Perivascular_myofib -> C1_CXCL12_RBP1",
    "C4_Perivascular_myofib\n→ C1_CXCL12_RBP1"
  )
  x <- stringr::str_replace_all(x, " - ", " → ")
  x
}

########################################################
## 3. helper function
########################################################
make_cellchat_bubble <- function(
    lr_table_file,
    outfile,
    title_text,
    width = 11.5
) {

  lr_df <- read.csv(
    lr_table_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"interaction_name" %in% colnames(lr_df)) {
    stop("LR table 里没有 interaction_name 列，无法用于 CellChat::netVisual_bubble。")
  }

  pair_use <- data.frame(
    interaction_name = unique(lr_df$interaction_name),
    stringsAsFactors = FALSE
  )

  n_pair <- length(unique(lr_df$interaction_name))

  ## 关键：根据受配体数量动态拉高 PDF
  fig_height <- max(9, 0.22 * n_pair + 4)

  p <- CellChat::netVisual_bubble(
    object = cellchat,
    sources.use = sources_use,
    targets.use = target_use,
    pairLR.use = pair_use,
    remove.isolate = FALSE,
    thresh = 0.05,
    color.heatmap = "Spectral",
    angle.x = 45
  )

  ## 只改排版，不改 CellChat 原始 visual style
  p <- p +
    scale_x_discrete(labels = x_label_fun) +
    labs(title = title_text) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14,
        color = "black"
      ),
      axis.text.x = element_text(
        angle = -45,
        hjust = 0,
        vjust = 1,
        size = 9,
        color = "black",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = ifelse(n_pair > 55, 5.8, ifelse(n_pair > 35, 6.8, 8)),
        color = "black"
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(8, 24, 8, 8)
    )

  ggsave(
    filename = outfile,
    plot = p,
    width = width,
    height = fig_height,
    device = cairo_pdf
  )

  message("Saved: ", outfile)
  message("n_pair = ", n_pair, "; height = ", round(fig_height, 2))
}

########################################################
## 4. regenerate FigS1 / FigS2
########################################################
make_cellchat_bubble(
  lr_table_file = all_file,
  outfile = file.path(
    out_dir,
    "FigS1_original_01_all_LR_bubble_CellChat_style_resized.pdf"
  ),
  title_text = "Original CellChat LR screening: all C2/C3/C4-to-C1 interactions",
  width = 11.8
)

make_cellchat_bubble(
  lr_table_file = focus_file,
  outfile = file.path(
    out_dir,
    "FigS2_original_02_focus_LR_bubble_CellChat_style_resized.pdf"
  ),
  title_text = "Original CellChat LR screening: migration/adhesion-focused interactions",
  width = 11.8
)






########################################################
## 空间转录组
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(stringr)
  library(scales)
  library(grid)
})

########################################################
## 1. 路径设置
########################################################
st_root <- file.path(hcc_dlc1_root(), "spatial_transcriptomics", "GSE238264")

out_dir <- file.path(st_root, "DLC1_spatial_HE_native_reanalysis_all7_noDetected_rightLegend")
plot_dir <- file.path(out_dir, "plots")
single_dir <- file.path(plot_dir, "per_sample")
table_dir <- file.path(out_dir, "tables")
rds_dir <- file.path(out_dir, "rds")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(single_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

########################################################
## 2. 自动寻找 7 个 Visium 样本
########################################################
h5_files <- list.files(
  st_root,
  pattern = "filtered_feature_bc_matrix.h5$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(h5_files) == 0) {
  stop("没有找到 filtered_feature_bc_matrix.h5，请检查 st_root 路径。")
}

sample_info <- data.frame(
  h5_file = h5_files,
  data_dir = dirname(h5_files),
  stringsAsFactors = FALSE
)

sample_info$sample_id <- basename(sample_info$data_dir)

sample_info$sample_id <- ifelse(
  sample_info$sample_id == "outs",
  basename(dirname(sample_info$data_dir)),
  sample_info$sample_id
)

sample_info <- sample_info %>%
  dplyr::arrange(sample_id)

write.csv(
  sample_info,
  file.path(table_dir, "00_detected_visium_samples.csv"),
  row.names = FALSE
)

cat("Detected Visium samples:\n")
print(sample_info)

########################################################
## 3. marker gene sets
########################################################
C1_genes <- c(
  "CXCL12", "RBP1", "IGFBP3", "CYGB", "FABP4", "DCN"
)

C3_genes <- c(
  "POSTN", "CTHRC1", "COL1A1", "COL1A2",
  "COL3A1", "COL6A1", "COL6A2", "FN1", "THBS2"
)

MIF_ligand_genes <- c("MIF")

MIF_receptor_genes <- c("CD74", "CXCR4")

ECM_ligand_genes <- c(
  "COL1A1", "COL1A2", "COL6A1", "COL6A2", "FN1", "SPP1"
)

Integrin_receptor_genes <- c(
  "ITGA1", "ITGA9", "ITGB1"
)

Adhesion_genes <- c(
  "ITGA1", "ITGA9", "ITGB1", "FN1",
  "COL1A1", "COL1A2", "COL6A1", "COL6A2"
)

########################################################
## 4. 工具函数
########################################################
save_pdf <- function(p, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = p,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

get_data_safe <- function(obj, assay = "Spatial") {
  mat <- tryCatch(
    suppressWarnings(GetAssayData(obj, assay = assay, layer = "data")),
    error = function(e) {
      suppressWarnings(GetAssayData(obj, assay = assay, slot = "data"))
    }
  )
  return(mat)
}

add_mean_expr_score <- function(obj, genes, score_name, assay = "Spatial") {

  mat <- get_data_safe(obj, assay = assay)
  genes_use <- intersect(genes, rownames(mat))

  if (length(genes_use) < 1) {
    obj[[score_name]] <- NA_real_
    obj[[paste0(score_name, "_z")]] <- NA_real_
    warning(score_name, " 没有可用基因。")
    return(obj)
  }

  score <- Matrix::colMeans(mat[genes_use, , drop = FALSE])
  score <- as.numeric(score)
  names(score) <- colnames(mat)

  obj[[score_name]] <- score[colnames(obj)]

  if (sd(score, na.rm = TRUE) > 0) {
    score_z <- as.numeric(scale(score))
    names(score_z) <- names(score)
    obj[[paste0(score_name, "_z")]] <- score_z[colnames(obj)]
  } else {
    obj[[paste0(score_name, "_z")]] <- 0
  }

  return(obj)
}

safe_rowmean <- function(df, cols) {
  cols <- intersect(cols, colnames(df))
  if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
  rowMeans(df[, cols, drop = FALSE], na.rm = TRUE)
}

add_context_scores <- function(obj) {

  meta <- obj@meta.data

  obj$MIF_axis_context_score <- safe_rowmean(
    meta,
    c("MIF_ligand_score_z", "MIF_receptor_score_z")
  )

  obj$ECM_integrin_context_score <- safe_rowmean(
    meta,
    c("ECM_ligand_score_z", "Integrin_receptor_score_z", "Adhesion_ECM_score_z")
  )

  return(obj)
}

legend_number_fmt <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(x)
  mx <- max(abs(x), na.rm = TRUE)
  acc <- ifelse(mx < 1, 0.01, 0.1)
  scales::number(x, accuracy = acc, trim = TRUE)
}

theme_sp_native <- function(base_size = 10) {
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = base_size + 1,
      color = "black"
    ),

    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box = "vertical",

    legend.title = element_text(
      face = "bold",
      size = base_size - 2,
      color = "black",
      hjust = 0
    ),
    legend.text = element_text(
      size = base_size - 3,
      color = "black"
    ),

    legend.key.height = unit(0.36, "cm"),
    legend.key.width  = unit(0.18, "cm"),
    legend.margin = margin(0, 0, 0, 2),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.03, "cm"),

    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),

    plot.margin = margin(2, 2, 2, 2)
  )
}

########################################################
## Seurat 原生 H&E overlay 画图函数
## 重点修改：
## 1. combine = FALSE，防止 Seurat 默认组合后 legend 位置不好控制
## 2. legend.position = right
## 3. breaks_pretty(n = 3)，减少色标数字重叠
########################################################
plot_spatial_native <- function(
    obj,
    feature,
    title,
    cols_use,
    pt_size = 1.65,
    image_alpha = 0.95,
    min_cutoff = "q02",
    max_cutoff = "q98",
    crop_use = TRUE,
    show_legend = TRUE
) {

  if (!feature %in% c(rownames(obj), colnames(obj@meta.data))) {
    warning("Feature not found: ", feature)
    obj[[feature]] <- NA_real_
  }

  p_list <- SpatialFeaturePlot(
    object = obj,
    features = feature,
    image.alpha = image_alpha,
    pt.size.factor = pt_size,
    alpha = c(0.03, 1),
    min.cutoff = min_cutoff,
    max.cutoff = max_cutoff,
    crop = crop_use,
    combine = FALSE
  )

  p <- p_list[[1]] +
    scale_fill_gradientn(
      colors = cols_use,
      na.value = "grey90",
      name = title,
      breaks = scales::breaks_pretty(n = 3),
      labels = legend_number_fmt
    ) +
    ggtitle(title) +
    theme_sp_native(base_size = 10) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0,
        label.position = "right",
        barheight = unit(1.55, "cm"),
        barwidth = unit(0.18, "cm"),
        ticks = TRUE
      )
    ) +
    theme(
      legend.position = "right"
    )

  if (!show_legend) {
    p <- p + NoLegend()
  }

  return(p)
}

plot_he_native <- function(obj, sid) {

  obj$HE_signal <- 1

  p_list <- SpatialFeaturePlot(
    object = obj,
    features = "HE_signal",
    image.alpha = 1,
    pt.size.factor = 0.35,
    alpha = c(0.02, 0.08),
    crop = TRUE,
    combine = FALSE
  )

  p <- p_list[[1]] +
    scale_fill_gradientn(
      colors = c("grey90", "grey90"),
      na.value = "grey90"
    ) +
    ggtitle(paste0(sid, " H&E")) +
    theme_sp_native(base_size = 10) +
    NoLegend()

  return(p)
}

########################################################
## 5. palettes
## 已删除 pal_detect
########################################################
pal_dlc1 <- c("#F7F7F7", "#FDD0A2", "#FB6A4A", "#A50F15")
pal_c1 <- c("#F7FCF5", "#C7E9C0", "#41AB5D", "#005A32")
pal_c3 <- c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#A50F15")
pal_mif_lig <- c("#F7F7F7", "#FDD0A2", "#FC9272", "#CB181D")
pal_mif_rec <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#08306B")
pal_ecm <- c("#F7F4F9", "#D4B9DA", "#8856A7", "#3F007D")

########################################################
## 6. loop all samples
########################################################
all_core4_list <- list()
all_full7_list <- list()
all_spot_tables <- list()
all_gene_check <- list()
all_summary <- list()
all_feature_check <- list()

for (i in seq_len(nrow(sample_info))) {

  sid <- sample_info$sample_id[i]
  data_dir_i <- sample_info$data_dir[i]
  h5_i <- basename(sample_info$h5_file[i])

  message("\n==============================")
  message("Processing sample: ", sid)
  message("data_dir: ", data_dir_i)
  message("==============================")

  obj <- Load10X_Spatial(
    data.dir = data_dir_i,
    filename = h5_i,
    assay = "Spatial",
    slice = sid
  )

  obj$sample_id <- sid

  DefaultAssay(obj) <- "Spatial"
  obj <- NormalizeData(obj, verbose = FALSE)

  data_mat <- get_data_safe(obj, assay = "Spatial")

  ######################################################
  ## 6.1 DLC1 expression
  ## 已删除 DLC1_detected_num
  ######################################################
  if ("DLC1" %in% rownames(data_mat)) {
    dlc1_vec <- as.numeric(data_mat["DLC1", colnames(obj), drop = TRUE])
    names(dlc1_vec) <- colnames(obj)
    obj$DLC1_expr <- dlc1_vec[colnames(obj)]
  } else {
    obj$DLC1_expr <- NA_real_
    warning("样本 ", sid, " 中没有 DLC1。")
  }

  ######################################################
  ## 6.2 module scores
  ######################################################
  obj <- add_mean_expr_score(obj, C1_genes, "C1_score")
  obj <- add_mean_expr_score(obj, C3_genes, "C3_score")
  obj <- add_mean_expr_score(obj, MIF_ligand_genes, "MIF_ligand_score")
  obj <- add_mean_expr_score(obj, MIF_receptor_genes, "MIF_receptor_score")
  obj <- add_mean_expr_score(obj, ECM_ligand_genes, "ECM_ligand_score")
  obj <- add_mean_expr_score(obj, Integrin_receptor_genes, "Integrin_receptor_score")
  obj <- add_mean_expr_score(obj, Adhesion_genes, "Adhesion_ECM_score")
  obj <- add_context_scores(obj)

  ######################################################
  ## 6.3 gene presence check
  ######################################################
  genes_check <- unique(c(
    "DLC1",
    C1_genes,
    C3_genes,
    MIF_ligand_genes,
    MIF_receptor_genes,
    ECM_ligand_genes,
    Integrin_receptor_genes,
    Adhesion_genes
  ))

  gene_check_i <- data.frame(
    sample_id = sid,
    gene = genes_check,
    present = genes_check %in% rownames(obj),
    stringsAsFactors = FALSE
  )

  all_gene_check[[sid]] <- gene_check_i

  ######################################################
  ## 6.4 feature value check
  ## 已删除 DLC1_detected_num
  ######################################################
  features_check <- c(
    "DLC1_expr",
    "C1_score",
    "C3_score",
    "MIF_ligand_score",
    "MIF_receptor_score",
    "ECM_ligand_score",
    "Integrin_receptor_score",
    "Adhesion_ECM_score",
    "MIF_axis_context_score",
    "ECM_integrin_context_score"
  )

  feature_check_i <- lapply(features_check, function(f) {
    x <- obj@meta.data[[f]]
    data.frame(
      sample_id = sid,
      feature = f,
      n_non_na = sum(!is.na(x)),
      n_positive = sum(x > 0, na.rm = TRUE),
      min_value = suppressWarnings(min(x, na.rm = TRUE)),
      median_value = suppressWarnings(median(x, na.rm = TRUE)),
      max_value = suppressWarnings(max(x, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }) %>%
    dplyr::bind_rows()

  all_feature_check[[sid]] <- feature_check_i

  ######################################################
  ## 6.5 spot score table
  ## 已删除 DLC1_detected_num
  ######################################################
  meta_i <- obj@meta.data

  spot_table_i <- meta_i %>%
    dplyr::mutate(
      sample_id = sid,
      spot = rownames(meta_i)
    ) %>%
    dplyr::select(
      sample_id,
      spot,
      nCount_Spatial,
      nFeature_Spatial,
      DLC1_expr,
      C1_score,
      C1_score_z,
      C3_score,
      C3_score_z,
      MIF_ligand_score,
      MIF_ligand_score_z,
      MIF_receptor_score,
      MIF_receptor_score_z,
      ECM_ligand_score,
      ECM_ligand_score_z,
      Integrin_receptor_score,
      Integrin_receptor_score_z,
      Adhesion_ECM_score,
      Adhesion_ECM_score_z,
      MIF_axis_context_score,
      ECM_integrin_context_score
    )

  all_spot_tables[[sid]] <- spot_table_i

  summary_i <- spot_table_i %>%
    dplyr::summarise(
      sample_id = sid,
      n_spots = dplyr::n(),
      DLC1_mean = mean(DLC1_expr, na.rm = TRUE),
      DLC1_median = median(DLC1_expr, na.rm = TRUE),
      C1_score_mean = mean(C1_score, na.rm = TRUE),
      C3_score_mean = mean(C3_score, na.rm = TRUE),
      MIF_ligand_mean = mean(MIF_ligand_score, na.rm = TRUE),
      MIF_receptor_mean = mean(MIF_receptor_score, na.rm = TRUE),
      ECM_integrin_context_mean = mean(ECM_integrin_context_score, na.rm = TRUE)
    )

  all_summary[[sid]] <- summary_i

  ######################################################
  ## 6.6 Native HE plots
  ######################################################
  p_he <- plot_he_native(obj, sid)

  p_dlc1 <- plot_spatial_native(
    obj,
    feature = "DLC1_expr",
    title = "DLC1",
    cols_use = pal_dlc1,
    pt_size = 1.65,
    image_alpha = 0.95,
    min_cutoff = "q00",
    max_cutoff = "q98"
  )

  p_c1 <- plot_spatial_native(
    obj,
    feature = "C1_score",
    title = "C1-like score",
    cols_use = pal_c1,
    pt_size = 1.65,
    image_alpha = 0.95,
    min_cutoff = "q02",
    max_cutoff = "q98"
  )

  p_c3 <- plot_spatial_native(
    obj,
    feature = "C3_score",
    title = "C3-like score",
    cols_use = pal_c3,
    pt_size = 1.65,
    image_alpha = 0.95,
    min_cutoff = "q02",
    max_cutoff = "q98"
  )

  p_mif_lig <- plot_spatial_native(
    obj,
    feature = "MIF_ligand_score",
    title = "MIF ligand",
    cols_use = pal_mif_lig,
    pt_size = 1.65,
    image_alpha = 0.95,
    min_cutoff = "q02",
    max_cutoff = "q98"
  )

  p_mif_rec <- plot_spatial_native(
    obj,
    feature = "MIF_receptor_score",
    title = "CD74/CXCR4",
    cols_use = pal_mif_rec,
    pt_size = 1.65,
    image_alpha = 0.95,
    min_cutoff = "q02",
    max_cutoff = "q98"
  )

  p_ecm_context <- plot_spatial_native(
    obj,
    feature = "ECM_integrin_context_score",
    title = "ECM-integrin context",
    cols_use = pal_ecm,
    pt_size = 1.65,
    image_alpha = 0.95,
    min_cutoff = "q02",
    max_cutoff = "q98"
  )

  ######################################################
  ## 6.7 Save individual panels
  ######################################################
  save_pdf(p_he, file.path(single_dir, paste0(sid, "_00_HE_native.pdf")), 6.4, 5.6)
  save_pdf(p_dlc1, file.path(single_dir, paste0(sid, "_01_DLC1_native.pdf")), 6.4, 5.6)
  save_pdf(p_c1, file.path(single_dir, paste0(sid, "_02_C1_score_native.pdf")), 6.4, 5.6)
  save_pdf(p_c3, file.path(single_dir, paste0(sid, "_03_C3_score_native.pdf")), 6.4, 5.6)
  save_pdf(p_mif_lig, file.path(single_dir, paste0(sid, "_04_MIF_ligand_native.pdf")), 6.4, 5.6)
  save_pdf(p_mif_rec, file.path(single_dir, paste0(sid, "_05_MIF_receptor_native.pdf")), 6.4, 5.6)
  save_pdf(p_ecm_context, file.path(single_dir, paste0(sid, "_06_ECM_integrin_context_native.pdf")), 6.4, 5.6)

  ######################################################
  ## 6.8 Core4: HE / DLC1 / C1 / C3
  ######################################################
  p_core4 <- (
    p_he | p_dlc1 | p_c1 | p_c3
  ) +
    plot_annotation(
      title = paste0(sid, " DLC1 / C1 / C3 spatial core"),
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 16,
          hjust = 0.5,
          color = "black"
        )
      )
    )

  save_pdf(
    p_core4,
    file.path(single_dir, paste0(sid, "_core4_native_HE_DLC1_C1_C3.pdf")),
    width = 19.2,
    height = 4.8
  )

  ######################################################
  ## 6.9 Full7
  ## 已删除 DLC1 detected
  ######################################################
  p_full7 <- (
    p_he | p_dlc1 | p_c1 | p_c3
  ) / (
    p_mif_lig | p_mif_rec | p_ecm_context | plot_spacer()
  ) +
    plot_annotation(
      title = paste0(sid, " spatial distribution on H&E"),
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 16,
          hjust = 0.5,
          color = "black"
        )
      )
    )

  save_pdf(
    p_full7,
    file.path(single_dir, paste0(sid, "_full7_native_HE_publication.pdf")),
    width = 19.2,
    height = 8.8
  )

  all_core4_list[[sid]] <- p_core4
  all_full7_list[[sid]] <- p_full7

  saveRDS(
    obj,
    file.path(rds_dir, paste0(sid, "_native_HE_scores.rds"))
  )
}

########################################################
## 7. 合并七个样本 PDF
########################################################
all7_core4_pdf <- file.path(plot_dir, "ALL7_core4_native_HE_DLC1_C1_C3.pdf")

cairo_pdf(
  filename = all7_core4_pdf,
  width = 19.2,
  height = 4.8,
  onefile = TRUE
)

for (sid in names(all_core4_list)) {
  print(all_core4_list[[sid]])
}

dev.off()

all7_full7_pdf <- file.path(plot_dir, "ALL7_full7_native_HE_publication.pdf")

cairo_pdf(
  filename = all7_full7_pdf,
  width = 19.2,
  height = 8.8,
  onefile = TRUE
)

for (sid in names(all_full7_list)) {
  print(all_full7_list[[sid]])
}

dev.off()

########################################################
## 8. 输出表格
########################################################
write.csv(
  dplyr::bind_rows(all_spot_tables),
  file.path(table_dir, "01_all7_spot_scores_native_HE.csv"),
  row.names = FALSE
)

write.csv(
  dplyr::bind_rows(all_gene_check),
  file.path(table_dir, "02_all7_gene_presence_native_HE.csv"),
  row.names = FALSE
)

write.csv(
  dplyr::bind_rows(all_summary),
  file.path(table_dir, "03_all7_sample_summary_native_HE.csv"),
  row.names = FALSE
)

write.csv(
  dplyr::bind_rows(all_feature_check),
  file.path(table_dir, "04_all7_feature_value_check_native_HE.csv"),
  row.names = FALSE
)

########################################################
## 9. README
########################################################
writeLines(
  c(
    "Native Seurat H&E overlay reanalysis for all seven samples.",
    "",
    "Modified version:",
    "1. DLC1 detected panel has been removed from all plots.",
    "2. Colorbars are placed on the right side of each panel.",
    "3. Colorbar tick labels are simplified to reduce overlap.",
    "",
    "Important:",
    "This script does not depend on old RDS files.",
    "It reloads raw Visium data from filtered_feature_bc_matrix.h5 and spatial folders.",
    "It uses Seurat native SpatialFeaturePlot H&E overlay style.",
    "",
    "Main outputs:",
    "plots/ALL7_core4_native_HE_DLC1_C1_C3.pdf",
    "plots/ALL7_full7_native_HE_publication.pdf",
    "plots/per_sample/*_native.pdf",
    "",
    "Core4 layout:",
    "1. H&E",
    "2. DLC1 expression",
    "3. C1-like score",
    "4. C3-like score",
    "",
    "Full7 layout:",
    "1. H&E",
    "2. DLC1 expression",
    "3. C1-like score",
    "4. C3-like score",
    "5. MIF ligand",
    "6. CD74/CXCR4 receptor score",
    "7. ECM-integrin context",
    "",
    "Tables:",
    "01_all7_spot_scores_native_HE.csv",
    "02_all7_gene_presence_native_HE.csv",
    "03_all7_sample_summary_native_HE.csv",
    "04_all7_feature_value_check_native_HE.csv"
  ),
  con = file.path(out_dir, "README_native_HE_noDetected_rightLegend.txt")
)

message("\n全部完成！")
message("输出目录：", out_dir)
message("Core4 PDF：", all7_core4_pdf)
message("Full7 PDF：", all7_full7_pdf)




########################################################
## 空间邻近统计分析完整脚本
########################################################
########################################################
## DLC1 空间转录组后续统计分析
## C1-high spots 是否靠近 C3-high spots
##
## 核心逻辑：
## 在 C1-high spots 中，根据是否靠近 C3-high spots 分组：
## 1. C1_near_C3_high
## 2. C1_not_near_C3_high
##
## 比较：
## DLC1 expression
## MIF receptor score
## Integrin receptor score
## ECM-integrin context score
## Adhesion/ECM score
########################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
set.seed(123)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(scales)
  library(grid)
})

########################################################
## 1. 路径设置
########################################################

st_root <- file.path(hcc_dlc1_root(), "spatial_transcriptomics", "GSE238264")

## 优先读取你最新修改后的 noDetected/rightLegend 版本
candidate_rds_dirs <- c(
  file.path(st_root, "DLC1_spatial_HE_native_reanalysis_all7_noDetected_rightLegend", "rds"),
  file.path(st_root, "DLC1_spatial_HE_native_reanalysis_all7", "rds")
)

candidate_rds_dirs <- candidate_rds_dirs[dir.exists(candidate_rds_dirs)]

if (length(candidate_rds_dirs) == 0) {
  stop("没有找到前一步生成的 rds 文件夹，请先运行 native H&E overlay 脚本。")
}

rds_dir <- candidate_rds_dirs[1]

out_dir <- file.path(st_root, "DLC1_spatial_C1_C3_neighbor_statistics")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
rds_out_dir <- file.path(out_dir, "rds_annotated")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_out_dir, showWarnings = FALSE, recursive = TRUE)

message("Using RDS directory: ", rds_dir)
message("Output directory: ", out_dir)

rds_files <- list.files(
  rds_dir,
  pattern = "_native_HE_scores\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("rds_dir 中没有找到 *_native_HE_scores.rds 文件。")
}

rds_files <- sort(rds_files)

sample_ids <- sub("_native_HE_scores\\.rds$", "", basename(rds_files))

sample_info <- data.frame(
  sample_id = sample_ids,
  rds_file = rds_files,
  stringsAsFactors = FALSE
)

write.csv(
  sample_info,
  file.path(table_dir, "00_input_rds_files.csv"),
  row.names = FALSE
)

print(sample_info)

########################################################
## 2. 参数设置
########################################################

## C1-high / C3-high 的定义：每个样本内部 top 25%
c1_high_quantile <- 0.75
c3_high_quantile <- 0.75

## 空间邻近半径
## 1.5 表示约一阶邻近；
## 2.5 可以理解为更宽松的邻近；
## 这里 primary 用 1.5，同时保留 sensitivity 分析。
primary_radius_multiplier <- 1.5
sensitivity_radius_multipliers <- c(1.5, 2.0, 2.5, 3.0)

## 至少多少个 C3-high spot 出现在邻近范围内，才认为 near
min_C3_neighbors <- 1

## 要比较的核心特征
compare_features <- c(
  "DLC1_expr",
  "C1_score",
  "C3_score",
  "MIF_ligand_score",
  "MIF_receptor_score",
  "Integrin_receptor_score",
  "ECM_integrin_context_score",
  "Adhesion_ECM_score"
)

feature_labels <- c(
  DLC1_expr = "DLC1 expression",
  C1_score = "C1-like score",
  C3_score = "C3-like score",
  MIF_ligand_score = "MIF ligand score",
  MIF_receptor_score = "CD74/CXCR4 receptor score",
  Integrin_receptor_score = "Integrin receptor score",
  ECM_integrin_context_score = "ECM-integrin context score",
  Adhesion_ECM_score = "Adhesion/ECM score"
)

########################################################
## 3. 工具函数
########################################################

save_pdf <- function(p, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = p,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
}

get_assay_data_safe <- function(obj, assay = "Spatial", layer = "data") {
  mat <- tryCatch(
    suppressWarnings(GetAssayData(obj, assay = assay, layer = layer)),
    error = function(e) {
      tryCatch(
        suppressWarnings(GetAssayData(obj, assay = assay, slot = layer)),
        error = function(e2) NULL
      )
    }
  )
  return(mat)
}

add_dlc1_raw_detected <- function(obj, assay = "Spatial") {

  counts_mat <- get_assay_data_safe(obj, assay = assay, layer = "counts")

  if (!is.null(counts_mat) && "DLC1" %in% rownames(counts_mat)) {
    dlc1_raw <- as.numeric(counts_mat["DLC1", colnames(obj), drop = TRUE])
    names(dlc1_raw) <- colnames(obj)
    obj$DLC1_raw_count <- dlc1_raw[colnames(obj)]
    obj$DLC1_raw_detected <- ifelse(obj$DLC1_raw_count > 0, 1, 0)
  } else {
    obj$DLC1_raw_count <- NA_real_
    obj$DLC1_raw_detected <- NA_real_
  }

  return(obj)
}

extract_spatial_coords <- function(obj) {

  image_name <- names(obj@images)[1]

  coords <- tryCatch(
    as.data.frame(GetTissueCoordinates(obj, image = image_name)),
    error = function(e) NULL
  )

  if (is.null(coords) || nrow(coords) == 0) {
    coords <- tryCatch(
      as.data.frame(obj@images[[image_name]]@coordinates),
      error = function(e) NULL
    )
  }

  if (is.null(coords) || nrow(coords) == 0) {
    stop("无法提取空间坐标。")
  }

  if ("cell" %in% colnames(coords)) {
    coords$spot <- as.character(coords$cell)
  } else if ("barcode" %in% colnames(coords)) {
    coords$spot <- as.character(coords$barcode)
  } else {
    coords$spot <- rownames(coords)
  }

  if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
    coords$x <- coords$imagecol
    coords$y <- coords$imagerow
  } else if (all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% colnames(coords))) {
    coords$x <- coords$pxl_col_in_fullres
    coords$y <- coords$pxl_row_in_fullres
  } else if (all(c("x", "y") %in% colnames(coords))) {
    coords$x <- coords$x
    coords$y <- coords$y
  } else if (all(c("col", "row") %in% colnames(coords))) {
    coords$x <- coords$col
    coords$y <- coords$row
  } else if (all(c("array_col", "array_row") %in% colnames(coords))) {
    coords$x <- coords$array_col
    coords$y <- coords$array_row
  } else {
    stop("坐标文件中没有识别到可用的 x/y 坐标列。")
  }

  coords <- coords %>%
    dplyr::select(spot, x, y) %>%
    dplyr::distinct(spot, .keep_all = TRUE)

  coords <- coords[coords$spot %in% colnames(obj), , drop = FALSE]

  if (nrow(coords) == 0) {
    stop("空间坐标 spot 名称和 Seurat object colnames 无法匹配。")
  }

  return(coords)
}

estimate_nn_distance <- function(coords_df, max_sample = 2500) {

  xy <- as.matrix(coords_df[, c("x", "y")])
  n <- nrow(xy)

  if (n < 2) return(NA_real_)

  if (n > max_sample) {
    set.seed(123)
    idx <- sample(seq_len(n), max_sample)
    xy_use <- xy[idx, , drop = FALSE]
  } else {
    xy_use <- xy
  }

  d <- as.matrix(dist(xy_use))
  diag(d) <- NA_real_

  nn <- apply(d, 1, min, na.rm = TRUE)
  nn <- nn[is.finite(nn)]

  if (length(nn) == 0) return(NA_real_)

  median(nn, na.rm = TRUE)
}

nearest_distance_to_targets <- function(query_df, target_df, exclude_self = FALSE) {

  if (nrow(query_df) == 0 || nrow(target_df) == 0) {
    return(rep(NA_real_, nrow(query_df)))
  }

  q_xy <- as.matrix(query_df[, c("x", "y")])
  t_xy <- as.matrix(target_df[, c("x", "y")])
  t_spot <- target_df$spot

  out <- numeric(nrow(query_df))

  for (i in seq_len(nrow(query_df))) {
    d <- sqrt((t_xy[, 1] - q_xy[i, 1])^2 + (t_xy[, 2] - q_xy[i, 2])^2)

    if (exclude_self) {
      d[t_spot == query_df$spot[i]] <- Inf
    }

    d_min <- suppressWarnings(min(d, na.rm = TRUE))

    if (!is.finite(d_min)) {
      out[i] <- NA_real_
    } else {
      out[i] <- d_min
    }
  }

  out
}

count_targets_within_radius <- function(query_df, target_df, radius, exclude_self = FALSE) {

  if (nrow(query_df) == 0 || nrow(target_df) == 0 || is.na(radius)) {
    return(rep(0, nrow(query_df)))
  }

  q_xy <- as.matrix(query_df[, c("x", "y")])
  t_xy <- as.matrix(target_df[, c("x", "y")])
  t_spot <- target_df$spot

  out <- integer(nrow(query_df))

  for (i in seq_len(nrow(query_df))) {
    d <- sqrt((t_xy[, 1] - q_xy[i, 1])^2 + (t_xy[, 2] - q_xy[i, 2])^2)

    if (exclude_self) {
      d[t_spot == query_df$spot[i]] <- Inf
    }

    out[i] <- sum(d <= radius, na.rm = TRUE)
  }

  out
}

safe_wilcox_compare <- function(df, group_col, feature, sample_id = "pooled", radius_multiplier = NA_real_) {

  group_levels <- c("C1_not_near_C3_high", "C1_near_C3_high")

  if (!all(c(group_col, feature) %in% colnames(df))) {
    return(data.frame(
      sample_id = sample_id,
      radius_multiplier = radius_multiplier,
      feature = feature,
      n_not_near = NA_integer_,
      n_near = NA_integer_,
      median_not_near = NA_real_,
      median_near = NA_real_,
      median_diff_near_minus_not = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  d <- df[, c(group_col, feature), drop = FALSE]
  colnames(d) <- c("group", "value")

  d <- d %>%
    dplyr::filter(group %in% group_levels) %>%
    dplyr::filter(!is.na(value))

  n_not <- sum(d$group == "C1_not_near_C3_high")
  n_near <- sum(d$group == "C1_near_C3_high")

  med_not <- ifelse(n_not > 0, median(d$value[d$group == "C1_not_near_C3_high"], na.rm = TRUE), NA_real_)
  med_near <- ifelse(n_near > 0, median(d$value[d$group == "C1_near_C3_high"], na.rm = TRUE), NA_real_)

  p_val <- NA_real_

  if (n_not >= 3 && n_near >= 3 && length(unique(d$value)) > 1) {
    p_val <- tryCatch(
      wilcox.test(value ~ group, data = d)$p.value,
      error = function(e) NA_real_
    )
  }

  data.frame(
    sample_id = sample_id,
    radius_multiplier = radius_multiplier,
    feature = feature,
    n_not_near = n_not,
    n_near = n_near,
    median_not_near = med_not,
    median_near = med_near,
    median_diff_near_minus_not = med_near - med_not,
    p_value = p_val,
    stringsAsFactors = FALSE
  )
}

p_to_star <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "ns")))
  )
}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

plot_spatial_neighbor_group <- function(obj, sid, group_col = "C1_C3_neighbor_group") {

  group_levels <- c(
    "Non_C1_high",
    "C1_not_near_C3_high",
    "C1_near_C3_high"
  )

  pal_group <- c(
    Non_C1_high = "#D9D9D9",
    C1_not_near_C3_high = "#2C7BB6",
    C1_near_C3_high = "#D7191C"
  )

  obj@meta.data[[group_col]] <- factor(
    as.character(obj@meta.data[[group_col]]),
    levels = group_levels
  )

  image_name <- names(obj@images)[1]

  p_list <- SpatialDimPlot(
    object = obj,
    group.by = group_col,
    images = image_name,
    image.alpha = 1,
    pt.size.factor = 1.65,
    crop = TRUE,
    cols = pal_group[group_levels],
    combine = FALSE
  )

  p <- p_list[[1]] +
    ggtitle(paste0(sid, ": C1-high spots near C3-high regions")) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 9, color = "black")
    )

  return(p)
}

########################################################
## 4. 主循环：每个样本做 C1/C3 空间邻近分析
########################################################

all_spatial_tables <- list()
all_c1_tables <- list()
all_sample_summary <- list()
all_per_sample_stats <- list()
all_sensitivity_stats <- list()
spatial_plots <- list()

for (i in seq_len(nrow(sample_info))) {

  sid <- sample_info$sample_id[i]
  rds_file <- sample_info$rds_file[i]

  message("\n==============================")
  message("Processing sample: ", sid)
  message("RDS: ", rds_file)
  message("==============================")

  obj <- readRDS(rds_file)
  DefaultAssay(obj) <- "Spatial"

  obj <- add_dlc1_raw_detected(obj, assay = "Spatial")

  meta <- obj@meta.data
  meta$spot <- rownames(meta)

  required_cols <- c("C1_score", "C3_score")

  if (!all(required_cols %in% colnames(meta))) {
    stop("样本 ", sid, " 缺少 C1_score 或 C3_score，请检查前一步脚本是否成功。")
  }

  coords <- extract_spatial_coords(obj)

  spatial_tbl <- meta %>%
    dplyr::inner_join(coords, by = "spot") %>%
    dplyr::mutate(sample_id = sid)

  ######################################################
  ## 4.1 定义 C1-high / C3-high
  ######################################################

  c1_thr <- as.numeric(quantile(
    spatial_tbl$C1_score,
    probs = c1_high_quantile,
    na.rm = TRUE
  ))

  c3_thr <- as.numeric(quantile(
    spatial_tbl$C3_score,
    probs = c3_high_quantile,
    na.rm = TRUE
  ))

  spatial_tbl <- spatial_tbl %>%
    dplyr::mutate(
      C1_high = ifelse(!is.na(C1_score) & C1_score >= c1_thr, TRUE, FALSE),
      C3_high = ifelse(!is.na(C3_score) & C3_score >= c3_thr, TRUE, FALSE)
    )

  ######################################################
  ## 4.2 估计 Visium spot 邻近距离
  ######################################################

  nn_dist <- estimate_nn_distance(spatial_tbl)
  primary_radius <- nn_dist * primary_radius_multiplier

  C3_high_tbl <- spatial_tbl %>%
    dplyr::filter(C3_high)

  ######################################################
  ## 4.3 primary 邻近分组：C1_near_C3_high vs C1_not_near_C3_high
  ######################################################

  spatial_tbl$nearest_C3_high_distance <- nearest_distance_to_targets(
    query_df = spatial_tbl,
    target_df = C3_high_tbl,
    exclude_self = FALSE
  )

  spatial_tbl$nearest_C3_high_distance_excluding_self <- nearest_distance_to_targets(
    query_df = spatial_tbl,
    target_df = C3_high_tbl,
    exclude_self = TRUE
  )

  spatial_tbl$n_C3_high_within_primary_radius <- count_targets_within_radius(
    query_df = spatial_tbl,
    target_df = C3_high_tbl,
    radius = primary_radius,
    exclude_self = FALSE
  )

  spatial_tbl <- spatial_tbl %>%
    dplyr::mutate(
      C1_near_C3_high = C1_high & n_C3_high_within_primary_radius >= min_C3_neighbors,
      C1_C3_neighbor_group = dplyr::case_when(
        !C1_high ~ "Non_C1_high",
        C1_high & C1_near_C3_high ~ "C1_near_C3_high",
        C1_high & !C1_near_C3_high ~ "C1_not_near_C3_high",
        TRUE ~ "Non_C1_high"
      ),
      C1_C3_neighbor_group = factor(
        C1_C3_neighbor_group,
        levels = c(
          "Non_C1_high",
          "C1_not_near_C3_high",
          "C1_near_C3_high"
        )
      )
    )

  ######################################################
  ## 4.4 sensitivity：不同邻近半径下重复分组和统计
  ######################################################

  sensitivity_stats_i <- list()

  for (rad_mul in sensitivity_radius_multipliers) {

    radius_i <- nn_dist * rad_mul
    group_col_i <- paste0("C1_C3_group_radius_", gsub("\\.", "_", rad_mul))
    count_col_i <- paste0("n_C3_high_radius_", gsub("\\.", "_", rad_mul))

    spatial_tbl[[count_col_i]] <- count_targets_within_radius(
      query_df = spatial_tbl,
      target_df = C3_high_tbl,
      radius = radius_i,
      exclude_self = FALSE
    )

    spatial_tbl[[group_col_i]] <- dplyr::case_when(
      !spatial_tbl$C1_high ~ "Non_C1_high",
      spatial_tbl$C1_high & spatial_tbl[[count_col_i]] >= min_C3_neighbors ~ "C1_near_C3_high",
      spatial_tbl$C1_high & spatial_tbl[[count_col_i]] < min_C3_neighbors ~ "C1_not_near_C3_high",
      TRUE ~ "Non_C1_high"
    )

    stats_rad_i <- lapply(compare_features, function(f) {
      safe_wilcox_compare(
        df = spatial_tbl,
        group_col = group_col_i,
        feature = f,
        sample_id = sid,
        radius_multiplier = rad_mul
      )
    }) %>%
      dplyr::bind_rows()

    sensitivity_stats_i[[as.character(rad_mul)]] <- stats_rad_i
  }

  all_sensitivity_stats[[sid]] <- dplyr::bind_rows(sensitivity_stats_i)

  ######################################################
  ## 4.5 样本 summary
  ######################################################

  sample_summary_i <- spatial_tbl %>%
    dplyr::summarise(
      sample_id = sid,
      n_spots = dplyr::n(),
      C1_threshold = c1_thr,
      C3_threshold = c3_thr,
      estimated_nn_distance = nn_dist,
      primary_radius_multiplier = primary_radius_multiplier,
      primary_radius = primary_radius,
      n_C1_high = sum(C1_high, na.rm = TRUE),
      n_C3_high = sum(C3_high, na.rm = TRUE),
      n_C1_near_C3_high = sum(C1_C3_neighbor_group == "C1_near_C3_high", na.rm = TRUE),
      n_C1_not_near_C3_high = sum(C1_C3_neighbor_group == "C1_not_near_C3_high", na.rm = TRUE),
      C1_near_C3_high_rate = n_C1_near_C3_high / n_C1_high
    )

  all_sample_summary[[sid]] <- sample_summary_i

  ######################################################
  ## 4.6 每个样本内统计比较
  ######################################################

  present_features <- intersect(compare_features, colnames(spatial_tbl))

  per_sample_stats_i <- lapply(present_features, function(f) {
    safe_wilcox_compare(
      df = spatial_tbl,
      group_col = "C1_C3_neighbor_group",
      feature = f,
      sample_id = sid,
      radius_multiplier = primary_radius_multiplier
    )
  }) %>%
    dplyr::bind_rows()

  all_per_sample_stats[[sid]] <- per_sample_stats_i

  ######################################################
  ## 4.7 写回 obj metadata，保存带分组注释的 RDS
  ######################################################

  add_cols <- spatial_tbl %>%
    dplyr::select(
      spot,
      C1_high,
      C3_high,
      nearest_C3_high_distance,
      nearest_C3_high_distance_excluding_self,
      n_C3_high_within_primary_radius,
      C1_near_C3_high,
      C1_C3_neighbor_group
    )

  for (cc in setdiff(colnames(add_cols), "spot")) {
    vv <- add_cols[[cc]]
    names(vv) <- add_cols$spot
    obj@meta.data[[cc]] <- vv[rownames(obj@meta.data)]
  }

  obj@meta.data$C1_C3_neighbor_group <- factor(
    as.character(obj@meta.data$C1_C3_neighbor_group),
    levels = c(
      "Non_C1_high",
      "C1_not_near_C3_high",
      "C1_near_C3_high"
    )
  )

  saveRDS(
    obj,
    file.path(rds_out_dir, paste0(sid, "_C1_C3_neighbor_annotated.rds"))
  )

  ######################################################
  ## 4.8 空间图：C1-high near/not near C3-high
  ######################################################

  p_spatial <- plot_spatial_neighbor_group(
    obj = obj,
    sid = sid,
    group_col = "C1_C3_neighbor_group"
  )

  save_pdf(
    p_spatial,
    file.path(plot_dir, paste0(sid, "_C1_high_near_C3_high_spatial_map.pdf")),
    width = 7.2,
    height = 5.8
  )

  spatial_plots[[sid]] <- p_spatial

  ######################################################
  ## 4.9 保存表格到 list
  ######################################################

  all_spatial_tables[[sid]] <- spatial_tbl

  all_c1_tables[[sid]] <- spatial_tbl %>%
    dplyr::filter(C1_high)
}

########################################################
## 5. 合并所有样本结果
########################################################

all_spatial_table <- dplyr::bind_rows(all_spatial_tables)
all_c1_table <- dplyr::bind_rows(all_c1_tables)
sample_summary <- dplyr::bind_rows(all_sample_summary)
per_sample_stats <- dplyr::bind_rows(all_per_sample_stats)
sensitivity_stats <- dplyr::bind_rows(all_sensitivity_stats)

per_sample_stats <- per_sample_stats %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significance = p_to_star(p_adj_BH)
  ) %>%
  dplyr::ungroup()

sensitivity_stats <- sensitivity_stats %>%
  dplyr::group_by(sample_id, radius_multiplier) %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significance = p_to_star(p_adj_BH)
  ) %>%
  dplyr::ungroup()

########################################################
## 6. pooled exploratory comparison
########################################################

present_features_all <- intersect(compare_features, colnames(all_c1_table))

pooled_stats <- lapply(present_features_all, function(f) {
  safe_wilcox_compare(
    df = all_c1_table,
    group_col = "C1_C3_neighbor_group",
    feature = f,
    sample_id = "ALL7_pooled_exploratory",
    radius_multiplier = primary_radius_multiplier
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significance = p_to_star(p_adj_BH)
  )

########################################################
## 7. sample-level effect test
## 更稳健：每个样本一个 median difference，再检验 7 个样本是否整体 > 0 或 != 0
########################################################

sample_effects <- per_sample_stats %>%
  dplyr::select(
    sample_id,
    feature,
    median_diff_near_minus_not
  ) %>%
  dplyr::filter(!is.na(median_diff_near_minus_not))

sample_level_tests <- sample_effects %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    median_effect = median(median_diff_near_minus_not, na.rm = TRUE),
    mean_effect = mean(median_diff_near_minus_not, na.rm = TRUE),
    n_positive_effect = sum(median_diff_near_minus_not > 0, na.rm = TRUE),
    p_value_two_sided = ifelse(
      n_samples >= 3 && length(unique(median_diff_near_minus_not)) > 1,
      tryCatch(
        wilcox.test(median_diff_near_minus_not, mu = 0, alternative = "two.sided")$p.value,
        error = function(e) NA_real_
      ),
      NA_real_
    ),
    p_value_greater = ifelse(
      n_samples >= 3 && length(unique(median_diff_near_minus_not)) > 1,
      tryCatch(
        wilcox.test(median_diff_near_minus_not, mu = 0, alternative = "greater")$p.value,
        error = function(e) NA_real_
      ),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj_BH_two_sided = p.adjust(p_value_two_sided, method = "BH"),
    p_adj_BH_greater = p.adjust(p_value_greater, method = "BH"),
    significance_two_sided = p_to_star(p_adj_BH_two_sided),
    significance_greater = p_to_star(p_adj_BH_greater)
  )

########################################################
## 8. 输出统计表格
########################################################

write.csv(
  sample_summary,
  file.path(table_dir, "01_sample_summary_C1_C3_neighbor.csv"),
  row.names = FALSE
)

write.csv(
  all_spatial_table,
  file.path(table_dir, "02_all_spots_C1_C3_neighbor_annotation.csv"),
  row.names = FALSE
)

write.csv(
  all_c1_table,
  file.path(table_dir, "03_C1_high_spots_near_C3_annotation.csv"),
  row.names = FALSE
)

write.csv(
  per_sample_stats,
  file.path(table_dir, "04_per_sample_wilcox_C1_near_vs_not_near.csv"),
  row.names = FALSE
)

write.csv(
  pooled_stats,
  file.path(table_dir, "05_pooled_exploratory_wilcox_C1_near_vs_not_near.csv"),
  row.names = FALSE
)

write.csv(
  sample_effects,
  file.path(table_dir, "06_sample_level_effects_near_minus_not.csv"),
  row.names = FALSE
)

write.csv(
  sample_level_tests,
  file.path(table_dir, "07_sample_level_effect_tests.csv"),
  row.names = FALSE
)

write.csv(
  sensitivity_stats,
  file.path(table_dir, "08_sensitivity_radius_wilcox_results.csv"),
  row.names = FALSE
)

########################################################
## 9. 作图：C1-high spots 两组比较
########################################################

plot_features <- intersect(
  c(
    "DLC1_expr",
    "MIF_receptor_score",
    "Integrin_receptor_score",
    "ECM_integrin_context_score",
    "Adhesion_ECM_score",
    "C3_score"
  ),
  colnames(all_c1_table)
)

all_c1_long <- all_c1_table %>%
  dplyr::filter(
    C1_C3_neighbor_group %in% c(
      "C1_not_near_C3_high",
      "C1_near_C3_high"
    )
  ) %>%
  dplyr::select(
    sample_id,
    spot,
    C1_C3_neighbor_group,
    dplyr::all_of(plot_features)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(plot_features),
    names_to = "feature",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    feature_label = ifelse(
      feature %in% names(feature_labels),
      feature_labels[feature],
      feature
    ),
    C1_C3_neighbor_group = factor(
      C1_C3_neighbor_group,
      levels = c(
        "C1_not_near_C3_high",
        "C1_near_C3_high"
      ),
      labels = c(
        "Not near C3-high",
        "Near C3-high"
      )
    )
  )


########################################################
## 9.2 绘图
########################################################

pal_compare <- c(
  "Not near C3-high" = "#2C7BB6",
  "Near C3-high" = "#D7191C"
)

p_violin_all <- ggplot(
  all_c1_long,
  aes(
    x = C1_C3_neighbor_group,
    y = value,
    fill = C1_C3_neighbor_group
  )
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    alpha = 0.75,
    color = "grey30",
    linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.95,
    color = "black",
    linewidth = 0.25
  ) +

  ## 显著性横线
  geom_segment(
    data = sig_df,
    aes(
      x = x_start,
      xend = x_end,
      y = y_bracket,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +

  ## 左侧短竖线
  geom_segment(
    data = sig_df,
    aes(
      x = x_start,
      xend = x_start,
      y = y_bracket - 0.03 * y_range,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +

  ## 右侧短竖线
  geom_segment(
    data = sig_df,
    aes(
      x = x_end,
      xend = x_end,
      y = y_bracket - 0.03 * y_range,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "black"
  ) +

  ## 显著性星号
  geom_text(
    data = sig_df,
    aes(
      x = x_mid,
      y = y_text,
      label = sig_label
    ),
    inherit.aes = FALSE,
    size = 4.2,
    fontface = "bold",
    color = "black"
  ) +

  facet_wrap(
    ~ feature_label,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(values = pal_compare) +
  labs(
    title = "C1-high spots grouped by proximity to C3-high regions",
    x = NULL,
    y = "Score / expression"
  ) +
  theme_pub(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      color = "black"
    )
  )

save_pdf(
  p_violin_all,
  file.path(plot_dir, "01_ALL7_C1_near_vs_not_near_C3_violin_box_with_significance.pdf"),
  width = 11.5,
  height = 8.2
)

########################################################
## 10. 作图：每个样本的分面箱线图
########################################################

p_box_sample <- ggplot(
  all_c1_long,
  aes(
    x = C1_C3_neighbor_group,
    y = value,
    fill = C1_C3_neighbor_group
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.85,
    linewidth = 0.25
  ) +
  facet_grid(feature_label ~ sample_id, scales = "free_y") +
  scale_fill_manual(values = pal_compare) +
  labs(
    title = "Per-sample comparison in C1-high spots",
    x = NULL,
    y = "Score / expression"
  ) +
  theme_pub(base_size = 8) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8)
  )

save_pdf(
  p_box_sample,
  file.path(plot_dir, "02_Per_sample_C1_near_vs_not_near_C3_boxplot.pdf"),
  width = 15.5,
  height = 10.5
)

########################################################
## 11. 作图：sample-level effect size
########################################################

sample_effects_plot <- sample_effects %>%
  dplyr::filter(feature %in% plot_features) %>%
  dplyr::mutate(
    feature_label = ifelse(
      feature %in% names(feature_labels),
      feature_labels[feature],
      feature
    )
  )

p_effect <- ggplot(
  sample_effects_plot,
  aes(
    x = feature_label,
    y = median_diff_near_minus_not
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey40") +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    fill = "grey85",
    color = "black",
    linewidth = 0.3
  ) +
  geom_point(
    aes(group = sample_id),
    position = position_jitter(width = 0.08, height = 0),
    size = 2.0,
    alpha = 0.85
  ) +
  labs(
    title = "Sample-level effect: C1_near_C3_high minus C1_not_near_C3_high",
    x = NULL,
    y = "Median difference"
  ) +
  theme_pub(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

save_pdf(
  p_effect,
  file.path(plot_dir, "03_Sample_level_effect_near_minus_not.pdf"),
  width = 10.5,
  height = 5.8
)

########################################################
## 12. 作图：C1-high 中 near C3-high 的比例
########################################################

p_rate <- ggplot(
  sample_summary,
  aes(
    x = sample_id,
    y = C1_near_C3_high_rate
  )
) +
  geom_col(
    width = 0.65,
    fill = "grey35",
    color = "black",
    linewidth = 0.25
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Proportion of C1-high spots located near C3-high regions",
    x = NULL,
    y = "C1 near C3-high rate"
  ) +
  theme_pub(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

save_pdf(
  p_rate,
  file.path(plot_dir, "04_C1_high_near_C3_high_rate_by_sample.pdf"),
  width = 7.5,
  height = 5.2
)

########################################################
## 13. 作图：空间图合并 PDF
########################################################

all_spatial_pdf <- file.path(plot_dir, "05_ALL7_C1_high_near_C3_high_spatial_maps.pdf")

cairo_pdf(
  filename = all_spatial_pdf,
  width = 7.2,
  height = 5.8,
  onefile = TRUE
)

for (sid in names(spatial_plots)) {
  print(spatial_plots[[sid]])
}

dev.off()

########################################################
## 14. 主图候选：HCC1R + HCC4R 空间邻近图
########################################################

rep_samples <- c("HCC1R", "HCC4R")
rep_samples <- rep_samples[rep_samples %in% names(spatial_plots)]

if (length(rep_samples) > 0) {

  p_rep <- wrap_plots(
    spatial_plots[rep_samples],
    ncol = length(rep_samples)
  ) +
    plot_annotation(
      title = "Representative spatial proximity of C1-high spots to C3-high regions",
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 15,
          hjust = 0.5
        )
      )
    )

  save_pdf(
    p_rep,
    file.path(plot_dir, "06_Representative_HCC1R_HCC4R_C1_C3_spatial_neighbor.pdf"),
    width = 13.5,
    height = 5.8
  )
}

########################################################
## 15. README
########################################################

writeLines(
  c(
    "DLC1 spatial transcriptomics C1/C3 neighbor analysis.",
    "",
    "Core logic:",
    "1. Define C1-high spots within each sample using the top 25% C1_score.",
    "2. Define C3-high spots within each sample using the top 25% C3_score.",
    "3. Estimate the native Visium nearest-neighbor spot distance.",
    "4. Classify C1-high spots as C1_near_C3_high if at least one C3-high spot is located within primary_radius = nearest_neighbor_distance * 1.5.",
    "5. Compare C1_near_C3_high versus C1_not_near_C3_high.",
    "",
    "Compared features:",
    paste(compare_features, collapse = ", "),
    "",
    "Main output tables:",
    "01_sample_summary_C1_C3_neighbor.csv",
    "02_all_spots_C1_C3_neighbor_annotation.csv",
    "03_C1_high_spots_near_C3_annotation.csv",
    "04_per_sample_wilcox_C1_near_vs_not_near.csv",
    "05_pooled_exploratory_wilcox_C1_near_vs_not_near.csv",
    "06_sample_level_effects_near_minus_not.csv",
    "07_sample_level_effect_tests.csv",
    "08_sensitivity_radius_wilcox_results.csv",
    "",
    "Main output plots:",
    "01_ALL7_C1_near_vs_not_near_C3_violin_box.pdf",
    "02_Per_sample_C1_near_vs_not_near_C3_boxplot.pdf",
    "03_Sample_level_effect_near_minus_not.pdf",
    "04_C1_high_near_C3_high_rate_by_sample.pdf",
    "05_ALL7_C1_high_near_C3_high_spatial_maps.pdf",
    "06_Representative_HCC1R_HCC4R_C1_C3_spatial_neighbor.pdf",
    "",
    "Interpretation note:",
    "If C1_near_C3_high spots show higher MIF_receptor_score, Integrin_receptor_score, ECM_integrin_context_score and Adhesion_ECM_score, this supports a spatial association between C3-like activated-matrix niches and C1-like fibroblast regions through MIF-CD74/CXCR4 and ECM-integrin-associated matrix context."
  ),
  con = file.path(out_dir, "README_C1_C3_neighbor_analysis.txt")
)

########################################################
## 16. 完成提示
########################################################

message("\n全部完成！")
message("输出目录：", out_dir)
message("主要统计表：", file.path(table_dir, "04_per_sample_wilcox_C1_near_vs_not_near.csv"))
message("合并统计表：", file.path(table_dir, "05_pooled_exploratory_wilcox_C1_near_vs_not_near.csv"))
message("样本层面效应表：", file.path(table_dir, "07_sample_level_effect_tests.csv"))
message("主图候选：", file.path(plot_dir, "06_Representative_HCC1R_HCC4R_C1_C3_spatial_neighbor.pdf"))
