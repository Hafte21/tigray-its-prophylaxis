###############################################################################
#
# ITS_V13.R - Single analytical pipeline for the Tigray ITS manuscript
#
# Manuscript: "Collapse of preventive HIV care during the Tigray war:
#              a 20-year interrupted time series"
#
# Author:  Hafte Kahsay Kebede (ORCID 0000-0002-1916-0819)
# Co-authors: Hailay Abrha Gesesew, Lillian Mwanri, Paul Ward
# License: Apache-2.0
#
# Reproducibility: R 4.6.0; primary CRAN packages listed in each part.
#
# This is the single comprehensive analytical script for the manuscript.
# It folds V13 (the originally submitted pipeline) and V13.1 (the
# reanalysis suite produced in response to internal review) into one
# file. No additional script is required to reproduce any number,
# table, or figure cited in the main text or supplementary appendix,
# with the single exception of the publication-quality plot
# generation, which is maintained separately as
# ITS_V13_Part2_Tables_Figures.R (see note at the end of this header).
#
# CONTENTS
# --------
# Part 1: V13 main analytical pipeline (data construction, ITS
#         modelling, and the V13 sensitivity suite). Produces all
#         primary tables and the master analytical workbook
#         ITS_V13_MASTER_RESULTS.xlsx.
#
# Part 2: V13.1 reanalysis (six targeted re-analyses produced in
#         response to internal cross-check of V13: negative binomial
#         under-recording sensitivity; Cox model n and event
#         extraction; per-period percentage of raw records resolving
#         to initiation events; drop-all-influential-quarters joint
#         sensitivity for BOTH cotrimoxazole and TPT; Table 1 cohort
#         filter forensics with the period-attribution bug fix; and
#         consolidated session info). Produces
#         ITS_V13_REANALYSIS_RESULTS.xlsx.
#
# Part 3: V13.1 proper-survival Cox analysis. Produces
#         ITS_V13_PROPER_SURVIVAL.xlsx with the Table S10a entries.
#
# DATA-CONSTRUCTION DIFFERENCE FROM V13
# -------------------------------------
# V13 parsed sex, DateOfBirth, and date_hiv_confirmed per row. If a
# patient demographic was blank on a given row, V13 treated it as
# missing for that row even when the same field was filled on another
# of the same patient's visits. V13.1 carries demographics forward
# within patient_id before de-duplication: for each unique patient,
# the first non-NA value seen in any of their records is propagated
# to all of their rows. A field is treated as truly missing only when
# it is blank on every row for the patient.
#
# RUN ORDER
# ---------
# 1. Edit data_dir at the top of Part 1 to your local path.
# 2. Source the file in a fresh R session:
#       source("ITS_V13.R")
#    Total runtime on a standard workstation (16 GB RAM, 4 cores) is
#    approximately 50 to 60 minutes for the complete sequence.
#
# NOTE ON V13 PART 2 (TABLES AND FIGURES GENERATION)
# --------------------------------------------------
# The V13 Part 2 script (Tables and Figures generation, approximately
# 4,100 additional lines) is maintained as a separate file
# (ITS_V13_Part2_Tables_Figures.R) in the GitHub repository and is
# sourced separately to keep this combined file under the 10,000-line
# threshold for editor convenience. Part 2 reads the checkpoint saved
# at the end of Part 1 and produces the publication-quality figures.
#
###############################################################################


###############################################################################
###############################################################################
###                                                                         ###
###                            PART 1: V13 MAIN PIPELINE                    ###
###                                                                         ###
###############################################################################
###############################################################################

################################################################################
#
#  ITS V13 - PART 1 of 2: Data construction, ITS modelling, sensitivity
#  Cotrimoxazole prophylaxis & TB preventive therapy during the Tigray conflict
#
#  Nature Communications | R 4.6.0 | CRAN-only packages
#
#  HOW TO RUN
#  ----------
#  Source Part 1 first, then Part 2:
#      source("ITS_V13_Part1_Analysis.R")
#      source("ITS_V13_Part2_Tables_Figures.R")
#  Part 1 saves a checkpoint at the end (v13_part1_checkpoint.RData) that
#  Part 2 loads if it cannot find the in-memory objects, so Part 2 is
#  restartable in a fresh R session.
#
#  FIGURE STANDARDS (Nature Communications style)
#  ------------------------------------------------------------
#     * dpi    = 600 (enforced via save_fig / save_fig_pub helpers)
#     * bg     = "white"
#     * units  = "mm"
#     * Lancet palette (#00468B, #AD002A, #925E9F, #42B540 ...)
#     * Sans-serif fonts (theme_lancet_pub)
#     * Main figures additionally exported as TIFF (LZW) via save_fig_pub
#
#  CHANGES IN V13 vs V12
#  ---------------------
#  Context: data were collected directly from each facility, including paper
#  records made during the conflict that were later transcribed to SmartCare.
#  Recording during the war may have been incomplete. V13 therefore adds
#  recording-completeness sensitivity, NOT EMR right-censoring analyses.
#
#  V13 adds five focused analyses without disturbing any V12 result:
#
#    Section 4B   Recording-pattern timeline (records, % initiated,
#                 ID/sex/age completeness, date-source distribution).
#    Section 13Z  TPT regimen-type filter sensitivity (refit without
#                 `tb_prophylaxis_type %in% c(1)` filter).
#    Section 13AA Under-recording sensitivity scenarios (inflate war-period
#                 observed counts by 10/25/50/75% and report war-effect
#                 threshold at which significance is lost).
#    Section 13BB ART-denominator adjustment for TPT (rate per 100 ART
#                 starts).
#    Section 13CC Year-level recording-quality summary.
#
#  Output folder bumped to Final_ITS_V13_Results so V12 results are preserved.
#
################################################################################
# ==============================================================================
# 0. SETUP
# ==============================================================================
rm(list = ls()); gc(); set.seed(42)
options(scipen = 999, digits = 4, warn = 1, stringsAsFactors = FALSE)
# Drain any OUTPUT sinks leftover from a prior failed run.
while (sink.number() > 0) sink()
flush.console()
# Script-wide timer
.script_t0 <- Sys.time()
tic <- function(label = "") {
  e <- as.numeric(difftime(Sys.time(), .script_t0, units = "secs"))
  cat(sprintf("[%6.1fs] %s\n", e, label))
  flush.console()
}
# Pre-emptive detach of packages that mask dplyr verbs
.detach_if_attached <- function(pkg) {
  pos <- paste0("package:", pkg)
  if (pos %in% search()) {
    suppressWarnings(suppressMessages(
      try(detach(pos, character.only = TRUE, unload = TRUE, force = TRUE),
          silent = TRUE)
    ))
  }
}
for (.p in c("MASS", "plyr", "raster", "Hmisc")) .detach_if_attached(.p)
rm(.p, .detach_if_attached)
cat("\n############################################################\n")
cat("#  ITS V13 PART 1 - Tigray HIV care continuum              #\n")
cat("############################################################\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "\n\n")
flush.console()
required_packages <- c(
  "tidyverse", "lubridate", "data.table", "zoo", "readxl",
  "forecast", "tseries", "segmented", "strucchange",
  "survival", "survminer",
  "nlme", "lmtest", "sandwich", "car",
  "nortest", "moments",
  "ggplot2", "ggtext", "patchwork", "scales", "viridis", "ggrepel",
  "broom", "naniar",
  "CausalImpact", "mgcv", "urca", "quantreg",
  "glmmTMB", "trend",
  "openxlsx"
)
ns_only_packages <- c("MASS", "pscl")
tic("Installing / loading packages")
for (pkg in c(required_packages, ns_only_packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
  }
}
suppressPackageStartupMessages(invisible(
  lapply(required_packages, function(p) {
    tryCatch(library(p, character.only = TRUE),
             error = function(e) message("Could not attach: ", p))
  })
))
suppressPackageStartupMessages(library(dplyr))  # last, wins masking
.check_select <- function() {
  s <- tryCatch(get("select", envir = globalenv(), inherits = TRUE),
                error = function(e) NULL)
  !is.null(s) && identical(s, dplyr::select)
}
if (!.check_select()) {
  for (.p in c("MASS", "plyr", "raster", "Hmisc")) {
    pos <- paste0("package:", .p)
    if (pos %in% search()) suppressWarnings(suppressMessages(
      try(detach(pos, character.only = TRUE, unload = TRUE, force = TRUE),
          silent = TRUE)))
  }
  suppressPackageStartupMessages(library(dplyr))
  if (!.check_select()) stop("dplyr::select still masked. Restart R.")
  message("Auto-remediated: detached masking package and re-attached dplyr.")
}
rm(.check_select)
cat("Packages OK.\n\n")
# ---- Paths ------------------------------------------------------------------
data_dir       <- "D:/PhD project Data/New data set/FINAL_ANALYSIS_ALL/HIV_care _outcomes/TB/incidence/tpt/prophylactic"
ctx_file_stems <- c("crtEthiopiaARTvisit_1", "crtEthiopiaARTvisit_2")
tpt_file_stems <- c("vcrtEthiopiaARTVisit_TB_All")
output_root <- file.path(data_dir, "Final_ITS_V13_Results")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
for (d in c("figures/main", "figures/supplement", "figures/diagnostics",
            "figures/cotri", "figures/tpt", "figures/v13",
            "tables", "diagnostics", "sensitivity", "survival",
            "quality_assurance", "data", "decomposition", "subgroup",
            "models", "validation")) {
  dir.create(file.path(output_root, d), showWarnings = FALSE, recursive = TRUE)
}
fig_main <- file.path(output_root, "figures/main")
fig_supp <- file.path(output_root, "figures/supplement")
fig_diag <- file.path(output_root, "figures/diagnostics")
fig_v13  <- file.path(output_root, "figures/v13")
config <- list(
  study_start   = ymd("2005-01-01"), study_end = ymd("2025-06-30"),
  covid_start   = ymd("2020-01-01"), covid_end = ymd("2020-12-31"),
  war_start     = ymd("2021-01-01"), war_end   = ymd("2022-12-31"),
  postwar_start = ymd("2023-01-01")
)
flow <- list()
# ==============================================================================
# 1. HELPERS  (see also the printed walkthrough)
# ==============================================================================
parse_date_flex <- function(x) {
  if (inherits(x, "Date"))   return(as.Date(x))
  if (inherits(x, "POSIXt")) return(as.Date(x))
  x_chr <- trimws(as.character(x))
  na_tokens <- c("", "NA", "N/A", "n/a", ".", "-", "NULL", "null",
                 "1/1/1900", "01/01/1900", "1900-01-01", "1900/1/1",
                 "1900/01/01", "1-1-1900", "01-01-1900")
  x_chr[x_chr %in% na_tokens] <- NA_character_
  fmts <- c("%m/%d/%Y", "%Y-%m-%d", "%d/%m/%Y", "%Y/%m/%d",
            "%m-%d-%Y", "%d-%m-%Y", "%m/%d/%y", "%d/%m/%y",
            "%Y%m%d", "%d-%b-%Y", "%d-%B-%Y")
  res <- as.Date(rep(NA, length(x_chr)))
  for (fmt in fmts) {
    idx <- is.na(res) & !is.na(x_chr); if (!any(idx)) break
    res[idx] <- suppressWarnings(as.Date(x_chr[idx], format = fmt))
  }
  na <- is.na(res) & !is.na(x_chr)
  if (any(na)) res[na] <- suppressWarnings(mdy(x_chr[na]))
  na <- is.na(res) & !is.na(x_chr)
  if (any(na)) {
    res[na] <- suppressWarnings(as.Date(parse_date_time(
      x_chr[na], orders = c("mdy", "dmy", "ymd"), quiet = TRUE)))
  }
  num_idx <- is.na(res) & suppressWarnings(!is.na(as.numeric(x_chr)))
  if (any(num_idx)) {
    res[num_idx] <- as.Date(as.numeric(x_chr[num_idx]),
                            origin = "1899-12-30")
  }
  res[!is.na(res) & res <= as.Date("1900-12-31")] <- NA
  res
}
read_artvisit <- function(stem, dir) {
  candidates <- c(file.path(dir, paste0(stem, ".csv")),
                  file.path(dir, paste0(stem, ".xlsx")),
                  file.path(dir, paste0(stem, ".xls")),
                  file.path(dir, paste0(stem, ".rds")))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop(sprintf("Could not find '%s' in:\n  %s", stem, dir))
  size_mb <- file.info(hit)$size / 1024^2
  cat(sprintf("  Reading: %s  (%.1f MB) ...\n", basename(hit), size_mb))
  flush.console()
  t0 <- Sys.time()
  ext <- tolower(tools::file_ext(hit))
  df <- switch(
    ext,
    csv  = data.table::fread(hit, na.strings = c("", "NA", "N/A"),
                             showProgress = FALSE,
                             colClasses = "character") |> as.data.frame(),
    xlsx = as.data.frame(readxl::read_excel(hit, guess_max = 100000,
                                            col_types = "text")),
    xls  = as.data.frame(readxl::read_excel(hit, guess_max = 100000,
                                            col_types = "text")),
    rds  = as.data.frame(readRDS(hit))
  )
  df[] <- lapply(df, function(col)
    if (is.character(col)) col else as.character(col))
  cat(sprintf("    \u2192 %s rows \u00d7 %d cols  in %.1fs\n",
              format(nrow(df), big.mark = ","), ncol(df),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  flush.console()
  df
}
get_col <- function(df, target) {
  hit <- names(df)[tolower(names(df)) == tolower(target)]
  if (length(hit) == 0) NA_character_ else hit[1]
}
# Nature Communications palette
lancet_colors <- c("Pre-war"  = "#00468BFF", "COVID-19" = "#925E9FFF",
                   "War"      = "#AD002AFF", "Post-war" = "#42B540FF")
period6_colors <- c("Early (2005-09)"    = "#1B9E77",
                    "Scale-up (2010-14)" = "#D95F02",
                    "Plateau (2015-19)"  = "#7570B3",
                    "COVID-19 (2020)"    = "#E7298A",
                    "War (2021-22)"      = "#E6AB02",
                    "Post-war (2023-25)" = "#66A61E")
lancet_lines  <- c("Fitted GLS"            = "#AD002AFF",
                   "Counterfactual GLS"    = "#00468BFF",
                   "Fitted NB"             = "#ED7D31FF",
                   "Counterfactual NB"     = "#925E9FFF",
                   "Counterfactual capped" = "#42B540FF")
theme_lancet_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace% theme(
    plot.title = element_blank(), plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    axis.title = element_text(face = "bold", size = rel(1)),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(size = rel(0.9), color = "grey20"),
    axis.line = element_line(color = "grey40", linewidth = 0.4),
    axis.ticks = element_line(color = "grey40", linewidth = 0.25),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = rel(0.85)),
    legend.text  = element_text(size = rel(0.8)),
    strip.background = element_rect(fill = "grey96", color = NA),
    strip.text = element_text(face = "bold", size = rel(0.9)),
    plot.margin = margin(8, 10, 8, 10)
  )
}
annots <- list(
  annotate("rect", xmin = config$covid_start, xmax = config$covid_end,
           ymin = -Inf, ymax = Inf, fill = "#925E9F", alpha = 0.12),
  annotate("rect", xmin = config$war_start, xmax = config$war_end,
           ymin = -Inf, ymax = Inf, fill = "#AD002A", alpha = 0.10),
  geom_vline(xintercept = config$covid_start, linetype = "dotted",
             color = "#925E9FFF", linewidth = 0.5),
  geom_vline(xintercept = config$war_start, linetype = "dashed",
             color = "#AD002AFF", linewidth = 0.6),
  geom_vline(xintercept = config$postwar_start, linetype = "dashed",
             color = "#42B540FF", linewidth = 0.6)
)
bp_idx <- function(df, yr, qtr = 1) which(df$year == yr & df$quarter == qtr)
save_csv <- function(x, file, ...) {
  utils::write.csv(x, file, ...)
  cat(sprintf("    saved: %s\n", normalizePath(file, mustWork = FALSE)))
  flush.console(); invisible(file)
}
save_rds <- function(object, file, ...) {
  base::saveRDS(object, file, ...)
  cat(sprintf("    saved: %s\n", normalizePath(file, mustWork = FALSE)))
  flush.console(); invisible(file)
}
# Nature Communications defaults: dpi=600, bg=white, units=mm
save_fig <- function(filename, plot, width = 180, height = 130,
                     dpi = 600, bg = "white", units = "mm", ...) {
  ggplot2::ggsave(filename, plot, width = width, height = height,
                  dpi = dpi, bg = bg, units = units, ...)
  cat(sprintf("    saved: %s\n", normalizePath(filename, mustWork = FALSE)))
  flush.console(); invisible(filename)
}
# Main-text figure saver: emits PNG + TIFF (LZW) at 600 dpi
save_fig_pub <- function(stem_no_ext, plot, width = 180, height = 130,
                         units = "mm", dpi = 600, also_pdf = FALSE) {
  png_file  <- paste0(stem_no_ext, ".png")
  tiff_file <- paste0(stem_no_ext, ".tiff")
  ggplot2::ggsave(png_file, plot, width = width, height = height,
                  dpi = dpi, bg = "white", units = units)
  ggplot2::ggsave(tiff_file, plot, width = width, height = height,
                  dpi = dpi, bg = "white", units = units,
                  device = "tiff", compression = "lzw")
  cat(sprintf("    saved: %s\n", normalizePath(png_file,  mustWork = FALSE)))
  cat(sprintf("    saved: %s\n", normalizePath(tiff_file, mustWork = FALSE)))
  if (also_pdf) {
    pdf_file <- paste0(stem_no_ext, ".pdf")
    tryCatch(
      ggplot2::ggsave(pdf_file, plot, width = width, height = height,
                      units = units, device = cairo_pdf, bg = "white"),
      error = function(e)
        ggplot2::ggsave(pdf_file, plot, width = width, height = height,
                        units = units, device = "pdf", bg = "white"))
    cat(sprintf("    saved: %s\n", normalizePath(pdf_file, mustWork = FALSE)))
  }
  flush.console(); invisible(stem_no_ext)
}
safe_print_summary <- function(m, label) {
  cat("\n═══", label, "═══\n\n")
  s <- tryCatch(summary(m), error = function(e) NULL,
                warning = function(w) NULL)
  if (is.null(s)) { cat("(summary() failed)\n"); return(invisible()) }
  tryCatch(
    print(s),
    error = function(e) {
      cat("(print(summary(.)) failed:", conditionMessage(e), ")\n")
      cat("Falling back to coefficients only:\n")
      tryCatch(print(coef(m)), error = function(e2)
        cat("(coef() also failed:", conditionMessage(e2), ")\n"))
    },
    warning = function(w) {
      cat("(print warning:", conditionMessage(w), ")\n")
      tryCatch(print(coef(m)), error = function(e) NULL)
    }
  )
  invisible()
}
run_subgroup_its <- function(event_data, date_col, quarter_grid,
                             subgroup_name, min_obs = 20) {
  sub_ts <- event_data %>%
    mutate(yr = year(.data[[date_col]]),
           qtr = quarter(.data[[date_col]])) %>%
    count(yr, qtr, name = "n_events") %>%
    right_join(quarter_grid, by = c("yr" = "year", "qtr" = "quarter")) %>%
    mutate(n_events = replace_na(n_events, 0)) %>% arrange(yr, qtr) %>%
    mutate(date = ymd(paste0(yr, "-", qtr * 3 - 2, "-01")),
           time_index = row_number(), season = factor(qtr),
           covid = as.numeric(yr == 2020),
           war_onset = as.numeric(date >= config$war_start),
           time_since_war = ifelse(war_onset == 1, as.numeric(difftime(
             date, config$war_start, units = "days")) / 91.25, 0),
           postwar = as.numeric(date >= config$postwar_start),
           time_since_postwar = ifelse(postwar == 1, as.numeric(difftime(
             date, config$postwar_start, units = "days")) / 91.25, 0))
  if (nrow(sub_ts) < min_obs || sum(sub_ts$war_onset == 1) < 4) return(NULL)
  m <- suppressWarnings(tryCatch(
    MASS::glm.nb(n_events ~ time_index + season + covid + war_onset +
                   time_since_war + postwar + time_since_postwar,
                 data = sub_ts,
                 control = glm.control(maxit = 300, epsilon = 1e-8)),
    error = function(e) NULL))
  model_type <- "NegBin"
  if (is.null(m) || (!is.null(m$converged) && !m$converged) ||
      any(is.na(coef(m))) || any(!is.finite(coef(m)))) {
    m <- suppressWarnings(tryCatch(
      gls(n_events ~ time_index + covid + war_onset + time_since_war +
            postwar + time_since_postwar, data = sub_ts,
          correlation = corAR1(form = ~time_index), method = "ML"),
      error = function(e) NULL))
    model_type <- "GLS-AR1"
  }
  if (is.null(m)) return(NULL)
  co <- tryCatch(coef(summary(m)), error = function(e) NULL)
  if (is.null(co) || !"war_onset" %in% rownames(co)) return(NULL)
  ec <- if ("Value" %in% colnames(co)) "Value" else "Estimate"
  sc <- if ("Std.Error" %in% colnames(co)) "Std.Error" else "Std. Error"
  ew <- co["war_onset", ec]; sw <- co["war_onset", sc]
  ep <- if ("postwar" %in% rownames(co)) co["postwar", ec] else NA
  sp <- if ("postwar" %in% rownames(co)) co["postwar", sc] else NA
  if (any(!is.finite(c(ew, sw)))) return(NULL)
  data.frame(subgroup = subgroup_name, war_estimate = ew, war_se = sw,
             war_ci_lower = ew - 1.96 * sw, war_ci_upper = ew + 1.96 * sw,
             postwar_estimate = ep, postwar_se = sp,
             postwar_ci_lower = ep - 1.96 * sp,
             postwar_ci_upper = ep + 1.96 * sp,
             n_quarters = nrow(sub_ts), total_events = sum(sub_ts$n_events),
             model_type = model_type)
}
# ==============================================================================
# 2. RAW DATA
# ==============================================================================
tic("2: DATA LOADING ")
cat("══ SECTION 2: DATA LOADING ══════════════════════════════════════════════\n")
# 2A. CTX SOURCE
ctx_raw_list <- lapply(ctx_file_stems, read_artvisit, dir = data_dir)
names(ctx_raw_list) <- ctx_file_stems
ctx_raw_list <- lapply(ctx_raw_list, function(d) { names(d) <- tolower(names(d)); d })
for (i in seq_along(ctx_raw_list)) {
  save_csv(data.frame(column = sort(names(ctx_raw_list[[i]]))),
           file.path(output_root, "quality_assurance",
                     paste0("columns_", ctx_file_stems[i], ".csv")),
           row.names = FALSE)
}
ctx_cols_needed <- c("uniqueartnumber", "sex", "dateofbirth", "date_hiv_confirmed",
                     "facility_level", "facility_ownership",
                     "cotri_dose_days_number", "cotrimoxazolestartdate",
                     "cortimoxazole_stop_date")
ctx_raw_list <- lapply(ctx_raw_list, function(d) {
  keep <- intersect(ctx_cols_needed, names(d))
  miss <- setdiff(ctx_cols_needed, keep)
  out  <- d[, keep, drop = FALSE]
  for (m in miss) out[[m]] <- NA_character_
  out[, ctx_cols_needed, drop = FALSE]
})
ctx_combined <- bind_rows(
  ctx_raw_list[[1]] %>% mutate(.source_file = ctx_file_stems[1]),
  ctx_raw_list[[2]] %>% mutate(.source_file = ctx_file_stems[2])
)
cat(sprintf("  CTX combined rows: %s (file 1: %s, file 2: %s)\n",
            format(nrow(ctx_combined), big.mark = ","),
            format(nrow(ctx_raw_list[[1]]), big.mark = ","),
            format(nrow(ctx_raw_list[[2]]), big.mark = ",")))
flow$ctx_file1 <- nrow(ctx_raw_list[[1]])
flow$ctx_file2 <- nrow(ctx_raw_list[[2]])
flow$ctx_combined <- nrow(ctx_combined)
n_before <- nrow(ctx_combined)
ctx_combined <- ctx_combined %>% distinct(across(-.source_file), .keep_all = TRUE)
cat(sprintf("  CTX: removed %s perfect-duplicate rows after merge.\n",
            format(n_before - nrow(ctx_combined), big.mark = ",")))
flow$ctx_dedup_removed     <- n_before - nrow(ctx_combined)
flow$ctx_after_merge_dedup <- nrow(ctx_combined)
# 2B. TPT SOURCE
tpt_raw_list <- lapply(tpt_file_stems, read_artvisit, dir = data_dir)
names(tpt_raw_list) <- tpt_file_stems
tpt_raw_list <- lapply(tpt_raw_list, function(d) { names(d) <- tolower(names(d)); d })
for (i in seq_along(tpt_raw_list)) {
  save_csv(data.frame(column = sort(names(tpt_raw_list[[i]]))),
           file.path(output_root, "quality_assurance",
                     paste0("columns_", tpt_file_stems[i], ".csv")),
           row.names = FALSE)
}
tpt_cols_needed <- c("uano", "sex", "dateofbirth", "art_start_date",
                     "registration_date",
                     "facility_level", "ownership",
                     "tb_prophylaxis_type", "inh_prophylaxisdurationmnth",
                     "inhprophylaxis_started_date",
                     "inhprophylaxisdiscontinueddate",
                     "inhprophylaxiscompleteddate")
tpt_raw_list <- lapply(tpt_raw_list, function(d) {
  keep <- intersect(tpt_cols_needed, names(d))
  miss <- setdiff(tpt_cols_needed, keep)
  out  <- d[, keep, drop = FALSE]
  for (m in miss) out[[m]] <- NA_character_
  out[, tpt_cols_needed, drop = FALSE]
})
tpt_raw <- tpt_raw_list[[1]]
cat(sprintf("  TPT raw rows: %s\n", format(nrow(tpt_raw), big.mark = ",")))
flow$tpt_raw <- nrow(tpt_raw)
n_before <- nrow(tpt_raw)
tpt_raw <- tpt_raw %>% distinct(.keep_all = TRUE)
cat(sprintf("  TPT: removed %s perfect-duplicate rows.\n\n",
            format(n_before - nrow(tpt_raw), big.mark = ",")))
flow$tpt_dedup_removed     <- n_before - nrow(tpt_raw)
flow$tpt_after_merge_dedup <- nrow(tpt_raw)
ctx_pid_col <- get_col(ctx_combined, "uniqueartnumber")
tpt_pid_col <- get_col(tpt_raw,      "uano")
if (is.na(ctx_pid_col)) warning("UniqueArtNumber not found in CTX files.")
if (is.na(tpt_pid_col)) warning("UANo not found in TPT file.")
cat(sprintf("  CTX patient ID: '%s'   TPT patient ID: '%s'\n\n",
            ctx_pid_col, tpt_pid_col))
c_pid <- if (!is.na(ctx_pid_col) || !is.na(tpt_pid_col)) "patient_id" else NA_character_
# ==============================================================================
# 3. BUILD PER-OUTCOME DATASETS
# ==============================================================================
tic("3: BUILD PER-OUTCOME DATASETS ")
cat("══ SECTION 3: BUILD PER-OUTCOME DATASETS ════════════════════════════════\n")
add_period_year_quarter <- function(df, date_col) {
  df %>% mutate(
    year = year(.data[[date_col]]),
    quarter = quarter(.data[[date_col]]),
    month = month(.data[[date_col]]),
    period = factor(case_when(
      year <= 2019 ~ "Pre-war", year == 2020 ~ "COVID-19",
      year <= 2022 ~ "War",     TRUE         ~ "Post-war"),
      levels = c("Pre-war", "COVID-19", "War", "Post-war")),
    period6 = factor(case_when(
      year <= 2009 ~ "Early (2005-09)",
      year <= 2014 ~ "Scale-up (2010-14)",
      year <= 2019 ~ "Plateau (2015-19)",
      year == 2020 ~ "COVID-19 (2020)",
      year <= 2022 ~ "War (2021-22)",
      TRUE         ~ "Post-war (2023-25)"),
      levels = c("Early (2005-09)", "Scale-up (2010-14)", "Plateau (2015-19)",
                 "COVID-19 (2020)", "War (2021-22)", "Post-war (2023-25)"))
  )
}
# 3A. COTRIMOXAZOLE
ctx_c_pid       <- get_col(ctx_combined, "uniqueartnumber")
ctx_c_sex       <- get_col(ctx_combined, "sex")
ctx_c_dob       <- get_col(ctx_combined, "dateofbirth")
ctx_c_hivconf   <- get_col(ctx_combined, "date_hiv_confirmed")
ctx_c_facility  <- get_col(ctx_combined, "facility_level")
ctx_c_ownership <- get_col(ctx_combined, "facility_ownership")
ctx_c_ctx_dose  <- get_col(ctx_combined, "cotri_dose_days_number")
ctx_c_ctx_start <- get_col(ctx_combined, "cotrimoxazolestartdate")
ctx_c_ctx_stop  <- get_col(ctx_combined, "cortimoxazole_stop_date")
cotri <- ctx_combined %>%
  mutate(
    patient_id = if (!is.na(ctx_c_pid)) trimws(as.character(.data[[ctx_c_pid]])) else NA_character_,
    sex_raw    = if (!is.na(ctx_c_sex)) toupper(trimws(.data[[ctx_c_sex]])) else NA_character_,
    sex = case_when(sex_raw == "F" ~ "Female",
                    sex_raw == "M" ~ "Male", TRUE ~ NA_character_),
    DateOfBirth        = if (!is.na(ctx_c_dob))     parse_date_flex(.data[[ctx_c_dob]])     else as.Date(NA),
    date_hiv_confirmed = if (!is.na(ctx_c_hivconf)) parse_date_flex(.data[[ctx_c_hivconf]]) else as.Date(NA),
    age = as.numeric(difftime(date_hiv_confirmed, DateOfBirth, units = "days")) / 365.25,
    age = ifelse(age < 0 | age > 110, NA_real_, age),
    facility_level     = if (!is.na(ctx_c_facility))  as.character(.data[[ctx_c_facility]])  else NA_character_,
    facility_ownership = if (!is.na(ctx_c_ownership)) as.character(.data[[ctx_c_ownership]]) else NA_character_,
    age_cat = cut(age, breaks = c(0, 14, 24, 34, 44, 54, Inf),
                  labels = c("<15", "15-24", "25-34", "35-44", "45-54", "55+"),
                  right = TRUE, include.lowest = TRUE),
    cotri_dose_days_number = if (!is.na(ctx_c_ctx_dose))
      suppressWarnings(as.numeric(.data[[ctx_c_ctx_dose]])) else NA_real_,
    CotrimoxazoleStartDate  = if (!is.na(ctx_c_ctx_start)) parse_date_flex(.data[[ctx_c_ctx_start]]) else as.Date(NA),
    cortimoxazole_stop_date = if (!is.na(ctx_c_ctx_stop))  parse_date_flex(.data[[ctx_c_ctx_stop]])  else as.Date(NA)
  ) %>%
  mutate(
    cotri_start_date = dplyr::coalesce(CotrimoxazoleStartDate, cortimoxazole_stop_date),
    date_source = case_when(
      !is.na(CotrimoxazoleStartDate)  ~ "initiation",
      !is.na(cortimoxazole_stop_date) ~ "stop_date",
      TRUE ~ "not_initiated"),
    initiated = !is.na(cotri_start_date) &
      cotri_start_date >= config$study_start &
      cotri_start_date <= config$study_end
  ) %>%
  mutate(cotri_start_date = ifelse(initiated, cotri_start_date, NA) %>%
           as.Date(origin = "1970-01-01")) %>%
  mutate(
    age = as.numeric(difftime(
      dplyr::coalesce(cotri_start_date, date_hiv_confirmed),
      DateOfBirth, units = "days")) / 365.25,
    age = ifelse(age < 0 | age > 110, NA_real_, age),
    age_cat = cut(age, breaks = c(0, 14, 24, 34, 44, 54, Inf),
                  labels = c("<15", "15-24", "25-34", "35-44", "45-54", "55+"),
                  right = TRUE, include.lowest = TRUE),
    time_to_init_days = as.numeric(cotri_start_date - date_hiv_confirmed),
    time_to_init_days = ifelse(time_to_init_days < 0, NA_real_, time_to_init_days)
  ) %>%
  add_period_year_quarter("cotri_start_date")
# V13 NEW: record_date anchor for Section 4B (recording-pattern timeline)
cotri$record_date <- cotri$date_hiv_confirmed
cat(sprintf("  CTX date sources: initiation=%s, stop_date=%s, not_initiated=%s\n",
            format(sum(cotri$date_source == "initiation",    na.rm = TRUE), big.mark = ","),
            format(sum(cotri$date_source == "stop_date",     na.rm = TRUE), big.mark = ","),
            format(sum(cotri$date_source == "not_initiated", na.rm = TRUE), big.mark = ",")))
flow$ctx_after_date_filter <- nrow(cotri)
cotri_has_id <- cotri %>% filter(!is.na(patient_id) & patient_id != "" & patient_id != "NA")
cotri_no_id  <- cotri %>% filter(is.na(patient_id) | patient_id == "" | patient_id == "NA")
n_has_id_pre <- nrow(cotri_has_id)
cotri_has_id <- cotri_has_id %>% distinct(patient_id, cotri_start_date, .keep_all = TRUE)
cotri <- bind_rows(cotri_has_id, cotri_no_id)
flow$ctx_n_initiated     <- sum(cotri$initiated)
flow$ctx_n_not_initiated <- sum(!cotri$initiated)
flow$ctx_blank_count     <- nrow(cotri_no_id)
flow$ctx_patdate_removed <- n_has_id_pre - nrow(cotri_has_id)
flow$ctx_final           <- nrow(cotri)
flow$ctx_patients        <- n_distinct(cotri_has_id$patient_id[cotri_has_id$initiated])
cat(sprintf("  CTX FINAL: %s initiated events (%s patients + %s unlinked) + %s not initiated = %s total\n",
            format(flow$ctx_n_initiated, big.mark = ","),
            format(flow$ctx_patients, big.mark = ","),
            format(sum(cotri_no_id$initiated), big.mark = ","),
            format(flow$ctx_n_not_initiated, big.mark = ","),
            format(nrow(cotri), big.mark = ",")))
save_csv(cotri, file.path(output_root, "data", "cotri_event_level.csv"),
         row.names = FALSE)
save_rds(cotri, file.path(output_root, "data", "cotri_event_level.rds"))
# 3B. TPT
tpt_c_pid       <- get_col(tpt_raw, "uano")
tpt_c_sex       <- get_col(tpt_raw, "sex")
tpt_c_dob       <- get_col(tpt_raw, "dateofbirth")
tpt_c_artstart  <- get_col(tpt_raw, "art_start_date")
tpt_c_regdate   <- get_col(tpt_raw, "registration_date")
tpt_c_facility  <- get_col(tpt_raw, "facility_level")
tpt_c_ownership <- get_col(tpt_raw, "ownership")
tpt_c_tb_type   <- get_col(tpt_raw, "tb_prophylaxis_type")
tpt_c_inh_dur   <- get_col(tpt_raw, "inh_prophylaxisdurationmnth")
tpt_c_inh_start <- get_col(tpt_raw, "inhprophylaxis_started_date")
tpt_c_inh_disc  <- get_col(tpt_raw, "inhprophylaxisdiscontinueddate")
tpt_c_inh_compl <- get_col(tpt_raw, "inhprophylaxiscompleteddate")
tpt <- tpt_raw %>%
  mutate(
    patient_id = if (!is.na(tpt_c_pid)) trimws(as.character(.data[[tpt_c_pid]])) else NA_character_,
    sex_raw    = if (!is.na(tpt_c_sex)) toupper(trimws(.data[[tpt_c_sex]])) else NA_character_,
    sex = case_when(sex_raw == "F" ~ "Female",
                    sex_raw == "M" ~ "Male", TRUE ~ NA_character_),
    DateOfBirth       = if (!is.na(tpt_c_dob))      parse_date_flex(.data[[tpt_c_dob]])      else as.Date(NA),
    art_start_date    = if (!is.na(tpt_c_artstart)) parse_date_flex(.data[[tpt_c_artstart]]) else as.Date(NA),
    registration_date = if (!is.na(tpt_c_regdate))  parse_date_flex(.data[[tpt_c_regdate]])  else as.Date(NA),
    age = as.numeric(difftime(registration_date, DateOfBirth, units = "days")) / 365.25,
    age = ifelse(age < 0 | age > 110, NA_real_, age),
    facility_level     = if (!is.na(tpt_c_facility))  as.character(.data[[tpt_c_facility]])  else NA_character_,
    facility_ownership = if (!is.na(tpt_c_ownership)) as.character(.data[[tpt_c_ownership]]) else NA_character_,
    age_cat = cut(age, breaks = c(0, 14, 24, 34, 44, 54, Inf),
                  labels = c("<15", "15-24", "25-34", "35-44", "45-54", "55+"),
                  right = TRUE, include.lowest = TRUE),
    tb_prophylaxis_type = if (!is.na(tpt_c_tb_type))
      suppressWarnings(as.numeric(.data[[tpt_c_tb_type]])) else NA_real_,
    inh_ProphylaxisDurationMnth = if (!is.na(tpt_c_inh_dur))
      suppressWarnings(as.numeric(.data[[tpt_c_inh_dur]])) else NA_real_,
    inhprophylaxis_started_date    = if (!is.na(tpt_c_inh_start)) parse_date_flex(.data[[tpt_c_inh_start]]) else as.Date(NA),
    InhprophylaxisDiscontinuedDate = if (!is.na(tpt_c_inh_disc))  parse_date_flex(.data[[tpt_c_inh_disc]])  else as.Date(NA),
    InhprophylaxisCompletedDate    = if (!is.na(tpt_c_inh_compl)) parse_date_flex(.data[[tpt_c_inh_compl]]) else as.Date(NA)
  ) %>%
  mutate(
    tpt_start_date = dplyr::coalesce(inhprophylaxis_started_date,
                                     InhprophylaxisDiscontinuedDate,
                                     InhprophylaxisCompletedDate),
    date_source = case_when(
      !is.na(inhprophylaxis_started_date)    ~ "initiation",
      !is.na(InhprophylaxisDiscontinuedDate) ~ "discontinued",
      !is.na(InhprophylaxisCompletedDate)    ~ "completed",
      TRUE ~ "not_initiated"),
    initiated = !is.na(tpt_start_date) &
      tpt_start_date >= config$study_start &
      tpt_start_date <= config$study_end &
      (is.na(tb_prophylaxis_type) | tb_prophylaxis_type %in% c(1))
  ) %>%
  mutate(tpt_start_date = ifelse(initiated, tpt_start_date, NA) %>%
           as.Date(origin = "1970-01-01")) %>%
  mutate(
    age = as.numeric(difftime(
      dplyr::coalesce(tpt_start_date, registration_date),
      DateOfBirth, units = "days")) / 365.25,
    age = ifelse(age < 0 | age > 110, NA_real_, age),
    age_cat = cut(age, breaks = c(0, 14, 24, 34, 44, 54, Inf),
                  labels = c("<15", "15-24", "25-34", "35-44", "45-54", "55+"),
                  right = TRUE, include.lowest = TRUE),
    time_to_init_days = as.numeric(tpt_start_date - art_start_date),
    time_to_init_days = ifelse(time_to_init_days < 0, NA_real_, time_to_init_days)
  ) %>%
  add_period_year_quarter("tpt_start_date")
# V13 NEW: record_date anchor for Section 4B
tpt$record_date <- tpt$registration_date
cat(sprintf("  TPT date sources: initiation=%s, discontinued=%s, completed=%s, not_initiated=%s\n",
            format(sum(tpt$date_source == "initiation",    na.rm = TRUE), big.mark = ","),
            format(sum(tpt$date_source == "discontinued",  na.rm = TRUE), big.mark = ","),
            format(sum(tpt$date_source == "completed",     na.rm = TRUE), big.mark = ","),
            format(sum(tpt$date_source == "not_initiated", na.rm = TRUE), big.mark = ",")))
flow$tpt_after_date_filter <- nrow(tpt)
tpt_has_id <- tpt %>% filter(!is.na(patient_id) & patient_id != "" & patient_id != "NA")
tpt_no_id  <- tpt %>% filter(is.na(patient_id) | patient_id == "" | patient_id == "NA")
n_has_id_pre <- nrow(tpt_has_id)
tpt_has_id <- tpt_has_id %>% distinct(patient_id, tpt_start_date, .keep_all = TRUE)
tpt <- bind_rows(tpt_has_id, tpt_no_id)
flow$tpt_n_initiated     <- sum(tpt$initiated)
flow$tpt_n_not_initiated <- sum(!tpt$initiated)
flow$tpt_blank_count     <- nrow(tpt_no_id)
flow$tpt_patdate_removed <- n_has_id_pre - nrow(tpt_has_id)
flow$tpt_final           <- nrow(tpt)
flow$tpt_patients        <- n_distinct(tpt_has_id$patient_id[tpt_has_id$initiated])
cat(sprintf("  TPT FINAL: %s initiated events (%s patients + %s unlinked) + %s not initiated = %s total\n\n",
            format(flow$tpt_n_initiated, big.mark = ","),
            format(flow$tpt_patients, big.mark = ","),
            format(sum(tpt_no_id$initiated), big.mark = ","),
            format(flow$tpt_n_not_initiated, big.mark = ","),
            format(nrow(tpt), big.mark = ",")))
save_csv(tpt, file.path(output_root, "data", "tpt_event_level.csv"),
         row.names = FALSE)
save_rds(tpt, file.path(output_root, "data", "tpt_event_level.rds"))
rm(ctx_combined, ctx_raw_list, tpt_raw, tpt_raw_list); invisible(gc())
# ==============================================================================
# 4. DATA QUALITY TIMELINE & TABLE 1
# ==============================================================================
tic("4: QUALITY TIMELINE & TABLE 1 ")
cat("══ SECTION 4: QUALITY TIMELINE & TABLE 1 ════════════════════════════════\n")
make_dq <- function(data, date_col, label) {
  data %>%
    mutate(yr = year(.data[[date_col]]), qtr = quarter(.data[[date_col]])) %>%
    group_by(yr, qtr) %>%
    summarise(n_events = n(),
              n_patients = if (!is.na(c_pid)) n_distinct(patient_id) else NA_integer_,
              pct_age_miss = mean(is.na(age)) * 100,
              pct_sex_miss = mean(is.na(sex)) * 100,
              pct_fac_miss = mean(is.na(facility_level)) * 100,
              .groups = "drop") %>%
    mutate(outcome = label)
}
save_csv(bind_rows(make_dq(cotri, "cotri_start_date", "Cotrimoxazole"),
                   make_dq(tpt,   "tpt_start_date",   "TPT")),
         file.path(output_root, "quality_assurance",
                   "data_quality_timeline.csv"), row.names = FALSE)
cotri_for_t1 <- if (!is.na(c_pid))
  cotri %>% filter(!is.na(patient_id) & patient_id != "" & patient_id != "NA") %>%
  arrange(patient_id, cotri_start_date) %>%
  group_by(patient_id) %>% slice(1) %>% ungroup() else cotri
tpt_for_t1   <- if (!is.na(c_pid))
  tpt %>% filter(!is.na(patient_id) & patient_id != "" & patient_id != "NA") %>%
  arrange(patient_id, tpt_start_date) %>%
  group_by(patient_id) %>% slice(1) %>% ungroup() else tpt
make_t1 <- function(d) {
  d %>% group_by(period) %>% group_modify(~ tibble(
    Variable = c("N", "Age, mean (SD)", "Age, median (IQR)",
                 "Female, n (%)", "Male, n (%)"),
    Value = c(
      as.character(nrow(.x)),
      sprintf("%.1f (%.1f)", mean(.x$age, na.rm = TRUE), sd(.x$age, na.rm = TRUE)),
      sprintf("%.1f (%.1f\u2013%.1f)", median(.x$age, na.rm = TRUE),
              quantile(.x$age, .25, na.rm = TRUE),
              quantile(.x$age, .75, na.rm = TRUE)),
      sprintf("%d (%.1f%%)", sum(.x$sex == "Female", na.rm = TRUE),
              sum(.x$sex == "Female", na.rm = TRUE) / nrow(.x) * 100),
      sprintf("%d (%.1f%%)", sum(.x$sex == "Male", na.rm = TRUE),
              sum(.x$sex == "Male", na.rm = TRUE) / nrow(.x) * 100))
  )) %>% ungroup() %>%
    pivot_wider(names_from = period, values_from = Value)
}
save_csv(make_t1(cotri_for_t1),
         file.path(output_root, "tables", "table1_cotri.csv"), row.names = FALSE)
save_csv(make_t1(tpt_for_t1),
         file.path(output_root, "tables", "table1_tpt.csv"), row.names = FALSE)
cat("  Table 1 saved.\n\n")
# ==============================================================================
# 4B. RECORDING-PATTERN TIMELINE  (V13 NEW)
# ==============================================================================
tic("4B: RECORDING-PATTERN TIMELINE (V13) ")
cat("══ SECTION 4B: RECORDING-PATTERN TIMELINE (V13 NEW) ═════════════════════\n")
records_timeline <- function(d, outcome_label) {
  has_id <- !is.na(d$patient_id) & d$patient_id != "" & d$patient_id != "NA"
  d$has_id  <- has_id
  d$has_sex <- !is.na(d$sex)
  d$has_age <- !is.na(d$age)
  d$stratum <- paste(ifelse(is.na(d$facility_level), "NA",
                            d$facility_level), "::",
                     ifelse(is.na(d$facility_ownership), "NA",
                            d$facility_ownership))
  d %>%
    filter(!is.na(record_date),
           record_date >= config$study_start,
           record_date <= config$study_end) %>%
    mutate(yr  = year(record_date),
           qtr = quarter(record_date)) %>%
    group_by(yr, qtr) %>%
    summarise(
      n_records_total = n(),
      n_initiated     = sum(initiated, na.rm = TRUE),
      pct_initiated   = 100 * sum(initiated, na.rm = TRUE) / n(),
      pct_with_id     = 100 * mean(has_id),
      pct_with_sex    = 100 * mean(has_sex),
      pct_with_age    = 100 * mean(has_age),
      n_strata_active = n_distinct(stratum),
      ds_initiation   = sum(date_source == "initiation"),
      ds_stop         = sum(date_source == "stop_date"),
      ds_disc         = sum(date_source == "discontinued"),
      ds_completed    = sum(date_source == "completed"),
      ds_not_init     = sum(date_source == "not_initiated"),
      .groups = "drop") %>%
    mutate(date    = ymd(paste0(yr, "-", qtr * 3 - 2, "-01")),
           period  = factor(case_when(
             yr <= 2019 ~ "Pre-war", yr == 2020 ~ "COVID-19",
             yr <= 2022 ~ "War", TRUE ~ "Post-war"),
             levels = c("Pre-war", "COVID-19", "War", "Post-war")),
           outcome = outcome_label) %>%
    arrange(yr, qtr)
}
rt_c <- records_timeline(cotri, "Cotrimoxazole")
rt_t <- records_timeline(tpt,   "TPT")
rt_all <- bind_rows(rt_c, rt_t)
save_csv(rt_all, file.path(output_root, "data",
                           "recording_timeline_quarterly.csv"),
         row.names = FALSE)
rt_period <- rt_all %>%
  group_by(outcome, period) %>%
  summarise(n_quarters         = n(),
            mean_records       = round(mean(n_records_total), 1),
            mean_initiated     = round(mean(n_initiated), 1),
            mean_pct_initiated = round(mean(pct_initiated), 1),
            mean_pct_id        = round(mean(pct_with_id), 1),
            mean_pct_sex       = round(mean(pct_with_sex), 1),
            mean_pct_age       = round(mean(pct_with_age), 1),
            mean_strata_active = round(mean(n_strata_active), 1),
            .groups = "drop")
save_csv(rt_period, file.path(output_root, "tables",
                              "table27_recording_quality_by_period.csv"),
         row.names = FALSE)
print(rt_period, row.names = FALSE)
ds_long_year <- function(rt, outcome_label) {
  rt %>%
    group_by(yr) %>%
    summarise(initiation   = sum(ds_initiation),
              stop         = sum(ds_stop),
              discontinued = sum(ds_disc),
              completed    = sum(ds_completed),
              not_init     = sum(ds_not_init),
              .groups = "drop") %>%
    pivot_longer(-yr, names_to = "source", values_to = "n") %>%
    mutate(outcome = outcome_label)
}
ds_all <- bind_rows(ds_long_year(rt_c, "Cotrimoxazole"),
                    ds_long_year(rt_t, "TPT")) %>%
  filter(n > 0)
save_csv(ds_all, file.path(output_root, "data",
                           "date_source_by_year.csv"),
         row.names = FALSE)
# ==============================================================================
# 5. QUARTERLY TIME SERIES
# ==============================================================================
tic("5: QUARTERLY TIME SERIES ")
cat("══ SECTION 5: QUARTERLY TIME SERIES ═════════════════════════════════════\n")
full_grid <- expand.grid(year = 2005:2025, quarter = 1:4) %>%
  filter(!(year == 2025 & quarter > 2)) %>% arrange(year, quarter)
build_ts <- function(data, date_col, label) {
  q <- data %>%
    filter(initiated == TRUE) %>%
    count(year, quarter, name = "n_events") %>%
    right_join(full_grid, by = c("year", "quarter")) %>%
    mutate(n_events = replace_na(n_events, 0)) %>% arrange(year, quarter) %>%
    mutate(date = ymd(paste0(year, "-", quarter * 3 - 2, "-01")),
           time_index = row_number(), season = factor(quarter),
           period = factor(case_when(
             year <= 2019 ~ "Pre-war", year == 2020 ~ "COVID-19",
             year <= 2022 ~ "War", TRUE ~ "Post-war"),
             levels = c("Pre-war", "COVID-19", "War", "Post-war")),
           period6 = factor(case_when(
             year <= 2009 ~ "Early (2005-09)",
             year <= 2014 ~ "Scale-up (2010-14)",
             year <= 2019 ~ "Plateau (2015-19)",
             year == 2020 ~ "COVID-19 (2020)",
             year <= 2022 ~ "War (2021-22)",
             TRUE         ~ "Post-war (2023-25)"),
             levels = c("Early (2005-09)", "Scale-up (2010-14)", "Plateau (2015-19)",
                        "COVID-19 (2020)", "War (2021-22)", "Post-war (2023-25)")),
           outcome = label)
  if (!is.na(c_pid)) {
    pat <- data %>%
      filter(initiated == TRUE) %>%
      mutate(yr = year(.data[[date_col]]),
             qtr = quarter(.data[[date_col]])) %>%
      group_by(yr, qtr) %>%
      summarise(n_patients = n_distinct(patient_id), .groups = "drop")
    q <- q %>% left_join(pat, by = c("year" = "yr", "quarter" = "qtr")) %>%
      mutate(n_patients = replace_na(n_patients, 0))
  } else q$n_patients <- NA_integer_
  bp_w <- bp_idx(q, 2021); bp_p <- bp_idx(q, 2023)
  q %>% mutate(
    trend_pre = time_index, covid = as.numeric(year == 2020),
    level_war = as.numeric(time_index >= bp_w),
    trend_war = ifelse(level_war == 1, time_index - bp_w + 1, 0),
    level_postwar = as.numeric(time_index >= bp_p),
    trend_postwar = ifelse(level_postwar == 1, time_index - bp_p + 1, 0),
    sin_q = sin(2 * pi * quarter / 4),
    cos_q = cos(2 * pi * quarter / 4))
}
cotri_q <- build_ts(cotri, "cotri_start_date", "Cotrimoxazole")
tpt_q   <- build_ts(tpt,   "tpt_start_date",   "TPT")
cat(sprintf("  CTX TS: %d quarters, %s initiated events\n",
            nrow(cotri_q), format(sum(cotri_q$n_events), big.mark = ",")))
cat(sprintf("  TPT TS: %d quarters, %s initiated events\n\n",
            nrow(tpt_q), format(sum(tpt_q$n_events), big.mark = ",")))
save_csv(cotri_q, file.path(output_root, "data", "cotri_quarterly_ts.csv"),
         row.names = FALSE)
save_csv(tpt_q,   file.path(output_root, "data", "tpt_quarterly_ts.csv"),
         row.names = FALSE)
# ==============================================================================
# 6. STL DECOMPOSITION
# ==============================================================================
tic("6: STL DECOMPOSITION ")
cat("══ SECTION 6: STL DECOMPOSITION ═════════════════════════════════════════\n")
do_stl <- function(ts_df, label, prefix) {
  stl_r <- stl(ts(ts_df$n_events, start = c(2005, 1), frequency = 4),
               s.window = "periodic", robust = TRUE)
  ts_df$trend     <- as.numeric(stl_r$time.series[, "trend"])
  ts_df$seasonal  <- as.numeric(stl_r$time.series[, "seasonal"])
  ts_df$remainder <- as.numeric(stl_r$time.series[, "remainder"])
  vr <- var(ts_df$remainder)
  str_df <- data.frame(
    Outcome = label,
    Trend_Strength    = round(max(0, 1 - vr / var(ts_df$trend + ts_df$remainder)), 3),
    Seasonal_Strength = round(max(0, 1 - vr / var(ts_df$seasonal + ts_df$remainder)), 3))
  save_csv(ts_df %>% dplyr::select(date, year, quarter, n_events,
                                   trend, seasonal, remainder),
           file.path(output_root, "decomposition", paste0(prefix, "_stl.csv")),
           row.names = FALSE)
  list(data = ts_df, strength = str_df)
}
dc <- do_stl(cotri_q, "Cotrimoxazole", "cotri"); cotri_q <- dc$data
dt <- do_stl(tpt_q,   "TPT",           "tpt");   tpt_q   <- dt$data
save_csv(bind_rows(dc$strength, dt$strength),
         file.path(output_root, "decomposition", "stl_strength.csv"),
         row.names = FALSE)
cat("  STL components extracted.\n\n")
# ==============================================================================
# 6B. SEASONALITY ANALYSIS
# ==============================================================================
tic("6B: SEASONALITY ANALYSIS ")
cat("══ SECTION 6B: SEASONALITY ANALYSIS ═════════════════════════════════════\n")
seasonal_indices <- function(ts_df, label) {
  overall <- ts_df %>%
    group_by(quarter) %>%
    summarise(mean_events = mean(n_events), sd_events = sd(n_events),
              median_events = median(n_events),
              ci_low = mean(n_events) - 1.96 * sd(n_events) / sqrt(n()),
              ci_high = mean(n_events) + 1.96 * sd(n_events) / sqrt(n()),
              .groups = "drop") %>%
    mutate(period = "Overall", outcome = label)
  by_period <- ts_df %>%
    group_by(period, quarter) %>%
    summarise(mean_events = mean(n_events), sd_events = sd(n_events),
              median_events = median(n_events),
              ci_low = mean(n_events) - 1.96 * sd(n_events) / sqrt(n()),
              ci_high = mean(n_events) + 1.96 * sd(n_events) / sqrt(n()),
              .groups = "drop") %>%
    mutate(outcome = label)
  grand_mean <- mean(ts_df$n_events)
  seasonal_idx <- ts_df %>%
    group_by(quarter) %>%
    summarise(seasonal_index = mean(n_events) / grand_mean, .groups = "drop") %>%
    mutate(outcome = label, grand_mean = grand_mean)
  list(overall = overall, by_period = by_period, indices = seasonal_idx)
}
si_c <- seasonal_indices(cotri_q, "Cotrimoxazole")
si_t <- seasonal_indices(tpt_q,   "TPT")
save_csv(bind_rows(si_c$overall, si_t$overall),
         file.path(output_root, "decomposition", "seasonal_means_overall.csv"),
         row.names = FALSE)
save_csv(bind_rows(si_c$by_period, si_t$by_period),
         file.path(output_root, "decomposition", "seasonal_means_by_period.csv"),
         row.names = FALSE)
seasonal_indices_6p <- function(ts_df, label) {
  ts_df %>%
    group_by(period6, quarter) %>%
    summarise(mean_events = mean(n_events), sd_events = sd(n_events),
              ci_low = mean(n_events) - 1.96 * sd(n_events) / sqrt(n()),
              ci_high = mean(n_events) + 1.96 * sd(n_events) / sqrt(n()),
              n_quarters = n(), .groups = "drop") %>%
    mutate(outcome = label)
}
si6_all <- bind_rows(seasonal_indices_6p(cotri_q, "Cotrimoxazole"),
                     seasonal_indices_6p(tpt_q,   "TPT"))
save_csv(si6_all,
         file.path(output_root, "decomposition", "seasonal_means_by_period6.csv"),
         row.names = FALSE)
save_csv(bind_rows(si_c$indices, si_t$indices),
         file.path(output_root, "decomposition", "seasonal_indices.csv"),
         row.names = FALSE)
season_tests <- function(ts_df, label) {
  kw <- tryCatch({ t <- kruskal.test(n_events ~ season, data = ts_df)
  data.frame(Test = "Kruskal-Wallis", Stat = t$statistic, P = t$p.value,
             DF = t$parameter) }, error = function(e) NULL)
  pre <- ts_df %>% filter(period == "Pre-war")
  fw <- tryCatch({ t <- kruskal.test(n_events ~ season, data = pre)
  data.frame(Test = "KW (pre-war only)", Stat = t$statistic, P = t$p.value,
             DF = t$parameter) }, error = function(e) NULL)
  lm_noseas <- lm(n_events ~ time_index, data = ts_df)
  lm_seas   <- lm(n_events ~ time_index + season, data = ts_df)
  ft <- tryCatch({ a <- anova(lm_noseas, lm_seas)
  data.frame(Test = "F-test (seasonal dummies)", Stat = a$F[2],
             P = a$`Pr(>F)`[2], DF = a$Df[2]) }, error = function(e) NULL)
  bind_rows(kw, fw, ft) %>% mutate(Outcome = label)
}
stest_all <- bind_rows(season_tests(cotri_q, "Cotrimoxazole"),
                       season_tests(tpt_q,   "TPT"))
save_csv(stest_all,
         file.path(output_root, "decomposition", "seasonality_tests.csv"),
         row.names = FALSE)
periodogram_data <- function(ts_df, label) {
  spec <- spectrum(ts(ts_df$n_events, frequency = 4), plot = FALSE, spans = c(3, 3))
  data.frame(Frequency = spec$freq, Period_Quarters = 1 / spec$freq,
             Spectral_Density = spec$spec, Outcome = label)
}
pg_all <- bind_rows(periodogram_data(cotri_q, "Cotrimoxazole"),
                    periodogram_data(tpt_q,   "TPT"))
save_csv(pg_all, file.path(output_root, "decomposition", "periodogram.csv"),
         row.names = FALSE)
season_war_interaction <- function(ts_df, label) {
  tryCatch({
    m <- suppressWarnings(gls(
      n_events ~ time_index + season + covid + level_war + trend_war +
        level_postwar + trend_postwar + season:level_war + season:level_postwar,
      data = ts_df, correlation = corAR1(form = ~time_index), method = "ML"))
    co <- coef(summary(m))
    data.frame(Parameter = rownames(co), Estimate = co[, "Value"],
               SE = co[, "Std.Error"],
               CI_Low = co[, "Value"] - 1.96 * co[, "Std.Error"],
               CI_High = co[, "Value"] + 1.96 * co[, "Std.Error"],
               P_Value = co[, "p-value"],
               Outcome = label, row.names = NULL)
  }, error = function(e) NULL)
}
swi_all <- bind_rows(season_war_interaction(cotri_q, "Cotrimoxazole"),
                     season_war_interaction(tpt_q,   "TPT"))
if (!is.null(swi_all) && nrow(swi_all) > 0) {
  save_csv(swi_all,
           file.path(output_root, "tables", "season_war_interaction.csv"),
           row.names = FALSE)
}
cotri_q$deseasonalised <- cotri_q$n_events - cotri_q$seasonal
tpt_q$deseasonalised   <- tpt_q$n_events   - tpt_q$seasonal
cat("  Seasonality analysis complete.\n\n")
# ==============================================================================
# 7. ALL ITS MODELS
# ==============================================================================
tic("7: ITS MODELS ")
cat("══ SECTION 7: ITS MODELS ════════════════════════════════════════════════\n")
fit_all_models <- function(ts_df, label, prefix) {
  ts_obj <- ts(ts_df$n_events, start = c(2005, 1), frequency = 4)
  xreg   <- as.matrix(ts_df[, c("trend_pre", "covid", "level_war", "trend_war",
                                "level_postwar", "trend_postwar")])
  models <- list()
  cat("  ", label, " · GLS-AR1\n")
  models$gls_ar1 <- suppressWarnings(tryCatch(
    gls(n_events ~ time_index + covid + level_war + trend_war +
          level_postwar + trend_postwar,
        data = ts_df, correlation = corAR1(form = ~time_index), method = "ML"),
    error = function(e) NULL))
  cat("  ", label, " · GLS-AR1 + season\n")
  models$gls_season <- suppressWarnings(tryCatch(
    gls(n_events ~ time_index + season + covid + level_war + trend_war +
          level_postwar + trend_postwar,
        data = ts_df, correlation = corAR1(form = ~time_index), method = "ML"),
    error = function(e) NULL))
  cat("  ", label, " · GLS-AR1 + season × war\n")
  models$gls_season_war <- suppressWarnings(tryCatch(
    gls(n_events ~ time_index + season + covid + level_war + trend_war +
          level_postwar + trend_postwar + season:level_war,
        data = ts_df, correlation = corAR1(form = ~time_index), method = "ML"),
    error = function(e) NULL))
  cat("  ", label, " · ARIMA(auto) + xreg\n")
  models$arima <- suppressWarnings(tryCatch(
    auto.arima(ts_obj, xreg = xreg, seasonal = TRUE,
               stepwise = FALSE, approximation = FALSE, trace = FALSE),
    error = function(e) NULL))
  cat("  ", label, " · Negative Binomial\n")
  models$nb <- suppressWarnings(tryCatch(
    MASS::glm.nb(n_events ~ time_index + season + covid + level_war + trend_war +
                   level_postwar + trend_postwar,
                 data = ts_df, control = glm.control(maxit = 300, epsilon = 1e-8)),
    error = function(e) NULL))
  if (!is.null(models$nb) && !isTRUE(models$nb$converged)) {
    cat("     NB did not converge → dropped\n"); models$nb <- NULL
  }
  cat("  ", label, " · Quasi-Poisson\n")
  models$qpois <- suppressWarnings(tryCatch(
    glm(n_events ~ time_index + season + covid + level_war + trend_war +
          level_postwar + trend_postwar,
        data = ts_df, family = quasipoisson),
    error = function(e) NULL))
  if (sum(ts_df$n_events == 0) > 0) {
    cat("  ", label, " · Zero-Inflated NB\n")
    models$zinb <- suppressWarnings(tryCatch(
      pscl::zeroinfl(n_events ~ time_index + covid + level_war + trend_war +
                       level_postwar + trend_postwar | 1,
                     data = ts_df, dist = "negbin",
                     control = pscl::zeroinfl.control(maxit = 500)),
      error = function(e) NULL))
  }
  cat("  ", label, " · Poisson (overdispersion check)\n")
  pois_check <- suppressWarnings(tryCatch(
    glm(n_events ~ time_index + covid + level_war + trend_war +
          level_postwar + trend_postwar, data = ts_df, family = poisson),
    error = function(e) NULL))
  models$poisson <- pois_check
  disp <- NA_real_
  if (!is.null(pois_check)) {
    disp <- sum(residuals(pois_check, type = "pearson")^2) / pois_check$df.residual
    cat(sprintf("     Dispersion ratio: %.2f → %s\n", disp,
                ifelse(disp > 1.5, "NB/QP warranted", "Gaussian acceptable")))
  }
  comp <- data.frame(
    Model = c("GLS-AR1", "GLS-AR1+season", "GLS-AR1+season×war",
              "ARIMA(auto)", "NegBin", "Quasi-Poisson", "ZINB", "Poisson"),
    AIC = sapply(list(models$gls_ar1, models$gls_season, models$gls_season_war,
                      models$arima, models$nb, models$qpois, models$zinb,
                      models$poisson),
                 function(m) if (!is.null(m)) tryCatch(AIC(m), error = function(e) NA) else NA),
    BIC = sapply(list(models$gls_ar1, models$gls_season, models$gls_season_war,
                      models$arima, models$nb, models$qpois, models$zinb,
                      models$poisson),
                 function(m) if (!is.null(m)) tryCatch(BIC(m), error = function(e) NA) else NA),
    Outcome = label, Dispersion = disp
  )
  save_csv(comp, file.path(output_root, "models",
                           paste0(prefix, "_model_comparison.csv")),
           row.names = FALSE)
  extract_war <- function(m, mname) {
    if (is.null(m)) return(NULL)
    co <- tryCatch(coef(summary(m)), error = function(e) NULL)
    if (is.null(co)) return(NULL)
    ec <- if ("Value" %in% colnames(co)) "Value"
    else if ("Estimate" %in% colnames(co)) "Estimate" else return(NULL)
    sc <- if ("Std.Error" %in% colnames(co)) "Std.Error"
    else if ("Std. Error" %in% colnames(co)) "Std. Error" else return(NULL)
    if (!"level_war" %in% rownames(co)) return(NULL)
    est <- co["level_war", ec]; se <- co["level_war", sc]
    if (any(!is.finite(c(est, se)))) return(NULL)
    pw <- if ("level_postwar" %in% rownames(co)) co["level_postwar", ec] else NA
    data.frame(Model = mname, War_Estimate = est, War_SE = se,
               War_CI_Low = est - 1.96 * se, War_CI_High = est + 1.96 * se,
               Postwar_Estimate = pw, Outcome = label)
  }
  war_effects <- bind_rows(
    extract_war(models$gls_ar1,        "GLS-AR1"),
    extract_war(models$gls_season,     "GLS-AR1+season"),
    extract_war(models$gls_season_war, "GLS-AR1+season×war"),
    extract_war(models$nb,             "NegBin"),
    extract_war(models$qpois,          "Quasi-Poisson"),
    extract_war(models$zinb,           "ZINB"),
    extract_war(models$poisson,        "Poisson")
  )
  if (!is.null(models$arima)) {
    ac <- coef(models$arima)
    as_ <- tryCatch(sqrt(diag(vcov(models$arima))), error = function(e) NULL)
    if (!is.null(as_) && "level_war" %in% names(ac)) {
      e <- ac["level_war"]; s <- as_["level_war"]
      if (all(is.finite(c(e, s)))) {
        war_effects <- bind_rows(war_effects, data.frame(
          Model = "ARIMA", War_Estimate = e, War_SE = s,
          War_CI_Low = e - 1.96 * s, War_CI_High = e + 1.96 * s,
          Postwar_Estimate = if ("level_postwar" %in% names(ac)) ac["level_postwar"] else NA,
          Outcome = label))
      }
    }
  }
  save_csv(war_effects, file.path(output_root, "models",
                                  paste0(prefix, "_war_effects.csv")),
           row.names = FALSE)
  sink_path <- file.path(output_root, "models", paste0(prefix, "_summaries.txt"))
  sink(sink_path)
  on.exit(if (sink.number() > 0) sink(), add = TRUE)
  for (nm in names(models)) {
    if (!is.null(models[[nm]])) safe_print_summary(models[[nm]], toupper(nm))
  }
  sink()
  on.exit()
  list(models = models, comparison = comp, war_effects = war_effects)
}
all_cotri <- fit_all_models(cotri_q, "Cotrimoxazole", "cotri")
all_tpt   <- fit_all_models(tpt_q,   "TPT",           "tpt")
m_cotri      <- all_cotri$models$gls_ar1
m_tpt        <- all_tpt$models$gls_ar1
nb_cotri     <- all_cotri$models$nb
nb_tpt       <- all_tpt$models$nb
arima_cotri  <- all_cotri$models$arima
arima_tpt    <- all_tpt$models$arima
stopifnot("Primary GLS-AR1 model failed for CTX" = !is.null(m_cotri))
stopifnot("Primary GLS-AR1 model failed for TPT" = !is.null(m_tpt))
cat("\n  All models fitted.\n\n")
# ==============================================================================
# 8. COUNTERFACTUAL & IMPACT
# ==============================================================================
tic("8: COUNTERFACTUALS & IMPACT ")
cat("══ SECTION 8: COUNTERFACTUALS & IMPACT ══════════════════════════════════\n")
calc_cf_gls <- function(ts_df, model) {
  ts_df$fitted_gls     <- as.numeric(fitted(model))
  ts_df$resid_norm_gls <- as.numeric(residuals(model, type = "normalized"))
  cf_data <- ts_df %>% mutate(covid = 0, level_war = 0, trend_war = 0,
                              level_postwar = 0, trend_postwar = 0)
  ts_df$counterfactual_gls <- predict(model, newdata = cf_data)
  ts_df
}
cotri_q <- calc_cf_gls(cotri_q, m_cotri)
tpt_q   <- calc_cf_gls(tpt_q,   m_tpt)
calc_cf_nb <- function(ts_df, model) {
  if (is.null(model)) {
    ts_df$counterfactual_nb <- NA_real_; ts_df$fitted_nb <- NA_real_; return(ts_df)
  }
  ts_df$fitted_nb <- predict(model, type = "response")
  cf_data <- ts_df %>% mutate(covid = 0, level_war = 0, trend_war = 0,
                              level_postwar = 0, trend_postwar = 0)
  ts_df$counterfactual_nb <- predict(model, newdata = cf_data, type = "response")
  ts_df
}
cotri_q <- calc_cf_nb(cotri_q, nb_cotri)
tpt_q   <- calc_cf_nb(tpt_q,   nb_tpt)
calc_cf_capped <- function(ts_df) {
  pre_anchor <- ts_df %>% filter(year %in% 2018:2019) %>%
    summarise(mu = mean(n_events)) %>% pull(mu)
  ts_df$counterfactual_capped <- ifelse(ts_df$year >= 2020, pre_anchor, NA_real_)
  ts_df
}
cotri_q <- calc_cf_capped(cotri_q)
tpt_q   <- calc_cf_capped(tpt_q)
calc_impact <- function(ts_df, label) {
  war_idx  <- which(ts_df$year %in% 2021:2022)
  post_idx <- which(ts_df$year >= 2023)
  obs_w  <- sum(ts_df$n_events[war_idx])
  obs_p  <- sum(ts_df$n_events[post_idx])
  exp_g_w <- sum(ts_df$counterfactual_gls[war_idx],  na.rm = TRUE)
  exp_g_p <- sum(ts_df$counterfactual_gls[post_idx], na.rm = TRUE)
  exp_n_w <- sum(ts_df$counterfactual_nb[war_idx],   na.rm = TRUE)
  exp_n_p <- sum(ts_df$counterfactual_nb[post_idx],  na.rm = TRUE)
  exp_c_w <- sum(ts_df$counterfactual_capped[war_idx],  na.rm = TRUE)
  exp_c_p <- sum(ts_df$counterfactual_capped[post_idx], na.rm = TRUE)
  tibble(
    Outcome = label,
    Period  = rep(c("War (2021-2022)", "Post-war (2023-2025)"), 3),
    CF_Source = rep(c("GLS-AR1", "NegBin", "Capped trend (2018-19 mean)"), each = 2),
    Observed = c(obs_w, obs_p, obs_w, obs_p, obs_w, obs_p),
    Expected = round(c(exp_g_w, exp_g_p, exp_n_w, exp_n_p, exp_c_w, exp_c_p)),
    Missed   = pmax(0, round(c(exp_g_w, exp_g_p, exp_n_w, exp_n_p,
                               exp_c_w, exp_c_p) -
                               c(obs_w, obs_p, obs_w, obs_p, obs_w, obs_p))),
    Pct_Change = round((c(obs_w, obs_p, obs_w, obs_p, obs_w, obs_p) -
                          c(exp_g_w, exp_g_p, exp_n_w, exp_n_p, exp_c_w, exp_c_p)) /
                         c(exp_g_w, exp_g_p, exp_n_w, exp_n_p, exp_c_w, exp_c_p) * 100, 1)
  )
}
impact <- bind_rows(calc_impact(cotri_q, "Cotrimoxazole"),
                    calc_impact(tpt_q,   "TPT"))
print(impact)
save_csv(impact, file.path(output_root, "tables", "impact_estimates.csv"),
         row.names = FALSE)
cat("\n")
# ==============================================================================
# 9. DIAGNOSTICS
# ==============================================================================
tic("9: DIAGNOSTICS ")
cat("══ SECTION 9: DIAGNOSTICS ═══════════════════════════════════════════════\n")
run_diag <- function(rv, label, prefix) {
  rv <- as.numeric(rv); rv <- rv[is.finite(rv)]
  if (length(rv) < 8) { cat("  ", label, ": too few residuals.\n"); return(NULL) }
  tests <- bind_rows(
    tryCatch({lb <- Box.test(rv, 20, "Ljung-Box")
    data.frame(Test = "Ljung-Box(20)", Stat = lb$statistic, P = lb$p.value)},
    error = function(e) NULL),
    tryCatch({sw <- shapiro.test(rv)
    data.frame(Test = "Shapiro-Wilk", Stat = sw$statistic, P = sw$p.value)},
    error = function(e) NULL),
    tryCatch({a <- tseries::adf.test(rv)
    data.frame(Test = "ADF", Stat = a$statistic, P = a$p.value)},
    error = function(e) NULL),
    suppressWarnings(tryCatch({k <- tseries::kpss.test(rv)
    data.frame(Test = "KPSS", Stat = k$statistic, P = k$p.value)},
    error = function(e) NULL)),
    suppressWarnings(tryCatch({p <- tseries::pp.test(rv)
    data.frame(Test = "Phillips-Perron", Stat = p$statistic, P = p$p.value)},
    error = function(e) NULL)),
    tryCatch({a <- nortest::ad.test(rv)
    data.frame(Test = "Anderson-Darling", Stat = a$statistic, P = a$p.value)},
    error = function(e) NULL),
    tryCatch({j <- moments::jarque.test(rv)
    data.frame(Test = "Jarque-Bera", Stat = j$statistic, P = j$p.value)},
    error = function(e) NULL),
    tryCatch({dw <- sum(diff(rv)^2) / sum(rv^2)
    data.frame(Test = "Durbin-Watson", Stat = dw, P = NA)},
    error = function(e) NULL)
  )
  tests$Outcome <- label
  save_csv(tests, file.path(output_root, "diagnostics",
                            paste0(prefix, "_tests.csv")), row.names = FALSE)
  # Lancet-styled ggplot diagnostics composite (replaces the previous
  # base-R png() + par(mfrow=c(3,3)) block, which rendered blank PNGs at
  # 600 dpi because per-panel line-based margins collapsed the plot area).
  color_pick <- if (grepl("^CTX", label)) "#00468BFF" else "#AD002AFF"
  diag_panels <- list(
    tryCatch({
      d <- data.frame(t = seq_along(rv), r = rv)
      ggplot(d, aes(t, r)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "#AD002A",
                   linewidth = 0.4) +
        geom_line(color = color_pick, linewidth = 0.45) +
        geom_point(color = color_pick, size = 0.9) +
        labs(x = "Observation index", y = "Residual",
             subtitle = "Residuals over time") +
        theme_lancet_pub() +
        theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
    }, error = function(e) NULL),
    tryCatch({
      d <- data.frame(r = rv); m <- mean(rv); s <- sd(rv)
      xs <- seq(min(rv) - s, max(rv) + s, length.out = 200)
      ndf <- data.frame(x = xs, y = dnorm(xs, m, s))
      ggplot(d, aes(r)) +
        geom_histogram(aes(y = after_stat(density)), bins = 25,
                       fill = color_pick, color = "white", alpha = 0.7) +
        geom_line(data = ndf, aes(x, y), color = "#AD002A", linewidth = 0.6) +
        labs(x = "Residual", y = "Density",
             subtitle = "Histogram with normal overlay") +
        theme_lancet_pub() +
        theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
    }, error = function(e) NULL),
    tryCatch({
      qq <- qqnorm(rv, plot.it = FALSE)
      d <- data.frame(theo = qq$x, samp = qq$y)
      ggplot(d, aes(theo, samp)) +
        geom_qq_line(aes(sample = samp), color = "#AD002A", linewidth = 0.5) +
        geom_point(color = color_pick, size = 1.1, alpha = 0.8) +
        labs(x = "Theoretical quantile", y = "Sample quantile",
             subtitle = "Normal Q-Q") +
        theme_lancet_pub() +
        theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
    }, error = function(e) NULL),
    tryCatch({
      ac <- acf(rv, lag.max = 20, plot = FALSE)
      d <- data.frame(lag = as.numeric(ac$lag), acf = as.numeric(ac$acf))
      d <- d[d$lag > 0, ]
      ci <- qnorm(0.975) / sqrt(length(rv))
      ggplot(d, aes(lag, acf)) +
        geom_hline(yintercept = c(-ci, ci), linetype = "dashed",
                   color = "#AD002A", linewidth = 0.35) +
        geom_hline(yintercept = 0, color = "grey40", linewidth = 0.25) +
        geom_segment(aes(xend = lag, yend = 0), color = color_pick,
                     linewidth = 0.55) +
        geom_point(color = color_pick, size = 1.4) +
        labs(x = "Lag (quarters)", y = "ACF", subtitle = "ACF") +
        theme_lancet_pub() +
        theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
    }, error = function(e) NULL),
    tryCatch({
      pc <- pacf(rv, lag.max = 20, plot = FALSE)
      d <- data.frame(lag = as.numeric(pc$lag), acf = as.numeric(pc$acf))
      ci <- qnorm(0.975) / sqrt(length(rv))
      ggplot(d, aes(lag, acf)) +
        geom_hline(yintercept = c(-ci, ci), linetype = "dashed",
                   color = "#AD002A", linewidth = 0.35) +
        geom_hline(yintercept = 0, color = "grey40", linewidth = 0.25) +
        geom_segment(aes(xend = lag, yend = 0), color = color_pick,
                     linewidth = 0.55) +
        geom_point(color = color_pick, size = 1.4) +
        labs(x = "Lag (quarters)", y = "PACF", subtitle = "PACF") +
        theme_lancet_pub() +
        theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
    }, error = function(e) NULL),
    tryCatch({
      d <- data.frame(t = seq_along(rv), r2 = rv^2)
      ggplot(d, aes(t, r2)) +
        geom_line(color = color_pick, linewidth = 0.45) +
        labs(x = "Observation index", y = expression(Residual^2),
             subtitle = "Squared residuals (ARCH check)") +
        theme_lancet_pub() +
        theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
    }, error = function(e) NULL),
    tryCatch({
      s <- cumsum(rv) / sqrt(sum(rv^2))
      d <- data.frame(t = seq_along(s), c = s)
      ggplot(d, aes(t, c)) +
        geom_hline(yintercept = c(-1.96, 1.96) / sqrt(length(rv)),
                   linetype = "dashed", color = "#AD002A", linewidth = 0.4) +
        geom_hline(yintercept = 0, color = "grey40", linewidth = 0.25) +
        geom_line(color = color_pick, linewidth = 0.55) +
        labs(x = "Observation index", y = "CUSUM",
             subtitle = "Standardised CUSUM") +
        theme_lancet_pub() +
        theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
    }, error = function(e) NULL),
    tryCatch({
      lags <- c(4, 8, 12, 16, 20)
      pvals <- sapply(lags, function(L) tryCatch(
        Box.test(rv, L, "Ljung-Box")$p.value, error = function(e) NA))
      d <- data.frame(lag = lags, p = pvals)
      ggplot(d, aes(lag, p)) +
        geom_hline(yintercept = 0.05, linetype = "dashed",
                   color = "#AD002A", linewidth = 0.4) +
        geom_line(color = color_pick, linewidth = 0.5) +
        geom_point(color = color_pick, size = 2.2) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(x = "Lag (quarters)", y = "p-value",
             subtitle = "Ljung-Box p-values") +
        theme_lancet_pub() +
        theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
    }, error = function(e) NULL))
  diag_panels <- diag_panels[!vapply(diag_panels, is.null, logical(1))]
  if (length(diag_panels) > 0) {
    fig <- Reduce(`+`, diag_panels) +
      patchwork::plot_layout(ncol = 3) +
      patchwork::plot_annotation(
        title = paste0(label, " - residual diagnostics"),
        theme = theme(plot.title = element_text(
          face = "bold", size = 12, color = color_pick,
          margin = margin(b = 6))))
    save_fig_pub(file.path(output_root, "diagnostics",
                           paste0(prefix, "_residuals")),
                 fig, width = 200, height = 200)
  }
  tests
}
diag_c_gls <- run_diag(residuals(m_cotri, type = "normalized"), "CTX GLS-AR1", "cotri_gls")
diag_t_gls <- run_diag(residuals(m_tpt,   type = "normalized"), "TPT GLS-AR1", "tpt_gls")
diag_c_nb  <- if (!is.null(nb_cotri))
  run_diag(residuals(nb_cotri, type = "deviance"), "CTX NegBin", "cotri_nb")
diag_t_nb  <- if (!is.null(nb_tpt))
  run_diag(residuals(nb_tpt,   type = "deviance"), "TPT NegBin", "tpt_nb")
save_csv(bind_rows(diag_c_gls, diag_t_gls, diag_c_nb, diag_t_nb),
         file.path(output_root, "diagnostics", "all_tests.csv"), row.names = FALSE)
cat("  Chow tests at war onset (2021 Q1):\n")
chow_results <- list()
for (nm in list(list(cotri_q, "Cotrimoxazole"), list(tpt_q, "TPT"))) {
  d <- nm[[1]]
  res <- tryCatch({
    sc <- strucchange::sctest(n_events ~ time_index, data = d,
                              type = "Chow", point = bp_idx(d, 2021))
    data.frame(Outcome = nm[[2]], Statistic = sc$statistic, P_Value = sc$p.value)
  }, error = function(e)
    data.frame(Outcome = nm[[2]], Statistic = NA, P_Value = NA))
  cat(sprintf("    %-15s  F = %.2f  p = %.4f\n",
              nm[[2]], res$Statistic, res$P_Value))
  chow_results[[nm[[2]]]] <- res
}
save_csv(bind_rows(chow_results),
         file.path(output_root, "diagnostics", "chow_tests.csv"), row.names = FALSE)
run_advanced_diag <- function(model, ts_df, label, prefix) {
  results <- list()
  lm_fit <- lm(n_events ~ time_index + covid + level_war + trend_war +
                 level_postwar + trend_postwar, data = ts_df)
  bg <- tryCatch({ t <- lmtest::bgtest(lm_fit, order = 4)
  data.frame(Test = "Breusch-Godfrey(4)", Stat = t$statistic, P = t$p.value)
  }, error = function(e) NULL)
  results <- c(results, list(bg))
  bp <- tryCatch({ t <- lmtest::bptest(lm_fit)
  data.frame(Test = "Breusch-Pagan", Stat = t$statistic, P = t$p.value)
  }, error = function(e) NULL)
  results <- c(results, list(bp))
  reset <- tryCatch({ t <- lmtest::resettest(lm_fit, power = 2:3)
  data.frame(Test = "Ramsey RESET", Stat = t$statistic, P = t$p.value)
  }, error = function(e) NULL)
  results <- c(results, list(reset))
  hac_comp <- tryCatch({
    ols_se <- coef(summary(lm_fit))[, "Std. Error"]
    hac_se <- sqrt(diag(sandwich::NeweyWest(lm_fit, lag = 4)))
    data.frame(Parameter = names(ols_se),
               OLS_SE = ols_se, HAC_SE = hac_se, Ratio = hac_se / ols_se,
               row.names = NULL)
  }, error = function(e) NULL)
  if (!is.null(hac_comp)) {
    save_csv(hac_comp,
             file.path(output_root, "diagnostics",
                       paste0(prefix, "_hac_vs_ols.csv")),
             row.names = FALSE)
  }
  cooks <- tryCatch({
    cd <- cooks.distance(lm_fit)
    influential <- which(cd > 4 / nrow(ts_df))
    png(file.path(output_root, "diagnostics",
                  paste0(prefix, "_influence.png")),
        3600, 2400, res = 600)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    plot(lm_fit, which = c(1, 2, 4, 5), pch = 19, cex = 0.6,
         col = ifelse(seq_along(cd) %in% influential, "red", "grey50"))
    dev.off()
    data.frame(Test = "Cook's D (n influential)",
               Stat = length(influential), P = NA)
  }, error = function(e) NULL)
  results <- c(results, list(cooks))
  adv_tests <- bind_rows(results); adv_tests$Outcome <- label
  save_csv(adv_tests,
           file.path(output_root, "diagnostics",
                     paste0(prefix, "_advanced_tests.csv")),
           row.names = FALSE)
  adv_tests
}
adv_c <- run_advanced_diag(m_cotri, cotri_q, "Cotrimoxazole", "cotri")
adv_t <- run_advanced_diag(m_tpt,   tpt_q,   "TPT",           "tpt")
save_csv(bind_rows(adv_c, adv_t),
         file.path(output_root, "diagnostics", "advanced_tests.csv"),
         row.names = FALSE)
cat("\n")
# ==============================================================================
# 10. SEGMENTED REGRESSION (6 BP forced)
# ==============================================================================
tic("10: SEGMENTED REGRESSION ")
cat("══ SECTION 10: SEGMENTED REGRESSION ═════════════════════════════════════\n")
fit_seg <- function(ts_df, label) {
  bps_full <- c(bp_idx(ts_df, 2010), bp_idx(ts_df, 2014), bp_idx(ts_df, 2018),
                bp_idx(ts_df, 2020), bp_idx(ts_df, 2021), bp_idx(ts_df, 2023))
  bps_full <- sort(unique(bps_full[!is.na(bps_full) & bps_full > 1 &
                                     bps_full < nrow(ts_df)]))
  lm0 <- lm(n_events ~ time_index, data = ts_df)
  seg <- suppressWarnings(tryCatch(
    segmented(lm0, seg.Z = ~time_index, psi = bps_full,
              control = seg.control(display = FALSE, it.max = 500,
                                    n.boot = 200, K = 4)),
    error = function(e) NULL))
  if (!is.null(seg) && inherits(seg, "segmented")) {
    chain_used <- "segmented (6 BP, iterative)"
    ts_df$seg_fitted <- as.numeric(predict(seg))
  } else {
    cat(sprintf("  %s: segmented(6 BP) failed; using fixed piecewise lm\n", label))
    pw_data <- ts_df
    for (i in seq_along(bps_full)) {
      bp <- bps_full[i]
      pw_data[[paste0("post_bp", i)]] <-
        ifelse(pw_data$time_index >= bp, pw_data$time_index - bp, 0)
    }
    pw_formula <- as.formula(paste0(
      "n_events ~ time_index + ",
      paste(paste0("post_bp", seq_along(bps_full)), collapse = " + ")))
    seg <- lm(pw_formula, data = pw_data)
    chain_used <- "fixed piecewise lm (6 BP)"
    ts_df$seg_fitted <- as.numeric(predict(seg))
  }
  list(model = seg, data = ts_df, chain = chain_used,
       breakpoint_indices = bps_full,
       breakpoint_dates = ts_df$date[bps_full])
}
sg_c <- fit_seg(cotri_q, "Cotrimoxazole"); cotri_q <- sg_c$data
sg_t <- fit_seg(tpt_q,   "TPT");           tpt_q   <- sg_t$data
save_csv(data.frame(Outcome = c("Cotrimoxazole", "TPT"),
                    Chain_Used = c(sg_c$chain, sg_t$chain)),
         file.path(output_root, "diagnostics", "segmented_chains.csv"),
         row.names = FALSE)
cat("\n")
# ==============================================================================
# 11. SURVIVAL
# ==============================================================================
tic("11: SURVIVAL ")
cat("══ SECTION 11: SURVIVAL ═════════════════════════════════════════════════\n")
run_surv <- function(data, time_col, label, prefix) {
  if (!time_col %in% names(data) || all(is.na(data[[time_col]]))) {
    cat(sprintf("  %s: '%s' unavailable. Skipped.\n", label, time_col)); return(NULL)
  }
  if (!is.na(c_pid)) {
    data <- data %>%
      filter(!is.na(patient_id) & patient_id != "" & patient_id != "NA") %>%
      arrange(patient_id, .data[[time_col]]) %>%
      group_by(patient_id) %>% slice(1) %>% ungroup()
    cat(sprintf("  %s: first-per-patient (n = %s)\n",
                label, format(nrow(data), big.mark = ",")))
  }
  sdf <- data %>% filter(!is.na(.data[[time_col]]),
                         .data[[time_col]] >= 0,
                         .data[[time_col]] <= 3650,
                         !is.na(period), !is.na(sex)) %>%
    mutate(event = 1, time_months = .data[[time_col]] / 30.44)
  if (nrow(sdf) < 50) { cat("    too few - skipped\n"); return(NULL) }
  km  <- survfit(Surv(time_months, event) ~ period, data = sdf)
  lr  <- survdiff(Surv(time_months, event) ~ period, data = sdf)
  cox <- coxph(Surv(time_months, event) ~ period + sex + age, data = sdf)
  ph  <- tryCatch(cox.zph(cox), error = function(e) NULL)
  ct  <- broom::tidy(cox, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(Outcome = label)
  save_csv(ct, file.path(output_root, "survival",
                         paste0(prefix, "_cox.csv")), row.names = FALSE)
  cox_fr <- tryCatch({
    if (!is.null(sdf$facility_level) &&
        sum(!is.na(sdf$facility_level)) > 100 &&
        length(unique(na.omit(sdf$facility_level))) >= 3) {
      coxph(Surv(time_months, event) ~ period + sex + age +
              frailty(facility_level, distribution = "gamma"),
            data = sdf %>% filter(!is.na(facility_level)))
    } else NULL
  }, error = function(e) NULL)
  if (!is.null(cox_fr)) {
    fr_summary <- tryCatch({
      s <- summary(cox_fr)
      data.frame(Term = rownames(s$coefficients),
                 Estimate = s$coefficients[, 1],
                 SE = s$coefficients[, "se(coef)"],
                 z = if ("z" %in% colnames(s$coefficients))
                   s$coefficients[, "z"] else NA_real_,
                 P = s$coefficients[, ncol(s$coefficients)],
                 HR = exp(s$coefficients[, 1]),
                 Outcome = label, row.names = NULL,
                 stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    if (!is.null(fr_summary))
      save_csv(fr_summary,
               file.path(output_root, "survival",
                         paste0(prefix, "_cox_frailty.csv")),
               row.names = FALSE)
  }
  if (!is.null(ph)) {
    tryCatch({
      png(file.path(output_root, "survival",
                    paste0(prefix, "_schoenfeld.png")),
          width = 4800, height = 3600, res = 600)
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); plot(ph); dev.off()
    }, error = function(e) cat("    Schoenfeld plot failed\n"))
  }
  km_p <- tryCatch(
    ggsurvplot(km, data = sdf, pval = TRUE, pval.method = TRUE,
               risk.table = TRUE, conf.int = TRUE,
               xlab = "Time (months)", ylab = "Probability of not initiating",
               legend.labs = levels(sdf$period),
               palette = unname(lancet_colors[levels(sdf$period)]),
               ggtheme = theme_lancet_pub(), title = ""),
    error = function(e) NULL)
  if (!is.null(km_p)) {
    save_fig(file.path(output_root, "survival", paste0(prefix, "_km.png")),
             km_p$plot, width = 200, height = 160)
  }
  cat(sprintf("    Log-rank χ² = %.2f\n", lr$chisq))
  if (!is.null(ph)) cat(sprintf("    cox.zph global p = %.4f\n",
                                ph$table[nrow(ph$table), "p"]))
  list(km = km, lr = lr, cox = cox, ph = ph, n = nrow(sdf))
}
sv_c <- run_surv(cotri, "time_to_init_days", "Cotrimoxazole", "cotri")
sv_t <- run_surv(tpt,   "time_to_init_days", "TPT",           "tpt")
cat("\n")
# ==============================================================================
# 12. SUBGROUPS
# ==============================================================================
tic("12: SUBGROUP ANALYSES ")
cat("══ SECTION 12: SUBGROUP ANALYSES ════════════════════════════════════════\n")
run_all_sub <- function(ed, dc, gmod, gdata, label) {
  ed <- ed %>% filter(initiated == TRUE)
  res <- list()
  co <- coef(summary(gmod))
  ec <- if ("Value" %in% colnames(co)) "Value" else "Estimate"
  sc <- if ("Std.Error" %in% colnames(co)) "Std.Error" else "Std. Error"
  ew <- co["level_war", ec]; sw <- co["level_war", sc]
  ep <- co["level_postwar", ec]; sp <- co["level_postwar", sc]
  res[[1]] <- data.frame(subgroup = "Overall", war_estimate = ew, war_se = sw,
                         war_ci_lower = ew - 1.96 * sw, war_ci_upper = ew + 1.96 * sw,
                         postwar_estimate = ep, postwar_se = sp,
                         postwar_ci_lower = ep - 1.96 * sp,
                         postwar_ci_upper = ep + 1.96 * sp,
                         n_quarters = nrow(gdata),
                         total_events = sum(gdata$n_events),
                         model_type = "GLS-AR1 (primary)")
  for (s in c("Male", "Female")) {
    r <- run_subgroup_its(ed %>% filter(sex == s), dc, full_grid, paste("Sex:", s))
    if (!is.null(r)) res[[length(res) + 1]] <- r
  }
  for (ac in levels(ed$age_cat)) {
    d <- ed %>% filter(age_cat == ac)
    if (nrow(d) >= 100) {
      r <- run_subgroup_its(d, dc, full_grid, paste("Age:", ac))
      if (!is.null(r)) res[[length(res) + 1]] <- r
    }
  }
  if (!all(is.na(ed$facility_level))) {
    for (fl in unique(na.omit(ed$facility_level))) {
      d <- ed %>% filter(facility_level == fl)
      if (nrow(d) >= 100) {
        r <- run_subgroup_its(d, dc, full_grid, paste("Facility:", fl))
        if (!is.null(r)) res[[length(res) + 1]] <- r
      }
    }
  }
  if (!all(is.na(ed$facility_ownership))) {
    for (fo in unique(na.omit(ed$facility_ownership))) {
      d <- ed %>% filter(facility_ownership == fo)
      if (nrow(d) >= 100) {
        r <- run_subgroup_its(d, dc, full_grid, paste("Ownership:", fo))
        if (!is.null(r)) res[[length(res) + 1]] <- r
      }
    }
  }
  bind_rows(res) %>% mutate(outcome = label)
}
sub_c <- run_all_sub(cotri, "cotri_start_date", m_cotri, cotri_q, "Cotrimoxazole")
sub_t <- run_all_sub(tpt,   "tpt_start_date",   m_tpt,   tpt_q,   "TPT")
save_csv(bind_rows(sub_c, sub_t),
         file.path(output_root, "subgroup", "subgroup_results.csv"),
         row.names = FALSE)
cat("  Subgroup analyses complete.\n\n")
# ==============================================================================
# 13. SENSITIVITY ANALYSES (13A-13Y from V12 + 13Z-13CC NEW V13)
# ==============================================================================
tic("13: SENSITIVITY ")
cat("══ SECTION 13: SENSITIVITY ══════════════════════════════════════════════\n")
cotri_ts_obj <- ts(cotri_q$n_events, start = c(2005, 1), frequency = 4)
tpt_ts_obj   <- ts(tpt_q$n_events,   start = c(2005, 1), frequency = 4)
xreg_cotri_mat <- as.matrix(cotri_q[, c("trend_pre", "covid", "level_war",
                                        "trend_war", "level_postwar", "trend_postwar")])
xreg_tpt_mat   <- as.matrix(tpt_q[,   c("trend_pre", "covid", "level_war",
                                        "trend_war", "level_postwar", "trend_postwar")])
# --- 13A. Alternative ARIMA orders ------------------------------------------
cat("  13A. Alternative ARIMA orders\n")
arima_specs <- list(c(0,1,1), c(1,0,0), c(1,0,1), c(1,1,0),
                    c(1,1,1), c(2,1,1), c(2,1,2))
arima_sens <- function(ts_obj, xreg_mat, label) {
  bind_rows(lapply(arima_specs, function(spec) suppressWarnings(tryCatch({
    fit <- Arima(ts_obj, order = spec, xreg = xreg_mat)
    co <- coef(fit)
    data.frame(Outcome = label,
               Order = paste0("(", paste(spec, collapse = ","), ")"),
               AIC = as.numeric(AIC(fit)), BIC = as.numeric(BIC(fit)),
               War_Effect = if ("level_war" %in% names(co))
                 as.numeric(co["level_war"]) else NA_real_,
               Converged = TRUE, stringsAsFactors = FALSE)
  }, error = function(e)
    data.frame(Outcome = label,
               Order = paste0("(", paste(spec, collapse = ","), ")"),
               AIC = NA_real_, BIC = NA_real_, War_Effect = NA_real_,
               Converged = FALSE, stringsAsFactors = FALSE))))) %>%
    filter(Converged) %>% arrange(AIC)
}
sa_all <- bind_rows(arima_sens(cotri_ts_obj, xreg_cotri_mat, "Cotrimoxazole"),
                    arima_sens(tpt_ts_obj,   xreg_tpt_mat,   "TPT"))
save_csv(sa_all, file.path(output_root, "sensitivity",
                           "arima_specifications.csv"), row.names = FALSE)
# --- 13B. GLS correlation structures ----------------------------------------
cat("  13B. GLS correlation structures\n")
corr_structures <- list(
  "None (OLS)" = NULL,
  "AR(1)"      = corAR1(form = ~time_index),
  "AR(2)"      = corARMA(form = ~time_index, p = 2, q = 0),
  "ARMA(1,1)"  = corARMA(form = ~time_index, p = 1, q = 1)
)
gls_formula <- n_events ~ time_index + season + covid + level_war +
  trend_war + level_postwar + trend_postwar
test_corr_struct <- function(ts_data, label) {
  bind_rows(lapply(names(corr_structures), function(nm)
    suppressWarnings(tryCatch({
      if (is.null(corr_structures[[nm]])) fit <- lm(gls_formula, data = ts_data)
      else fit <- gls(gls_formula, data = ts_data,
                      correlation = corr_structures[[nm]], method = "ML")
      co <- coef(summary(fit))
      ec <- if ("Value" %in% colnames(co)) "Value" else "Estimate"
      sc <- if ("Std.Error" %in% colnames(co)) "Std.Error" else "Std. Error"
      data.frame(Outcome = label, Structure = nm,
                 War_Estimate = as.numeric(co["level_war", ec]),
                 War_SE = as.numeric(co["level_war", sc]),
                 AIC = as.numeric(AIC(fit)), Converged = TRUE,
                 stringsAsFactors = FALSE)
    }, error = function(e)
      data.frame(Outcome = label, Structure = nm,
                 War_Estimate = NA_real_, War_SE = NA_real_,
                 AIC = NA_real_, Converged = FALSE,
                 stringsAsFactors = FALSE))))) %>%
    filter(Converged) %>% arrange(AIC)
}
sens_corr <- bind_rows(test_corr_struct(cotri_q, "Cotrimoxazole"),
                       test_corr_struct(tpt_q,   "TPT"))
save_csv(sens_corr, file.path(output_root, "sensitivity",
                              "correlation_structures.csv"), row.names = FALSE)
# --- 13C. Breakpoint timing (±2 quarters) -----------------------------------
cat("  13C. Breakpoint timing (±2 quarters)\n")
test_breakpoints <- function(ts_obj, xreg_base, bp_war_base, bp_post_base,
                             n_obs, label) {
  offsets <- c(-2, -1, 0, 1, 2); results <- list()
  for (o1 in offsets) for (o2 in offsets) {
    b1 <- bp_war_base + o1; b2 <- bp_post_base + o2
    if (b1 > 3 && b2 > b1 && b2 < n_obs - 2) {
      tryCatch({
        xr <- xreg_base
        xr[, "level_war"]     <- as.numeric(seq_len(n_obs) >= b1)
        xr[, "trend_war"]     <- ifelse(seq_len(n_obs) >= b1,
                                        seq_len(n_obs) - b1 + 1, 0)
        xr[, "level_postwar"] <- as.numeric(seq_len(n_obs) >= b2)
        xr[, "trend_postwar"] <- ifelse(seq_len(n_obs) >= b2,
                                        seq_len(n_obs) - b2 + 1, 0)
        fit <- suppressWarnings(Arima(ts_obj, order = c(1, 1, 1), xreg = xr))
        results[[length(results) + 1]] <- data.frame(
          Outcome = label, War_Offset = o1, Postwar_Offset = o2,
          War_BP_Q = paste0("2021 Q", 1 + o1),
          Postwar_BP_Q = paste0("2023 Q", 1 + o2),
          AIC = as.numeric(AIC(fit)), BIC = as.numeric(BIC(fit)),
          stringsAsFactors = FALSE)
      }, error = function(e) NULL)
    }
  }
  bind_rows(results) %>% arrange(AIC)
}
bp_sens <- bind_rows(
  test_breakpoints(cotri_ts_obj, xreg_cotri_mat,
                   bp_idx(cotri_q, 2021), bp_idx(cotri_q, 2023),
                   nrow(cotri_q), "Cotrimoxazole"),
  test_breakpoints(tpt_ts_obj,   xreg_tpt_mat,
                   bp_idx(tpt_q, 2021),   bp_idx(tpt_q, 2023),
                   nrow(tpt_q),   "TPT")
)
save_csv(bp_sens, file.path(output_root, "sensitivity",
                            "breakpoint_sensitivity.csv"), row.names = FALSE)
# --- 13D. COVID indicator on/off --------------------------------------------
cat("  13D. COVID indicator on/off\n")
fit_gls_covid <- function(data, include_covid, label) {
  suppressWarnings(tryCatch({
    f <- if (include_covid)
      n_events ~ time_index + season + covid + level_war + trend_war +
      level_postwar + trend_postwar
    else n_events ~ time_index + season + level_war + trend_war +
      level_postwar + trend_postwar
    fit <- gls(f, data = data, correlation = corAR1(form = ~time_index), method = "ML")
    co <- coef(summary(fit))
    ec <- if ("Value" %in% colnames(co)) "Value" else "Estimate"
    data.frame(Outcome = label,
               Model = if (include_covid) "GLS with COVID" else "GLS without COVID",
               War_Effect = as.numeric(co["level_war", ec]),
               War_IRR = NA_real_, AIC = as.numeric(AIC(fit)),
               Converged = TRUE, stringsAsFactors = FALSE)
  }, error = function(e)
    data.frame(Outcome = label,
               Model = if (include_covid) "GLS with COVID" else "GLS without COVID",
               War_Effect = NA_real_, War_IRR = NA_real_, AIC = NA_real_,
               Converged = FALSE, stringsAsFactors = FALSE)))
}
fit_nb_covid <- function(data, include_covid, label) {
  suppressWarnings(tryCatch({
    fit <- if (include_covid)
      MASS::glm.nb(n_events ~ time_index + covid + level_war + trend_war +
                     level_postwar + trend_postwar, data = data)
    else MASS::glm.nb(n_events ~ time_index + level_war + trend_war +
                        level_postwar + trend_postwar, data = data)
    data.frame(Outcome = label,
               Model = if (include_covid) "NB with COVID" else "NB without COVID",
               War_Effect = NA_real_,
               War_IRR = as.numeric(exp(coef(fit)["level_war"])),
               AIC = as.numeric(AIC(fit)), Converged = TRUE,
               stringsAsFactors = FALSE)
  }, error = function(e)
    data.frame(Outcome = label,
               Model = if (include_covid) "NB with COVID" else "NB without COVID",
               War_Effect = NA_real_, War_IRR = NA_real_, AIC = NA_real_,
               Converged = FALSE, stringsAsFactors = FALSE)))
}
covid_all <- bind_rows(
  fit_gls_covid(cotri_q, TRUE,  "Cotrimoxazole"),
  fit_gls_covid(cotri_q, FALSE, "Cotrimoxazole"),
  fit_gls_covid(tpt_q,   TRUE,  "TPT"),
  fit_gls_covid(tpt_q,   FALSE, "TPT"),
  fit_nb_covid (cotri_q, TRUE,  "Cotrimoxazole"),
  fit_nb_covid (cotri_q, FALSE, "Cotrimoxazole"),
  fit_nb_covid (tpt_q,   TRUE,  "TPT"),
  fit_nb_covid (tpt_q,   FALSE, "TPT")
)
save_csv(covid_all,
         file.path(output_root, "sensitivity", "covid_indicator_sensitivity.csv"),
         row.names = FALSE)
# --- 13E. Outlier sensitivity -----------------------------------------------
cat("  13E. Outlier sensitivity (>3 SD)\n")
outlier_sens <- function(data, ts_obj, xreg_mat, nb_model, arima_model, label) {
  thr <- mean(data$n_events, na.rm = TRUE) + 3 * sd(data$n_events, na.rm = TRUE)
  outliers <- which(data$n_events > thr)
  if (length(outliers) == 0) {
    return(data.frame(Outcome = label, N_Outliers = 0, Threshold = round(thr, 1),
                      Outlier_Quarters = "",
                      Original_NB_War_IRR = NA_real_, Clean_NB_War_IRR = NA_real_,
                      Original_ARIMA_War = NA_real_, Clean_ARIMA_War = NA_real_,
                      stringsAsFactors = FALSE))
  }
  orig_nb_irr <- if (!is.null(nb_model) && "level_war" %in% names(coef(nb_model)))
    as.numeric(exp(coef(nb_model)["level_war"])) else NA_real_
  orig_arima_war <- if (!is.null(arima_model) && "level_war" %in% names(coef(arima_model)))
    as.numeric(coef(arima_model)["level_war"]) else NA_real_
  data_clean <- data[-outliers, ]
  clean_nb_irr <- NA_real_
  if (!is.null(nb_model)) {
    nb_clean <- suppressWarnings(tryCatch(MASS::glm.nb(
      n_events ~ time_index + covid + level_war + trend_war +
        level_postwar + trend_postwar, data = data_clean),
      error = function(e) NULL))
    if (!is.null(nb_clean))
      clean_nb_irr <- as.numeric(exp(coef(nb_clean)["level_war"]))
  }
  clean_arima_war <- NA_real_
  if (!is.null(arima_model)) {
    xreg_clean <- xreg_mat[-outliers, , drop = FALSE]
    ts_clean <- ts(data_clean$n_events, start = c(2005, 1), frequency = 4)
    arima_clean <- suppressWarnings(tryCatch(
      Arima(ts_clean, order = c(1, 1, 1), xreg = xreg_clean),
      error = function(e) NULL))
    if (!is.null(arima_clean))
      clean_arima_war <- as.numeric(coef(arima_clean)["level_war"])
  }
  data.frame(Outcome = label, N_Outliers = length(outliers),
             Threshold = round(thr, 1),
             Outlier_Quarters = paste(outliers, collapse = ";"),
             Original_NB_War_IRR = orig_nb_irr, Clean_NB_War_IRR = clean_nb_irr,
             Original_ARIMA_War = orig_arima_war, Clean_ARIMA_War = clean_arima_war,
             stringsAsFactors = FALSE)
}
outlier_results <- bind_rows(
  outlier_sens(cotri_q, cotri_ts_obj, xreg_cotri_mat,
               nb_cotri, arima_cotri, "Cotrimoxazole"),
  outlier_sens(tpt_q,   tpt_ts_obj,   xreg_tpt_mat,
               nb_tpt,   arima_tpt,   "TPT"))
save_csv(outlier_results,
         file.path(output_root, "sensitivity", "outlier_sensitivity.csv"),
         row.names = FALSE)
# --- 13F. Newey-West HAC ----------------------------------------------------
cat("  13F. Newey-West HAC robust inference\n")
hac_war <- function(ts_df, label) {
  lm_fit <- lm(n_events ~ time_index + covid + level_war + trend_war +
                 level_postwar + trend_postwar, data = ts_df)
  co <- coef(lm_fit)
  hac_vcov <- sandwich::NeweyWest(lm_fit, lag = 4)
  hac_se <- sqrt(diag(hac_vcov))
  ols_se <- coef(summary(lm_fit))[, "Std. Error"]
  data.frame(Outcome = label, Parameter = names(co),
             Estimate = co, OLS_SE = ols_se, HAC_SE = hac_se,
             HAC_CI_Low = co - 1.96 * hac_se, HAC_CI_High = co + 1.96 * hac_se,
             row.names = NULL, stringsAsFactors = FALSE)
}
hac_all <- bind_rows(hac_war(cotri_q, "Cotrimoxazole"), hac_war(tpt_q, "TPT"))
save_csv(hac_all, file.path(output_root, "sensitivity", "newey_west_hac.csv"),
         row.names = FALSE)
# --- 13G. Pre-war trend validation ------------------------------------------
cat("  13G. Pre-war trend validation\n")
prewar_validation <- function(ts_df, label) {
  tryCatch({
    pre  <- ts_df %>% filter(year <= 2019)
    post <- ts_df %>% filter(year >= 2020)
    m_pre <- lm(n_events ~ time_index, data = pre)
    pred  <- predict(m_pre, newdata = post, interval = "prediction", level = 0.95)
    post$predicted <- pred[, "fit"]; post$pred_lwr <- pred[, "lwr"]; post$pred_upr <- pred[, "upr"]
    post$in_band  <- post$n_events >= post$pred_lwr & post$n_events <= post$pred_upr
    data.frame(Outcome = label,
               RMSE_PreWar = round(sqrt(mean(residuals(m_pre)^2)), 1),
               MAE_PostWar = round(mean(abs(post$n_events - post$predicted)), 1),
               Pct_In_95PI = round(mean(post$in_band) * 100, 1),
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}
preval <- bind_rows(prewar_validation(cotri_q, "Cotrimoxazole"),
                    prewar_validation(tpt_q,   "TPT"))
save_csv(preval, file.path(output_root, "sensitivity", "prewar_validation.csv"),
         row.names = FALSE)
# --- 13H. Cross-correlation -------------------------------------------------
cat("  13H. Cross-correlation (CTX × TPT)\n")
tryCatch({
  ccf_obj <- ccf(cotri_q$n_events, tpt_q$n_events, lag.max = 12, plot = FALSE)
  save_csv(data.frame(Lag = as.numeric(ccf_obj$lag), CCF = as.numeric(ccf_obj$acf)),
           file.path(output_root, "sensitivity", "cross_correlation.csv"),
           row.names = FALSE)
  png(file.path(output_root, "diagnostics", "cross_correlation_ctx_tpt.png"),
      3600, 2400, res = 600)
  par(mar = c(5, 5, 3, 1))
  plot(ccf_obj, main = "", xlab = "Lag (quarters)", ylab = "Cross-correlation",
       col = "#AD002AFF", lwd = 2, ci.col = "#00468B40")
  mtext("CTX vs TPT quarterly events", side = 3, line = 0.5, font = 2, cex = 1.1)
  dev.off()
}, error = function(e) cat("    cross-correlation failed\n"))
# --- 13I. Leave-one-quarter-out ---------------------------------------------
cat("  13I. Leave-one-quarter-out (LOQO)\n")
loqo_check <- function(ts_df, label) {
  war_est <- numeric(nrow(ts_df))
  for (i in seq_len(nrow(ts_df))) {
    d_loo <- ts_df[-i, ]
    m_loo <- suppressWarnings(tryCatch(
      gls(n_events ~ time_index + covid + level_war + trend_war +
            level_postwar + trend_postwar,
          data = d_loo, correlation = corAR1(form = ~time_index), method = "ML"),
      error = function(e) NULL))
    if (!is.null(m_loo)) {
      war_est[i] <- coef(summary(m_loo))["level_war", "Value"]
    } else war_est[i] <- NA
  }
  data.frame(Outcome = label, Quarter_Removed = ts_df$date,
             War_Estimate = war_est, stringsAsFactors = FALSE) %>%
    filter(!is.na(War_Estimate))
}
loqo_all <- bind_rows(loqo_check(cotri_q, "Cotrimoxazole"),
                      loqo_check(tpt_q,   "TPT"))
save_csv(loqo_all, file.path(output_root, "sensitivity", "loqo_robustness.csv"),
         row.names = FALSE)
# --- 13J. CausalImpact ------------------------------------------------------
cat("  13J. CausalImpact (Bayesian structural time series)\n")
ci_run <- function(ts_df, label, prefix) {
  tryCatch({
    if (!requireNamespace("CausalImpact", quietly = TRUE)) return(NULL)
    suppressPackageStartupMessages(library(CausalImpact))
    z <- zoo::zoo(ts_df$n_events, order.by = ts_df$date)
    pre_per  <- range(ts_df$date[ts_df$year <= 2020])
    post_per <- range(ts_df$date[ts_df$year >= 2021])
    ci <- CausalImpact::CausalImpact(
      z, pre.period = pre_per, post.period = post_per,
      model.args = list(niter = 2000, nseasons = 4))
    sm <- ci$summary
    out <- data.frame(
      Outcome = label, Metric = rownames(sm),
      Average_Actual = sm$Actual, Average_Predicted = sm$Pred,
      Pred_CI_Low = sm$Pred.lower, Pred_CI_High = sm$Pred.upper,
      Abs_Effect = sm$AbsEffect,
      Abs_Effect_Low = sm$AbsEffect.lower,
      Abs_Effect_High = sm$AbsEffect.upper,
      Rel_Effect_Pct = sm$RelEffect * 100,
      Rel_CI_Low_Pct = sm$RelEffect.lower * 100,
      Rel_CI_High_Pct = sm$RelEffect.upper * 100,
      P_Value = ci$summary$p[1],
      row.names = NULL, stringsAsFactors = FALSE)
    save_csv(out, file.path(output_root, "sensitivity",
                            paste0(prefix, "_causalimpact.csv")),
             row.names = FALSE)
    saveRDS(ci, file.path(output_root, "models",
                          paste0(prefix, "_causalimpact.rds")))
    ci
  }, error = function(e) {
    cat("    CausalImpact failed:", conditionMessage(e), "\n"); NULL
  })
}
ci_cotri <- ci_run(cotri_q, "Cotrimoxazole", "cotri")
ci_tpt   <- ci_run(tpt_q,   "TPT",           "tpt")
# --- 13K. Block-bootstrap impact CIs ----------------------------------------
cat("  13K. Block-bootstrap cumulative impact CIs\n")
boot_impact <- function(ts_df, label, B = 1000, block_len = 4) {
  tryCatch({
    war_idx  <- which(ts_df$year %in% 2021:2022)
    post_idx <- which(ts_df$year >= 2023)
    if (length(war_idx) < block_len) return(NULL)
    obs_w  <- ts_df$n_events[war_idx]
    cf_w   <- ts_df$counterfactual_gls[war_idx]
    obs_p  <- ts_df$n_events[post_idx]
    cf_p   <- ts_df$counterfactual_gls[post_idx]
    if (any(is.na(cf_w)) || any(is.na(cf_p))) return(NULL)
    block_resample <- function(x) {
      if (length(x) < block_len) return(x)
      n_blocks <- ceiling(length(x) / block_len)
      starts   <- sample(seq_len(length(x) - block_len + 1),
                         n_blocks, replace = TRUE)
      out <- unlist(lapply(starts, function(s) x[s:(s + block_len - 1)]))
      out[seq_along(x)]
    }
    set.seed(42)
    boots_w <- replicate(B, sum(block_resample(obs_w) - cf_w))
    boots_p <- replicate(B, sum(block_resample(obs_p) - cf_p))
    data.frame(
      Outcome = label,
      Period  = c("War (2021-2022)", "Post-war (2023-2025)"),
      Observed = c(sum(obs_w), sum(obs_p)),
      Expected_GLS = round(c(sum(cf_w), sum(cf_p))),
      Missed_Mean = round(c(mean(boots_w), mean(boots_p))),
      Missed_CI_Low  = round(c(quantile(boots_w, 0.025),
                               quantile(boots_p, 0.025))),
      Missed_CI_High = round(c(quantile(boots_w, 0.975),
                               quantile(boots_p, 0.975))),
      B_replicates = B, Block_length_q = block_len,
      stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}
bs_all <- bind_rows(boot_impact(cotri_q, "Cotrimoxazole"),
                    boot_impact(tpt_q,   "TPT"))
if (!is.null(bs_all) && nrow(bs_all) > 0) {
  save_csv(bs_all, file.path(output_root, "sensitivity",
                             "bootstrap_impact_ci.csv"), row.names = FALSE)
}
# --- 13L. Placebo test ------------------------------------------------------
cat("  13L. Placebo / permutation test\n")
placebo_test <- function(ts_df, actual_war_est, label) {
  tryCatch({
    sham_qs <- which(ts_df$year %in% 2010:2019)
    if (length(sham_qs) < 5) return(NULL)
    estimates <- sapply(sham_qs, function(bp) {
      d <- ts_df
      d$level_war <- as.numeric(d$time_index >= bp)
      d$trend_war <- ifelse(d$level_war == 1, d$time_index - bp + 1, 0)
      d$level_postwar <- 0; d$trend_postwar <- 0
      m <- suppressWarnings(tryCatch(
        gls(n_events ~ time_index + covid + level_war + trend_war,
            data = d, correlation = corAR1(form = ~time_index), method = "ML"),
        error = function(e) NULL))
      if (is.null(m)) return(NA_real_)
      coef(summary(m))["level_war", "Value"]
    })
    estimates <- estimates[is.finite(estimates)]
    if (length(estimates) < 5) return(NULL)
    p_emp <- mean(abs(estimates) >= abs(actual_war_est))
    save_rds(estimates,
             file.path(output_root, "sensitivity",
                       paste0(tolower(substr(label, 1, 3)),
                              "_placebo_dist.rds")))
    data.frame(Outcome = label, Actual_War_Estimate = actual_war_est,
               N_Placebos = length(estimates),
               Placebo_Mean = mean(estimates), Placebo_SD = sd(estimates),
               Placebo_CI_Low = quantile(estimates, 0.025),
               Placebo_CI_High = quantile(estimates, 0.975),
               Empirical_P = p_emp, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}
actual_war_c <- coef(summary(m_cotri))["level_war", "Value"]
actual_war_t <- coef(summary(m_tpt))["level_war",   "Value"]
placebo_all <- bind_rows(
  placebo_test(cotri_q, actual_war_c, "Cotrimoxazole"),
  placebo_test(tpt_q,   actual_war_t, "TPT"))
if (!is.null(placebo_all) && nrow(placebo_all) > 0) {
  save_csv(placebo_all, file.path(output_root, "sensitivity",
                                  "placebo_test.csv"), row.names = FALSE)
}
# --- 13M. GAM smooth-trend ITS ---------------------------------------------
cat("  13M. GAM smooth-trend ITS\n")
gam_its <- function(ts_df, label) {
  tryCatch({
    m <- mgcv::gam(n_events ~ s(time_index, k = 12, bs = "cr") +
                     covid + level_war + trend_war +
                     level_postwar + trend_postwar,
                   family = mgcv::nb(), data = ts_df, method = "REML")
    co <- summary(m)$p.coeff
    se <- summary(m)$se[names(co)]
    pv <- summary(m)$p.pv[names(co)]
    list(model = m,
         table = data.frame(Outcome = label, Parameter = names(co),
                            Log_Estimate = co, SE = se,
                            IRR = exp(co),
                            IRR_CI_Lower = exp(co - 1.96 * se),
                            IRR_CI_Upper = exp(co + 1.96 * se),
                            P_Value = pv, row.names = NULL,
                            stringsAsFactors = FALSE))
  }, error = function(e) NULL)
}
gam_c <- gam_its(cotri_q, "Cotrimoxazole")
gam_t <- gam_its(tpt_q,   "TPT")
gam_all <- bind_rows(if (!is.null(gam_c)) gam_c$table,
                     if (!is.null(gam_t)) gam_t$table)
if (nrow(gam_all) > 0) {
  save_csv(gam_all, file.path(output_root, "sensitivity",
                              "gam_smooth_trend_its.csv"), row.names = FALSE)
}
# --- 13N. Granger causality -------------------------------------------------
cat("  13N. Granger causality\n")
granger_run <- function(y_ts, x_ts, y_lab, x_lab, lags = 1:4) {
  bind_rows(lapply(lags, function(L) {
    g <- tryCatch(lmtest::grangertest(y_ts ~ x_ts, order = L),
                  error = function(e) NULL)
    if (is.null(g)) return(NULL)
    data.frame(Direction = sprintf("%s -> %s", x_lab, y_lab),
               Lag_Quarters = L, F_Stat = g$F[2], DF1 = g$Df[2],
               DF2 = g$Res.Df[2], P_Value = g$`Pr(>F)`[2],
               row.names = NULL, stringsAsFactors = FALSE)
  }))
}
granger_all <- bind_rows(
  granger_run(cotri_q$n_events, tpt_q$n_events, "CTX", "TPT"),
  granger_run(tpt_q$n_events,   cotri_q$n_events, "TPT", "CTX"))
if (nrow(granger_all) > 0) {
  save_csv(granger_all, file.path(output_root, "sensitivity",
                                  "granger_causality.csv"), row.names = FALSE)
}
# --- 13O. Engle-Granger cointegration ---------------------------------------
cat("  13O. Engle-Granger cointegration test\n")
tryCatch({
  fit_eg <- lm(cotri_q$n_events ~ tpt_q$n_events)
  res_eg <- residuals(fit_eg)
  adf_res <- suppressWarnings(tseries::adf.test(res_eg))
  save_csv(data.frame(
    Step = c("Cointegrating regression", "ADF on residuals (Engle-Granger)"),
    Statistic = c(coef(fit_eg)[2], adf_res$statistic),
    SE_or_Lag = c(summary(fit_eg)$coefficients[2, 2], adf_res$parameter),
    P_or_Pval = c(summary(fit_eg)$coefficients[2, 4], adf_res$p.value),
    row.names = NULL, stringsAsFactors = FALSE),
    file.path(output_root, "sensitivity", "cointegration_eg.csv"),
    row.names = FALSE)
}, error = function(e) cat("    cointegration failed\n"))
# --- 13P. Bai-Perron --------------------------------------------------------
cat("  13P. Bai-Perron multiple structural breakpoints\n")
bp_run <- function(ts_df, label) {
  tryCatch({
    bp <- strucchange::breakpoints(n_events ~ time_index, data = ts_df,
                                   h = 0.10, breaks = 5)
    valid_bps <- if (length(bp$breakpoints) > 0 &&
                     !any(is.na(bp$breakpoints))) bp$breakpoints else integer(0)
    bdates <- if (length(valid_bps) > 0) ts_df$date[valid_bps] else as.Date(NA)
    rss_val <- tryCatch({
      rss_t <- bp$RSS.table
      if (is.null(rss_t)) NA_real_
      else if (is.matrix(rss_t) && "RSS" %in% colnames(rss_t))
        as.numeric(rss_t[nrow(rss_t), "RSS"])
      else if (is.matrix(rss_t) && ncol(rss_t) >= 2)
        as.numeric(rss_t[nrow(rss_t), 2])
      else NA_real_
    }, error = function(e) NA_real_)
    bic_val <- tryCatch(as.numeric(AIC(bp, k = log(nrow(ts_df)))),
                        error = function(e) NA_real_)
    list(table = data.frame(
      Outcome = label, N_Breaks = length(valid_bps),
      Break_Dates = if (length(valid_bps) > 0)
        paste(paste0(year(bdates), " Q", quarter(bdates)), collapse = "; ")
      else "(none detected)",
      RSS = rss_val, BIC = bic_val,
      row.names = NULL, stringsAsFactors = FALSE),
      model = bp)
  }, error = function(e) NULL)
}
bp_c <- bp_run(cotri_q, "Cotrimoxazole")
bp_t <- bp_run(tpt_q,   "TPT")
bp_table <- bind_rows(if (!is.null(bp_c)) bp_c$table,
                      if (!is.null(bp_t)) bp_t$table)
if (nrow(bp_table) > 0) {
  save_csv(bp_table, file.path(output_root, "sensitivity",
                               "baiperron_breakpoints.csv"), row.names = FALSE)
}
# --- 13Q. Rolling-origin tsCV -----------------------------------------------
cat("  13Q. Rolling-origin tsCV\n")
tscv_run <- function(ts_df, label) {
  tryCatch({
    pre_q <- ts_df %>% filter(year <= 2020)
    if (nrow(pre_q) < 20) return(NULL)
    ts_obj <- ts(pre_q$n_events, start = c(2005, 1), frequency = 4)
    fcast_fn <- function(x, h) {
      m <- suppressWarnings(forecast::auto.arima(x, seasonal = TRUE, stepwise = TRUE))
      forecast::forecast(m, h = h)
    }
    e1 <- forecast::tsCV(ts_obj, fcast_fn, h = 1, initial = 12)
    e4 <- forecast::tsCV(ts_obj, fcast_fn, h = 4, initial = 12)
    actual <- as.numeric(ts_obj)
    metrics <- function(e) {
      e <- e[is.finite(e)]
      data.frame(N = length(e), RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
                 MAPE = mean(abs(e) / pmax(actual[1:length(e)], 1)) * 100)
    }
    bind_rows(
      cbind(Outcome = label, Horizon = "h=1 quarter",  metrics(as.numeric(e1))),
      cbind(Outcome = label, Horizon = "h=4 quarters", metrics(as.numeric(e4[, 4]))))
  }, error = function(e) NULL)
}
tscv_all <- bind_rows(tscv_run(cotri_q, "Cotrimoxazole"),
                      tscv_run(tpt_q,   "TPT"))
if (nrow(tscv_all) > 0) {
  save_csv(tscv_all, file.path(output_root, "sensitivity",
                               "tscv_accuracy.csv"), row.names = FALSE)
}
# --- 13R. Zivot-Andrews -----------------------------------------------------
cat("  13R. Zivot-Andrews\n")
za_run <- function(ts_df, label) {
  tryCatch({
    za <- urca::ur.za(ts_df$n_events, model = "both", lag = 4)
    cv <- za@cval
    cv_get <- function(name) {
      if (is.matrix(cv)) {
        if (name %in% colnames(cv)) as.numeric(cv[1, name]) else NA_real_
      } else if (is.numeric(cv) && !is.null(names(cv))) {
        if (name %in% names(cv)) as.numeric(cv[name]) else NA_real_
      } else NA_real_
    }
    cv5 <- cv_get("5pct"); bp_idx_za <- za@bpoint
    data.frame(Outcome = label, ZA_Statistic = as.numeric(za@teststat),
               Break_Index = bp_idx_za,
               Break_Date = if (!is.na(bp_idx_za) && bp_idx_za > 0 &&
                                bp_idx_za <= nrow(ts_df))
                 as.character(ts_df$date[bp_idx_za]) else NA_character_,
               CV_1pct = cv_get("1pct"), CV_5pct = cv5,
               CV_10pct = cv_get("10pct"),
               Reject_unit_root_5pct =
                 if (!is.na(cv5)) ifelse(as.numeric(za@teststat) < cv5, "Yes", "No")
               else NA_character_,
               row.names = NULL, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}
za_all <- bind_rows(za_run(cotri_q, "Cotrimoxazole"),
                    za_run(tpt_q,   "TPT"))
if (nrow(za_all) > 0) {
  save_csv(za_all, file.path(output_root, "sensitivity",
                             "zivot_andrews.csv"), row.names = FALSE)
}
# --- 13S. Quantile ITS ------------------------------------------------------
cat("  13S. Quantile ITS regression\n")
qreg_run <- function(ts_df, label) {
  tryCatch({
    fit <- quantreg::rq(n_events ~ time_index + covid + level_war + trend_war +
                          level_postwar + trend_postwar,
                        tau = c(0.25, 0.5, 0.75), data = ts_df)
    sm <- summary(fit, se = "boot", R = 500)
    bind_rows(lapply(seq_along(sm), function(i) {
      co <- sm[[i]]$coefficients
      data.frame(Outcome = label, Tau = sm[[i]]$tau,
                 Parameter = rownames(co),
                 Estimate = co[, "Value"], SE = co[, "Std. Error"],
                 t_or_z = co[, ncol(co) - 1], P_Value = co[, ncol(co)],
                 CI_Low  = co[, "Value"] - 1.96 * co[, "Std. Error"],
                 CI_High = co[, "Value"] + 1.96 * co[, "Std. Error"],
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
  }, error = function(e) NULL)
}
qreg_all <- bind_rows(qreg_run(cotri_q, "Cotrimoxazole"),
                      qreg_run(tpt_q,   "TPT"))
if (nrow(qreg_all) > 0) {
  save_csv(qreg_all, file.path(output_root, "sensitivity",
                               "quantile_its.csv"), row.names = FALSE)
}
# --- 13T. Distributed-lag NB-ITS --------------------------------------------
cat("  13T. Distributed-lag NB-ITS\n")
distributed_lag_nb <- function(ts_df, label) {
  tryCatch({
    d <- ts_df %>% mutate(
      level_war_lag0 = level_war,
      level_war_lag1 = dplyr::lag(level_war, 1, default = 0),
      level_war_lag2 = dplyr::lag(level_war, 2, default = 0))
    m <- MASS::glm.nb(n_events ~ time_index + covid +
                        level_war_lag0 + level_war_lag1 + level_war_lag2 +
                        trend_war + level_postwar + trend_postwar,
                      data = d,
                      control = glm.control(maxit = 300, epsilon = 1e-8))
    if (!isTRUE(m$converged)) return(NULL)
    co <- coef(summary(m))
    keep <- intersect(c("level_war_lag0", "level_war_lag1", "level_war_lag2"),
                      rownames(co))
    est <- co[keep, "Estimate"]; se <- co[keep, "Std. Error"]
    pv  <- co[keep, "Pr(>|z|)"]
    data.frame(Outcome = label,
               Lag_Quarters = as.integer(sub("level_war_lag", "", keep)),
               Log_Estimate = est, SE = se,
               IRR = exp(est),
               IRR_CI_Lower = exp(est - 1.96 * se),
               IRR_CI_Upper = exp(est + 1.96 * se),
               P_Value = pv,
               row.names = NULL, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}
dl_all <- bind_rows(distributed_lag_nb(cotri_q, "Cotrimoxazole"),
                    distributed_lag_nb(tpt_q,   "TPT"))
if (!is.null(dl_all) && nrow(dl_all) > 0) {
  save_csv(dl_all, file.path(output_root, "sensitivity",
                             "distributed_lag.csv"), row.names = FALSE)
}
# --- 13U. Robust ITS via rlm ------------------------------------------------
cat("  13U. Robust ITS (rlm Huber)\n")
robust_its <- function(ts_df, label) {
  tryCatch({
    m <- MASS::rlm(n_events ~ time_index + covid + level_war + trend_war +
                     level_postwar + trend_postwar,
                   data = ts_df, psi = MASS::psi.huber, maxit = 100)
    co <- coef(summary(m))
    data.frame(Outcome = label, Parameter = rownames(co),
               Estimate = co[, "Value"], SE = co[, "Std. Error"],
               t_value = co[, "t value"],
               CI_Low  = co[, "Value"] - 1.96 * co[, "Std. Error"],
               CI_High = co[, "Value"] + 1.96 * co[, "Std. Error"],
               row.names = NULL, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}
rlm_all <- bind_rows(robust_its(cotri_q, "Cotrimoxazole"),
                     robust_its(tpt_q,   "TPT"))
if (!is.null(rlm_all) && nrow(rlm_all) > 0) {
  save_csv(rlm_all, file.path(output_root, "sensitivity",
                              "robust_its_rlm.csv"), row.names = FALSE)
}
# --- 13V. Mixed-effects NB via glmmTMB --------------------------------------
cat("  13V. Mixed-effects NB-ITS (random year)\n")
mixed_effects_its <- function(ts_df, label) {
  tryCatch({
    m <- glmmTMB::glmmTMB(
      n_events ~ time_index + covid + level_war + trend_war +
        level_postwar + trend_postwar + (1 | year),
      data = ts_df, family = glmmTMB::nbinom2)
    if (!is.null(m$fit$convergence) && m$fit$convergence != 0) return(NULL)
    co <- summary(m)$coefficients$cond
    re_var <- as.numeric(VarCorr(m)$cond$year)
    data.frame(Outcome = label, Parameter = rownames(co),
               Log_Estimate = co[, "Estimate"], SE = co[, "Std. Error"],
               IRR = exp(co[, "Estimate"]),
               IRR_CI_Lower = exp(co[, "Estimate"] - 1.96 * co[, "Std. Error"]),
               IRR_CI_Upper = exp(co[, "Estimate"] + 1.96 * co[, "Std. Error"]),
               Z = co[, "z value"], P_Value = co[, "Pr(>|z|)"],
               Year_RE_Variance = re_var,
               row.names = NULL, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}
me_all <- bind_rows(mixed_effects_its(cotri_q, "Cotrimoxazole"),
                    mixed_effects_its(tpt_q,   "TPT"))
if (!is.null(me_all) && nrow(me_all) > 0) {
  save_csv(me_all, file.path(output_root, "sensitivity",
                             "mixed_effects_its.csv"), row.names = FALSE)
}
# --- 13W. Mann-Kendall + Theil-Sen ------------------------------------------
cat("  13W. Mann-Kendall + Theil-Sen\n")
mk_sen <- function(ts_df, label) {
  tryCatch({
    pre  <- ts_df %>% filter(year <= 2019)
    sets <- list("Pre-war (2005-2019)" = pre$n_events,
                 "Full series (2005-2025 Q2)" = ts_df$n_events)
    bind_rows(lapply(names(sets), function(s) {
      x <- sets[[s]]
      mk <- trend::mk.test(x); ss <- trend::sens.slope(x)
      data.frame(Outcome = label, Window = s, N = length(x),
                 MK_Tau = as.numeric(mk$estimates["tau"]),
                 MK_S = as.numeric(mk$estimates["S"]),
                 MK_P_Value = mk$p.value,
                 Sen_Slope = as.numeric(ss$estimates),
                 Sen_CI_Low = as.numeric(ss$conf.int[1]),
                 Sen_CI_High = as.numeric(ss$conf.int[2]),
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
  }, error = function(e) NULL)
}
mk_all <- bind_rows(mk_sen(cotri_q, "Cotrimoxazole"),
                    mk_sen(tpt_q,   "TPT"))
if (!is.null(mk_all) && nrow(mk_all) > 0) {
  save_csv(mk_all, file.path(output_root, "sensitivity",
                             "mann_kendall_sen.csv"), row.names = FALSE)
}
# --- 13X. Forecast benchmarks (ETS/TBATS/NNAR/ARIMA) ------------------------
cat("  13X. ETS / TBATS / NNAR / ARIMA forecast benchmarks\n")
forecast_benchmarks <- function(ts_df, label) {
  tryCatch({
    pre_q  <- ts_df %>% filter(year <= 2020)
    post_q <- ts_df %>% filter(year >= 2021)
    if (nrow(pre_q) < 24 || nrow(post_q) < 4) return(NULL)
    ts_pre <- ts(pre_q$n_events, start = c(2005, 1), frequency = 4)
    h <- nrow(post_q); actual <- post_q$n_events
    fits <- list(
      ARIMA = function() forecast::auto.arima(ts_pre, seasonal = TRUE, stepwise = TRUE),
      ETS   = function() forecast::ets(ts_pre),
      TBATS = function() forecast::tbats(ts_pre),
      NNAR  = function() forecast::nnetar(ts_pre))
    bind_rows(lapply(names(fits), function(nm) {
      fit <- suppressWarnings(tryCatch(fits[[nm]](), error = function(e) NULL))
      if (is.null(fit)) return(NULL)
      fc <- suppressWarnings(tryCatch(forecast::forecast(fit, h = h),
                                      error = function(e) NULL))
      if (is.null(fc)) return(NULL)
      pred <- as.numeric(fc$mean); e <- actual - pred
      data.frame(Outcome = label, Model = nm,
                 AIC = if ("aic" %in% names(fit) && is.numeric(fit$aic))
                   fit$aic else NA_real_,
                 Forecast_RMSE = sqrt(mean(e^2, na.rm = TRUE)),
                 Forecast_MAE  = mean(abs(e), na.rm = TRUE),
                 Forecast_MAPE_pct = mean(abs(e) / pmax(actual, 1), na.rm = TRUE) * 100,
                 Sum_Predicted = sum(pred), Sum_Observed = sum(actual),
                 Pct_Gap = (sum(actual) - sum(pred)) / sum(pred) * 100,
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
  }, error = function(e) NULL)
}
fb_all <- bind_rows(forecast_benchmarks(cotri_q, "Cotrimoxazole"),
                    forecast_benchmarks(tpt_q,   "TPT"))
if (!is.null(fb_all) && nrow(fb_all) > 0) {
  save_csv(fb_all, file.path(output_root, "sensitivity",
                             "forecast_benchmarks.csv"), row.names = FALSE)
}
# --- 13Y. Bayesian NB-ITS via rstanarm (optional) ---------------------------
cat("  13Y. Bayesian NB-ITS (rstanarm; optional)\n")
bayesian_nb_its <- function(ts_df, label) {
  tryCatch({
    if (!requireNamespace("rstanarm", quietly = TRUE)) {
      cat("    rstanarm not installed; skipped.\n"); return(NULL)
    }
    suppressPackageStartupMessages(library(rstanarm))
    m <- rstanarm::stan_glm.nb(
      n_events ~ time_index + covid + level_war + trend_war +
        level_postwar + trend_postwar,
      data = ts_df, prior = normal(0, 5), prior_intercept = normal(0, 10),
      chains = 4, iter = 2000, seed = 42, refresh = 0)
    post <- as.matrix(m)
    keep <- intersect(c("(Intercept)", "time_index", "covid", "level_war",
                        "trend_war", "level_postwar", "trend_postwar"),
                      colnames(post))
    bind_rows(lapply(keep, function(p) {
      x <- post[, p]
      data.frame(Outcome = label, Parameter = p,
                 Posterior_Mean = mean(x), Posterior_SD = sd(x),
                 IRR_Mean = mean(exp(x)),
                 Cred_2_5 = quantile(x, 0.025), Cred_50 = quantile(x, 0.5),
                 Cred_97_5 = quantile(x, 0.975),
                 IRR_Cred_2_5 = quantile(exp(x), 0.025),
                 IRR_Cred_97_5 = quantile(exp(x), 0.975),
                 Prob_Negative = mean(x < 0),
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
  }, error = function(e) NULL)
}
bay_all <- bind_rows(bayesian_nb_its(cotri_q, "Cotrimoxazole"),
                     bayesian_nb_its(tpt_q,   "TPT"))
if (!is.null(bay_all) && nrow(bay_all) > 0) {
  save_csv(bay_all, file.path(output_root, "sensitivity",
                              "bayesian_nb_its.csv"), row.names = FALSE)
}
cat("  Advanced sensitivity (13J-13Y) complete.\n\n")
# ==============================================================================
# 13Z. TPT REGIMEN-TYPE FILTER SENSITIVITY  (V13 NEW)
# ==============================================================================
tic("13Z: TPT TYPE FILTER (V13) ")
cat("══ SECTION 13Z: TPT TYPE FILTER SENSITIVITY (V13 NEW) ═══════════════════\n")
if ("tb_prophylaxis_type" %in% names(tpt)) {
  # Year × type tabulation
  type_yearly <- tpt %>%
    filter(!is.na(tpt_start_date) | !is.na(registration_date)) %>%
    mutate(yr = year(dplyr::coalesce(tpt_start_date, registration_date)),
           type_label = ifelse(is.na(tb_prophylaxis_type), "NA / blank",
                               as.character(tb_prophylaxis_type))) %>%
    filter(!is.na(yr), yr >= 2005, yr <= 2025) %>%
    count(yr, type_label, name = "n") %>%
    arrange(yr, type_label)
  save_csv(type_yearly, file.path(output_root, "tables",
                                  "table28_tpt_type_by_year.csv"),
           row.names = FALSE)
  save_csv(type_yearly, file.path(output_root, "sensitivity",
                                  "tpt_type_by_year.csv"),
           row.names = FALSE)
  # Rebuild quarterly series WITHOUT the type filter
  tpt_all <- tpt %>%
    mutate(
      tpt_start_date_any = dplyr::coalesce(inhprophylaxis_started_date,
                                           InhprophylaxisDiscontinuedDate,
                                           InhprophylaxisCompletedDate),
      initiated_any = !is.na(tpt_start_date_any) &
        tpt_start_date_any >= config$study_start &
        tpt_start_date_any <= config$study_end)
  build_q_no_filter <- function(d) {
    q <- d %>%
      filter(initiated_any) %>%
      mutate(year = year(tpt_start_date_any),
             quarter = quarter(tpt_start_date_any)) %>%
      count(year, quarter, name = "n_events") %>%
      right_join(full_grid, by = c("year", "quarter")) %>%
      mutate(n_events = replace_na(n_events, 0)) %>% arrange(year, quarter) %>%
      mutate(date = ymd(paste0(year, "-", quarter * 3 - 2, "-01")),
             time_index = row_number(), season = factor(quarter),
             covid = as.numeric(year == 2020))
    bp_w <- which(q$year == 2021 & q$quarter == 1)
    bp_p <- which(q$year == 2023 & q$quarter == 1)
    q %>% mutate(
      level_war = as.numeric(time_index >= bp_w),
      trend_war = ifelse(level_war == 1, time_index - bp_w + 1, 0),
      level_postwar = as.numeric(time_index >= bp_p),
      trend_postwar = ifelse(level_postwar == 1, time_index - bp_p + 1, 0))
  }
  tpt_q_unfiltered <- build_q_no_filter(tpt_all)
  save_csv(tpt_q_unfiltered,
           file.path(output_root, "data",
                     "tpt_quarterly_ts_no_type_filter.csv"),
           row.names = FALSE)
  m_tpt_unf <- suppressWarnings(tryCatch(
    gls(n_events ~ time_index + covid + level_war + trend_war +
          level_postwar + trend_postwar,
        data = tpt_q_unfiltered,
        correlation = corAR1(form = ~time_index), method = "ML"),
    error = function(e) NULL))
  extract_war_gls <- function(m, label) {
    if (is.null(m)) return(NULL)
    co <- coef(summary(m))
    e <- as.numeric(co["level_war", "Value"])
    s <- as.numeric(co["level_war", "Std.Error"])
    ep <- as.numeric(co["level_postwar", "Value"])
    sp <- as.numeric(co["level_postwar", "Std.Error"])
    data.frame(Model = label,
               War_Estimate = e, War_SE = s,
               War_CI_Low = e - 1.96 * s, War_CI_High = e + 1.96 * s,
               Postwar_Estimate = ep, Postwar_SE = sp,
               Postwar_CI_Low = ep - 1.96 * sp,
               Postwar_CI_High = ep + 1.96 * sp,
               stringsAsFactors = FALSE)
  }
  type_compare <- bind_rows(
    extract_war_gls(m_tpt,     "V12 (type == 1 OR NA)"),
    extract_war_gls(m_tpt_unf, "V13 (no type filter)"))
  save_csv(type_compare,
           file.path(output_root, "tables",
                     "table29_tpt_type_filter_war_effects.csv"),
           row.names = FALSE)
  save_csv(type_compare,
           file.path(output_root, "sensitivity",
                     "tpt_type_filter_war_effects.csv"),
           row.names = FALSE)
  print(type_compare, row.names = FALSE)
} else {
  cat("  tb_prophylaxis_type column not found - skipped.\n")
  tpt_q_unfiltered <- NULL
}
# ==============================================================================
# 13AA. UNDER-RECORDING SENSITIVITY SCENARIOS  (V13 NEW)
# ==============================================================================
tic("13AA: UNDER-RECORDING SENSITIVITY (V13) ")
cat("══ SECTION 13AA: UNDER-RECORDING SENSITIVITY (V13 NEW) ══════════════════\n")
inflate_fit <- function(ts_df, label, scenarios = c(0, 0.10, 0.25, 0.50, 0.75),
                        target = c("war", "war_and_postwar")) {
  target <- match.arg(target)
  war_q   <- ts_df$year %in% 2021:2022
  post_q  <- ts_df$year >= 2023
  bind_rows(lapply(scenarios, function(u) {
    d <- ts_df
    if (u > 0) {
      d$n_events[war_q] <- round(d$n_events[war_q] / (1 - u))
      if (target == "war_and_postwar") {
        d$n_events[post_q] <- round(d$n_events[post_q] / (1 - u))
      }
    }
    m <- suppressWarnings(tryCatch(
      gls(n_events ~ time_index + covid + level_war + trend_war +
            level_postwar + trend_postwar,
          data = d, correlation = corAR1(form = ~time_index), method = "ML"),
      error = function(e) NULL))
    if (is.null(m)) return(NULL)
    co <- coef(summary(m))
    e  <- as.numeric(co["level_war", "Value"])
    s  <- as.numeric(co["level_war", "Std.Error"])
    ep <- as.numeric(co["level_postwar", "Value"])
    sp <- as.numeric(co["level_postwar", "Std.Error"])
    data.frame(
      Outcome = label,
      Scenario = paste0(round(u * 100), "% under-recording (",
                        ifelse(target == "war", "war only", "war + post-war"), ")"),
      U = u, Target = target,
      War_Estimate = e, War_SE = s,
      War_CI_Low = e - 1.96 * s, War_CI_High = e + 1.96 * s,
      War_Significant = (e - 1.96 * s > 0) | (e + 1.96 * s < 0),
      Postwar_Estimate = ep, Postwar_CI_Low = ep - 1.96 * sp,
      Postwar_CI_High = ep + 1.96 * sp,
      Sum_War_Obs_Inflated = sum(d$n_events[war_q]),
      Sum_Post_Obs_Inflated = sum(d$n_events[post_q]),
      stringsAsFactors = FALSE)
  }))
}
ur_results <- bind_rows(
  inflate_fit(cotri_q, "Cotrimoxazole", target = "war"),
  inflate_fit(cotri_q, "Cotrimoxazole", target = "war_and_postwar"),
  inflate_fit(tpt_q,   "TPT",           target = "war"),
  inflate_fit(tpt_q,   "TPT",           target = "war_and_postwar")
)
save_csv(ur_results, file.path(output_root, "tables",
                               "table30_under_recording_scenarios.csv"),
         row.names = FALSE)
save_csv(ur_results, file.path(output_root, "sensitivity",
                               "under_recording_scenarios.csv"),
         row.names = FALSE)
threshold_table <- ur_results %>%
  filter(Target == "war") %>%
  group_by(Outcome) %>%
  summarise(
    first_nonsig_U = ifelse(any(!War_Significant),
                            min(U[!War_Significant]) * 100, NA_real_),
    max_tested_U   = max(U) * 100,
    interpretation = ifelse(is.na(first_nonsig_U),
                            paste0("War effect remains significant even at ",
                                   max_tested_U, "% war-only under-recording"),
                            paste0("War effect becomes non-significant at ~",
                                   round(first_nonsig_U), "% under-recording")),
    .groups = "drop")
save_csv(threshold_table, file.path(output_root, "tables",
                                    "table31_under_recording_threshold.csv"),
         row.names = FALSE)
cat("\n  Under-recording threshold:\n")
print(threshold_table, row.names = FALSE)
# ==============================================================================
# 13BB. ART-DENOMINATOR ADJUSTMENT FOR TPT  (V13 NEW)
# ==============================================================================
tic("13BB: ART-DENOMINATOR ADJUSTMENT (V13) ")
cat("══ SECTION 13BB: ART-DENOMINATOR ADJUSTMENT (V13 NEW) ═══════════════════\n")
art_starts_q <- function(d) {
  d %>%
    filter(!is.na(art_start_date),
           art_start_date >= config$study_start,
           art_start_date <= config$study_end) %>%
    distinct(patient_id, art_start_date) %>%   # dedup BEFORE mutate
    mutate(year    = year(art_start_date),
           quarter = quarter(art_start_date)) %>%
    count(year, quarter, name = "n_art_starts")
}
art_q <- art_starts_q(tpt) %>%
  right_join(full_grid, by = c("year", "quarter")) %>%
  mutate(n_art_starts = replace_na(n_art_starts, 0)) %>%
  arrange(year, quarter)
tpt_rate <- tpt_q %>%
  dplyr::select(year, quarter, date, n_events_tpt = n_events,
                time_index, season, covid,
                level_war, trend_war, level_postwar, trend_postwar) %>%
  left_join(art_q, by = c("year", "quarter")) %>%
  mutate(tpt_per_100_art = ifelse(n_art_starts > 0,
                                  100 * n_events_tpt / n_art_starts, NA_real_))
save_csv(tpt_rate, file.path(output_root, "data",
                             "tpt_per_100_art_quarterly.csv"),
         row.names = FALSE)
tpt_rate_fit <- tpt_rate %>% filter(!is.na(tpt_per_100_art),
                                    is.finite(tpt_per_100_art))
m_rate <- NULL
if (nrow(tpt_rate_fit) >= 20 && sum(tpt_rate_fit$level_war) >= 4) {
  m_rate <- suppressWarnings(tryCatch(
    gls(tpt_per_100_art ~ time_index + covid + level_war + trend_war +
          level_postwar + trend_postwar,
        data = tpt_rate_fit,
        correlation = corAR1(form = ~time_index), method = "ML"),
    error = function(e) NULL))
  if (!is.null(m_rate)) {
    co <- coef(summary(m_rate))
    rate_table <- data.frame(
      Parameter = rownames(co),
      Estimate = co[, "Value"], SE = co[, "Std.Error"],
      CI_Low = co[, "Value"] - 1.96 * co[, "Std.Error"],
      CI_High = co[, "Value"] + 1.96 * co[, "Std.Error"],
      P_Value = co[, "p-value"],
      row.names = NULL, stringsAsFactors = FALSE)
    save_csv(rate_table, file.path(output_root, "tables",
                                   "table32_tpt_rate_per_100_art_gls.csv"),
             row.names = FALSE)
    save_csv(rate_table, file.path(output_root, "sensitivity",
                                   "tpt_rate_per_100_art_gls.csv"),
             row.names = FALSE)
    cat("\n  GLS-AR1 on TPT rate per 100 ART starts:\n")
    print(rate_table, row.names = FALSE)
  }
}
# ==============================================================================
# 13CC. YEAR-LEVEL RECORDING-QUALITY SUMMARY  (V13 NEW)
# ==============================================================================
tic("13CC: YEAR RECORDING QUALITY (V13) ")
cat("══ SECTION 13CC: YEAR-LEVEL RECORDING QUALITY (V13 NEW) ═════════════════\n")
year_quality <- function(d, anchor, outcome_label) {
  has_id <- !is.na(d$patient_id) & d$patient_id != "" & d$patient_id != "NA"
  d$has_id <- has_id; d$has_sex <- !is.na(d$sex); d$has_age <- !is.na(d$age)
  d %>%
    filter(!is.na(.data[[anchor]]),
           .data[[anchor]] >= config$study_start,
           .data[[anchor]] <= config$study_end) %>%
    mutate(yr = year(.data[[anchor]])) %>%
    group_by(yr) %>%
    summarise(n_records = n(),
              n_initiated = sum(initiated, na.rm = TRUE),
              pct_initiated = round(100 * mean(initiated, na.rm = TRUE), 1),
              pct_with_id = round(100 * mean(has_id), 1),
              pct_with_sex = round(100 * mean(has_sex), 1),
              pct_with_age = round(100 * mean(has_age), 1),
              n_unique_pat = ifelse(any(has_id),
                                    n_distinct(patient_id[has_id]), 0L),
              .groups = "drop") %>%
    mutate(outcome = outcome_label) %>% arrange(yr)
}
yq_all <- bind_rows(year_quality(cotri, "record_date", "Cotrimoxazole"),
                    year_quality(tpt,   "record_date", "TPT"))
save_csv(yq_all, file.path(output_root, "tables",
                           "table33_year_recording_quality.csv"),
         row.names = FALSE)
cat("\n  Year-level recording quality table saved.\n\n")
# ==============================================================================
# CHECKPOINT - save all objects Part 2 will need
# ==============================================================================
tic("CHECKPOINT")
cat("══ CHECKPOINT: saving objects for Part 2 ════════════════════════════════\n")
checkpoint_file <- file.path(output_root, "v13_part1_checkpoint.RData")
save(
  # config / paths / counters
  data_dir, output_root, fig_main, fig_supp, fig_diag, fig_v13,
  config, flow, full_grid, c_pid,
  # styling
  lancet_colors, period6_colors, lancet_lines, annots,
  # event-level + quarterly
  cotri, tpt, cotri_q, tpt_q,
  cotri_for_t1, tpt_for_t1,
  # primary models
  m_cotri, m_tpt, nb_cotri, nb_tpt, arima_cotri, arima_tpt,
  all_cotri, all_tpt,
  # segmented / survival / subgroups
  sg_c, sg_t, sv_c, sv_t, sub_c, sub_t,
  # decomposition / seasonality
  si_c, si_t, si6_all,
  # V12 sensitivity results
  sa_all, sens_corr, bp_sens, covid_all, outlier_results,
  hac_all, preval, loqo_all,
  ci_cotri, ci_tpt, bs_all, placebo_all,
  gam_c, gam_t, gam_all, granger_all, bp_c, bp_t, bp_table,
  tscv_all, za_all, qreg_all, dl_all, rlm_all, me_all, mk_all,
  fb_all, bay_all,
  # diagnostics
  diag_c_gls, diag_t_gls, diag_c_nb, diag_t_nb, adv_c, adv_t,
  chow_results, impact,
  # V13 new results
  rt_c, rt_t, rt_all, rt_period, ds_all,
  tpt_q_unfiltered, type_compare, type_yearly, m_tpt_unf,
  ur_results, threshold_table,
  tpt_rate, m_rate, yq_all,
  # actual war effects for placebo plots
  actual_war_c, actual_war_t,
  # I/O helpers
  save_csv, save_rds, save_fig, save_fig_pub, theme_lancet_pub,
  bp_idx, tic, .script_t0,
  file = checkpoint_file
)
cat(sprintf("    saved: %s  (%.1f MB)\n",
            normalizePath(checkpoint_file, mustWork = FALSE),
            file.info(checkpoint_file)$size / 1024^2))
cat("\n############################################################\n")
cat("#  V13 PART 1 complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("#  Now run Part 2 to build tables, figures, and the master\n")
cat("#  workbook.\n")
cat("############################################################\n")

################################################################################
#
#  ITS V13 - PART 2 of 2: Tables, figures, session info, master workbook
#  Cotrimoxazole prophylaxis & TB preventive therapy during the Tigray conflict
#
#  Nature Communications | R 4.6.0 | CRAN-only
#
#  HOW TO RUN
#  ----------
#  Either source after Part 1 in the same session:
#      source("ITS_V13_Part1_Analysis.R")
#      source("ITS_V13_Part2_Tables_Figures.R")
#  Or run Part 2 standalone - it loads v13_part1_checkpoint.RData automatically.
#
#  FIGURE STANDARDS
#  ----------------
#  Every figure dpi = 600, bg = "white", units = "mm" (enforced by save_fig /
#  save_fig_pub). Main figures (Fig0-Fig3b, FigV13_1-8) additionally export
#  as TIFF (LZW). Fig0 also exports PDF via cairo.
#
#  CONTENTS
#  --------
#    Section 14   Tables 1-33 (consolidation; Part 1 saved the CSVs)
#    Section 15   Figures: Fig0 (flow), Fig1-Fig3b (main), FigS1-FigS19,
#                  FigV13_1 ... FigV13_8 (V13 recording / under-recording / rate)
#    Section 16   sessionInfo()
#    Section 17   File manifest
#    Section 18   Master workbook ITS_V13_MASTER_RESULTS.xlsx
#
################################################################################
# ==============================================================================
# 14A. CHECKPOINT LOAD (restartable in fresh R session)
# ==============================================================================
if (!exists("cotri_q") || !exists("tpt_q") || !exists("output_root")) {
  candidates <- c(
    file.path("D:/PhD project Data/New data set/FINAL_ANALYSIS_ALL",
              "HIV_care _outcomes/TB/incidence/tpt/prophylactic",
              "Final_ITS_V13_Results", "v13_part1_checkpoint.RData"),
    file.path(getwd(), "Final_ITS_V13_Results", "v13_part1_checkpoint.RData"),
    file.path(getwd(), "v13_part1_checkpoint.RData"))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit))
    stop("Could not find v13_part1_checkpoint.RData. Run Part 1 first.")
  cat("Loading checkpoint:", hit, "\n"); flush.console()
  load(hit, envir = globalenv())
}
# Ensure dplyr verbs are not masked
suppressPackageStartupMessages({
  for (.p in c("MASS", "plyr", "raster", "Hmisc")) {
    pos <- paste0("package:", .p)
    if (pos %in% search()) try(detach(pos, character.only = TRUE,
                                      unload = TRUE, force = TRUE),
                               silent = TRUE)
  }
  library(tidyverse); library(lubridate); library(scales); library(patchwork)
  library(ggtext); library(ggrepel); library(survminer); library(openxlsx)
})
cat("\n############################################################\n")
cat("#  ITS V13 PART 2 - Tables, figures, master workbook       #\n")
cat("############################################################\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
flush.console()
# Re-establish convenient folder shorthands
fig_main <- file.path(output_root, "figures/main")
fig_supp <- file.path(output_root, "figures/supplement")
fig_diag <- file.path(output_root, "figures/diagnostics")
fig_v13  <- file.path(output_root, "figures/v13")
fig_ctx  <- file.path(output_root, "figures/cotri")
fig_tpt  <- file.path(output_root, "figures/tpt")
for (d in c(fig_main, fig_supp, fig_diag, fig_v13, fig_ctx, fig_tpt))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
# ==============================================================================
# 14. TABLES (consolidation of Part 1 CSV outputs)
# ==============================================================================
cat("══ SECTION 14: TABLE CONSOLIDATION ══════════════════════════════════════\n")
tables_dir <- file.path(output_root, "tables")
# Tables 1-26 from V12 (CSV files Part 1 already saved):
#   table1_cotri.csv / table1_tpt.csv             (Section 4)
#   impact_estimates.csv                          (Section 8)
#   season_war_interaction.csv                    (Section 6B)
# plus 23 sensitivity tables in /sensitivity/.
# V13 NEW tables 27-33 already saved by Part 1:
#   table27_recording_quality_by_period.csv       (Section 4B)
#   table28_tpt_type_by_year.csv                  (Section 13Z)
#   table29_tpt_type_filter_war_effects.csv       (Section 13Z)
#   table30_under_recording_scenarios.csv         (Section 13AA)
#   table31_under_recording_threshold.csv         (Section 13AA)
#   table32_tpt_rate_per_100_art_gls.csv          (Section 13BB)
#   table33_year_recording_quality.csv            (Section 13CC)
# Build a master "Table X" consolidated summary CSV combining war/post-war
# effects from primary GLS-AR1 + NB across both outcomes.
build_primary_summary <- function() {
  rows <- list()
  for (lab in c("Cotrimoxazole", "TPT")) {
    m_gls <- if (lab == "Cotrimoxazole") m_cotri else m_tpt
    m_nb  <- if (lab == "Cotrimoxazole") nb_cotri else nb_tpt
    if (!is.null(m_gls)) {
      co <- coef(summary(m_gls))
      e_w <- co["level_war",     "Value"];     s_w <- co["level_war",     "Std.Error"]
      e_p <- co["level_postwar", "Value"];     s_p <- co["level_postwar", "Std.Error"]
      rows[[length(rows) + 1]] <- data.frame(
        Outcome = lab, Model = "GLS-AR1 (primary)",
        Scale = "Events per quarter",
        War_Estimate = e_w, War_CI_Low = e_w - 1.96 * s_w,
        War_CI_High = e_w + 1.96 * s_w, War_P = co["level_war", "p-value"],
        Postwar_Estimate = e_p, Postwar_CI_Low = e_p - 1.96 * s_p,
        Postwar_CI_High = e_p + 1.96 * s_p, Postwar_P = co["level_postwar", "p-value"])
    }
    if (!is.null(m_nb)) {
      co <- coef(summary(m_nb))
      e_w <- co["level_war",     "Estimate"]; s_w <- co["level_war",     "Std. Error"]
      e_p <- co["level_postwar", "Estimate"]; s_p <- co["level_postwar", "Std. Error"]
      rows[[length(rows) + 1]] <- data.frame(
        Outcome = lab, Model = "Negative Binomial",
        Scale = "Log-IRR",
        War_Estimate = e_w, War_CI_Low = e_w - 1.96 * s_w,
        War_CI_High = e_w + 1.96 * s_w, War_P = co["level_war", "Pr(>|z|)"],
        Postwar_Estimate = e_p, Postwar_CI_Low = e_p - 1.96 * s_p,
        Postwar_CI_High = e_p + 1.96 * s_p, Postwar_P = co["level_postwar", "Pr(>|z|)"])
    }
  }
  bind_rows(rows)
}
primary_summary <- build_primary_summary()
save_csv(primary_summary,
         file.path(tables_dir, "table2_primary_war_postwar_effects.csv"),
         row.names = FALSE)
# Subgroup summary table
if (exists("sub_c") && exists("sub_t")) {
  sub_all <- bind_rows(sub_c, sub_t) %>%
    mutate(
      War_Effect = sprintf("%.2f (%.2f to %.2f)",
                           war_estimate, war_ci_lower, war_ci_upper),
      Postwar_Effect = sprintf("%.2f (%.2f to %.2f)",
                               postwar_estimate, postwar_ci_lower, postwar_ci_upper)) %>%
    dplyr::select(outcome, subgroup, n_quarters, total_events, model_type,
                  War_Effect, Postwar_Effect)
  save_csv(sub_all, file.path(tables_dir, "table3_subgroup_effects.csv"),
           row.names = FALSE)
}
# Model comparison summary
mc_c <- tryCatch(read.csv(file.path(output_root, "models",
                                    "cotri_model_comparison.csv")),
                 error = function(e) NULL)
mc_t <- tryCatch(read.csv(file.path(output_root, "models",
                                    "tpt_model_comparison.csv")),
                 error = function(e) NULL)
if (!is.null(mc_c) && !is.null(mc_t)) {
  save_csv(bind_rows(mc_c, mc_t),
           file.path(tables_dir, "table4_model_comparison.csv"),
           row.names = FALSE)
}
cat("  All consolidated tables saved.\n\n")
# ==============================================================================
# 15. FIGURES
# ==============================================================================
cat("══ SECTION 15: FIGURES ══════════════════════════════════════════════════\n")
# ---- FIG 0: Study flow diagram (CONSORT-style) -----------------------------
cat("  Fig0: Study flow diagram (CONSORT-style)\n")
make_flow_panel <- function(title_text, n_raw, n_dedup_removed,
                            n_after_dedup, n_not_init, n_blank,
                            n_patdate_removed, n_initiated, n_patients,
                            accent_color, badge_color) {
  fmt <- function(x) format(x, big.mark = ",")
  stages <- data.frame(
    stage = 1:6,
    y = seq(6.5, 0.5, length.out = 6),
    text = c(
      paste0("**Source records**<br><span style='font-size:9pt;'>",
             fmt(n_raw), " rows</span>"),
      paste0("**After de-duplication**<br><span style='font-size:9pt;'>",
             fmt(n_after_dedup), " rows<br>(removed: ",
             fmt(n_dedup_removed), " perfect duplicates)</span>"),
      paste0("**Events meeting initiation criteria**<br><span style='font-size:9pt;'>",
             fmt(n_raw - n_not_init), " of ", fmt(n_raw), "<br>",
             "(excluded: ", fmt(n_not_init), " not initiated /<br>",
             "invalid date)</span>"),
      paste0("**Records with patient ID**<br><span style='font-size:9pt;'>",
             fmt(n_initiated + n_patdate_removed), " linkable +<br>",
             fmt(n_blank), " unlinked records</span>"),
      paste0("**After patient × date de-duplication**<br><span style='font-size:9pt;'>",
             fmt(n_initiated), " unique events<br>(removed: ",
             fmt(n_patdate_removed), " same-patient/same-date)</span>"),
      paste0("**Analytic sample**<br><span style='font-size:11pt; color:",
             accent_color, ";'>**", fmt(n_initiated), " initiation events**</span><br>",
             "<span style='font-size:9pt;'>(", fmt(n_patients),
             " unique patients)</span>"))
  )
  ggplot(stages, aes(x = 0.5, y = y)) +
    geom_rect(aes(xmin = 0.06, xmax = 0.94,
                  ymin = y - 0.42, ymax = y + 0.42),
              fill = "white", color = accent_color, linewidth = 0.6) +
    geom_richtext(aes(label = text), fill = NA, label.color = NA,
                  size = 3.0, lineheight = 1.15) +
    geom_richtext(aes(x = 0.03, label = paste0("<span style='color:white;'>**",
                                               stage, "**</span>")),
                  fill = badge_color, label.color = NA, label.r = unit(0.18, "lines"),
                  label.padding = unit(c(0.10, 0.20, 0.10, 0.20), "lines"),
                  size = 3.2) +
    geom_segment(data = stages[-nrow(stages), ],
                 aes(x = 0.5, xend = 0.5,
                     y = y - 0.42, yend = stages$y[-1][seq_len(nrow(stages) - 1)] + 0.42),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 color = "grey55", linewidth = 0.45) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 7), expand = c(0, 0)) +
    labs(title = title_text) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 11,
                                    color = accent_color, hjust = 0.5,
                                    margin = margin(b = 6)),
          plot.margin = margin(4, 4, 4, 4))
}
# Null/NA-coalesce helper (must be defined BEFORE the flow-number block)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a
# Build flow numbers from the `flow` list saved in Part 1 (best-effort)
fctx_raw         <- flow$ctx_combined         %||% NA_integer_
fctx_dedup_rm    <- flow$ctx_dedup_removed    %||% 0L
fctx_after_dedup <- flow$ctx_after_merge_dedup %||% (fctx_raw - fctx_dedup_rm)
fctx_not_init    <- flow$ctx_n_not_initiated  %||% 0L
fctx_blank       <- flow$ctx_blank_count      %||% 0L
fctx_patdate_rm  <- flow$ctx_patdate_removed  %||% 0L
fctx_init        <- flow$ctx_n_initiated      %||% sum(cotri_q$n_events)
fctx_pat         <- flow$ctx_patients         %||% NA_integer_
ftpt_raw         <- flow$tpt_raw              %||% NA_integer_
ftpt_dedup_rm    <- flow$tpt_dedup_removed    %||% 0L
ftpt_after_dedup <- flow$tpt_after_merge_dedup %||% (ftpt_raw - ftpt_dedup_rm)
ftpt_not_init    <- flow$tpt_n_not_initiated  %||% 0L
ftpt_blank       <- flow$tpt_blank_count      %||% 0L
ftpt_patdate_rm  <- flow$tpt_patdate_removed  %||% 0L
ftpt_init        <- flow$tpt_n_initiated      %||% sum(tpt_q$n_events)
ftpt_pat         <- flow$tpt_patients         %||% NA_integer_
p_ctx <- make_flow_panel(
  "(A) Cotrimoxazole prophylaxis",
  fctx_raw, fctx_dedup_rm, fctx_after_dedup, fctx_not_init,
  fctx_blank, fctx_patdate_rm, fctx_init, fctx_pat,
  accent_color = "#00468B", badge_color = "#00468B")
p_tpt <- make_flow_panel(
  "(B) TB preventive therapy (TPT)",
  ftpt_raw, ftpt_dedup_rm, ftpt_after_dedup, ftpt_not_init,
  ftpt_blank, ftpt_patdate_rm, ftpt_init, ftpt_pat,
  accent_color = "#AD002A", badge_color = "#AD002A")
# Bottom analytical framework panel
p_frame <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.05, ymax = 0.95,
           fill = "grey97", color = "grey70", linewidth = 0.3) +
  annotate("richtext",
           x = 0.5, y = 0.55,
           label = paste0(
             "**(C) Analytical framework**<br>",
             "<span style='font-size:9pt;'>",
             "Quarterly counts (Q1 2005 - Q2 2025) &middot; GLS-AR(1) primary &middot; ",
             "Negative-binomial sensitivity &middot; Capped-trend counterfactual<br>",
             "Pre-war &middot; <span style='color:#925E9F;'>COVID-19 (2020)</span> &middot; ",
             "<span style='color:#AD002A;'>**War (2021-2022)**</span> &middot; ",
             "<span style='color:#42B540;'>Post-war (2023-2025 Q2)</span></span>"),
           fill = NA, label.color = NA, size = 3.1, hjust = 0.5) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(2, 4, 2, 4))
fig0 <- (p_ctx | p_tpt) / p_frame + plot_layout(heights = c(1.55, 0.22))
save_fig_pub(file.path(fig_main, "Fig0_study_flow"), fig0,
             width = 200, height = 240, also_pdf = TRUE)
# ---- FIG 1: Raw quarterly time series --------------------------------------
cat("  Fig1: Quarterly time series\n")
build_ts_panel <- function(ts_df, label, color, ylab) {
  ggplot(ts_df, aes(x = date, y = n_events)) +
    annots +
    geom_line(color = color, linewidth = 0.65) +
    geom_point(color = color, size = 1.6) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = ylab, subtitle = label) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 11,
                                       hjust = 0, color = color))
}
p1_ctx <- build_ts_panel(cotri_q, "(A) Cotrimoxazole",
                         "#00468BFF", "Initiations per quarter")
p1_tpt <- build_ts_panel(tpt_q, "(B) TPT",
                         "#AD002AFF", "Initiations per quarter")
fig1 <- p1_ctx / p1_tpt + plot_layout(heights = c(1, 1))
save_fig_pub(file.path(fig_main, "Fig1_quarterly_time_series"),
             fig1, width = 180, height = 160)
# ---- FIG 2: ITS with fitted & counterfactual -------------------------------
cat("  Fig2: ITS counterfactual\n")
build_its_panel <- function(ts_df, label, color) {
  plot_df <- ts_df %>%
    transmute(date, n_events,
              fitted_gls = fitted_gls,
              counterfactual_gls = counterfactual_gls,
              counterfactual_capped = counterfactual_capped)
  ggplot(plot_df, aes(x = date)) +
    annots +
    geom_line(aes(y = n_events), color = "grey45", linewidth = 0.4) +
    geom_point(aes(y = n_events), color = "grey45", size = 1.2, alpha = 0.6) +
    geom_line(aes(y = fitted_gls, color = "Fitted (GLS-AR1)"), linewidth = 0.8) +
    geom_line(aes(y = counterfactual_gls, color = "Counterfactual (GLS)"),
              linetype = "dashed", linewidth = 0.8) +
    geom_line(aes(y = counterfactual_capped, color = "Counterfactual (capped)"),
              linetype = "dotted", linewidth = 0.7,
              data = subset(plot_df, !is.na(counterfactual_capped))) +
    scale_color_manual(name = "",
                       values = c("Fitted (GLS-AR1)" = "#AD002AFF",
                                  "Counterfactual (GLS)" = "#00468BFF",
                                  "Counterfactual (capped)" = "#42B540FF")) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Initiations per quarter", subtitle = label) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 11,
                                       hjust = 0, color = color))
}
p2_ctx <- build_its_panel(cotri_q, "(A) Cotrimoxazole", "#00468BFF")
p2_tpt <- build_its_panel(tpt_q,   "(B) TPT",            "#AD002AFF")
fig2 <- p2_ctx / p2_tpt + plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")
save_fig_pub(file.path(fig_main, "Fig2_its_counterfactual"),
             fig2, width = 180, height = 170)
# ---- FIG 3: Subgroup forest (GLS-AR1, scale-consistent) --------------------
cat("  Fig3: Subgroup forest\n")
# To enforce scale consistency, refit every subgroup with GLS-AR1 specifically
# (not NB fallback). This means re-running subgroup ITS with GLS forced.
refit_gls_subgroup <- function(event_data, date_col, quarter_grid, label) {
  sub_ts <- event_data %>%
    mutate(yr  = year(.data[[date_col]]),
           qtr = quarter(.data[[date_col]])) %>%
    count(yr, qtr, name = "n_events") %>%
    right_join(quarter_grid, by = c("yr" = "year", "qtr" = "quarter")) %>%
    mutate(n_events = replace_na(n_events, 0)) %>% arrange(yr, qtr) %>%
    mutate(date = ymd(paste0(yr, "-", qtr * 3 - 2, "-01")),
           time_index = row_number(),
           covid = as.numeric(yr == 2020),
           level_war = as.numeric(date >= config$war_start),
           trend_war = ifelse(level_war == 1,
                              as.numeric(difftime(date, config$war_start,
                                                  units = "days")) / 91.25, 0),
           level_postwar = as.numeric(date >= config$postwar_start),
           trend_postwar = ifelse(level_postwar == 1,
                                  as.numeric(difftime(date, config$postwar_start,
                                                      units = "days")) / 91.25, 0))
  if (nrow(sub_ts) < 20 || sum(sub_ts$level_war == 1) < 4) return(NULL)
  m <- suppressWarnings(tryCatch(
    gls(n_events ~ time_index + covid + level_war + trend_war +
          level_postwar + trend_postwar,
        data = sub_ts, correlation = corAR1(form = ~time_index), method = "ML"),
    error = function(e) NULL))
  if (is.null(m)) return(NULL)
  co <- coef(summary(m))
  if (!"level_war" %in% rownames(co)) return(NULL)
  ew <- co["level_war",     "Value"];     sw <- co["level_war",     "Std.Error"]
  ep <- co["level_postwar", "Value"];     sp <- co["level_postwar", "Std.Error"]
  data.frame(subgroup = label, war_estimate = ew, war_se = sw,
             war_ci_lower = ew - 1.96 * sw, war_ci_upper = ew + 1.96 * sw,
             postwar_estimate = ep, postwar_se = sp,
             postwar_ci_lower = ep - 1.96 * sp,
             postwar_ci_upper = ep + 1.96 * sp,
             n_quarters = nrow(sub_ts), total_events = sum(sub_ts$n_events))
}
build_gls_forest_data <- function(event_data, date_col, label) {
  event_data <- event_data %>% filter(initiated == TRUE)
  res <- list()
  for (s in c("Male", "Female")) {
    r <- refit_gls_subgroup(event_data %>% filter(sex == s),
                            date_col, full_grid, paste("Sex:", s))
    if (!is.null(r)) res[[length(res) + 1]] <- r
  }
  for (ac in levels(event_data$age_cat)) {
    d <- event_data %>% filter(age_cat == ac)
    if (nrow(d) >= 100) {
      r <- refit_gls_subgroup(d, date_col, full_grid, paste("Age:", ac))
      if (!is.null(r)) res[[length(res) + 1]] <- r
    }
  }
  if (!all(is.na(event_data$facility_level))) {
    for (fl in unique(na.omit(event_data$facility_level))) {
      d <- event_data %>% filter(facility_level == fl)
      if (nrow(d) >= 100) {
        r <- refit_gls_subgroup(d, date_col, full_grid, paste("Facility:", fl))
        if (!is.null(r)) res[[length(res) + 1]] <- r
      }
    }
  }
  if (!all(is.na(event_data$facility_ownership))) {
    for (fo in unique(na.omit(event_data$facility_ownership))) {
      d <- event_data %>% filter(facility_ownership == fo)
      if (nrow(d) >= 100) {
        r <- refit_gls_subgroup(d, date_col, full_grid, paste("Ownership:", fo))
        if (!is.null(r)) res[[length(res) + 1]] <- r
      }
    }
  }
  if (length(res) == 0) return(NULL)
  bind_rows(res) %>% mutate(outcome = label)
}
forest_c <- build_gls_forest_data(cotri, "cotri_start_date", "Cotrimoxazole")
forest_t <- build_gls_forest_data(tpt,   "tpt_start_date",   "TPT")
# Also add Overall row from the primary model
make_overall_row <- function(m, label) {
  co <- coef(summary(m))
  ew <- co["level_war",     "Value"];     sw <- co["level_war",     "Std.Error"]
  ep <- co["level_postwar", "Value"];     sp <- co["level_postwar", "Std.Error"]
  data.frame(subgroup = "Overall",
             war_estimate = ew, war_se = sw,
             war_ci_lower = ew - 1.96 * sw, war_ci_upper = ew + 1.96 * sw,
             postwar_estimate = ep, postwar_se = sp,
             postwar_ci_lower = ep - 1.96 * sp,
             postwar_ci_upper = ep + 1.96 * sp,
             n_quarters = NA, total_events = NA,
             outcome = label)
}
forest_all <- bind_rows(
  make_overall_row(m_cotri, "Cotrimoxazole"),
  if (!is.null(forest_c)) forest_c,
  make_overall_row(m_tpt, "TPT"),
  if (!is.null(forest_t)) forest_t
) %>% mutate(sig = (war_ci_lower > 0) | (war_ci_upper < 0))
save_csv(forest_all, file.path(tables_dir,
                               "table5_forest_subgroup_gls.csv"),
         row.names = FALSE)
# Build the forest plot
make_forest <- function(d, outcome_label, color_main) {
  d <- d %>% filter(outcome == outcome_label) %>%
    mutate(subgroup = factor(subgroup, levels = rev(subgroup)),
           sig_color = ifelse(sig, "#AD002A", "#00468B"))
  xrng <- range(c(d$war_ci_lower, d$war_ci_upper), na.rm = TRUE)
  pad  <- diff(xrng) * 0.30
  ggplot(d, aes(x = war_estimate, y = subgroup, color = sig_color)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = war_ci_lower, xmax = war_ci_upper),
                  orientation = "y", width = 0.20, linewidth = 0.55) +
    geom_point(size = 2.4, shape = 18) +
    geom_text(aes(x = xrng[2] + pad,
                  label = sprintf("%.1f (%.1f, %.1f)",
                                  war_estimate, war_ci_lower, war_ci_upper)),
              hjust = 0, size = 2.7, color = "grey20") +
    scale_color_identity() +
    scale_x_continuous(limits = c(xrng[1] - pad * 0.2, xrng[2] + pad * 3),
                       expand = c(0, 0)) +
    labs(x = "War-onset effect (GLS-AR1; events per quarter)", y = NULL,
         subtitle = outcome_label) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 11,
                                       color = color_main, hjust = 0),
          panel.grid.major.y = element_blank())
}
fig3 <- make_forest(forest_all, "Cotrimoxazole", "#00468BFF") /
  make_forest(forest_all, "TPT",            "#AD002AFF")
save_fig_pub(file.path(fig_main, "Fig3_subgroup_forest"),
             fig3, width = 180, height = 200)
# ---- FIG 3b: Segmented regression colored by period6 -----------------------
cat("  Fig3b: Segmented regression\n")
make_seg_panel <- function(ts_df, label, color) {
  ts_df$period6 <- factor(case_when(
    ts_df$year <= 2009 ~ "Early (2005-09)",
    ts_df$year <= 2014 ~ "Scale-up (2010-14)",
    ts_df$year <= 2019 ~ "Plateau (2015-19)",
    ts_df$year == 2020 ~ "COVID-19 (2020)",
    ts_df$year <= 2022 ~ "War (2021-22)",
    TRUE               ~ "Post-war (2023-25)"),
    levels = names(period6_colors))
  ggplot(ts_df, aes(x = date)) +
    annots +
    geom_point(aes(y = n_events, color = period6), size = 1.6) +
    geom_line(aes(y = seg_fitted), color = color, linewidth = 0.7) +
    scale_color_manual(values = period6_colors, name = NULL) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Initiations per quarter", subtitle = label) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 11,
                                       color = color, hjust = 0),
          legend.position = "bottom")
}
fig3b <- make_seg_panel(cotri_q, "(A) Cotrimoxazole", "#00468BFF") /
  make_seg_panel(tpt_q,   "(B) TPT",            "#AD002AFF") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig_pub(file.path(fig_main, "Fig3b_segmented_regression"),
             fig3b, width = 180, height = 200)
# ---- FigS1: STL decomposition ----------------------------------------------
cat("  FigS1: STL decomposition\n")
build_stl_panel <- function(ts_df, label, color) {
  stl_long <- ts_df %>%
    transmute(date,
              Observed   = n_events,
              Trend      = trend,
              Seasonal   = seasonal,
              Remainder  = remainder) %>%
    pivot_longer(-date, names_to = "component", values_to = "value") %>%
    mutate(component = factor(component,
                              levels = c("Observed", "Trend",
                                         "Seasonal", "Remainder")))
  ggplot(stl_long, aes(x = date, y = value)) +
    geom_line(color = color, linewidth = 0.5) +
    facet_wrap(~ component, ncol = 1, scales = "free_y") +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    labs(x = NULL, y = NULL, subtitle = label) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", color = color,
                                       size = 11, hjust = 0))
}
figS1 <- build_stl_panel(cotri_q, "(A) Cotrimoxazole", "#00468BFF") |
  build_stl_panel(tpt_q,   "(B) TPT",            "#AD002AFF")
save_fig(file.path(fig_supp, "FigS1_stl_decomposition.png"),
         figS1, width = 200, height = 220)
# ---- FigS2: Seasonal pattern by period6 ------------------------------------
cat("  FigS2: Seasonality by period\n")
if (exists("si6_all") && nrow(si6_all) > 0) {
  figS2 <- ggplot(si6_all, aes(x = factor(quarter), y = mean_events,
                               color = period6, group = period6)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                  linewidth = 0.4) +
    geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
    facet_wrap(~ outcome, scales = "free_y", ncol = 2) +
    scale_color_manual(values = period6_colors, name = NULL) +
    labs(x = "Quarter", y = "Mean initiations") +
    theme_lancet_pub() + theme(legend.position = "bottom")
  save_fig(file.path(fig_supp, "FigS2_seasonality_by_period6.png"),
           figS2, width = 200, height = 130)
}
# ---- FigS3: Residual diagnostics -------------------------------------------
cat("  FigS3: residual diagnostics (PNG files already in /diagnostics/)\n")
# Already saved by Section 9 in Part 1 as multi-panel PNG.
# ---- FigS4: Kaplan-Meier survival ------------------------------------------
cat("  FigS4: Kaplan-Meier survival (PNG already in /survival/)\n")
# Already saved by Section 11 in Part 1.
# ---- FigS5: Bootstrap impact distributions ---------------------------------
cat("  FigS5: Bootstrap impact CI summary\n")
if (exists("bs_all") && !is.null(bs_all) && nrow(bs_all) > 0) {
  figS5 <- ggplot(bs_all, aes(x = Period, y = Missed_Mean,
                              color = Outcome, group = Outcome)) +
    geom_errorbar(aes(ymin = Missed_CI_Low, ymax = Missed_CI_High),
                  width = 0.15, linewidth = 0.6,
                  position = position_dodge(width = 0.3)) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    scale_color_manual(values = c("Cotrimoxazole" = "#00468BFF",
                                  "TPT"           = "#AD002AFF")) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Missed initiations (block-bootstrap mean & 95% CI)",
         color = NULL) +
    theme_lancet_pub() + theme(legend.position = "bottom")
  save_fig(file.path(fig_supp, "FigS5_bootstrap_impact.png"),
           figS5, width = 180, height = 130)
}
# ---- FigS6: Placebo distribution -------------------------------------------
cat("  FigS6: Placebo / permutation null\n")
make_placebo_panel <- function(prefix, actual, label, color) {
  rds_path <- file.path(output_root, "sensitivity",
                        paste0(prefix, "_placebo_dist.rds"))
  if (!file.exists(rds_path)) return(NULL)
  d <- data.frame(est = readRDS(rds_path))
  ggplot(d, aes(x = est)) +
    geom_histogram(bins = 30, fill = color, color = "white", alpha = 0.7) +
    geom_vline(xintercept = actual, color = "#AD002A", linetype = "dashed",
               linewidth = 0.8) +
    annotate("text", x = actual, y = Inf, vjust = 1.5, hjust = -0.05,
             label = sprintf("Actual = %.1f", actual),
             color = "#AD002A", size = 3) +
    labs(x = "Placebo war-effect estimate", y = "Frequency", subtitle = label) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", color = color,
                                       size = 11, hjust = 0))
}
plc_c <- make_placebo_panel("cot", actual_war_c, "(A) Cotrimoxazole", "#00468BFF")
plc_t <- make_placebo_panel("tpt", actual_war_t, "(B) TPT",            "#AD002AFF")
if (!is.null(plc_c) && !is.null(plc_t)) {
  figS6 <- plc_c | plc_t
  save_fig(file.path(fig_supp, "FigS6_placebo_distribution.png"),
           figS6, width = 200, height = 100)
}
# ---- FigS7: LOQO sensitivity -----------------------------------------------
cat("  FigS7: LOQO\n")
if (exists("loqo_all") && nrow(loqo_all) > 0) {
  figS7 <- ggplot(loqo_all, aes(x = as.Date(Quarter_Removed), y = War_Estimate,
                                color = Outcome)) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Cotrimoxazole" = "#00468BFF",
                                  "TPT"           = "#AD002AFF")) +
    facet_wrap(~ Outcome, ncol = 1, scales = "free_y") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(x = "Quarter removed", y = "War-onset estimate (re-fit)",
         color = NULL) +
    theme_lancet_pub() + theme(legend.position = "none")
  save_fig(file.path(fig_supp, "FigS7_loqo_robustness.png"),
           figS7, width = 180, height = 140)
}
# ---- FigS8: ARIMA-order sensitivity heatmap --------------------------------
cat("  FigS8: ARIMA-order sensitivity\n")
if (exists("sa_all") && nrow(sa_all) > 0) {
  figS8 <- ggplot(sa_all, aes(x = Order, y = Outcome, fill = War_Effect)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f", War_Effect)),
              color = "grey20", size = 3) +
    scale_fill_gradient2(low = "#00468B", mid = "white", high = "#AD002A",
                         midpoint = 0, name = "War effect") +
    labs(x = "ARIMA order (p,d,q)", y = NULL) +
    theme_lancet_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(file.path(fig_supp, "FigS8_arima_orders.png"),
           figS8, width = 180, height = 90)
}
# ---- FigS9: GAM smooth ITS -------------------------------------------------
cat("  FigS9: GAM smooth ITS\n")
if (!is.null(gam_c) && !is.null(gam_t)) {
  build_gam_panel <- function(gres, ts_df, label, color) {
    ts_df$gam_fit <- predict(gres$model, type = "response")
    ggplot(ts_df, aes(x = date)) +
      annots +
      geom_point(aes(y = n_events), color = "grey45", size = 1.2, alpha = 0.6) +
      geom_line(aes(y = gam_fit), color = color, linewidth = 0.8) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      scale_y_continuous(labels = scales::comma) +
      labs(x = NULL, y = "Initiations per quarter", subtitle = label) +
      theme_lancet_pub() +
      theme(plot.subtitle = element_text(face = "bold", color = color,
                                         size = 11, hjust = 0))
  }
  figS9 <- build_gam_panel(gam_c, cotri_q, "(A) Cotrimoxazole - GAM NB",
                           "#00468BFF") /
    build_gam_panel(gam_t, tpt_q,   "(B) TPT - GAM NB",
                    "#AD002AFF")
  save_fig(file.path(fig_supp, "FigS9_gam_smooth.png"),
           figS9, width = 180, height = 170)
}
# ---- FigS10: Distributed-lag effects ---------------------------------------
cat("  FigS10: Distributed lag\n")
if (exists("dl_all") && !is.null(dl_all) && nrow(dl_all) > 0) {
  figS10 <- ggplot(dl_all, aes(x = factor(Lag_Quarters), y = IRR,
                               color = Outcome, group = Outcome)) +
    geom_errorbar(aes(ymin = IRR_CI_Lower, ymax = IRR_CI_Upper), width = 0.15,
                  position = position_dodge(width = 0.3), linewidth = 0.6) +
    geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Cotrimoxazole" = "#00468BFF",
                                  "TPT"           = "#AD002AFF")) +
    scale_y_log10() +
    labs(x = "Lag (quarters from war onset)", y = "IRR (95% CI)", color = NULL) +
    theme_lancet_pub() + theme(legend.position = "bottom")
  save_fig(file.path(fig_supp, "FigS10_distributed_lag.png"),
           figS10, width = 160, height = 120)
}
# ---- FigS11: Cross-correlation (PNG saved in /diagnostics/) ----------------
# Already saved during 13H
# ---- FigS12: Bai-Perron breakpoints ----------------------------------------
cat("  FigS12: Bai-Perron breakpoints\n")
make_bp_panel <- function(bp_res, ts_df, label, color) {
  if (is.null(bp_res)) return(NULL)
  bps <- bp_res$model$breakpoints
  if (is.null(bps) || any(is.na(bps))) bps <- integer(0)
  bp_dates <- if (length(bps) > 0) ts_df$date[bps] else as.Date(integer(0))
  ggplot(ts_df, aes(x = date, y = n_events)) +
    geom_line(color = color, linewidth = 0.5) +
    geom_point(color = color, size = 1.4) +
    { if (length(bp_dates) > 0)
      geom_vline(xintercept = bp_dates, linetype = "dashed",
                 color = "#AD002A", linewidth = 0.7) } +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(x = NULL, y = "Initiations per quarter", subtitle = label) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", color = color,
                                       size = 11, hjust = 0))
}
bp_panel_c <- make_bp_panel(bp_c, cotri_q, "(A) Cotrimoxazole", "#00468BFF")
bp_panel_t <- make_bp_panel(bp_t, tpt_q,   "(B) TPT",            "#AD002AFF")
if (!is.null(bp_panel_c) && !is.null(bp_panel_t)) {
  figS12 <- bp_panel_c / bp_panel_t
  save_fig(file.path(fig_supp, "FigS12_baiperron_breakpoints.png"),
           figS12, width = 180, height = 170)
}
# ---- FigS13: Zivot-Andrews breakpoint --------------------------------------
cat("  FigS13: Zivot-Andrews\n")
if (exists("za_all") && nrow(za_all) > 0) {
  za_plot_df <- za_all %>%
    mutate(Break_Date = as.Date(Break_Date))
  figS13 <- ggplot(za_plot_df, aes(x = Outcome, y = ZA_Statistic, fill = Outcome)) +
    geom_col(width = 0.5) +
    geom_hline(aes(yintercept = CV_5pct), linetype = "dashed",
               color = "#AD002A", linewidth = 0.5) +
    geom_text(aes(label = paste("Break:", Break_Date)),
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("Cotrimoxazole" = "#00468BFF",
                                 "TPT"           = "#AD002AFF")) +
    labs(x = NULL, y = "Zivot-Andrews statistic (dashed = 5% CV)") +
    theme_lancet_pub() + theme(legend.position = "none")
  save_fig(file.path(fig_supp, "FigS13_zivot_andrews.png"),
           figS13, width = 160, height = 110)
}
# ---- FigS14: Forecast benchmarks -------------------------------------------
cat("  FigS14: Forecast benchmarks\n")
if (exists("fb_all") && !is.null(fb_all) && nrow(fb_all) > 0) {
  figS14 <- ggplot(fb_all, aes(x = Model, y = Pct_Gap, fill = Outcome)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = c("Cotrimoxazole" = "#00468BFF",
                                 "TPT"           = "#AD002AFF")) +
    labs(x = NULL, y = "Observed − predicted (% of predicted)", fill = NULL) +
    theme_lancet_pub() + theme(legend.position = "bottom")
  save_fig(file.path(fig_supp, "FigS14_forecast_benchmarks.png"),
           figS14, width = 180, height = 110)
}
# ---- FigS15: Quantile ITS coefficients -------------------------------------
cat("  FigS15: Quantile ITS\n")
if (exists("qreg_all") && !is.null(qreg_all) && nrow(qreg_all) > 0) {
  qr_war <- qreg_all %>% filter(Parameter == "level_war")
  figS15 <- ggplot(qr_war, aes(x = factor(Tau), y = Estimate,
                               color = Outcome, group = Outcome)) +
    geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.12,
                  position = position_dodge(0.3), linewidth = 0.6) +
    geom_point(size = 2.4, position = position_dodge(0.3)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Cotrimoxazole" = "#00468BFF",
                                  "TPT"           = "#AD002AFF")) +
    labs(x = "Quantile (τ)", y = "War-onset effect", color = NULL) +
    theme_lancet_pub() + theme(legend.position = "bottom")
  save_fig(file.path(fig_supp, "FigS15_quantile_its.png"),
           figS15, width = 160, height = 110)
}
# ---- FigS16: Bayesian posterior --------------------------------------------
cat("  FigS16: Bayesian posterior (if rstanarm ran)\n")
if (exists("bay_all") && !is.null(bay_all) && nrow(bay_all) > 0) {
  bay_war <- bay_all %>% filter(Parameter == "level_war")
  figS16 <- ggplot(bay_war, aes(x = Outcome, y = Posterior_Mean,
                                color = Outcome)) +
    geom_errorbar(aes(ymin = Cred_2_5, ymax = Cred_97_5), width = 0.15,
                  linewidth = 0.7) +
    geom_point(size = 3.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Cotrimoxazole" = "#00468BFF",
                                  "TPT"           = "#AD002AFF")) +
    labs(x = NULL, y = "Posterior war effect (log-IRR; 95% credible)") +
    theme_lancet_pub() + theme(legend.position = "none")
  save_fig(file.path(fig_supp, "FigS16_bayesian_posterior.png"),
           figS16, width = 140, height = 100)
}
# ---- FigS17: CausalImpact plots --------------------------------------------
# Note: do NOT use plot(ci_obj) + theme_*(...) - CausalImpact's plot method
# returns a faceted ggplot whose strip/facet structure breaks when overridden
# externally (renders to a near-blank canvas). Build the three panels manually
# from ci_obj$series instead.
cat("  FigS17: CausalImpact (if it ran)\n")
build_ci_three_panel <- function(ci_obj, outcome_label, accent_color) {
  if (is.null(ci_obj) || is.null(ci_obj$series)) return(NULL)
  s <- as.data.frame(ci_obj$series); s$date <- zoo::index(ci_obj$series)
  if (!inherits(s$date, "Date")) s$date <- as.Date(s$date)
  int_date <- tryCatch({
    pp <- ci_obj$model$post.period
    if (!is.null(pp)) as.Date(pp[1]) else config$war_start
  }, error = function(e) config$war_start)
  pA <- ggplot(s, aes(x = date)) +
    annotate("rect", xmin = config$war_start, xmax = config$war_end,
             ymin = -Inf, ymax = Inf, fill = "#AD002A", alpha = 0.10) +
    geom_ribbon(aes(ymin = point.pred.lower, ymax = point.pred.upper),
                fill = accent_color, alpha = 0.20) +
    geom_line(aes(y = point.pred), color = accent_color,
              linetype = "dashed", linewidth = 0.6) +
    geom_line(aes(y = response), color = "grey25", linewidth = 0.55) +
    geom_vline(xintercept = int_date, linetype = "dashed",
               color = "#AD002A", linewidth = 0.6) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Initiations per quarter",
         subtitle = paste0(outcome_label, " - observed vs counterfactual")) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 10,
                                       color = accent_color, hjust = 0))
  pB <- ggplot(s, aes(x = date)) +
    annotate("rect", xmin = config$war_start, xmax = config$war_end,
             ymin = -Inf, ymax = Inf, fill = "#AD002A", alpha = 0.10) +
    geom_ribbon(aes(ymin = point.effect.lower, ymax = point.effect.upper),
                fill = accent_color, alpha = 0.20) +
    geom_line(aes(y = point.effect), color = accent_color, linewidth = 0.55) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = int_date, linetype = "dashed",
               color = "#AD002A", linewidth = 0.6) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Pointwise effect",
         subtitle = "Pointwise effect (observed − counterfactual)") +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 10, hjust = 0))
  pC <- ggplot(s, aes(x = date)) +
    annotate("rect", xmin = config$war_start, xmax = config$war_end,
             ymin = -Inf, ymax = Inf, fill = "#AD002A", alpha = 0.10) +
    geom_ribbon(aes(ymin = cum.effect.lower, ymax = cum.effect.upper),
                fill = accent_color, alpha = 0.20) +
    geom_line(aes(y = cum.effect), color = accent_color, linewidth = 0.55) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = int_date, linetype = "dashed",
               color = "#AD002A", linewidth = 0.6) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Cumulative effect",
         subtitle = "Cumulative effect over time") +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 10, hjust = 0))
  pA / pB / pC + plot_layout(heights = c(1, 1, 1))
}
if (exists("ci_cotri") && !is.null(ci_cotri)) {
  fig_ci_c <- build_ci_three_panel(ci_cotri, "(A) Cotrimoxazole", "#00468BFF")
  if (!is.null(fig_ci_c))
    save_fig_pub(file.path(fig_supp, "FigS17_causalimpact_cotri"),
                 fig_ci_c, width = 180, height = 220)
}
if (exists("ci_tpt") && !is.null(ci_tpt)) {
  fig_ci_t <- build_ci_three_panel(ci_tpt, "(B) TPT", "#AD002AFF")
  if (!is.null(fig_ci_t))
    save_fig_pub(file.path(fig_supp, "FigS17_causalimpact_tpt"),
                 fig_ci_t, width = 180, height = 220)
}
# ---- FigS18: Mixed-effects NB ITS coefficients -----------------------------
cat("  FigS18: Mixed-effects ITS\n")
if (exists("me_all") && !is.null(me_all) && nrow(me_all) > 0) {
  me_keep <- me_all %>% filter(Parameter %in%
                                 c("level_war", "trend_war", "level_postwar"))
  figS18 <- ggplot(me_keep, aes(x = Parameter, y = IRR, color = Outcome)) +
    geom_errorbar(aes(ymin = IRR_CI_Lower, ymax = IRR_CI_Upper), width = 0.15,
                  position = position_dodge(width = 0.4), linewidth = 0.6) +
    geom_point(size = 2.4, position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Cotrimoxazole" = "#00468BFF",
                                  "TPT"           = "#AD002AFF")) +
    scale_y_log10() +
    labs(x = NULL, y = "IRR (95% CI; random year)", color = NULL) +
    theme_lancet_pub() + theme(legend.position = "bottom")
  save_fig(file.path(fig_supp, "FigS18_mixed_effects.png"),
           figS18, width = 160, height = 110)
}
# ---- FigS19: Periodogram ---------------------------------------------------
cat("  FigS19: Periodogram\n")
pg_path <- file.path(output_root, "decomposition", "periodogram.csv")
if (file.exists(pg_path)) {
  pg_df <- read.csv(pg_path)
  figS19 <- ggplot(pg_df, aes(x = Period_Quarters, y = Spectral_Density,
                              color = Outcome)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = c("Cotrimoxazole" = "#00468BFF",
                                  "TPT"           = "#AD002AFF")) +
    scale_x_continuous(limits = c(0, 12)) +
    labs(x = "Period (quarters)", y = "Spectral density", color = NULL) +
    theme_lancet_pub() + theme(legend.position = "bottom")
  save_fig(file.path(fig_supp, "FigS19_periodogram.png"),
           figS19, width = 160, height = 110)
}
# ==============================================================================
# 15B. V13 NEW FIGURES
# ==============================================================================
cat("\n══ V13 NEW FIGURES ══════════════════════════════════════════════════════\n")
# ---- FigV13_1: Recording timeline CTX --------------------------------------
cat("  FigV13_1: Recording-pattern timeline (CTX)\n")
if (exists("rt_c") && nrow(rt_c) > 0) {
  rt_c_long <- rt_c %>%
    transmute(date,
              `n records (total)`     = n_records_total,
              `n initiated`           = n_initiated,
              `% initiated`           = pct_initiated,
              `% with patient ID`     = pct_with_id,
              `% with sex`            = pct_with_sex,
              `% with age`            = pct_with_age,
              `n strata active`       = n_strata_active) %>%
    pivot_longer(-date, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c(
      "n records (total)", "n initiated", "% initiated",
      "% with patient ID", "% with sex", "% with age", "n strata active")))
  figV13_1 <- ggplot(rt_c_long, aes(x = date, y = value)) +
    annots +
    geom_line(color = "#00468BFF", linewidth = 0.55) +
    geom_point(color = "#00468BFF", size = 1.1) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    labs(x = NULL, y = NULL) + theme_lancet_pub()
  save_fig_pub(file.path(fig_v13, "FigV13_1_recording_timeline_cotri"),
               figV13_1, width = 200, height = 230)
}
# ---- FigV13_2: Recording timeline TPT --------------------------------------
cat("  FigV13_2: Recording-pattern timeline (TPT)\n")
if (exists("rt_t") && nrow(rt_t) > 0) {
  rt_t_long <- rt_t %>%
    transmute(date,
              `n records (total)`     = n_records_total,
              `n initiated`           = n_initiated,
              `% initiated`           = pct_initiated,
              `% with patient ID`     = pct_with_id,
              `% with sex`            = pct_with_sex,
              `% with age`            = pct_with_age,
              `n strata active`       = n_strata_active) %>%
    pivot_longer(-date, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c(
      "n records (total)", "n initiated", "% initiated",
      "% with patient ID", "% with sex", "% with age", "n strata active")))
  figV13_2 <- ggplot(rt_t_long, aes(x = date, y = value)) +
    annots +
    geom_line(color = "#AD002AFF", linewidth = 0.55) +
    geom_point(color = "#AD002AFF", size = 1.1) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    labs(x = NULL, y = NULL) + theme_lancet_pub()
  save_fig_pub(file.path(fig_v13, "FigV13_2_recording_timeline_tpt"),
               figV13_2, width = 200, height = 230)
}
# ---- FigV13_3: Date source distribution stacked ----------------------------
cat("  FigV13_3: Date-source distribution\n")
if (exists("ds_all") && nrow(ds_all) > 0) {
  ds_colors <- c("initiation"   = "#00468BFF",
                 "stop"         = "#925E9F",
                 "discontinued" = "#E7298A",
                 "completed"    = "#42B540",
                 "not_init"     = "grey70")
  figV13_3 <- ggplot(ds_all, aes(x = yr, y = n, fill = source)) +
    geom_col(width = 0.85) +
    facet_wrap(~ outcome, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = ds_colors, name = "Date source") +
    scale_x_continuous(breaks = seq(2005, 2025, 2)) +
    labs(x = NULL, y = "Records per year") +
    theme_lancet_pub() + theme(legend.position = "bottom")
  save_fig_pub(file.path(fig_v13, "FigV13_3_date_source_stacked"),
               figV13_3, width = 200, height = 170)
}
# ---- FigV13_4: TPT type × year heatmap -------------------------------------
cat("  FigV13_4: TPT regimen type × year heatmap\n")
if (exists("type_yearly") && nrow(type_yearly) > 0) {
  figV13_4 <- ggplot(type_yearly, aes(x = yr, y = type_label, fill = n)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(n > 0, format(n, big.mark = ","), "")),
              color = "grey20", size = 2.6) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08306B", name = "n records",
                        trans = "sqrt", labels = scales::comma) +
    scale_x_continuous(breaks = seq(2005, 2025, 2)) +
    labs(x = NULL, y = "tb_prophylaxis_type code") +
    theme_lancet_pub()
  save_fig_pub(file.path(fig_v13, "FigV13_4_tpt_type_heatmap"),
               figV13_4, width = 200, height = 130)
}
# ---- FigV13_5: TPT type-filter sensitivity comparison ----------------------
cat("  FigV13_5: TPT type-filter sensitivity\n")
if (exists("type_compare") && nrow(type_compare) > 0) {
  tc_long <- bind_rows(
    type_compare %>% transmute(Model, Period = "War",
                               Estimate = War_Estimate,
                               CI_Low = War_CI_Low, CI_High = War_CI_High),
    type_compare %>% transmute(Model, Period = "Post-war",
                               Estimate = Postwar_Estimate,
                               CI_Low = Postwar_CI_Low, CI_High = Postwar_CI_High)
  )
  figV13_5 <- ggplot(tc_long, aes(x = Estimate, y = Model, color = Period)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = CI_Low, xmax = CI_High),
                  orientation = "y", width = 0.15,
                  position = position_dodge(width = 0.4), linewidth = 0.55) +
    geom_point(size = 2.6, position = position_dodge(width = 0.4)) +
    scale_color_manual(values = c("War" = "#AD002AFF",
                                  "Post-war" = "#42B540FF")) +
    labs(x = "GLS-AR1 estimate (events per quarter; 95% CI)",
         y = NULL, color = NULL) +
    theme_lancet_pub() + theme(legend.position = "bottom")
  save_fig_pub(file.path(fig_v13, "FigV13_5_tpt_type_filter"),
               figV13_5, width = 180, height = 110)
}
# ---- FigV13_6 / FigV13_7: Under-recording forest plots ---------------------
cat("  FigV13_6/7: Under-recording sensitivity forest plots\n")
make_underrec_forest <- function(target_lab, file_stem) {
  if (!exists("ur_results") || nrow(ur_results) == 0) return(NULL)
  d <- ur_results %>% filter(Target == target_lab) %>%
    mutate(label = sprintf("%2d%% inflation", round(U * 100)),
           label = factor(label, levels = rev(unique(label))),
           sig_color = ifelse(War_Significant, "#AD002A", "#00468B"))
  ggplot(d, aes(x = War_Estimate, y = label, color = sig_color)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = War_CI_Low, xmax = War_CI_High),
                  orientation = "y", width = 0.15, linewidth = 0.55) +
    geom_point(size = 2.4, shape = 18) +
    scale_color_identity() +
    facet_wrap(~ Outcome, scales = "free_x") +
    labs(x = "War-onset effect (GLS-AR1; 95% CI)",
         y = NULL,
         subtitle = ifelse(target_lab == "war",
                           "(A) War-only inflation",
                           "(B) War + post-war inflation")) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 11, hjust = 0))
}
figV13_6 <- make_underrec_forest("war",             "war_only")
figV13_7 <- make_underrec_forest("war_and_postwar", "war_and_postwar")
if (!is.null(figV13_6))
  save_fig_pub(file.path(fig_v13, "FigV13_6_underrecording_war_only"),
               figV13_6, width = 200, height = 130)
if (!is.null(figV13_7))
  save_fig_pub(file.path(fig_v13, "FigV13_7_underrecording_war_postwar"),
               figV13_7, width = 200, height = 130)
# ---- FigV13_8: TPT per 100 ART starts --------------------------------------
cat("  FigV13_8: TPT per 100 ART starts\n")
if (exists("tpt_rate") && nrow(tpt_rate) > 0) {
  tpt_rate_long <- tpt_rate %>%
    transmute(date,
              `ART starts per quarter`     = n_art_starts,
              `TPT initiations per quarter`= n_events_tpt,
              `TPT per 100 ART starts`     = tpt_per_100_art) %>%
    pivot_longer(-date, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c(
      "ART starts per quarter", "TPT initiations per quarter",
      "TPT per 100 ART starts")))
  figV13_8 <- ggplot(tpt_rate_long, aes(x = date, y = value)) +
    annots +
    geom_line(color = "#AD002AFF", linewidth = 0.55) +
    geom_point(color = "#AD002AFF", size = 1.1) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = NULL) + theme_lancet_pub()
  save_fig_pub(file.path(fig_v13, "FigV13_8_tpt_per_100_art_starts"),
               figV13_8, width = 180, height = 180)
}
cat("\n  All figures saved.\n\n")
# ==============================================================================
# 16. SESSION INFO
# ==============================================================================
cat("══ SECTION 16: SESSION INFO ═════════════════════════════════════════════\n")
session_path <- file.path(output_root, "validation", "session_info.txt")
sink(session_path)
cat("ITS V13 Session Information\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n\n")
print(sessionInfo())
sink()
cat("    saved:", normalizePath(session_path, mustWork = FALSE), "\n\n")
# ==============================================================================
# 17. FILE MANIFEST
# ==============================================================================
cat("══ SECTION 17: FILE MANIFEST ════════════════════════════════════════════\n")
all_files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
manifest <- data.frame(
  Folder = dirname(gsub(paste0(output_root, "/?"), "", all_files)),
  File   = basename(all_files),
  Size_KB = round(file.info(all_files)$size / 1024, 1),
  Modified = format(file.info(all_files)$mtime, "%Y-%m-%d %H:%M"),
  stringsAsFactors = FALSE
) %>% arrange(Folder, File)
manifest_path <- file.path(output_root, "FILE_MANIFEST.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
cat(sprintf("    saved: %s   (%s files)\n",
            normalizePath(manifest_path, mustWork = FALSE),
            format(nrow(manifest), big.mark = ",")))
# ==============================================================================
# 18. MASTER WORKBOOK - ITS_V13_MASTER_RESULTS.xlsx
# ==============================================================================
cat("\n══ SECTION 18: MASTER WORKBOOK ══════════════════════════════════════════\n")
excel_path <- file.path(output_root, "ITS_V13_MASTER_RESULTS.xlsx")
# Gather all CSVs except FILE_MANIFEST.csv
csv_files <- list.files(output_root, pattern = "\\.csv$", recursive = TRUE,
                        full.names = TRUE)
csv_files <- csv_files[!grepl("FILE_MANIFEST\\.csv$", csv_files)]
# Build sheet-name candidate (folder_file, strip illegal chars, 31-char limit)
sanitize_sheet <- function(folder, file) {
  base <- paste0(folder, "__", tools::file_path_sans_ext(file))
  base <- gsub("[\\/?*\\[\\]:'\"]", "_", base)
  base <- gsub("[[:space:]]+", "_", base)
  base <- substr(base, 1, 31)
  base
}
wb <- createWorkbook()
# 00_TOC sheet (built first - but data populated last after we know sheet names)
addWorksheet(wb, "00_TOC", tabColour = "#00468B")
# Style definitions
header_style <- createStyle(textDecoration = "bold", fontColour = "white",
                            fgFill = "#00468B", halign = "center",
                            border = "bottom", borderColour = "#00468B")
toc_header_style <- createStyle(textDecoration = "bold", fontSize = 12,
                                fontColour = "white", fgFill = "#00468B",
                                halign = "center")
# Process each CSV
toc_entries <- list()
used_names  <- character(0)
for (f in csv_files) {
  rel <- gsub(paste0(output_root, "/?"), "", f)
  folder <- dirname(rel); file <- basename(rel)
  sn <- sanitize_sheet(folder, file)
  # Ensure unique
  base_sn <- sn; ctr <- 1
  while (sn %in% used_names) {
    suffix <- paste0("_", ctr)
    sn <- paste0(substr(base_sn, 1, 31 - nchar(suffix)), suffix)
    ctr <- ctr + 1
  }
  used_names <- c(used_names, sn)
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) next
  addWorksheet(wb, sn)
  writeData(wb, sn, df, headerStyle = header_style)
  freezePane(wb, sn, firstActiveRow = 2)
  setColWidths(wb, sn, cols = seq_len(ncol(df)), widths = "auto")
  toc_entries[[length(toc_entries) + 1]] <- data.frame(
    Sheet = sn, Folder = folder, File = file,
    Rows = nrow(df), Columns = ncol(df),
    stringsAsFactors = FALSE)
}
toc_df <- if (length(toc_entries) > 0) bind_rows(toc_entries) %>%
  arrange(Folder, File) else data.frame()
# Write TOC
writeData(wb, "00_TOC",
          data.frame(Title = "ITS V13 Master Results - Tigray HIV care continuum",
                     Generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     Outcomes = "Cotrimoxazole prophylaxis & TB preventive therapy",
                     N_Sheets = nrow(toc_df) + 1L,
                     stringsAsFactors = FALSE),
          startRow = 1, headerStyle = toc_header_style)
writeData(wb, "00_TOC", toc_df, startRow = 7, headerStyle = header_style)
freezePane(wb, "00_TOC", firstActiveRow = 8)
setColWidths(wb, "00_TOC", cols = 1:5, widths = "auto")
# Add hyperlinks from TOC to each sheet
if (nrow(toc_df) > 0) {
  for (i in seq_len(nrow(toc_df))) {
    writeFormula(wb, "00_TOC",
                 x = makeHyperlinkString(sheet = toc_df$Sheet[i], row = 1,
                                         col = 1, text = toc_df$Sheet[i]),
                 startRow = 7 + i, startCol = 1)
  }
}
# Save
saveWorkbook(wb, excel_path, overwrite = TRUE)
cat(sprintf("    saved: %s   (%s sheets, %.1f MB)\n",
            normalizePath(excel_path, mustWork = FALSE),
            format(nrow(toc_df) + 1L, big.mark = ","),
            file.info(excel_path)$size / 1024^2))
# ==============================================================================
# DONE
# ==============================================================================
cat("\n############################################################\n")
cat("#  V13 PART 2 complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("#  Master workbook: ITS_V13_MASTER_RESULTS.xlsx\n")
cat("#  Manifest: FILE_MANIFEST.csv\n")
cat("############################################################\n")

################################################################################
#  STANDALONE DIAGNOSTICS FIGURE BUILDER
#  -------------------------------------
#  * Loads v13_part1_checkpoint.RData if objects are not already in memory.
#  * Writes residual-diagnostic figures to figures/supplement/ as FigS3a-d.
#  * "Not-produced check": before building, scans the target folder and only
#    builds figures whose PNG file does NOT already exist (or is empty / a
#    blank stub). Use FORCE = TRUE to rebuild everything anyway.
#  * After building, prints a manifest of what was found, what was built,
#    and what is still missing.
#
#  Run this file on its own:
#      source("ITS_V13_diagnostics_standalone.R")
################################################################################
FORCE <- FALSE   # set TRUE to rebuild even if file already exists & non-empty
# ---- 1. Bring everything we need into scope -------------------------------
# Either we already have objects in memory (e.g. continuing a session), or
# we load the Part-1 checkpoint.
need <- c("m_cotri", "m_tpt", "output_root",
          "save_fig_pub", "theme_lancet_pub")
missing_objs <- need[!vapply(need, exists, logical(1), where = globalenv())]
if (length(missing_objs) > 0) {
  candidates <- c(
    file.path("D:/PhD project Data/New data set/FINAL_ANALYSIS_ALL",
              "HIV_care _outcomes/TB/incidence/tpt/prophylactic",
              "Final_ITS_V13_Results", "v13_part1_checkpoint.RData"),
    file.path(getwd(), "Final_ITS_V13_Results", "v13_part1_checkpoint.RData"),
    file.path(getwd(), "v13_part1_checkpoint.RData"))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit))
    stop("Cannot find v13_part1_checkpoint.RData. Run Part 1 first.")
  cat("Loading checkpoint:", hit, "\n"); flush.console()
  load(hit, envir = globalenv())
}
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
})
# ---- 2. Target paths & not-produced check ---------------------------------
fig_supp <- file.path(output_root, "figures/supplement")
dir.create(fig_supp, showWarnings = FALSE, recursive = TRUE)
specs <- tibble::tribble(
  ~slot,                              ~label,                       ~color,        ~model_obj_name, ~resid_type,
  "FigS3a_diagnostics_cotri_gls",     "CTX GLS-AR1 (primary)",      "#00468BFF",   "m_cotri",       "normalized",
  "FigS3b_diagnostics_tpt_gls",       "TPT GLS-AR1 (primary)",      "#AD002AFF",   "m_tpt",         "normalized",
  "FigS3c_diagnostics_cotri_nb",      "CTX Negative Binomial",      "#00468BFF",   "nb_cotri",      "deviance",
  "FigS3d_diagnostics_tpt_nb",        "TPT Negative Binomial",      "#AD002AFF",   "nb_tpt",        "deviance"
)
specs$png_path  <- file.path(fig_supp, paste0(specs$slot, ".png"))
specs$exists    <- file.exists(specs$png_path)
specs$size_kb   <- ifelse(specs$exists,
                          round(file.info(specs$png_path)$size / 1024, 1),
                          NA_real_)
# Treat <50 KB at 600 dpi over 200×200 mm as a blank/stub file (real
# diagnostic figures are typically 400 KB - 2 MB).
specs$is_blank  <- ifelse(specs$exists, specs$size_kb < 50, NA)
specs$to_build  <- FORCE | !specs$exists | (specs$exists & specs$is_blank %in% TRUE)
cat("\n══ NOT-PRODUCED CHECK ══════════════════════════════════════════════════\n")
print(specs %>% dplyr::select(slot, exists, size_kb, is_blank, to_build),
      row.names = FALSE)
cat("\n")
if (!any(specs$to_build)) {
  cat("All four diagnostic figures already produced and non-blank.\n")
  cat("Set FORCE <- TRUE at the top of the script to rebuild anyway.\n")
  return(invisible())
}
# ---- 3. Panel builders -----------------------------------------------------
ggdiag_residuals <- function(rv, color) {
  d <- data.frame(t = seq_along(rv), r = rv)
  ggplot(d, aes(t, r)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#AD002A",
               linewidth = 0.4) +
    geom_line(color = color, linewidth = 0.45) +
    geom_point(color = color, size = 0.9) +
    labs(x = "Observation index", y = "Residual",
         subtitle = "Residuals over time") +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
}
ggdiag_hist <- function(rv, color) {
  d <- data.frame(r = rv); m <- mean(rv); s <- sd(rv)
  xs <- seq(min(rv) - s, max(rv) + s, length.out = 200)
  ndf <- data.frame(x = xs, y = dnorm(xs, m, s))
  ggplot(d, aes(r)) +
    geom_histogram(aes(y = after_stat(density)), bins = 25,
                   fill = color, color = "white", alpha = 0.7) +
    geom_line(data = ndf, aes(x, y), color = "#AD002A", linewidth = 0.6) +
    labs(x = "Residual", y = "Density",
         subtitle = "Histogram + normal overlay") +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
}
ggdiag_qq <- function(rv, color) {
  qq <- qqnorm(rv, plot.it = FALSE)
  d <- data.frame(theo = qq$x, samp = qq$y)
  ggplot(d, aes(theo, samp)) +
    geom_qq_line(aes(sample = samp), color = "#AD002A", linewidth = 0.5) +
    geom_point(color = color, size = 1.1, alpha = 0.8) +
    labs(x = "Theoretical quantile", y = "Sample quantile",
         subtitle = "Normal Q-Q") +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
}
ggdiag_acf_like <- function(rv, color, type = c("ACF", "PACF"), lag.max = 20) {
  type <- match.arg(type)
  ac <- if (type == "ACF") acf(rv, lag.max = lag.max, plot = FALSE)
  else pacf(rv, lag.max = lag.max, plot = FALSE)
  d <- data.frame(lag = as.numeric(ac$lag), acf = as.numeric(ac$acf))
  if (type == "ACF") d <- d[d$lag > 0, ]
  ci <- qnorm(0.975) / sqrt(length(rv))
  ggplot(d, aes(lag, acf)) +
    geom_hline(yintercept = c(-ci, ci), linetype = "dashed",
               color = "#AD002A", linewidth = 0.35) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.25) +
    geom_segment(aes(xend = lag, yend = 0), color = color, linewidth = 0.55) +
    geom_point(color = color, size = 1.4) +
    labs(x = "Lag (quarters)", y = type, subtitle = type) +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
}
ggdiag_squared <- function(rv, color) {
  d <- data.frame(t = seq_along(rv), r2 = rv^2)
  ggplot(d, aes(t, r2)) +
    geom_line(color = color, linewidth = 0.45) +
    labs(x = "Observation index", y = expression(Residual^2),
         subtitle = "Squared residuals (ARCH check)") +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
}
ggdiag_cusum <- function(rv, color) {
  s <- cumsum(rv) / sqrt(sum(rv^2))
  d <- data.frame(t = seq_along(s), c = s)
  ggplot(d, aes(t, c)) +
    geom_hline(yintercept = c(-1.96, 1.96) / sqrt(length(rv)),
               linetype = "dashed", color = "#AD002A", linewidth = 0.4) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.25) +
    geom_line(color = color, linewidth = 0.55) +
    labs(x = "Observation index", y = "CUSUM",
         subtitle = "Standardised CUSUM") +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
}
ggdiag_ljungbox <- function(rv, color) {
  lags <- c(4, 8, 12, 16, 20)
  pvals <- sapply(lags, function(L) tryCatch(
    Box.test(rv, L, "Ljung-Box")$p.value, error = function(e) NA))
  d <- data.frame(lag = lags, p = pvals)
  ggplot(d, aes(lag, p)) +
    geom_hline(yintercept = 0.05, linetype = "dashed",
               color = "#AD002A", linewidth = 0.4) +
    geom_line(color = color, linewidth = 0.5) +
    geom_point(color = color, size = 2.2) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Lag (quarters)", y = "p-value",
         subtitle = "Ljung-Box p-values") +
    theme_lancet_pub() +
    theme(plot.subtitle = element_text(face = "bold", size = 9, hjust = 0))
}
build_diagnostics_figure <- function(rv, label, color) {
  rv <- as.numeric(rv); rv <- rv[is.finite(rv)]
  if (length(rv) < 8) {
    message("Too few residuals for ", label); return(NULL)
  }
  safely <- function(fn, ...) tryCatch(fn(...), error = function(e) {
    message(label, " panel failed: ", conditionMessage(e)); NULL
  })
  panels <- list(
    safely(ggdiag_residuals, rv, color),
    safely(ggdiag_hist,      rv, color),
    safely(ggdiag_qq,        rv, color),
    safely(ggdiag_acf_like,  rv, color, "ACF"),
    safely(ggdiag_acf_like,  rv, color, "PACF"),
    safely(ggdiag_squared,   rv, color),
    safely(ggdiag_cusum,     rv, color),
    safely(ggdiag_ljungbox,  rv, color))
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (length(panels) == 0) return(NULL)
  Reduce(`+`, panels) +
    plot_layout(ncol = 3) +
    plot_annotation(title = paste0(label, " - residual diagnostics"),
                    theme = theme(plot.title = element_text(
                      face = "bold", size = 12, color = color,
                      margin = margin(b = 6))))
}
# ---- 4. Build the ones that need building ---------------------------------
cat("══ BUILDING DIAGNOSTIC FIGURES ══════════════════════════════════════════\n")
built <- character(0); skipped <- character(0); failed <- character(0)
for (i in seq_len(nrow(specs))) {
  row <- specs[i, ]
  if (!row$to_build) { skipped <- c(skipped, row$slot); next }
  m_obj <- tryCatch(get(row$model_obj_name, envir = globalenv()),
                    error = function(e) NULL)
  if (is.null(m_obj)) {
    cat(sprintf("  [SKIP] %s - model object '%s' not in memory.\n",
                row$slot, row$model_obj_name))
    failed <- c(failed, row$slot); next
  }
  rv <- tryCatch(residuals(m_obj, type = row$resid_type),
                 error = function(e) NULL)
  if (is.null(rv)) {
    cat(sprintf("  [SKIP] %s - residuals(%s, type='%s') failed.\n",
                row$slot, row$model_obj_name, row$resid_type))
    failed <- c(failed, row$slot); next
  }
  fig <- build_diagnostics_figure(rv, row$label, row$color)
  if (is.null(fig)) { failed <- c(failed, row$slot); next }
  save_fig_pub(file.path(fig_supp, row$slot),
               fig, width = 200, height = 200)
  built <- c(built, row$slot)
}
# ---- 5. Post-build manifest -----------------------------------------------
cat("\n══ POST-BUILD MANIFEST ══════════════════════════════════════════════════\n")
final <- specs %>%
  mutate(png_path = file.path(fig_supp, paste0(slot, ".png"))) %>%
  mutate(now_exists = file.exists(png_path),
         now_size_kb = ifelse(now_exists,
                              round(file.info(png_path)$size / 1024, 1),
                              NA_real_),
         status = case_when(
           slot %in% built   ~ "BUILT",
           slot %in% skipped ~ "SKIPPED (already non-blank)",
           slot %in% failed  ~ "FAILED (see message)",
           TRUE              ~ "?"))
print(final %>% dplyr::select(slot, status, now_exists, now_size_kb),
      row.names = FALSE)
cat(sprintf("\n  Built: %d   Skipped: %d   Failed: %d\n",
            length(built), length(skipped), length(failed)))
cat("  Folder:", normalizePath(fig_supp, mustWork = FALSE), "\n")
invisible(final)



###############################################################################
###############################################################################
###                                                                         ###
###                  PART 2: V13.1 REANALYSIS (SIX SECTIONS)                ###
###                                                                         ###
###############################################################################
###############################################################################

###############################################################################
# ITS V13.1 reanalysis (consolidated, standalone)
#
# Companion to: ITS_V13_Part1_Analysis.R and ITS_V13_Part2_Tables_Figures.R
# Manuscript: "Collapse of preventive HIV care during the Tigray war: a
#             20-year interrupted time series"
#
# Author:  Hafte Kahsay Kebede
# License: Apache-2.0
# Reproducibility:
#   R 4.6.0; nlme 3.1.169; MASS 7.3.65; survival 3.8.6; openxlsx 4.2.8.1;
#   dplyr 1.2.1; tidyr 1.3.2; lubridate 1.9.5; data.table 1.18.4
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Six targeted re-analyses requested during post-V13 manuscript review,
# bundled into one standalone script that reads the raw V13 inputs and
# produces a single Excel workbook plus four CSV tables. The script is
# self-contained: it does not require V13 Part 1 to have run, although
# it will use the V13 workspace objects if they are present.
#
#   1. NB under-recording sensitivity. Refits the V13 segmented negative
#      binomial model with war and post-war counts inflated by 0, 10, 25,
#      50, and 75 percent. Reports the war-onset incidence rate ratio
#      and 95 percent confidence interval under each scenario.
#
#   2. Cox model n and nevent extraction. Reports the V13-style Cox
#      specification (ever-initiators only) and a proper survival
#      specification with administrative censoring at the study end.
#
#   3. Per-period percentage of raw records resolving to initiation
#      events. Anchored on the appropriate record date column for each
#      cohort.
#
#   4. Drop-all-influential-quarters sensitivity. Refits the primary
#      GLS-AR(1) and NB models after jointly removing the quarters
#      identified as influential by the V13 leave-one-quarter-out check.
#
#   5. Table 1 cohort filter forensics and bug-fix. Reproduces V13's
#      published Table 1 (which contained a case_when fallthrough bug
#      that assigned all non-initiators to the post-war column), then
#      builds a corrected version.
#
#   6. Cox refit with correct row selection. The v2 patch used
#      distinct(patient_id, .keep_all=TRUE) which picked rows at random
#      per patient; this refit uses arrange(patient_id, desc(initiated),
#      start_date) so the earliest initiation row is preferred.
#
# DATA-CONSTRUCTION DIFFERENCE FROM V13
# -------------------------------------
# V13 parsed sex, DateOfBirth, and date_hiv_confirmed per row. If a
# patient demographic was blank on a given row, V13 treated it as
# missing for that row even when the same field was filled on another
# of the same patient's visits. V13.1 carries demographics forward
# within patient_id before de-duplication: for each unique patient,
# the first non-NA value seen in any of their records is propagated
# to all of their rows. A field is treated as truly missing only when
# it is blank on every row for the patient. Rows with no patient_id
# are not grouped and retain whatever they had on the specific row.
#
# INPUTS
# ------
#   data_dir/crtEthiopiaARTvisit_1.csv   (cotrimoxazole, part 1)
#   data_dir/crtEthiopiaARTvisit_2.csv   (cotrimoxazole, part 2)
#   data_dir/vcrtEthiopiaARTVisit_TB_All.csv  (tuberculosis preventive therapy)
#
# OUTPUTS
# -------
#   data_dir/Final_ITS_V13_Results/reanalysis/
#     ITS_V13_REANALYSIS_RESULTS.xlsx     (all six sections, 12 sheets)
#     table1_cotri_buggy_reproduction.csv
#     table1_cotri_corrected.csv
#     table1_tpt_buggy_reproduction.csv
#     table1_tpt_corrected.csv
#     reanalysis_console_summary.txt
###############################################################################


###############################################################################
# Section 0. Setup, paths, packages                                           #
###############################################################################

data_dir       <- "D:/PhD project Data/New data set/FINAL_ANALYSIS_ALL/HIV_care _outcomes/TB/incidence/tpt/prophylactic"
ctx_file_stems <- c("crtEthiopiaARTvisit_1", "crtEthiopiaARTvisit_2")
tpt_file_stems <- c("vcrtEthiopiaARTVisit_TB_All")

config <- list(
  study_start   = as.Date("2005-01-01"),
  study_end     = as.Date("2025-05-31"),
  war_start     = as.Date("2021-01-01"),
  postwar_start = as.Date("2023-01-01")
)

out_dir <- file.path(data_dir, "Final_ITS_V13_Results", "reanalysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
xlsx_path <- file.path(out_dir, "ITS_V13_REANALYSIS_RESULTS.xlsx")
log_path  <- file.path(out_dir, "reanalysis_console_summary.txt")

pkgs <- c("nlme", "MASS", "survival", "openxlsx", "dplyr", "readxl",
          "tidyr", "lubridate", "data.table", "tools")
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cran.r-project.org")
}
invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE,
                 warn.conflicts = FALSE))

set.seed(42)

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)

cat("ITS V13.1 reanalysis. Run started:", format(Sys.time()), "\n")
cat(strrep("=", 75), "\n\n")
cat("Working directory: ", getwd(), "\n")
cat("V13 data_dir:      ", data_dir, "\n")
cat("Output directory:  ", out_dir, "\n\n")

wb <- createWorkbook()


###############################################################################
# Section 0.5. Data construction (skipped if workspace has V13 objects)       #
###############################################################################

parse_date_flex <- function(x) {
  if (inherits(x, "Date")) {
    res <- x
  } else if (inherits(x, "POSIXt")) {
    res <- as.Date(x)
  } else {
    s <- as.character(x); s[s == "" | s == "NA"] <- NA
    res <- as.Date(rep(NA, length(s)))
    fmts <- c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%Y/%m/%d",
              "%d-%m-%Y", "%m-%d-%Y", "%d.%m.%Y")
    for (f in fmts) {
      idx <- is.na(res) & !is.na(s)
      if (!any(idx)) break
      res[idx] <- as.Date(s[idx], format = f)
    }
    idx <- is.na(res) & !is.na(s) & grepl("^[0-9]+(\\.[0-9]+)?$", s)
    if (any(idx)) {
      n <- suppressWarnings(as.numeric(s[idx]))
      res[idx] <- as.Date(n, origin = "1899-12-30")
    }
  }
  res[!is.na(res) & res <= as.Date("1900-12-31")] <- NA
  res
}

read_artvisit <- function(stem, dir) {
  candidates <- c(file.path(dir, paste0(stem, ".csv")),
                  file.path(dir, paste0(stem, ".xlsx")),
                  file.path(dir, paste0(stem, ".xls")),
                  file.path(dir, paste0(stem, ".rds")))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop(sprintf("Could not find '%s' in:\n  %s", stem, dir))
  cat(sprintf("  Reading: %s\n", basename(hit)))
  flush.console()
  ext <- tolower(tools::file_ext(hit))
  df <- switch(
    ext,
    csv  = data.table::fread(hit, na.strings = c("", "NA", "N/A"),
                             showProgress = FALSE,
                             colClasses = "character") |> as.data.frame(),
    xlsx = as.data.frame(readxl::read_excel(hit, guess_max = 100000,
                                            col_types = "text")),
    xls  = as.data.frame(readxl::read_excel(hit, guess_max = 100000,
                                            col_types = "text")),
    rds  = as.data.frame(readRDS(hit))
  )
  df[] <- lapply(df, function(col)
    if (is.character(col)) col else as.character(col))
  cat(sprintf("    %s rows, %d cols\n",
              format(nrow(df), big.mark = ","), ncol(df)))
  df
}

get_col <- function(df, target) {
  hits <- which(tolower(names(df)) == tolower(target))
  if (length(hits) == 0) return(NA_character_)
  names(df)[hits[1]]
}

# V13.1 helper: take the first non-NA value within a group.
.first_nonNA <- function(x) {
  v <- x[!is.na(x)]
  if (length(v) == 0) return(x[NA_integer_][seq_along(x)])
  out <- x
  out[] <- v[1]
  out
}

add_period_year_quarter <- function(df, date_col) {
  df %>% mutate(
    year = lubridate::year(.data[[date_col]]),
    quarter = lubridate::quarter(.data[[date_col]]),
    month = lubridate::month(.data[[date_col]]),
    period = factor(case_when(
      year <= 2019 ~ "Pre-war",   year == 2020 ~ "COVID-19",
      year <= 2022 ~ "War",       TRUE         ~ "Post-war"),
      levels = c("Pre-war", "COVID-19", "War", "Post-war"))
  )
}

# Period assigner without NA fallthrough (used for Table 1 correction)
assign_period_fixed <- function(date_vec) {
  yr <- lubridate::year(date_vec)
  out <- case_when(
    is.na(yr)  ~ NA_character_,
    yr <= 2019 ~ "Pre-war",
    yr == 2020 ~ "COVID-19",
    yr <= 2022 ~ "War",
    yr >= 2023 ~ "Post-war"
  )
  factor(out, levels = c("Pre-war", "COVID-19", "War", "Post-war"))
}

have_workspace_data <- all(c("cotri_q", "tpt_q", "cotri", "tpt") %in%
                             ls(envir = .GlobalEnv))

if (have_workspace_data) {
  cat("V13 workspace objects found. Skipping data construction.\n\n")
  ctx_combined_raw <- if (exists("ctx_combined", envir = .GlobalEnv))
    get("ctx_combined", envir = .GlobalEnv) else NULL
  tpt_raw_full <- if (exists("tpt_raw", envir = .GlobalEnv))
    get("tpt_raw", envir = .GlobalEnv) else NULL
} else {
  cat("V13 workspace objects not found. Rebuilding from raw inputs.\n\n")
  
  cat("Reading cotrimoxazole raw files.\n")
  ctx_raw_list <- lapply(ctx_file_stems, read_artvisit, dir = data_dir)
  ctx_raw_list <- lapply(ctx_raw_list, function(d) {
    names(d) <- tolower(names(d)); d
  })
  ctx_combined <- bind_rows(
    lapply(seq_along(ctx_raw_list), function(i) {
      d <- ctx_raw_list[[i]]; d$.source_file <- ctx_file_stems[i]; d
    })
  )
  ctx_combined <- ctx_combined %>%
    distinct(across(-.source_file), .keep_all = TRUE)
  
  cat("Reading tuberculosis preventive therapy raw file.\n")
  tpt_raw_list <- lapply(tpt_file_stems, read_artvisit, dir = data_dir)
  tpt_raw_list <- lapply(tpt_raw_list, function(d) {
    names(d) <- tolower(names(d)); d
  })
  tpt_raw <- tpt_raw_list[[1]] %>% distinct(.keep_all = TRUE)
  
  # Build cotrimoxazole event-level
  c_pid       <- get_col(ctx_combined, "uniqueartnumber")
  c_sex       <- get_col(ctx_combined, "sex")
  c_dob       <- get_col(ctx_combined, "dateofbirth")
  c_hivconf   <- get_col(ctx_combined, "date_hiv_confirmed")
  c_ctx_start <- get_col(ctx_combined, "cotrimoxazolestartdate")
  c_ctx_stop  <- get_col(ctx_combined, "cortimoxazole_stop_date")
  
  cotri <- ctx_combined %>%
    mutate(
      patient_id = if (!is.na(c_pid))
        trimws(as.character(.data[[c_pid]])) else NA_character_,
      patient_id = ifelse(is.na(patient_id) | patient_id == "" |
                            patient_id == "NA", NA_character_, patient_id),
      sex_raw = if (!is.na(c_sex))
        toupper(trimws(.data[[c_sex]])) else NA_character_,
      sex = case_when(sex_raw == "F" ~ "Female",
                      sex_raw == "M" ~ "Male", TRUE ~ NA_character_),
      DateOfBirth = if (!is.na(c_dob))
        parse_date_flex(.data[[c_dob]]) else as.Date(NA),
      date_hiv_confirmed = if (!is.na(c_hivconf))
        parse_date_flex(.data[[c_hivconf]]) else as.Date(NA),
      CotrimoxazoleStartDate = if (!is.na(c_ctx_start))
        parse_date_flex(.data[[c_ctx_start]]) else as.Date(NA),
      cortimoxazole_stop_date = if (!is.na(c_ctx_stop))
        parse_date_flex(.data[[c_ctx_stop]]) else as.Date(NA)
    )
  
  # V13.1 demographic carry-forward
  has_pid <- !is.na(cotri$patient_id)
  pid_part <- cotri[has_pid, ] %>%
    group_by(patient_id) %>%
    mutate(across(c(sex_raw, sex, DateOfBirth, date_hiv_confirmed),
                  .first_nonNA)) %>%
    ungroup()
  cotri <- bind_rows(pid_part, cotri[!has_pid, ])
  
  cotri <- cotri %>%
    mutate(
      cotri_start_date = dplyr::coalesce(CotrimoxazoleStartDate,
                                         cortimoxazole_stop_date),
      initiated = !is.na(cotri_start_date) &
        cotri_start_date >= config$study_start &
        cotri_start_date <= config$study_end
    ) %>%
    mutate(cotri_start_date = ifelse(initiated, cotri_start_date, NA) %>%
             as.Date(origin = "1970-01-01")) %>%
    mutate(
      age = as.numeric(difftime(
        dplyr::coalesce(cotri_start_date, date_hiv_confirmed),
        DateOfBirth, units = "days")) / 365.25,
      age = ifelse(age < 0 | age > 110, NA_real_, age),
      time_to_init_days = as.numeric(cotri_start_date - date_hiv_confirmed),
      time_to_init_days = ifelse(time_to_init_days < 0, NA_real_,
                                 time_to_init_days),
      record_date = date_hiv_confirmed
    ) %>%
    add_period_year_quarter("cotri_start_date")
  
  cotri_has <- cotri %>%
    filter(!is.na(patient_id)) %>%
    distinct(patient_id, cotri_start_date, .keep_all = TRUE)
  cotri_no <- cotri %>% filter(is.na(patient_id))
  cotri <- bind_rows(cotri_has, cotri_no)
  
  cat(sprintf("Cotrimoxazole: %s rows, %s initiations, %s unique initiating patients\n",
              format(nrow(cotri), big.mark = ","),
              format(sum(cotri$initiated, na.rm = TRUE), big.mark = ","),
              format(n_distinct(cotri$patient_id[cotri$initiated]),
                     big.mark = ",")))
  
  # Build TPT event-level
  t_pid       <- get_col(tpt_raw, "uano")
  t_sex       <- get_col(tpt_raw, "sex")
  t_dob       <- get_col(tpt_raw, "dateofbirth")
  t_artstart  <- get_col(tpt_raw, "art_start_date")
  t_regdate   <- get_col(tpt_raw, "registration_date")
  t_tbtype    <- get_col(tpt_raw, "tb_prophylaxis_type")
  t_inh_start <- get_col(tpt_raw, "inhprophylaxis_started_date")
  t_inh_disc  <- get_col(tpt_raw, "inhprophylaxisdiscontinueddate")
  t_inh_compl <- get_col(tpt_raw, "inhprophylaxiscompleteddate")
  
  tpt <- tpt_raw %>%
    mutate(
      patient_id = if (!is.na(t_pid))
        trimws(as.character(.data[[t_pid]])) else NA_character_,
      patient_id = ifelse(is.na(patient_id) | patient_id == "" |
                            patient_id == "NA", NA_character_, patient_id),
      sex_raw = if (!is.na(t_sex))
        toupper(trimws(.data[[t_sex]])) else NA_character_,
      sex = case_when(sex_raw == "F" ~ "Female",
                      sex_raw == "M" ~ "Male", TRUE ~ NA_character_),
      DateOfBirth = if (!is.na(t_dob))
        parse_date_flex(.data[[t_dob]]) else as.Date(NA),
      art_start_date = if (!is.na(t_artstart))
        parse_date_flex(.data[[t_artstart]]) else as.Date(NA),
      registration_date = if (!is.na(t_regdate))
        parse_date_flex(.data[[t_regdate]]) else as.Date(NA),
      tb_prophylaxis_type = if (!is.na(t_tbtype))
        suppressWarnings(as.numeric(.data[[t_tbtype]])) else NA_real_,
      inhprophylaxis_started_date = if (!is.na(t_inh_start))
        parse_date_flex(.data[[t_inh_start]]) else as.Date(NA),
      InhprophylaxisDiscontinuedDate = if (!is.na(t_inh_disc))
        parse_date_flex(.data[[t_inh_disc]]) else as.Date(NA),
      InhprophylaxisCompletedDate = if (!is.na(t_inh_compl))
        parse_date_flex(.data[[t_inh_compl]]) else as.Date(NA)
    )
  
  has_pid <- !is.na(tpt$patient_id)
  pid_part <- tpt[has_pid, ] %>%
    group_by(patient_id) %>%
    mutate(across(c(sex_raw, sex, DateOfBirth, art_start_date,
                    registration_date), .first_nonNA)) %>%
    ungroup()
  tpt <- bind_rows(pid_part, tpt[!has_pid, ])
  
  tpt <- tpt %>%
    filter(is.na(tb_prophylaxis_type) | tb_prophylaxis_type == 1) %>%
    mutate(
      tpt_start_date = dplyr::coalesce(inhprophylaxis_started_date,
                                       InhprophylaxisDiscontinuedDate,
                                       InhprophylaxisCompletedDate),
      initiated = !is.na(tpt_start_date) &
        tpt_start_date >= config$study_start &
        tpt_start_date <= config$study_end &
        (is.na(art_start_date) | tpt_start_date >= art_start_date)
    ) %>%
    mutate(tpt_start_date = ifelse(initiated, tpt_start_date, NA) %>%
             as.Date(origin = "1970-01-01")) %>%
    mutate(
      age = as.numeric(difftime(
        dplyr::coalesce(registration_date, art_start_date),
        DateOfBirth, units = "days")) / 365.25,
      age = ifelse(age < 0 | age > 110, NA_real_, age),
      time_to_init_days = as.numeric(tpt_start_date - art_start_date),
      time_to_init_days = ifelse(time_to_init_days < 0, NA_real_,
                                 time_to_init_days),
      record_date = dplyr::coalesce(art_start_date, registration_date)
    ) %>%
    add_period_year_quarter("tpt_start_date")
  
  tpt_has <- tpt %>%
    filter(!is.na(patient_id)) %>%
    distinct(patient_id, tpt_start_date, .keep_all = TRUE)
  tpt_no <- tpt %>% filter(is.na(patient_id))
  tpt <- bind_rows(tpt_has, tpt_no)
  
  cat(sprintf("Tuberculosis preventive therapy: %s rows, %s initiations, %s unique initiating patients\n\n",
              format(nrow(tpt), big.mark = ","),
              format(sum(tpt$initiated, na.rm = TRUE), big.mark = ","),
              format(n_distinct(tpt$patient_id[tpt$initiated]),
                     big.mark = ",")))
  
  # Quarterly time series
  full_grid <- expand.grid(year = 2005:2025, quarter = 1:4) %>%
    filter(!(year == 2025 & quarter > 2)) %>% arrange(year, quarter)
  bp_idx <- function(q, start_year)
    which(q$year == start_year & q$quarter == 1)[1]
  
  build_ts <- function(data, date_col, label) {
    q <- data %>%
      filter(initiated == TRUE) %>%
      count(year, quarter, name = "n_events") %>%
      right_join(full_grid, by = c("year", "quarter")) %>%
      mutate(n_events = replace_na(n_events, 0)) %>%
      arrange(year, quarter) %>%
      mutate(
        date = ymd(paste0(year, "-", quarter * 3 - 2, "-01")),
        time_index = row_number(), season = factor(quarter),
        period = factor(case_when(
          year <= 2019 ~ "Pre-war", year == 2020 ~ "COVID-19",
          year <= 2022 ~ "War", TRUE ~ "Post-war"),
          levels = c("Pre-war", "COVID-19", "War", "Post-war")),
        outcome = label
      )
    pat <- data %>%
      filter(initiated == TRUE) %>%
      mutate(yr = lubridate::year(.data[[date_col]]),
             qtr = lubridate::quarter(.data[[date_col]])) %>%
      group_by(yr, qtr) %>%
      summarise(n_patients = n_distinct(patient_id), .groups = "drop")
    q <- q %>%
      left_join(pat, by = c("year" = "yr", "quarter" = "qtr")) %>%
      mutate(n_patients = replace_na(n_patients, 0))
    bp_w <- bp_idx(q, 2021); bp_p <- bp_idx(q, 2023)
    q %>% mutate(
      trend_pre = time_index,
      covid = as.numeric(year == 2020),
      level_war = as.numeric(time_index >= bp_w),
      trend_war = ifelse(level_war == 1, time_index - bp_w + 1, 0),
      level_postwar = as.numeric(time_index >= bp_p),
      trend_postwar = ifelse(level_postwar == 1, time_index - bp_p + 1, 0)
    )
  }
  cotri_q <- build_ts(cotri, "cotri_start_date", "Cotrimoxazole")
  tpt_q   <- build_ts(tpt,   "tpt_start_date",   "TPT")
  
  ctx_combined_raw <- ctx_combined
  tpt_raw_full     <- tpt_raw
  assign("cotri",            cotri,            envir = .GlobalEnv)
  assign("tpt",              tpt,              envir = .GlobalEnv)
  assign("cotri_q",          cotri_q,          envir = .GlobalEnv)
  assign("tpt_q",            tpt_q,            envir = .GlobalEnv)
  assign("ctx_combined",     ctx_combined,     envir = .GlobalEnv)
  assign("tpt_raw",          tpt_raw,          envir = .GlobalEnv)
}

cotri    <- get("cotri",    envir = .GlobalEnv)
tpt      <- get("tpt",      envir = .GlobalEnv)
cotri_q  <- get("cotri_q",  envir = .GlobalEnv)
tpt_q    <- get("tpt_q",    envir = .GlobalEnv)


###############################################################################
# Section 1. Negative binomial under-recording sensitivity                    #
###############################################################################

cat(strrep("=", 75), "\n", sep = "")
cat("Section 1. Negative binomial under-recording sensitivity\n")
cat(strrep("=", 75), "\n\n", sep = "")

NB_FORMULA <- n_events ~ time_index + season + covid + level_war + trend_war +
  level_postwar + trend_postwar
GLS_FORMULA <- n_events ~ time_index + covid + level_war + trend_war +
  level_postwar + trend_postwar

inflate_counts <- function(qdata, inflation_pct) {
  if (inflation_pct == 0) return(qdata)
  m <- 1 / (1 - inflation_pct / 100)
  out <- qdata
  idx <- out$period %in% c("War", "Post-war")
  out$n_events[idx] <- round(out$n_events[idx] * m)
  out
}

fit_nb_war <- function(qdata) {
  fit <- tryCatch(
    suppressWarnings(MASS::glm.nb(NB_FORMULA, data = qdata,
                                  control = glm.control(maxit = 300,
                                                        epsilon = 1e-8))),
    error = function(e) list(err = conditionMessage(e))
  )
  if (!is.null(fit$err)) return(data.frame(IRR=NA, IRR_lower=NA, IRR_upper=NA,
                                           p_value=NA, converged=FALSE,
                                           note=fit$err))
  co <- summary(fit)$coefficients
  if (!"level_war" %in% rownames(co))
    return(data.frame(IRR=NA, IRR_lower=NA, IRR_upper=NA, p_value=NA,
                      converged=FALSE, note="level_war absent"))
  est <- co["level_war","Estimate"]; se <- co["level_war","Std. Error"]
  p   <- co["level_war","Pr(>|z|)"]
  ci  <- tryCatch(suppressMessages(confint(fit,"level_war")),
                  error = function(e) est + c(-1,1) * 1.96 * se)
  data.frame(IRR=exp(est), IRR_lower=exp(ci[1]), IRR_upper=exp(ci[2]),
             p_value=p, converged=isTRUE(fit$converged), note="")
}

run_nb_sensitivity <- function(qdata, lbl) {
  res <- list()
  for (pct in c(0, 10, 25, 50, 75)) {
    qinf <- inflate_counts(qdata, pct)
    r <- fit_nb_war(qinf)
    r$inflation_pct <- pct
    r$sum_war_obs <- sum(qinf$n_events[qinf$period == "War"])
    r$sum_postwar_obs <- sum(qinf$n_events[qinf$period == "Post-war"])
    cat(sprintf("  %-30s inflation=%2d%% : IRR=%.4f (%.4f, %.4f), p=%.4g\n",
                lbl, pct, r$IRR, r$IRR_lower, r$IRR_upper, r$p_value))
    res[[length(res) + 1]] <- r
  }
  d <- do.call(rbind, res)
  d$outcome <- lbl
  d$war_significant <- with(d, IRR_upper < 1 & !is.na(IRR_upper))
  d[, c("outcome", "inflation_pct", "sum_war_obs", "sum_postwar_obs",
        "IRR", "IRR_lower", "IRR_upper", "p_value", "war_significant",
        "converged", "note")]
}

nb_cotri <- run_nb_sensitivity(cotri_q, "Cotrimoxazole")
cat("\n")
nb_tpt   <- run_nb_sensitivity(tpt_q,   "Tuberculosis preventive therapy")
addWorksheet(wb, "NB_under_recording_cotri")
writeData(wb, "NB_under_recording_cotri", nb_cotri)
addWorksheet(wb, "NB_under_recording_tpt")
writeData(wb, "NB_under_recording_tpt", nb_tpt)


###############################################################################
# Section 2. Cox model with proper censoring                                  #
###############################################################################

cat("\n", strrep("=", 75), "\n", sep = "")
cat("Section 2. Cox model with proper censoring\n")
cat(strrep("=", 75), "\n\n", sep = "")

build_cox_data <- function(event_data, time_origin_col, start_date_col) {
  df <- event_data %>%
    filter(!is.na(patient_id) & patient_id != "" & patient_id != "NA")
  # Prefer the initiation row per patient (fixes the v2 patch row-selection bug)
  df <- df %>%
    arrange(patient_id, desc(initiated), .data[[start_date_col]]) %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup()
  df %>%
    mutate(
      time_origin = .data[[time_origin_col]],
      end_date = case_when(
        initiated ~ .data[[start_date_col]],
        TRUE      ~ config$study_end
      ),
      time_days = as.numeric(end_date - time_origin),
      time_months = time_days / 30.44,
      event = as.numeric(initiated),
      period_origin = factor(case_when(
        is.na(lubridate::year(time_origin)) ~ NA_character_,
        lubridate::year(time_origin) <= 2019 ~ "Pre-war",
        lubridate::year(time_origin) == 2020 ~ "COVID-19",
        lubridate::year(time_origin) <= 2022 ~ "War",
        lubridate::year(time_origin) >= 2023 ~ "Post-war"
      ), levels = c("Pre-war", "COVID-19", "War", "Post-war"))
    ) %>%
    filter(!is.na(time_months), !is.na(sex), !is.na(age), !is.na(period_origin),
           time_months >= 0, time_months <= 240)
}

fit_cox_and_report <- function(cox_df, label) {
  fit <- coxph(Surv(time_months, event) ~ period_origin + sex + age,
               data = cox_df, ties = "efron")
  zph_p <- tryCatch(cox.zph(fit)$table["GLOBAL", "p"],
                    error = function(e) NA_real_)
  cat(sprintf("  %s: n=%s patients; events=%s initiations; Schoenfeld p=%.3g\n",
              label,
              format(fit$n,      big.mark = ","),
              format(fit$nevent, big.mark = ","),
              zph_p))
  data.frame(
    outcome = label,
    n_patients = fit$n,
    n_events = fit$nevent,
    event_rate_pct = round(100 * fit$nevent / fit$n, 1),
    schoenfeld_global_p = zph_p
  )
}

cox_cotri_df <- build_cox_data(cotri, "date_hiv_confirmed", "cotri_start_date")
cox_tpt_df   <- build_cox_data(tpt,   "art_start_date",     "tpt_start_date")
cox_results <- rbind(
  fit_cox_and_report(cox_cotri_df, "Cotrimoxazole"),
  fit_cox_and_report(cox_tpt_df,   "Tuberculosis preventive therapy")
)
addWorksheet(wb, "cox_fit_summary")
writeData(wb, "cox_fit_summary", cox_results)


###############################################################################
# Section 3. Percentage of raw records resolving to initiation events         #
###############################################################################

cat("\n", strrep("=", 75), "\n", sep = "")
cat("Section 3. Percentage of raw records resolving to initiation events\n")
cat(strrep("=", 75), "\n\n", sep = "")

assign_period_str <- function(d) {
  ifelse(is.na(d), NA_character_,
         ifelse(d <= as.Date("2019-12-31"), "Pre-war",
                ifelse(d <= as.Date("2020-12-31"), "COVID-19",
                       ifelse(d <= as.Date("2022-12-31"), "War", "Post-war"))))
}

compute_pct <- function(raw_records, event_level, anchor_cols, label) {
  if (is.null(raw_records)) return(NULL)
  hit_col_lc <- anchor_cols[anchor_cols %in% tolower(names(raw_records))][1]
  if (is.na(hit_col_lc)) return(NULL)
  hit_col_real <- names(raw_records)[tolower(names(raw_records)) == hit_col_lc][1]
  rec_date <- parse_date_flex(raw_records[[hit_col_real]])
  per <- assign_period_str(rec_date)
  per[rec_date < config$study_start | rec_date > config$study_end] <-
    NA_character_
  rec_counts <- as.data.frame(table(period = per, useNA = "no"))
  names(rec_counts)[2] <- "n_records"
  init_counts <- event_level %>% filter(initiated == TRUE) %>%
    count(period, name = "n_initiations")
  out <- left_join(rec_counts, init_counts, by = "period") %>%
    mutate(outcome = label,
           pct_initiations = round(100 * n_initiations / n_records, 2)) %>%
    mutate(period = factor(period, levels =
                             c("Pre-war", "COVID-19", "War", "Post-war"))) %>%
    arrange(period) %>%
    select(outcome, period, n_records, n_initiations, pct_initiations)
  for (i in seq_len(nrow(out))) {
    cat(sprintf("  %-30s %-9s records=%s initiations=%s pct=%.2f%%\n",
                label, out$period[i],
                format(out$n_records[i],     big.mark = ","),
                format(out$n_initiations[i], big.mark = ","),
                out$pct_initiations[i]))
  }
  out
}

pct_cotri <- compute_pct(ctx_combined_raw, cotri,
                         c("date_hiv_confirmed"), "Cotrimoxazole")
pct_tpt   <- compute_pct(tpt_raw_full, tpt,
                         c("art_start_date", "registration_date"),
                         "Tuberculosis preventive therapy")
pct_all <- bind_rows(pct_cotri, pct_tpt)
addWorksheet(wb, "pct_records_initiation")
writeData(wb, "pct_records_initiation", as.data.frame(pct_all))


###############################################################################
# Section 4. Drop-all-influential-quarters sensitivity                        #
###############################################################################

cat("\n", strrep("=", 75), "\n", sep = "")
cat("Section 4. Drop-all-influential-quarters sensitivity\n")
cat(strrep("=", 75), "\n\n", sep = "")

fit_gls_war <- function(qdata) {
  fit <- tryCatch(
    nlme::gls(GLS_FORMULA, data = qdata,
              correlation = nlme::corAR1(form = ~ time_index),
              method = "REML"),
    error = function(e) list(err = conditionMessage(e))
  )
  if (!is.null(fit$err))
    return(data.frame(coef=NA, ci_lower=NA, ci_upper=NA, p_value=NA,
                      note=fit$err))
  s <- summary(fit)$tTable["level_war", ]
  data.frame(coef = s["Value"],
             ci_lower = s["Value"] - 1.96 * s["Std.Error"],
             ci_upper = s["Value"] + 1.96 * s["Std.Error"],
             p_value = s["p-value"], note = "")
}

influential_quarters <- list(
  cotri = list(years = c(2019, 2020, 2020),       qtrs = c(4, 1, 4)),
  tpt   = list(years = c(2018, 2019, 2020, 2020), qtrs = c(3, 4, 1, 4))
)

drop_and_refit <- function(qdata, inf_q, label) {
  drop_idx <- which(paste0(qdata$year, "_", qdata$quarter) %in%
                      paste0(inf_q$years, "_", inf_q$qtrs))
  cat(sprintf("  %s: dropping %d quarters\n", label, length(drop_idx)))
  qd <- qdata[-drop_idx, ]
  gls_full <- fit_gls_war(qdata); gls_d <- fit_gls_war(qd)
  nb_full  <- fit_nb_war(qdata);  nb_d  <- fit_nb_war(qd)
  cat(sprintf("    GLS coef:  full=%.2f  dropped=%.2f\n",
              gls_full$coef, gls_d$coef))
  cat(sprintf("    NB IRR:    full=%.4f  dropped=%.4f\n",
              nb_full$IRR,   nb_d$IRR))
  rbind(
    data.frame(outcome=label, fit="primary (full series)",
               n_quarters=nrow(qdata), model="GLS-AR1",
               coef_or_IRR=gls_full$coef,
               lower=gls_full$ci_lower, upper=gls_full$ci_upper,
               p_value=gls_full$p_value),
    data.frame(outcome=label, fit="drop influential quarters",
               n_quarters=nrow(qd), model="GLS-AR1",
               coef_or_IRR=gls_d$coef,
               lower=gls_d$ci_lower, upper=gls_d$ci_upper,
               p_value=gls_d$p_value),
    data.frame(outcome=label, fit="primary (full series)",
               n_quarters=nrow(qdata), model="NB IRR",
               coef_or_IRR=nb_full$IRR,
               lower=nb_full$IRR_lower, upper=nb_full$IRR_upper,
               p_value=nb_full$p_value),
    data.frame(outcome=label, fit="drop influential quarters",
               n_quarters=nrow(qd), model="NB IRR",
               coef_or_IRR=nb_d$IRR,
               lower=nb_d$IRR_lower, upper=nb_d$IRR_upper,
               p_value=nb_d$p_value)
  )
}

drop_combined <- rbind(
  drop_and_refit(cotri_q, influential_quarters$cotri, "Cotrimoxazole"),
  drop_and_refit(tpt_q,   influential_quarters$tpt,
                 "Tuberculosis preventive therapy")
)
addWorksheet(wb, "drop_influential_quarters")
writeData(wb, "drop_influential_quarters", drop_combined)


###############################################################################
# Section 5. Table 1 cohort: buggy reproduction and corrected version         #
###############################################################################

cat("\n", strrep("=", 75), "\n", sep = "")
cat("Section 5. Table 1 cohort with bug fix\n")
cat(strrep("=", 75), "\n\n", sep = "")

# Buggy reproduction: V13's case_when fallthrough puts NA periods in Post-war.
assign_period_buggy <- function(date_vec) {
  yr <- lubridate::year(date_vec)
  factor(case_when(
    yr <= 2019 ~ "Pre-war",
    yr == 2020 ~ "COVID-19",
    yr <= 2022 ~ "War",
    TRUE       ~ "Post-war"
  ), levels = c("Pre-war", "COVID-19", "War", "Post-war"))
}

make_table1 <- function(d, date_col, period_fn, label) {
  d <- d %>%
    filter(!is.na(patient_id) & patient_id != "" & patient_id != "NA")
  d$.period_this <- period_fn(d[[date_col]])
  d_first <- d %>%
    arrange(patient_id, .data[[date_col]]) %>%
    group_by(patient_id) %>% slice(1) %>% ungroup() %>%
    filter(!is.na(.period_this))
  d_first %>%
    group_by(.period_this) %>%
    group_modify(~ tibble(
      Variable = c("N", "Age, mean (SD)", "Age, median (IQR)",
                   "Female, n (%)", "Male, n (%)"),
      Value = c(
        as.character(nrow(.x)),
        sprintf("%.1f (%.1f)",
                mean(.x$age, na.rm = TRUE),
                sd(.x$age, na.rm = TRUE)),
        sprintf("%.1f (%.1f to %.1f)",
                median(.x$age, na.rm = TRUE),
                quantile(.x$age, 0.25, na.rm = TRUE),
                quantile(.x$age, 0.75, na.rm = TRUE)),
        sprintf("%d (%.1f%%)",
                sum(.x$sex == "Female", na.rm = TRUE),
                100 * sum(.x$sex == "Female", na.rm = TRUE) / nrow(.x)),
        sprintf("%d (%.1f%%)",
                sum(.x$sex == "Male", na.rm = TRUE),
                100 * sum(.x$sex == "Male", na.rm = TRUE) / nrow(.x))
      )
    )) %>%
    ungroup() %>%
    rename(period = .period_this) %>%
    pivot_wider(names_from = period, values_from = Value) %>%
    mutate(Outcome = label) %>%
    select(Outcome, Variable,
           any_of(c("Pre-war", "COVID-19", "War", "Post-war")))
}

t1_cotri_buggy <- make_table1(cotri, "cotri_start_date",
                              assign_period_buggy, "Cotrimoxazole")
t1_tpt_buggy   <- make_table1(tpt,   "tpt_start_date",
                              assign_period_buggy,
                              "Tuberculosis preventive therapy")
t1_cotri_fixed <- make_table1(cotri, "cotri_start_date",
                              assign_period_fixed, "Cotrimoxazole")
t1_tpt_fixed   <- make_table1(tpt,   "tpt_start_date",
                              assign_period_fixed,
                              "Tuberculosis preventive therapy")

cat("V13 buggy reproduction (case_when fallthrough places non-initiators in Post-war):\n")
print(t1_cotri_buggy); cat("\n")
print(t1_tpt_buggy);   cat("\n")
cat("Corrected version (non-initiators excluded from named periods):\n")
print(t1_cotri_fixed); cat("\n")
print(t1_tpt_fixed);   cat("\n")

addWorksheet(wb, "table1_cotri_buggy")
writeData(wb, "table1_cotri_buggy", as.data.frame(t1_cotri_buggy))
addWorksheet(wb, "table1_tpt_buggy")
writeData(wb, "table1_tpt_buggy", as.data.frame(t1_tpt_buggy))
addWorksheet(wb, "table1_cotri_corrected")
writeData(wb, "table1_cotri_corrected", as.data.frame(t1_cotri_fixed))
addWorksheet(wb, "table1_tpt_corrected")
writeData(wb, "table1_tpt_corrected", as.data.frame(t1_tpt_fixed))

# CSV copies for direct table use
write.csv(t1_cotri_buggy,
          file.path(out_dir, "table1_cotri_buggy_reproduction.csv"),
          row.names = FALSE, na = "")
write.csv(t1_tpt_buggy,
          file.path(out_dir, "table1_tpt_buggy_reproduction.csv"),
          row.names = FALSE, na = "")
write.csv(t1_cotri_fixed,
          file.path(out_dir, "table1_cotri_corrected.csv"),
          row.names = FALSE, na = "")
write.csv(t1_tpt_fixed,
          file.path(out_dir, "table1_tpt_corrected.csv"),
          row.names = FALSE, na = "")


###############################################################################
# Section 6. Save workbook and session info                                   #
###############################################################################

si <- data.frame(
  field = c("R version", "platform", "run_finished", "seed", "data_dir",
            paste0("pkg_", pkgs)),
  value = c(R.version.string, R.version$platform, format(Sys.time()), "42",
            data_dir,
            sapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)),
                                              error = function(e) "not loaded")))
)
addWorksheet(wb, "session_info")
writeData(wb, "session_info", si)
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

cat("\n", strrep("=", 75), "\n", sep = "")
cat("Reanalysis complete.\n")
cat(strrep("=", 75), "\n\n", sep = "")
cat("Workbook: ", xlsx_path, "\n")
cat("Log:      ", log_path, "\n")
cat("CSV tables in ", out_dir, "\n")
cat("Run finished:", format(Sys.time()), "\n")




###############################################################################
###############################################################################
###                                                                         ###
###             PART 3: V13.1 PROPER-SURVIVAL COX (STANDALONE)              ###
###                                                                         ###
###############################################################################
###############################################################################

###############################################################################
# ITS V13.1 PROPER SURVIVAL ANALYSIS, STANDALONE
#
# Manuscript: "Collapse of preventive HIV care during the Tigray war: a
#             20-year interrupted time series"
#
# Author:  Hafte Kahsay Kebede
# License: Apache-2.0
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Runs a full proper-survival Cox analysis that complements the
# ever-initiators-only specification in the V13 main pipeline.
#
# This script is standalone. It does not require any other script to
# have run first. It reads V13's raw inputs from data_dir and builds
# the event-level cotrimoxazole and tuberculosis preventive therapy
# datasets using V13's exact column-name conventions plus the V13.1
# demographic carry-forward correction.
#
# Cohort definition:
#   For each unique patient with a non-missing time origin (date of HIV
#   confirmation for cotrimoxazole; ART start date for tuberculosis
#   preventive therapy), one row is created.
#     - event = 1 if the patient initiated the preventive drug during
#       the study window; time = days from time origin to initiation
#     - event = 0 if the patient never initiated; time = days from time
#       origin to study end (administrative censoring at 31 May 2025)
#
# Exposure: the period in which the patient entered the cohort (the
# period of their time origin). This answers the question: did patients
# who entered care during the war get initiated at the same rate as
# patients who entered care before the war?
#
# IMPORTANT INTERPRETATION NOTE
# -----------------------------
# The proper-survival hazard ratios estimate the rate of initiation
# among patients who reached care, not the volume of care delivered.
# They can move in the opposite direction from the interrupted time
# series counts. The time series measures total events per quarter
# across the whole catchment; the Cox HRs measure per-patient timing
# conditional on reaching care. A higher Cox HR during the war
# indicates that the smaller cohort of patients who did reach care
# during the war was initiated faster than the pre-war cohort, which
# is consistent with selection of sicker patients into facility care
# combined with shorter administrative follow-up. Both findings
# together describe an access-driven collapse: total volume dropped,
# while per-patient processing within facilities was preserved or
# accelerated.
#
# OUTPUT
# ------
#   data_dir/Final_ITS_V13_Results/reanalysis/
#     ITS_V13_PROPER_SURVIVAL.xlsx
#       cox_results
#       schoenfeld_tests
#       km_summary
#       table_S10_ready
#       session_info
#     table_S10_proper_survival_ready.csv
#     proper_survival_console.txt
###############################################################################


###############################################################################
# Section 0. Setup, paths, packages                                           #
###############################################################################

data_dir       <- "D:/PhD project Data/New data set/FINAL_ANALYSIS_ALL/HIV_care _outcomes/TB/incidence/tpt/prophylactic"
ctx_file_stems <- c("crtEthiopiaARTvisit_1", "crtEthiopiaARTvisit_2")
tpt_file_stems <- c("vcrtEthiopiaARTVisit_TB_All")

config <- list(
  study_start   = as.Date("2005-01-01"),
  study_end     = as.Date("2025-05-31"),
  war_start     = as.Date("2021-01-01"),
  postwar_start = as.Date("2023-01-01")
)

out_dir <- file.path(data_dir, "Final_ITS_V13_Results", "reanalysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
xlsx_path <- file.path(out_dir, "ITS_V13_PROPER_SURVIVAL.xlsx")
log_path  <- file.path(out_dir, "proper_survival_console.txt")
csv_path  <- file.path(out_dir, "table_S10_proper_survival_ready.csv")

pkgs <- c("dplyr", "tidyr", "lubridate", "data.table", "tools",
          "survival", "openxlsx", "readxl")
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cran.r-project.org")
}
invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE,
                 warn.conflicts = FALSE))

set.seed(42)

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)

cat("ITS V13.1 proper-survival Cox analysis (standalone). Run started:",
    format(Sys.time()), "\n")
cat(strrep("=", 75), "\n\n")
cat("V13 data_dir:    ", data_dir, "\n")
cat("Output directory:", out_dir, "\n\n")

wb <- createWorkbook()


###############################################################################
# Section 0.5. Data construction                                              #
#                                                                             #
# Uses workspace objects if cotri and tpt2 (or cotri and tpt with > 0         #
# initiations) are present. Otherwise reads V13 raw inputs and rebuilds.      #
###############################################################################

cat(strrep("=", 75), "\n", sep = "")
cat("Section 0.5. Data construction\n")
cat(strrep("=", 75), "\n\n", sep = "")

# --- V13 helpers (verbatim) -------------------------------------------------

parse_date_flex <- function(x) {
  if (inherits(x, "Date")) {
    res <- x
  } else if (inherits(x, "POSIXt")) {
    res <- as.Date(x)
  } else {
    s <- as.character(x); s[s == "" | s == "NA"] <- NA
    res <- as.Date(rep(NA, length(s)))
    fmts <- c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%Y/%m/%d",
              "%d-%m-%Y", "%m-%d-%Y", "%d.%m.%Y")
    for (f in fmts) {
      idx <- is.na(res) & !is.na(s)
      if (!any(idx)) break
      res[idx] <- as.Date(s[idx], format = f)
    }
    idx <- is.na(res) & !is.na(s) & grepl("^[0-9]+(\\.[0-9]+)?$", s)
    if (any(idx)) {
      n <- suppressWarnings(as.numeric(s[idx]))
      res[idx] <- as.Date(n, origin = "1899-12-30")
    }
  }
  res[!is.na(res) & res <= as.Date("1900-12-31")] <- NA
  res
}

read_artvisit <- function(stem, dir) {
  candidates <- c(file.path(dir, paste0(stem, ".csv")),
                  file.path(dir, paste0(stem, ".xlsx")),
                  file.path(dir, paste0(stem, ".xls")),
                  file.path(dir, paste0(stem, ".rds")))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop(sprintf("Could not find '%s' in:\n  %s", stem, dir))
  cat(sprintf("  Reading: %s\n", basename(hit)))
  flush.console()
  ext <- tolower(tools::file_ext(hit))
  df <- switch(
    ext,
    csv  = data.table::fread(hit, na.strings = c("", "NA", "N/A"),
                             showProgress = FALSE,
                             colClasses = "character") |> as.data.frame(),
    xlsx = as.data.frame(readxl::read_excel(hit, guess_max = 100000,
                                            col_types = "text")),
    xls  = as.data.frame(readxl::read_excel(hit, guess_max = 100000,
                                            col_types = "text")),
    rds  = as.data.frame(readRDS(hit))
  )
  df[] <- lapply(df, function(col)
    if (is.character(col)) col else as.character(col))
  cat(sprintf("    %s rows, %d cols\n",
              format(nrow(df), big.mark = ","), ncol(df)))
  df
}

get_col <- function(df, target) {
  hits <- which(tolower(names(df)) == tolower(target))
  if (length(hits) == 0) return(NA_character_)
  names(df)[hits[1]]
}

.first_nonNA <- function(x) {
  v <- x[!is.na(x)]
  if (length(v) == 0) return(x[NA_integer_][seq_along(x)])
  out <- x
  out[] <- v[1]
  out
}

add_period_year_quarter <- function(df, date_col) {
  df %>% mutate(
    year = lubridate::year(.data[[date_col]]),
    quarter = lubridate::quarter(.data[[date_col]]),
    month = lubridate::month(.data[[date_col]]),
    period = factor(case_when(
      year <= 2019 ~ "Pre-war",   year == 2020 ~ "COVID-19",
      year <= 2022 ~ "War",       TRUE         ~ "Post-war"),
      levels = c("Pre-war", "COVID-19", "War", "Post-war"))
  )
}

assign_period_fixed <- function(date_vec) {
  yr <- lubridate::year(date_vec)
  out <- case_when(
    is.na(yr)  ~ NA_character_,
    yr <= 2019 ~ "Pre-war",
    yr == 2020 ~ "COVID-19",
    yr <= 2022 ~ "War",
    yr >= 2023 ~ "Post-war"
  )
  factor(out, levels = c("Pre-war", "COVID-19", "War", "Post-war"))
}

# --- Decide whether to use workspace or rebuild -----------------------------

workspace_cotri_ok <- exists("cotri", envir = .GlobalEnv) &&
  is.data.frame(get("cotri", envir = .GlobalEnv)) &&
  sum(get("cotri", envir = .GlobalEnv)$initiated, na.rm = TRUE) > 1000

workspace_tpt_name <- NULL
if (exists("tpt2", envir = .GlobalEnv) &&
    is.data.frame(get("tpt2", envir = .GlobalEnv)) &&
    sum(get("tpt2", envir = .GlobalEnv)$initiated, na.rm = TRUE) > 1000) {
  workspace_tpt_name <- "tpt2"
} else if (exists("tpt", envir = .GlobalEnv) &&
           is.data.frame(get("tpt", envir = .GlobalEnv)) &&
           sum(get("tpt", envir = .GlobalEnv)$initiated, na.rm = TRUE) > 1000) {
  workspace_tpt_name <- "tpt"
}

# --- Build cotri --------------------------------------------------------

if (workspace_cotri_ok) {
  cotri <- get("cotri", envir = .GlobalEnv)
  cat(sprintf("Using cotri from workspace: %s rows, %s initiations\n",
              format(nrow(cotri), big.mark = ","),
              format(sum(cotri$initiated, na.rm = TRUE), big.mark = ",")))
} else {
  cat("Building cotri from raw input files.\n")
  ctx_raw_list <- lapply(ctx_file_stems, read_artvisit, dir = data_dir)
  ctx_raw_list <- lapply(ctx_raw_list, function(d) {
    names(d) <- tolower(names(d)); d
  })
  ctx_combined <- bind_rows(
    lapply(seq_along(ctx_raw_list), function(i) {
      d <- ctx_raw_list[[i]]; d$.source_file <- ctx_file_stems[i]; d
    })
  )
  ctx_combined <- ctx_combined %>%
    distinct(across(-.source_file), .keep_all = TRUE)
  
  c_pid       <- get_col(ctx_combined, "uniqueartnumber")
  c_sex       <- get_col(ctx_combined, "sex")
  c_dob       <- get_col(ctx_combined, "dateofbirth")
  c_hivconf   <- get_col(ctx_combined, "date_hiv_confirmed")
  c_ctx_start <- get_col(ctx_combined, "cotrimoxazolestartdate")
  c_ctx_stop  <- get_col(ctx_combined, "cortimoxazole_stop_date")
  
  cotri <- ctx_combined %>%
    mutate(
      patient_id = if (!is.na(c_pid))
        trimws(as.character(.data[[c_pid]])) else NA_character_,
      patient_id = ifelse(is.na(patient_id) | patient_id == "" |
                            patient_id == "NA", NA_character_, patient_id),
      sex_raw = if (!is.na(c_sex))
        toupper(trimws(.data[[c_sex]])) else NA_character_,
      sex = case_when(sex_raw == "F" ~ "Female",
                      sex_raw == "M" ~ "Male", TRUE ~ NA_character_),
      DateOfBirth = if (!is.na(c_dob))
        parse_date_flex(.data[[c_dob]]) else as.Date(NA),
      date_hiv_confirmed = if (!is.na(c_hivconf))
        parse_date_flex(.data[[c_hivconf]]) else as.Date(NA),
      CotrimoxazoleStartDate = if (!is.na(c_ctx_start))
        parse_date_flex(.data[[c_ctx_start]]) else as.Date(NA),
      cortimoxazole_stop_date = if (!is.na(c_ctx_stop))
        parse_date_flex(.data[[c_ctx_stop]]) else as.Date(NA)
    )
  
  # V13.1 demographic carry-forward
  has_pid <- !is.na(cotri$patient_id)
  pid_part <- cotri[has_pid, ] %>%
    group_by(patient_id) %>%
    mutate(across(c(sex_raw, sex, DateOfBirth, date_hiv_confirmed),
                  .first_nonNA)) %>%
    ungroup()
  cotri <- bind_rows(pid_part, cotri[!has_pid, ])
  
  cotri <- cotri %>%
    mutate(
      cotri_start_date = dplyr::coalesce(CotrimoxazoleStartDate,
                                         cortimoxazole_stop_date),
      initiated = !is.na(cotri_start_date) &
        cotri_start_date >= config$study_start &
        cotri_start_date <= config$study_end
    ) %>%
    mutate(cotri_start_date = ifelse(initiated, cotri_start_date, NA) %>%
             as.Date(origin = "1970-01-01")) %>%
    mutate(
      age = as.numeric(difftime(
        dplyr::coalesce(cotri_start_date, date_hiv_confirmed),
        DateOfBirth, units = "days")) / 365.25,
      age = ifelse(age < 0 | age > 110, NA_real_, age),
      time_to_init_days = as.numeric(cotri_start_date - date_hiv_confirmed),
      time_to_init_days = ifelse(time_to_init_days < 0, NA_real_,
                                 time_to_init_days),
      record_date = date_hiv_confirmed
    ) %>%
    add_period_year_quarter("cotri_start_date")
  
  cotri_has <- cotri %>%
    filter(!is.na(patient_id)) %>%
    distinct(patient_id, cotri_start_date, .keep_all = TRUE)
  cotri_no <- cotri %>% filter(is.na(patient_id))
  cotri <- bind_rows(cotri_has, cotri_no)
  cat(sprintf("  cotri built: %s rows, %s initiations, %s unique initiating patients\n",
              format(nrow(cotri), big.mark = ","),
              format(sum(cotri$initiated, na.rm = TRUE), big.mark = ","),
              format(n_distinct(cotri$patient_id[cotri$initiated]),
                     big.mark = ",")))
  assign("cotri", cotri, envir = .GlobalEnv)
}

# --- Build tpt with the CORRECT V13 column names ----------------------------

if (!is.null(workspace_tpt_name)) {
  tpt <- get(workspace_tpt_name, envir = .GlobalEnv)
  cat(sprintf("Using %s from workspace: %s rows, %s initiations\n",
              workspace_tpt_name,
              format(nrow(tpt), big.mark = ","),
              format(sum(tpt$initiated, na.rm = TRUE), big.mark = ",")))
} else {
  cat("Building tpt from raw input file.\n")
  tpt_raw_list <- lapply(tpt_file_stems, read_artvisit, dir = data_dir)
  tpt_raw_list <- lapply(tpt_raw_list, function(d) {
    names(d) <- tolower(names(d)); d
  })
  tpt_raw <- tpt_raw_list[[1]] %>% distinct(.keep_all = TRUE)
  
  # V13's actual TPT column names
  t_pid       <- get_col(tpt_raw, "uano")
  t_sex       <- get_col(tpt_raw, "sex")
  t_dob       <- get_col(tpt_raw, "dateofbirth")
  t_artstart  <- get_col(tpt_raw, "art_start_date")
  t_regdate   <- get_col(tpt_raw, "registration_date")
  t_tbtype    <- get_col(tpt_raw, "tb_prophylaxis_type")
  t_inh_start <- get_col(tpt_raw, "inhprophylaxis_started_date")
  t_inh_disc  <- get_col(tpt_raw, "inhprophylaxisdiscontinueddate")
  t_inh_compl <- get_col(tpt_raw, "inhprophylaxiscompleteddate")
  
  tpt <- tpt_raw %>%
    mutate(
      patient_id = if (!is.na(t_pid))
        trimws(as.character(.data[[t_pid]])) else NA_character_,
      patient_id = ifelse(is.na(patient_id) | patient_id == "" |
                            patient_id == "NA", NA_character_, patient_id),
      sex_raw = if (!is.na(t_sex))
        toupper(trimws(.data[[t_sex]])) else NA_character_,
      sex = case_when(sex_raw == "F" ~ "Female",
                      sex_raw == "M" ~ "Male", TRUE ~ NA_character_),
      DateOfBirth = if (!is.na(t_dob))
        parse_date_flex(.data[[t_dob]]) else as.Date(NA),
      art_start_date = if (!is.na(t_artstart))
        parse_date_flex(.data[[t_artstart]]) else as.Date(NA),
      registration_date = if (!is.na(t_regdate))
        parse_date_flex(.data[[t_regdate]]) else as.Date(NA),
      tb_prophylaxis_type = if (!is.na(t_tbtype))
        suppressWarnings(as.numeric(.data[[t_tbtype]])) else NA_real_,
      inhprophylaxis_started_date = if (!is.na(t_inh_start))
        parse_date_flex(.data[[t_inh_start]]) else as.Date(NA),
      InhprophylaxisDiscontinuedDate = if (!is.na(t_inh_disc))
        parse_date_flex(.data[[t_inh_disc]]) else as.Date(NA),
      InhprophylaxisCompletedDate = if (!is.na(t_inh_compl))
        parse_date_flex(.data[[t_inh_compl]]) else as.Date(NA)
    )
  
  # V13.1 demographic carry-forward
  has_pid <- !is.na(tpt$patient_id)
  pid_part <- tpt[has_pid, ] %>%
    group_by(patient_id) %>%
    mutate(across(c(sex_raw, sex, DateOfBirth, art_start_date,
                    registration_date), .first_nonNA)) %>%
    ungroup()
  tpt <- bind_rows(pid_part, tpt[!has_pid, ])
  
  tpt <- tpt %>%
    filter(is.na(tb_prophylaxis_type) | tb_prophylaxis_type == 1) %>%
    mutate(
      tpt_start_date = dplyr::coalesce(inhprophylaxis_started_date,
                                       InhprophylaxisDiscontinuedDate,
                                       InhprophylaxisCompletedDate),
      initiated = !is.na(tpt_start_date) &
        tpt_start_date >= config$study_start &
        tpt_start_date <= config$study_end &
        (is.na(art_start_date) | tpt_start_date >= art_start_date)
    ) %>%
    mutate(tpt_start_date = ifelse(initiated, tpt_start_date, NA) %>%
             as.Date(origin = "1970-01-01")) %>%
    mutate(
      age = as.numeric(difftime(
        dplyr::coalesce(registration_date, art_start_date),
        DateOfBirth, units = "days")) / 365.25,
      age = ifelse(age < 0 | age > 110, NA_real_, age),
      time_to_init_days = as.numeric(tpt_start_date - art_start_date),
      time_to_init_days = ifelse(time_to_init_days < 0, NA_real_,
                                 time_to_init_days),
      record_date = dplyr::coalesce(art_start_date, registration_date)
    ) %>%
    add_period_year_quarter("tpt_start_date")
  
  tpt_has <- tpt %>%
    filter(!is.na(patient_id)) %>%
    distinct(patient_id, tpt_start_date, .keep_all = TRUE)
  tpt_no <- tpt %>% filter(is.na(patient_id))
  tpt <- bind_rows(tpt_has, tpt_no)
  cat(sprintf("  tpt built: %s rows, %s initiations, %s unique initiating patients\n",
              format(nrow(tpt), big.mark = ","),
              format(sum(tpt$initiated, na.rm = TRUE), big.mark = ","),
              format(n_distinct(tpt$patient_id[tpt$initiated]),
                     big.mark = ",")))
  assign("tpt2", tpt, envir = .GlobalEnv)
}

cat("\nFinal datasets in scope:\n")
cat(sprintf("  cotri: %s rows, %s initiated, %s unique initiating patients\n",
            format(nrow(cotri), big.mark = ","),
            format(sum(cotri$initiated, na.rm = TRUE), big.mark = ","),
            format(n_distinct(cotri$patient_id[cotri$initiated]),
                   big.mark = ",")))
cat(sprintf("  tpt:   %s rows, %s initiated, %s unique initiating patients\n\n",
            format(nrow(tpt), big.mark = ","),
            format(sum(tpt$initiated, na.rm = TRUE), big.mark = ","),
            format(n_distinct(tpt$patient_id[tpt$initiated]),
                   big.mark = ",")))


###############################################################################
# Section 1. Build the proper-survival Cox dataset                            #
###############################################################################

cat(strrep("=", 75), "\n", sep = "")
cat("Section 1. Building proper-survival Cox datasets\n")
cat(strrep("=", 75), "\n\n", sep = "")

build_proper_cox_data <- function(event_data, time_origin_col,
                                  start_date_col, label) {
  cat(sprintf("Building proper-survival dataset for %s.\n", label))
  df <- event_data %>%
    filter(!is.na(patient_id) & patient_id != "" & patient_id != "NA")
  df <- df %>%
    arrange(patient_id, desc(initiated), .data[[start_date_col]]) %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup()
  df <- df %>%
    mutate(
      time_origin = .data[[time_origin_col]],
      end_date = case_when(
        initiated ~ .data[[start_date_col]],
        TRUE      ~ config$study_end
      ),
      time_days = as.numeric(end_date - time_origin),
      time_months = time_days / 30.44,
      event = as.numeric(initiated),
      period_origin = assign_period_fixed(time_origin)
    )
  before <- nrow(df)
  df <- df %>%
    filter(!is.na(time_months), !is.na(sex), !is.na(age),
           !is.na(period_origin),
           time_months >= 0, time_months <= 240)
  after <- nrow(df)
  cat(sprintf("  %s rows after filters (lost %s rows for missing time, sex, age, or period)\n",
              format(after, big.mark = ","),
              format(before - after, big.mark = ",")))
  cat(sprintf("  %s events; %s censored; event rate %.1f%%\n",
              format(sum(df$event == 1), big.mark = ","),
              format(sum(df$event == 0), big.mark = ","),
              100 * mean(df$event)))
  ev_by_per <- table(df$period_origin[df$event == 1])
  cat(sprintf("  Events by period: %s\n\n",
              paste(names(ev_by_per), ev_by_per, sep = "=",
                    collapse = ", ")))
  df
}

cox_cotri <- build_proper_cox_data(cotri, "date_hiv_confirmed",
                                   "cotri_start_date", "Cotrimoxazole")
cox_tpt   <- build_proper_cox_data(tpt,   "art_start_date",
                                   "tpt_start_date",
                                   "Tuberculosis preventive therapy")


###############################################################################
# Section 2. Cox models                                                       #
###############################################################################

cat(strrep("=", 75), "\n", sep = "")
cat("Section 2. Cox models\n")
cat(strrep("=", 75), "\n\n", sep = "")

format_hr <- function(hr, lo, hi, p) {
  sprintf("%.3f (%.3f to %.3f); p=%s",
          hr, lo, hi, formatC(p, format = "g", digits = 3))
}

fit_cox_models <- function(cox_df, label) {
  cat(sprintf("\nCox models for %s\n", label))
  cat(strrep("-", 50), "\n", sep = "")
  if (sum(cox_df$event) == 0) {
    cat("  No events. Skipping Cox fit.\n")
    return(list(rows = data.frame(), zph = NULL, fit_a = NULL, fit_s = NULL))
  }
  
  fit_u <- coxph(Surv(time_months, event) ~ period_origin,
                 data = cox_df, ties = "efron")
  u_smry <- summary(fit_u)
  cat(sprintf("Unadjusted Cox (period only). n=%s, events=%s\n",
              format(fit_u$n,      big.mark = ","),
              format(fit_u$nevent, big.mark = ",")))
  for (term in rownames(u_smry$conf.int)) {
    cat(sprintf("  %-24s HR %s\n",
                sub("period_origin", "", term),
                format_hr(u_smry$conf.int[term, "exp(coef)"],
                          u_smry$conf.int[term, "lower .95"],
                          u_smry$conf.int[term, "upper .95"],
                          u_smry$coefficients[term, "Pr(>|z|)"])))
  }
  
  fit_a <- coxph(Surv(time_months, event) ~ period_origin + sex + age,
                 data = cox_df, ties = "efron")
  a_smry <- summary(fit_a)
  cat(sprintf("\nAdjusted Cox (period + sex + age). n=%s, events=%s\n",
              format(fit_a$n,      big.mark = ","),
              format(fit_a$nevent, big.mark = ",")))
  for (term in rownames(a_smry$conf.int)) {
    cat(sprintf("  %-24s HR %s\n",
                sub("period_origin", "", term),
                format_hr(a_smry$conf.int[term, "exp(coef)"],
                          a_smry$conf.int[term, "lower .95"],
                          a_smry$conf.int[term, "upper .95"],
                          a_smry$coefficients[term, "Pr(>|z|)"])))
  }
  
  zph <- tryCatch(cox.zph(fit_a), error = function(e) NULL)
  ph_violated <- FALSE
  if (!is.null(zph)) {
    cat("\nSchoenfeld proportional hazards test:\n")
    print(zph$table, digits = 3)
    ph_violated <- isTRUE(zph$table["GLOBAL", "p"] < 0.05)
  }
  
  fit_s <- NULL
  if (ph_violated) {
    cat("\nProportional hazards violated. Refitting stratified by period.\n")
    fit_s <- coxph(Surv(time_months, event) ~ sex + age +
                     strata(period_origin),
                   data = cox_df, ties = "efron")
    s_smry <- summary(fit_s)
    cat(sprintf("Stratified Cox (strata = period). n=%s, events=%s\n",
                format(fit_s$n,      big.mark = ","),
                format(fit_s$nevent, big.mark = ",")))
    for (term in rownames(s_smry$conf.int)) {
      cat(sprintf("  %-24s HR %s\n", term,
                  format_hr(s_smry$conf.int[term, "exp(coef)"],
                            s_smry$conf.int[term, "lower .95"],
                            s_smry$conf.int[term, "upper .95"],
                            s_smry$coefficients[term, "Pr(>|z|)"])))
    }
  }
  
  rows <- list()
  for (term in rownames(u_smry$conf.int)) {
    rows[[length(rows) + 1]] <- data.frame(
      outcome = label, model = "Unadjusted",
      term = sub("period_origin", "", term),
      n_patients = fit_u$n, n_events = fit_u$nevent,
      HR = u_smry$conf.int[term, "exp(coef)"],
      HR_lower = u_smry$conf.int[term, "lower .95"],
      HR_upper = u_smry$conf.int[term, "upper .95"],
      p_value = u_smry$coefficients[term, "Pr(>|z|)"]
    )
  }
  for (term in rownames(a_smry$conf.int)) {
    rows[[length(rows) + 1]] <- data.frame(
      outcome = label, model = "Adjusted (period + sex + age)",
      term = sub("period_origin", "", term),
      n_patients = fit_a$n, n_events = fit_a$nevent,
      HR = a_smry$conf.int[term, "exp(coef)"],
      HR_lower = a_smry$conf.int[term, "lower .95"],
      HR_upper = a_smry$conf.int[term, "upper .95"],
      p_value = a_smry$coefficients[term, "Pr(>|z|)"]
    )
  }
  if (!is.null(fit_s)) {
    s_smry <- summary(fit_s)
    for (term in rownames(s_smry$conf.int)) {
      rows[[length(rows) + 1]] <- data.frame(
        outcome = label, model = "Stratified by period (sex + age)",
        term = term,
        n_patients = fit_s$n, n_events = fit_s$nevent,
        HR = s_smry$conf.int[term, "exp(coef)"],
        HR_lower = s_smry$conf.int[term, "lower .95"],
        HR_upper = s_smry$conf.int[term, "upper .95"],
        p_value = s_smry$coefficients[term, "Pr(>|z|)"]
      )
    }
  }
  
  zph_tbl <- if (!is.null(zph)) {
    as.data.frame(zph$table) %>%
      mutate(outcome = label, term = rownames(zph$table)) %>%
      select(outcome, term, everything())
  } else NULL
  
  list(rows = do.call(rbind, rows), zph = zph_tbl,
       fit_a = fit_a, fit_s = fit_s)
}

res_cotri <- fit_cox_models(cox_cotri, "Cotrimoxazole")
res_tpt   <- fit_cox_models(cox_tpt,   "Tuberculosis preventive therapy")


###############################################################################
# Section 3. Kaplan-Meier summary by period                                   #
###############################################################################

cat("\n", strrep("=", 75), "\n", sep = "")
cat("Section 3. Kaplan-Meier survival summary by period\n")
cat(strrep("=", 75), "\n\n", sep = "")

km_summary <- function(cox_df, label) {
  cat(sprintf("\n%s\n", label))
  cat(strrep("-", 50), "\n", sep = "")
  if (sum(cox_df$event) == 0) {
    cat("  No events. Skipping KM summary.\n")
    return(data.frame())
  }
  rows <- list()
  for (per in c("Pre-war", "COVID-19", "War", "Post-war")) {
    sub <- cox_df %>% filter(period_origin == per)
    if (nrow(sub) < 5 || sum(sub$event) == 0) next
    sf <- survfit(Surv(time_months, event) ~ 1, data = sub)
    qs <- quantile(sf, probs = c(0.25, 0.5, 0.75))$quantile
    s12 <- summary(sf, times = 12)$surv
    s24 <- summary(sf, times = 24)$surv
    s60 <- summary(sf, times = 60)$surv
    s12 <- if (length(s12) == 0) NA else 1 - s12
    s24 <- if (length(s24) == 0) NA else 1 - s24
    s60 <- if (length(s60) == 0) NA else 1 - s60
    cat(sprintf("  %-9s n=%5d events=%5d  median(months)=%s  init by 24m=%s%%\n",
                per, sf$n, sum(sub$event),
                ifelse(is.na(qs[2]), "  NA (>50%% censored)",
                       sprintf("%6.1f", qs[2])),
                ifelse(is.na(s24), "  NA", sprintf("%5.1f", 100 * s24))))
    rows[[length(rows) + 1]] <- data.frame(
      outcome = label, period = per,
      n_patients = sf$n, n_events = sum(sub$event),
      median_months_to_init = unname(qs[2]),
      q25_months = unname(qs[1]),
      q75_months = unname(qs[3]),
      pct_initiated_by_12m = 100 * s12,
      pct_initiated_by_24m = 100 * s24,
      pct_initiated_by_60m = 100 * s60
    )
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

km_cotri <- km_summary(cox_cotri, "Cotrimoxazole")
km_tpt   <- km_summary(cox_tpt,   "Tuberculosis preventive therapy")
km_all   <- if (nrow(km_cotri) + nrow(km_tpt) > 0)
  rbind(km_cotri, km_tpt) else data.frame()


###############################################################################
# Section 4. Ready-to-paste Table S10 block                                   #
###############################################################################

cat("\n", strrep("=", 75), "\n", sep = "")
cat("Section 4. Composing ready-to-paste Table S10 entries\n")
cat(strrep("=", 75), "\n\n", sep = "")

format_table_s10 <- function(rows, label) {
  if (nrow(rows) == 0) {
    return(data.frame(Outcome = label,
                      Specification = "Adjusted Cox (period + sex + age), proper survival with censoring at study end",
                      n_patients = "not estimable",
                      n_events = "not estimable",
                      HR_COVID_vs_Pre_war   = "not estimable",
                      HR_War_vs_Pre_war     = "not estimable",
                      HR_Post_war_vs_Pre_war = "not estimable"))
  }
  adj <- rows %>%
    filter(model == "Adjusted (period + sex + age)" &
             grepl("Pre-war|COVID-19|War|Post-war", term))
  out <- data.frame(
    Outcome = label,
    Specification = "Adjusted Cox (period + sex + age), proper survival with censoring at study end",
    n_patients = if (nrow(adj) > 0) format(adj$n_patients[1], big.mark = ",") else NA,
    n_events   = if (nrow(adj) > 0) format(adj$n_events[1],   big.mark = ",") else NA
  )
  for (per in c("COVID-19", "War", "Post-war")) {
    row <- adj %>% filter(term == per)
    val <- if (nrow(row) > 0) {
      sprintf("%.3f (%.3f to %.3f); p=%s",
              row$HR, row$HR_lower, row$HR_upper,
              formatC(row$p_value, format = "g", digits = 3))
    } else "not estimable"
    colname <- sprintf("HR_%s_vs_Pre_war",
                       gsub("[^A-Za-z0-9]", "_", per))
    out[[colname]] <- val
  }
  out
}

ts10_cotri <- format_table_s10(res_cotri$rows, "Cotrimoxazole")
ts10_tpt   <- format_table_s10(res_tpt$rows,   "Tuberculosis preventive therapy")
ts10 <- rbind(ts10_cotri, ts10_tpt)

cat("Table S10 ready-to-paste block:\n\n")
print(ts10, row.names = FALSE)

cat("\nNarrative phrasing for the supplementary text:\n\n")
for (out_lbl in c("Cotrimoxazole", "Tuberculosis preventive therapy")) {
  r <- if (out_lbl == "Cotrimoxazole") res_cotri$rows else res_tpt$rows
  if (nrow(r) == 0) next
  adj <- r %>% filter(model == "Adjusted (period + sex + age)")
  war_row <- adj %>% filter(term == "War")
  pw_row  <- adj %>% filter(term == "Post-war")
  if (nrow(war_row) == 1 && nrow(pw_row) == 1) {
    cat(sprintf(
      paste0("  %s: proper-survival analysis (n=%s patients with %s ",
             "initiations) showed a war-period adjusted hazard ratio of ",
             "%.3f (95%% CI %.3f to %.3f, p=%s) and a post-war adjusted ",
             "hazard ratio of %.3f (95%% CI %.3f to %.3f, p=%s), with ",
             "censoring at the study end for patients who did not initiate.\n\n"),
      out_lbl,
      format(war_row$n_patients, big.mark = ","),
      format(war_row$n_events,   big.mark = ","),
      war_row$HR, war_row$HR_lower, war_row$HR_upper,
      formatC(war_row$p_value, format = "g", digits = 3),
      pw_row$HR,  pw_row$HR_lower,  pw_row$HR_upper,
      formatC(pw_row$p_value,  format = "g", digits = 3)
    ))
  }
}

cat("INTERPRETATION NOTE:\n")
cat("  The proper-survival hazard ratios above are conditional on reaching\n")
cat("  care. They estimate per-patient timing of initiation among patients\n")
cat("  who reached a facility, not the total volume of care delivered. They\n")
cat("  can move in the opposite direction from the time-series counts. The\n")
cat("  interrupted time series measures total events per quarter across the\n")
cat("  catchment; the Cox HRs measure per-patient timing. A higher Cox HR\n")
cat("  during the war indicates that the smaller cohort of patients who did\n")
cat("  reach care during the war was processed faster than the pre-war\n")
cat("  cohort, consistent with selection of sicker patients and shorter\n")
cat("  administrative follow-up. Both findings together describe an access-\n")
cat("  driven collapse: total volume dropped while per-patient processing\n")
cat("  within the surviving facilities was preserved or accelerated.\n\n")


###############################################################################
# Section 5. Save outputs                                                     #
###############################################################################

all_rows <- if (nrow(res_cotri$rows) + nrow(res_tpt$rows) > 0)
  rbind(res_cotri$rows, res_tpt$rows) else data.frame()
zph_rows <- rbind(res_cotri$zph, res_tpt$zph)

addWorksheet(wb, "cox_results")
writeData(wb, "cox_results", all_rows)
addWorksheet(wb, "schoenfeld_tests")
writeData(wb, "schoenfeld_tests",
          if (!is.null(zph_rows)) zph_rows else data.frame())
addWorksheet(wb, "km_summary")
writeData(wb, "km_summary", km_all)
addWorksheet(wb, "table_S10_ready")
writeData(wb, "table_S10_ready", ts10)

si <- data.frame(
  field = c("R version", "platform", "run_finished", "seed", "data_dir",
            "cohort_definition",
            paste0("pkg_", c("survival", "dplyr", "tidyr", "lubridate",
                             "openxlsx"))),
  value = c(R.version.string, R.version$platform, format(Sys.time()), "42",
            data_dir,
            "One row per patient. Time origin: date_hiv_confirmed (cotri), art_start_date (tpt). Censoring: study_end if not initiated.",
            sapply(c("survival","dplyr","tidyr","lubridate","openxlsx"),
                   function(p) as.character(packageVersion(p))))
)
addWorksheet(wb, "session_info")
writeData(wb, "session_info", si)

write.csv(ts10, csv_path, row.names = FALSE, na = "")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

cat(strrep("=", 75), "\n", sep = "")
cat("Proper-survival analysis complete.\n")
cat(strrep("=", 75), "\n\n", sep = "")
cat("Workbook: ", xlsx_path, "\n")
cat("CSV:      ", csv_path, "\n")
cat("Log:      ", log_path, "\n")
cat("Run finished:", format(Sys.time()), "\n")
