# CroptimizR calibration of the FACE decks, the AgMIP procedure of Yaron (calibrate_croptimizr.R, 29 Oct 2026) applied
# to the 9 phenology control decks and the 4 cultivar decks. Criterion crit_ols, optimizer Brent, one parameter and one
# target per step.
#   phase "phen" (per control deck): step 1 Anthesis -> observed anthesis day; step 2 GrainFill_TDU -> observed maturity
#     day. Bounds 0.7-1.3 x deck Anthesis, 0.5-1.5 x deck GrainFill_TDU (as AgMIP).
#   phase "cult" (per cultivar deck, after phase phen was written to the decks): [Blade_L_max_cm -> observed peak LAI,
#     bounds = literature range] then fruiting_efficiency -> observed grain yield, bounds 50-200 (model range, as AgMIP).
#     Decks listed in NO_FE keep their measured fruiting_efficiency (only Blade_L_max_cm is calibrated); decks in KN_DECKS
#     calibrate fruiting_efficiency to observed kernel number; ONLY_DECKS limits the phase to the listed decks.
# Usage: Rscript calibrate_face.R <phase> <targets.csv> <out_dir>
suppressPackageStartupMessages(library(CroptimizR))
args <- commandArgs(trailingOnly = TRUE)
phase <- args[1]; tg <- read.csv(args[2], stringsAsFactors = FALSE); out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
runner <- file.path(Sys.getenv("FCAL_HOME"), "face_run.py")
KN_DECKS <- strsplit(Sys.getenv("KN_DECKS"), ",")[[1]]
NO_FE <- strsplit(Sys.getenv("NO_FE"), ",")[[1]]   # decks whose fruiting efficiency is MEASURED (not calibrated)
ONLY <- Sys.getenv("ONLY_DECKS")
LEAF <- list(WDF25_Ruta_aCO2_WW = c(17, 26), WDF25_Zahir_aCO2_WW = c(17, 26), MARICOPA_906 = c(16, 25))

wrapper <- function(param_values = NULL, situation, model_options, var = NULL, ...) {
  res <- lapply(situation, function(sit) {
    pv <- if (is.data.frame(param_values)) as.list(param_values[1, setdiff(names(param_values), "situation")]) else
      as.list(param_values)
    pj <- tempfile(fileext = ".json", tmpdir = model_options$work); oc <- tempfile(fileext = ".csv", tmpdir = model_options$work)
    writeLines(paste0("{", paste(sprintf('"%s": %.10g', names(pv), unlist(pv)), collapse = ", "), "}"), pj)
    st <- system2("python3", c("-I", shQuote(runner), "--deck", sit, "--params", paste0("@", pj), "--out", shQuote(oc),
                               "--work", shQuote(model_options$work), "--log", shQuote(model_options$log)))
    unlink(pj)
    if (st != 0 || !file.exists(oc)) return(NULL)
    df <- read.csv(oc, stringsAsFactors = FALSE); unlink(oc)
    df$Date <- as.POSIXct(df$Date, format = "%Y-%m-%d", tz = "UTC"); df
  })
  names(res) <- situation
  if (any(vapply(res, is.null, logical(1)))) return(list(sim_list = NULL, error = TRUE))
  if (!is.null(var)) res <- lapply(res, function(df) df[, c("Date", intersect(var, names(df))), drop = FALSE])
  attr(res, "class") <- "cropr_simulation"
  list(sim_list = res, error = FALSE)
}
brent <- function(tol) list(method = "Brent", nb_rep = 1, control = list(reltol = tol))

calibrate <- function(t) {
  d <- t$deck; sdir <- file.path(out_dir, d); dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  opts <- list(work = file.path(sdir, "tmp"), log = file.path(sdir, "evals.csv")); dir.create(opts$work, showWarnings = FALSE)
  unlink(opts$log)
  date <- as.POSIXct(t$start_date, tz = "UTC")
  if (phase == "phen") {
    obs <- setNames(list(data.frame(Date = date, anthesis_dayf = t$anthesis_dayf, maturity_dayf = t$maturity_dayf)), d)
    pi <- list(lb = c(Anthesis = 0.7 * t$A0, GrainFill_TDU = 0.5 * t$G0), ub = c(Anthesis = 1.3 * t$A0, GrainFill_TDU = 1.5 * t$G0),
               init_values = data.frame(Anthesis = t$A0, GrainFill_TDU = t$G0))
    step <- list(list(major_param = "Anthesis", obs_var = "anthesis_dayf", optim_options = brent(0.5)),
                 list(major_param = "GrainFill_TDU", obs_var = "maturity_dayf", optim_options = brent(0.5)))
  } else {
    obs <- setNames(list(data.frame(Date = date, yield_final = t$yield_final, lai_peak = t$lai_peak, kernel_number = t$kernel_number)), d)
    fe_var <- if (d %in% KN_DECKS) "kernel_number" else "yield_final"
    fe0 <- min(max(t$FE0, 50), 200)
    lb <- c(fruiting_efficiency = 50); ub <- c(fruiting_efficiency = 200); iv <- data.frame(fruiting_efficiency = fe0)
    step <- list(list(major_param = "fruiting_efficiency", obs_var = fe_var, optim_options = brent(0.02)))
    if (d %in% NO_FE) { lb <- c(); ub <- c(); iv <- NULL; step <- list() }
    if (!is.null(LEAF[[d]])) {
      b <- LEAF[[d]]; l0 <- min(max(t$L0, b[1]), b[2])
      lb <- c(Blade_L_max_cm = b[1], lb); ub <- c(Blade_L_max_cm = b[2], ub)
      iv <- if (is.null(iv)) data.frame(Blade_L_max_cm = l0) else cbind(data.frame(Blade_L_max_cm = l0), iv)
      step <- c(list(list(major_param = "Blade_L_max_cm", obs_var = "lai_peak", optim_options = brent(0.02))), step)
    }
    pi <- list(lb = lb, ub = ub, init_values = iv)
  }
  res <- estim_param(obs_list = obs, model_function = wrapper, model_options = opts, crit_function = crit_ols,
                     optim_method = "optim", optim_options = brent(0.5), param_info = pi, step = step,
                     out_dir = file.path(sdir, "croptimizr"))
  unlink(opts$work, recursive = TRUE)
  fv <- res$final_values
  cat(d, "done:", paste(names(fv), signif(fv, 7), collapse = ", "), "\n")
  data.frame(deck = d, t(fv), runs = nrow(read.csv(opts$log)))
}
rows <- tg[tg$kind == phase, ]
if (nzchar(ONLY)) rows <- rows[rows$deck %in% strsplit(ONLY, ",")[[1]], ]
res <- parallel::mclapply(seq_len(nrow(rows)), function(i) tryCatch(calibrate(rows[i, ]), error = function(e) conditionMessage(e)),
                          mc.cores = nrow(rows), mc.preschedule = FALSE)
bad <- vapply(res, function(r) !is.data.frame(r), logical(1))
if (any(bad)) for (i in which(bad)) cat("FAILED", rows$deck[i], ":", format(res[[i]]), "\n")
tab <- do.call(plyr::rbind.fill, res[!bad])
write.csv(tab, file.path(out_dir, paste0("calibrated_", phase, ".csv")), row.names = FALSE)
print(tab)
