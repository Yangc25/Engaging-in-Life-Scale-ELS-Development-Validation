# =========================================================
# Engaging in Life Scale (ELS) – Development & Validation
# Pipeline for OSF
# Stages: Stage 1 (Item analysis + EFA) →
#         Stage 2 (CFA + MIRT + validity) →
#         Stage 3 (ACO short form → CFA/MIRT + validity)
# =========================================================

options(stringsAsFactors = FALSE)

# -------------------- Reproducibility --------------------
SEED_SPLIT  <- 20250822   # train/test split
SEED_ACO    <- 20250822   # ACO short-form search
set.seed(SEED_SPLIT)

# -------------------- Working directory -----------------
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("--file=", "", file_arg))
    setwd(dirname(script_path))
  }
}

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
sink("outputs/sessionInfo.txt"); print(sessionInfo()); sink()

# -------------------- Dependencies ----------------------
pkg_needed <- c(
  "readxl","psych","GPArotation","mirt","lavaan","semPlot","dplyr",
  "semTools","effsize","purrr","ggplot2","stringr"
)
missing <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(missing)) {
  stop(sprintf(
    "Missing packages: %s. Install with: install.packages(c(%s))",
    paste(missing, collapse = ", "),
    paste(sprintf('"%s"', missing), collapse = ", ")
  ))
}
suppressPackageStartupMessages({
  library(readxl); library(psych); library(GPArotation); library(mirt)
  library(lavaan); library(semPlot); library(dplyr); library(semTools)
  library(effsize); library(purrr); library(ggplot2); library(stringr)
})

# -------------------- Helpers ---------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
ensure_parent <- function(path){ d <- dirname(path); if (!dir.exists(d)) dir.create(d, recursive = TRUE) }
to_numeric_safe <- function(x){
  if (is.factor(x)) x <- as.character(x)
  zh <- c("从不"=1,"很少"=2,"有时"=3,"经常"=4,"总是"=5)
  if (is.character(x)) {
    if (all(na.omit(x) %in% names(zh))) return(as.numeric(unname(zh[x])))
    if (all(grepl("^\\s*[0-9]+\\s*$", na.omit(x)))) return(as.numeric(trimws(x)))
  }
  suppressWarnings(as.numeric(x))
}
reverse_by_range <- function(v){
  if (all(is.na(v))) return(v)
  mn <- suppressWarnings(min(v, na.rm=TRUE)); mx <- suppressWarnings(max(v, na.rm=TRUE))
  if (!is.finite(mn) || !is.finite(mx) || mn == mx) return(v)
  ifelse(is.na(v), NA_real_, (mn + mx) - v)
}
check_cols    <- function(df, cols) intersect(cols, colnames(df))
row_mean_safe <- function(df, cols){ cols <- check_cols(df, cols); if (!length(cols)) return(rep(NA_real_, nrow(df))); suppressWarnings(rowMeans(df[, cols, drop=FALSE], na.rm=TRUE)) }
row_sum_safe  <- function(df, cols){ cols <- check_cols(df, cols); if (!length(cols)) return(rep(NA_real_, nrow(df))); suppressWarnings(rowSums(df[, cols, drop=FALSE], na.rm=TRUE)) }
fmt_p         <- function(p) ifelse(is.na(p), NA_character_, ifelse(p < .001, "< .001", sprintf("= %.3f", round(p, 3))))
pair_row <- function(df, v1, v2, method="pearson"){
  x <- df[[v1]]; y <- df[[v2]]
  if (is.null(x) || is.null(y)) return(data.frame(var1=v1,var2=v2,r=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p=NA_real_))
  ct <- try(suppressWarnings(cor.test(x,y,method=method,use="pairwise.complete.obs")), silent=TRUE)
  if (inherits(ct,"try-error")) return(data.frame(var1=v1,var2=v2,r=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p=NA_real_))
  ci <- try(as.numeric(ct$conf.int), silent=TRUE); if (inherits(ci,"try-error")||length(ci)!=2) ci <- c(NA_real_,NA_real_)
  data.frame(var1=v1,var2=v2,r=as.numeric(ct$estimate),ci_low=ci[1],ci_high=ci[2],p=as.numeric(ct$p.value))
}
ordify_like_short <- function(df){
  out <- df
  out[] <- lapply(out, function(x){
    if (is.ordered(x)) return(x)
    if (is.factor(x))  return(ordered(x))
    if (is.numeric(x)) return(ordered(as.integer(x)))
    ordered(x)
  })
  out
}
cr_ave_table <- function(fit){
  std_loadings <- lavaan::inspect(fit, "std")$lambda
  fac_names <- colnames(std_loadings)
  out_list <- lapply(fac_names, function(fac) {
    loadings <- std_loadings[, fac]; loadings <- loadings[loadings != 0]
    cr  <- sum(loadings)^2 / (sum(loadings)^2 + sum(1 - loadings^2))
    ave <- sum(loadings^2) / length(loadings)
    data.frame(Factor = fac, CR = round(cr,3), AVE = round(ave,3), N_items = length(loadings))
  })
  do.call(rbind, out_list)
}
htmt_safely <- function(fit, model_syntax = NULL, data = NULL) {
  lv_cor <- try(lavaan::lavInspect(fit, "cor.lv"), silent=TRUE)
  if (!inherits(lv_cor,"try-error") && is.matrix(lv_cor)) return(as.data.frame(lv_cor))
  if (!is.null(model_syntax) && !is.null(data)) {
    data_num <- as.data.frame(lapply(data, function(x) as.numeric(as.character(x))))
    fit_ml <- lavaan::sem(model_syntax, data = data_num, std.lv = TRUE, estimator = "MLR")
    out <- try(semTools::htmt(fit_ml), silent = TRUE)
    if (!inherits(out,"try-error")) return(out)
  }
  warning("HTMT unavailable for this model."); NULL
}

# -------------------- Privacy ---------------------------
# Demographic columns to keep (whitelist).
demo_vars_canonical <- c(
  "Gender","Age","Education","Ethnic","Occupation",
  "Social status","Emotional disorders","Psychotherapy"
)
# Optional name mapping: list("sex"="Gender", ...)
demo_vars_alias <- list()

drop_identifier_cols <- function(df) {
  if (length(demo_vars_alias)) {
    for (nm in names(demo_vars_alias)) {
      if (nm %in% colnames(df)) names(df)[names(df)==nm] <- demo_vars_alias[[nm]]
    }
  }
  canonical_present <- intersect(demo_vars_canonical, colnames(df))
  keep_protect <- rep(FALSE, ncol(df)); names(keep_protect) <- colnames(df)
  keep_protect[canonical_present] <- TRUE
  patt <- paste(
    c("(^|_)id$", "name", "full.?name", "real.?name", "phone", "mobile",
      "email", "mail", "wechat", "qq", "whats", "address", "contact",
      "ip", "geo", "lat", "lon", "longitude", "latitude",
      "timestamp", "time_created", "created_at", "submitted_at", "uid",
      "national.?id", "passport", "student.?id", "employee.?id"),
    collapse = "|"
  )
  is_identifier <- grepl(patt, tolower(colnames(df)))
  drop_mask <- is_identifier & !keep_protect
  if (any(drop_mask)) {
    message("Dropping identifier-like columns: ", paste(colnames(df)[drop_mask], collapse = ", "))
  }
  df[, !drop_mask, drop = FALSE]
}

# -------------------- Analysis thresholds ----------------
CITC_MIN       <- 0.40
TTEST_P_MAX    <- 0.05
ALPHA_GAIN_MIN <- 0.01
LOAD_MIN       <- 0.32
CROSS_GAP_MIN  <- 0.20
H2_MIN         <- 0.20
N_FACTORS_EFA  <- 4

# -------------------- Data input ------------------------
infile <- "data/ELS_raw.xlsx"
if (!file.exists(infile)) stop("Input file not found: ", infile)

dat_raw <- readxl::read_xlsx(infile)

dir.create("outputs/descriptives", showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(Column = colnames(dat_raw)), "outputs/descriptives/columns_present.csv", row.names = FALSE)

dat_raw <- drop_identifier_cols(dat_raw)
demopres <- intersect(demo_vars_canonical, colnames(dat_raw))
write.csv(data.frame(DemographicColumn = demopres),
          "outputs/descriptives/demographics_columns.csv", row.names = FALSE)

# -------------------- Stage 1: Item analysis + EFA ------
dir.create("outputs/stage1", showWarnings = FALSE, recursive = TRUE)

item_names    <- paste0("ELS", 1:24)
reverse_items <- c("ELS1","ELS2","ELS3","ELS5","ELS7","ELS8","ELS9","ELS11",
                   "ELS13","ELS14","ELS15","ELS17","ELS19","ELS20","ELS21","ELS23")

dat <- dat_raw
for (it in intersect(reverse_items, colnames(dat))) {
  dat[[it]] <- 6 - to_numeric_safe(dat[[it]])
}
dat$total_score_all24 <- rowSums(data.frame(lapply(dat[, intersect(item_names, colnames(dat)), drop=FALSE], to_numeric_safe)), na.rm = TRUE)

n_total  <- nrow(dat)
idx_all  <- sample(seq_len(n_total))
half     <- floor(n_total/2)
train_idx <- idx_all[1:half]
test_idx  <- idx_all[(half+1):n_total]
dat_train <- dat[train_idx, , drop=FALSE]
dat_test  <- dat[test_idx,  , drop=FALSE]

# Demographic exports (Total/Train/Test)
demo_vars_raw <- intersect(demo_vars_canonical, colnames(dat))
if (length(demo_vars_raw)) {
  demo_total <- dat[,       demo_vars_raw, drop = FALSE]; demo_total$.set <- "Total"
  demo_train <- dat_train[, demo_vars_raw, drop = FALSE]; demo_train$.set <- "Train"
  demo_test  <- dat_test[,  demo_vars_raw, drop = FALSE]; demo_test$.set  <- "Test"
  demo <- dplyr::bind_rows(demo_total, demo_train, demo_test)
  
  chisq_or_fisher <- function(tab){
    if (any(tab < 5, na.rm = TRUE)) {
      list(p = tryCatch(fisher.test(tab)$p.value, error=function(e) NA_real_), method = "Fisher")
    } else {
      list(p = tryCatch(chisq.test(tab, correct = FALSE)$p.value, error=function(e) NA_real_), method = "Chi-square")
    }
  }
  cramers_v <- function(tab){
    k <- min(nrow(tab), ncol(tab))
    cs <- tryCatch(chisq.test(tab, correct = FALSE)$statistic, error=function(e) NA_real_)
    n  <- sum(tab)
    if (!is.finite(cs) || !is.finite(n) || k <= 1) return(NA_real_)
    as.numeric(sqrt(cs / (n * (k - 1))))
  }
  
  # Categorical
  cat_vars <- setdiff(demo_vars_raw, c("Age","Social status"))
  cat_rows <- list()
  for (v in cat_vars) {
    dv <- demo[, c(".set", v), drop = FALSE]
    names(dv)[2] <- "Level"
    dv$Level <- as.factor(dv$Level)
    tab <- with(dv, table(.set, Level, useNA = "no"))
    prop <- prop.table(tab, margin = 1) * 100
    
    mk_df <- function(i_set){
      if (!(i_set %in% rownames(tab))) return(data.frame(Level = colnames(tab), val = NA_character_))
      n   <- as.numeric(tab[i_set, ])
      pct <- as.numeric(prop[i_set, ])
      data.frame(Level = colnames(tab),
                 val   = ifelse(is.na(n), NA_character_, sprintf("%d (%.1f%%)", n, pct)))
    }
    w_total <- mk_df("Total"); names(w_total)[2] <- "Total"
    w_train <- mk_df("Train"); names(w_train)[2] <- "Train"
    w_test  <- mk_df("Test");  names(w_test)[2]  <- "Test"
    
    out <- Reduce(function(x,y) dplyr::left_join(x,y, by="Level"),
                  list(w_total, w_train, w_test))
    out$Variable <- v
    out <- out[, c("Variable","Level","Total","Train","Test")]
    
    tab_tt <- tab[intersect(rownames(tab), c("Train","Test")), , drop = FALSE]
    tst <- chisq_or_fisher(tab_tt); cv <- cramers_v(tab_tt)
    out$Test_method <- c(tst$method, rep(NA_character_, nrow(out)-1))
    out$p_value     <- c(round(tst$p, 4), rep(NA_real_, nrow(out)-1))
    out$CramersV    <- c(round(cv, 3),    rep(NA_real_, nrow(out)-1))
    
    cat_rows[[length(cat_rows)+1]] <- out
  }
  demo_cat <- if (length(cat_rows)) dplyr::bind_rows(cat_rows) else
    data.frame(Variable=character(),Level=character(),Total=character(),Train=character(),Test=character(),
               Test_method=character(),p_value=numeric(),CramersV=numeric(), stringsAsFactors = FALSE)
  write.csv(demo_cat, "outputs/descriptives/demo_categorical.csv", row.names = FALSE)
  
  # Numeric (Age, Social status)
  num_vars <- intersect(c("Age","Social status"), demo_vars_raw)
  num_rows <- list()
  for (v in num_vars) {
    x_total <- to_numeric_safe(demo_total[[v]])
    x_train <- to_numeric_safe(demo_train[[v]])
    x_test  <- to_numeric_safe(demo_test[[v]])
    m_sd <- function(x){ if (all(is.na(x))) return(NA_character_); sprintf("%.2f (%.2f)", mean(x,na.rm=TRUE), stats::sd(x,na.rm=TRUE)) }
    p_t <- tryCatch(t.test(x_train, x_test, var.equal = FALSE)$p.value, error=function(e) NA_real_)
    d_es <- tryCatch(as.numeric(effsize::cohen.d(x_train, x_test, hedges.correction = TRUE)$estimate), error=function(e) NA_real_)
    num_rows[[length(num_rows)+1]] <- data.frame(
      Variable = v, Total = m_sd(x_total), Train = m_sd(x_train), Test = m_sd(x_test),
      t_p_value = round(p_t, 4), Hedges_g  = round(d_es, 3), stringsAsFactors = FALSE
    )
  }
  demo_num <- if (length(num_rows)) dplyr::bind_rows(num_rows) else
    data.frame(Variable=character(),Total=character(),Train=character(),Test=character(),
               t_p_value=numeric(),Hedges_g=numeric(), stringsAsFactors = FALSE)
  write.csv(demo_num, "outputs/descriptives/demo_numeric.csv", row.names = FALSE)
  
  demo_all <- dplyr::bind_rows(
    dplyr::mutate(demo_cat, Type = "Categorical"),
    dplyr::mutate(demo_num, Type = "Numeric")
  )
  write.csv(demo_all, "outputs/descriptives/demo_all.csv", row.names = FALSE)
}

# Item analysis (train)
item_in_train <- intersect(item_names, colnames(dat_train))
dat_tr_items  <- data.frame(lapply(dat_train[, item_in_train, drop=FALSE], to_numeric_safe))

n_tr <- nrow(dat_tr_items)
top_n <- floor(n_tr * 0.27); bottom_n <- floor(n_tr * 0.27)
ord   <- order(dat_train$total_score_all24)
low_g <- dat_tr_items[ord[1:bottom_n], , drop=FALSE]
high_g<- dat_tr_items[ord[(n_tr - top_n + 1):n_tr], , drop=FALSE]

safe_t <- function(xh, xl){ tryCatch(t.test(xh, xl, var.equal = FALSE)$p.value, error = function(e) NA_real_) }
safe_t_stat <- function(xh, xl){ tryCatch(as.numeric(t.test(xh, xl, var.equal = FALSE)$statistic), error = function(e) NA_real_) }

t_results <- sapply(colnames(dat_tr_items), function(item) safe_t(high_g[[item]], low_g[[item]]))
t_stats   <- sapply(colnames(dat_tr_items), function(item) safe_t_stat(high_g[[item]], low_g[[item]]))
item_ttest <- data.frame(item = colnames(dat_tr_items), t_stat = as.numeric(t_stats), p_value = as.numeric(t_results))

citc <- sapply(colnames(dat_tr_items), function(item){
  total_wo_item <- rowSums(dat_tr_items[, setdiff(colnames(dat_tr_items), item), drop=FALSE], na.rm = TRUE)
  suppressWarnings(cor(dat_tr_items[[item]], total_wo_item, use="pairwise.complete.obs"))
})
citc_df <- data.frame(item = names(citc), CITC = as.numeric(citc))

alpha_all <- psych::alpha(dat_tr_items, check.keys = FALSE)
alpha_drop_df <- as.data.frame(alpha_all$alpha.drop)
alpha_if_deleted <- data.frame(item = rownames(alpha_drop_df),
                               alpha_if_deleted = alpha_drop_df$`std.alpha`, row.names = NULL)
alpha_total <- alpha_all$total$`std.alpha`

item_analysis_summary <- item_ttest %>%
  left_join(citc_df, by="item") %>%
  left_join(alpha_if_deleted, by="item") %>%
  mutate(
    flag_t_nonsig   = p_value >= TTEST_P_MAX,
    flag_citc_fail  = CITC < CITC_MIN,
    flag_alpha_gain = (alpha_if_deleted - alpha_total) >= ALPHA_GAIN_MIN,
    remove_stage1_suggest = flag_t_nonsig | flag_citc_fail | flag_alpha_gain,
    reason_stage1 = mapply(function(a,b,c){
      paste(c(
        if (a) "t-test non-significant (p>=.05)" else NULL,
        if (b) "CITC<.40" else NULL,
        if (c) "alpha increases by >= .01 after deletion" else NULL
      ), collapse = "; ")
    }, flag_t_nonsig, flag_citc_fail, flag_alpha_gain)
  )

write.csv(item_analysis_summary, "outputs/stage1/ItemAnalysis_full.csv", row.names = FALSE)
remove_stage1        <- item_analysis_summary$item[item_analysis_summary$remove_stage1_suggest]
items_after_stage1   <- setdiff(item_in_train, remove_stage1)
write.csv(data.frame(item = item_in_train,
                     keep = !(item_in_train %in% remove_stage1),
                     reason = ifelse(item_in_train %in% remove_stage1,
                                     item_analysis_summary$reason_stage1[match(item_in_train, item_analysis_summary$item)], "")),
          "outputs/stage1/ItemAnalysis_decisions.csv", row.names = FALSE)

# EFA (train)
efa_data <- dat_tr_items[, items_after_stage1, drop=FALSE]
Rpoly    <- psych::polychoric(efa_data)$rho

sink("outputs/stage1/EFA_KMO_Bartlett.txt")
print(psych::KMO(Rpoly))
print(psych::cortest.bartlett(Rpoly, n = nrow(efa_data)))
sink()

png("outputs/stage1/EFA_parallel.png", width=1400, height=1000, res=200)
pa_out <- psych::fa.parallel(Rpoly, n.obs = nrow(efa_data), fm="pa", fa="fa", error.bars = FALSE, main="Parallel Analysis (polychoric)")
dev.off()
cat("Parallel-suggested factors: ", pa_out$nfact %||% NA, "\n", file="outputs/stage1/EFA_parallel_suggest.txt")

efa_result <- psych::fa(Rpoly, nfactors = N_FACTORS_EFA, fm = "pa", rotate = "oblimin", n.obs = nrow(efa_data))
write.csv(round(unclass(efa_result$loadings), 3), "outputs/stage1/EFA_pattern_matrix.csv")
if (!is.null(efa_result$Phi)) write.csv(round(efa_result$Phi, 3), "outputs/stage1/EFA_phi_cor.csv")
write.csv(round(efa_result$Vaccounted, 3), "outputs/stage1/EFA_variance_accounted.csv")

L1 <- as.matrix(efa_result$loadings); absL1 <- abs(L1)
max1 <- apply(absL1, 1, max)
sec1 <- apply(absL1, 1, function(x) sort(x, decreasing=TRUE)[min(2,length(x))])
gap1 <- max1 - sec1
comm_efa <- if (!is.null(efa_result$communality)) efa_result$communality else rowSums(as.matrix(unclass(efa_result$loadings))^2)
names(comm_efa) <- rownames(L1)

efa_items_df <- data.frame(
  item = rownames(L1),
  primary_loading = max1,
  second_loading  = sec1,
  gap             = gap1,
  h2              = as.numeric(comm_efa),
  drop_low_loading  = max1 < LOAD_MIN,
  drop_crossloading = gap1 < CROSS_GAP_MIN,
  drop_low_h2       = as.numeric(comm_efa) < H2_MIN
)
efa_items_df$drop_any <- with(efa_items_df, drop_low_loading | drop_crossloading | drop_low_h2)
efa_items_df$reason <- apply(efa_items_df, 1, function(r){
  rs <- c()
  if (as.logical(r["drop_low_loading"]))  rs <- c(rs, paste0("primary<", LOAD_MIN))
  if (as.logical(r["drop_crossloading"])) rs <- c(rs, paste0("gap<", CROSS_GAP_MIN))
  if (as.logical(r["drop_low_h2"]))       rs <- c(rs, paste0("h2<", H2_MIN))
  paste(rs, collapse = "; ")
})
write.csv(efa_items_df, "outputs/stage1/EFA_drop_rules.csv", row.names = FALSE)

items_after_efa <- setdiff(items_after_stage1, efa_items_df$item[efa_items_df$drop_any])

efa_data2  <- dat_tr_items[, items_after_efa, drop=FALSE]
Rpoly2     <- psych::polychoric(efa_data2)$rho
efa_result2<- psych::fa(Rpoly2, nfactors = N_FACTORS_EFA, fm="pa", rotate="oblimin", n.obs = nrow(efa_data2))

write.csv(data.frame(item = items_after_efa),                   "outputs/stage1/Items_after_EFA.csv", row.names = FALSE)
write.csv(round(unclass(efa_result2$loadings), 3),              "outputs/stage1/EFA2_pattern_matrix.csv")
if (!is.null(efa_result2$Phi)) write.csv(round(efa_result2$Phi, 3), "outputs/stage1/EFA2_phi_cor.csv")
write.csv(round(efa_result2$Vaccounted, 3),                     "outputs/stage1/EFA2_variance_accounted.csv")

L2 <- as.matrix(efa_result2$loadings); absL2 <- abs(L2)
prim2 <- apply(absL2, 1, which.max)
max2  <- apply(absL2, 1, max)
sec2  <- apply(absL2, 1, function(x) sort(x, decreasing = TRUE)[min(2, length(x))])
gap2  <- max2 - sec2
h2_2  <- if (!is.null(efa_result2$communality)) efa_result2$communality else rowSums(as.matrix(unclass(efa_result2$loadings))^2)
efa2_items_table <- data.frame(
  Item            = rownames(L2),
  Factor          = paste0("F", prim2),
  PrimaryLoading  = round(max2, 3),
  SecondLoading   = round(sec2, 3),
  Gap             = round(gap2, 3),
  h2              = round(as.numeric(h2_2), 3),
  check.names = FALSE
)
write.csv(efa2_items_table, "outputs/stage1/EFA2_items_table.csv", row.names = FALSE)

saveRDS(list(
  items_initial        = item_names,
  items_after_stage1   = items_after_stage1,
  items_after_efa      = items_after_efa,
  efa_result_final     = efa_result2,
  Rpoly_final          = Rpoly2,
  dat_train            = dat_train,
  dat_test             = dat_test,
  efa2_items_table     = efa2_items_table
), file = "outputs/stage1/stage1_artifacts.rds")

# -------------------- Stage 2: CFA/MIRT/Validity --------
dir.create("outputs/stage2/long", showWarnings = FALSE, recursive = TRUE)

L2 <- as.matrix(efa_result2$loadings); colnames(L2) <- paste0("F", seq_len(ncol(L2)))
assign_df2 <- data.frame(
  item    = rownames(L2),
  factor  = colnames(L2)[ max.col(abs(L2), ties.method = "first") ],
  loading = apply(L2, 1, function(x) x[which.max(abs(x))]),
  stringsAsFactors = FALSE
)
mk_line <- function(f){
  its <- assign_df2$item[assign_df2$factor == f]
  if (length(its)>=2) paste0(f," =~ ", paste(its, collapse=" + ")) else NA_character_
}
model_lines2 <- na.omit(vapply(unique(assign_df2$factor), mk_line, character(1)))
model_cfa <- paste(model_lines2, collapse = "\n")

items_for_cfa <- rownames(L2)
cfa_data <- dat_test[, items_for_cfa]
cfa_data[] <- lapply(cfa_data, ordered)
fit_cfa <- lavaan::cfa(model_cfa, data=cfa_data, ordered=colnames(cfa_data), estimator="WLSMV", std.lv=TRUE)
fit_vals <- lavaan::fitMeasures(fit_cfa, c("chisq","df","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr"))
fit_tab_long <- data.frame(Model = sprintf("Long (%d items) CFA (test)", length(items_for_cfa)), t(as.data.frame(fit_vals)), check.names = FALSE)
write.csv(fit_tab_long, "outputs/stage2/long/fit_table_long.csv", row.names = FALSE)

omega_long <- try(semTools::reliability(fit_cfa), silent = TRUE)
if (!inherits(omega_long, "try-error")) write.csv(omega_long, "outputs/stage2/long/omega_long_validation.csv")

if ("Gender" %in% colnames(dat)) {
  mg_data <- rbind(
    cbind(dat_train[, items_for_cfa, drop = FALSE], Gender = dat_train$Gender),
    cbind(dat_test[,  items_for_cfa, drop = FALSE], Gender = dat_test$Gender)
  )
  mg_data[items_for_cfa] <- lapply(mg_data[items_for_cfa], ordered)
  mg_data$Gender <- factor(mg_data$Gender)
  if (nlevels(mg_data$Gender) >= 2) {
    get_measures <- function(fit) {
      fm <- lavaan::fitMeasures(fit, c("cfi","tli","rmsea","srmr","chisq","df","pvalue"))
      data.frame(t(unclass(fm)), check.names = FALSE)
    }
    fit_config <- lavaan::cfa(model_cfa, data=mg_data, group="Gender", ordered=items_for_cfa, estimator="WLSMV", parameterization="delta", std.lv=TRUE)
    fit_metric <- lavaan::cfa(model_cfa, data=mg_data, group="Gender", group.equal=c("loadings"), ordered=items_for_cfa, estimator="WLSMV", parameterization="delta", std.lv=TRUE)
    fit_scalar <- lavaan::cfa(model_cfa, data=mg_data, group="Gender", group.equal=c("loadings","thresholds"), ordered=items_for_cfa, estimator="WLSMV", parameterization="delta", std.lv=TRUE)
    rows <- list(
      cbind(Level="Configural",                   get_measures(fit_config)),
      cbind(Level="Metric (loadings)",            get_measures(fit_metric)),
      cbind(Level="Scalar (loadings+thresholds)", get_measures(fit_scalar))
    )
    tab <- do.call(rbind, rows); row.names(tab) <- NULL
    tab$DeltaCFI   <- c(NA, diff(tab$cfi)); tab$DeltaRMSEA <- c(NA, diff(tab$rmsea))
    write.csv(tab, "outputs/stage2/long/mi_gender_summary.csv", row.names = FALSE)
  }
}

mirt_base <- rbind(dat_train[, items_for_cfa], dat_test[, items_for_cfa])
mirt_base[] <- lapply(mirt_base, as.integer)
items_by_factor <- split(assign_df2$item, assign_df2$factor)
items_by_factor <- lapply(items_by_factor, function(v) intersect(v, items_for_cfa))
items_by_factor <- items_by_factor[sapply(items_by_factor, function(v) length(v) >= 2)]
grid <- seq(-3, 3, length.out = 121)
testinfo_by_dim <- list(); mods_uni <- list()
for (f in names(items_by_factor)) {
  subdat <- mirt_base[, items_by_factor[[f]], drop=FALSE]
  mod_uni <- mirt::mirt(subdat, 1, itemtype="graded", method="QMCEM", verbose=FALSE)
  mods_uni[[f]] <- mod_uni
  info <- mirt::testinfo(mod_uni, Theta = matrix(grid, ncol = 1))
  testinfo_by_dim[[f]] <- as.numeric(info)
  dfp <- data.frame(theta = grid, info = as.numeric(info))
  p <- ggplot2::ggplot(dfp, ggplot2::aes(theta, info)) +
    ggplot2::geom_line() + ggplot2::geom_hline(yintercept=0, linetype=2) +
    ggplot2::labs(title = paste0("Test Information - ", f), x = expression(theta), y = "Information") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path("outputs/stage2/long", paste0("test_info_", f, ".png")), p, width = 6, height = 4, dpi = 600)
}
marginal_rel <- data.frame(Factor = names(items_by_factor), Reliability = NA_real_)
for (i in seq_along(items_by_factor)) {
  f <- names(items_by_factor)[i]
  fs <- mirt::fscores(mods_uni[[f]], method="EAP", full.scores=TRUE, full.scores.SE=TRUE, QMC=TRUE)
  theta_hat <- as.numeric(fs); se <- attr(fs,"SE"); if (is.null(se)) se <- attr(fs,"SE2")
  if (!is.null(se)) {
    marginal_rel$Reliability[i] <- 1 - mean(se^2, na.rm=TRUE)/stats::var(theta_hat, na.rm=TRUE)
  } else {
    info <- pmax(testinfo_by_dim[[f]], 1e-6)
    w <- dnorm(grid); w <- w/sum(w)
    marginal_rel$Reliability[i] <- 1 - sum((1/info)*w)
  }
}
write.csv(marginal_rel, "outputs/stage2/long/mirt_marginal_reliability_long.csv", row.names = FALSE)

cr_ave_long <- cr_ave_table(fit_cfa)
write.csv(cbind(Factor=rownames(cr_ave_long), cr_ave_long), "outputs/stage2/long/cr_ave_long.csv", row.names=FALSE)
Xlong_for_htmt <- dat_test[, items_for_cfa, drop = FALSE]
htmt_long <- htmt_safely(fit_cfa, model_cfa, Xlong_for_htmt)
if (!is.null(htmt_long)) write.csv(htmt_long, "outputs/stage2/long/htmt_long.csv", row.names = FALSE)

saveRDS(list(
  items_final_cfa = items_for_cfa,
  model_cfa       = model_cfa,
  assign_df_pruned= assign_df2,
  dat_train       = dat_train,
  dat_test        = dat_test
), file = "outputs/stage2/long/stage2_artifacts.rds")

# Validity (long form)
LET_items  <- paste0("LET",  1:6)
SWLS_items <- paste0("SWLS", 1:5)
WHO5_items <- paste0("WHO",  1:5)
K10_items  <- paste0("K",    1:10)
CPSS_items <- paste0("CPSS", 1:4)
LET_rev    <- c("LET1","LET3","LET5")

X_test <- dat_test
need_cols <- unique(c(items_for_cfa, LET_items, SWLS_items, WHO5_items, K10_items, CPSS_items))
need_cols <- intersect(need_cols, colnames(X_test))
X_test[need_cols] <- lapply(X_test[need_cols], to_numeric_safe)
if (length(intersect(LET_rev, colnames(X_test)))) {
  X_test[intersect(LET_rev, colnames(X_test))] <- lapply(X_test[intersect(LET_rev, colnames(X_test))], reverse_by_range)
}

ELS_mean_long <- row_mean_safe(X_test, items_for_cfa)
items_by_factor_val <- split(assign_df2$item, assign_df2$factor)
items_by_factor_val <- lapply(items_by_factor_val, function(v) intersect(v, items_for_cfa))
items_by_factor_val <- items_by_factor_val[sapply(items_by_factor_val, function(v) length(v) >= 2)]
factor_means <- lapply(items_by_factor_val, function(cols) row_mean_safe(X_test, cols))
names(factor_means) <- paste0(names(factor_means), "_mean")

LET_mean  <- row_mean_safe(X_test, intersect(LET_items,  colnames(X_test)))
SWLS_mean <- row_mean_safe(X_test, intersect(SWLS_items, colnames(X_test)))
WHO5_mean <- row_mean_safe(X_test, intersect(WHO5_items, colnames(X_test)))
K10_sum   <- row_sum_safe (X_test, intersect(K10_items,  colnames(X_test)))
cpss_total <- function(df, items){
  cols <- intersect(items, colnames(df)); if (!length(cols)) return(rep(NA_real_, nrow(df)))
  Z <- df[, cols, drop = FALSE]; Z[] <- lapply(Z, to_numeric_safe)
  rev_cols <- intersect(c("CPSS2","CPSS3"), colnames(Z))
  if (length(rev_cols)) Z[rev_cols] <- lapply(Z[rev_cols], reverse_by_range)
  rowSums(Z, na.rm = TRUE) - 4
}
CPSS_tot  <- cpss_total(X_test, CPSS_items)

valid_df_long <- data.frame(ELS_mean = ELS_mean_long, factor_means, LET_mean, SWLS_mean, WHO5_mean, K10_sum, CPSS_tot)

conv_pairs <- list(c("ELS_mean","LET_mean"), c("ELS_mean","SWLS_mean"), c("ELS_mean","WHO5_mean"))
for (f in names(factor_means)) conv_pairs <- c(conv_pairs, list(c(f,"LET_mean")), list(c(f,"SWLS_mean")), list(c(f,"WHO5_mean")))
conv_long <- dplyr::bind_rows(lapply(conv_pairs, function(p) pair_row(valid_df_long, p[1], p[2], "pearson")))
conv_long$p_fdr <- p.adjust(conv_long$p, method = "BH")
conv_print <- within(conv_long, { r <- round(r,3); ci_low <- round(ci_low,3); ci_high <- round(ci_high,3); p_txt <- fmt_p(p); p_fdr_txt <- fmt_p(p_fdr) })
write.csv(conv_print[, c("var1","var2","r","ci_low","ci_high","p_txt","p_fdr_txt")], "outputs/stage2/long/convergent_validity_long.csv", row.names = FALSE)

disc_pairs <- list(c("ELS_mean","K10_sum"), c("ELS_mean","CPSS_tot"))
for (f in names(factor_means)) disc_pairs <- c(disc_pairs, list(c(f,"K10_sum")), list(c(f,"CPSS_tot")))
disc_long <- dplyr::bind_rows(lapply(disc_pairs, function(p) pair_row(valid_df_long, p[1], p[2], "pearson")))
disc_long$p_fdr <- p.adjust(disc_long$p, method = "BH")
disc_print <- within(disc_long, { r <- round(r,3); ci_low <- round(ci_low,3); ci_high <- round(ci_high,3); p_txt <- fmt_p(p); p_fdr_txt <- fmt_p(p_fdr) })
write.csv(disc_print[, c("var1","var2","r","ci_low","ci_high","p_txt","p_fdr_txt")], "outputs/stage2/long/discriminant_validity_long.csv", row.names = FALSE)

# -------------------- Stage 3: ACO short form -----------
dir.create("outputs/stage3/short", showWarnings = FALSE, recursive = TRUE)
art <- readRDS("outputs/stage2/long/stage2_artifacts.rds")
items_final_cfa  <- art$items_final_cfa
model_cfa_long   <- art$model_cfa
assign_df_pruned <- art$assign_df_pruned
dat_train        <- art$dat_train
dat_test         <- art$dat_test

auto <- function(pat, df) grep(pat, colnames(df), value = TRUE)
LET_items_auto  <- auto("^LET\\d+$",  dat_test); SWLS_items_auto <- auto("^SWLS\\d+$", dat_test)
WHO5_items_auto <- auto("^WHO\\d+$",  dat_test); K10_items_auto  <- auto("^K\\d+$",    dat_test)
CPSS_items_auto <- auto("^CPSS\\d+$", dat_test)
LET_items  <- if (length(LET_items_auto))  LET_items_auto  else LET_items
SWLS_items <- if (length(SWLS_items_auto)) SWLS_items_auto else SWLS_items
WHO5_items <- if (length(WHO5_items_auto)) WHO5_items_auto else WHO5_items
K10_items  <- if (length(K10_items_auto))  K10_items_auto  else K10_items
CPSS_items <- if (length(CPSS_items_auto)) CPSS_items_auto else CPSS_items
LET_rev    <- intersect(c("LET1","LET3","LET5"), colnames(dat_train))

if (!requireNamespace("ShortForm", quietly = TRUE)) {
  stop("Package 'ShortForm' is required for Stage 3. Install with install.packages('ShortForm').")
}
library(ShortForm)

RNGkind("L'Ecuyer-CMRG"); set.seed(SEED_ACO)
ANTS  <- 40; STEPS <- 60; EVAP <- 0.80
CANDIDATE_ITEMS <- items_final_cfa
map_df <- subset(assign_df_pruned, item %in% CANDIDATE_ITEMS)
factors <- unique(map_df$factor); n_by_factor <- table(map_df$factor)
if (!all(n_by_factor >= 2L)) stop("Each factor must have at least 2 candidate items.")
i_per_f <- rep(2L, length(n_by_factor)); names(i_per_f) <- names(n_by_factor)
list_items <- split(map_df$item, map_df$factor); list_items <- list_items[factors]

ordify_items <- function(df, items){
  out <- df
  out[items] <- lapply(out[items], function(x){
    if (is.ordered(x)) return(x)
    if (is.factor(x))  return(ordered(x))
    if (is.numeric(x)) return(ordered(as.integer(x)))
    ordered(x)
  })
  out
}
dat_train_ord <- ordify_items(dat_train, CANDIDATE_ITEMS)

score_scale <- function(df, items, reverse = NULL) {
  z <- df[, items, drop = FALSE]
  z[] <- lapply(z, function(x) if (is.ordered(x) || is.factor(x)) as.numeric(x) else to_numeric_safe(x))
  if (length(reverse)) z[reverse] <- lapply(z[reverse], function(x) max(x,na.rm=T)+min(x,na.rm=T)-x)
  rowMeans(z, na.rm = TRUE)
}
dat_train$LET_tot  <- if (length(LET_items))  score_scale(dat_train, LET_items, LET_rev) else NA
dat_train$SWLS_tot <- if (length(SWLS_items)) score_scale(dat_train, SWLS_items) else NA
dat_train$WHO5_tot <- if (length(WHO5_items)) score_scale(dat_train, WHO5_items) else NA
dat_train$K10_tot  <- if (length(K10_items))  score_scale(dat_train, K10_items)  else NA
cpss_total <- function(df, items){
  cols <- intersect(items, colnames(df)); if (!length(cols)) return(rep(NA_real_, nrow(df)))
  Z <- df[, cols, drop = FALSE]; Z[] <- lapply(Z, to_numeric_safe)
  rev_cols <- intersect(c("CPSS2","CPSS3"), colnames(Z))
  if (length(rev_cols)) Z[rev_cols] <- lapply(Z[rev_cols], reverse_by_range)
  rowSums(Z, na.rm = TRUE) - 4
}
dat_train$CPSS_tot <- if (length(CPSS_items)) cpss_total(dat_train, CPSS_items) else NA

build_factor_line <- function(f, its) paste0(f, " =~ ", paste(its, collapse = " + "))
factor_lines <- vapply(factors, function(f) build_factor_line(f, list_items[[f]]), character(1))
reg_lines <- c(
  paste0("LET_tot  ~ ", paste(factors, collapse=" + ")),
  paste0("SWLS_tot ~ ", paste(factors, collapse=" + ")),
  paste0("WHO5_tot ~ ", paste(factors, collapse=" + ")),
  paste0("K10_tot  ~ ", paste(factors, collapse=" + ")),
  paste0("CPSS_tot ~ ", paste(factors, collapse=" + "))
)
antModel <- paste(c(factor_lines, reg_lines), collapse = "\n")
dat_train_ord_ext <- cbind(dat_train_ord, dat_train[, c("LET_tot","SWLS_tot","WHO5_tot","K10_tot","CPSS_tot"), drop=FALSE])

SEEDS <- SEED_ACO + 0:9
dir.create("outputs/stage3/short/multiseed", showWarnings = FALSE, recursive = TRUE)

get_best_syntax_from_S4 <- function(obj) {
  syn <- try(methods::slot(obj, "best_syntax"), silent = TRUE)
  if (!inherits(syn,"try-error") && is.character(syn) && any(grepl("=~", syn))) return(paste(syn, collapse="\n"))
  bm <- try(methods::slot(obj, "best_model"), silent = TRUE)
  if (!inherits(bm,"try-error") && inherits(bm, "lavaan")) {
    PT <- lavaan::parTable(bm); PT <- PT[PT$op == "=~" & is.finite(PT$est), ]
    if (nrow(PT)) {
      syn <- PT |>
        split(.$lhs) |>
        lapply(function(df) paste0(unique(df$lhs), " =~ ", paste(df$rhs, collapse=" + "))) |>
        unlist() |> paste(collapse="\n")
      if (grepl("=~", syn)) return(syn)
    }
  }
  stop("Unable to extract best_syntax from ShortForm object.")
}
parse_items_by_factor <- function(syn){
  lines <- unlist(strsplit(syn, "\n"))
  flines <- lines[grepl("=~", lines, fixed=TRUE)]
  out <- lapply(flines, function(ln){
    sp <- strsplit(ln, "=~", fixed=TRUE)[[1]]
    f  <- stringr::str_trim(sp[1])
    its<- stringr::str_trim(unlist(strsplit(sp[2], "\\+")))
    list(factor=f, items=its)
  })
  setNames(lapply(out, `[[`, "items"), vapply(out, `[[`, "", "factor"))
}

build_aco <- function(seed){
  set.seed(seed)
  ShortForm::antcolony.lavaan(
    data        = dat_train_ord_ext,
    ants        = ANTS,
    evaporation = EVAP,
    antModel    = antModel,
    list.items  = unname(list_items),
    full        = length(CANDIDATE_ITEMS),
    i.per.f     = as.integer(i_per_f),
    factors     = factors,
    steps       = STEPS,
    lavaan.model.specs = list(
      model.type = "sem", estimator = "WLSMV",
      ordered    = CANDIDATE_ITEMS, std.lv = TRUE,
      parameterization = "delta", auto.th = TRUE, auto.delta = TRUE
    ),
    pheromone.calculation = "variance",
    fit.indices           = c("cfi","tli","rmsea"),
    fit.statistics.test   = "(cfi > 0.95)&(tli > 0.95)&(rmsea < 0.06)",
    summaryfile  = file.path("outputs/stage3/short/multiseed", paste0("ACO_summary_seed", seed, ".txt")),
    feedbackfile = file.path("outputs/stage3/short/multiseed", paste0("ACO_iteration_seed", seed, ".html")),
    max.run      = ANTS*STEPS,
    parallel     = FALSE
  )
}
multi_runs <- lapply(SEEDS, function(s){
  obj <- build_aco(s)
  syn <- get_best_syntax_from_S4(obj)
  cho <- parse_items_by_factor(syn)
  writeLines(syn, file.path("outputs/stage3/short/multiseed", paste0("short_syntax_seed", s, ".lav")))
  list(seed=s, syntax=syn, choice=cho)
})

vote_rows <- list()
for (rr in multi_runs) {
  for (f in names(rr$choice)) {
    its <- rr$choice[[f]]
    vote_rows[[length(vote_rows)+1]] <- data.frame(factor=f, item=its, seed=rr$seed)
  }
}
vote_tbl <- do.call(rbind, vote_rows)
vote_tbl <- as.data.frame(table(vote_tbl$factor, vote_tbl$item), stringsAsFactors = FALSE)
colnames(vote_tbl) <- c("factor","item","freq")

aux2 <- assign_df_pruned[, c("item","factor","loading"), drop=FALSE]
aux2$loading <- suppressWarnings(as.numeric(aux2$loading))
vote_tbl <- merge(vote_tbl, aux2, by = c("item","factor"), all.x = TRUE)
vote_tbl$loading[!is.finite(vote_tbl$loading)] <- -Inf
write.csv(vote_tbl, "outputs/stage3/short/aco_vote_table.csv", row.names = FALSE)

pick_two_per_factor <- function(df){ df <- df[order(-df$freq, -df$loading, df$item), ]; head(df$item, 2L) }
short_choice_voted <- tapply(seq_len(nrow(vote_tbl)), vote_tbl$factor, function(idx){
  df <- vote_tbl[idx, c("factor","item","freq","loading")]
  data.frame(factor = unique(df$factor), items = I(list(pick_two_per_factor(df))))
})
short_choice_voted <- setNames(lapply(short_choice_voted, function(x) x$items[[1]]),
                               nm = vapply(short_choice_voted, function(x) as.character(x$factor[1]), ""))

model_short <- paste(vapply(names(short_choice_voted),
                            function(f) paste0(f, " =~ ", paste(short_choice_voted[[f]], collapse = " + ")),
                            character(1)),
                     collapse = "\n")
writeLines(model_short, "outputs/stage3/short/short_cfa_model.lav")

short_items <- unique(unlist(short_choice_voted, use.names = FALSE))
fit_tr <- lavaan::cfa(model_short, data = ordify_like_short(dat_train[, short_items, drop=FALSE]),
                      estimator="WLSMV", ordered=short_items, std.lv=TRUE, parameterization="delta")
fit_te <- lavaan::cfa(model_short, data = ordify_like_short(dat_test[,  short_items, drop=FALSE]),
                      estimator="WLSMV", ordered=short_items, std.lv=TRUE, parameterization="delta")

fit_vals_short_te <- lavaan::fitMeasures(fit_te, c("chisq","df","pvalue","cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr"))
fit_tab_short <- data.frame(Model = sprintf("Short (%d items) CFA (test)", length(short_items)),
                            t(as.data.frame(fit_vals_short_te)), check.names = FALSE)
write.csv(fit_tab_short, "outputs/stage3/short/fit_table_short.csv", row.names = FALSE)

cr_ave_short <- cr_ave_table(fit_te)
write.csv(cbind(Factor=rownames(cr_ave_short), cr_ave_short), "outputs/stage3/short/cr_ave_short.csv", row.names=FALSE)
htmt_short <- htmt_safely(fit_te, model_short, dat_test[, short_items, drop = FALSE])
if (!is.null(htmt_short)) write.csv(htmt_short, "outputs/stage3/short/htmt_short.csv", row.names = FALSE)

assign_short <- do.call(rbind, lapply(names(short_choice_voted), function(f){
  data.frame(item = short_choice_voted[[f]], factor = f, stringsAsFactors = FALSE)
}))
grid <- seq(-3, 3, length.out = 121); testinfo_by_dim <- list(); mods_uni <- list()
mirt_base_short <- rbind(dat_train[, short_items], dat_test[, short_items]); mirt_base_short[] <- lapply(mirt_base_short, as.integer)
for (f in unique(assign_short$factor)) {
  its <- assign_short$item[assign_short$factor == f]
  mod_uni <- mirt::mirt(mirt_base_short[, its, drop=FALSE], 1, itemtype="graded", method="QMCEM", verbose=FALSE); mods_uni[[f]] <- mod_uni
  info <- mirt::testinfo(mod_uni, Theta = matrix(grid, ncol=1)); testinfo_by_dim[[f]] <- as.numeric(info)
  dfp <- data.frame(theta=grid, info=as.numeric(info))
  p <- ggplot2::ggplot(dfp, ggplot2::aes(theta, info)) + ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept=0, linetype=2) + ggplot2::labs(title=paste0("Test Information - ", f, " (Short)"), x=expression(theta), y="Information") +
    ggplot2::theme_minimal(base_size=12)
  ggplot2::ggsave(file.path("outputs/stage3/short", paste0("test_info_short_", f, ".png")), p, width=6, height=4, dpi=600)
}
marginal_rel_short <- data.frame(Factor = unique(assign_short$factor), Reliability = NA_real_)
for (i in seq_len(nrow(marginal_rel_short))) {
  f <- marginal_rel_short$Factor[i]
  fs <- mirt::fscores(mods_uni[[f]], method="EAP", full.scores=TRUE, full.scores.SE=TRUE, QMC=TRUE)
  theta_hat <- as.numeric(fs); se <- attr(fs,"SE"); if (is.null(se)) se <- attr(fs,"SE2")
  if (!is.null(se)) {
    marginal_rel_short$Reliability[i] <- 1 - mean(se^2, na.rm=TRUE)/stats::var(theta_hat, na.rm=TRUE)
  } else {
    info <- pmax(testinfo_by_dim[[f]], 1e-6); w <- dnorm(grid); w <- w/sum(w)
    marginal_rel_short$Reliability[i] <- 1 - sum((1/info) * w)
  }
}
write.csv(marginal_rel_short, "outputs/stage3/short/mirt_marginal_reliability_short.csv", row.names = FALSE)

# Validity (short)
X_test_short <- dat_test
need_cols2 <- unique(c(short_items, LET_items, SWLS_items, WHO5_items, K10_items, CPSS_items))
need_cols2 <- intersect(need_cols2, colnames(X_test_short))
X_test_short[need_cols2] <- lapply(X_test_short[need_cols2], to_numeric_safe)
if (length(intersect(LET_rev, colnames(X_test_short)))) {
  X_test_short[intersect(LET_rev, colnames(X_test_short))] <- lapply(X_test_short[intersect(LET_rev, colnames(X_test_short))], reverse_by_range)
}
ELS8_mean <- row_mean_safe(X_test_short, short_items)
items_by_factor_short <- split(assign_short$item, assign_short$factor)
items_by_factor_short <- items_by_factor_short[sapply(items_by_factor_short, length) >= 2]
factor_means_short <- lapply(items_by_factor_short, function(cols) row_mean_safe(X_test_short, cols))
names(factor_means_short) <- paste0(names(factor_means_short), "_mean")
LET_mean  <- row_mean_safe(X_test_short, intersect(LET_items,  colnames(X_test_short)))
SWLS_mean <- row_mean_safe(X_test_short, intersect(SWLS_items, colnames(X_test_short)))
WHO5_mean <- row_mean_safe(X_test_short, intersect(WHO5_items, colnames(X_test_short)))
K10_sum   <- row_sum_safe (X_test_short, intersect(K10_items,  colnames(X_test_short)))
CPSS_tot  <- cpss_total   (X_test_short, CPSS_items)

valid_df_short <- data.frame(ELS8_mean = ELS8_mean, factor_means_short, LET_mean, SWLS_mean, WHO5_mean, K10_sum, CPSS_tot)
conv_pairs <- list(c("ELS8_mean","LET_mean"), c("ELS8_mean","SWLS_mean"), c("ELS8_mean","WHO5_mean"))
for (f in names(factor_means_short)) conv_pairs <- c(conv_pairs, list(c(f,"LET_mean")), list(c(f,"SWLS_mean")), list(c(f,"WHO5_mean")))
conv_short <- dplyr::bind_rows(lapply(conv_pairs, function(p) pair_row(valid_df_short, p[1], p[2], "pearson"))); conv_short$p_fdr <- p.adjust(conv_short$p, method = "BH")
conv_short_print <- within(conv_short, { r<-round(r,3); ci_low<-round(ci_low,3); ci_high<-round(ci_high,3); p_txt<-fmt_p(p); p_fdr_txt<-fmt_p(p_fdr) })
write.csv(conv_short_print[, c("var1","var2","r","ci_low","ci_high","p_txt","p_fdr_txt")], "outputs/stage3/short/convergent_validity_short.csv", row.names = FALSE)

disc_pairs <- list(c("ELS8_mean","K10_sum"), c("ELS8_mean","CPSS_tot"))
for (f in names(factor_means_short)) disc_pairs <- c(disc_pairs, list(c(f,"K10_sum")), list(c(f,"CPSS_tot")))
disc_short <- dplyr::bind_rows(lapply(disc_pairs, function(p) pair_row(valid_df_short, p[1], p[2], "pearson"))); disc_short$p_fdr <- p.adjust(disc_short$p, method = "BH")
disc_short_print <- within(disc_short, { r<-round(r,3); ci_low<-round(ci_low,3); ci_high<-round(ci_high,3); p_txt<-fmt_p(p); p_fdr_txt<-fmt_p(p_fdr) })
write.csv(disc_short_print[, c("var1","var2","r","ci_low","ci_high","p_txt","p_fdr_txt")], "outputs/stage3/short/discriminant_validity_short.csv", row.names = FALSE)

saveRDS(list(
  short_choice = short_choice_voted,
  short_items  = short_items,
  model_short  = model_short,
  fit_tr       = fit_tr
), file = "outputs/stage3/short/stage3_short_artifacts.rds")

message("Pipeline completed. See outputs/stage1, stage2, stage3.")
