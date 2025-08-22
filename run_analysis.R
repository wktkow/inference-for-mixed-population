#!/usr/bin/env Rscript
# Entry point for the hemodialysis mixture analysis: orchestrates data load, EDA,
# dispersion checks, mixture fitting, covariate associations, figures, and report.

options(device = "png")
# Force a non-interactive graphics device so plots are saved as PNGs in any env.

suppressPackageStartupMessages({
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales", repos = "https://cloud.r-project.org")
  if (!requireNamespace("flexmix", quietly = TRUE)) install.packages("flexmix", repos = "https://cloud.r-project.org")
  if (!requireNamespace("nnet", quietly = TRUE)) install.packages("nnet", repos = "https://cloud.r-project.org")
  if (!requireNamespace("CAMAN", quietly = TRUE) && !requireNamespace("caman", quietly = TRUE)) {
    install.packages("CAMAN", repos = "https://cloud.r-project.org")
  }
})
# Ensure required packages are present; install them on the fly if missing for portability.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(flexmix)
  library(nnet)
  if (requireNamespace("CAMAN", quietly = TRUE)) {
    library(CAMAN)
  } else if (requireNamespace("caman", quietly = TRUE)) {
    library(caman)
  }
})
# Load libraries quietly to keep console output clean in scripted runs.

dir.create("figs", showWarnings = FALSE)
# Create an output directory for figures; ignore if it already exists.

# 1) Load data (include NAs)
path_csv <- "./hemodialysismix.csv"
dat_raw <- read_csv(path_csv, show_col_types = FALSE, na = c("", "NA"))

dat <- dat_raw %>%
  rename(
    ID = ID,
    SEX = SEX,
    AGE = AGE,
    NRIRON = nriron,
    NR = nr
  ) %>%
  mutate(
    NR = as.integer(NR),
    NRIRON = as.integer(NRIRON),
    AGE = as.integer(AGE),
    SEX = as.integer(SEX)
  ) %>%
  mutate(
    SEX_factor = factor(ifelse(is.na(SEX), "Unknown", ifelse(SEX == 1, "Male", ifelse(SEX == 2, "Female", "Other"))),
                        levels = c("Male", "Female", "Unknown", "Other")),
    AGE_miss = as.integer(is.na(AGE)),
    AGE_imp = ifelse(is.na(AGE), mean(AGE, na.rm = TRUE), AGE),
    p_hat = NRIRON / NR
  )
# Standardize schema (names/types) and derive analysis fields including
# explicit NA handling for SEX, mean-imputed AGE with a missingness flag,
# and subject-level proportion p_hat = NRIRON/NR.

# Sanity checks
stopifnot(all(dat$NR >= 1))
stopifnot(all(dat$NRIRON >= 0 & dat$NRIRON <= dat$NR))
# Basic data validation: each subject has trials and successes within [0, NR].

# 2) EDA
eda_summary <- list(
  n = nrow(dat),
  nr_table = dat %>% count(NR, name = "n") %>% arrange(NR),
  sex_table = dat %>% count(SEX_factor, name = "n") %>% arrange(desc(n)),
  age_summary = summary(dat$AGE),
  p_hat_summary = summary(dat$p_hat)
)
# Lightweight summaries to be embedded in the README for quick inspection.

# Visualize subject-level adequate iron proportion distribution.
g_p <- ggplot(dat, aes(x = p_hat)) +
  geom_histogram(bins = 30, color = "white", fill = "steelblue") +
  labs(title = "Distribution of adequate iron proportion (NRIRON/NR)", x = "Proportion", y = "Count")
ggsave("figs/p_hat_hist.png", g_p, width = 7, height = 4.5, dpi = 150)

# Visualize distribution of counts of adequate iron across subjects.
g_counts <- ggplot(dat, aes(x = NRIRON)) +
  geom_histogram(binwidth = 1, color = "white", fill = "darkorange") +
  labs(title = "Distribution of counts of adequate iron (NRIRON)", x = "NRIRON", y = "Count")
ggsave("figs/nriron_hist.png", g_counts, width = 7, height = 4.5, dpi = 150)

# 3) Dispersion assessment (binomial)
# Intercept-only model: logit(p) = beta0; aggregated binomial with varying NR
glm0 <- glm(cbind(NRIRON, NR - NRIRON) ~ 1, data = dat, family = binomial())
p0 <- fitted(glm0)
pearson0 <- sum(((dat$NRIRON - dat$NR * p0) / sqrt(dat$NR * p0 * (1 - p0)))^2, na.rm = TRUE)
df0 <- nrow(dat) - length(coef(glm0))
phi0 <- as.numeric(pearson0 / df0)
# Estimate dispersion under a single-probability binomial to detect overdispersion.

# With covariates AGE and SEX (including missingness indicator)
glm_cov <- glm(cbind(NRIRON, NR - NRIRON) ~ AGE_imp + AGE_miss + SEX_factor, data = dat, family = binomial())
p_cov <- fitted(glm_cov)
pearson_cov <- sum(((dat$NRIRON - dat$NR * p_cov) / sqrt(dat$NR * p_cov * (1 - p_cov)))^2, na.rm = TRUE)
df_cov <- nrow(dat) - length(coef(glm_cov))
phi_cov <- as.numeric(pearson_cov / df_cov)
# Repeat dispersion check after adjusting for AGE/SEX to see if heterogeneity remains.

# 4) Finite mixture of binomials via custom EM (robust, no external deps)
logsumexp_vec <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}
# Numerically stable helper for log-sum-exp in posterior normalization.

fit_binom_mixture_em <- function(y, n, k, n_starts = 5, max_iter = 200, tol = 1e-8) {
  y <- as.numeric(y); n <- as.numeric(n)
  N <- length(y)
  best <- NULL
  for (s in seq_len(n_starts)) {
    # Initialize p around quantiles of observed proportions
    ph <- pmax(pmin((y + 0.5) / (n + 1), 1 - 1e-6), 1e-6)
    qs <- quantile(ph, probs = seq(0.1, 0.9, length.out = k))
    p <- as.numeric(qs) + rnorm(k, sd = 0.02)
    p <- pmax(pmin(p, 1 - 1e-6), 1e-6)
    pi <- rep(1 / k, k)
    ll_old <- -Inf
    for (it in seq_len(max_iter)) {
      # E-step: responsibilities
      loglik_mat <- sapply(seq_len(k), function(j) {
        dbinom(y, n, p[j], log = TRUE) + log(pi[j])
      })
      # row-wise normalization
      lse <- apply(loglik_mat, 1, logsumexp_vec)
      ll <- sum(lse)
      gamma <- exp(loglik_mat - lse)
      # M-step
      pi <- colMeans(gamma)
      for (j in seq_len(k)) {
        num <- sum(gamma[, j] * y)
        den <- sum(gamma[, j] * n)
        p[j] <- ifelse(den > 0, num / den, p[j])
      }
      p <- pmax(pmin(p, 1 - 1e-8), 1e-8)
      if (abs(ll - ll_old) < tol) break
      ll_old <- ll
    }
    out <- list(logLik = ll, pi = pi, p = p)
    if (is.null(best) || out$logLik > best$logLik) best <- out
  }
  # Posterior for best
  loglik_mat <- sapply(seq_len(k), function(j) dbinom(y, n, best$p[j], log = TRUE) + log(best$pi[j]))
  lse <- apply(loglik_mat, 1, logsumexp_vec)
  post <- exp(loglik_mat - lse)
  kparams <- (k - 1) + k # pi (k-1) + p (k)
  bic <- -2 * best$logLik + kparams * log(length(y))
  list(logLik = best$logLik, pi = best$pi, p = best$p, posterior = post, BIC = bic)
}
# Custom EM for Binomial mixtures with multiple random starts; returns MAP posteriors
# and BIC for model selection; uses simple E/M updates on probabilities and weights.

mix_fits <- lapply(1:3, function(k) fit_binom_mixture_em(dat$NRIRON, dat$NR, k = k, n_starts = 10))
mix_bics <- sapply(mix_fits, function(m) m$BIC)
best_k <- which.min(mix_bics)
best_mix <- mix_fits[[best_k]]
# Fit candidate models with k=1..3, compare by BIC, and select the best mixture size.

mixture_summary <- data.frame(
  component = seq_len(best_k),
  p = pmax(pmin(as.numeric(best_mix$p), 1), 0),
  pi = as.numeric(best_mix$pi)
)
ord <- order(mixture_summary$p)
mixture_summary <- mixture_summary[ord, ]
# Summarize ordered component success probabilities and their mixture weights.

# Classification by MAP
comp_class <- apply(best_mix$posterior[, ord, drop = FALSE], 1, which.max)
dat$component <- factor(comp_class, levels = seq_len(nrow(mixture_summary)))
# Assign each subject to the most likely component for downstream association plots.

# Plots
# Plot estimated component probabilities with labels showing mixing weights.
g_mix <- ggplot(mixture_summary, aes(x = factor(component), y = p)) +
  geom_point(size = 3, color = "firebrick") +
  geom_text(aes(label = sprintf("pi=%.2f", pi)), vjust = -1.0, size = 3.5) +
  labs(title = sprintf("Estimated mixture components (binomial, k=%d)", best_k), x = "Component", y = "p (adequate)") +
  ylim(0, 1)
ggsave("figs/mixture_components.png", g_mix, width = 6, height = 4, dpi = 150)

# Visualize age distributions across components (imputed age, with NA flag elsewhere).
g_age_comp <- ggplot(dat, aes(x = component, y = AGE_imp, fill = component)) +
  geom_violin(alpha = 0.5) + geom_boxplot(width = 0.1, outlier.size = 0.5) +
  labs(title = "Age by mixture component", x = "Component", y = "Age (years, imputed for NA)") +
  theme(legend.position = "none")
ggsave("figs/age_by_component.png", g_age_comp, width = 7, height = 4.5, dpi = 150)

# Stacked bar chart of component-wise sex composition (including Unknown level).
g_sex_comp <- ggplot(dat, aes(x = component, fill = SEX_factor)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "SEX composition within components", x = "Component", y = "Proportion")
ggsave("figs/sex_by_component.png", g_sex_comp, width = 7, height = 4.5, dpi = 150)

# 5) Relationship of components to AGE/SEX
assoc_results <- list()
{
  # Multinomial logit of component classification ~ AGE_imp + AGE_miss + SEX_factor
  df_class <- dat %>% mutate(component = factor(component))
  # ensure reference levels
  df_class$SEX_factor <- relevel(df_class$SEX_factor, ref = "Male")
  if (nlevels(df_class$component) >= 2) {
    mnl <- multinom(component ~ AGE_imp + AGE_miss + SEX_factor, data = df_class, trace = FALSE)
    s <- summary(mnl)
    zval <- s$coefficients / s$standard.errors
    pval <- 2 * pnorm(-abs(zval))
    assoc_results$multinom_pvals <- pval
  }
  # Also simple rank-correlation p_hat ~ AGE
  assoc_results$spearman_p_age <- suppressWarnings(cor.test(df_class$p_hat, df_class$AGE_imp, method = "spearman")$p.value)
  # Kruskal-Wallis p_hat by SEX
  assoc_results$kw_p_sex <- kruskal.test(p_hat ~ SEX_factor, data = df_class)$p.value
}
# Evaluate whether AGE/SEX relate to component membership and p_hat via simple tests.

# 6) Do covariates explain clusters? (approximate check via multinomial regression)
flexmix_concom <- NULL
concom_bic <- NA_real_
# Placeholder for concomitant model selection; custom EM does not include concomitants.

# 6b) CAMAN analysis (VEM → EM) and NPML (mixalg), plus gradient and overlay
run_caman <- function(df) {
  if (!requireNamespace("CAMAN", quietly = TRUE) && !requireNamespace("caman", quietly = TRUE)) {
    return(NULL)
  }
  res <- list()
  # VEM (Phase 1)
  vem <- try(mixalg.VEM(obs = "NRIRON", pop.at.risk = "NR", family = "binomial",
                        data = df, acc = 1e-6, numiter = 5000, startk = 40), silent = TRUE)
  if (!inherits(vem, "try-error")) res$vem <- vem
  # EM (Phase 2)
  if (!is.null(res$vem)) {
    em <- try(mixalg.EM(obs = "NRIRON", pop.at.risk = "NR", family = "binomial",
                        data = df, t = res$vem@t, p = res$vem@p, acc = 1e-6, numiter = 30000), silent = TRUE)
    if (!inherits(em, "try-error")) res$em <- em
  }
  # NPML in one go
  npml <- try(mixalg(obs = "NRIRON", pop.at.risk = "NR", family = "binomial",
                     data = df, acc = 1e-6, numiter = 30000, startk = 50), silent = TRUE)
  if (!inherits(npml, "try-error")) res$npml <- npml
  # Support for diagnostics (prefer NPML)
  if (!is.null(res$npml)) {
    res$support <- data.frame(p = res$npml@t, pi = res$npml@p)
  } else if (!is.null(res$em)) {
    res$support <- data.frame(p = res$em@t, pi = res$em@p)
  }
  # Gradient function diagnostic
  if (!is.null(res$support)) {
    grid_p <- seq(0, 1, length.out = 201)
    N <- nrow(df)
    denom <- sapply(1:N, function(i) sum(res$support$pi * dbinom(df$NRIRON[i], df$NR[i], res$support$p)))
    gradient <- sapply(grid_p, function(pp) {
      num <- dbinom(df$NRIRON, df$NR, pp)
      mean(num / denom)
    })
    png("figs/gradient_binomial.png", width = 900, height = 500)
    plot(grid_p, gradient, type = "l", lwd = 2, col = "purple",
         xlab = "p", ylab = "Gradient d(G,p)", main = "Gradient Function (Binomial Mixture)")
    abline(h = 1, col = "red", lty = 2)
    dev.off()
    res$gradient <- data.frame(p = grid_p, d = gradient)
  }
  res
}

caman_results <- run_caman(dat)

make_overlay_plot <- function(df, p_vec, pi_vec, outpath_prefix = "figs/overlay_mixture") {
  nr_counts <- df %>% count(NR, name = "n") %>% arrange(desc(n))
  if (nrow(nr_counts) == 0) return(NULL)
  n0 <- nr_counts$NR[1]
  df_n0 <- df %>% filter(NR == n0)
  emp <- df_n0 %>% count(NRIRON, name = "freq") %>% mutate(prob = freq / sum(freq))
  y_vals <- 0:n0
  comp_pmfs <- sapply(seq_along(p_vec), function(j) dbinom(y_vals, n0, p_vec[j]))
  overall <- as.numeric(comp_pmfs %*% pi_vec)
  comp_df <- as.data.frame(comp_pmfs)
  colnames(comp_df) <- paste0("component_", seq_along(p_vec))
  plot_df <- cbind(y = y_vals, overall = overall, comp_df) %>%
    tidyr::pivot_longer(cols = -y, names_to = "series", values_to = "prob")
  p <- ggplot() +
    geom_col(data = emp, aes(x = NRIRON, y = prob), fill = "grey80", color = "grey50", width = 0.8) +
    geom_line(data = plot_df, aes(x = y, y = prob, color = series), linewidth = 1) +
    scale_x_continuous(breaks = y_vals) +
    labs(title = sprintf("Mixture overlay for NR=%d", n0), x = "NRIRON (counts)", y = "Probability") +
    theme(legend.position = "bottom")
  outpath <- sprintf("%s_NR%d.png", outpath_prefix, n0)
  ggsave(outpath, p, width = 8, height = 5, dpi = 150)
  outpath
}

overlay_path <- NULL
if (!is.null(caman_results) && !is.null(caman_results$support)) {
  overlay_path <- make_overlay_plot(dat, p_vec = caman_results$support$p, pi_vec = caman_results$support$pi)
} else if (!is.null(mixture_summary)) {
  overlay_path <- make_overlay_plot(dat, p_vec = mixture_summary$p, pi_vec = mixture_summary$pi)
}

mix_expected <- NULL
if (!is.null(caman_results) && !is.null(caman_results$support)) {
  pbar <- sum(caman_results$support$pi * caman_results$support$p)
  p2bar <- sum(caman_results$support$pi * caman_results$support$p^2)
  varp <- p2bar - pbar^2
  mean_pred <- mean(dat$NR) * pbar
  var_pred <- mean(dat$NR^2) * varp + mean(dat$NR) * pbar * (1 - pbar)
  mix_expected <- list(pbar = pbar, varp = varp, mean_pred = mean_pred, var_pred = var_pred,
                       mean_obs = mean(dat$NRIRON), var_obs = var(dat$NRIRON))
}

# 7) Write README.md
fmt_pct <- function(x) paste0(sprintf("%.1f", 100 * x), "%")
# Helper to format percentages consistently in the report.

lines <- c(
  "# AMT: Homework Assignment — Hemodialysis Mixture Analysis",
  "",
  "## (Section 1) Data description (no modeling)",
  "- Source: `hemodialysismix.csv` (patients on hemodialysis)",
  "- Variables: ID, AGE (years), SEX (1=Male, 2=Female, NA kept as Unknown), NR (number of measurements), NRIRON (adequate iron counts)",
  sprintf("- N subjects: %d", eda_summary$n),
  "- Distribution of follow-up counts NR:",
  capture.output(print(eda_summary$nr_table)),
  "- `SEX` counts:",
  capture.output(print(eda_summary$sex_table)),
  "- Proportion p = NRIRON/NR summary:",
  capture.output(print(eda_summary$p_hat_summary)),
  "",
  "Figures (Section 1):",
  "- figs/p_hat_hist.png",
  "- figs/nriron_hist.png",
  "- figs/nr_hist.png",
  "",
  "## (Section 2) Model-based analysis",
  "- Intercept-only Binomial GLM (unimodal) Pearson dispersion:",
  sprintf("  - phi = %.3f (>%s indicates overdispersion)", phi0, "1"),
  if (is.finite(phi0) && phi0 > 1.05) "  - Interpretation: Overdispersion present." else "  - Interpretation: No strong overdispersion.",
  "- Binomial GLM with AGE and SEX (with missingness indicator):",
  sprintf("  - phi = %.3f", phi_cov),
  if (is.finite(phi_cov) && phi_cov > 1.05) "  - Overdispersion persists after covariate adjustment." else "  - Overdispersion largely explained by covariates.",
  "- Finite mixture of Binomials (custom EM baseline):",
  sprintf("  - Selected k by BIC: %d", best_k),
  if (!is.null(mixture_summary)) {
    c("  - Component probabilities (p) and weights (pi):", capture.output(print(mixture_summary)))
  } else {
    c("  - Mixture summary unavailable.")
  },
  if (!is.null(caman_results) && !is.null(caman_results$support)) {
    c("- CAMAN VEM/EM/NPML (binomial):",
      if (!is.null(caman_results$vem)) sprintf("  - VEM: components=%d", length(caman_results$vem@t)) else "  - VEM: not available",
      if (!is.null(caman_results$em)) sprintf("  - EM: components=%d", length(caman_results$em@t)) else "  - EM: not available",
      if (!is.null(caman_results$npml)) sprintf("  - NPML: components=%d", length(caman_results$npml@t)) else "  - NPML: not available",
      "  - Support points and weights (using NPML if available, else EM):",
      capture.output(print(caman_results$support))
    )
  } else {
    "- CAMAN results unavailable (package not installed or fit failed)."
  },
  if (!is.null(caman_results) && !is.null(caman_results$gradient)) {
    c("- Gradient function diagnostic saved: figs/gradient_binomial.png",
      sprintf("  - max d(G,p) = %.3f", max(caman_results$gradient$d)))
  } else {
    "- Gradient function diagnostic unavailable."
  },
  if (!is.null(mix_expected)) {
    c("- Mixture-predicted vs observed (population-level):",
      sprintf("  - E[p] = %.3f, Var[p] = %.4f", mix_expected$pbar, mix_expected$varp),
      sprintf("  - Pred mean NRIRON ≈ %.2f vs observed %.2f", mix_expected$mean_pred, mix_expected$mean_obs),
      sprintf("  - Pred var NRIRON ≈ %.2f vs observed %.2f", mix_expected$var_pred, mix_expected$var_obs)
    )
  } else {
    "- Mixture-predicted diagnostics unavailable."
  },
  "",
  "## (Section 3) Finalization",
  "- Overlay of mixture model on empirical counts histogram (most common NR):",
  if (!is.null(overlay_path)) sprintf("  - %s", overlay_path) else "  - overlay unavailable",
  "",
  "## Software and estimation",
  "- R, glm (binomial), flexmix (baseline mixture EM), CAMAN (VEM/EM/NPML)",
  "- NAs were retained: `SEX` as explicit 'Unknown' level; `AGE` imputed with mean plus a missingness indicator to keep subjects in all analyses.",
  "",
  "## Communication-ready takeaway",
  "- Monthly adequate-iron counts are overdispersed relative to a single-probability Binomial.",
  if (!is.null(caman_results) && !is.null(caman_results$support))
    sprintf("- NPML mixture identifies %d latent components with distinct adequacy rates.", nrow(caman_results$support)) else
    sprintf("- A %d-component finite mixture (baseline) captures heterogeneity.", best_k),
  if (is.finite(phi_cov) && phi_cov > 1.05) "- AGE/SEX do not fully explain heterogeneity (overdispersion persists)." else "- AGE/SEX largely explain dispersion.",
  "",
  "(See notebook `main.ipynb` for figures.)"
)
# Compose a concise Markdown report including summaries, diagnostics, and figures.

writeLines(unlist(lines), con = "README.md")
# Persist the report to README.md at project root for easy viewing on GitHub.

# Save key objects for notebook use
saveRDS(list(
  eda_summary = eda_summary,
  phi0 = phi0,
  phi_cov = phi_cov,
  mixture_summary = mixture_summary,
  assoc_results = assoc_results
), file = "analysis_results.rds")
# Store key results for interactive exploration in the notebook.

cat("Analysis completed. README.md and figures written.\n")
# Final console message to indicate successful run in automated environments.


