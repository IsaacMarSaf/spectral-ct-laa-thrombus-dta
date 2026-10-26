# Post hoc analyses (labeled post hoc in the article and in Online Resource 1); run after 12_dta_meta_analysis.R.
#
# 1. Comparacion pareada (S1) con un efecto prueba-por-estudio: la ventaja espectral puede variar entre estudios.
# 2. S1 sin Hur 2012 (comparador convencional reconstruido por los revisores) y sin Hur 2012 ni Yu 2024 (lectura espectral
#    anidada en positivos convencionales).
# 3. Analisis primario sin Schlett 2020 (intervalo de prediccion de especificidad), sin Li 2022 y sin Yu 2024; con los
#    analisis de sensibilidad planeados (solo un gate: sin Hur 2012; sin verificacion parcial: sin Rodriguez-Granillo
#    2021), cada estudio primario sale una vez.
# 4. Modelo primario con priors vagos reajustado con adapt_delta 0.9999 para eliminar transiciones divergentes.
# 5. Probabilidad post-test y consecuencias por 1000 defectos de llenado con las muestras del modelo S2.
# 6. AUC binormal equivalente a las SMD agrupadas (Phi(|d| / sqrt(2)), varianzas iguales).
#
# Usage: Rscript code/12b_posthoc_analyses.R <repository root>   -> results/posthoc_analyses.csv
suppressPackageStartupMessages({
  library(brms); library(posterior); library(dplyr); library(readr); library(tidyr); library(tibble)
})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
options(brms.backend = "cmdstanr", mc.cores = 4)
SEED <- 20261025
fit_dir <- file.path(root, "fits")
d <- read_csv(file.path(root, "data", "analysis", "dta_2x2.csv"), show_col_types = FALSE)
to_long <- function(x) bind_rows(
  x %>% transmute(study, measure = "sens", y = tp, n = tp + fn),
  x %>% transmute(study, measure = "spec", y = tn, n = tn + fp)) %>% mutate(measure = factor(measure, c("sens", "spec")))
PRI <- c(prior(normal(0, 2), class = b), prior(normal(0, 1), class = sd), prior(lkj(2), class = cor))
CTRL <- list(adapt_delta = 0.995, max_treedepth = 12)
q <- function(v) unname(quantile(v, c(.5, .025, .975)))
f1 <- function(v, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), v[1], v[2], v[3])
out <- list()
add <- function(analysis, quantity, est, note = "") out[[length(out) + 1]] <<- tibble(analysis, quantity, estimate = est, note)
ndiv <- function(fit) { np <- nuts_params(fit); sum(np$Value[np$Parameter == "divergent__"]) }
# convergencia con el mismo criterio que diag_fit() del paso 12: R-hat maximo y ESS bulk minimo de efectos fijos y aleatorios
conv <- function(fit, tag) {
  s <- summary(fit)
  rh <- c(s$fixed$Rhat, unlist(lapply(s$random, function(r) r$Rhat)))
  es <- c(s$fixed$Bulk_ESS, unlist(lapply(s$random, function(r) r$Bulk_ESS)))
  add(tag, "Maximum R-hat", sprintf("%.4f", max(rh, na.rm = TRUE)))
  add(tag, "Minimum bulk ESS", sprintf("%.0f", min(es, na.rm = TRUE)))
  add(tag, "Divergent transitions", as.character(ndiv(fit)))
}

prim <- d %>% filter(set == "PRIMARY"); s1c <- d %>% filter(set == "S1"); s2 <- d %>% filter(set == "S2")
s1s <- prim %>% filter(paper_id %in% s1c$paper_id)
paired <- function(ids) bind_rows(to_long(s1s %>% filter(paper_id %in% ids)) %>% mutate(test = "spectral"),
                                  to_long(s1c %>% filter(paper_id %in% ids)) %>% mutate(test = "conventional")) %>%
  mutate(test = factor(test, c("spectral", "conventional")),
         xse = as.numeric(measure == "sens" & test == "spectral"), xsp = as.numeric(measure == "spec" & test == "spectral"))
summ_cmp <- function(fit, tag) {
  dr <- as_draws_df(fit)
  se_s <- plogis(dr$`b_measuresens:testspectral`); se_c <- plogis(dr$`b_measuresens:testconventional`)
  sp_s <- plogis(dr$`b_measurespec:testspectral`); sp_c <- plogis(dr$`b_measurespec:testconventional`)
  add(tag, "Sensitivity difference, pp", f1(100 * q(se_s - se_c)), sprintf("P(spectral > conventional) = %.2f", mean(se_s > se_c)))
  add(tag, "Specificity difference, pp", f1(100 * q(sp_s - sp_c)), sprintf("P(spectral > conventional) = %.2f", mean(sp_s > sp_c)))
  add(tag, "Specificity, spectral %", f1(100 * q(sp_s))); add(tag, "Specificity, conventional %", f1(100 * q(sp_c)))
  conv(fit, tag)
  dr
}
all4 <- unique(s1c$paper_id)

# ---- 1. efecto prueba-por-estudio: u_study (sens, spec) comun a ambas pruebas + v_study (sens, spec) solo en la espectral
L <- paired(all4)
f_var <- brm(y | trials(n) ~ 0 + measure:test + (0 + measure | study) + (0 + xse + xsp || study), data = L,
             family = binomial(), prior = PRI, chains = 4, iter = 6000, warmup = 2000, seed = SEED, refresh = 0, silent = 2,
             control = CTRL, file = file.path(fit_dir, "cmp_s1_studyvarying"), file_refit = "on_change")
dr <- summ_cmp(f_var, "S1, study-varying spectral effect")
set.seed(SEED)
nd <- nrow(dr)
z <- matrix(rnorm(nd * 3), nd)
rho <- dr$cor_study__measuresens__measurespec
u_sp <- dr$sd_study__measurespec * (rho * z[, 1] + sqrt(1 - rho^2) * z[, 2])
v_sp <- dr$sd_study__xsp * z[, 3]
new_diff <- plogis(dr$`b_measurespec:testspectral` + u_sp + v_sp) - plogis(dr$`b_measurespec:testconventional` + u_sp)
add("S1, study-varying spectral effect", "Specificity difference in a new study, pp (95% PI)", f1(100 * q(new_diff)),
    sprintf("P(new study favors conventional) = %.2f", mean(new_diff < 0)))
add("S1, study-varying spectral effect", "SD of spectral effect on logit specificity", f1(q(dr$sd_study__xsp), 2))

# ---- 2. S1 dejando fuera estudios (modelo original, efectos de estudio compartidos)
fit_cmp <- function(ids, tag, file) {
  f <- brm(y | trials(n) ~ 0 + measure:test + (0 + measure | study), data = paired(ids), family = binomial(), prior = PRI,
           chains = 4, iter = 6000, warmup = 2000, seed = SEED, refresh = 0, silent = 2, control = CTRL,
           file = file.path(fit_dir, file), file_refit = "on_change")
  summ_cmp(f, tag)
}
fit_cmp(setdiff(all4, "Hur_2012"), "S1 without Hur et al. 2012", "cmp_s1_noHur")
fit_cmp(setdiff(all4, c("Hur_2012", "Yu_2024")), "S1 without Hur et al. 2012 and Yu et al. 2024", "cmp_s1_noHurYu")

# ---- 3. primario sin Schlett 2020
summ_biv <- function(fit, tag) {
  dr <- as_draws_df(fit)
  ms <- dr$b_measuresens; mp <- dr$b_measurespec
  ss <- dr$sd_study__measuresens; sp <- dr$sd_study__measurespec; rho <- dr$cor_study__measuresens__measurespec
  set.seed(SEED)
  z1 <- rnorm(length(ms)); z2 <- rnorm(length(ms))
  ns <- ms + ss * z1; np <- mp + sp * (rho * z1 + sqrt(1 - rho^2) * z2)
  add(tag, "Sensitivity %", f1(100 * q(plogis(ms)))); add(tag, "Specificity %", f1(100 * q(plogis(mp))))
  add(tag, "Sensitivity 95% PI %", f1(100 * q(plogis(ns)))); add(tag, "Specificity 95% PI %", f1(100 * q(plogis(np))))
  conv(fit, tag)
  list(se = plogis(ms), sp = plogis(mp))
}
f_noS <- brm(y | trials(n) ~ 0 + measure + (0 + measure | study), data = to_long(prim %>% filter(paper_id != "Schlett_2020")),
             family = binomial(), prior = PRI, chains = 4, iter = 6000, warmup = 2000, seed = SEED, refresh = 0, silent = 2,
             control = CTRL, file = file.path(fit_dir, "biv_sa_noSchlett_main"), file_refit = "on_change")
summ_biv(f_noS, "Primary without Schlett et al. 2020")
for (x in list(c("Li_2022", "Primary without Li et al. 2022", "biv_sa_noLi2022_main"),
               c("Yu_2024", "Primary without Yu et al. 2024", "biv_sa_noYu_main"))) {
  f_x <- brm(y | trials(n) ~ 0 + measure + (0 + measure | study), data = to_long(prim %>% filter(paper_id != x[1])),
             family = binomial(), prior = PRI, chains = 4, iter = 6000, warmup = 2000, seed = SEED, refresh = 0, silent = 2,
             control = CTRL, file = file.path(fit_dir, x[3]), file_refit = "on_change")
  summ_biv(f_x, x[2])
}

# ---- 4. priors vagos con adapt_delta mas alto
PRI_V <- c(prior(normal(0, 5), class = b), prior(cauchy(0, 1), class = sd), prior(lkj(1), class = cor))
f_vague <- brm(y | trials(n) ~ 0 + measure + (0 + measure | study), data = to_long(prim), family = binomial(), prior = PRI_V,
               chains = 4, iter = 6000, warmup = 2000, seed = SEED, refresh = 0, silent = 2,
               control = list(adapt_delta = 0.9999, max_treedepth = 15),
               file = file.path(fit_dir, "biv_primary_vague_ad9999"), file_refit = "on_change")
summ_biv(f_vague, "Primary, vague priors, adapt_delta 0.9999")

# ---- 5. S2: post-test y consecuencias por 1000 defectos
f_s2 <- readRDS(file.path(fit_dir, "biv_s2_main.rds"))
dr2 <- as_draws_df(f_s2)
se2 <- plogis(dr2$b_measuresens); sp2 <- plogis(dr2$b_measurespec)
for (pre in c(0.27, 0.31, 0.36)) {
  odds <- pre / (1 - pre)
  lrn <- (1 - se2) / sp2; lrp <- se2 / (1 - sp2)
  pn <- odds * lrn / (1 + odds * lrn); pp <- odds * lrp / (1 + odds * lrp)
  tag <- sprintf("S2 at %.0f%% thrombus among defects", 100 * pre)
  add(tag, "Posttest probability after negative spectral read, %", f1(100 * q(pn)))
  add(tag, "Posttest probability after positive spectral read, %", f1(100 * q(pp)))
  add(tag, "Per 1000 defects: TP", f1(q(1000 * pre * se2), 0)); add(tag, "Per 1000 defects: FN", f1(q(1000 * pre * (1 - se2)), 0))
  add(tag, "Per 1000 defects: TN", f1(q(1000 * (1 - pre) * sp2), 0)); add(tag, "Per 1000 defects: FP", f1(q(1000 * (1 - pre) * (1 - sp2)), 0))
}
add("S2", "LR+", f1(q(se2 / (1 - sp2)), 1)); add("S2", "LR-", f1(q((1 - se2) / sp2), 2))

# ---- 6. AUC binormal de las SMD
for (x in list(c("Thrombus vs stasis", -2.09, -4.24, 0.07), c("Thrombus vs normal LAA", -4.04, -7.14, -0.94))) {
  v <- as.numeric(x[2:4])  # SMD (thrombus - comparator), lower and upper CI limits
  a <- pnorm(-v / sqrt(2))  # AUC = P(thrombus value < comparator value) under equal-variance normality
  add("S3, binormal AUC", x[1], sprintf("%.2f (%.2f to %.2f)", a[1], a[3], a[2]), "equal-variance binormal; CI bounds from the SMD CI")
}

res <- bind_rows(out)
write_csv(res, file.path(root, "results", "posthoc_analyses.csv"))
print(res, n = 100, width = 200)
