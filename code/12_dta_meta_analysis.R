# Paso 11: meta-analisis de exactitud diagnostica (DTA) y figuras.
#
# Modelo primario: bivariado binomial-normal (Chu y Cole 2006; Reitsma 2005) estimado en forma bayesiana con brms/Stan,
# verosimilitud binomial exacta (sin correccion por continuidad) y priors debilmente informativos:
#   logit(Se), logit(Sp) ~ Normal(0, 2); DE entre estudios ~ Half-Normal(0, 1); correlacion ~ LKJ(2).
# Sensibilidad: priors alternativos, GLMM bivariado frecuentista (lme4), Reitsma (mada, correccion 0.5) y modelos
# univariados GLMM (meta). Comparativo (S1): misma estructura con la prueba como covariable y efectos de estudio
# compartidos (diseno pareado). Continuos (S3): SMD (Hedges g), REML con Knapp-Hartung (metafor).
#
# Uso: Rscript scripts/12_dta_meta_analysis.R <carpeta Mi--Review>
suppressPackageStartupMessages({
  library(brms); library(lme4); library(mada); library(meta); library(metafor); library(posterior)
  library(dplyr); library(tidyr); library(readr); library(purrr); library(tibble)
  library(ggplot2); library(patchwork); library(ggdist); library(distributional); library(tidybayes); library(grid); library(ragg)
  library(ggrepel); library(systemfonts)
  library(jsonlite); library(writexl)
})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
dat_dir <- file.path(root, "data", "analysis")
res_dir <- file.path(root, "results")
# FIG_DIR, FIG_TAG y FIG_DPI permiten una segunda tanda de figuras para otra revista (p. ej. IJCI: paneles a, b y 822 dpi
# a 5 pulgadas, que equivalen a 600 dpi a 174 mm); sin variables de entorno se generan las figuras del AJR.
fig_dir <- Sys.getenv("FIG_DIR", file.path(root, "figures"))
FIG_TAG <- Sys.getenv("FIG_TAG", "A")
FIG_DPI <- as.numeric(Sys.getenv("FIG_DPI", "600"))
fit_dir <- file.path(root, "fits")
for (d in c(res_dir, fig_dir, fit_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
options(brms.backend = "cmdstanr", mc.cores = 4)
SEED <- 20261025
set.seed(SEED)

# ------------------------------------------------------------------ estilo (paleta validada, ver dataviz) ----------
COL <- c(spectral = "#2a78d6", conventional = "#eb6834", pcd = "#1baf7a")
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e0"
STATUS <- c(low = "#0ca30c", insufficient = "#fab219", high = "#d03b3b")
FONT <- "Arial"
theme_ajr <- function(base = 7) {
  theme_minimal(base_size = base, base_family = FONT) %+replace%
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.25),
          axis.text = element_text(colour = INK2, size = base - 0.5), axis.title = element_text(colour = INK),
          axis.line = element_line(colour = INK2, linewidth = 0.25), axis.ticks = element_line(colour = INK2, linewidth = 0.25),
          legend.text = element_text(colour = INK, size = base - 0.5), legend.title = element_text(colour = INK, size = base - 0.5),
          plot.tag = element_text(face = "bold", size = base + 2, colour = INK), plot.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(colour = INK, face = "bold", size = base))
}
save_gg <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, units = "in", dpi = FIG_DPI, device = ragg::agg_png, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, "_rgb300.tif")), p, width = w, height = h, units = "in", dpi = 300,
         device = ragg::agg_tiff, bg = "white", compression = "lzw")
}
save_grob <- function(g, name, w, h) {
  cairo_pdf(file.path(fig_dir, paste0(name, ".pdf")), width = w, height = h, family = FONT); grid.newpage(); grid.draw(g); dev.off()
  ragg::agg_png(file.path(fig_dir, paste0(name, ".png")), width = w, height = h, units = "in", res = FIG_DPI, background = "white")
  grid.newpage(); grid.draw(g); dev.off()
  ragg::agg_tiff(file.path(fig_dir, paste0(name, "_rgb300.tif")), width = w, height = h, units = "in", res = 300,
                 background = "white", compression = "lzw")
  grid.newpage(); grid.draw(g); dev.off()
}
fmt3 <- function(q, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), q[1], q[2], q[3])
pct3 <- function(q) sprintf("%.1f%% (%.1f–%.1f)", 100 * q[1], 100 * q[2], 100 * q[3])

# ------------------------------------------------------------------ datos -------------------------------------------
d <- read_csv(file.path(dat_dir, "dta_2x2.csv"), show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
cont <- read_csv(file.path(dat_dir, "continuous.csv"), show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
q3 <- read_csv(file.path(dat_dir, "quadas3.csv"), show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
prisma <- fromJSON(file.path(dat_dir, "prisma.json"))

to_long <- function(x, extra = character()) {
  bind_rows(
    x %>% transmute(study, across(all_of(extra)), measure = "sens", y = tp, n = tp + fn),
    x %>% transmute(study, across(all_of(extra)), measure = "spec", y = tn, n = tn + fp)
  ) %>% mutate(measure = factor(measure, levels = c("sens", "spec")))
}
exact_ci <- function(x, n) { b <- binom.test(x, n); c(x / n, b$conf.int[1], b$conf.int[2]) }

# ------------------------------------------------------------------ modelos -----------------------------------------
PRI <- list(
  main  = c(prior(normal(0, 2), class = b), prior(normal(0, 1), class = sd), prior(lkj(2), class = cor)),
  vague = c(prior(normal(0, 5), class = b), prior(cauchy(0, 1), class = sd), prior(lkj(1), class = cor)),
  tight = c(prior(normal(0, 2), class = b), prior(normal(0, 0.5), class = sd), prior(lkj(2), class = cor))
)
fit_biv <- function(x, prior = "main", tag) {
  brm(y | trials(n) ~ 0 + measure + (0 + measure | study), data = to_long(x), family = binomial(),
      prior = PRI[[prior]], chains = 4, iter = 6000, warmup = 2000, seed = SEED, refresh = 0, silent = 2,
      control = list(adapt_delta = 0.995, max_treedepth = 12),
      file = file.path(fit_dir, paste0("biv_", tag, "_", prior)), file_refit = "on_change")
}
diag_fit <- function(fit, tag) {
  s <- summary(fit)
  rh <- c(s$fixed$Rhat, unlist(lapply(s$random, function(r) r$Rhat)))
  es <- c(s$fixed$Bulk_ESS, unlist(lapply(s$random, function(r) r$Bulk_ESS)))
  np <- nuts_params(fit)
  tibble(model = tag, max_rhat = max(rh, na.rm = TRUE), min_bulk_ess = min(es, na.rm = TRUE),
         divergences = sum(np$Value[np$Parameter == "divergent__"]))
}
summ_biv <- function(fit) {
  dr <- as_draws_df(fit)
  ms <- dr$b_measuresens; mp <- dr$b_measurespec
  ss <- dr$sd_study__measuresens; sp <- dr$sd_study__measurespec; rho <- dr$cor_study__measuresens__measurespec
  set.seed(SEED)  # la prediccion para un estudio nuevo no depende del orden en que se cargan los modelos
  z1 <- rnorm(length(ms)); z2 <- rnorm(length(ms))
  ns <- ms + ss * z1; np <- mp + sp * (rho * z1 + sqrt(1 - rho^2) * z2)
  se <- plogis(ms); spc <- plogis(mp)
  lrp <- se / (1 - spc); lrn <- (1 - se) / spc
  q <- function(v) unname(quantile(v, c(.5, .025, .975)))
  list(sens = q(se), spec = q(spc), sens_pi = q(plogis(ns)), spec_pi = q(plogis(np)), lrp = q(lrp), lrn = q(lrn),
       dor = q(lrp / lrn), tau_sens = q(ss), tau_spec = q(sp), rho = q(rho),
       draws = tibble(ms = ms, mp = mp, ns = ns, np = np, se = se, sp = spc))
}
fit_glmm <- function(x) {
  L <- to_long(x)
  m <- tryCatch(suppressWarnings(glmer(cbind(y, n - y) ~ 0 + measure + (0 + measure | study), data = L, family = binomial,
                                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- fixef(m); s <- sqrt(diag(as.matrix(vcov(m))))
  list(sens = plogis(b[1] + c(0, -1.96, 1.96) * s[1]), spec = plogis(b[2] + c(0, -1.96, 1.96) * s[2]),
       note = paste0(if (isSingular(m)) "singular fit" else "regular fit",
                     if (length(m@optinfo$conv$lme4$messages)) paste0("; ", paste(m@optinfo$conv$lme4$messages, collapse = "; ")) else ""))
}
fit_reitsma <- function(x) {
  r <- tryCatch(reitsma(data.frame(TP = x$tp, FN = x$fn, FP = x$fp, TN = x$tn), correction = 0.5, correction.control = "all"),
                error = function(e) NULL)
  if (is.null(r)) return(NULL)
  cf <- summary(r)$coefficients
  list(sens = unname(cf["sensitivity", c(1, 5, 6)]), spec = 1 - unname(cf["false pos. rate", c(1, 6, 5)]))
}
fit_uni <- function(e, n) {
  m <- tryCatch(metaprop(event = e, n = n, method = "GLMM", sm = "PLOGIT", common = FALSE, random = TRUE, prediction = TRUE),
                error = function(err) NULL)
  if (is.null(m)) return(rep(NA, 5))
  plogis(c(m$TE.random, m$lower.random, m$upper.random, m$lower.predict, m$upper.predict))
}
row_of <- function(label, x, s, method) {
  tibble(analysis = label, method = method, k = nrow(x), studies = paste(x$study, collapse = "; "), patients = sum(x$n),
         thrombi = sum(x$tp + x$fn), sens = pct3(s$sens), spec = pct3(s$spec),
         sens_pi = if (!is.null(s$sens_pi)) pct3(s$sens_pi) else "", spec_pi = if (!is.null(s$spec_pi)) pct3(s$spec_pi) else "",
         lr_pos = if (!is.null(s$lrp)) fmt3(s$lrp, 1) else "", lr_neg = if (!is.null(s$lrn)) fmt3(s$lrn, 2) else "",
         dor = if (!is.null(s$dor)) sprintf("%.0f (%.0f–%.0f)", s$dor[1], s$dor[2], s$dor[3]) else "",
         tau_sens = if (!is.null(s$tau_sens)) fmt3(s$tau_sens) else "", tau_spec = if (!is.null(s$tau_spec)) fmt3(s$tau_spec) else "",
         rho = if (!is.null(s$rho)) fmt3(s$rho) else "", sens_med = s$sens[1], sens_lo = s$sens[2], sens_hi = s$sens[3],
         spec_med = s$spec[1], spec_lo = s$spec[2], spec_hi = s$spec[3])
}

prim <- d %>% filter(set == "PRIMARY")
alt <- d %>% filter(set == "PRIMARY_ALT")
s1c <- d %>% filter(set == "S1"); s2 <- d %>% filter(set == "S2"); s4 <- d %>% filter(set == "S4")

fits <- list(); summ <- list(); rows <- list(); diags <- list()
run <- function(tag, label, x, prior = "main") {
  f <- fit_biv(x, prior, tag); fits[[paste(tag, prior)]] <<- f; s <- summ_biv(f); summ[[paste(tag, prior)]] <<- s
  diags[[paste(tag, prior)]] <<- diag_fit(f, paste(tag, prior))
  rows[[length(rows) + 1]] <<- row_of(label, x, s, paste0("Bayesian bivariate (", prior, " priors)"))
  s
}
message("PRIMARY")
sP <- run("primary", "Primary: spectral CT, all patients", prim)
run("primary", "Primary: spectral CT, all patients", prim, "vague")
run("primary", "Primary: spectral CT, all patients", prim, "tight")
g <- fit_glmm(prim); if (!is.null(g)) rows[[length(rows) + 1]] <- row_of("Primary: spectral CT, all patients", prim, g, paste("Frequentist bivariate GLMM,", g$note))
r <- fit_reitsma(prim); if (!is.null(r)) rows[[length(rows) + 1]] <- row_of("Primary: spectral CT, all patients", prim, r, "Reitsma bivariate (continuity correction 0.5)")
u_se <- fit_uni(prim$tp, prim$tp + prim$fn); u_sp <- fit_uni(prim$tn, prim$tn + prim$fp)
rows[[length(rows) + 1]] <- row_of("Primary: spectral CT, all patients", prim,
                                   list(sens = u_se[1:3], spec = u_sp[1:3], sens_pi = u_se[c(1, 4, 5)], spec_pi = u_sp[c(1, 4, 5)]),
                                   "Univariate GLMM (separate random-effects logit models)")
message("SENSITIVITY ANALYSES")
run("sa_li2019", "Sensitivity: Li et al. 2019 instead of Li et al. 2022", bind_rows(prim %>% filter(paper_id != "Li_2022"), alt))
run("sa_singlegate", "Sensitivity: single-gate studies only", prim %>% filter(gate == "single-gate"))
run("sa_noRG", "Sensitivity: excluding derived LAA counts with partial verification", prim %>% filter(paper_id != "Rodriguez-Granillo_2021"))
run("sa_af", "Sensitivity: AF pre-procedure cohorts only", prim %>% filter(paper_id %in% c("Li_2022", "Schlett_2020", "Yu_2024")))
message("S2")
sS2 <- run("s2", "S2: thrombus vs SEC in filling defects", s2)
run("s2_li2019", "S2 sensitivity: Li et al. 2019 instead of Li et al. 2022", bind_rows(s2 %>% filter(paper_id != "Li_2022"), d %>% filter(set == "S2_ALT")))
run("s2_zhangIC", "S2 sensitivity: Zhang et al. 2024 with absolute iodine concentration", bind_rows(s2 %>% filter(paper_id != "Zhang_2024"), d %>% filter(set == "S2_IC")))

# ---- S1: comparacion pareada espectral vs convencional ----
message("S1")
s1s <- prim %>% filter(paper_id %in% s1c$paper_id)
L1 <- bind_rows(to_long(s1s) %>% mutate(test = "spectral"), to_long(s1c) %>% mutate(test = "conventional")) %>%
  mutate(test = factor(test, levels = c("spectral", "conventional")))
fit_cmp <- brm(y | trials(n) ~ 0 + measure:test + (0 + measure | study), data = L1, family = binomial(), prior = PRI$main,
               chains = 4, iter = 6000, warmup = 2000, seed = SEED, refresh = 0, silent = 2,
               control = list(adapt_delta = 0.995, max_treedepth = 12), file = file.path(fit_dir, "cmp_s1_main"), file_refit = "on_change")
diags[["s1 comparative"]] <- diag_fit(fit_cmp, "S1 comparative")
dc <- as_draws_df(fit_cmp)
qq <- function(v) unname(quantile(v, c(.5, .025, .975)))
se_s <- plogis(dc$`b_measuresens:testspectral`); se_c <- plogis(dc$`b_measuresens:testconventional`)
sp_s <- plogis(dc$`b_measurespec:testspectral`); sp_c <- plogis(dc$`b_measurespec:testconventional`)
s1_tab <- tibble(
  quantity = c("Sensitivity, spectral", "Sensitivity, conventional", "Difference in sensitivity (spectral − conventional)",
               "Specificity, spectral", "Specificity, conventional", "Difference in specificity (spectral − conventional)",
               "Relative specificity (spectral / conventional)"),
  estimate = c(pct3(qq(se_s)), pct3(qq(se_c)), sprintf("%.1f (%.1f to %.1f) percentage points", 100 * qq(se_s - se_c)[1], 100 * qq(se_s - se_c)[2], 100 * qq(se_s - se_c)[3]),
               pct3(qq(sp_s)), pct3(qq(sp_c)), sprintf("%.1f (%.1f to %.1f) percentage points", 100 * qq(sp_s - sp_c)[1], 100 * qq(sp_s - sp_c)[2], 100 * qq(sp_s - sp_c)[3]),
               fmt3(qq(sp_s / sp_c), 3)),
  posterior_probability = c(NA, NA, mean(se_s > se_c), NA, NA, mean(sp_s > sp_c), mean(sp_s > sp_c)))
s1_tab$k <- nrow(s1c); s1_tab$studies <- paste(s1c$study, collapse = "; ")

# ---- S3: datos continuos (SMD) ----
message("S3")
s3_one <- function(cmp) {
  x <- cont %>% filter(comparison == cmp)
  es <- escalc("SMD", m1i = m1, sd1i = sd1, n1i = n1, m2i = m2, sd2i = sd2, n2i = n2, data = x)
  m <- rma(yi, vi, data = es, method = "REML", test = "knha")
  p <- predict(m)
  list(es = es %>% mutate(lo = yi - 1.96 * sqrt(vi), hi = yi + 1.96 * sqrt(vi), w = weights(m)), m = m, p = p)
}
s3a <- s3_one("thrombus_vs_SEC"); s3b <- s3_one("thrombus_vs_normal")
s3_tab <- bind_rows(lapply(list(`Thrombus vs SEC` = s3a, `Thrombus vs normal LAA` = s3b), function(o)
  tibble(k = o$m$k, smd = sprintf("%.2f (%.2f to %.2f)", o$p$pred, o$p$ci.lb, o$p$ci.ub),
         prediction_interval = sprintf("%.2f to %.2f", o$p$pi.lb, o$p$pi.ub), tau2 = round(o$m$tau2, 3), I2 = round(o$m$I2, 1))),
  .id = "comparison")

# ---- S4: un solo estudio ----
s4_tab <- s4 %>% rowwise() %>% mutate(sens = pct3(exact_ci(tp, tp + fn)), spec = pct3(exact_ci(tn, tn + fp))) %>% ungroup() %>%
  select(study, test_label, tp, fp, fn, tn, sens, spec)

# ---- consecuencias por 1000 (GRADE) ----
prev_obs <- prim %>% filter(gate == "single-gate") %>% pull(prevalence)
prevs <- c(min(prev_obs), max(prev_obs))
cons <- bind_rows(lapply(prevs, function(pv) {
  dr <- sP$draws; q <- function(v) sprintf("%.0f (%.0f–%.0f)", quantile(v, .5), quantile(v, .025), quantile(v, .975))
  tibble(prevalence = sprintf("%.1f%%", 100 * pv), per_1000_with_thrombus = round(1000 * pv),
         true_positives = q(1000 * pv * dr$se), false_negatives = q(1000 * pv * (1 - dr$se)),
         true_negatives = q(1000 * (1 - pv) * dr$sp), false_positives = q(1000 * (1 - pv) * (1 - dr$sp)))
}))

# ------------------------------------------------------------------ tablas de resultados ----------------------------
pooled <- bind_rows(rows)
diagnostics <- bind_rows(diags)
per_study <- d %>% filter(set %in% c("PRIMARY", "PRIMARY_ALT", "S1", "S1_ALT", "S1_FD", "S2", "S2_ALT", "S2_IC", "S4")) %>% rowwise() %>%
  mutate(sensitivity = sprintf("%d/%d (%.1f%%)", tp, tp + fn, 100 * tp / (tp + fn)),
         specificity = sprintf("%d/%d (%.1f%%)", tn, tn + fp, 100 * tn / (tn + fp)),
         sens_ci = pct3(exact_ci(tp, tp + fn)), spec_ci = pct3(exact_ci(tn, tn + fp))) %>% ungroup()
write_csv(pooled, file.path(res_dir, "pooled_estimates.csv"))
write_csv(diagnostics, file.path(res_dir, "model_diagnostics.csv"))
write_csv(s1_tab, file.path(res_dir, "s1_spectral_vs_conventional.csv"))
write_csv(s3_tab, file.path(res_dir, "s3_continuous_smd.csv"))
write_csv(s4_tab, file.path(res_dir, "s4_photon_counting.csv"))
write_csv(cons, file.path(res_dir, "consequences_per_1000.csv"))
write_csv(per_study, file.path(res_dir, "per_study_accuracy.csv"))
write_xlsx(list(pooled = pooled, diagnostics = diagnostics, s1 = s1_tab, s3 = s3_tab, s4 = s4_tab, consequences = cons,
                per_study = per_study, s3_effects_vs_SEC = as.data.frame(s3a$es), s3_effects_vs_normal = as.data.frame(s3b$es)),
           file.path(res_dir, "meta_analysis_results.xlsx"))

# ------------------------------------------------------------------ forest plots (ggdist) ---------------------------
message("FIGURES")
# Forest plots de lienzo unico: tabla y paneles en un solo ggplot con el eje x en pulgadas (ancho 5 in), asi las bandas
# alternas, los encabezados agrupados y las reglas son continuos. Cada fila lleva un "eye" de ggdist (densidad
# simetrica): en DTA, la posterior especifica del estudio en el modelo bivariado (con shrinkage); en SMD, la
# distribucion muestral normal. Encima van la estimacion observada y su IC. El resumen es la posterior del modelo (DTA)
# o la t de Knapp-Hartung con k - 1 gl (SMD); la ultima fila es la distribucion predictiva para un estudio nuevo.
PT_BODY <- 5.6; PT_SUB <- 4.7; PT_AX <- 5.0
ZEBRA <- "#f3f2ee"
EYE <- 0.62                                   # altura total del eye (unidades de fila)
fs <- function(pt) pt / ggplot2::.pt
sw <- function(s, pt = PT_BODY, bold = FALSE) {
  s <- as.character(s); s[is.na(s)] <- ""
  vapply(s, function(z) if (!nzchar(z)) 0 else as.numeric(systemfonts::string_width(z, family = FONT, size = pt, weight = if (bold) "bold" else "normal", res = 72)) / 72, numeric(1))
}
map_x <- function(v, lim, reg) reg[1] + (v - lim[1]) / diff(lim) * diff(reg)
minus <- function(s) gsub("-", "\u2212", s, fixed = TRUE)
alpha_hex <- function(col, a) paste0(col, sprintf("%02X", round(255 * a)))
FILL_STUDY <- function(col) alpha_hex(col, 0.42)
FILL_SUMM <- function(col) alpha_hex(col, 0.95)
FILL_PRED <- alpha_hex(INK2, 0.22)

zebra_layer <- function(rows, W) {
  zb <- rows %>% filter(kind == "study") %>% group_by(block) %>% arrange(desc(row), .by_group = TRUE) %>%
    filter(row_number() %% 2 == 1) %>% ungroup()
  if (!nrow(zb)) return(NULL)
  annotate("rect", xmin = 0, xmax = W, ymin = zb$row - 0.5, ymax = zb$row + 0.5, fill = ZEBRA, colour = NA)
}
rules_layer <- function(rows, W) {
  top <- max(rows$row); sm <- rows %>% filter(kind == "summary")
  list(annotate("segment", x = 0, xend = W, y = top + 0.5, yend = top + 0.5, colour = INK, linewidth = 0.4),
       if (nrow(sm)) annotate("segment", x = 0, xend = W, y = sm$row + 0.5, yend = sm$row + 0.5, colour = INK2, linewidth = 0.2))
}
axis_layer <- function(reg, lim, breaks, labels, top, y0 = 0.5, title = NULL) {
  bx <- map_x(breaks, lim, reg)
  list(annotate("segment", x = bx, xend = bx, y = y0, yend = top + 0.5, colour = GRID, linewidth = 0.25),
       annotate("segment", x = reg[1], xend = reg[2], y = y0, yend = y0, colour = INK2, linewidth = 0.25),
       annotate("segment", x = bx, xend = bx, y = y0, yend = y0 - 0.09, colour = INK2, linewidth = 0.25),
       annotate("text", x = bx, y = y0 - 0.15, label = labels, size = fs(PT_AX), colour = INK2, family = FONT, vjust = 1),
       if (!is.null(title)) annotate("text", x = mean(reg), y = y0 - 0.62, label = title, size = fs(PT_AX), colour = INK, family = FONT, vjust = 1))
}
# intervalos con flecha si salen de la region
int_layer <- function(dat, reg, lw, lty = "solid") {
  if (is.null(dat) || !nrow(dat)) return(NULL)
  dat <- dat %>% mutate(cl = lo < reg[1], ch = hi > reg[2], lo = pmax(lo, reg[1]), hi = pmin(hi, reg[2]))
  a <- arrow(length = unit(0.032, "in"), type = "closed", angle = 25)
  list(geom_segment(data = dat, aes(x = lo, xend = hi, y = row, yend = row), colour = INK, linewidth = lw, linetype = lty),
       if (any(dat$cl)) geom_segment(data = dat %>% filter(cl), aes(x = lo + 0.05, xend = lo, y = row, yend = row), colour = INK, linewidth = lw, arrow = a),
       if (any(dat$ch)) geom_segment(data = dat %>% filter(ch), aes(x = hi - 0.05, xend = hi, y = row, yend = row), colour = INK, linewidth = lw, arrow = a))
}
obs_pts <- function(dat) geom_point(data = dat, aes(est, row), shape = 21, fill = "white", colour = INK, size = 1.15, stroke = 0.4)
summ_pts <- function(dat) geom_point(data = dat, aes(est, row), shape = 21, fill = INK, colour = "white", size = 1.5, stroke = 0.35)
txt <- function(x, y, label, hjust = 0, pt = PT_BODY, face = "plain", colour = INK, vjust = 0.5, lineheight = 0.9)
  annotate("text", x = x, y = y, label = label, hjust = hjust, vjust = vjust, size = fs(pt), fontface = face, colour = colour,
           family = FONT, lineheight = lineheight)
# leyenda dibujada (glifos iguales a las marcas)
key_layer <- function(items, color, W, y) {
  gw <- 0.2; gap <- 0.05; sep <- 0.24
  lw <- sw(items$label, PT_AX)
  tot <- sum(gw + gap + lw) + sep * (nrow(items) - 1)
  x0 <- (W - tot) / 2 + cumsum(c(0, head(gw + gap + lw + sep, -1)))
  t <- seq(0, 1, length.out = 61); h <- dbeta(t, 5, 2.4); h <- h / max(h) * EYE / 2
  out <- list()
  for (i in seq_len(nrow(items))) {
    xs <- x0[i] + t * gw; ty <- items$type[i]
    poly <- tibble(x = c(xs, rev(xs)), y = c(y + h, rev(y - h)))
    if (ty == "study") out <- c(out, list(geom_polygon(data = poly, aes(x, y), fill = FILL_STUDY(color), colour = NA)))
    if (ty == "summary") out <- c(out, list(geom_polygon(data = poly, aes(x, y), fill = FILL_SUMM(color), colour = NA),
                                            annotate("segment", x = x0[i] + 0.02, xend = x0[i] + gw - 0.01, y = y, yend = y, colour = INK, linewidth = 0.55),
                                            annotate("point", x = x0[i] + gw * 0.7, y = y, shape = 21, fill = INK, colour = "white", size = 1.5, stroke = 0.35)))
    if (ty == "pred") out <- c(out, list(geom_polygon(data = poly, aes(x, y), fill = FILL_PRED, colour = NA),
                                         annotate("segment", x = x0[i], xend = x0[i] + gw, y = y, yend = y, colour = INK, linewidth = 0.4, linetype = "22")))
    if (ty == "obs") out <- c(out, list(annotate("segment", x = x0[i] + 0.01, xend = x0[i] + gw - 0.01, y = y, yend = y, colour = INK, linewidth = 0.32),
                                        annotate("point", x = x0[i] + gw / 2, y = y, shape = 21, fill = "white", colour = INK, size = 1.15, stroke = 0.4)))
    out <- c(out, list(txt(x0[i] + gw + gap, y, items$label[i], pt = PT_AX, colour = INK2)))
  }
  out
}
canvas <- function(W, ylim) list(coord_cartesian(xlim = c(0, W), ylim = ylim, expand = FALSE, clip = "off"),
                                 scale_fill_identity(), scale_colour_identity(), theme_void(base_family = FONT),
                                 theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = "white", colour = NA)))

pc0 <- function(v) sprintf("%.0f", 100 * v)
pc1 <- function(v) ifelse(round(100 * v, 1) >= 100, "100", sprintf("%.1f", 100 * v))
clean_test <- function(s) gsub("\u2264 ?", "\u2264 ", s)
# rotulos de prueba comunes con las Tablas 2, 3 y S4 (scripts/14_manuscript_tables.py, TEST)
TEST_LAB <- c("IC \u22641.74 mg/mL" = "IC \u2264 1.74 mg/mL", "IC \u22648.9 mg/cm\u00b3" = "IC \u2264 8.9 mg/cm\u00b3",
              "Phase III + IC" = "Phase III defect plus IC", "Iodine map / low-keV reading" = "Iodine maps and low-keV VMI, visual",
              "NIC \u22641.98%" = "NIC (LAA/DA) \u2264 1.98%")
test_lab <- function(s) unname(ifelse(s %in% names(TEST_LAB), TEST_LAB[s], clean_test(s)))
# adquisicion en la que se midio el parametro espectral (igual que la Tabla 2)
ACQ <- c(Hur_2012 = "dual-bolus, delayed", Li_2019 = "dual-bolus", Li_2022 = "dual-bolus", Schlett_2020 = "dual-bolus",
         `Rodriguez-Granillo_2021` = "delayed", Yu_2024 = "delayed", Zhang_2024 = "single-bolus early",
         He_2019 = "single-bolus early", Kirschfink_2026 = "single-bolus early")

# ---------------------------------------------------------------------------------------- forest pareado (DTA) ----
dta_forest <- function(x, fit, s, name, color = COL["spectral"], W = 5, ROW = 0.285) {
  x <- x %>% mutate(year = as.integer(sub(".* ", "", study))) %>% arrange(year, study)
  k <- nrow(x); top <- k + 2
  st <- x %>% rowwise() %>% mutate(se = list(exact_ci(tp, tp + fn)), sp = list(exact_ci(tn, tn + fp))) %>% ungroup() %>%
    mutate(row = top:3, sub = paste0(test_lab(test_label), " \u00b7 ", ACQ[paper_id]))
  rows <- tibble(row = c(st$row, 2, 1), kind = c(rep("study", k), "summary", "pi"), block = 1)
  obs <- function(m) tibble(row = st$row, est = map_dbl(st[[m]], 1), lo = map_dbl(st[[m]], 2), hi = map_dbl(st[[m]], 3))
  qs <- list(se = list(s = s$sens, p = s$sens_pi), sp = list(s = s$spec, p = s$spec_pi))
  lab <- function(m) { o <- obs(m); q <- qs[[m]]
    c(sprintf("%s (%s\u2013%s)", pc1(o$est), pc1(o$lo), pc1(o$hi)), sprintf("%s (%s\u2013%s)", pc1(q$s[1]), pc1(q$s[2]), pc1(q$s[3])),
      sprintf("%s\u2013%s", pc1(q$p[2]), pc1(q$p[3]))) }
  L <- list(se = lab("se"), sp = lab("sp"))
  # ---- columnas (pulgadas), calculadas desde el ancho real del texto
  sum_lab <- c("Summary estimate", "95% prediction interval")
  w_study <- max(sw(st$study), sw(st$sub, PT_SUB), sw(sum_lab[1], bold = TRUE), sw(sum_lab[2]))
  w_num <- max(sw(c(st$tp, st$fp, st$fn, st$tn)), sw("TN", bold = TRUE)) + 0.07
  w_txt <- max(sw(L$se[1:k]), sw(L$sp[1:k]), sw(L$se[k + 1], bold = TRUE), sw(L$sp[k + 1], bold = TRUE), sw(L$se[k + 2]), sw("% (95% CI)")) + 0.08
  x_study <- 0.03
  x_num <- x_study + w_study + 0.08 + w_num * (1:4)          # bordes derechos de TP, FP, FN, TN
  x_txt1 <- x_num[4] + 0.06 + w_txt / 2
  w_reg <- (W - 0.1 - (x_num[4] + 0.06 + 2 * w_txt + 0.04 + 0.08)) / 2
  reg1 <- x_txt1 + w_txt / 2 + c(0.02, 0.02 + w_reg)
  x_txt2 <- reg1[2] + 0.08 + w_txt / 2
  reg2 <- x_txt2 + w_txt / 2 + c(0.02, 0.02 + w_reg)
  # ---- posteriores por estudio (con shrinkage) y resumen/prediccion
  lp <- add_linpred_draws(to_long(x), fit, transform = TRUE) %>% ungroup() %>% select(study, measure, value = .linpred) %>%
    left_join(st %>% select(study, row), by = "study")
  dr <- s$draws
  slab_dat <- function(m, reg) {
    meas <- if (m == "se") "sens" else "spec"
    bind_rows(lp %>% filter(measure == meas) %>% transmute(row, value, fill = FILL_STUDY(color), outline = NA),
              tibble(row = 2, value = if (m == "se") dr$se else dr$sp, fill = FILL_SUMM(color), outline = NA),
              tibble(row = 1, value = plogis(if (m == "se") dr$ns else dr$np), fill = FILL_PRED, outline = NA)) %>%
      mutate(x = map_x(value, c(0, 1), reg))
  }
  eye <- function(m, reg) stat_slab(data = slab_dat(m, reg), aes(x = x, y = row, group = row, fill = fill, colour = outline),
                                     side = "both", scale = EYE, normalize = "groups", linewidth = 0.15,
                                     density = density_bounded(bounds = reg), inherit.aes = FALSE)
  marks <- function(m, reg) {
    o <- obs(m) %>% mutate(across(c(est, lo, hi), ~ map_x(.x, c(0, 1), reg)))
    q <- qs[[m]]
    sm <- tibble(row = 2, est = q$s[1], lo = q$s[2], hi = q$s[3]) %>% mutate(across(c(est, lo, hi), ~ map_x(.x, c(0, 1), reg)))
    pr <- tibble(row = 1, est = q$p[1], lo = q$p[2], hi = q$p[3]) %>% mutate(across(c(est, lo, hi), ~ map_x(.x, c(0, 1), reg)))
    list(annotate("segment", x = sm$est, xend = sm$est, y = 2.5, yend = top + 0.5, colour = INK2, linewidth = 0.25, linetype = "22"),
         eye(m, reg), int_layer(o, reg, 0.32), obs_pts(o), int_layer(pr, reg, 0.4, "22"), int_layer(sm, reg, 0.55), summ_pts(sm))
  }
  keys <- tibble(type = c("study", "obs", "summary", "pred"),
                 label = c("Study-level posterior", "Observed (exact 95% CI)", "Summary (median, 95% CrI)", "New study (95% PI)"))
  y_head <- top + 0.78; y_span <- top + 1.22
  p <- ggplot() + zebra_layer(rows, W) +
    axis_layer(reg1, c(0, 1), seq(0, 1, 0.25), c("0", "", "50", "", "100"), top) +
    axis_layer(reg2, c(0, 1), seq(0, 1, 0.25), c("0", "", "50", "", "100"), top) +
    rules_layer(rows, W) + marks("se", reg1) + marks("sp", reg2) +
    # tabla
    txt(x_study, st$row + 0.13, st$study) + txt(x_study, st$row - 0.18, st$sub, pt = PT_SUB, colour = INK2) +
    txt(x_study, c(2, 1), sum_lab, face = c("bold", "plain")) +
    txt(rep(x_num, each = k), rep(st$row, 4), c(st$tp, st$fp, st$fn, st$tn), hjust = 1) +
    txt(x_txt1, rows$row, L$se, hjust = 0.5, face = c(rep("plain", k), "bold", "plain")) +
    txt(x_txt2, rows$row, L$sp, hjust = 0.5, face = c(rep("plain", k), "bold", "plain")) +
    # encabezados
    txt(c(x_study, x_num), y_head, c("Study (index test \u00b7 acquisition)", "TP", "FP", "FN", "TN"), hjust = c(0, 1, 1, 1, 1), face = "bold") +
    txt(c(x_txt1, x_txt2), y_head, "% (95% CI)", hjust = 0.5, colour = INK2) +
    txt(c(mean(c(x_txt1 - w_txt / 2, reg1[2])), mean(c(x_txt2 - w_txt / 2, reg2[2]))), y_span, c("Sensitivity", "Specificity"),
        hjust = 0.5, face = "bold") +
    annotate("segment", x = c(x_txt1 - w_txt / 2, x_txt2 - w_txt / 2) + 0.02, xend = c(reg1[2], reg2[2]), y = y_span - 0.22,
             yend = y_span - 0.22, colour = INK2, linewidth = 0.3) +
    key_layer(keys, color, W, -0.62) + canvas(W, c(-0.95, top + 1.45))
  save_gg(p, name, W, (top + 2.4) * ROW)
  invisible(p)
}
dta_forest(prim, fits[["primary main"]], sP, "Figure_03_forest_primary")
dta_forest(s2, fits[["s2 main"]], sS2, "Figure_S1_forest_S2")

# ------------------------------------------------------------------ figura 3: SROC ----------------------------------
ell <- function(mx, my, S, level = 0.95, npts = 200) {
  r <- sqrt(qchisq(level, 2)); t <- seq(0, 2 * pi, length.out = npts)
  e <- t(chol(S)) %*% rbind(cos(t), sin(t)) * r
  tibble(lx = mx + e[1, ], ly = my + e[2, ])
}
region <- function(dr, pred = FALSE) {
  # espacio logit: x = logit(FPR) = -logit(Sp), y = logit(Se)
  if (pred) m <- cbind(-dr$np, dr$ns) else m <- cbind(-dr$mp, dr$ms)
  e <- ell(mean(m[, 1]), mean(m[, 2]), cov(m))
  tibble(fpr = plogis(e$lx), sens = plogis(e$ly))
}
sroc_points <- function(x) x %>% mutate(fpr = fp / (fp + tn), sens = tp / (tp + fn))
pA_dat <- sroc_points(prim)
cr <- region(sP$draws); pr <- region(sP$draws, TRUE)
regs <- bind_rows(cr %>% mutate(region = "95% credible region"), pr %>% mutate(region = "95% prediction region"))
# Etiquetas fijas: tres estudios caen en (0, 1) y ggrepel no logra separarlos
LAB_POS <- tribble(~study, ~lx, ~ly,
                   "Hur et al. 2012", 0.25, 0.895, "Li et al. 2022", 0.25, 0.835, "Yu et al. 2024", 0.25, 0.775,
                   "Rodríguez-Granillo et al. 2021", 0.07, 0.47, "Schlett et al. 2020", 0.47, 0.955)
labA <- pA_dat %>% select(study, fpr, sens) %>% inner_join(LAB_POS, by = "study")
pA <- ggplot() +
  geom_abline(slope = 1, intercept = 0, colour = GRID, linewidth = 0.3) +
  geom_path(data = regs, aes(fpr, sens, linetype = region), colour = COL["spectral"], linewidth = 0.45) +
  geom_point(data = pA_dat, aes(fpr, sens, size = n), shape = 21, fill = COL["spectral"], colour = "white", stroke = 0.5, alpha = 0.85) +
  geom_point(aes(x = 1 - sP$spec[1], y = sP$sens[1]), shape = 23, size = 2.4, fill = INK, colour = "white", stroke = 0.5) +
  geom_segment(data = labA, aes(x = fpr, y = sens, xend = lx - 0.01, yend = ly), colour = INK2, linewidth = 0.2) +
  geom_text(data = labA, aes(lx, ly, label = study), family = FONT, size = 1.8, colour = INK2, hjust = 0) +
  scale_linetype_manual(values = c(`95% credible region` = "solid", `95% prediction region` = "22"), name = NULL) +
  scale_size_area(max_size = 4.5, breaks = c(25, 100, 300), name = "Patients") +
  scale_x_continuous("1 − Specificity", limits = c(0, 1), breaks = seq(0, 1, .25), expand = c(0.01, 0.01)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, .25), expand = c(0.01, 0.01)) +
  guides(linetype = guide_legend(order = 1, override.aes = list(linewidth = 0.45)), size = guide_legend(order = 2)) +
  coord_equal() + theme_ajr() + theme(legend.position = c(0.72, 0.24), legend.key.size = unit(0.22, "cm"),
                                      legend.spacing.y = unit(0.02, "cm"), legend.box.just = "left")
# B: espectral vs convencional, mismos pacientes
dc_reg <- function(sens_col, spec_col) {
  m <- cbind(-dc[[spec_col]], dc[[sens_col]]); e <- ell(mean(m[, 1]), mean(m[, 2]), cov(m))
  tibble(fpr = plogis(e$lx), sens = plogis(e$ly))
}
reg_s <- dc_reg("b_measuresens:testspectral", "b_measurespec:testspectral") %>% mutate(test = "Spectral CT")
reg_c <- dc_reg("b_measuresens:testconventional", "b_measurespec:testconventional") %>% mutate(test = "Conventional CT")
pts_s1 <- bind_rows(sroc_points(s1s) %>% mutate(test = "Spectral CT"), sroc_points(s1c) %>% mutate(test = "Conventional CT"))
seg <- pts_s1 %>% select(study, test, fpr, sens) %>% pivot_wider(names_from = test, values_from = c(fpr, sens))
sum_pts <- tibble(test = c("Spectral CT", "Conventional CT"), fpr = c(1 - median(sp_s), 1 - median(sp_c)), sens = c(median(se_s), median(se_c)))
tcol <- c(`Spectral CT` = unname(COL["spectral"]), `Conventional CT` = unname(COL["conventional"]))
pB <- ggplot() +
  geom_abline(slope = 1, intercept = 0, colour = GRID, linewidth = 0.3) +
  geom_segment(data = seg, aes(x = `fpr_Conventional CT`, y = `sens_Conventional CT`, xend = `fpr_Spectral CT`, yend = `sens_Spectral CT`),
               colour = INK2, linewidth = 0.25, alpha = 0.7) +
  geom_path(data = bind_rows(reg_s, reg_c), aes(fpr, sens, colour = test), linewidth = 0.5) +
  geom_point(data = pts_s1, aes(fpr, sens, fill = test, size = n), shape = 21, colour = "white", stroke = 0.5, alpha = 0.9) +
  geom_point(data = sum_pts, aes(fpr, sens, fill = test), shape = 23, size = 2.4, colour = INK, stroke = 0.4, show.legend = FALSE) +
  annotate("text", x = c(0.09, 0.27), y = c(0.71, 0.62), label = c("Spectral CT", "Conventional CT"), family = FONT, size = 1.9,
           colour = INK, hjust = 0) +
  scale_colour_manual(values = tcol, guide = "none") + scale_fill_manual(values = tcol, name = NULL) +
  scale_size_area(max_size = 4.5, guide = "none") +
  scale_x_continuous("1 − Specificity", limits = c(0, 1), breaks = seq(0, 1, .25), expand = c(0.01, 0.01)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, .25), expand = c(0.01, 0.01)) +
  guides(fill = guide_legend(override.aes = list(size = 2.4, shape = 21, colour = "white"))) +
  coord_equal() + theme_ajr() + theme(legend.position = c(0.72, 0.2), legend.key.size = unit(0.22, "cm"))
save_gg(pA + pB + plot_annotation(tag_levels = FIG_TAG), "Figure_04_sroc", 5, 2.7)

# ------------------------------------------------------------------ figura 2: QUADAS-3 ------------------------------
# una fila por reporte: la estimacion que entra en su sintesis principal (Zhang: defectos de llenado, E3), en orden cronologico
MAIN_EST <- c(Zhang_2024 = "E3")
chron <- unique(q3$study); chron <- chron[order(as.integer(sub(".* ", "", chron)), chron)]
qe1 <- q3 %>% filter(estimate == ifelse(paper_id %in% names(MAIN_EST), MAIN_EST[paper_id], "E1")) %>%
  mutate(study = factor(study, levels = rev(chron)))
dom <- c(rob_participants = "Participants", rob_index = "Index test", rob_target = "Target condition", rob_analysis = "Analysis",
         overall_risk = "Overall", app_participants = "Participants ", app_index = "Index test ", app_target = "Target condition ",
         overall_applicability = "Overall ")
tl <- qe1 %>% select(study, all_of(names(dom))) %>% pivot_longer(-study, names_to = "domain", values_to = "j") %>%
  mutate(j = recode(tolower(j), "insufficient information" = "insufficient"), lab = recode(j, low = "+", high = "−", insufficient = "?"),
         domain = factor(dom[domain], levels = unname(dom)), block = ifelse(grepl(" $", as.character(domain)), "Applicability concerns", "Risk of bias"))
p4a <- ggplot(tl, aes(domain, study)) +
  geom_point(aes(fill = j), shape = 21, size = 3.6, colour = "white", stroke = 0.6) +
  geom_text(aes(label = lab), colour = "white", size = 2.3, fontface = "bold", family = FONT) +
  facet_grid(. ~ factor(block, levels = c("Risk of bias", "Applicability concerns")), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = STATUS, breaks = intersect(c("low", "insufficient", "high"), unique(tl$j)),
                    labels = c(low = "Low (+)", insufficient = "Insufficient information (?)", high = "High (−)")[intersect(c("low", "insufficient", "high"), unique(tl$j))],
                    name = NULL) +
  scale_x_discrete(position = "top", labels = function(x) trimws(x)) + labs(x = NULL, y = NULL) + theme_ajr(6.5) +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle = 35, hjust = 0, vjust = 0, colour = INK), legend.position = "bottom")
sb <- tl %>% filter(!grepl("Overall", domain)) %>% count(block, domain, j) %>% group_by(block, domain) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
  mutate(j = factor(j, levels = c("high", "insufficient", "low")))
bar_block <- function(b, xlab) {
  ggplot(sb %>% filter(block == b) %>% droplevels(), aes(p, domain, fill = j)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(p >= 0.1, paste0(round(100 * p), "%"), "")), position = position_stack(vjust = 0.5), colour = "white",
              size = 1.8, family = FONT, fontface = "bold") +
    scale_fill_manual(values = STATUS, guide = "none") +
    scale_x_continuous(xlab, labels = if (is.null(xlab)) NULL else scales::percent, expand = c(0, 0)) +
    scale_y_discrete(labels = function(x) trimws(x), limits = rev) + labs(y = NULL, title = b) + theme_ajr(6.5) +
    theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold", size = 6.5, margin = margin(b = 2)),
          plot.title.position = "plot")
}
p4b <- bar_block("Risk of bias", NULL) / bar_block("Applicability concerns", "Proportion of reports") + plot_layout(heights = c(4, 3))
save_gg(wrap_elements(p4a) / wrap_elements(p4b) + plot_layout(heights = c(2.0, 1.75)) + plot_annotation(tag_levels = FIG_TAG) &
          theme(plot.tag = element_text(face = "bold", size = 9, family = FONT, colour = INK)),
        "Figure_02_quadas3", 5, 5)

# ------------------------------------------------------------------ figura 5: SMD (S3) ------------------------------
smd_forest <- function(name, W = 5, ROW = 0.262, XL = c(-9, 3), color = COL["spectral"]) {
  blk <- function(o, grp, b) {
    e <- as.data.frame(o$es)
    e <- e[order(as.integer(sub(".* ", "", e$study)), e$study), ]  # orden cronologico, como en las tablas
    unit <- ifelse(e$unit == "ratio", "ratio", e$unit)
    bind_rows(tibble(kind = "header", block = b, label = grp, sub = "", nn = "", wt = "", est = NA, lo = NA, hi = NA, sdv = NA, df = NA),
              tibble(kind = "study", block = b, label = e$study,
                     sub = paste0(e$measure, ", ", unit, ifelse(!is.na(e$conversion) & e$conversion != "", "*", "")),
                     nn = paste0(e$n1, "/", e$n2), wt = sprintf("%.1f", e$w), est = e$yi, lo = e$lo, hi = e$hi, sdv = sqrt(e$vi), df = NA),
              tibble(kind = c("summary", "pi"), block = b, label = c("Summary (REML, Hartung-Knapp)", "95% prediction interval"), sub = "",
                     nn = "", wt = c("100", ""), est = o$p$pred, lo = c(o$p$ci.lb, o$p$pi.lb), hi = c(o$p$ci.ub, o$p$pi.ub),
                     sdv = c(o$p$se, sqrt(o$p$se^2 + o$m$tau2)), df = o$m$ddf))
  }
  sm <- bind_rows(blk(s3a, "Thrombus vs stasis", 1), blk(s3b, "Thrombus vs normal LAA", 2)) %>%
    mutate(row = rev(seq_len(n())), face = ifelse(kind %in% c("header", "summary"), "bold", "plain"),
           lab = minus(case_when(kind == "header" ~ "", kind == "pi" ~ sprintf("%.2f to %.2f", lo, hi),
                                 TRUE ~ sprintf("%.2f (%.2f to %.2f)", est, lo, hi))))
  top <- max(sm$row)
  stt <- sm %>% filter(kind == "study"); oth <- sm %>% filter(kind != "study")
  w_study <- max(sw(stt$label), sw(stt$sub, PT_SUB), sw(oth$label[oth$kind != "pi"], bold = TRUE), sw("95% prediction interval"))
  w_nn <- max(sw(sm$nn), sw("n (T/C)", bold = TRUE)) + 0.08
  w_wt <- max(sw(sm$wt), sw("Weight, %", bold = TRUE)) + 0.08
  w_txt <- max(sw(sm$lab[sm$kind %in% c("study", "pi")]), sw(sm$lab[sm$kind == "summary"], bold = TRUE)) + 0.1
  x_study <- 0.03; x_nn <- x_study + w_study + 0.06 + w_nn; x_wt <- x_nn + w_wt
  x_txt <- x_wt + 0.06 + w_txt / 2
  reg <- c(x_txt + w_txt / 2 + 0.06, W - 0.12)
  mx <- function(v) map_x(v, XL, reg); ms <- diff(reg) / diff(XL)
  sl <- bind_rows(sm %>% filter(kind == "study") %>% mutate(dist = dist_normal(mx(est), sdv * ms)),
                  sm %>% filter(kind %in% c("summary", "pi")) %>% mutate(dist = dist_student_t(df, mx(est), sdv * ms))) %>%
    mutate(fill = case_when(kind == "study" ~ FILL_STUDY(color), kind == "summary" ~ FILL_SUMM(color), TRUE ~ FILL_PRED),
           outline = NA)
  mk <- function(k) sm %>% filter(kind == k) %>% transmute(row, est = mx(est), lo = mx(lo), hi = mx(hi))
  keys <- tibble(type = c("study", "obs", "summary", "pred"),
                 label = c("Sampling distribution", "Observed SMD (95% CI)", "Summary (t, k \u2212 1 df; 95% CI)", "New study (95% PI)"))
  y_head <- top + 0.92
  p <- ggplot() + zebra_layer(sm, W) +
    axis_layer(reg, XL, seq(-8, 2, 2), minus(as.character(seq(-8, 2, 2))), top, title = "SMD (Hedges g), thrombus \u2212 comparator") +
    annotate("segment", x = mx(0), xend = mx(0), y = 0.5, yend = top + 0.5, colour = INK2, linewidth = 0.3) +
    rules_layer(sm, W) +
    stat_slab(data = sl, aes(y = row, xdist = dist, group = row, fill = fill, colour = outline), side = "both", scale = EYE,
              normalize = "groups", linewidth = 0.15, limits = reg, p_limits = c(0.001, 0.999), inherit.aes = FALSE) +
    int_layer(mk("study"), reg, 0.32) + obs_pts(mk("study")) + int_layer(mk("pi"), reg, 0.4, "22") + int_layer(mk("summary"), reg, 0.55) +
    summ_pts(mk("summary")) +
    txt(x_study, stt$row + 0.13, stt$label) + txt(x_study, stt$row - 0.18, stt$sub, pt = PT_SUB, colour = INK2) +
    txt(x_study, oth$row, oth$label, face = oth$face) +
    txt(x_nn, sm$row, sm$nn, hjust = 1) + txt(x_wt, sm$row, sm$wt, hjust = 1, face = sm$face) +
    txt(x_txt, sm$row, sm$lab, hjust = 0.5, face = sm$face) +
    txt(c(x_study, x_nn, x_wt, x_txt), y_head, c("Comparison and study", "n (T/C)", "Weight, %", "SMD (95% CI)"),
        hjust = c(0, 1, 1, 0.5), face = "bold") +
    txt(reg[1] + 0.02, -0.55, "\u2190 Lower iodine in thrombus", pt = PT_SUB, colour = INK2) +
    key_layer(keys, color, W, -1.25) + canvas(W, c(-1.6, top + 1.3))
  save_gg(p, name, W, (top + 2.9) * ROW)
}
smd_forest("Figure_05_smd_iodine")

# ------------------------------------------------------------------ figura 6: probabilidad post-test ----------------
# A: todos los pacientes (modelo primario). B: riesgo residual tras una lectura espectral negativa dentro de los defectos
# de llenado (modelo S2), con lo observado en Zhang 2024 (unica cohorte temprana de defectos) y la prevalencia de trombo
# en anticoagulados que hoy van a TEE (2.73%, Lurie 2021) como referencia.
post_curve <- function(dr, pre) bind_rows(lapply(pre, function(pt) {
  odds <- pt / (1 - pt); lrp <- dr$se / (1 - dr$sp); lrn <- (1 - dr$se) / dr$sp
  pp <- odds * lrp / (1 + odds * lrp); pn <- odds * lrn / (1 + odds * lrn)
  tibble(pre = pt, result = c("Positive", "Negative"),
         med = c(median(pp), median(pn)), lo = c(quantile(pp, .025), quantile(pn, .025)), hi = c(quantile(pp, .975), quantile(pn, .975)))
}))
post <- post_curve(sP$draws, seq(0.005, 0.5, by = 0.005))
p6a <- ggplot(post, aes(pre, med, group = result)) +
  annotate("rect", xmin = prevs[1], xmax = prevs[2], ymin = -Inf, ymax = Inf, fill = GRID, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.25, linetype = "dotted") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = COL["spectral"], alpha = 0.18) +
  geom_line(colour = COL["spectral"], linewidth = 0.6) +
  annotate("text", x = c(0.30, 0.30), y = c(0.70, 0.13), label = c("Positive read", "Negative read"), family = FONT, size = 1.95,
           colour = INK, hjust = 0) +
  annotate("text", x = mean(prevs), y = 0.99, label = "Prevalence range,
single-gate studies", family = FONT, size = 1.75,
           colour = INK2, vjust = 1, lineheight = 0.9) +
  scale_x_continuous("Pretest probability, all patients", labels = scales::percent, limits = c(0, 0.5), expand = c(0, 0)) +
  scale_y_continuous("Posttest probability", labels = scales::percent, limits = c(0, 1), expand = c(0, 0.01)) +
  labs(title = "All patients (primary model)") + theme_ajr(6.5) +
  theme(plot.title = element_text(size = 6.5, face = "bold", colour = INK, margin = margin(b = 3)), plot.title.position = "plot")
post_b <- post_curve(sS2$draws, seq(0.05, 0.6, by = 0.005)) %>% filter(result == "Negative")
zb <- binom.test(5, 43)
zh <- tibble(pre = 18 / 66, est = 5 / 43, lo = zb$conf.int[1], hi = zb$conf.int[2])
p6b <- ggplot(post_b, aes(pre, med)) +
  geom_hline(yintercept = 0.0273, colour = INK2, linewidth = 0.3, linetype = "22") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = COL["spectral"], alpha = 0.18) +
  geom_line(colour = COL["spectral"], linewidth = 0.6) +
  geom_errorbar(data = zh, aes(x = pre, ymin = lo, ymax = hi), inherit.aes = FALSE, width = 0.014, colour = INK, linewidth = 0.35) +
  geom_point(data = zh, aes(pre, est), inherit.aes = FALSE, shape = 21, fill = "white", colour = INK, size = 1.6, stroke = 0.45) +
  annotate("text", x = zh$pre + 0.025, y = zh$est + 0.01, label = "Zhang et al. 2024,
observed (5 of 43)", hjust = 0, size = 1.8,
           family = FONT, colour = INK, lineheight = 0.9) +
  annotate("text", x = 0.06, y = 0.0273 + 0.004, label = "2.7%: anticoagulated\npatients sent to TEE", hjust = 0, vjust = 0,
           size = 1.7, family = FONT, colour = INK2, lineheight = 0.9) +
  annotate("text", x = 0.32, y = 0.265, label = "Filling-defect model,
median and 95% CrI", hjust = 0, size = 1.8, family = FONT,
           colour = INK2, lineheight = 0.9) +
  scale_x_continuous("Pretest probability, filling defects", labels = scales::percent, limits = c(0.05, 0.6), expand = c(0, 0)) +
  scale_y_continuous("Posttest probability after a negative read", labels = scales::percent, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 0.3)) +
  labs(title = "Residual risk after a negative read") + theme_ajr(6.5) +
  theme(plot.title = element_text(size = 6.5, face = "bold", colour = INK, margin = margin(b = 3)), plot.title.position = "plot")
save_gg(p6a + p6b + plot_annotation(tag_levels = FIG_TAG) & theme(plot.tag = element_text(face = "bold", size = 9, family = FONT, colour = INK)),
        "Figure_06_posttest", 5, 2.8)

# ------------------------------------------------------------------ figura 1: PRISMA 2020 ---------------------------
# Paquete PRISMA2020 (Haddaway et al., Campbell Syst Rev 2022) con sus colores por defecto: encabezado Goldenrod1 y rotulos
# laterales LightSteelBlue2. Solo bases de datos; la caja de reportes excluidos lleva los motivos.
# El SVG se exporta con DiagrammeRsvg (PRISMA_save pide pandoc) y los rotulos laterales se escriben aqui, como hace el paquete.
tpl <- read.csv(system.file("extdata", "PRISMA.csv", package = "PRISMA2020"), stringsAsFactors = FALSE)
tpl$n <- NA
setn <- function(key, value) tpl$n[tpl$data %in% key] <<- as.character(value)
setx <- function(key, text, col = "data") tpl$boxtext[tpl[[col]] %in% key] <<- text
setn("database_results", prisma$identified)
setn("database_specific_results", sprintf("PubMed, %d; Embase, %d; Web of Science, %d", prisma$db$PubMed, prisma$db$Embase, prisma$db$WoS))
setn("duplicates", prisma$duplicates)
setn("records_screened", prisma$screened); setn("records_excluded", prisma$excluded_tiab)
setn("dbr_sought_reports", prisma$sought); setn("dbr_notretrieved_reports", prisma$not_retrieved)
setn("dbr_assessed", prisma$assessed)
setn("dbr_excluded", "No accuracy analysis against a reference, 3; Accuracy only for findings other than LAA thrombus, 1")
setn("new_studies", prisma$included_studies); setn("new_reports", prisma$included_reports)
setx("newstud", "Identification of studies via databases", col = "box")
setx("duplicates", "Duplicate records removed")
setx("dbr_excluded", sprintf("Reports excluded (n = %d)", prisma$excluded_ft))
setx("new_studies", "Studies included in review"); setx("new_reports", "Reports of included studies")
fd <- PRISMA2020::PRISMA_flowdiagram(PRISMA2020::PRISMA_data(tpl), previous = FALSE, other = FALSE, detail_databases = TRUE,
                                     fontsize = 12, font = FONT)
px <- xml2::read_xml(DiagrammeRsvg::export_svg(fd)); ns <- xml2::xml_ns(px)
side <- c(identification = "Identification", screening = "Screening", included = "Included")
for (k in names(side)) {
  tx <- xml2::xml_find_first(px, sprintf("//d1:g[d1:title = '%s']//d1:text", k), ns)
  x0 <- as.numeric(xml2::xml_attr(tx, "x")); y0 <- as.numeric(xml2::xml_attr(tx, "y"))
  xml2::xml_attr(tx, "x") <- -y0; xml2::xml_attr(tx, "y") <- x0 + 4
  xml2::xml_attr(tx, "transform") <- "rotate(-90)"; xml2::xml_attr(tx, "font-weight") <- "bold"
  xml2::xml_text(tx) <- side[[k]]
}
hd <- xml2::xml_find_first(px, "//d1:g[d1:title = '3']//d1:text", ns); xml2::xml_attr(hd, "font-weight") <- "bold"
svg_f <- tempfile(fileext = ".svg"); xml2::write_xml(px, svg_f)
dm <- dim(rsvg::rsvg(svg_f)); asp <- dm[1] / dm[2]                  # alto / ancho
w_in <- min(5, 5 / asp); h_in <- w_in * asp                          # AJR: 5 pulgadas como maximo por lado
rsvg::rsvg_pdf(svg_f, file.path(fig_dir, "Figure_01_prisma.pdf"), width = w_in * 72, height = h_in * 72)
img <- magick::image_read(rsvg::rsvg_png(svg_f, width = round(w_in * FIG_DPI), height = round(h_in * FIG_DPI)))
img <- magick::image_convert(magick::image_flatten(magick::image_background(img, "white")), type = "TrueColor", matte = FALSE)
magick::image_write(img, file.path(fig_dir, "Figure_01_prisma.png"), format = "png", density = "600x600")
magick::image_write(magick::image_resize(img, sprintf("%dx%d!", round(w_in * 300), round(h_in * 300))),
                    file.path(fig_dir, "Figure_01_prisma_rgb300.tif"), format = "tiff", density = "300x300", compression = "LZW")

writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo.txt"))
message("DONE")
print(pooled %>% select(analysis, method, k, sens, spec, sens_pi, spec_pi), n = 50, width = 250)
print(diagnostics)
print(s1_tab); print(s3_tab); print(cons)
