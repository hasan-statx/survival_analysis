# =============================================================================
# 02_survival_analysis.R
# Advanced Survival Analysis — Neonatal Mortality
# Methods: KM, Log-rank, Cox PH, Schoenfeld residuals,
#          Competing risks (cmprsk/tidycmprsk), Frailty models
#
# Author : Hasan Mahmud Sujan
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(broom)
  library(gtsummary)
  library(gt)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(ggtext)
  library(cmprsk)        # competing risks
  library(kableExtra)
  library(pammtools)     # geom_stepribbon
})

# ── Theme & Palette ────────────────────────────────────────────────────────────
BLUE_PALETTE <- c("#1F4E79","#2E75B6","#70AD47","#ED7D31","#FFC000")
JiVitA_COLORS <- c(
  "Control"     = "#C00000",
  "BEP_supplement" = "#2E75B6",
  "MNP"         = "#70AD47",
  "BEP_plus_MNP"= "#7030A0"
)
ANAEMIA_COLORS <- c("None"="#2E75B6","Mild"="#70AD47",
                    "Moderate"="#ED7D31","Severe"="#C00000")

pub_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 1,
                                   colour = "#1F4E79"),
      plot.subtitle = element_text(colour = "grey40", size = base_size - 1),
      plot.caption  = element_text(colour = "grey50", size = 8, hjust = 0),
      axis.title    = element_text(face = "bold", colour = "#1F4E79"),
      legend.title  = element_text(face = "bold"),
      legend.background = element_rect(fill = "white", colour = "grey85"),
      strip.background = element_rect(fill = "#D6E4F0", colour = NA),
      strip.text    = element_text(face = "bold", colour = "#1F4E79"),
      panel.grid.major.y = element_line(colour = "grey92"),
    )
}
theme_set(pub_theme())

# ── Load data ─────────────────────────────────────────────────────────────────
cat("=================================================================\n")
cat("  ADVANCED SURVIVAL ANALYSIS — NEONATAL MORTALITY COHORT\n")
cat("=================================================================\n\n")

df <- read_csv("data/cohort_baseline.csv", show_col_types = FALSE) %>%
  mutate(
    intervention = factor(intervention,
                          levels = c("Control","BEP_supplement","MNP","BEP_plus_MNP"),
                          labels = c("Control","BEP Supplement","MNP","BEP + MNP")),
    anaemia_enrol = factor(anaemia_enrol,
                           levels = c("None","Mild","Moderate","Severe")),
    wealth_cat    = factor(wealth_index,
                           labels = c("Poorest","Poor","Middle","Wealthy","Wealthiest")),
    lbw_factor    = factor(low_birthweight, labels = c("Normal","LBW")),
    preterm_factor= factor(preterm, labels = c("Term","Preterm")),
    age_cat       = cut(age_enrol, breaks = c(14,19,24,29,35,46),
                        labels = c("<20","20–24","25–29","30–34","35+")),
    hb_cat        = cut(haemoglobin_enrol,
                        breaks = c(0, 7, 10, 11, 20),
                        labels = c("Severe (<7)","Moderate (7–10)",
                                   "Mild (10–11)","None (≥11)")),
    # Surv object
    t = surv_time_days,
    e = neonatal_death
  )

cat(sprintf("[Data] N = %d | Events = %d (%.1f per 1,000)\n\n",
            nrow(df), sum(df$e), mean(df$e)*1000))

dir.create("outputs", showWarnings = FALSE)

# =============================================================================
# SECTION 1 — KAPLAN-MEIER CURVES
# =============================================================================
cat("── Section 1: Kaplan-Meier Estimates ──\n")

fit_km_overall <- survfit(Surv(t, e) ~ 1, data = df)
fit_km_int     <- survfit(Surv(t, e) ~ intervention, data = df)
fit_km_anaemia <- survfit(Surv(t, e) ~ anaemia_enrol, data = df)
fit_km_lbw     <- survfit(Surv(t, e) ~ lbw_factor, data = df)

# ── 1a: KM by intervention ────────────────────────────────────────────────────
p_km_int <- ggsurvplot(
  fit_km_int,
  data          = df,
  fun           = "event",       # cumulative incidence
  palette       = unname(JiVitA_COLORS),
  conf.int      = TRUE,
  conf.int.alpha= 0.12,
  risk.table    = TRUE,
  risk.table.height = 0.28,
  risk.table.title  = "Number at risk",
  risk.table.fontsize = 3.5,
  pval          = TRUE,
  pval.method   = TRUE,
  log.rank.weights = "S1",       # Peto-Peto weights
  legend.title  = "Intervention",
  legend.labs   = levels(df$intervention),
  xlab          = "Days after birth",
  ylab          = "Cumulative neonatal mortality",
  title         = "Cumulative Neonatal Mortality by Intervention Arm",
  subtitle      = "Cluster-RCT cohort, Bangladesh 2019–2023 (N = 5,000)",
  caption       = "Shaded bands: 95% CI. P-value: Peto-Peto weighted log-rank test.",
  surv.scale    = "percent",
  ggtheme       = pub_theme(),
  break.time.by = 7,
  xlim          = c(0, 28)
)

# ── 1b: KM by anaemia status ──────────────────────────────────────────────────
p_km_anaemia <- ggsurvplot(
  fit_km_anaemia,
  data          = df,
  fun           = "event",
  palette       = unname(ANAEMIA_COLORS),
  conf.int      = TRUE,
  conf.int.alpha= 0.12,
  risk.table    = TRUE,
  risk.table.height = 0.28,
  pval          = TRUE,
  pval.method   = TRUE,
  legend.title  = "Maternal Anaemia\n(Hb at enrolment)",
  xlab          = "Days after birth",
  ylab          = "Cumulative neonatal mortality",
  title         = "Neonatal Mortality by Maternal Anaemia Status at Enrolment",
  ggtheme       = pub_theme(),
  break.time.by = 7,
  xlim          = c(0, 28)
)

# Save KM plots
ggsave("outputs/fig1a_km_intervention.png",
       print(p_km_int)$plot, width = 9, height = 7, dpi = 300)
ggsave("outputs/fig1b_km_anaemia.png",
       print(p_km_anaemia)$plot, width = 9, height = 7, dpi = 300)
cat("  → KM plots saved\n")

# =============================================================================
# SECTION 2 — LOG-RANK & WEIGHTED TESTS
# =============================================================================
cat("\n── Section 2: Log-rank Tests ──\n")

lr_int     <- survdiff(Surv(t, e) ~ intervention,   data = df)
lr_anaemia <- survdiff(Surv(t, e) ~ anaemia_enrol,  data = df)
lr_lbw     <- survdiff(Surv(t, e) ~ lbw_factor,     data = df)

p_chi <- function(x) round(pchisq(x$chisq, df = length(x$n)-1, lower.tail=FALSE), 4)
cat(sprintf("  Log-rank p: Intervention=%.4f | Anaemia=%.4f | LBW=%.4f\n",
            p_chi(lr_int), p_chi(lr_anaemia), p_chi(lr_lbw)))

# =============================================================================
# SECTION 3 — COX PROPORTIONAL HAZARDS REGRESSION
# =============================================================================
cat("\n── Section 3: Cox PH Regression ──\n")

# 3a: Unadjusted
cox_unadj <- coxph(
  Surv(t, e) ~ intervention,
  data = df, ties = "efron"
)

# 3b: Fully adjusted
cox_adj <- coxph(
  Surv(t, e) ~ intervention
    + low_birthweight + preterm
    + anaemia_enrol
    + age_enrol + parity
    + edu_years
    + wealth_index
    + food_insecurity
    + skilled_birth_attendant
    + anc_visits
    + child_sex
    + season,
  data = df, ties = "efron"
)

# 3c: Frailty model (cluster random effect)
cox_frailty <- coxph(
  Surv(t, e) ~ intervention
    + low_birthweight + preterm
    + anaemia_enrol
    + age_enrol + parity
    + edu_years + wealth_index
    + food_insecurity
    + skilled_birth_attendant
    + child_sex
    + frailty(cluster_id, distribution = "gamma"),
  data = df, ties = "efron"
)

print(summary(cox_adj))

# ── Forest plot of adjusted HR ────────────────────────────────────────────────
tidy_adj <- tidy(cox_adj, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(!str_detect(term, "Intercept")) %>%
  mutate(
    term_clean = case_when(
      term == "interventionBEP Supplement" ~ "BEP Supplement vs Control",
      term == "interventionMNP"            ~ "MNP vs Control",
      term == "interventionBEP + MNP"      ~ "BEP + MNP vs Control",
      term == "low_birthweightTRUE"        ~ "Low birthweight",
      term == "pretermTRUE"               ~ "Preterm birth",
      term == "anaemia_enrolMild"         ~ "Anaemia: Mild vs None",
      term == "anaemia_enrolModerate"     ~ "Anaemia: Moderate vs None",
      term == "anaemia_enrolSevere"       ~ "Anaemia: Severe vs None",
      term == "age_enrol"                 ~ "Age at enrolment (years)",
      term == "parity"                    ~ "Parity",
      term == "edu_years"                 ~ "Education (years)",
      term == "wealth_index"              ~ "Wealth index",
      term == "food_insecurity"           ~ "Food insecurity score",
      term == "skilled_birth_attendant"   ~ "Skilled birth attendant",
      term == "anc_visits"                ~ "ANC visits",
      term == "child_sexMale"             ~ "Child sex: Male",
      term == "seasonPreMonsoon"          ~ "Season: Pre-monsoon",
      term == "seasonMonsoon"             ~ "Season: Monsoon",
      term == "seasonPostMonsoon"         ~ "Season: Post-monsoon",
      TRUE                               ~ term
    ),
    group = case_when(
      str_detect(term_clean, "BEP|MNP")           ~ "Intervention",
      str_detect(term_clean, "birthweight|Preterm")~ "Birth outcomes",
      str_detect(term_clean, "Anaemia")            ~ "Maternal nutrition",
      str_detect(term_clean, "Skilled|ANC")        ~ "Healthcare access",
      TRUE                                         ~ "Socioeconomic"
    ),
    sig = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",
      TRUE            ~ "p ≥ 0.05"
    ),
    direction = if_else(estimate < 1, "Protective", "Risk")
  )

p_forest <- ggplot(tidy_adj,
  aes(x = estimate, y = fct_reorder(term_clean, estimate),
      colour = group, shape = sig)) +
  geom_vline(xintercept = 1, lty = 2, colour = "grey50", lwd = 0.7) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.3, lwd = 0.7) +
  geom_point(size = 3.5) +
  scale_colour_brewer(palette = "Set1", name = "Variable group") +
  scale_shape_manual(
    values = c("p < 0.001"=16,"p < 0.01"=17,"p < 0.05"=15,"p ≥ 0.05"=1),
    name   = "Significance"
  ) +
  scale_x_log10(breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5),
                labels = c("0.25","0.50","0.75","1.0","1.5","2.0","3.0","5.0")) +
  labs(
    x       = "Adjusted Hazard Ratio (95% CI) — log scale",
    y       = NULL,
    title   = "Forest Plot: Cox PH Regression — Neonatal Mortality",
    subtitle= "Fully adjusted model (N = 5,000). Reference: Control arm, No anaemia, Female child, Winter.",
    caption = "HR < 1: lower hazard (protective). HR > 1: higher hazard. Tie-handling: Efron."
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.y = element_text(size = 9)
  )

ggsave("outputs/fig2_forest_cox.png", p_forest,
       width = 11, height = 10, dpi = 300)
cat("  → Forest plot saved\n")

# ── Publication table: gtsummary ──────────────────────────────────────────────
cox_tbl <- tbl_regression(
  cox_adj,
  exponentiate = TRUE,
  label = list(
    intervention ~ "Intervention arm",
    low_birthweight ~ "Low birthweight (<2,500g)",
    preterm ~ "Preterm birth (<37 wk)",
    anaemia_enrol ~ "Maternal anaemia (Hb at enrolment)",
    age_enrol ~ "Maternal age at enrolment (years)",
    parity ~ "Parity",
    edu_years ~ "Maternal education (years)",
    wealth_index ~ "Household wealth index",
    food_insecurity ~ "Household food insecurity score",
    skilled_birth_attendant ~ "Skilled birth attendant (yes)",
    anc_visits ~ "ANC visits (count)",
    child_sex ~ "Child sex",
    season ~ "Season of enrolment"
  )
) %>%
  bold_p(t = 0.05) %>%
  add_significance_stars() %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("**Table 1. Cox PH regression: adjusted hazard ratios for neonatal mortality**")

# Save as HTML
cox_tbl %>%
  as_gt() %>%
  gt::gtsave("outputs/table1_cox_regression.html")
cat("  → Regression table saved\n")

# =============================================================================
# SECTION 4 — PROPORTIONAL HAZARDS ASSUMPTION TEST
# =============================================================================
cat("\n── Section 4: PH Assumption (Schoenfeld Residuals) ──\n")

ph_test <- cox.zph(cox_adj)
print(ph_test)

p_zph <- ggcoxzph(ph_test, point.size = 1, point.alpha = 0.4,
                  ggtheme = pub_theme(),
                  font.main = c(11, "bold", "#1F4E79"))
# Save first 4 panels as representative
ggsave("outputs/fig3_schoenfeld_residuals.png",
       ggpubr::ggarrange(plotlist = p_zph[1:4], ncol = 2, nrow = 2),
       width = 12, height = 9, dpi = 300)
cat("  → Schoenfeld residual plots saved\n")

# =============================================================================
# SECTION 5 — COMPETING RISKS (Fine-Gray subdistribution HR)
# =============================================================================
cat("\n── Section 5: Competing Risks Analysis ──\n")

# Competing events: neonatal death (1) vs censored/alive (0)
# Treat post-neonatal death in 1st week as competing risk conceptually
df_cr <- df %>%
  mutate(
    cr_status = case_when(
      neonatal_death == 1 ~ 1L,   # event of interest
      TRUE                ~ 0L    # censored
    )
  )

# Cumulative Incidence Function by intervention (cmprsk)
cif_int <- cuminc(
  ftime   = df_cr$t,
  fstatus = df_cr$cr_status,
  group   = df_cr$intervention
)

# Fine-Gray model (intervention, LBW, preterm, anaemia)
fg_covars <- model.matrix(
  ~ intervention + low_birthweight + preterm + anaemia_enrol
    + edu_years + wealth_index + skilled_birth_attendant,
  data = df_cr
)[, -1]   # drop intercept

fg_model <- crr(
  ftime   = df_cr$t,
  fstatus = df_cr$cr_status,
  cov1    = fg_covars
)

fg_summary <- summary(fg_model)
cat("\n  Fine-Gray Model (neonatal mortality):\n")
print(fg_summary$coef[,c("coef","exp(coef)","se(coef)","p-value")] %>%
      as_tibble(rownames="term") %>%
      mutate(across(where(is.numeric), ~round(.,4))), n=20)

# CIF plot using tidycmprsk-style with ggplot2
cif_df <- imap_dfr(cif_int[str_detect(names(cif_int), "^[^\\s]+ 1$")], function(x, nm) {
  arm <- str_remove(nm, " 1$")
  tibble(time = x$time, est = x$est, arm = arm)
})

p_cif <- ggplot(cif_df, aes(x = time, y = est * 100, colour = arm)) +
  geom_step(lwd = 1.1) +
  scale_colour_manual(values = unname(JiVitA_COLORS), name = "Intervention") +
  scale_x_continuous(breaks = seq(0,28,7)) +
  labs(
    x       = "Days after birth",
    y       = "Cumulative Incidence (%)",
    title   = "Cumulative Incidence of Neonatal Mortality by Intervention",
    subtitle= "Competing risks framework (Fine-Gray method)",
    caption = "Competing event: discharged alive at day 28."
  )

ggsave("outputs/fig4_competing_risks_cif.png", p_cif,
       width = 9, height = 6, dpi = 300)
cat("  → Competing risks plot saved\n")

# =============================================================================
# SECTION 6 — RESTRICTED MEAN SURVIVAL TIME (RMST)
# =============================================================================
cat("\n── Section 6: Restricted Mean Survival Time ──\n")

rmst_res <- survRM2::rmst2(
  time   = df$t,
  status = df$e,
  arm    = as.integer(df$intervention == "BEP Supplement"),
  tau    = 28
)
cat("\n  RMST (BEP Supplement vs Control, tau=28 days):\n")
print(rmst_res$unadjusted.result)

# =============================================================================
# SECTION 7 — HAZARD FUNCTION PLOTS (smoothed)
# =============================================================================
cat("\n── Section 7: Smoothed Hazard Functions ──\n")

library(muhaz)

hz_plots <- map(levels(df$intervention), function(arm) {
  sub_df <- filter(df, intervention == arm)
  hz  <- muhaz(sub_df$t, sub_df$e, min.time = 0.1, max.time = 27,
               bw.method = "local", b.cor = "left")
  tibble(time = hz$est.grid, hazard = hz$haz.est, intervention = arm)
}) %>% bind_rows()

p_hazard <- ggplot(hz_plots, aes(x = time, y = hazard * 1000,
                                  colour = intervention)) +
  geom_line(lwd = 1.2) +
  scale_colour_manual(values = unname(JiVitA_COLORS)) +
  scale_x_continuous(breaks = seq(0, 28, 7)) +
  labs(
    x       = "Days after birth",
    y       = "Estimated hazard (per 1,000 person-days)",
    colour  = "Intervention",
    title   = "Smoothed Hazard Function — Neonatal Mortality",
    subtitle= "Kernel smoothing (local bandwidth, Müller-Wang)",
    caption = "Bandwidth selected using local bandwidth correction."
  )

ggsave("outputs/fig5_hazard_function.png", p_hazard,
       width = 9, height = 6, dpi = 300)
cat("  → Hazard function plot saved\n")

# =============================================================================
# SECTION 8 — LONGITUDINAL Hb TRAJECTORY PLOT
# =============================================================================
cat("\n── Section 8: Longitudinal Haemoglobin Trajectories ──\n")

hb <- read_csv("data/haemoglobin_longitudinal.csv", show_col_types = FALSE) %>%
  left_join(df %>% select(participant_id, intervention, neonatal_death,
                           anaemia_enrol, wealth_cat), by = "participant_id") %>%
  mutate(
    timepoint = factor(timepoint, levels = c(
      "Enrolment (~20wk)","28 weeks","34 weeks",
      "Delivery","6wk PP","6mo PP","12mo PP"
    ))
  )

hb_summary <- hb %>%
  group_by(intervention, timepoint) %>%
  summarise(
    mean_hb = mean(haemoglobin_g_dl, na.rm=TRUE),
    se_hb   = sd(haemoglobin_g_dl, na.rm=TRUE) / sqrt(n()),
    n       = n(),
    .groups = "drop"
  )

p_hb <- ggplot(hb_summary, aes(x = timepoint, y = mean_hb,
                                 colour = intervention, group = intervention)) +
  geom_line(lwd = 1.0) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_hb - 1.96*se_hb,
                    ymax = mean_hb + 1.96*se_hb), width = 0.2, lwd = 0.7) +
  geom_hline(yintercept = 11, lty = 2, colour = "grey50") +
  annotate("text", x = 0.6, y = 11.12, label = "Hb 11 g/dL threshold",
           size = 3, colour = "grey40", hjust = 0) +
  scale_colour_manual(values = unname(JiVitA_COLORS), name = "Intervention") +
  scale_y_continuous(limits = c(9.5, 12.5), breaks = seq(9.5, 12.5, 0.5)) +
  labs(
    x       = "Measurement timepoint",
    y       = "Mean haemoglobin (g/dL) ± 95% CI",
    title   = "Longitudinal Haemoglobin Trajectories by Intervention Arm",
    subtitle= "From enrolment (~20 weeks gestation) to 12 months postpartum",
    caption = "Dashed line: WHO anaemia threshold (11 g/dL) for pregnant/postpartum women."
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("outputs/fig6_hb_trajectories.png", p_hb,
       width = 10, height = 6, dpi = 300)
cat("  → Hb trajectory plot saved\n")

# =============================================================================
# SECTION 9 — COMPOSITE FIGURE (Publication-grade multi-panel)
# =============================================================================
cat("\n── Section 9: Composite Figure ──\n")

# Rebuild simpler versions for patchwork
km_df <- surv_summary(fit_km_int, data = df) %>%
  as_tibble() %>%
  mutate(
    intervention = if ("intervention" %in% names(.))
                     intervention
                   else str_remove(strata, "intervention="),
    surv_pct = (1 - surv) * 100
  )

p_km_gg <- ggplot(km_df, aes(x = time, y = surv_pct,
                              colour = intervention, fill = intervention)) +
  geom_step(lwd = 1.0) +
  geom_stepribbon(aes(ymin = (1-upper)*100, ymax = (1-lower)*100),
                  alpha = 0.12, colour = NA) +
  scale_colour_manual(values = unname(JiVitA_COLORS), name = "Intervention") +
  scale_fill_manual(  values = unname(JiVitA_COLORS), name = "Intervention") +
  scale_x_continuous(breaks = seq(0,28,7)) +
  labs(x = "Days after birth", y = "Cumulative mortality (%)",
       title = "A) Kaplan-Meier Curves") +
  theme(legend.position = "none")

p_hb_small <- p_hb + labs(title = "B) Haemoglobin Trajectories") +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4,"cm"),
        axis.text.x = element_text(angle=30, hjust=1, size=8))

p_hz_small <- p_hazard + labs(title = "C) Smoothed Hazard Functions") +
  theme(legend.position = "none")

p_forest_small <- p_forest +
  labs(title = "D) Adjusted Cox Hazard Ratios") +
  theme(legend.position = "none",
        axis.text.y = element_text(size=7.5),
        strip.text = element_text(size=7.5))

composite <- (p_km_gg | p_hz_small) /
             (p_hb_small | p_forest_small) +
  plot_annotation(
    title    = "Neonatal Mortality: Survival Analysis in a Cluster-RCT Cohort, Bangladesh",
    subtitle = "N = 5,000 pregnant women | 4 intervention arms | 2019–2023",
    caption  = "Methods: Kaplan-Meier, kernel-smoothed hazard, Cox PH with Efron ties, and Fine-Gray competing risks.\nAuthor: Hasan Mahmud Sujan",
    theme    = theme(
      plot.title    = element_text(face="bold", size=14, colour="#1F4E79"),
      plot.subtitle = element_text(colour="grey40"),
      plot.caption  = element_text(colour="grey50", size=8)
    )
  )

ggsave("outputs/fig_composite_publication.png", composite,
       width = 16, height = 14, dpi = 300, bg = "white")
cat("  → Composite publication figure saved\n")

# =============================================================================
# SECTION 10 — MODEL PERFORMANCE & DIAGNOSTICS SUMMARY
# =============================================================================
cat("\n── Section 10: Model diagnostics ──\n")

diag_df <- bind_rows(
  tidy(cox_unadj, exponentiate=TRUE, conf.int=TRUE) %>% mutate(model="Unadjusted Cox"),
  tidy(cox_adj,   exponentiate=TRUE, conf.int=TRUE) %>% mutate(model="Adjusted Cox"),
  tidy(cox_frailty, exponentiate=TRUE, conf.int=TRUE) %>% mutate(model="Frailty Cox")
) %>%
  filter(str_detect(term,"intervention")) %>%
  mutate(
    arm = str_remove(term,"intervention"),
    label = sprintf("HR %.2f (%.2f–%.2f)", estimate, conf.low, conf.high)
  )

p_hr_compare <- ggplot(diag_df,
  aes(x = estimate, y = arm, colour = model, shape = model)) +
  geom_vline(xintercept=1, lty=2, colour="grey60") +
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high),
                  position = position_dodge(0.5), size=0.7) +
  scale_x_log10(breaks=c(0.3,0.5,0.7,1.0,1.5),
                labels=c("0.3","0.5","0.7","1.0","1.5")) +
  scale_colour_brewer(palette="Dark2") +
  labs(x="Hazard Ratio (95% CI)", y="Intervention vs Control",
       colour="Model", shape="Model",
       title="E) HR Sensitivity: Unadjusted vs Adjusted vs Frailty Cox",
       caption="Frailty model accounts for cluster-level (village) random effects.") +
  theme(legend.position="bottom")

ggsave("outputs/fig7_hr_model_comparison.png", p_hr_compare,
       width=9, height=5, dpi=300)
cat("  → Model comparison plot saved\n")

cat("\n=================================================================\n")
cat("  ANALYSIS COMPLETE — all outputs in ./outputs/\n")
cat("=================================================================\n")
cat("\n  Files generated:\n")
fs::dir_ls("outputs") %>% walk(~cat(sprintf("  • %s\n", .x)))
