# =============================================================================
# 01_generate_data.R
# Synthetic Longitudinal Maternal & Child Health Cohort
# Inspired by icddr,b / JiVitA-style surveillance systems
#
# Study Design:
#   A cluster-randomised prospective cohort of pregnant women followed
#   from enrolment (~20 weeks gestation) to 12 months postpartum.
#   Primary outcome : neonatal mortality (death within 28 days of birth)
#   Secondary outcomes: preterm birth, low birthweight, maternal anaemia
#
# Author : Hasan Mahmud Sujan
# =============================================================================

set.seed(2024)
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(simsurv)      # parametric survival data generation
})

cat("=================================================================\n")
cat("  GENERATING SYNTHETIC MATERNAL & CHILD HEALTH COHORT\n")
cat("=================================================================\n\n")

# ── Parameters ────────────────────────────────────────────────────────────────
N          <- 5000        # total participants
N_CLUSTERS <- 50          # villages/clusters
ENROL_START <- as.Date("2019-01-01")
ENROL_END   <- as.Date("2022-12-31")

DISTRICTS <- c("Gaibandha","Nilphamari","Lalmonirhat",
                "Kurigram","Rangpur","Sirajganj")

INTERVENTIONS <- c("Control","BEP_supplement","MNP","BEP_plus_MNP")

# ── 1. Cluster-level data ─────────────────────────────────────────────────────
clusters <- tibble(
  cluster_id    = sprintf("C%03d", 1:N_CLUSTERS),
  district      = sample(DISTRICTS, N_CLUSTERS, replace = TRUE),
  intervention  = rep(INTERVENTIONS, length.out = N_CLUSTERS),
  cluster_size  = sample(80:120, N_CLUSTERS, replace = TRUE),
  rural         = rbinom(N_CLUSTERS, 1, 0.82),
  health_facility_dist_km = rexp(N_CLUSTERS, rate = 0.15) + 0.5,
  cluster_food_insecurity = rbeta(N_CLUSTERS, 2, 5),  # cluster-level FI score
  cluster_re    = rnorm(N_CLUSTERS, 0, 0.3)            # random effect
)

# ── 2. Individual-level baseline data ─────────────────────────────────────────
cat("[1] Generating baseline participant data...\n")

indiv <- tibble(
  participant_id = sprintf("JVT-%05d", 1:N),
  cluster_id     = sample(clusters$cluster_id, N,
                          prob = clusters$cluster_size / sum(clusters$cluster_size),
                          replace = TRUE)
) %>%
  left_join(clusters, by = "cluster_id") %>%
  mutate(
    # Sociodemographic
    age_enrol        = round(rnorm(N, 24.5, 4.8), 1) %>% pmax(15) %>% pmin(45),
    parity           = sample(0:6, N, prob = c(.28,.30,.22,.12,.05,.02,.01), replace = TRUE),
    edu_years        = pmax(0, round(rnorm(N, 5.2, 3.8))),
    husband_edu_yrs  = pmax(0, round(rnorm(N, 6.1, 4.0))),
    wealth_index     = sample(1:5, N, prob = c(.28,.25,.22,.15,.10), replace = TRUE),
    food_insecurity  = rbeta(N, 2, 5) + cluster_food_insecurity * 0.3,
    food_insecurity  = pmin(food_insecurity, 1),

    # Anthropometric at enrolment (~20wk gestation)
    muac_mm          = round(rnorm(N, 228, 22)),   # mid-upper arm circumference
    height_cm        = round(rnorm(N, 150.5, 5.8), 1),
    weight_kg_enrol  = round(rnorm(N, 48.5 + wealth_index * 0.8, 7.2), 1),
    bmi_enrol        = round(weight_kg_enrol / (height_cm / 100)^2, 1),

    # Clinical
    haemoglobin_enrol = round(rnorm(N, 10.8, 1.9), 1) %>% pmax(5) %>% pmin(17),
    anaemia_enrol    = case_when(
      haemoglobin_enrol < 7  ~ "Severe",
      haemoglobin_enrol < 10 ~ "Moderate",
      haemoglobin_enrol < 11 ~ "Mild",
      TRUE                   ~ "None"
    ),
    systolic_bp_enrol = round(rnorm(N, 108, 12)),
    gravida           = parity + 1L,

    # Gestation & timing
    gestational_age_enrol_wk = round(rnorm(N, 20, 3)) %>% pmax(13) %>% pmin(28),
    enrol_date       = sample(seq(ENROL_START, ENROL_END, by = "day"), N, replace = TRUE),
    birth_date       = enrol_date +
                         days(round((40 - gestational_age_enrol_wk) * 7 +
                                    rnorm(N, 0, 7))),

    # Compute gestation at birth
    gestational_age_birth_wk = round(40 + rnorm(N, 0, 2.1)),
    preterm          = gestational_age_birth_wk < 37,

    # Intervention assignment (cluster-level)
    intervention     = factor(intervention, levels = INTERVENTIONS),
    compliance       = case_when(
      intervention == "Control"     ~ 0,
      TRUE ~ rbeta(N, 5, 2)        # compliance 0-1
    ),

    # Season of enrolment
    enrol_month = month(enrol_date),
    season = case_when(
      enrol_month %in% 6:9  ~ "Monsoon",
      enrol_month %in% 10:11 ~ "PostMonsoon",
      enrol_month %in% 12:2  ~ "Winter",
      TRUE                   ~ "PreMonsoon"
    ) %>% factor(levels = c("Winter","PreMonsoon","Monsoon","PostMonsoon")),

    # ANC visits
    anc_visits = pmin(8, rpois(N, lambda = 3.2 + edu_years * 0.08 +
                                  wealth_index * 0.25)),
    skilled_birth_attendant = rbinom(N, 1,
      prob = plogis(-1.2 + 0.08 * edu_years + 0.25 * wealth_index +
                    0.3 * (rural == 0) + 0.4 * (anc_visits >= 4)))
  )

# ── 3. Birth outcomes ──────────────────────────────────────────────────────────
cat("[2] Simulating birth outcomes...\n")

indiv <- indiv %>%
  mutate(
    # Birthweight (g) — affected by intervention, anaemia, nutrition, etc.
    birthweight_g = round(
      3000
      - 50  * (anaemia_enrol == "Severe")
      - 25  * (anaemia_enrol == "Moderate")
      - 60  * (food_insecurity > 0.6)
      + 80  * (intervention == "BEP_supplement")
      + 50  * (intervention == "MNP")
      + 120 * (intervention == "BEP_plus_MNP")
      - 120 * preterm
      + 20  * compliance
      + rnorm(N, 0, 320)
      + cluster_re * 60
    ) %>% pmax(500) %>% pmin(5500),

    low_birthweight  = birthweight_g < 2500,
    very_lbw         = birthweight_g < 1500,

    # Sex of child
    child_sex = sample(c("Male","Female"), N, replace = TRUE, prob = c(.514,.486)),

    # APGAR score
    apgar_5min = pmin(10, pmax(0, round(rnorm(N, 8.2, 1.5))))
  )

# ── 4. Survival outcome: neonatal mortality (days 0–28) ─────────────────────
cat("[3] Simulating neonatal survival times...\n")

# True log-hazard model:
#   ln h(t) = ln h0(t)
#             + β1*low_birthweight + β2*preterm + β3*anaemia_severe
#             + β4*food_insecurity + β5*skilled_birth_attendant
#             + β6*(intervention==BEP) + β7*(sex==Male) + random_effect

covdf <- indiv %>%
  transmute(
    lbw          = as.numeric(low_birthweight),
    vlbw         = as.numeric(very_lbw),
    preterm      = as.numeric(preterm),
    anemia_sev   = as.numeric(anaemia_enrol == "Severe"),
    anemia_mod   = as.numeric(anaemia_enrol == "Moderate"),
    food_ins_hi  = as.numeric(food_insecurity > 0.6),
    sba          = as.numeric(skilled_birth_attendant),
    bep          = as.numeric(intervention %in% c("BEP_supplement","BEP_plus_MNP")),
    mnp          = as.numeric(intervention %in% c("MNP","BEP_plus_MNP")),
    male         = as.numeric(child_sex == "Male"),
    wealth       = (wealth_index - 3) / 2,
    re           = indiv$cluster_re
  )

surv_dat <- simsurv(
  dist      = "weibull",
  lambdas   = 0.003,
  gammas    = 0.75,
  betas     = c(lbw=1.20, vlbw=0.60, preterm=0.90, anemia_sev=0.55,
                anemia_mod=0.25, food_ins_hi=0.35, sba=-0.70,
                bep=-0.50, mnp=-0.25, male=0.15, wealth=-0.20, re=0.40),
  x         = covdf,
  maxt      = 28,
  interval  = c(1e-8, 200)
)

indiv <- indiv %>%
  bind_cols(
    surv_dat %>% select(eventtime, status) %>%
      rename(surv_time_days = eventtime, neonatal_death = status)
  ) %>%
  mutate(
    surv_time_days  = round(surv_time_days, 2),
    neonatal_death  = as.integer(neonatal_death),
    # Postneonatal follow-up (days 29–365)
    postneonatal_death = case_when(
      neonatal_death == 1 ~ 0L,
      TRUE ~ rbinom(N, 1,
        prob = plogis(-4.5 + 0.5*low_birthweight + 0.3*food_insecurity
                      - 0.4*sba + cluster_re*0.3))
    ),
    postneonatal_death_day = if_else(
      postneonatal_death == 1,
      round(runif(N, 29, 365)), NA_real_
    ),
    # Combined under-5 event
    any_death     = pmax(neonatal_death, postneonatal_death),
    death_day     = case_when(
      neonatal_death == 1 ~ surv_time_days,
      postneonatal_death == 1 ~ postneonatal_death_day,
      TRUE ~ 365
    )
  )

# ── 5. Time-varying: repeated Hb measurements ─────────────────────────────────
cat("[4] Generating time-varying haemoglobin measurements...\n")

timepoints_label <- c("Enrolment (~20wk)","28 weeks","34 weeks","Delivery","6wk PP","6mo PP","12mo PP")
timepoints_day   <- c(-140, -84, -42, 0, 42, 182, 365)  # days relative to birth

hb_long <- map2_dfr(timepoints_label, timepoints_day, function(tp, day) {
  indiv %>%
    filter(!is.na(haemoglobin_enrol)) %>%
    transmute(
      participant_id,
      timepoint       = tp,
      day_relative_birth = day,
      haemoglobin_g_dl = case_when(
        tp == "Enrolment (~20wk)" ~ haemoglobin_enrol,
        tp == "28 weeks"  ~ haemoglobin_enrol + rnorm(N, 0.1, 0.5),
        tp == "34 weeks"  ~ haemoglobin_enrol + rnorm(N, 0.3, 0.6)
                           + 0.3*(intervention %in% c("MNP","BEP_plus_MNP")),
        tp == "Delivery"  ~ haemoglobin_enrol + rnorm(N, -0.2, 0.7),
        tp == "6wk PP"    ~ haemoglobin_enrol + rnorm(N, 0.5, 0.8)
                           + 0.4*(intervention %in% c("MNP","BEP_plus_MNP")),
        tp == "6mo PP"    ~ haemoglobin_enrol + rnorm(N, 0.8, 0.9)
                           + 0.6*(intervention %in% c("MNP","BEP_plus_MNP")),
        tp == "12mo PP"   ~ haemoglobin_enrol + rnorm(N, 1.0, 1.0)
                           + 0.7*(intervention %in% c("MNP","BEP_plus_MNP")),
      ) %>% round(1) %>% pmax(5) %>% pmin(17),
      measured = rbinom(N, 1, prob = 0.88)  # 88% measurement rate
    ) %>%
    filter(measured == 1) %>%
    select(-measured)
})

# ── 6. Save ────────────────────────────────────────────────────────────────────
cat("[5] Saving datasets...\n")

# Main cohort
write_csv(indiv,    "data/cohort_baseline.csv")

# Haemoglobin longitudinal
write_csv(hb_long,  "data/haemoglobin_longitudinal.csv")

# Cluster data
write_csv(clusters, "data/cluster_data.csv")

cat(sprintf("  cohort_baseline.csv      : %d rows × %d cols\n", nrow(indiv), ncol(indiv)))
cat(sprintf("  haemoglobin_longitudinal : %d rows × %d cols\n", nrow(hb_long), ncol(hb_long)))
cat(sprintf("  cluster_data             : %d rows × %d cols\n", nrow(clusters), ncol(clusters)))
cat(sprintf("\n  Neonatal mortality rate  : %.1f per 1,000 live births\n",
            mean(indiv$neonatal_death) * 1000))
cat(sprintf("  Low birthweight rate     : %.1f%%\n", mean(indiv$low_birthweight) * 100))
cat(sprintf("  Preterm birth rate       : %.1f%%\n", mean(indiv$preterm) * 100))
cat("\n✔ Data generation complete.\n")
