# =============================================================================
# shiny_app/app.R
# Interactive Survival Analysis Dashboard
# Neonatal Mortality — Cluster-RCT Cohort, Bangladesh
#
# Author: Hasan Mahmud Sujan
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(tidyverse)
  library(survival)
  library(survminer)
  library(gtsummary)
  library(gt)
  library(plotly)
  library(DT)
  library(broom)
})

# ── Load data (relative path — run from project root) ─────────────────────────
df_raw <- read_csv("../data/cohort_baseline.csv", show_col_types = FALSE)
hb_raw <- read_csv("../data/haemoglobin_longitudinal.csv", show_col_types = FALSE)

df <- df_raw %>%
  mutate(
    intervention = factor(intervention,
      levels = c("Control","BEP_supplement","MNP","BEP_plus_MNP"),
      labels = c("Control","BEP Supplement","MNP","BEP + MNP")),
    anaemia_enrol = factor(anaemia_enrol, levels=c("None","Mild","Moderate","Severe")),
    wealth_cat    = factor(wealth_index, labels=c("Poorest","Poor","Middle","Wealthy","Wealthiest")),
    lbw_factor    = factor(low_birthweight, labels=c("Normal BW","Low BW")),
    preterm_factor= factor(preterm, labels=c("Term","Preterm")),
    season        = factor(season, levels=c("Winter","PreMonsoon","Monsoon","PostMonsoon"))
  )

INT_COLORS <- c("Control"="#C00000","BEP Supplement"="#2E75B6",
                "MNP"="#70AD47","BEP + MNP"="#7030A0")

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════
ui <- page_navbar(
  title = div(
    img(src="", height="30px"),
    strong("Neonatal Survival Dashboard"),
    span(" — Cluster-RCT Bangladesh 2019–2023", style="font-size:0.85em;color:#aaa;")
  ),
  theme = bs_theme(
    bootswatch  = "flatly",
    primary     = "#1F4E79",
    secondary   = "#2E75B6",
    base_font   = font_google("Inter"),
    heading_font= font_google("Inter")
  ),
  bg = "#1F4E79",
  inverse = TRUE,

  # ── Tab 1: Overview ─────────────────────────────────────────────────────────
  nav_panel("📊 Overview",
    layout_columns(
      fill = FALSE,
      value_box("Total participants", nrow(df),       showcase=bsicons::bs_icon("people-fill"),     theme="primary"),
      value_box("Neonatal deaths",    sum(df$neonatal_death),
                showcase=bsicons::bs_icon("heart-pulse"), theme="danger"),
      value_box("Mortality rate",
                paste0(round(mean(df$neonatal_death)*1000,1)," per 1,000"),
                showcase=bsicons::bs_icon("graph-down"),   theme="warning"),
      value_box("Study clusters",     n_distinct(df$cluster_id),
                showcase=bsicons::bs_icon("pin-map-fill"), theme="success")
    ),
    layout_columns(
      col_widths = c(8, 4),
      card(
        card_header("Baseline characteristics by intervention arm"),
        gt_output("tbl_baseline")
      ),
      card(
        card_header("Deaths by intervention arm"),
        plotOutput("bar_deaths", height="280px"),
        card_footer("Error bars: 95% CI using Wilson method")
      )
    )
  ),

  # ── Tab 2: Kaplan-Meier ─────────────────────────────────────────────────────
  nav_panel("📈 Kaplan-Meier",
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        selectInput("km_strata", "Stratify by:",
          choices = c("Intervention"="intervention",
                      "Anaemia status"="anaemia_enrol",
                      "Birthweight"="lbw_factor",
                      "Preterm"="preterm_factor",
                      "Season"="season",
                      "Wealth quintile"="wealth_cat"),
          selected = "intervention"
        ),
        checkboxInput("km_ci",     "Show 95% CI",         value=TRUE),
        checkboxInput("km_risktbl","Show risk table",      value=TRUE),
        checkboxInput("km_pval",   "Show log-rank p-value",value=TRUE),
        hr(),
        sliderInput("km_tau","Restrict to (days):",min=1,max=28,value=28),
        hr(),
        p(em("Tip: stratify by intervention to assess treatment effect."),
          style="font-size:0.85em;color:grey;")
      ),
      card(
        card_header("Kaplan-Meier Survival Curve"),
        plotOutput("km_plot", height="520px")
      )
    )
  ),

  # ── Tab 3: Cox Regression ───────────────────────────────────────────────────
  nav_panel("🔬 Cox Regression",
    layout_columns(
      col_widths=c(4,8),
      card(
        card_header("Model covariates"),
        checkboxGroupInput("cox_vars","Include in model:",
          choices = c(
            "Intervention arm"    = "intervention",
            "Low birthweight"     = "low_birthweight",
            "Preterm birth"       = "preterm",
            "Maternal anaemia"    = "anaemia_enrol",
            "Maternal age"        = "age_enrol",
            "Parity"              = "parity",
            "Education (years)"   = "edu_years",
            "Wealth index"        = "wealth_index",
            "Food insecurity"     = "food_insecurity",
            "Skilled attendant"   = "skilled_birth_attendant",
            "ANC visits"          = "anc_visits",
            "Child sex"           = "child_sex",
            "Season"              = "season"
          ),
          selected = c("intervention","low_birthweight","anaemia_enrol",
                       "edu_years","wealth_index","skilled_birth_attendant")
        ),
        actionButton("run_cox","▶ Fit Model", class="btn-primary w-100"),
        hr(),
        verbatimTextOutput("cox_fit_info")
      ),
      card(
        card_header("Forest Plot — Adjusted Hazard Ratios"),
        plotOutput("forest_plot", height="500px")
      )
    ),
    card(
      card_header("Regression Table"),
      gt_output("cox_table")
    )
  ),

  # ── Tab 4: Longitudinal Hb ──────────────────────────────────────────────────
  nav_panel("🩸 Haemoglobin",
    layout_columns(
      col_widths=c(3,9),
      card(
        card_header("Filters"),
        checkboxGroupInput("hb_arms","Intervention arms:",
          choices = levels(df$intervention),
          selected= levels(df$intervention)
        ),
        checkboxGroupInput("hb_anaemia","Baseline anaemia:",
          choices = levels(df$anaemia_enrol),
          selected= levels(df$anaemia_enrol)
        ),
        radioButtons("hb_outcome","Overlay neonatal death:",
                     c("Both","Survived","Died"), inline=TRUE)
      ),
      card(
        card_header("Longitudinal Haemoglobin Trajectory"),
        plotlyOutput("hb_plot", height="420px"),
        card_footer("Points: group mean ± 95% CI. Dashed line: Hb 11 g/dL threshold.")
      )
    )
  ),

  # ── Tab 5: Data Explorer ────────────────────────────────────────────────────
  nav_panel("🗂️ Data Explorer",
    card(
      card_header("Cohort Data (first 1,000 rows)"),
      DTOutput("data_table")
    )
  ),

  # ── Footer nav ──────────────────────────────────────────────────────────────
  nav_spacer(),
  nav_item(tags$a(
    href="https://github.com/hasan-statx",
    target="_blank",
    bsicons::bs_icon("github"), "GitHub"
  ))
)

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Overview: baseline table ───────────────────────────────────────────────
  output$tbl_baseline <- render_gt({
    df %>%
      select(intervention, age_enrol, parity, edu_years, wealth_index,
             haemoglobin_enrol, birthweight_g, low_birthweight,
             preterm, skilled_birth_attendant, neonatal_death) %>%
      tbl_summary(
        by = intervention,
        statistic = list(
          all_continuous()  ~ "{mean} ({sd})",
          all_categorical() ~ "{n} ({p}%)"
        ),
        digits = all_continuous() ~ 1,
        label = list(
          age_enrol              ~ "Age at enrolment (years)",
          parity                 ~ "Parity",
          edu_years              ~ "Education (years)",
          wealth_index           ~ "Wealth index",
          haemoglobin_enrol      ~ "Haemoglobin at enrolment (g/dL)",
          birthweight_g          ~ "Birthweight (g)",
          low_birthweight        ~ "Low birthweight (<2,500g)",
          preterm                ~ "Preterm birth (<37 wk)",
          skilled_birth_attendant~ "Skilled birth attendant",
          neonatal_death         ~ "Neonatal death (day 0–28)"
        )
      ) %>%
      add_p() %>%
      bold_p(t=0.05) %>%
      as_gt() %>%
      tab_options(table.font.size=px(13))
  })

  # ── Overview: bar chart ────────────────────────────────────────────────────
  output$bar_deaths <- renderPlot({
    df %>%
      group_by(intervention) %>%
      summarise(
        n = n(), deaths = sum(neonatal_death),
        rate = deaths/n,
        se   = sqrt(rate*(1-rate)/n),
        lo   = pmax(0, rate - 1.96*se),
        hi   = rate + 1.96*se
      ) %>%
      ggplot(aes(x=intervention, y=rate*1000, fill=intervention)) +
      geom_col(width=0.6, alpha=0.9) +
      geom_errorbar(aes(ymin=lo*1000,ymax=hi*1000),width=0.2,lwd=0.8) +
      scale_fill_manual(values=INT_COLORS, guide="none") +
      labs(x=NULL, y="Deaths per 1,000 live births") +
      theme_minimal(base_size=12) +
      theme(axis.text.x=element_text(angle=20,hjust=1))
  })

  # ── KM plot ────────────────────────────────────────────────────────────────
  output$km_plot <- renderPlot({
    formula_str <- paste0("Surv(surv_time_days, neonatal_death) ~ ", input$km_strata)
    fit <- survfit(as.formula(formula_str), data = df)

    ggsurvplot(
      fit,
      data          = df,
      fun           = "event",
      conf.int      = input$km_ci,
      risk.table    = input$km_risktbl,
      risk.table.height = if(input$km_risktbl) 0.28 else 0,
      pval          = input$km_pval,
      surv.scale    = "percent",
      xlab          = "Days after birth",
      ylab          = "Cumulative neonatal mortality (%)",
      xlim          = c(0, input$km_tau),
      break.time.by = 7,
      ggtheme       = theme_classic(base_size=13),
      legend.title  = input$km_strata
    )
  }, res=110)

  # ── Cox model (reactive) ───────────────────────────────────────────────────
  cox_model <- eventReactive(input$run_cox, {
    req(length(input$cox_vars) >= 1)
    formula_str <- paste0(
      "Surv(surv_time_days, neonatal_death) ~ ",
      paste(input$cox_vars, collapse=" + ")
    )
    coxph(as.formula(formula_str), data=df, ties="efron")
  })

  output$cox_fit_info <- renderPrint({
    req(cox_model())
    m <- cox_model()
    cat(sprintf("Concordance : %.3f\n", summary(m)$concordance[1]))
    cat(sprintf("Likelihood ratio test p : %.4f\n",
                summary(m)$logtest["pvalue"]))
    cat(sprintf("n events : %d / %d\n", m$nevent, nrow(df)))
  })

  output$forest_plot <- renderPlot({
    req(cox_model())
    tidy(cox_model(), exponentiate=TRUE, conf.int=TRUE) %>%
      filter(!str_detect(term,"Intercept")) %>%
      mutate(term = str_replace_all(term, c(
        "intervention"="Int: ","low_birthweight"="LBW","preterm"="Preterm",
        "anaemia_enrol"="Anaemia: ","edu_years"="Education","wealth_index"="Wealth",
        "food_insecurity"="Food insecurity","skilled_birth_attendant"="SBA",
        "anc_visits"="ANC visits","child_sex"="Sex: ","season"="Season: "
      ))) %>%
      ggplot(aes(x=estimate, y=fct_reorder(term,estimate),
                 colour=p.value<0.05)) +
      geom_vline(xintercept=1,lty=2,colour="grey50") +
      geom_pointrange(aes(xmin=conf.low,xmax=conf.high),size=0.6) +
      scale_x_log10(breaks=c(0.25,0.5,1,2,4),
                    labels=c("0.25","0.50","1.0","2.0","4.0")) +
      scale_colour_manual(values=c("TRUE"="#1F4E79","FALSE"="grey60"),
                          labels=c("TRUE"="p<0.05","FALSE"="p≥0.05"),
                          name="") +
      labs(x="Hazard Ratio (95% CI) — log scale", y=NULL,
           caption="Solid: p<0.05 | Open: p≥0.05") +
      theme_classic(base_size=12) +
      theme(legend.position="bottom")
  })

  output$cox_table <- render_gt({
    req(cox_model())
    tbl_regression(cox_model(), exponentiate=TRUE) %>%
      bold_p(t=0.05) %>%
      add_significance_stars() %>%
      as_gt() %>%
      tab_options(table.font.size=px(13))
  })

  # ── Hb trajectory ──────────────────────────────────────────────────────────
  output$hb_plot <- renderPlotly({
    hb_filt <- hb_raw %>%
      left_join(df %>% select(participant_id, intervention,
                               anaemia_enrol, neonatal_death), by="participant_id") %>%
      filter(
        intervention %in% input$hb_arms,
        anaemia_enrol %in% input$hb_anaemia
      )

    if(input$hb_outcome == "Survived") hb_filt <- filter(hb_filt, neonatal_death==0)
    if(input$hb_outcome == "Died")     hb_filt <- filter(hb_filt, neonatal_death==1)

    tp_order <- c("Enrolment (~20wk)","28 weeks","34 weeks",
                  "Delivery","6wk PP","6mo PP","12mo PP")
    hb_filt <- hb_filt %>%
      mutate(timepoint = factor(timepoint, levels=tp_order))

    summ <- hb_filt %>%
      group_by(intervention, timepoint) %>%
      summarise(
        mean_hb = mean(haemoglobin_g_dl, na.rm=TRUE),
        se      = sd(haemoglobin_g_dl, na.rm=TRUE)/sqrt(n()),
        .groups = "drop"
      )

    p <- ggplot(summ, aes(x=timepoint, y=mean_hb,
                           colour=intervention, group=intervention,
                           text=sprintf("%s<br>Arm: %s<br>Hb: %.2f g/dL",
                                        timepoint, intervention, mean_hb))) +
      geom_line(lwd=1.1) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin=mean_hb-1.96*se, ymax=mean_hb+1.96*se),
                    width=0.2, lwd=0.7) +
      geom_hline(yintercept=11, lty=2, colour="grey60") +
      scale_colour_manual(values=INT_COLORS) +
      labs(x=NULL, y="Mean Hb (g/dL) ± 95% CI", colour="Intervention") +
      theme_minimal(base_size=12) +
      theme(axis.text.x=element_text(angle=25,hjust=1))

    ggplotly(p, tooltip="text")
  })

  # ── Data table ─────────────────────────────────────────────────────────────
  output$data_table <- renderDT({
    df %>%
      select(participant_id, cluster_id, intervention, district,
             age_enrol, parity, edu_years, wealth_index,
             haemoglobin_enrol, anaemia_enrol, birthweight_g,
             low_birthweight, preterm, skilled_birth_attendant,
             surv_time_days, neonatal_death) %>%
      slice_head(n=1000) %>%
      datatable(
        filter="top", rownames=FALSE,
        options=list(pageLength=15, scrollX=TRUE,
                     dom="Bfrtip", buttons=c("csv","excel")),
        extensions="Buttons"
      )
  })
}

# ── Run ────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
