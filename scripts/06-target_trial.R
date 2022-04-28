
# Libraries ---------------------------------------------------------------
library(tidyverse)
library(bernr)
library(here)
library(survminer)
library(survival)
library(ipw)
library(boot)
library(paletteer)
library(patchwork)
library(broom)
library(survey)

# Functions/Themes ----------------------------------------------------------
theme_set(theme_light(base_family = "Lato"))
theme_update(
  plot.title.position = "plot",
  axis.title.x = element_text(margin = margin(t = 1, unit = "line")),
  axis.title.y = element_text(margin = margin(r = 1, unit = "line")),
  plot.title = element_text(face = "bold", margin = margin(b = 1, unit = "line")),
  strip.background = element_rect(fill = "white"),
  strip.text = element_text(color = "black", face = "bold")
)



# Read data ---------------------------------------------------------------

df_analysis <- pro_read("05-analysis_data.rds")




df_analysis %>% 
  distinct(id, event, event_type, group) %>% 
  count(group, event, event_type)
  

df_analysis %>% 
  distinct(id, event_dicho, event_type, group) %>% 
  count(group, event_dicho) %>% 
  group_by(group) %>% 
  mutate(total = sum(n), 
         p = n / total)

df_analysis %>% 
  distinct(id, event_3, event_type, group) %>% 
  count(group, event_3) %>% 
  group_by(group) %>% 
  mutate(total = sum(n), 
         p = n / total)



# Kaplan Meier Plot -------------------------------------------------------


p <- ggsurvplot(fit = survfit(Surv(time, event_dicho) ~ group, data = df_analysis),
                ggtheme = theme_light(base_family = "Lato"),
                palette = paletteer_d(`"ggsci::default_jama"`),
                xlim = c(0, 7),
                risk.table = TRUE,
                legend = "right",
                legend.title = "",
                legend.labs = c("Current ART", "Switch to INSTI-based regimen"),
                censor = FALSE,
                break.time.by = 1,
                fontsize = 3,
                size = 0.7,
                risk.table.title = "",
                xlab = "\nYears since Baseline",
                ylab = "Viral Failure\n",
                tables.y.text.col = FALSE,
                fun = "event",
                conf.int = FALSE
)

kaplan_p <- p$plot + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 1, color = "grey90", size = 0.1),
        legend.position = c(0.3, 0.8)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

kaplan_t <- p$table + 
  labs(y = NULL, subtitle = "Individuals at risk") +
  theme(panel.grid = element_blank())

kaplan_p / kaplan_t + 
  plot_layout(heights = c(5, 1)) &
  theme(
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(
      face = "bold",
      margin = margin(b = 1.24, unit = "line")
    ),
    plot.subtitle = element_text(
      face = "bold",
      margin = margin(b = 1.24, unit = "line")
    ),
    axis.title.x = element_text(margin = margin(t = 1.24, unit = "line"))
  )


ggsave(here("graphs", "06-kaplan_meier.png"), 
       bg = "white", width = 7, height = 5, dpi = 300)
ggsave(here("graphs", "06-kaplan_meier.pdf"), 
       bg = "white", width = 7, height = 5, dpi = 300, device = cairo_pdf)





# Inverse Probability Weighting -------------------------------------------


weights <- ipwpoint(exposure = group_numeric, 
                    family = "binomial", 
                    link = "logit", 
                    numerator = ~1, # Stabilized weights
                    denominator = ~adherence_locf + history_VF + n_conmeds, 
                    data = as.data.frame(df_analysis))

summary(weights$ipw.weights)

df_analysis$weights <- weights$ipw.weights






# Fit weighted model ------------------------------------------------------

m <- coxph(Surv(time, event_dicho) ~ group, data = df_analysis, 
           weights = df_analysis$weights)

broom::tidy(m, exp = TRUE)


# Calculate bootstrap confidence intervals
boot_function <- function(data, indices) {
  d <- data[indices, ]
  weights <- ipwpoint(exposure = group_numeric, 
                      family = "binomial", 
                      link = "logit", 
                      numerator = ~1,
                      denominator = ~adherence_locf + history_VF + n_conmeds + 
                        adherence_locf + history_VF + n_conmeds, 
                      data = as.data.frame(d))
  d$weights <- weights$ipw.weights
  
  obj <- coxph(Surv(time, event_dicho) ~ group, data = d, 
               weights = d$weights)
  broom::tidy(obj)$estimate
} 


# BOOT starts here
set.seed(1)
# b <- boot(df_analysis, boot_function, R = 500, parallel = "multicore", ncpus = 12)
# write_rds(b, here("processed", "06-vf_model_bootstraps.rds"))


b <- pro_read("06-vf_model_bootstraps.rds")

tibble("HR" = b$t0,
       "LCI" = quantile(b$t, 0.025)[[1]],
       "UCI" = quantile(b$t, 0.975)[[1]]) %>%
  mutate_all(~scales::comma(exp(.x), accuracy = 0.01))



### DO THE SAME WITH any HIV > 300

# Kaplan Meier Plot -------------------------------------------------------


p <- ggsurvplot(fit = survfit(Surv(time, event_3) ~ group, data = df_analysis),
                ggtheme = theme_light(base_family = "Lato"),
                palette = paletteer_d(`"ggsci::default_jama"`),
                xlim = c(0, 7),
                risk.table = TRUE,
                legend = "right",
                legend.title = "",
                legend.labs = c("Current ART", "Switch to INSTI-based regimen"),
                censor = FALSE,
                break.time.by = 1,
                fontsize = 3,
                size = 0.7,
                risk.table.title = "",
                xlab = "\nYears since Baseline",
                ylab = "Viral Failure\n",
                tables.y.text.col = FALSE,
                # fun = "event",
                conf.int = FALSE
)

kaplan_p <- p$plot + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 1, color = "grey90", size = 0.1),
        legend.position = c(0.3, 0.8)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

kaplan_t <- p$table + 
  labs(y = NULL, subtitle = "Individuals at risk") +
  theme(panel.grid = element_blank())

kaplan_p / kaplan_t + 
  plot_layout(heights = c(5, 1)) &
  theme(
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(
      face = "bold",
      margin = margin(b = 1.24, unit = "line")
    ),
    plot.subtitle = element_text(
      face = "bold",
      margin = margin(b = 1.24, unit = "line")
    ),
    axis.title.x = element_text(margin = margin(t = 1.24, unit = "line"))
  )


ggsave(here("graphs", "06-kaplan_meier_any_detection.png"), 
       bg = "white", width = 7, height = 5, dpi = 300)
ggsave(here("graphs", "06-kaplan_meier_any_detection.pdf"), 
       bg = "white", width = 7, height = 5, dpi = 300, device = cairo_pdf)





# Inverse Probability Weighting -------------------------------------------


weights <- ipwpoint(exposure = group_numeric, 
                    family = "binomial", 
                    link = "logit", 
                    numerator = ~1, # Stabilized weights
                    denominator = ~female + adherence_locf + history_VF + n_conmeds, 
                    data = as.data.frame(df_analysis))

summary(weights$ipw.weights)

df_analysis$weights <- weights$ipw.weights


?ipwpoint



# Fit weighted model ------------------------------------------------------

m <- coxph(Surv(time, event_3) ~ group, data = df_analysis, 
           weights = df_analysis$weights)

broom::tidy(m, exp = TRUE)


# Calculate bootstrap confidence intervals
boot_function <- function(data, indices) {
  d <- data[indices, ]
  weights <- ipwpoint(exposure = group_numeric, 
                      family = "binomial", 
                      link = "logit", 
                      numerator = ~1,
                      denominator = ~adherence_locf + history_VF + n_conmeds + 
                        adherence_locf + history_VF + n_conmeds, 
                      data = as.data.frame(d))
  d$weights <- weights$ipw.weights
  
  obj <- coxph(Surv(time, event_3) ~ group, data = d, 
               weights = d$weights)
  broom::tidy(obj)$estimate
} 


# BOOT starts here
set.seed(1)
# b <- boot(df_analysis, boot_function, R = 500, parallel = "multicore", ncpus = 12)
# write_rds(b, here("processed", "06-vf_model_bootstraps_any_detection.rds"))

b <- pro_read("06-vf_model_bootstraps_any_detection.rds")

tibble("HR" = b$t0,
       "LCI" = quantile(b$t, 0.025)[[1]],
       "UCI" = quantile(b$t, 0.975)[[1]]) %>%
  mutate_all(~scales::comma(exp(.x), accuracy = 0.01))

svy_design <- svydesign(id = ~id, weights = ~weights, data = df_analysis)
svycoxph(Surv(time, event_dicho) ~ group,  design=svy_design) %>% 
  tidy(exp = TRUE, conf.int = TRUE) %>% 
  select(term, estimate, conf.low, conf.high, p.value)
coxph(Surv(time, event_dicho) ~ group,  data=df_analysis) %>% 
  tidy(exp = TRUE, conf.int = TRUE) %>% 
  select(term, estimate, conf.low, conf.high, p.value)
mod_wght <- svykm(Surv(time, event_dicho) ~ group, design=svy_design)

plot(mod_wght[[1]], col="black")
lines(mod_wght[[2]], col="red")

surv_curve <- tibble(time = c(mod_wght[[1]]$time, mod_wght[[2]]$time),
                     surv = c(mod_wght[[1]]$surv, mod_wght[[2]]$surv), 
                     treatment = rep(c("Current ART", "Switch to INSTI-based regimen"), 
                                     times =  c(length(mod_wght[[1]]$time), 
                                                length(mod_wght[[2]]$time))))

kaplan_p
ad_p <- surv_curve %>% 
  ggplot(aes(x = time, y = surv)) + 
  geom_line(aes(color = treatment)) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 1, color = "grey90", size = 0.1),
        legend.position = c(0.3, 0.8)) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, 1)) + 
  labs(color = NULL, 
       x = "Years since Baseline", 
       y = "Viral Failure")


kaplan_p + ad_p + 
  scale_color_paletteer_d("ggsci::default_jama") & 
  theme(legend.position = c(0.3, 0.5))

ggsave(here("graphs", "06-kaplan_adj_surv_any.png"), 
       dpi = 300, bg = "white", width = 10, height = 6)
