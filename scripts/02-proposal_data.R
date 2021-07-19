
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(bernr)
library(lubridate)
library(flextable)
library(tableone)
library(here)
library(data.table)
source("themes.R")



# Functions ---------------------------------------------------------------

# remove backbones from treatment string
drop_backbone <- function(string) {
  string <- str_remove_all(string, "ABC|3TC|ETC|TAF|TDF|COB|RTV|AZT|DDI|D4T")
  str_squish(string)
}



# Data --------------------------------------------------------------------


modif_art <- pro_read("modif_art.rds")
fup <- pro_read("fup.rds")
tail <- pro_read("tail.rds")
pat <- pro_read("pat.rds")
lab <- pro_read("lab.rds")


# Preparation -------------------------------------------------------------


# identify individuals with active follow-up after 1.11.13 (Stribild)
recent <- fup %>% 
  select(id, fupdate) %>% 
  arrange(id, desc(fupdate)) %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter(fupdate >= dmy("01.11.2013"))


# Calculate n of classes including NRTIs
art <- modif_art %>% 
  mutate(dur = as.numeric(enddate - moddate)) %>% 
  # Exclude time intervals where we have no data, but were less than 60 days
  filter(!(is.na(treatment) & enddate - moddate <= 60)) %>%
  arrange(id, moddate) %>% 
  group_by(id) %>% 
  mutate(current_art = last(treatment)) %>% 
  ungroup() %>% 
  mutate(num_pi = if_else(str_detect(treatment, "(RTV)"),
                          num_pi - 1, num_pi)) %>%
  mutate(num_others = if_else(str_detect(treatment, "(COB)"),
                              num_others - 1, num_others)) %>%
  mutate(classes_n = (num_nrti != 0 | num_ntrti != 0) + 
           (num_nnrti != 0) + 
           (num_pi != 0) + 
           (num_inti !=0) + 
           (num_others != 0) + 
           (num_fi != 0)) %>% 
  select(-(num_art:precision))




# Treated:
# Patients who switched from 2 to 1 anchor agent after 11/13, 
# and the regimen they switched to was one of the EACS recommended ones.

# Treated are those who switched from 3 to 2 classes
# Patients needed to be on the 2-class regimen for at least 30 days
treat <- art %>%
  mutate(flag = if_else(classes_n == 2 & 
                          lag(classes_n > 2) & 
                          moddate >= dmy("01.11.2013") & 
                          # at least 30 days or NA (means current)
                          (dur >= 30 | is.na(dur)),  
                        "X", "")) %>% 
  group_by(id) %>% 
  mutate(switch = any(flag == "X"), 
         switch_date = if_else(flag == "X", moddate, NA_real_)) %>% 
  filter(any(flag == "X")) %>% 
  # Check regimens that were before and after switch
  mutate(simple_reg = if_else(flag == "X", treatment, NA_character_), 
         comp_reg = if_else(lead(flag) == "X", treatment, NA_character_)) %>% 
  arrange(id, simple_reg, moddate) %>% 
  mutate(simple_reg = first(simple_reg)) %>% 
  arrange(id, comp_reg, moddate) %>% 
  mutate(comp_reg = first(comp_reg)) %>%
  arrange(id, switch_date) %>% 
  mutate(switch_date = first(switch_date)) %>% 
  mutate(comp_reg_s = drop_backbone(comp_reg), 
         simple_reg_s = drop_backbone(simple_reg)) %>% 
  mutate(across(contains("_reg_s"), ~str_squish(.x))) %>% 
  arrange(id, moddate) %>% 
  slice(1) %>%
  ungroup() %>% 
  select(id, last_art = current_art, switch, switch_date, 
         simple_reg, simple_reg_s, 
         comp_reg, comp_reg_s)




# Control group
# More than 2 anchor agents, after 11/13

control <- art %>% 
  ungroup() %>% 
  mutate(flag = if_else(classes_n > 1 & 
                        moddate >= dmy("01.11.2013"), 
                        "X", "")) %>% 
  group_by(id) %>% 
  filter(any(flag == "X")) %>% 
  mutate(switch = FALSE) %>% 
  arrange(id, moddate) %>% 
  ungroup()


# There are some individuals in the control group. This is because they simplified, 
# but later had an intensification of their treatment, and therefore the current
# ART regimen contains >= classes. They remain in the "treat" group. I therefore
# remove them from the controls

control <- control %>% 
  filter(!(id %in% treat$id))




# Baseline date -----------------------------------------------------------
# There is an error in ID 51171, it says 8020 instead fo 2020, correct it.
treat$switch_date[treat$id == 51171] <- treat$switch_date[treat$id == 51171] - years(6000)

# Same for ID 42485, it says 2917 instead of 2017
treat$switch_date[treat$id == 42485] <- treat$switch_date[treat$id == 42485] - years(900)


# I take a random sample of the switch dates and assign them to the control


set.seed(123)
switch_dates_ctrl <- tibble(id = unique(control$id), 
                            switch_date = sample(treat$switch_date, 
                                                 size = length(unique(control$id)), 
                                                 replace = TRUE))

# Now, add new baseline dates to control df, and then take closest treatment
# prior to baseline date

control_n1 <- control %>% 
  left_join(switch_dates_ctrl, by = "id") %>% 
  mutate(flag2 = as.numeric(switch_date - moddate)) %>%
  group_by(id, before = flag2 >= 0) %>% 
  mutate(flag3 = if_else(before == TRUE, min(flag2), NA_real_)) %>% 
  filter(flag2 == flag3) %>%
  ungroup() %>% 
  select(id, bl_treatment = treatment, bl_classes_n = classes_n, moddate, 
         switch, bl_date = switch_date, 
         days_started_bf_bl = flag2) %>% 
  mutate(bl_treatment_s = drop_backbone(bl_treatment)) %>% 
  filter(bl_classes_n > 1)
  

# Join treated and control together
full <- treat %>% 
  full_join(control_n1) %>% 
  relocate(bl_treatment_s, .after = bl_treatment)

# Add last fup to dataset
full <- full %>% 
  left_join(recent %>% select(id, fupdate)) %>% 
  rename(last_fup = fupdate)


# Assess the distribution of baseline dates
full <- full %>% 
  mutate(baseline_date = if_else(switch == TRUE, switch_date, bl_date))

full %>% 
  ggplot(aes(x = baseline_date)) + 
  geom_density(aes(fill = switch), alpha = 0.5, adjust = 0.5)




# add patient characteristics ---------------------------------------------

risk <- tail %>% 
  mutate(riskgroup = case_when(
    riskgroup %in% c("BLOOD", "PERINAT", "UNKNOWN") ~ "OTHER", 
    TRUE                                            ~ riskgroup), 
    riskgroup = factor(riskgroup, levels = c("MSM", "HET", "IDU", "OTHER"))) %>% 
  select(id, riskgroup)


# Last center and last source
source <- fup %>% 
  arrange(id, desc(fupdate)) %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(source = factor(source, 
                         levels = c(1, 2, 3), 
                         labels = c("Cohort center", "Hospital",
                                            "Physician")), 
         center = factor(center2, 
                          levels = c(10, 20, 30, 40, 50, 60, 70), 
                          labels = c("ZH", "BS", "BE", "GE", "LAU", 
                                     "LUG", "STG"))) %>% 
  select(id, center, source)


# Calculate number of changes in database and pre-treatment time
n_treatments <- modif_art %>% 
  left_join(full %>% select(id, baseline_date)) %>% 
  filter(id %in% full$id) %>% 
  filter(moddate < baseline_date) %>% 
  group_by(id) %>% 
  summarise(n_treatments = n(), 
            years_treatment = as.numeric(sum(enddate - moddate, 
                                            na.rm = TRUE)/365.25))


# Join together
pop_det <- pat %>% 
  mutate(male = if_else(sex == 1, 1, 0), 
         ethnicity = case_when(ethnicity == 1 ~ "White", 
                               ethnicity == 2 ~ "Black",
                               ethnicity %in% c(0, 3, 4) ~ "Other", 
                               TRUE ~ NA_character_)) %>%
  mutate(ethnicity = factor(ethnicity, levels = c("White", "Black", "Other"))) %>% 
  select(id, born, male, ethnicity) %>% 
  left_join(risk, by = "id") %>% 
  left_join(source, by = "id") %>% 
  left_join(n_treatments, by = "id")

full <- full %>% 
  left_join(pop_det, by = "id")



# Calculate age at baseline
full <- full %>% 
  mutate(age = year(baseline_date) - born) %>% 
  select(-born)


# Time since HIV diagnosis (if HIV posdate = NA, take regdate)
full <- full %>% 
  left_join(pat %>% select(id, hiv_posdate, regdate)) %>% 
  mutate(first_hiv = if_else(is.na(hiv_posdate), regdate, hiv_posdate)) %>% 
  select(-(hiv_posdate:regdate)) %>% 
  mutate(years_first_hiv = as.numeric(baseline_date - first_hiv) / 365.25)

# Baseline Treatment an N of anchor agents
full <- full %>% 
  mutate(baseline_treatment = if_else(switch == TRUE, comp_reg, bl_treatment),
         baseline_n_anchor = str_count(drop_backbone(baseline_treatment), "\\S+"))


# Baseline anchor agents lumped
full <- full %>% 
  mutate(baseline_anchors_fct = drop_backbone(baseline_treatment)) %>% 
  mutate(baseline_anchors_fct = fct_lump_n(baseline_anchors_fct, 10), 
         baseline_anchors_fct = fct_infreq(baseline_anchors_fct))

# Viral failure -----------------------------------------------------------

# Definition used from A. Scherrer in CID 2016:
# ---
# The high-risk group included patients who
# had ever experienced a virological failure or who were treated
# with single- or dual-NRTI therapy for >28 days. A virological
# failure was defined as either 2 consecutive viral loads >500
# HIV-1 RNA copies/mL or 1 viral load >500 HIV-1 RNA copies/mL followed 
# by a treatment change if the patient had experienced ≥180 days of 
# continuous ART or ≥90 days of ART if
# viral suppression was reached (<50 HIV-1 RNA copies/mL).



# Single or dual NRTI treatment in the past
nrti_mono <- modif_art %>% 
  filter(id %in% full$id) %>% 
  # add baseline date
  left_join(full %>% select(id, baseline_date), by = "id") %>% 
  group_by(id) %>% 
  mutate(no_third = (num_nnrti + num_pi + num_ntrti + 
                       num_fi + num_others + num_inti) == 0) %>% 
  summarise(nrti_mono = any(num_nrti <= 2 & 
                            no_third & 
                            !is.na(treatment) & 
                            moddate < baseline_date))



# Viral failure before baseline

# Prepare ART data
art <- modif_art %>% 
  filter(id %in% full$id) %>% 
  select(-(num_art:precision))

# Prepare RNA data (we don't have data for 2 individuals)
rna <- lab %>% 
  select(id, rna, labdate) %>% 
  drop_na() %>% 
  filter(id %in% full$id)


# Aggregate art data into periods
on_art <- art %>% 
  filter(id %in% full$id) %>%
  group_by(id) %>% 
  mutate(on_art = as.numeric(!is.na(treatment))) %>% 
  mutate(p = cumsum(on_art != lag(on_art, default = 2))) %>% 
  group_by(id, p) %>% 
  summarise(start = first(moddate), 
            stop = last(enddate), 
            on_art = max(on_art),
            .groups = "drop") %>% 
  mutate(stop = if_else(is.na(stop), dmy("01.01.2030"), stop))


# Join ART and RNA data
setDT(rna);setDT(on_art)

comb_art <- on_art[rna, 
                   on = .(id, start <= labdate, stop >= labdate), 
                   .(id, rna = i.rna, 
                     rna_date = i.labdate, 
                     art_period = x.p, 
                     on_art = x.on_art, 
                     period_start = x.start, 
                     period_stop = x.stop)] %>% 
  as_tibble() %>% 
  # add baseline date
  left_join(full %>% select(id, baseline_date), by = "id")


# Code failure as having 2x RNA > 500 while being on ART >= 180 days, 
# and only count failures PRIOR to baseline
failure_long <- comb_art %>% 
  group_by(id) %>% 
  mutate(flag = if_else(rna > 500 & lag(rna) > 500 & on_art == 1, 
                        "X", "")) %>% 
  group_by(id, art_period) %>% 
  mutate(t = as.numeric(rna_date - first(period_start)), 
         failure = if_else(flag == "X" & 
                           t >= 180 & 
                           rna_date <= baseline_date, "Y", "N"))

any_failure <- failure_long %>% 
  group_by(id) %>% 
  summarise(any_failure = any(failure == "Y")) %>% 
  mutate(any_failure = replace_na(any_failure, FALSE))
  


# Add NRTI mono and viral failure to full data set

full <- full %>% 
  left_join(nrti_mono, by = "id") %>% 
  left_join(any_failure, by = "id")



# Tables   -----------------------------------------------------------------


# Define variables for table
vars <- c("age", "male", "ethnicity", "riskgroup", "nrti_mono", "any_failure",
          "source", "n_treatments", "years_first_hiv", 
          "years_treatment", "baseline_n_anchor", "baseline_anchors_fct")
cat_vars <- c("male", "ethnicity", "riskgroup", "nrti_mono", "any_failure",
              "source", "baseline_anchors_fct")

# Create vector used later for indentation in flextable
cvars <- paste0(c(vars, "^n", "Prior ART changes"), collapse = "|^")


# Table 1
tab1 <- CreateTableOne(vars = vars, strata = "switch", factorVars = cat_vars, 
               data = full)
  
t_patchar <- (print(tab1, nonnormal = c("age", "n_treatments", "years_first_hiv",
                                        "years_treatment", "baseline_n_anchor"), 
                   printToggle = FALSE, contDigits = 1, dropEqual = TRUE) %>% 
  as_tibble(rownames = "Variable") %>% 
  select(-test) %>%
  mutate(Variable = if_else(str_detect(Variable, "n_treatments"), 
                            "Prior ART changes (median [IQR])", Variable)) %>% 
  rename("No simplification" = `FALSE`, 
         "Simplification" = `TRUE`) %>% 
  flextable() %>% 
  padding(i = ~ !str_detect(Variable, cvars), 
          j = 1, padding.left = 14) %>% 
  bold(i = ~ str_detect(Variable, cvars), 
       j = 1) %>% 
  f_theme_surial() %>% 
  align(j = 2:4, align = "center") %>% 
  autofit())


# ART regimes of those who had simplifications
t_switch <- full %>% 
  filter(switch == TRUE) %>% 
  count(comp_reg_s, simple_reg_s, sort = T) %>% 
  rename("Prior ART" = comp_reg_s, 
         "Simplified ART" = simple_reg_s) %>% 
  flextable() %>% 
  f_theme_surial()


# ART regimes of those who remained on complicated ART
t_control <- full %>% 
  filter(switch == FALSE) %>% 
  count(bl_treatment_s, sort = T) %>% 
  rename("Baseline ART" = bl_treatment_s) %>%
  flextable() %>% 
  f_theme_surial()


# Look at individuals who have DOR, as Gilles will be interested in those
modif_art %>% 
  filter(str_detect(treatment, "DOR")) %>% 
  arrange(id, desc(moddate)) %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(treatment, sort = T) %>% 
  flextable() %>% 
  f_theme_surial()

# Output ------------------------------------------------------------------


# Write tables out to word
t_patchar %>% 
  set_caption("Overview of Patients") %>% 
  save_as_docx(path = here("tables", "02-overview_of_patients.docx"))

t_switch %>% 
  set_caption("Anchor agents of patients who had simplification") %>% 
  save_as_docx(path = here("tables", "02-art_of_simplifications.docx"))

t_control %>% 
  set_caption("Anchor agents of patients who did not simplify") %>% 
  save_as_docx(path = here("tables", "02-art_of_complicated_art.docx"))



# Write dataframe as rds
write_rds(full, here("processed", "02-proposal_population.rds"))
