
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
  string <- str_remove_all(string, "ABC|3TC|ETC|TAF|TDF|COB|RTV|AZT")
  str_squish(string)
}



# Data --------------------------------------------------------------------


modif_art <- pro_read("modif_art.rds")
fup <- pro_read("fup.rds")
tail <- pro_read("tail.rds")
pat <- pro_read("pat.rds")
lab <- pro_read("lab.rds")


# Preparation -------------------------------------------------------------


# Restrict analysis to individuals with active follow-up after 1.1.20
recent <- fup %>% 
  select(id, fupdate) %>% 
  arrange(id, desc(fupdate)) %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter(fupdate >= dmy("01.01.2020"))


# Calculate number of classes
art <- modif_art %>% 
  mutate(dur = as.numeric(enddate - moddate)) %>% 
  # Exclude time intervals where we have no data, but were less than 60 days
  filter(!(is.na(treatment) & enddate - moddate <= 60)) %>%
  arrange(id, moddate) %>% 
  group_by(id) %>% 
  mutate(current_art = last(treatment)) %>% 
  select(-(num_art:precision)) %>% 
  # Compute number of classes
  mutate(classes = drop_backbone(treatment),
         classes_n = str_count(classes, "\\S+")) 


# Treated critera: 
# - at least 2 classes, of which 1 has to be a boosted regime
# - switched to either DTG/F/TAF, BIG/F/TAF, EVG/c/F/TAF or DTG/ABC/3TC

treat <- art %>% 
  # Now, apply complicated filter
  mutate(flag = if_else(classes_n == 1 & 
                        lag(classes_n > 1) & 
                        treatment %in% c("DTG ETC TAF", 
                                         "BIC ETC TAF", 
                                         "3TC ABC DTG", 
                                         "COB ETC EVG TAF") & 
                        str_detect(lag(treatment), "RTV|COB"),
                        "X", "")) %>% 
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




# Control group inclusion = 
control <- art %>% 
  mutate(flag = if_else(classes_n > 1 & 
                        str_detect(treatment, "RTV|COB") & 
                        moddate == last(moddate), 
                        "X", "")) %>% 
  filter(any(flag == "X")) %>% 
  mutate(switch = FALSE) %>% 
  arrange(id, desc(moddate)) %>% 
  select(id, last_art = current_art, switch) %>% 
  slice(1)


# There are some individuals in the control group. This is because they simplified, 
# but later had an intensification of their treatment, and therefore the current
# ART regimen contains >= classes. They remain in the "treat" group. I therefore
# remove them from the controls

control <- control %>% 
  filter(!(id %in% treat$id))


# I don't want individuals who have DDI or D4T or so in their final regime, this
# is not current practice so I exlcude them here.

control <- control %>% 
  filter(!str_detect(last_art, "DDI|D4T"))

# Join them together
full <- treat %>% 
  full_join(control) %>% 
  mutate(last_art_s = drop_backbone(last_art)) %>% 
  relocate(last_art_s, .after = last_art)


full <- full %>% 
  left_join(recent %>% select(id, fupdate)) %>% 
  rename(last_fup = fupdate)


full %>% 
  # filter(last_fup >= dmy("01.01.2020")) %>% 
  count(comp_reg_s, simple_reg_s, sort = T) %>% 
  flextable() %>% 
  autofit()




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


# Calculate number of changes in database
n_treatments <- modif_art %>% 
  group_by(id) %>% 
  summarise(n_treatments = n())


# Join together
pop_det <- pat %>% 
  mutate(male = if_else(sex == 1, 1, 0), 
         ethnicity = case_when(ethnicity == 1 ~ "White", 
                               ethnicity == 2 ~ "Black",
                               ethnicity %in% c(0, 3, 4) ~ "Other", 
                               TRUE ~ NA_character_), 
         age = 2021 - born) %>%
  mutate(ethnicity = factor(ethnicity, levels = c("White", "Black", "Other"))) %>% 
  select(id, age, male, ethnicity) %>% 
  left_join(risk, by = "id") %>% 
  left_join(source, by = "id") %>% 
  left_join(n_treatments, by = "id")

full <- full %>% 
  left_join(pop_det, by = "id")




# Baseline Date -----------------------------------------------------------

# I assign a random baseline for every control patient


# There is an error in ID 51171, it says 8020 instead fo 2020, correct it.
full$switch_date[full$id == 51171] <- full$switch_date[full$id == 51171] - years(6000)

# Same for ID 42485, it says 2917 instead of 2017
full$switch_date[full$id == 42485] <- full$switch_date[full$id == 42485] - years(900)

switch_dates <- full %>% 
  filter(switch == TRUE) %>% 
  pull(switch_date)


n_non_switchers <- full %>% 
  filter(switch == FALSE) %>% 
  nrow()


# Select random samples of switch_dates
set.seed(123)
bl_dates_non_switchers <- sample(switch_dates, n_non_switchers, replace = TRUE)

# Assign them to baseline_date
full$baseline_date <- full$switch_date
full$baseline_date[full$switch == FALSE] <- bl_dates_non_switchers

full$baseline_date %>% summary()

# Looks good
full %>% 
  ggplot(aes(x = baseline_date)) + 
  geom_density(aes(fill = switch), alpha = 0.5, adjust = 0.5)



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
  group_by(id) %>% 
  mutate(no_third = (num_nnrti + num_pi + num_ntrti + 
                       num_fi + num_others + num_inti) == 0) %>% 
  summarise(nrti_mono = any(num_nrti <= 2 & no_third & !is.na(treatment)))



# Viral failure
art <- modif_art %>% 
  filter(id %in% full$id) %>% 
  select(-(num_art:precision))

rna <- lab %>% 
  select(id, rna, labdate) %>% 
  drop_na() %>% 
  filter(id %in% full$id)


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


setDT(rna);setDT(on_art)

comb_art <- on_art[rna, 
                   on = .(id, start <= labdate, stop >= labdate), 
                   .(id, rna = i.rna, 
                     rna_date = i.labdate, 
                     art_period = x.p, 
                     on_art = x.on_art, 
                     period_start = x.start, 
                     period_stop = x.stop)] %>% 
  as_tibble()

set.seed(1235)
smpl <- sample(unique(full$id), 10)

comb_art %>% 
  group_by(id) %>% 
  mutate(flag = if_else(rna > 500 & lag(rna) > 500 & on_art == 1, 
                        "X", "")) %>% 
  group_by(id, art_period) %>% 
  mutate(t = as.numeric(rna_date - first(period_start)), 
         failure = if_else(flag == "X" & t >= 180, "Y", "N")) %>% 
  filter(id %in% smpl) %>% 
  View()




comb_art %>% 
  group_by(id, rna_date, art_period) %>% 
  filter(id %in% smpl) %>% 
  ggplot(aes(x = rna_date, y = rna)) +
  geom_line(aes(group = id)) +
  geom_point(aes(color = factor(on_art))) +
  facet_wrap(~id) +
  scale_y_log10()


comb_art %>% 
  filter(id %in% smpl) %>% 
  print(n = 220)

rna %>% 
  filter(id %in% smpl) %>% 
  group_by(id) %>% 
  # Flag virological failure
  mutate(flag = if_else(rna > 500 & lag(rna) <= 500, "X", "")) %>% 
  View("rna")







# Tables   -----------------------------------------------------------------


# Define variables for table
vars <- c("age", "male", "ethnicity", "riskgroup", 
          "source", "center", "n_treatments")
cat_vars <- c("male", "ethnicity", "riskgroup", "source", "center")

# Create vector used later for indentation in flextable
cvars <- paste0(c(vars, "^n"), collapse = "|^")


# Table 1
tab1 <- CreateTableOne(vars = vars, strata = "switch", factorVars = cat_vars, 
               data = full)
  
t_patchar <- (print(tab1, nonnormal = c("age", "n_treatments"), 
                   printToggle = FALSE, contDigits = 1, dropEqual = TRUE) %>% 
  as_tibble(rownames = "Variable") %>% 
  select(-test) %>% 
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
  count(last_art_s, sort = T) %>% 
  rename("Last ART" = last_art_s) %>%
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
