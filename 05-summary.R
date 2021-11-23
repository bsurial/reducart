library(tidyverse)
library(bernr)
library(gtsummary)
library(gt)
library(here)

# Themes ------------------------------------------------------------------

# Set GT Summary Theme
my_theme <-
  list(
    "tbl_summary-str:default_con_type" = "continuous",
    "tbl_summary-str:continuous_stat" = "{median} ({p25} to {p75})",
    "style_number-arg:big.mark" = "",
    "tbl_summary-arg:missing_text" = "(Missing)"
  )
set_gtsummary_theme(my_theme)


studypop <- pro_read("04-nested_trial_data.rds")
lab <- pro_read("lab.rds")

switchers_1 <- studypop %>% 
  filter(elig_switch == 1) %>% 
  distinct(id, riskgroup, female, ethnicity) %>% 
  mutate(group = "switched")

current_1 <- studypop %>% 
  filter(elig_current == 1) %>% 
  distinct(id, riskgroup, female, ethnicity) %>% 
  mutate(group = "current")


all_1 <- bind_rows(switchers_1, current_1)
  

max_rna <- lab %>% 
  select(id, rna, labdate) %>% 
  drop_na() %>% 
  group_by(id) %>% 
  summarise(max_rna = max(rna))

nadir_cd4 <- lab %>% 
  select(id, cd4, cd4date) %>% 
  drop_na() %>% 
  group_by(id) %>% 
  summarise(nadir_cd4 = min(cd4))

all_1 <- all_1 %>% 
  left_join(max_rna) %>% 
  left_join(nadir_cd4) %>% 
  mutate(max_rna_log = log(max_rna))


(t_1 <- all_1 %>% 
    select(-max_rna) %>% 
    labelled::set_variable_labels(
      female = "Female sex",
      ethnicity = "Ethnicity",
      riskgroup = "HIV transmission group",
      max_rna_log = "Pretreatment HIV viral load, log cp/mL (IQR)",
      nadir_cd4 = "Nadir CD4 count, cells/µL (IQR)"
    ) %>% 
    tbl_summary(by = group, 
                missing = "no",
                digits = list(all_continuous() ~ 1, 
                              nadir_cd4 ~ 0)) %>% 
    bold_labels())
  
gtsave(as_gt(t_1), here::here("tables", "05-patient_characteristics.png"))





# All trial participants --------------------------------------------------

studypop_filtered <- studypop %>% 
  filter(elig_current == 1 | elig_switch == 1) %>% 
  mutate(group = if_else(elig_current == 1, "current", 
                         if_else(elig_switch == 1, "switch", NA_character_)))


(t_2 <- studypop_filtered %>% 
  select(id, age, female, ethnicity, riskgroup, egfr,
         days_suppressed_actual, 
         ci_drug, time_good_adh_actual, group) %>% 
  left_join(max_rna) %>% 
  left_join(nadir_cd4) %>% 
  mutate(max_rna_log = log(max_rna)) %>% 
  select(-id, -max_rna) %>% 
  mutate(t_suppressed_actual = days_suppressed_actual / 365.25, 
         t_good_adh = time_good_adh_actual / 365.25) %>% 
  select(-days_suppressed_actual, -time_good_adh_actual) %>% 
  labelled::set_variable_labels(
    female = "Female sex",
    ethnicity = "Ethnicity",
    age = "Age at trial start, years (IQR)",
    riskgroup = "HIV transmission group",
    max_rna_log = "Pretreatment HIV viral load, log cp/mL (IQR)",
    nadir_cd4 = "Nadir CD4 count, cells/µL (IQR)",
    t_suppressed_actual = "Time with HIV <50 cp/mL, years (IQR)",
    t_good_adh = "Time with good adherence, years (IQR)",
    egfr = "eGFR at trial start, ml/min (IQR)",
    ci_drug = "Receipt of contra-indicated drug"
  ) %>%
  tbl_summary(by = group))
  
gtsave(as_gt(t_2), here::here("tables", "05-patient_characteristics_alltrials.png"))


# How many trials on average
studypop %>% 
  filter(elig_current == 1) %>% 
  group_by(id) %>% 
  count() %>% 
  ungroup() %>% 
  summarise(median = median(n), 
            p25 = quantile(n, 0.25),
            p75 = quantile(n, 0.75))



# Write id list for Roger
studypop %>% 
  distinct(id) %>% 
  write_csv(here("processed", "05-ids_for_resistance.csv"))

