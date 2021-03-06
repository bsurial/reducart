
# Setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(bernr)
library(lubridate)
library(nephro)
library(gtsummary)
library(gt)
library(here)
library(dtplyr)
set.seed(12)


# Load Datasets -----------------------------------------------------------

df <- pro_read("02-proposal_population.rds")
tail <- pro_read("tail.rds")
art <- pro_read("modif_art.rds")
lab <- pro_read("lab.rds")
lab2 <- pro_read("lab2.rds")
brand_dose <- pro_read("brand_dose.rds")
adhe <- pro_read("adhe.rds")
pat <- pro_read("pat.rds")
fup <- pro_read("fup.rds")


# Add contemporary treatments ---------------------------------------------

# TREATMENTS
insti <- c("DTG|RGV|EVG|BIC")
pi <- c("DRV|ATV|LPV|SQV|IDV|FAPV|NFV|TPV")
nnrti <- c("EFV|NVP|RPV|ETV|DOR")
booster <- c("RTV|COB")


# FUNCTIONS

# remove backbones from treatment string
drop_backbone <- function(string) {
  string <- str_remove_all(string, "ABC|3TC|ETC|TAF|TDF|COB|RTV|AZT|DDI|D4T")
  str_squish(string)
}


# THEMES
# Set GT Summary Theme
my_theme <-
  list(
    "tbl_summary-str:default_con_type" = "continuous",
    "tbl_summary-str:continuous_stat" = "{median} ({p25} to {p75})",
    "style_number-arg:big.mark" = "",
    "tbl_summary-arg:missing_text" = "(Missing)"
  )
set_gtsummary_theme(my_theme)


# eligible patients -------------------------------------------------------

# FUP after Jan 2013
last_fupdate <- fup %>% 
  select(id, fupdate) %>% 
  arrange(id, desc(fupdate)) %>% 
  group_by(id) %>% 
  summarise(last_fup_date = first(fupdate))


last_fup_sub <- last_fupdate %>% 
  filter(last_fup_date >= dmy("01.01.2013"))


# Ever had booster + 2 of NNRTI/INSTI/PI
df <- last_fup_sub %>% 
  left_join(art, by = "id") %>% 
  select(id, last_fup_date, treatment, moddate, enddate) %>% 
  mutate(booster = str_detect(treatment, booster),
         insti = str_detect(treatment, insti),
         nnrti = str_detect(treatment, nnrti), 
         pi = str_detect(treatment, pi)) %>%
  group_by(id) %>% 
  filter(any(booster == TRUE & (insti + nnrti + pi) >1)) %>% 
  ungroup() %>% 
  distinct(id, last_fup_date)





# Create Trials -----------------------------------------------------------


# Set-up the grid for the trials. One trial starts every month since 1. June
# 2013 (introduction of DTG in Switzerland)
trial_grid <- tibble(trial_start = seq(ymd('2013-06-01'),
                                       ymd('2020-01-01'), 
                                       by = '1 month')) %>% 
  mutate(trial_nr = row_number()) %>% 
  relocate(trial_nr)

# Add id to grid
trial_grid <- trial_grid %>% 
  crossing(id = df$id) %>% 
  arrange(id, trial_nr)




# Add info on treatment ---------------------------------------------------

df_sub <- df %>% 
  # select(-moddate) %>% 
  left_join(art %>% select(-(num_art:precision)))



# Since the last treatment is ongoin, I create a fictive "enddate in future"
df_sub <- df_sub %>% 
  mutate(enddate2 = if_else(is.na(enddate), dmy("01.01.2100"), enddate))

# Join current treatment to grid
setDT(trial_grid); setDT(df_sub)
df_step1 <- df_sub[trial_grid, 
           on = .(id, moddate <= trial_start, enddate2 >= trial_start),
           .(id, trial_start =  i.trial_start, trial_nr = i.trial_nr, 
             treatment = x.treatment)] %>% 
  as_tibble()


# Add HIV RNA -------------------------------------------------------------


# Now I add Info on HIV RNA. Eligible RNA values will be within 1 year prior
# to the start of the trial, and also 14 days after
df_step1 <- df_step1 %>% 
  mutate(rna_start = trial_start - 365, 
         rna_stop = trial_start + 14)


rna <- lab %>% 
  select(id, rna, labdate) %>% 
  drop_na()


rna_sub <- rna 



# Here I calculate the time since viral suppression for each lab visit
# See also https://stackoverflow.com/questions/69941208/count-time-since-last-condition-met-and-reset-to-0-when-not

time_supp_df <- rna_sub %>% 
  as_tibble() %>% 
  lazy_dt() %>% 
  # I allow blips below 200
  mutate(blip = if_else(rna >= 50 & rna < 200 & lag(rna) < 50 & lead(rna) < 50, 
                        TRUE, FALSE)) %>% 
  mutate(suppressed = rna < 50 | blip == TRUE) %>% 
  group_by(id) %>% 
  mutate(time_suppressed = case_when(
    labdate == first(labdate) ~ labdate - labdate,
    suppressed == FALSE ~ labdate - labdate, 
    suppressed == TRUE & lag(suppressed) == FALSE ~ labdate - labdate,
    suppressed == TRUE & lag(suppressed) == TRUE ~ labdate - lag(labdate)),
    time_suppressed = as.numeric(time_suppressed)) %>% 
  mutate(grp = cumsum(time_suppressed == 0)) %>% 
  group_by(id, grp) %>% 
  mutate(time_suppressed = cumsum(time_suppressed)) %>% 
  ungroup() %>% 
  select(-grp) %>% 
  as_tibble()


# Add time since viral suppression to dataset
setDT(df_step1); setDT(time_supp_df)

df_step2 <- time_supp_df[df_step1, 
        on = .(id, labdate >= rna_start, labdate <= rna_stop), 
        .(id, trial_start = i.trial_start, trial_nr = i.trial_nr, 
          treatment = i.treatment, rna = x.rna, rna_date = x.labdate, 
          days_suppressed = x.time_suppressed)] %>% 
  as_tibble() %>% 
  arrange(id, trial_nr, desc(rna_date)) %>% 
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  
  # KEY Assumption: If an id is suppressed at a given time point "A" for 30 days, 
  # then at the trial that start 10 days later, the id is suppressed for 40 days
  mutate(days_suppressed_actual = if_else(days_suppressed == 0, 0, 
                                          days_suppressed + 
                                            as.numeric(trial_start - rna_date)))


# Add contraindicated drugs -----------------------------------------------

# Add drugs
# Rifampicin
rifa_df <- brand_dose %>% 
  filter(str_detect(brand, "(J04AM0[2567])|(J04AB02)")) %>% 
  mutate(rifamp = TRUE) %>% 
  select(id, startdate, stopdate, rifamp, rifa_drug = var_desc) %>% 
  # stopdate far in future for join later
  mutate(stopdate = if_else(is.na(stopdate), dmy("01.01.2100"), stopdate))

# Carboxamide derivatives (Carbamazepine, Oxcarbazepine etc)
carbam_df <- brand_dose %>% 
  filter(str_detect(brand, "N03AF0[1234]")) %>% 
  mutate(carbam = TRUE) %>% 
  select(id, startdate, stopdate, carbam, carbam_drug = var_desc) %>% 
  # stopdate far in future for join later
  mutate(stopdate = if_else(is.na(stopdate), dmy("01.01.2100"), stopdate))

# Phenytoin and other hydantoin derivatives
phenytoin_df <- brand_dose %>% 
  filter(str_detect(brand, "N03AB0[12345]|N03AB5[24]")) %>% 
  mutate(phenytoin = TRUE) %>% 
  select(id, startdate, stopdate, phenytoin, phenytoin_drug = var_desc) %>% 
  # stopdate far in future for join later
  mutate(stopdate = if_else(is.na(stopdate), dmy("01.01.2100"), stopdate))

# Primidone 
primidone_df <- brand_dose %>% 
  filter(str_detect(brand, "N03AA030")) %>% 
  mutate(primidone = TRUE) %>% 
  select(id, startdate, stopdate, primidone, primidone_drug = var_desc) %>% 
  # stopdate far in future for join later
  mutate(stopdate = if_else(is.na(stopdate), dmy("01.01.2100"), stopdate))





setDT(df_step2) 
setDT(rifa_df); setDT(carbam_df); setDT(phenytoin_df); setDT(primidone_df) 

df_step3 <- rifa_df[df_step2, 
                         on = .(id, startdate <= trial_start, 
                                stopdate >= trial_start), 
                         .(id, trial_start = i.trial_start, trial_nr = i.trial_nr, 
                           treatment = i.treatment, rna = i.rna, 
                           rna_date = i.rna_date, 
                           days_suppressed = i.days_suppressed, 
                           days_suppressed_actual = i.days_suppressed_actual, 
                           rifamp = x.rifamp, rifa_drug = x.rifa_drug)] %>% 
  lazy_dt() %>% 
  mutate(rifamp = replace_na(rifamp, FALSE)) %>% 
  arrange(id, trial_nr, desc(rifamp)) %>% # This places any TRUE before FALSE
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  as_tibble()


setDT(df_step3) 
df_step4 <- carbam_df[df_step3, 
        on = .(id, startdate <= trial_start, 
               stopdate >= trial_start), 
        .(id, trial_start = i.trial_start, trial_nr = i.trial_nr, 
          treatment = i.treatment, rna = i.rna, 
          rna_date = i.rna_date, 
          days_suppressed = i.days_suppressed, 
          days_suppressed_actual = i.days_suppressed_actual, 
          rifamp = i.rifamp, rifa_drug = i.rifa_drug, 
          carbam = x.carbam, carbam_drug = x.carbam_drug)] %>% 
  lazy_dt() %>% 
  mutate(carbam = replace_na(carbam, FALSE)) %>% 
  arrange(id, trial_nr, desc(carbam)) %>% # This places any TRUE before FALSE
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  as_tibble()

setDT(df_step4)
df_step5 <- phenytoin_df[df_step4, 
                      on = .(id, startdate <= trial_start, 
                             stopdate >= trial_start), 
                      .(id, trial_start = i.trial_start, trial_nr = i.trial_nr, 
                        treatment = i.treatment, rna = i.rna, 
                        rna_date = i.rna_date, 
                        days_suppressed = i.days_suppressed, 
                        days_suppressed_actual = i.days_suppressed_actual, 
                        rifamp = i.rifamp, rifa_drug = i.rifa_drug, 
                        carbam = i.carbam, carbam_drug = i.carbam_drug, 
                        phenytoin = x.phenytoin, 
                        phenytoin_drug = x.phenytoin_drug)] %>% 
  lazy_dt() %>%
  mutate(phenytoin = replace_na(phenytoin, FALSE)) %>% 
  arrange(id, trial_nr, desc(phenytoin)) %>% # This places any TRUE before FALSE
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  as_tibble()

setDT(df_step5)
df_step6 <- primidone_df[df_step5, 
             on = .(id, startdate <= trial_start, 
                    stopdate >= trial_start), 
             .(id, trial_start = i.trial_start, trial_nr = i.trial_nr, 
               treatment = i.treatment, rna = i.rna, 
               rna_date = i.rna_date, 
               days_suppressed = i.days_suppressed, 
               days_suppressed_actual = i.days_suppressed_actual, 
               rifamp = i.rifamp, rifa_drug = i.rifa_drug, 
               carbam = i.carbam, carbam_drug = i.carbam_drug, 
               phenytoin = i.phenytoin, 
               phenytoin_drug = i.phenytoin_drug, 
               primidone = x.primidone, 
               primidone_drug = x.primidone_drug
               )] %>% 
  lazy_dt() %>%
  mutate(primidone = replace_na(primidone, FALSE)) %>% 
  arrange(id, trial_nr, desc(primidone)) %>% # This places any TRUE before FALSE
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  as_tibble()


df_step6 <- df_step6 %>% 
  mutate(ci_drug = (rifamp + carbam + phenytoin + primidone) > 0)


# Add adherence data ------------------------------------------------------

adhe_df <- adhe %>% 
  mutate(missed = if_else(missed == "Z", NA_character_, missed),
         missed = factor(missed, 
                         labels = c("every day", 
                                    "more than 1/week",
                                    "once a week", 
                                    "once every 2 weeks", 
                                    "once a month", 
                                    "never"))) %>% 
  arrange(id, ad_date) %>% 
  # I fill missing data using LOCF
  group_by(id) %>% 
  mutate(missed_locf = missed) %>% 
  fill(missed_locf, .direction = "downup") %>% 
  ungroup()

# Calculate Time with adherence once every 2 weeks or better
adhe_df <- adhe_df %>% 
  # filter(id == 12601) %>% 
  lazy_dt() %>% 
  mutate(good_adh = missed_locf %in% 
           c("once a month", "never", "once every 2 weeks")) %>% 
  group_by(id) %>% 
  mutate(time_good_adh = case_when(
    good_adh == FALSE ~ ad_date - ad_date, 
    good_adh == TRUE & lag(good_adh) == FALSE ~ ad_date -  ad_date,
    row_number() == 1 ~ ad_date - ad_date,
    good_adh = TRUE & lag(good_adh) == TRUE ~ ad_date - lag(ad_date)),
    time_good_adh = as.numeric(time_good_adh)) %>% 
  mutate(grp = cumsum(time_good_adh == 0)) %>% 
  group_by(id, grp) %>% 
  mutate(time_good_adh = cumsum(time_good_adh)) %>% 
  ungroup() %>% 
  as_tibble()

# Set consider all adherence data one year prior to and 14 days after the trial
# start day

df_step6 <- df_step6 %>% 
  mutate(adhe_low = trial_start - 365.25, 
         adhe_high = trial_start + 14)



setDT(df_step6); setDT(adhe_df)
df_step7 <- adhe_df[df_step6, 
                    on = .(id, ad_date >= adhe_low, 
                           ad_date <= adhe_high), 
                    .(id, trial_start = i.trial_start, trial_nr = i.trial_nr, 
                      treatment = i.treatment, rna = i.rna, 
                      rna_date = i.rna_date, 
                      days_suppressed = i.days_suppressed, 
                      days_suppressed_actual = i.days_suppressed_actual, 
                      rifamp = i.rifamp, rifa_drug = i.rifa_drug, 
                      carbam = i.carbam, carbam_drug = i.carbam_drug, 
                      phenytoin = i.phenytoin, 
                      phenytoin_drug = i.phenytoin_drug, 
                      primidone = i.primidone, 
                      primidone_drug = i.primidone_drug, 
                      ci_drug = i.ci_drug,
                      good_adh = x.good_adh,
                      time_good_adh = x.time_good_adh, 
                      ad_date = x.ad_date, 
                      adherence = x.missed,
                      adherence_locf = x.missed_locf
                    )] %>% 
  lazy_dt() %>% 
  arrange(id, trial_nr, abs(trial_start - ad_date)) %>% 
  group_by(id, trial_nr) %>% 
  slice(1) %>%
  ungroup() %>% 
  mutate(time_good_adh_actual = if_else(good_adh == FALSE & 
                                          time_good_adh == 0, 0, 
                                        time_good_adh + 
                                          as.numeric(trial_start - ad_date)),
         time_good_adh_actual = if_else(time_good_adh_actual < 0, 0, 
                                        time_good_adh_actual)) %>% 
  as_tibble()




# Add eGFR ----------------------------------------------------------------

crea_df <- lab2 %>% 
  filter(item == "CRE") %>% 
  select(id, crea = value, cre_date = labdate) %>% 
  drop_na() %>% 
  left_join(pat %>% select(id, born, sex, ethnicity), by = "id") %>% 
  mutate(male = as.numeric(sex == 1)) %>% 
  mutate(egfr = CKDEpi.creat(creatinine = crea / 88.42,
                             sex = male, 
                             age = year(cre_date) - born, 
                             ethnicity = ethnicity == 2)) %>% 
  select(-(sex:male))

# Window for Creatinine
df_step7 <- df_step7 %>% 
  mutate(crea_hi = trial_start + 14, 
         crea_lo = trial_start - 365.25)

setDT(df_step7); setDT(crea_df)

df_step8 <- crea_df[df_step7, 
        on = .(id, cre_date >= crea_lo, cre_date <= crea_hi),
        .(id, trial_start = i.trial_start, trial_nr = i.trial_nr, 
          treatment = i.treatment, rna = i.rna, 
          rna_date = i.rna_date, 
          days_suppressed = i.days_suppressed, 
          days_suppressed_actual = i.days_suppressed_actual, 
          rifamp = i.rifamp, rifa_drug = i.rifa_drug, 
          carbam = i.carbam, carbam_drug = i.carbam_drug, 
          phenytoin = i.phenytoin, 
          phenytoin_drug = i.phenytoin_drug, 
          primidone = i.primidone, 
          primidone_drug = i.primidone_drug, 
          ci_drug = i.ci_drug,
          time_good_adh = i.time_good_adh, 
          ad_date = i.ad_date, 
          adherence = i.adherence,
          adherence_locf = i.adherence_locf, 
          time_good_adh_actual = i.time_good_adh_actual, 
          crea = x.crea, 
          egfr = x.egfr, 
          cre_date = x.cre_date, 
          born = x.born)] %>% 
  lazy_dt() %>% 
  arrange(id, trial_nr, abs(trial_start - cre_date)) %>% 
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  as_tibble()




# Add patient characteristics -----------------------------------------------

patchars <- pat %>% 
  select(id, sex, ethnicity) %>% 
  left_join(tail %>% select(id, riskgroup)) %>% 
  mutate(riskgroup = case_when(riskgroup %in% c("BLOOD", "OTHER", "PERINAT", 
                                                "UNKNOWN") ~ "OTHER", 
                               TRUE ~ riskgroup)) %>% 
  mutate(female = as.numeric(sex == 2),
         ethnicity = factor(ethnicity, labels = c("Other", "White", 
                                                  "Black", "Other", 
                                                  "Other", "Unknown"))) %>% 
  select(id, female, ethnicity, riskgroup)


df_step9 <- df_step8 %>% 
  mutate(age = year(trial_start) - born) %>% 
  left_join(patchars, by = "id")





# See whether patient has eligible treatment ------------------------------

# Eligible means: 
# Booster + 2 or more drugs of either NNRTI, PI or INST = current regimen
# BIC or DTG based standrard therapy and just before current regimen = switch

df_step10 <- df_step9 %>% 
  mutate(booster = str_detect(treatment, booster),
         insti = str_detect(treatment, insti),
         nnrti = str_detect(treatment, nnrti), 
         pi = str_detect(treatment, pi)) %>% 
  mutate(elig_treatment = if_else(booster == TRUE & (insti + nnrti + pi) >1, 
                                  1, 0), 
         switch_treatment = if_else(booster == FALSE & 
                                    str_detect(treatment, "BIC|DTG") & 
                                    nnrti == FALSE & 
                                    pi == FALSE & 
                                    lag(elig_treatment == 1), 1, 0)) %>% 
  relocate(treatment, .before = booster)



# Add resistance data -----------------------------------------------------

resi_data <- read_csv(here("data", "2201_resistance", "2201_drm.csv")) %>% 
  select(-1) 

res_cum <- resi_data %>% 
  select(id = PATIENT, dat, ends_with("_c"))

df_step10_res <- df_step10 %>% 
  select(id, trial_start) %>% 
  left_join(res_cum) %>% 
  lazy_dt() %>% 
  filter(dat < trial_start | is.na(dat)) %>%
  arrange(id, trial_start, desc(dat)) %>% 
  group_by(id, trial_start) %>% 
  slice(1) %>% 
  ungroup() %>% 
  as_tibble()


tams <- "(M41[A-Z]?[A-Z]?)|(D67[A-Z]?[A-Z]?)|(K70[R])|(L210[A-Z]?[A-Z]?)|(T215[A-Z]?[A-Z]?)|(K219[A-Z]?[A-Z]?)"


df_step10_nrti <- df_step10 %>% 
  left_join(df_step10_res) %>% 
  select(id, trial_start, dat, NRTI_c) %>% 
  mutate(TAMS = str_extract_all(NRTI_c, tams)) %>% 
  mutate(TAMS_n = map(TAMS, ~length(.x)), 
         TAMS = map(TAMS, ~paste0(.x, collapse = "|"))) %>% 
  unnest(c(TAMS, TAMS_n)) %>% 
  mutate(K65R = str_detect(NRTI_c, "K65R"), 
         M184VI = str_detect(NRTI_c, "M184V?I?"), 
         T69_ins = str_detect(NRTI_c, "T69Insertion"), 
         Q151M = str_detect(NRTI_c, "Q151M[A-Z]?|Q151QKLM"))


df_step10_insti <- df_step10 %>% 
  left_join(df_step10_res) %>% 
  select(id, trial_start, dat, starts_with("INSTI.")) %>% 
  mutate(major_insti = str_detect(INSTI.Major_c, "Q148|R263|G118"), 
         n_insti = str_count(str_replace(INSTI.Major_c, ",", " "), "\\w+"), 
         major_insti = if_else(n_insti > 1, TRUE, major_insti))


df_step11 <- df_step10 %>% 
  left_join(df_step10_nrti) %>% 
  left_join(df_step10_insti) %>% 
  left_join(df_step10_res %>% select(-INSTI.Accessory_c))





# Assess eligibility for each trial ---------------------------------------

elig_data <- df_step11 %>% 
  mutate(elig_current = if_else(elig_treatment == 1 & 
                          days_suppressed_actual >= 180 & 
                          egfr >= 30 & 
                          time_good_adh_actual >= 180 &
                          ci_drug == FALSE & 
                          (major_insti == FALSE | is.na(major_insti)) &
                          (Q151M == FALSE | is.na(Q151M)), 
                        1, 0), 
         elig_switch = if_else(switch_treatment == 1 & 
                                 days_suppressed_actual >= 180 & 
                                 egfr >= 30 & 
                                 time_good_adh_actual >= 180 &
                                 ci_drug == FALSE &
                                 (major_insti == FALSE | is.na(major_insti)) &
                                 (Q151M == FALSE | is.na(Q151M)),
                               1, 0))






# Add history of VF and NRTI Mono ---------------------------------------------

### NEED TO INCLUDE IT UPSTREAM


## NRTI Mono one or two NRTIs, not 3
nrti_mono <- art %>% 
  group_by(id) %>% 
  mutate(no_third = (num_nnrti + num_pi + num_ntrti + 
                       num_fi + num_others + num_inti) == 0) %>% 
  mutate(nrti_mono = num_nrti <= 2 & 
           no_third & 
           !is.na(treatment)) %>% 
  select(id, treatment, moddate, enddate, no_third, nrti_mono) %>% 
  ungroup()


elig_data <- elig_data %>% 
  select(id, trial_start) %>% 
  left_join(nrti_mono, by = "id") %>% 
  group_by(id, trial_start) %>% 
  # Any use of NRTI_mono prior to trial start
  summarise(nrti_mono = any(nrti_mono & moddate < trial_start)) %>% 
  ungroup() %>% 
  right_join(elig_data, by = c("id", "trial_start")) %>% 
  to_last(nrti_mono)
  
  

# Add Virological failure:
# RNA > 200 after 180 days of ART or 2 x >200 after previously suppressed


on_treatment <- art %>% 
  select(-(num_art:precision)) %>% 
  group_by(id) %>% 
  mutate(on_art = as.numeric(!is.na(treatment))) %>% 
  mutate(p = cumsum(on_art != lag(on_art, default = 2))) %>% 
  group_by(id, p) %>% 
  summarise(start = first(moddate), 
            stop = last(enddate), 
            on_art = max(on_art),
            .groups = "drop") %>% 
  mutate(stop = if_else(is.na(stop), dmy("01.01.2100"), stop))


setDT(rna);setDT(on_treatment)
rna_art <- on_treatment[rna, 
                        on = .(id, start <= labdate, stop > labdate), 
                        .(id, rna = i.rna, 
                          rna_date = i.labdate, 
                          art_period = x.p, 
                          on_art = x.on_art, 
                          period_start = x.start, 
                          period_stop = x.stop)] %>% 
  as_tibble() %>% 
  mutate(days_on_art = if_else(on_art == 1, 
                               as.numeric(rna_date - period_start), 
                               0))




# Adherence data
ad_table <- adhe_df %>% 
  as_tibble() %>% 
  arrange(id, ad_date) %>% 
  select(-(amenddate:inputdate)) %>% 
  group_by(id) %>% 
  mutate(end_period = lead(ad_date)) %>% 
  ungroup() %>% 
  mutate(ad_date_windowed = ad_date - 21, 
         end_period_windowed = end_period - 21) %>% 
  mutate(end_period = replace_na(end_period, dmy("01.01.2100"))) %>%  
  select(id, missed_locf, missed, in_row, ad_date, end_period, 
         ad_date_windowed, end_period_windowed)


setDT(rna_art);setDT(ad_table)

rna_detailed_long <- ad_table[rna_art, 
                              on = .(id, ad_date_windowed <= rna_date, 
                                     end_period_windowed > rna_date), 
                              .(id, 
                                rna = i.rna, 
                                rna_date = i.rna_date, 
                                art_period = i.art_period, 
                                on_art = i.on_art, 
                                period_start = i.period_start, 
                                period_stop = i.period_stop, 
                                missed = x.missed_locf, 
                                in_row = x.in_row, 
                                ad_date = x.ad_date)] %>% 
  as_tibble()




failure_df <- rna_detailed_long %>% 
  lazy_dt() %>% 
  group_by(id, art_period) %>% 
  mutate(days_on_art = if_else(on_art == 1, 
                               as.numeric(rna_date - first(period_start)), 
                               0)) %>% 
  # add art_period so if we start a new ART, previoulsy suppressed resets to 0
  group_by(id, art_period) %>% 
  mutate(time_lapsed = as.numeric(rna_date - lag(rna_date))) %>% 
  mutate(previously_suppressed = cummax(rna < 200)) %>% 
  mutate(failure = case_when(
    rna >= 200 & days_on_art > 180  & previously_suppressed == 0 ~ 1, 
    rna >= 200 & lag(rna) >= 200 & 
      time_lapsed >= 14 & previously_suppressed == 1 ~ 1, 
    TRUE ~ 0
  )) %>% 
  as_tibble()


elig_data <- elig_data %>% 
  lazy_dt() %>% 
  select(id, trial_start) %>% 
  left_join(failure_df, by = "id") %>% 
  group_by(id, trial_start) %>% 
  summarise(history_VF = any(failure == 1 & rna_date < trial_start)) %>% 
  ungroup() %>% 
  right_join(elig_data, by = c("id", "trial_start")) %>% 
  to_last(history_VF) %>% 
  as_tibble()



# Add covariates ----------------------------------------------------------

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


elig_data <- elig_data %>% 
  left_join(max_rna, by = "id") %>% 
  left_join(nadir_cd4, by = "id") %>% 
  mutate(max_rna_log = log(max_rna))


# Add n_conmeds -----------------------------------------------------------

drug_small <- brand_dose %>% 
  select(id, brand, var_desc, startdate, stopdate) %>% 
  mutate(stopdate = if_else(is.na(stopdate), dmy("01.01.2100"), 
                            stopdate))

# Exclude HIV drugs and vaccines from drug_small vvv
HIV_drugs <- c("J05AR", "J05AE", "J05AX", "J05AG",
               "J05AF", "J05AJ", "_trialproduct_J05A_", 
               "J07")
HIV_drugs_pattern <- paste0(HIV_drugs, collapse = "|")
drug_nohiv <- drug_small %>% 
  filter(!str_detect(brand, HIV_drugs_pattern))

# Join them vvv
elig_small <- elig_data %>% 
  select(id, trial_start)

setDT(drug_nohiv); setDT(elig_small)

conmeds <- drug_nohiv[elig_small, 
                      on = .(id, startdate <= trial_start, 
                             stopdate >= trial_start), 
                      .(id, 
                        trial_start = i.trial_start,
                        brand = x.brand)] %>% 
  lazy_dt() %>% 
  distinct(id, trial_start, brand) %>% # remove duplicates per trial
  group_by(id, trial_start) %>% 
  summarise(n_conmeds = sum(!is.na(brand))) %>% 
  arrange(desc(n_conmeds)) %>% 
  as_tibble()


elig_data <- elig_data %>% 
  left_join(conmeds, by = c("id", "trial_start"))






# Summary -----------------------------------------------------------------


# Unique Patients
(t_switch_reg <- elig_data %>% 
  filter(elig_switch == 1) %>% 
  mutate(treatment = fct_infreq(treatment)) %>% 
  select(treatment) %>% 
  tbl_summary(label = treatment ~ "Switch Regimen") %>% 
  bold_labels())

gtsave(as_gt(t_switch_reg), here::here("tables", "04-switch_regimens.png"))


elig_data %>% 
  filter(elig_current == 1) %>% 
  distinct_ids()


(t_current_reg <- elig_data %>% 
  filter(elig_current == 1) %>% 
  mutate(treatment_lumped = fct_infreq(fct_lump(treatment, n = 15))) %>% 
  select(treatment_lumped) %>% 
  tbl_summary(label = treatment_lumped ~ "Current Regimen") %>% 
  bold_labels())

gtsave(as_gt(t_current_reg), here::here("tables", "04-current_regimens.png"))



switchers_1 <- elig_data %>% 
  filter(elig_switch == 1) %>% 
  distinct(id, female, ethnicity, riskgroup, max_rna_log, nadir_cd4) %>% 
  mutate(group = "switched")

current_1 <- elig_data %>% 
  filter(elig_current == 1) %>% 
  distinct(id, female, ethnicity, riskgroup, max_rna_log, nadir_cd4) %>% 
  mutate(group = "current")


all_1 <- bind_rows(switchers_1, current_1)


(t_1 <- all_1 %>% 
    select(-id) %>% 
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
    add_p() %>% 
    bold_labels())

gtsave(as_gt(t_1), here::here("tables", "04-patient_characteristics.png"))



# Trials

studypop_filtered <- elig_data %>% 
  filter(elig_current == 1 | elig_switch == 1) %>% 
  mutate(group = if_else(elig_current == 1, "current", 
                         if_else(elig_switch == 1, "switch", NA_character_)))


(t_2 <- studypop_filtered %>%
    mutate(t_suppressed_actual = days_suppressed_actual / 365.25, 
           t_good_adh = time_good_adh_actual / 365.25) %>% 
    select(group, female, ethnicity, age, riskgroup, max_rna_log, 
           nadir_cd4, t_suppressed_actual, t_good_adh, adherence_locf, 
           egfr, history_VF, nrti_mono, n_conmeds) %>% 
    labelled::set_variable_labels(
      female = "Female sex",
      ethnicity = "Ethnicity",
      age = "Age at trial start, years (IQR)",
      riskgroup = "HIV transmission group",
      max_rna_log = "Pretreatment HIV viral load, log cp/mL (IQR)",
      nadir_cd4 = "Nadir CD4 count, cells/µL (IQR)",
      t_suppressed_actual = "Time with HIV <50 cp/mL, years (IQR)",
      t_good_adh = "Time with good adherence, years (IQR)",
      adherence_locf = "Missed dose in the last 4 weeks",
      egfr = "eGFR at trial start, ml/min (IQR)",
      history_VF = "History of virological failure", 
      nrti_mono = "History of single/dual NRTI therapy",
      n_conmeds = "Nr. of non-HIV drugs"
    ) %>%
    tbl_summary(by = group, 
                digits = all_continuous() ~ 1) %>% 
    add_p() %>% 
    bold_labels())

gtsave(as_gt(t_2), here::here("tables", "04-patient_characteristics_alltrials.png"))


# Number of trials eligible per patient
elig_data %>% 
  filter(elig_current == 1) %>% 
  group_by(id) %>% 
  count() %>% 
  ungroup() %>% 
  summarise(median = median(n), 
            p25 = quantile(n, 0.25),
            p75 = quantile(n, 0.75))





# Write Data --------------------------------------------------------------

write_rds(elig_data, here::here("processed", "04-nested_trial_data.rds"))



