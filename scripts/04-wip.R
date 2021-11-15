
# Setup -------------------------------------------------------------------


library(data.table)
library(tidyverse)
library(bernr)
library(lubridate)
library(nephro)

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
  mutate(id = 31462)



# Add info on treatment ---------------------------------------------------

df_art <- df %>% 
  select(-moddate) %>% 
  left_join(art %>% select(-(num_art:precision)))


# I first work on 1 patient, and extend that later on
df_sub <- df_art %>% 
  filter(id == 31462) %>% 
  select(id, treatment, moddate, enddate, switch) 

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

# Only work with one ID
rna_sub <- rna %>% 
  filter(id == 31462) 



# Here I calculate the time since viral suppression for each lab visit
# See also https://stackoverflow.com/questions/69941208/count-time-since-last-condition-met-and-reset-to-0-when-not

time_supp_df <- rna_sub %>% 
  as_tibble() %>% 
  # I allow blips below 200
  mutate(blip = if_else(rna >= 50 & rna < 200 & lag(rna) < 50 & lead(rna) < 50, 
                        TRUE, FALSE)) %>% 
  mutate(suppressed = rna < 50 | blip == TRUE) %>% 
  group_by(id) %>% 
  mutate(time_suppressed = case_when(
    suppressed == FALSE ~ labdate - labdate, 
    suppressed == TRUE & lag(suppressed) == FALSE ~ labdate - labdate,
    suppressed == TRUE & lag(suppressed) == TRUE ~ labdate - lag(labdate)),
    time_suppressed = as.numeric(time_suppressed)) %>% 
  group_by(id, grp = cumsum(time_suppressed == 0)) %>% 
  mutate(time_suppressed = cumsum(time_suppressed)) %>% 
  ungroup() %>% 
  select(-grp) 


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
  as_tibble() %>%
  mutate(across(rifamp:rifa_drug, ~replace_na(.x, FALSE))) %>% 
  arrange(id, trial_nr, desc(rifamp)) %>% # This places any TRUE before FALSE
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup()


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
          carbam = x.carbam, carbam_drug = x.carbam)] %>% 
  as_tibble() %>%
  mutate(across(carbam:carbam_drug, ~replace_na(.x, FALSE))) %>% 
  arrange(id, trial_nr, desc(carbam)) %>% # This places any TRUE before FALSE
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup()

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
                        carbam = i.carbam, carbam_drug = i.carbam, 
                        phenytoin = x.phenytoin, 
                        phenytoin_drug = x.phenytoin_drug)] %>% 
  as_tibble() %>%
  mutate(across(phenytoin:phenytoin_drug, ~replace_na(.x, FALSE))) %>% 
  arrange(id, trial_nr, desc(phenytoin)) %>% # This places any TRUE before FALSE
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup()

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
               carbam = i.carbam, carbam_drug = i.carbam, 
               phenytoin = i.phenytoin, 
               phenytoin_drug = i.phenytoin_drug, 
               primidone = x.primidone, 
               primidone_drug = x.primidone_drug
               )] %>% 
  as_tibble() %>%
  mutate(across(primidone:primidone_drug, ~replace_na(.x, FALSE))) %>% 
  arrange(id, trial_nr, desc(primidone)) %>% # This places any TRUE before FALSE
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup()




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
  filter(id == 31462) %>% 
  mutate(good_adh = missed_locf %in% 
           c("once a month", "never", "once every 2 weeks")) %>% 
  group_by(id) %>% 
  mutate(time_good_adh = case_when(
    good_adh == FALSE ~ ad_date - ad_date, 
    good_adh == TRUE & lag(good_adh) == FALSE ~ ad_date -  ad_date,
    good_adh == TRUE & row_number() == 1 ~ ad_date - ad_date,
    good_adh = TRUE & lag(good_adh) == TRUE ~ ad_date - lag(ad_date)),
    time_good_adh = as.numeric(time_good_adh)) %>% 
  group_by(id, grp = cumsum(time_good_adh == 0)) %>% 
  mutate(time_good_adh = cumsum(time_good_adh)) %>% 
  ungroup()

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
                      carbam = i.carbam, carbam_drug = i.carbam, 
                      phenytoin = i.phenytoin, 
                      phenytoin_drug = i.phenytoin_drug, 
                      primidone = i.primidone, 
                      primidone_drug = i.primidone_drug, 
                      time_good_adh = x.time_good_adh, 
                      ad_date = x.ad_date, 
                      adherence = x.missed,
                      adherence_locf = x.missed_locf
                    )] %>% 
  as_tibble() %>% 
  arrange(id, trial_nr, abs(trial_start - ad_date)) %>% 
  group_by(id, trial_nr) %>% 
  slice(1) %>%
  ungroup() %>% 
  mutate(time_good_adh_actual = if_else(time_good_adh == 0, 0, 
                                        time_good_adh + 
                                          as.numeric(trial_start - ad_date))) 




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
  select(-(born:male))

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
          carbam = i.carbam, carbam_drug = i.carbam, 
          phenytoin = i.phenytoin, 
          phenytoin_drug = i.phenytoin_drug, 
          primidone = i.primidone, 
          primidone_drug = i.primidone_drug, 
          time_good_adh = i.time_good_adh, 
          ad_date = i.ad_date, 
          adherence = i.adherence,
          adherence_locf = i.adherence_locf, 
          time_good_adh_actual = i.time_good_adh_actual, 
          crea = x.crea, 
          egfr = x.egfr, 
          cre_date = x.cre_date)] %>% 
  as_tibble() %>% 
  arrange(id, trial_nr, abs(trial_start - cre_date)) %>% 
  group_by(id, trial_nr) %>% 
  slice(1) %>% 
  ungroup()

