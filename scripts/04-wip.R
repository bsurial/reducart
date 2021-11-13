library(data.table)
library(tidyverse)
library(bernr)
library(lubridate)


set.seed(12)

df <- pro_read("02-proposal_population.rds")


tail <- pro_read("tail.rds")
art <- pro_read("modif_art.rds")
lab <- pro_read("lab.rds")
brand_dose <- pro_read("brand_dose.rds")

rna <- lab %>% 
  select(id, rna, labdate) %>% 
  drop_na()

df_art <- df %>% 
  select(-moddate) %>% 
  left_join(art %>% select(-(num_art:precision)))


# I first work on 1 patient, and extend that later on
df_sub <- df_art %>% 
  filter(id == 31462) %>% 
  select(id, treatment, moddate, enddate, switch) 


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



# Add drugs
# Rifampicin
rifa_df <- brand_dose %>% 
  filter(str_detect(brand, "(J04AM0[2567])|(J04AB02)")) %>% 
  mutate(rifamp = TRUE) %>% 
  select(id, startdate, stopdate, rifamp, rifa_drug = var_desc)

# Carboxamide derivatives (Carbamazepine, Oxcarbazepine etc)
carbam_df <- brand_dose %>% 
  filter(str_detect(brand, "N03AF0[1234]")) %>% 
  mutate(carbam = TRUE) %>% 
  select(id, startdate, stopdate, carbam, carbam_drug = var_desc)

# Phenytoin and other hydantoin derivatives
phenytoin_df <- brand_dose %>% 
  filter(str_detect(brand, "N03AB0[12345]|N03AB5[24]")) %>% 
  mutate(phenytoin = TRUE) %>% 
  select(id, startdate, stopdate, phenytoin, phyentoin_drug = var_desc)

# Primidone 
primidone_df <- brand_dose %>% 
  filter(str_detect(brand, "N03AA030")) %>% 
  mutate(primidone = TRUE) %>% 
  select(id, startdate, stopdate, primidone, primidone_drug = var_desc)
