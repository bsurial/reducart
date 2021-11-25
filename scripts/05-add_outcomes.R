library(data.table)
library(tidyverse)
library(here)
library(bernr)
library(lubridate)


# Load data ---------------------------------------------------------------

elig_data <- pro_read("04-nested_trial_data.rds")
stop <- pro_read("stop.rds")
lab <- pro_read("lab.rds")
art <- pro_read("modif_art.rds")
pat <- pro_read("pat.rds")


# Add death and loss-to-fup -----------------------------------------------


# deaths
deaths <- stop %>% 
  filter(stop == 0) %>% 
  select(id, death_date = exitdate)

# loss-to-follow-up
loss_to_fup <- stop %>% 
  filter(stop %in% 1:6) %>% 
  select(id, loss_to_fup_date = stopdate)

# Join
df_death_fup <- elig_data %>% 
  select(id, trial_start, trial_nr) %>% 
  left_join(deaths, by = "id") %>% 
  left_join(loss_to_fup, by = "id")


# Add virologic data ------------------------------------------------------

rna <- lab %>% 
  select(id, rna, labdate) %>% 
  drop_na()




rna_filtered <- rna %>% 
  filter(id %in% df_death_fup$id) %>% 
  filter(labdate >= min(df_death_fup$trial_start)) %>% 
  arrange(id, labdate)


art_select <- art %>% 
  select(id, treatment, moddate, enddate) %>% 
  mutate(enddate = if_else(is.na(enddate), dmy("01.01.2100"), 
                           enddate))

setDT(rna_filtered); setDT(art_select)

rna_art <- art_select[rna_filtered,
           on = .(id, moddate <= labdate, 
                  enddate > labdate),
           .(id, 
            rna = i.rna, 
            labdate = i.labdate,
            treatment = x.treatment)] %>% 
  as_tibble()


# Time to confirmed virological failure, defined as 
# 2 consecutive HIV viral loads > 200 copies/mL (FAILURE 1) or 
# 1 viral load > 200 copies/mL followed by a treatment change (FAILURE 2)

failures <- rna_art %>% 
  group_by(id) %>% 
  mutate(failure = (rna > 200 & lag(rna) > 200), 
         failure2 = (rna > 200 & lag(treatment) != treatment))

failure_1 <- failures %>% 
  filter(failure == TRUE) %>% 
  select(id, failure_1_date = labdate) %>% 
  # Only take the first failure
  arrange(id, failure_1_date) %>% 
  group_by(id) %>% 
  summarise(failure_1_date = first(failure_1_date)) %>% 
  ungroup()

failure_2 <- failures %>% 
  filter(failure2 == TRUE) %>% 
  select(id, failure_2_date = labdate) %>% 
  # Only take the first failure
  arrange(id, failure_2_date) %>% 
  group_by(id) %>% 
  summarise(failure_2_date = first(failure_2_date)) %>% 
  ungroup()




df_death_fup_failure <- df_death_fup %>% 
  left_join(failure_1, by = "id") %>% 
  left_join(failure_2, by = "id")



# add Events and Censoring -----------------------------------------------------

event_data <- df_death_fup_failure %>% 
  # Add regdate 
  left_join(pat %>% select(id, regdate)) %>% 
  # failure 1 before failure 2
  mutate(event_type = case_when(
    !is.na(failure_1_date) & failure_1_date > trial_start ~ "failure 1",
    !is.na(failure_2_date) & failure_2_date > trial_start ~ "failure 2",
    TRUE ~ "no failure")) %>%  
  # loss_to_fup before death (some are loss to fup and found later to be dead)
  mutate(censor_reason = case_when(
    loss_to_fup_date > trial_start & !is.na(loss_to_fup_date) ~ "loss to fup", 
    death_date > trial_start & !is.na(death_date) ~ "death", 
    TRUE ~ "none")) %>% 
  # code event as 1 = failure 1, 2 = failure 2 and add censor dates
  mutate(event = if_else(event_type == "failure 1", 1, 
                         if_else(event_type == "failure 2", 2, 0)), 
         censor = if_else(censor_reason %in% c("death", "loss to fup"), 
                          1, 0),
         event_date = if_else(event_type == "failure 1", failure_1_date, NA_Date_), 
         censor_date =  case_when(
           censor_reason == "loss to fup" ~ loss_to_fup_date,
           censor_reason == "death" ~ death_date)) %>% 
  # Clean up
  select(id, trial_start, trial_nr, regdate, event, event_date, event_type, 
         censor, censor_date, censor_reason)
