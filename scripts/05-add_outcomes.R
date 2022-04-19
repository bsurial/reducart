
# Load packages -----------------------------------------------------------

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
fup <- pro_read("fup.rds")




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
# 1 viral load > 200 copies/mL followed by a treatment change (FAILURE 2) or
# any HIV viral load > 200 copies/mL

failures <- rna_art %>% 
  group_by(id) %>% 
  mutate(failure = (rna > 200 & lag(rna) > 200), 
         failure2 = (rna > 200 & lag(treatment) != treatment), 
         failure3 = rna > 200, 
         blip = rna > 50)

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

failure_3 <- failures %>% 
  filter(failure3 == TRUE) %>% 
  select(id, failure_3_date = labdate) %>% 
  # Only take the first failure
  arrange(id, failure_3_date) %>% 
  group_by(id) %>% 
  summarise(failure_3_date = first(failure_3_date)) %>% 
  ungroup()




df_death_fup_failure <- df_death_fup %>% 
  left_join(failure_1, by = "id") %>% 
  left_join(failure_2, by = "id") %>% 
  left_join(failure_3, by = "id")


# add Last FUP -----------------------------------------------------

last_fup <- fup %>% 
  select(id, fupdate) %>% 
  group_by(id) %>% 
  summarise(last_fupdate = last(fupdate, order_by = fupdate)) %>% 
  ungroup()

# add Events and Censoring -----------------------------------------------------

switch_date <- elig_data %>% 
  filter(switch_treatment == 1) %>% 
  select(id, trial_start, switch_treatment) %>% 
  arrange(id, trial_start) %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(id, switch_date = trial_start)

event_data <- df_death_fup_failure %>% 
  # add switch_date
  left_join(switch_date) %>% 
  # Add regdate 
  left_join(pat %>% select(id, regdate)) %>% 
  # add last fup 
  left_join(last_fup) %>% 
  mutate(event_type = case_when(
    !is.na(failure_1_date) & failure_1_date > trial_start ~ "failure 1",
    !is.na(failure_2_date) & failure_2_date > trial_start ~ "failure 2",
    !is.na(failure_3_date) & failure_3_date > trial_start ~ "failure 3",
    TRUE ~ "no failure")) %>% 
  # loss_to_fup before death (some are loss to fup and found later to be dead)
  filter(trial_start < death_date | is.na(death_date), 
         trial_start < loss_to_fup_date |  is.na(loss_to_fup_date)) %>% 
  mutate(time_date = pmin(failure_1_date, failure_2_date, 
                          switch_date, death_date, loss_to_fup_date, 
                          last_fupdate, na.rm = TRUE),
         time3_date = pmin(failure_1_date, failure_2_date, failure_3_date, 
                           switch_date, death_date, loss_to_fup_date, 
                           last_fupdate, na.rm = TRUE)) %>% 
  mutate(event_type = case_when(time_date == failure_1_date ~ "failure 1", 
                                time_date == failure_2_date ~ "failure 2", 
                                time3_date == failure_3_date ~ "failure 3", 
                                TRUE ~ "no failure"),
         censor_reason = case_when(time_date == death_date ~ "death", 
                                   time_date == loss_to_fup_date ~ "loss to fup",
                                   time_date == switch_date ~ "switched", 
                                   time_date == last_fupdate ~ "last fup date",
                                   TRUE ~ "none")) %>% 
  mutate(event = as.numeric(event_type %in% 
                              c("failure 1", "failure 2")), 
         event3 = as.numeric(event_type %in% 
                               c("failure 1", "failure 2", "failure 3"))) %>% 
  select(id:trial_nr, time_date, time3_date, event_type, censor_reason, event, event3)




# Add time variable -------------------------------------------------------


# Full data
full_df <- elig_data %>% 
  left_join(event_data, by = c("id", "trial_start", "trial_nr"))


current <- full_df %>% 
  # Filter only those who are eligible and where trial starts before time
  filter(elig_current == 1) 

switch <- full_df %>% 
  filter(elig_switch == 1) %>% 
  arrange(id, trial_start) %>% 
  group_by(id) %>% 
  slice(1) %>% # Some have 2 switch episodes. Remove them here
  ungroup()

aset <- bind_rows(current, switch) %>% 
  filter(trial_start < time_date)




# Add outcome variables
aset <- aset %>% 
  mutate(time = as.numeric((time_date - trial_start) / 365.25), 
         time_3 = as.numeric((time3_date - trial_start) / 365.25)) %>% 
  mutate(group = if_else(elig_switch == 1, "Switch", "Current"), 
         group_numeric = as.numeric(group == "Switch"))





write_rds(aset, here("processed", "05-analysis_data.rds"))



