library(tidyverse)
library(here)

resi_data <- read_csv(here("data", "2201_resistance", "2201_drm.csv")) %>% 
  select(-1) 

# Explore the 40 most common NRTI Resistance Mutations
resi_data %>% 
  separate_rows(NRTI, sep = ",") %>% 
  filter(!is.na(NRTI), NRTI != "None") %>% 
  mutate(NRTI = fct_lump(NRTI, n = 40),
         NRTI = fct_infreq(NRTI)) %>%
  count(NRTI, sort = T) %>% 
  ggplot(aes(y = NRTI, x = n)) + 
  geom_point() + 
  scale_x_log10() + 
  theme_minimal(base_family = "IBM Plex Sans")

# Explore the 40 most common INSTI Mutations (Major)
resi_data %>% 
  select(INSTI.Major_c) %>% 
  separate_rows(INSTI.Major_c, sep = ",") %>% 
  filter(!is.na(INSTI.Major_c), INSTI.Major_c != "None") %>% 
  filter(str_detect(INSTI.Major_c, "(Q148.?[HKR].?)|R263K|G118")) %>% 
  mutate(INSTI.Major_c = fct_lump(INSTI.Major_c, n = 40),
         INSTI.Major_c = fct_infreq(INSTI.Major_c)) %>%
  count(INSTI.Major_c, sort = T) %>% 
  ggplot(aes(y = INSTI.Major_c, x = n)) + 
  geom_point() + 
  scale_x_log10() + 
  theme_minimal(base_family = "IBM Plex Sans")

resi_data %>% 
  select(PATIENT, dat, NRTI_c, INSTI.Major_c, INSTI.Accessory_c) %>% 
  mutate(INSTI = str_detect(INSTI.Major_c, "(Q148.?[HKR].?)|R263K|G118")) %>% 
  filter(INSTI)

         