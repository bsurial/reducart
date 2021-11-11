library(tidyverse)  # ggplot, dplyr, %>%, and friends
library(ggdag)  # Make DAGs with ggplot
library(dagitty)  # Do basic DAG math
library(broom)  # For converting model output to data frames


node_details <- tribble(
  ~name, ~label, ~x, ~y,
  "switch", "Switch", 2, 1,
  "viral_failure", "Virologic Failure", 4, 1
)

set.seed(1)
dag_data <- dagify(VF ~ switch + DDI,
                   VF ~ resistance + adherence + history_VF + nadir_CD4 + max_RNA,
                   history_VF ~ adherence + nadir_CD4 + max_RNA,
                   switch ~ resistance + adherence + comorbidities + history_VF + DDI,
                   resistance ~ adherence + time_infected + history_VF,
                   switch ~ time_infected,
                   DDI ~ comorbidities,
                   labels = c("VF" = "Viral Failure (VF)", 
                              "switch" = "Switch",
                              "DDI" = "DDI",
                              "resistance" = "NRTI/NNRTI Resistance",
                              "adherence" = "Adherence",
                              "history_VF" = "Past VF",
                              "nadir_CD4" = "Nadir CD4",
                              "max_RNA" = "Max. HIV-VL",
                              "comorbidities" = "Co-morbidities", 
                              "time_infected" = "Time on ART"),
       exposure = "switch",
       outcome = "VF")

tidy_dag_data <- tidy_dagitty(dag_data, layout = "nicely")


ggdag_adjustment_set(tidy_dag_data, text_col = "black", shadow = TRUE, node_size = 20) + 
  theme_dag() + 
  theme(legend.position = "None")

ggsave(here::here("graphs", "03-adjusted_dag.png"), dpi = 300, bg = "white", 
       width = 7, height = 6)
