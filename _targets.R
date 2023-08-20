##### Set Configs #####
library(targets)
library(tarchetypes) 
library(tidyverse)
library(tidymodels)
library(pins)

tar_option_set(packages = c("tidyverse", "tidymodels", "qs", "future", "pins", "foreach"), format = "qs", deployment = "main", storage = "main", envir = globalenv(), tidy_eval = T, error = "continue", iteration = "list")

future::plan(future.callr::callr, workers = 16)
tar_source()
tar_source("../Helpers/SetGlobals.R")
tar_source("../Preprocessing/R/preprocessing.R")
tar_source("../Retrospective/R")

#### Run pipeline
list(
  
  # Read preprocessed data from LIS extract board.
  tar_target(preprocessed_bmp_inputs, list(BJH = read_csv(here::here("../../Data/BJH_inpatient_bmps_anonymized_sim_input.csv")))),
  tar_target(train_cohort, c("BJH")),
  tar_target(ref_ranges, read_csv(here::here("../../Data/ref_ranges_anonymized.csv"))),
  tar_target(bmp_no_NA, preprocessed_bmp_inputs %>% drop_na(any_of(lab_strings_bmp), ends_with("prior"), ends_with("post")), pattern = map(preprocessed_bmp_inputs)),
  
  # Remove contaminated and code specimens.
  tar_target(bmp_no_comment, bmp_no_NA %>% filter(!code_comment & !contam_comment & !unlikely_comment), pattern = map(bmp_no_NA)),
  #tar_target(bmp_uncontaminated, filterAnomalyWithResolution(bmp_no_comment, analyte_threshold = 3), pattern = map(bmp_no_comment)),
  tar_target(bmp_uncontaminated, bmp_no_comment),
  
  # Simulate across all mixture ratios.
  tar_target(fluids_tar, fluids),
  tar_target(fluid_names_tar, fluid_names),
  tar_target(fluid_analyte_map, getFluidMap()),
  tar_target(analytes_to_show, list(NS = c("calcium", "potassium_plas"), LR = c("calcium", "co2_totl"), D5W = c("sodium", "glucose"))),
  tar_target(contam_sim_all, makeFullContaminationTibble(bmp_uncontaminated, mix_ratios = rep(seq(0, 1, by = 0.01), each = 10000), fluid = fluids_tar, fluid_name = fluid_names_tar), pattern = cross(bmp_uncontaminated, map(fluids_tar, fluid_names_tar)), iteration = "list"),
  
  # Load wet-bench contamination experiments.
  tar_target(simulation_validation, data_board %>% pin_read("contamination_simulation_validation")),
  tar_target(val_long, convertContamValToLong(simulation_validation)),
  tar_target(val_nest, nestContamVal(val_long)),
  tar_target(val_reg, makeContamValRegModels(val_nest)),
  tar_target(val_figures, makeContamValRegFigures(val_long)),
 
  # Analyze the effects of simulated contamination.
  tar_target(sim_effect_figures, makeSimEffectFigures(contam_sim_all %>% bind_rows())),
  tar_target(sim_TAE_tables, makeTAEthresholdTables(contam_sim_all %>% bind_rows())),
  tar_target(sim_RCV_tables, makeRCVthresholdTables(contam_sim_all %>% bind_rows())),
  tar_target(sim_abnormal_flags, simulateRefRangeEffect(contam_sim_all %>% bind_rows(), ref_ranges)),
  tar_target(sim_flag_change_heatmaps, makeAbnormalFlagChangeHeatmap(contam_sim_all %>% bind_rows() %>% filter(mix_ratio == 0.10, label == fluid_names_tar) %>% left_join(ref_ranges)), pattern = map(fluid_names_tar)),
   
  # Calculate fluid distances.
  tar_target(normalizers, qs::qread("../Retrospective/_targets/objects/normalizers")),
  tar_target(prior_distances, getDistanceFromFluid(contam_sim_all %>% select(matches("_prior$"), -matches("delta")) %>% setNames(str_replace_all(names(.),  "_prior", "")), fluid_name = fluid_names_tar, label = paste0("dist_", fluid_names_tar, "_prior"), normalizers = normalizers), pattern = map(contam_sim_all, fluid_names_tar), iteration = "list"),
  tar_target(current_distances, getDistanceFromFluid(contam_sim_all %>% select(any_of(lab_strings)) %>% setNames(str_replace_all(names(.),  "_current", "")), fluid_name = fluid_names_tar, label = paste0("dist_", fluid_names_tar, "_current"), normalizers = normalizers), pattern = map(contam_sim_all, fluid_names_tar), iteration = "list"),
  tar_target(real_distances, getDistanceFromFluid(contam_sim_all %>% select(matches("_real")) %>% setNames(str_replace_all(names(.),  "_real", "")), fluid_name = fluid_names_tar, label = paste0("dist_", fluid_names_tar, "_real"), normalizers = normalizers), pattern = map(contam_sim_all, fluid_names_tar), iteration = "list"),
  tar_target(post_distances, getDistanceFromFluid(contam_sim_all %>% select(matches("_post$"), -matches("delta")) %>% setNames(str_replace_all(names(.),  "_post", "")), fluid_name = fluid_names_tar, label = paste0("dist_", fluid_names_tar, "_post"), normalizers = normalizers), pattern = map(contam_sim_all, fluid_names_tar), iteration = "list"),
  
  # Make figures for paper.
  tar_target(figure2, makeFigure2(contam_sim_all %>% bind_rows(), fluids = fluid_names, analytes = analytes_to_show)),
  tar_target(figure4, makeFigure4(contam_sim_all %>% bind_rows() %>% filter(mix_ratio == 0.10 & label %in% fluid_names) %>% left_join(ref_ranges), analytes = analytes_to_show))
  
)

