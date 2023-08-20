simulateContaminationRow <- function(input, mix_ratio, fluid){
  
  cols <- names(fluid)[which(names(fluid) %in% names(input))]
  
  output <- input %>%
    dplyr::mutate(across(all_of(cols), ~(1 - mix_ratio) * . + fluid[cur_column()] * mix_ratio)) %>%
    select(all_of(cols))
  
  output %>% 
    mutate(across(c("sodium", "chloride", "co2_totl", "bun", "glucose"), ~round(.))) %>%
    mutate(across(c("potassium_plas", "calcium"), ~round(., 1))) %>%
    mutate(creatinine = round(creatinine, 2)) %>% 
    mutate(anion_gap = sodium - chloride - co2_totl) %>%
    mutate(mix_ratio = mix_ratio)
  
}

getFluidMap <- function(){
  
  list(sodium = c("NS", "D5halfNS", "LR", "hyperNS"), 
       chloride = c("NS", "D5halfNS", "LR", "hyperNS"), 
       potassium_plas = c("NS", "LR", "D5halfNSwK"),
       co2_totl = c("NS"), 
       bun = c("NS"),
       creatinine = c("NS"),
       calcium = c("NS", "LR"),
       glucose = c("NS", "D5NS"))
  
}

makeFluidTable <- function(fluids){
  
  input = fluids %>% map(~.x[lab_strings_bmp]) %>% map(~.x[1:8])
  
  library(gt)
  gg_input = 
    input %>%
      bind_rows() %>% 
      select(sodium, chloride, potassium_plas, co2_totl, creatinine, bun, calcium, glucose) %>%
      setNames(analyte_labels[names(.)]) %>%
      mutate(fluid = names(fluids))
  
  gg_input %>%
      gt(rowname_col = "fluid") %>%
        tab_header(
          title = md("***Fluid Compositions***"),
          subtitle = "Six most commonly ordered fluids at BJH") %>%
        tab_footnote("Na, Cl, K in mEq/L, ") %>%
          gtsave("../../Figures/Simulation/fluid_table.pdf")
  
}

makeContaminationUnivariatePlots <- function(input = contam_sim_all, analytes = lab_strings_bmp[1:8], fluid_analyte_map = fluid_analyte_map, TAE = TAE_CLIA){
  
  library(foreach)
  library(geomtextpath)
  gg_list <- 
    foreach(analyte = analytes) %do%{
    
    tmp = input[which(fluid_names %in% fluid_analyte_map[[analyte]])] %>% bind_rows()
                                                      
    gg_wide = tmp %>% filter(mix_ratio > 0) %>% mutate(mix_ratio = factor(mix_ratio)) %>% select(label, mix_ratio, result = !!analyte)

    gg_input = 
      gg_wide %>%
        group_by(label, mix_ratio) %>% 
        summarise(p05 = quantile(result, probs = 0.05), 
                  p25 = quantile(result, probs = 0.25), 
                  median = quantile(result, probs = 0.5), 
                  p75 = quantile(result, probs = 0.75), 
                  p95 = quantile(result, probs = 0.95)) %>%
        mutate(mix_ratio = as.numeric(as.character(mix_ratio)))
    
    ggplot() +
      geom_smooth(data = gg_input, aes(x = mix_ratio, y = p05, color = label), alpha = 0.5, linewidth = 0.4, linetype = "dotted", se = F) +
      geom_smooth(data = gg_input, aes(x = mix_ratio, y = p25, color = label), alpha = 0.8, linewidth = 0.5, linetype = "dashed", se = F) +
      geom_smooth(data = gg_input, aes(x = mix_ratio, y = median, color = label), linewidth = 1, se = F) + 
      geom_smooth(data = gg_input, aes(x = mix_ratio, y = p75, color = label), alpha = 0.8, linewidth = 0.5, linetype = "dashed", se = F) +
      geom_smooth(data = gg_input, aes(x = mix_ratio, y = p95, color = label), alpha = 0.5, linewidth = 0.4, linetype = "dotted", se = F) + 
      scale_color_viridis_d(begin = 0.1, end = 0.9, option = "F") +
      ggtitle(analyte_labels[[analyte]]) +
      scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
      coord_cartesian(ylim = c(analyte_ranges[[analyte]][["min"]], analyte_ranges[[analyte]][["max"]])) +
      ylab("Result Value") + xlab("Mixture Ratio") + 
      theme(legend.position = "none", axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank(), title = element_text(size = 12, face = "bold.italic"))
    
    }
  

  library(ggpubr)
  library(grid)
  ggarrange(plotlist = gg_list, nrow = 2, ncol = 4, align = "hv") %>%  
    annotate_figure(left = textGrob("Result Value", rot = 90, vjust = 1, gp = gpar(cex = 1, fontface = "bold")),
                    bottom = textGrob("Mixture Ratio", gp = gpar(cex = 1, fontface = "bold")))
  ggsave("../../Figures/Simulation/Ranges/contam_sim_results_all_inpatient_mixtures.pdf", width = 8.5, height = 5)
  
}

makeTEaPlots <- function(input = contam_sim_all, analytes = lab_strings_bmp[1:8], fluid_analyte_map = fluid_analyte_map, TAE = TAE_CLIA){
  
  library(foreach)
  library(geomtextpath)
  gg_list <- 
    foreach(analyte = analytes) %do%{
      
      tmp = input[which(fluid_names %in% fluid_analyte_map[[analyte]])] %>% bind_rows()
      
      gg_wide = tmp %>% filter(mix_ratio > 0) %>% mutate(mix_ratio = factor(mix_ratio)) %>% select(label, mix_ratio, analyte, paste0(analyte, "_real"))
      gg_wide$error = abs(gg_wide[,3] - gg_wide[,4])
      
      gg_wide = gg_wide %>% mutate(exceeds_TEa = error > TAE[[analyte]][["threshold"]])
      
      gg_input = 
        gg_wide %>%
        group_by(label, mix_ratio, .drop = F) %>% 
        summarise(n = sum(exceeds_TEa), prop = n/10000) %>% 
        mutate(x = as.numeric(as.character(mix_ratio)))
      
      gg <- 
        gg_input %>% 
        ggplot(aes(x, prop, color = fct_rev(label), label = label)) +
        geom_textline(hjust = 0.2, straight = T, show.legend = F, size = 2, fontface = "bold", vjust = 1.25) +
        scale_color_viridis_d(begin = 0.1, end = 0.9, option = "F") +
        ggtitle(analyte_labels[[analyte]]) +
        scale_y_continuous(breaks = c(0, 0.5, 1), label = c(0, 0.5, 1)) + scale_x_continuous(breaks = c(0, 0.5, 1), label = c(0, 0.5, 1)) + 
        theme(axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank(), title = element_text(size = 12, face = "bold.italic"))
      ggsave(paste0("../../Figures/Simulation/Significant_Mixture_Thresholds/TEa_Proportion_By_Analyte/", analyte, "_TEa_cdf.pdf"), gg, width = 8, height = 6)  
      
      gg
      
    }
  
  library(ggpubr)
  library(grid)
  ggarrange(plotlist = gg_list, nrow = 2, ncol = 4) %>%  
    annotate_figure(left = textGrob("Proportion Exceeding TEa", rot = 90, vjust = 1, gp = gpar(cex = 1, fontface = "bold")),
                    bottom = textGrob("Mixture Ratio", gp = gpar(cex = 1, fontface = "bold")))
  ggsave("../../Figures/Simulation/Significant_Mixture_Thresholds/TEa_Proportion_By_Analyte/TEa_all.pdf", width = 8.5, height = 3)
  
}

makeFullContaminationTibble <- function(input, mix_ratios, fluid, fluid_name){
  
  input_prep = 
    input %>% drop_na(any_of(lab_strings)) %>%
      slice_sample(n = length(mix_ratios), replace = T) %>% 
      mutate(across(any_of(lab_strings), ~.x, .names = "{col}_real"), mix_ratio = mix_ratios, label = fluid_name)
  
  tmp = simulateContaminationRow(input_prep, input_prep$mix_ratio, fluid)
  
  input_prep[,names(tmp)] <- tmp
  
  input_prep %>% 
    mutate(label = fluid_name) %>%
    select(patient_id, specimen_id, drawn_dt_tm, label, mix_ratio, any_of(lab_strings), matches("_real"), matches("_prior"), matches("_post")) %>%
    addPrePostDeltas()
  
}

convertContamValToLong <- function(contam_validation){
  
  control_mean <- 
    contam_validation %>% 
      filter(Fluid == "Control") %>% 
      summarise(across(!matches("ID|Fluid|Mixture", ignore.case = F), ~mean(.x))) %>%
      mutate(across(c("sodium", "chloride", "co2_totl", "bun", "glucose"), ~round(.))) %>%
      mutate(across(c("potassium_plas", "calcium"), ~round(., 1))) %>%
      mutate(creatinine = round(creatinine, 2)) %>% 
      mutate(anion_gap = sodium - chloride - co2_totl) %>% select(any_of(lab_strings_bmp))
    
  
  controls <- tibble(control_mean, Fluid = c("NS", "LR", "D5NS", "D5LR")) %>% mutate(ID = 0, Mixture = 0)
  
  control_long <- tibble(Analyte = names(control_mean), Control_Result = unlist(control_mean))
  
  validation_long <- 
    contam_validation %>% 
      select(ID, Fluid, Mixture, any_of(lab_strings_bmp)) %>%
      bind_rows(controls) %>%
      filter(Fluid != "Control") %>%
      pivot_longer(!matches("Fluid|Mixture|ID", ignore.case = F), names_to = "Analyte", values_to = "Result") %>% 
      left_join(control_long) %>%
      mutate(Delta = Result - Control_Result, Prop = Delta/Control_Result)
  
  validation_long
  
}

nestContamVal <- function(validation_long){
  
  validation_long %>% 
    nest(data = everything(), .by = c(Fluid, Analyte))
  
}

makeContamValRegModels <- function(val_nest){
  
  model_list = 
    val_nest %>% 
      mutate(model = map(data, ~lm(Result ~ Mixture, data = .)) %>% setNames(paste(val_nest$Fluid, val_nest$Analyte, sep = "_")), 
             tidy = map(model, ~tidy(.x)) %>% setNames(paste(val_nest$Fluid, val_nest$Analyte, sep = "_")))
  
  model_list$intercept = model_list %>% pluck("tidy") %>% map(~filter(.x, term == "(Intercept)") %>% select("estimate") %>% pluck(1)) %>% unlist()
  model_list$slope = model_list %>% pluck("tidy") %>% map(~filter(.x, term == "Mixture") %>% select("estimate") %>% pluck(1)) %>% unlist()
  model_list = model_list %>% mutate(estimated_fluid_concentration = intercept + slope)
  
  bmp_models = model_list 
  
  bmp_models = 
    bmp_models %>% 
      mutate(sim_concentration = map2(bmp_models$Fluid, bmp_models$Analyte, ~fluids[[.x]][[.y]]) %>% unlist(), 
             sim_diff = sim_concentration - estimated_fluid_concentration, 
             sim_diff_prop = sim_diff/intercept)
  
  diff_summary = bmp_models %>% group_by(Analyte) %>% summarise(mean_diff = mean(sim_diff), mean_prop_diff = mean(sim_diff_prop))
  
  bmp_models
  
}

makeContamValRegFigures <- function(val_long){
  
  Analyte = map(val_long$Analyte, ~analyte_labels[[.x]]) %>% unlist()
  val_long$Analyte = map(val_long$Analyte, ~analyte_labels[[.x]]) %>% unlist()
  
  library(ggpmisc)
  ggplot(val_long %>% filter(Analyte != "Anion Gap"), aes(Mixture, Result)) + 
    geom_point() + 
    stat_poly_eq() + 
    stat_poly_line() + 
    facet_grid(Analyte ~ Fluid, scales = "free") + 
    theme(strip.text = element_text(size = 14, face = "bold.italic"))
  ggsave("../../Figures/Simulation/contam_sim_linear_concordance.pdf", width = 10, height = 6)
  
}

makeTAEthresholdTables <- function(input, label = ""){
  
  out_dir = paste0(here::here("../../Figures/Simulation/Significant_Mixture_Thresholds/"), label)
  
  if (!file.exists(out_dir)){
    dir.create(out_dir)
  }
  
  sim_error <- 
    input %>% 
      group_by(label, mix_ratio) %>%
      reframe(across(any_of(lab_strings) & !matches("anion_gap"), ~ abs(. - get(paste0(cur_column(), "_real")))))
  
  TAE_absolute_flags <- 
    sim_error %>% 
    mutate(across(any_of(lab_strings), ~ . > TAE_CLIA[[cur_column()]][["threshold"]])) %>% 
    pivot_longer(cols = any_of(lab_strings), names_to = "task_assay", values_to = "TAE_absolute_flag")
  
  TAE_absolute_counts <- 
    TAE_absolute_flags %>% 
    group_by(label, mix_ratio, task_assay, TAE_absolute_flag) %>% 
    count() %>% 
    group_by(label, mix_ratio, task_assay) %>%
    mutate(prop = n/sum(n))
  
  TAE_absolute_95th <- 
    TAE_absolute_counts %>% 
    ungroup() %>%
    filter(prop >= 0.95 & TAE_absolute_flag) %>% 
    group_by(label, task_assay) %>% 
    slice_head(n = 1) %>% 
    pivot_wider(id_cols = "label", names_from = "task_assay", values_from = "mix_ratio") %>% 
    left_join(tibble(label = fluid_names), .)
  
  TAE_absolute_median <- 
    TAE_absolute_counts %>% 
    ungroup() %>%
    filter(prop >= 0.5 & TAE_absolute_flag) %>% 
    group_by(label, task_assay) %>% 
    slice_head(n = 1) %>% 
    pivot_wider(id_cols = "label", names_from = "task_assay", values_from = "mix_ratio") %>% 
    left_join(tibble(label = fluid_names), .)
  
  TAE_absolute_5th <- 
    TAE_absolute_counts %>% 
    ungroup() %>%
    filter(prop >= 0.05 & TAE_absolute_flag) %>% 
    group_by(label, task_assay) %>% 
    slice_head(n = 1) %>% 
    pivot_wider(id_cols = "label", names_from = "task_assay", values_from = "mix_ratio") %>% 
    left_join(tibble(label = fluid_names), .)
  
  
  library(gt)
  
  TAE_absolute_95th %>%
    select(Fluid = label, any_of(lab_strings_bmp)) %>% 
    rename(any_of(relabeller)) %>%
    gt() %>%
    tab_header(
      title = md("***Minimum Significant Mixtures***"),
      subtitle = "95% of results contaminated at this ratio will exceed TAE thresholds."
    ) %>% 
    sub_missing(missing_text = "-") %>%
    tab_footnote(footnote = "TAE = Total Allowable Error (absolute)") %>% 
    tab_options(
      column_labels.padding = px(8),
      table.align = "center",
      table_body.vlines.color = "#EDEDED",
      table_body.vlines.width = px(1),
      table.border.bottom.color = "#EDEDED",
      heading.title.font.size = px(24),
      heading.subtitle.font.size = px(18)
    ) %>%
    cols_align("center") %>%
    tab_style(
      cells_column_labels(columns = everything()),
      style = list(cell_text(weight = "bold", size = "large"))
    ) %>%
    tab_style(
      cells_body(columns = Fluid),
      style = list(cell_text(style = "italic"))
    ) %>%
    gtsave(paste0(out_dir, "minimum_mixture_threshold_95_TAE_absolute_table.html"))
  
  TAE_absolute_median %>%
    select(` ` = label, any_of(lab_strings_bmp)) %>% 
    rename(any_of(relabeller)) %>%
    gt() %>%
    tab_header(
      title = md("***Minimum Significant Mixtures***"),
      subtitle = "50% of results contaminated at this ratio will exceed TEa thresholds."
    ) %>% 
    sub_missing(missing_text = "-") %>%
    tab_footnote(footnote = "TEa = Total Allowable Error. Absolute thresholds reported by CLIA or RPCA (chloride and CO2).") %>% 
    tab_options(
      column_labels.padding = px(8),
      table.align = "center",
      table_body.vlines.color = "#EDEDED",
      table_body.vlines.width = px(1),
      table.border.bottom.color = "#EDEDED",
      heading.title.font.size = px(24),
      heading.subtitle.font.size = px(18)
    ) %>%
    cols_align("center") %>%
    tab_style(
      cells_column_labels(columns = everything()),
      style = list(cell_text(weight = "bold", size = "large"))
    ) %>%
    tab_style(
      cells_body(columns = " "),
      style = list(cell_text(style = "italic", weight = "bold", align = "right"))
    ) %>%
    data_color(
      columns = !matches(" "), alpha = 0.5,
      fn = scales::col_bin(palette = c("#d19999", "#e8cccc"), domain = c(0, 1), bins = c(0, 0.1, 0.25), right = T, na.color = "white")
    ) %>%
    gtsave(paste0(out_dir, "minimum_mixture_threshold_50_TAE_CLIA_table.html"))

  return(NULL)
  
}

makeTAEat10Tables <- function(input, label = ""){
  
  out_dir = paste0(here::here("../../Figures/Simulation/Significant_Mixture_Thresholds/"), label)
  
  if (!file.exists(out_dir)){
    dir.create(out_dir)
  }
  
  sim_error <- 
    input %>% 
      group_by(label, mix_ratio) %>%
      reframe(across(any_of(lab_strings) & !matches("anion_gap"), ~ abs(. - get(paste0(cur_column(), "_real")))))
  
  TAE_flags <- 
    sim_error %>% 
      mutate(across(any_of(lab_strings), ~ . > TAE_CLIA[[cur_column()]][["threshold"]])) %>% 
      pivot_longer(cols = any_of(lab_strings), names_to = "task_assay", values_to = "tae_flag")
  
  gg_input <- 
    TAE_flags %>% 
      group_by(label, task_assay) %>% 
      count(tae_flag) %>% 
      mutate(prop = n/sum(n)) %>% 
      filter(tae_flag) %>%
      pivot_wider(id_cols = label, values_from = "prop", names_from = "task_assay", values_fill = 0) %>%
      select(Fluid = label, any_of(lab_strings_bmp)) %>% 
      rename(any_of(relabeller)) %>% 
      ungroup()
  
  library(gt)
  gg_input %>%
    gt() %>%
    tab_header(
      title = md("***Minimum Significant Mixtures***"),
      subtitle = "Proportion of Results Exceeded TEa at 10% Mixture"
    ) %>% 
    fmt_number() %>%
    sub_missing(missing_text = "-") %>%
    tab_footnote(footnote = "RCV = Reference Change Value (https://biologicalvariation.eu/)") %>% 
    tab_options(
      column_labels.padding = px(8),
      table.align = "center",
      table_body.vlines.color = "#EDEDED",
      table_body.vlines.width = px(1),
      table.border.bottom.color = "#EDEDED",
      heading.title.font.size = px(24),
      heading.subtitle.font.size = px(18)
    ) %>%
    cols_align("center") %>%
    tab_style(
      cells_column_labels(columns = everything()),
      style = list(cell_text(weight = "bold", size = "large"))
    ) %>%
    tab_style(
      cells_body(columns = Fluid),
      style = list(cell_text(style = "italic", weight = "bold", align = "right"))
    ) %>%
    gtsave(paste0(out_dir, "minimum_mixture_threshold_95_RCV_table.html"))
  

  }

makeRCVthresholdTables <- function(input, label = ""){
  
  out_dir = paste0(here::here("../../Figures/Simulation/Significant_Mixture_Thresholds/"), label)
  if (!file.exists(out_dir)){
    dir.create(out_dir)
  }
  setwd(out_dir)
  
  sim_error_prop <- 
    input %>% 
      group_by(label, mix_ratio) %>%
      reframe(across(any_of(lab_strings) & !matches("anion_gap"), ~ (. - get(paste0(cur_column(), "_real")))/get(paste0(cur_column(), "_real"))))
  
  RCV_percent_flags <- 
    sim_error_prop %>% 
      mutate(across(any_of(lab_strings), ~ (. > RCV_increase[[cur_column()]]) | (. < -RCV_decrease[[cur_column()]]))) %>% 
      pivot_longer(cols = any_of(lab_strings), names_to = "task_assay", values_to = "RCV_percent_flag")
  
  RCV_percent_counts <- 
    RCV_percent_flags %>% 
      group_by(label, mix_ratio, task_assay, RCV_percent_flag) %>% 
      count() %>% 
      group_by(label, mix_ratio, task_assay) %>%
      mutate(prop = n/sum(n))
  
  RCV_percent_95th <- 
    RCV_percent_counts %>% 
    ungroup() %>%
    filter(prop >= 0.95 & RCV_percent_flag) %>% 
    group_by(label, task_assay) %>% 
    slice_head(n = 1) %>% 
    pivot_wider(id_cols = "label", names_from = "task_assay", values_from = "mix_ratio") %>% 
    left_join(tibble(label = fluid_names), .)
  
  RCV_percent_median <- 
    RCV_percent_counts %>% 
      ungroup() %>%
      filter(prop >= 0.5 & RCV_percent_flag) %>% 
      group_by(label, task_assay) %>% 
      slice_head(n = 1) %>% 
      pivot_wider(id_cols = "label", names_from = "task_assay", values_from = "mix_ratio") %>% 
      left_join(tibble(label = fluid_names), .)
  
  RCV_percent_5th <- 
    RCV_percent_counts %>% 
    ungroup() %>%
    filter(prop >= 0.05 & RCV_percent_flag) %>% 
    group_by(label, task_assay) %>% 
    slice_head(n = 1) %>% 
    pivot_wider(id_cols = "label", names_from = "task_assay", values_from = "mix_ratio") %>% 
    left_join(tibble(label = fluid_names), .)
  
  
  library(gt)
  
  RCV_percent_95th %>%
    select(Fluid = label, any_of(lab_strings_bmp)) %>% 
    rename(any_of(relabeller)) %>%
    gt() %>%
    tab_header(
      title = md("***Minimum Significant Mixtures***"),
      subtitle = "95% of results contaminated at this ratio will exceed RCV thresholds."
    ) %>% 
    sub_missing(missing_text = "-") %>%
    tab_footnote(footnote = "RCV = Reference Change Value (https://biologicalvariation.eu/)") %>% 
    tab_options(
      column_labels.padding = px(8),
      table.align = "center",
      table_body.vlines.color = "#EDEDED",
      table_body.vlines.width = px(1),
      table.border.bottom.color = "#EDEDED",
      heading.title.font.size = px(24),
      heading.subtitle.font.size = px(18)
    ) %>%
    cols_align("center") %>%
    tab_style(
      cells_column_labels(columns = everything()),
      style = list(cell_text(weight = "bold", size = "large"))
    ) %>%
    tab_style(
      cells_body(columns = "Fluid"),
      style = list(cell_text(style = "italic", weight = "bold", align = "right"))
    ) %>%
    data_color(
      columns = !matches("Fluid"), alpha = 0.5,
      fn = scales::col_bin(palette = c("#d19999", "#e8cccc"), domain = c(0, 1), bins = c(0, 0.1, 0.25), right = T, na.color = "white")
    ) %>%
    gtsave(paste0(out_dir, "minimum_mixture_threshold_95_RCV_table.html"))
  
  RCV_percent_median %>%
    select(Fluid = label, any_of(lab_strings_bmp)) %>% 
    rename(any_of(relabeller)) %>%
    gt() %>%
    tab_header(
      title = md("***Minimum Significant Mixtures***"),
      subtitle = "50% of results contaminated at this ratio will exceed RCV thresholds."
    ) %>% 
    sub_missing(missing_text = "-") %>%
    tab_footnote(footnote = "RCV = Reference Change Value (https://biologicalvariation.eu/)") %>% 
    tab_options(
      column_labels.padding = px(8),
      table.align = "center",
      table_body.vlines.color = "#EDEDED",
      table_body.vlines.width = px(1),
      table.border.bottom.color = "#EDEDED",
      heading.title.font.size = px(24),
      heading.subtitle.font.size = px(18)
    ) %>%
    cols_align("center") %>%
    tab_style(
      cells_column_labels(columns = everything()),
      style = list(cell_text(weight = "bold", size = "large"))
    ) %>%
    tab_style(
      cells_body(columns = "Fluid"),
      style = list(cell_text(style = "italic", weight = "bold", align = "right"))
    ) %>%
    data_color(
      columns = !matches("Fluid"), alpha = 0.5,
      fn = scales::col_bin(palette = c("#d19999", "#e8cccc"), domain = c(0, 1), bins = c(0, 0.1, 0.25), right = T, na.color = "white")
    ) %>%
    gtsave(paste0(out_dir, "minimum_mixture_threshold_50_RCV_percentage_table.html"))
  
  
  RCV_percent_5th %>%
    select(Fluid = label, any_of(lab_strings_bmp)) %>% 
    rename(any_of(relabeller)) %>%
    gt() %>%
    tab_header(
      title = md("***Minimum Significant Mixtures***"),
      subtitle = "5% of results contaminated at this ratio will exceed RCV thresholds."
    ) %>% 
    sub_missing(missing_text = "-") %>%
    tab_footnote(footnote = "RCV = Reference Change Value (https://biologicalvariation.eu/)") %>% 
    tab_options(
      column_labels.padding = px(8),
      table.align = "center",
      table_body.vlines.color = "#EDEDED",
      table_body.vlines.width = px(1),
      table.border.bottom.color = "#EDEDED",
      heading.title.font.size = px(24),
      heading.subtitle.font.size = px(18)
    ) %>%
    cols_align("center") %>%
    tab_style(
      cells_column_labels(columns = everything()),
      style = list(cell_text(weight = "bold", size = "large"))
    ) %>%
    tab_style(
      cells_body(columns = "Fluid"),
      style = list(cell_text(style = "italic", weight = "bold", align = "right"))
    ) %>%
    data_color(
      columns = !matches("Fluid"), alpha = 0.5,
      fn = scales::col_bin(palette = c("#d19999", "#e8cccc"), domain = c(0, 1), bins = c(0, 0.1, 0.25), right = T, na.color = "white")
    ) %>%
    gtsave(paste0(out_dir, "minimum_mixture_threshold_5_RCV_percentage_table.html"))
  
  
  return(NULL)
  
}

simulateRefRangeEffect <- function(contam_sim_full, ref_ranges){
  
  input = left_join(contam_sim_full %>% bind_rows(), ref_ranges)
  
  library(foreach)
  gg_list = foreach(analyte = c("sodium", "chloride", "potassium_plas", "calcium", "glucose")) %do%{
    
    tmp = input %>% select(label, mix_ratio, matches(!!analyte)) %>% setNames(str_replace_all(names(.), eval(analyte), "analyte"))
    tmp <- tmp %>% mutate(real_flag = 
                            case_when(
                              analyte_real > analyte_normal_high & analyte_real < analyte_critical_high ~ "High", 
                              analyte_real >= analyte_critical_high ~ "Critical High",
                              analyte_real < analyte_normal_low & analyte_real > analyte_critical_low ~ "Low", 
                              analyte_real < analyte_critical_low ~ "Critical Low",
                              T ~ "Normal"))
    
    tmp <- tmp %>% mutate(sim_flag = 
                            case_when(
                              analyte > analyte_normal_high & analyte < analyte_critical_high ~ "High", 
                              analyte >= analyte_critical_high ~ "Critical High",
                              analyte < analyte_normal_low & analyte > analyte_critical_low ~ "Low", 
                              analyte <= analyte_critical_low ~ "Critical Low",
                              T ~ "Normal"))
    
    out <- tmp %>% mutate(analyte = !!analyte) %>% select(label, analyte, mix_ratio, real_flag, sim_flag)
    
    out
    
  } %>% bind_rows()
  
  gg_input <- 
    gg_list %>% 
    mutate(
      sim_flag = fct_rev(factor(sim_flag, levels = c("Critical Low", "Low", "Normal", "High", "Critical High"))),
      real_flag = fct_rev(factor(real_flag, levels = c("Critical Low", "Low", "Normal", "High", "Critical High"))),
      change = paste(real_flag, "to", sim_flag, sep = "-")) %>% 
    mutate(change = fct_rev(factor(change, levels = c("Critical Low-to-Critical Low", "Critical Low-to-Low", "Critical Low-to-Normal", "Critical Low-to-High", "Critical Low-to-Critical High", 
                                                      "Low-to-Critical Low", "Low-to-Low", "Low-to-Normal", "Low-to-High", "Low-to-Critical High",
                                                      "Normal-to-Critical Low", "Normal-to-Low", "Normal-to-Normal", "Normal-to-High", "Normal-to-Critical High",
                                                      "High-to-Critical Low", "High-to-Low", "High-to-Normal", "High-to-High", "High-to-Critical High",
                                                      "Critical High-to-Critical Low", "Critical High-to-Low", "Critical High-to-Normal", "Critical High-to-High", "Critical High-to-Critical High")))) %>%
    complete(label, analyte, change, fill = list(count = 0)) %>%
    group_by(label, analyte, mix_ratio) %>%
    count(sim_flag) %>% 
    mutate(prop = n/sum(n), 
           label = factor(label, levels = fluid_names), 
           analyte = factor(analyte, levels = c("sodium", "chloride", "potassium_plas", "calcium", "glucose"), labels = c("sodium" = "Sodium", "chloride" = "Chloride", "potassium_plas" = "Potassium", "calcium" = "Calcium",  "glucose" = "Glucose")))
  
  ggplot(gg_input %>% filter(label %in% c("NS", "LR", "D5NS", "D5LR", "D5W")), aes(x = mix_ratio, y = prop, fill = sim_flag, alpha = sim_flag)) + 
    geom_area(position = "stack") + 
    scale_fill_manual(values = c("#FA7F5EFF", "#FA7F5EFF", "grey60", "#20114BFF", "#20114BFF")) +
    scale_alpha_manual(values = c(0.9, 0.6, 0.3, 0.6, 0.9)) + 
    xlab("Mixture Ratio") + ylab("Proportion of Flags") + 
    scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) + 
    facet_grid(label ~ analyte) + 
    theme(legend.background = element_blank(), legend.position = "bottom", legend.title = element_blank())
  ggsave(here::here("../../Figures/Simulation/Significant_Mixture_Thresholds/abnormal_flag_proportions.pdf"), width = 8, height = 5)
  
}

getSimFlagDiffs <- function(input = contam_sim_all %>% bind_rows() %>% left_join(ref_ranges), analytes = analytes_to_show){
  
  library(foreach)
  gg_list = foreach(analyte = unique(unlist(analytes))) %do%{
    
    tmp = input %>% select(label, mix_ratio, matches(!!analyte)) %>% setNames(str_replace_all(names(.), eval(analyte), "analyte"))
    tmp <- tmp %>% mutate(real_flag = 
                            case_when(
                              analyte_real > analyte_normal_high & analyte_real < analyte_critical_high ~ "High", 
                              analyte_real >= analyte_critical_high ~ "Critical High",
                              analyte_real < analyte_normal_low & analyte_real > analyte_critical_low ~ "Low", 
                              analyte_real <= analyte_critical_low ~ "Critical Low",
                              T ~ "Normal"))
    
    tmp <- tmp %>% mutate(sim_flag = 
                            case_when(
                              analyte > analyte_normal_high & analyte < analyte_critical_high ~ "High", 
                              analyte >= analyte_critical_high ~ "Critical High",
                              analyte < analyte_normal_low & analyte > analyte_critical_low ~ "Low", 
                              analyte <= analyte_critical_low ~ "Critical Low",
                              T ~ "Normal"))
    
    out <- tmp %>% mutate(analyte = !!analyte) %>% select(label, analyte, mix_ratio, real_flag, sim_flag)
    
    out
    
  } %>% bind_rows()
  
  gg_list
  
}

makeAbnormalFlagChangeHeatmap <- function(input = contam_sim_all %>% bind_rows() %>% filter(mix_ratio == 0.10, label == fluid_names_tar) %>% left_join(ref_ranges), fluid = "D5NS", ratio = 0.10){
  
  gg_list = foreach(analyte = c("sodium", "chloride", "potassium_plas", "calcium", "glucose")) %do%{
    
    tmp = input %>% select(label, mix_ratio, matches(!!analyte)) %>% setNames(str_replace_all(names(.), eval(analyte), "analyte"))
    tmp <- tmp %>% mutate(real_flag = 
                            case_when(
                              analyte_real > analyte_normal_high & analyte_real < analyte_critical_high ~ "High", 
                              analyte_real >= analyte_critical_high ~ "Critical High",
                              analyte_real < analyte_normal_low & analyte_real > analyte_critical_low ~ "Low", 
                              analyte_real < analyte_critical_low ~ "Critical Low",
                              T ~ "Normal"))
    
    tmp <- tmp %>% mutate(sim_flag = 
                            case_when(
                              analyte > analyte_normal_high & analyte < analyte_critical_high ~ "High", 
                              analyte >= analyte_critical_high ~ "Critical High",
                              analyte < analyte_normal_low & analyte > analyte_critical_low ~ "Low", 
                              analyte <= analyte_critical_low ~ "Critical Low",
                              T ~ "Normal"))
    
    out <- tmp %>% mutate(analyte = !!analyte) %>% select(label, analyte, mix_ratio, real_flag, sim_flag)
    
    out
    
  } %>% bind_rows()
  
  gg_list %>% 
    mutate(index = 1, real_flag = factor(real_flag, levels = c("Critical Low", "Low", "Normal", "High", "Critical High")), sim_flag = factor(sim_flag, levels = c("Critical Low", "Low", "Normal", "High", "Critical High"))) %>% 
    group_by(analyte, real_flag, sim_flag) %>% 
    summarise(n = sum(index)) %>%
    group_by(analyte) %>%
    mutate(prop = round(n/sum(n) * 100, digits = 2)) %>%
    mutate(text_color = ifelse(n > 1000, "white", "white")) %>%
    ggplot(aes(real_flag, sim_flag, fill = n)) +
      geom_tile(show.legend = F) + 
      geom_text(aes(label = prop, color = text_color), fontface = "bold", show.legend = F) + 
      scale_color_identity() + 
      scale_fill_viridis_c(option = "B", begin = 0.1, end = 0.9, trans = scales::pseudo_log_trans(sigma = 0.001)) +
      facet_wrap(~ analyte, labeller = as_labeller(analyte_labels), nrow = 1) + 
      xlab("Ground Truth") + ylab("Contaminated Result Flag") + theme(plot.background = element_rect(fill = NULL))
  ggsave(here::here(paste0("../../Figures/Simulation/Significant_Mixture_Thresholds/Heatmaps/abnormal_flags_heatmaps_", fluid, "_", ratio, ".pdf")), device = "pdf", width = 18, height = 4)
  
}

makeSimEffectFigures <- function(input = contam_sim_all %>% bind_rows()){
  
  tmp_input <- 
    input %>% 
      select(label, mix_ratio, any_of(lab_strings)) %>% 
      pivot_longer(any_of(lab_strings), names_to = "Analyte", values_to = "Result") %>% 
      group_by(label, mix_ratio, Analyte) %>% 
      summarise(p05 = quantile(Result, probs = 0.05), 
                p25 = quantile(Result, probs = 0.25), 
                median = quantile(Result, probs = 0.5), 
                p75 = quantile(Result, probs = 0.75), 
                p95 = quantile(Result, probs = 0.95))
  
  library(foreach)
  gg_input = map(fluid_names, ~tmp_input %>% filter(label == .x) %>% mutate(Analyte = factor(Analyte, levels = lab_strings_bmp)))
  
  plots = 
    map2(gg_input, fluid_names, 
         ~ggplot() +
           geom_rect(data = .x, aes(xmin = mix_ratio, xmax = mix_ratio + 0.01, ymin = p05, ymax = p95, fill = label), alpha = 0.5) + 
           geom_rect(data = .x, aes(xmin = mix_ratio, xmax = mix_ratio + 0.01, ymin = p25, ymax = p75, fill = label), alpha = 0.8) + 
           geom_smooth(data = .x, aes(mix_ratio, median), se = F, color = "black", linetype = "dashed", alpha = 0.5) + 
           scale_fill_manual(values = color_map_global) + 
           scale_color_manual(values = color_map_global) + 
           scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
           facet_wrap(~ Analyte, scales = "free", labeller = as_labeller(analyte_labels), nrow = 1) + 
           ylab(.y) + xlab("Mixture Ratio") + 
           theme(legend.position = "none", axis.title.x = element_blank()))
  
  library(ggpubr)
  ggarrange(plots[[1]], plots[[2]], plots[[3]] + theme(axis.title.x = element_text(size = 18, face = "bold.italic", margin = margin(4,4,4,4))), nrow = length(plots), heights = c(rep(0.3, times = 5), 0.4), align = "hv")
  ggsave(paste0("../../Figures/Simulation/Ranges/contam_sim_results_all_mixture.pdf"), width = 18, height = 12)
         
}

makeFigure2 <- function(input = contam_sim_all %>% bind_rows(), fluids = c("NS", "LR", "D5W"), analytes = analytes_to_show){
  
  input = input %>% filter(label %in% fluids)
  
  tmp_input <- 
    input %>% 
      select(label, mix_ratio, any_of(lab_strings)) %>% 
      pivot_longer(any_of(lab_strings), names_to = "Analyte", values_to = "Result") %>% 
      group_by(label, mix_ratio, Analyte) %>% 
      summarise(p05 = quantile(Result, probs = 0.05), 
                p25 = quantile(Result, probs = 0.25), 
                median = quantile(Result, probs = 0.5), 
                p75 = quantile(Result, probs = 0.75), 
                p95 = quantile(Result, probs = 0.95))
  
  library(foreach)
  gg_input = map2(analytes, names(analytes), ~tmp_input %>% filter(label == .y & Analyte %in% .x) %>% mutate(Analyte = factor(Analyte, levels = .x)))
  
  plots = 
    map2(gg_input, names(gg_input), 
      ~ggplot() +
        geom_rect(data = .x, aes(xmin = mix_ratio, xmax = mix_ratio + 0.01, ymin = p05, ymax = p95, fill = label), alpha = 0.5) + 
        geom_rect(data = .x, aes(xmin = mix_ratio, xmax = mix_ratio + 0.01, ymin = p25, ymax = p75, fill = label), alpha = 0.8) + 
        geom_smooth(data = .x, aes(mix_ratio, median), se = F, color = "black", linetype = "dashed", alpha = 0.5) + 
        scale_fill_manual(values = color_map_global) + 
        scale_color_manual(values = color_map_global) + 
        scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
        facet_wrap(~ Analyte, scales = "free", labeller = as_labeller(analyte_labels)) + 
        ylab(.y) + xlab("Mixture Ratio") + 
        theme(legend.position = "none", axis.title.x = element_blank()))
  
  library(ggpubr)
  ggarrange(plots[[1]], plots[[2]], plots[[3]] + theme(axis.title.x = element_text(size = 18, face = "bold.italic", margin = margin(4,4,4,4))), nrow = length(plots), heights = c(0.35, 0.35, 0.4), align = "hv")
  ggsave("../../Figures/Simulation/Ranges/figure2.pdf", width = 10, height = 8)
  
}

makeFigure4 <- function(input = contam_sim_all %>% bind_rows() %>% filter(mix_ratio == 0.10 & label %in% c("NS", "LR", "D5W")) %>% left_join(ref_ranges), analytes = analytes_to_show){
  
  gg_list = foreach(analyte = unique(unlist(analytes))) %do%{
    
    tmp = input %>% select(label, mix_ratio, matches(!!analyte)) %>% setNames(str_replace_all(names(.), eval(analyte), "analyte"))
    tmp <- tmp %>% mutate(real_flag = 
                            case_when(
                              analyte_real > analyte_normal_high & analyte_real < analyte_critical_high ~ "High", 
                              analyte_real >= analyte_critical_high ~ "Critical High",
                              analyte_real < analyte_normal_low & analyte_real > analyte_critical_low ~ "Low", 
                              analyte_real <= analyte_critical_low ~ "Critical Low",
                              T ~ "Normal"))
    
    tmp <- tmp %>% mutate(sim_flag = 
                            case_when(
                              analyte > analyte_normal_high & analyte < analyte_critical_high ~ "High", 
                              analyte >= analyte_critical_high ~ "Critical High",
                              analyte < analyte_normal_low & analyte > analyte_critical_low ~ "Low", 
                              analyte <= analyte_critical_low ~ "Critical Low",
                              T ~ "Normal"))
    
    out <- tmp %>% mutate(analyte = !!analyte) %>% select(label, analyte, mix_ratio, real_flag, sim_flag)
    
    out
    
  } %>% bind_rows()
  
  props =   
    gg_list %>% 
      mutate(index = 1, real_flag = factor(real_flag, levels = c("Critical Low", "Low", "Normal", "High", "Critical High")), sim_flag = factor(sim_flag, levels = c("Critical Low", "Low", "Normal", "High", "Critical High"))) %>% 
      group_by(label, analyte, real_flag, sim_flag) %>% 
      summarise(n = sum(index)) %>%
      group_by(label, analyte) %>%
      mutate(prop = round(n/sum(n) * 100, digits = 2)) %>%
      mutate(text_color = ifelse(n > 1000, "white", "white"))
  
  grid = 
    expand_grid(real_flag = levels(props[[3]]), sim_flag = levels(props[[3]])) %>%
    as_tibble() %>%
    mutate(real_flag = factor(real_flag, levels = c("Critical Low", "Low", "Normal", "High", "Critical High")), 
           sim_flag = factor(sim_flag, levels = c("Critical Low", "Low", "Normal", "High", "Critical High"))) 
  
  correct_grid = grid %>% filter(real_flag == sim_flag)
  
  gg_input = map2(analytes, names(analytes), ~props %>% filter(label == .y & analyte %in% .x))
  
  plots = 
    map2(gg_input, names(gg_input),
      ~ggplot(.x, aes(real_flag, sim_flag, fill = n)) +
        geom_tile(data = grid, aes(real_flag, sim_flag), color = "grey50", fill = NA) +
        geom_tile(data = correct_grid, aes(real_flag, sim_flag), color = "grey50", fill = "grey90") +
        geom_tile(show.legend = F) + 
        geom_text(aes(label = prop, color = text_color), fontface = "bold", show.legend = F) + 
        scale_color_identity() + 
        scale_fill_viridis_c(option = "B", begin = 0.1, end = 0.9, trans = scales::pseudo_log_trans(sigma = 0.001)) +
        facet_wrap(~ analyte, labeller = as_labeller(analyte_labels), nrow = 1) + 
        ylab("Contaminated Result Flag") + xlab("Ground Truth Result Flag") + 
        theme(axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank(), panel.grid.minor = element_line()))
  
  library(ggpubr)
  ggarrange(plots[[1]], plots[[2]] + theme(axis.title.y.left = element_text(size = 18, face = "bold.italic", hjust = -2, margin = margin(4,4,4,4), angle = 90)), plots[[3]] + theme(axis.title.x.bottom = element_text(size = 18, face = "bold.italic", margin = margin(4,4,4,4))), nrow = length(plots), align = "hv")
  ggsave(here::here("../../Figures/Simulation/Significant_Mixture_Thresholds/Heatmaps/Figure4.svg"), width = 10, height = 10)
  
}

plotSimulatedAnomalyWithResolutionCounts <- function(input = contam_sim_all){
  
  gg_input <- 
    bind_cols(input[[1]], getAnomalyResolutionCount(input[[1]], "NS")) %>% 
    bind_rows() %>% 
    group_by(mix_ratio) %>% 
    reframe(p05 = quantile(AwR_count, probs = 0.05),
            p25 = quantile(AwR_count, probs = 0.25),
            median = quantile(AwR_count, probs = 0.50), 
            p75 = quantile(AwR_count, probs = 0.75),
            p95 = quantile(AwR_count, probs = 0.95))
  
  library(geomtextpath)
  ggplot(gg_input) + 
    geom_rect(xmin = -1, xmax = 0.10, ymin = -1, ymax = 4, alpha = 0.01, fill = NA, linewidth = 1.5, color = "darkred", linetype = "dashed") +
    geom_textpath(aes(x = mix_ratio, y = p05), label = "5%", alpha = 0.5, hjust = 0.7, size = 3) +
    geom_textpath(aes(x = mix_ratio, y = median), label = "Median", alpha = 1, hjust = 0.7, linewidth = 2, size = 4) +
    geom_textpath(aes(x = mix_ratio, y = p95), label = "95%", alpha = 0.5, hjust = 0.7, size = 3) +
    scale_x_continuous(limits = c(0, 0.51), breaks = seq(0, 0.5, by = 0.1), labels = seq(0, 0.5, by = 0.1), expand = c(0, 0)) + 
    scale_y_continuous(breaks = seq(0, 8), labels = seq(0, 8), expand = c(0, 0.1)) +
    ylab("Analytes Showing AwR Pattern") + xlab("Mixture Ratio") + 
    theme(panel.grid = element_blank(), panel.border = element_blank(), axis.title.x.bottom = element_text(size = 14), axis.title.y.left = element_text(size = 14), axis.text = element_text(size = 12, face = "bold"))
  ggsave("../../Figures/Simulation/Significant_Mixture_Thresholds/number_of_analytes_showing_anomaly_resolution_NS_sim.pdf", width = 8.5, height = 4)  
  
}

plotFlagChangesByMixture <- function(input = getSimFlagDiffs(contam_sim_all %>% bind_rows() %>% filter(label %in% c("NS", "LR", "D5W")) %>% left_join(ref_ranges), analytes = lab_strings_bmp_no_gap)){
  
  counts = 
    input %>%
      mutate(incorrect_flag = (real_flag != sim_flag), incorrect_critical = (incorrect_flag & (grepl("Critical", real_flag) | grepl("Critical", sim_flag)))) %>%
      group_by(label, analyte, mix_ratio)
  
  flag_counts <- counts %>% group_by(label, analyte, mix_ratio) %>% count(incorrect_flag)
  critical_counts <- counts %>% group_by(label, analyte, mix_ratio) %>% count(incorrect_critical)
  
  majority_flag = flag_counts %>% filter(incorrect_flag, n > 5000) %>% ungroup() %>% slice_head(n = 1, by = c(label, analyte))
  majority_critical = critical_counts %>% filter(incorrect_critical, n > 5000) %>% ungroup() %>% slice_head(n = 1, by = c(label, analyte))
  
  library(ggh4x)
  library(geomtextpath)
    ggplot() + 
      stat_smooth(data = flag_counts %>% filter(incorrect_flag), aes(mix_ratio, n), geom = "area", span = 0.05,
                fill = scico::scico(1, palette = "lajolla", begin = 0.75), alpha = 0.5, fullrange = F) +
      stat_smooth(data = critical_counts %>% filter(incorrect_critical), aes(mix_ratio, n), geom = "area", span = 0.05,
                  fill = scico::scico(1, palette = "lajolla", begin = 0.75), alpha = 0.75, fullrange = F) + 
      scale_y_continuous(breaks = c(0, 5000, 10000), labels = c("0%", "50%", "100%")) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = c("0", rep("", 4), "0.5", rep("", 4), "1")) +
      geom_segment(data = majority_flag, aes(x = mix_ratio, xend = mix_ratio, y = 0, yend = 5000), linetype = "dashed", alpha = 0.5) +
      geom_segment(data = majority_critical, aes(x = mix_ratio, xend = mix_ratio, y = 0, yend = 5000), linetype = "dashed", alpha = 0.5) +
      geom_textsmooth(data = flag_counts %>% filter(incorrect_flag & label == "NS" & analyte == "calcium"), aes(mix_ratio, n, label = "Abnormal"), span = 0.25, size = 3,
                      hjust = 0.2, vjust = -0.1, color = scico::scico(1, palette = "lajolla", begin = 0.75), alpha = 0.5, linewidth = NA, fontface = "bold") + 
      geom_textsmooth(data = critical_counts %>% filter(incorrect_critical & label == "NS" & analyte == "calcium"), aes(mix_ratio, n, label = "Critical"), span = 0.25,
                      hjust = 0.98, vjust = 1, color = scico::scico(1, palette = "lajolla", begin = 0.8), alpha = 1, linewidth = NA, fontface = "bold", size = 3) + 
      facet_grid2(fct_rev(label) ~ analyte, scales = "free", axes = "x", remove_labels = "y", labeller = as_labeller(c(analyte_labels, c("D5W" = "D5W", "NS" = "NS", "LR" = "LR"))), ) +
      xlab("Mixture Ratio") + ylab("Percentage Incorrectly Flagged") + 
      theme(axis.ticks.x = element_line())
  ggsave("../../Figures/Simulation/Ranges/percent_incorrectly_flagged_by_mixture.pdf", width = 8.5, height = 4)
    
}

plotYaleFlags <- function(input = preprocessed_bmp_inputs[[1]]){
  
  flags <- bind_cols(input, getYaleDeltaCheckRules(input))
  
  gg_input <- 
    flags %>% 
      transmute(yale_NS, yale_LR, yale_D5W, yale_any = (yale_NS | yale_LR | yale_D5W)) %>%
      drop_na()

  ggplot(gg_input, aes(yale_any)) +
    geom_bar()
  
}
