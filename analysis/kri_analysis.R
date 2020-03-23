# read in metadata from homeira's demographic information file

library(readxl)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(lme4)
library(stringr)

library(lubridate)

# custom function to import data from Homeira's data frames
parseCellDateTimes = function(cell_mappings, cell_dates, layer_name = 'L5', recorder_name = 'Homeira'){
  
  elem_df = lapply(1:length(cell_dates), function(ind){
    cell_inds = cell_mappings[!is.na(cell_mappings[, ind]), ind] %>% unlist
    date = cell_dates[[ind]] %>% as.character.Date()
    elem_df = data.frame(cell_inds = as.character(cell_inds), expt_date =  date)
    # print(elem_df)
    # elem_df$cell_ids = cell_inds
    # elem_df$expt_date = ephys_meta_dates[ind] %>% unlist
  }) %>% bind_rows()
  
  matching_inds = !str_detect(elem_df$cell_inds, '\\.abf')
  elem_df[matching_inds, 'cell_inds'] = paste0(elem_df[matching_inds,'cell_inds'], '.abf')
  
  elem_df$layer_name = layer_name
  elem_df$recorder_name = recorder_name
  return(elem_df)
  
}

### read in cell metadata information - expt dates, layers, recorder name, etc.

# read in mapping information from excel sheets in demographic file that associates experiment dates (subjects, usually) to cell .abf files
# note that information split up across two experimenters, Homeira and Lihua
demo_file_name = 'raw-data/final_kri_data/Demographic information Final version_March 12, 2020.xlsx'
l5_sheet_name = 'Layer 5- cells'
l23_sheet_name = 'L23-cells'
l3_sheet_name = 'Layer 3'

# layer 5
ephys_meta = read_excel(path = demo_file_name, sheet = l5_sheet_name, skip = 3, col_names = F, n_max = 12)
ephys_meta_dates = read_excel(path =demo_file_name , sheet = l5_sheet_name, skip = 2, col_names = F, n_max = 1)
l5_hom = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L5', recorder_name = 'Homeira')

ephys_meta = read_excel(path = demo_file_name, sheet = l5_sheet_name, skip = 18, col_names = F, n_max = 8)
ephys_meta_dates = read_excel(path = demo_file_name, sheet = l5_sheet_name, skip = 17, col_names = F, n_max = 1)
l5_lihua = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L5', recorder_name = 'Lihua')

# layer 2-3
ephys_meta = read_excel(path = demo_file_name, sheet = l23_sheet_name, range = 'B6:H11', col_names = F)
ephys_meta_dates = read_excel(path = demo_file_name, sheet = l23_sheet_name, skip = 3, col_names = F, n_max = 1)
l23_hom = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L2.3', recorder_name = 'Homeira')

ephys_meta = read_excel(path = demo_file_name, sheet = l23_sheet_name, range = 'B20:K26', col_names = F)
ephys_meta_dates = read_excel(path = demo_file_name, sheet = l23_sheet_name, skip = 17, col_names = F, n_max = 1)
l23_lihua = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L2.3', recorder_name = 'Lihua')

# layer 3c - just recordings made by homeira
ephys_meta = read_excel(path = demo_file_name, sheet = l3_sheet_name, range = 'A3:E10', col_names = F)
ephys_meta_dates = read_excel(path = demo_file_name, sheet = l3_sheet_name, skip = 1, col_names = F, n_max = 1)
l3c_hom = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L3c', recorder_name = 'Homeira')


# a final dataframe with all cell metadata
cell_meta_merged = bind_rows(l3c_hom, l23_hom, l5_hom, l5_lihua, l23_lihua) %>% as_tibble()
cell_meta_merged$expt_date = ymd(cell_meta_merged$expt_date)

cell_meta_merged = cell_meta_merged %>% distinct(cell_inds, .keep_all = T)

cell_meta_merged = cell_meta_merged %>% mutate(layer_name = factor(layer_name),
                             recorder_name = factor(recorder_name)) 


### read in patient chart info

# read in demographic information from more unstructured sheet
patient_meta = read_excel(path = demo_file_name,
                          range = 'A1:H61',
                          col_names = T,
                          col_types = c('date', 'guess', 'guess', 'guess', 'guess', 'guess', 'guess', 'guess'),
                          sheet = 'Demographic information-paper') %>%
  rename(expt_date = `Date of surgery`, age = `Age (Years)`, 
         sex = Sex, seizure_duration = "Years of seizures history", drugs = "Antiepileptic drugs",
         Diagnosis_long = Diagnosis, Resection_location_long = `Resection location`)

patient_meta$expt_date = ymd(patient_meta$expt_date)
patient_meta$subject_id = factor(make.names(patient_meta$expt_date), levels = unique(make.names(patient_meta$expt_date)))

# add a field to tell if subjects are unique by their date - some expt dates had two cases and not clear which cells from which patient
patient_meta$unique_subject = T
ambiguous_subject_ids = patient_meta %>% group_by(subject_id) %>% add_count() %>% filter(n > 1) %>% select(subject_id) %>% unlist()
patient_meta[patient_meta$subject_id %in% ambiguous_subject_ids, 'unique_subject'] = F

# load in demographic info from structured csv - contains normalized diagnosis and resection location
homeira_meta_structured = read.csv(file = 'raw-data/final_kri_data/Demographic information Final version_March 12, 2020_structured.csv') %>%
  select(-X) %>%
  rename(age = Age..Years., sex = Sex, seizure_duration = Years.of.seizure.history, drugs = Antiepileptic.drugs, diagnosis = Diagnosis,
         resection_location = Resection.location)

# merge patient metadata across two sources
patient_meta = left_join(patient_meta %>% select(-drugs), homeira_meta_structured, by = c('age', 'sex', 'seizure_duration')) %>% 
  arrange(expt_date) %>%
  select(subject_id, expt_date, age, sex, seizure_duration, unique_subject, diagnosis, resection_location, drugs, Diagnosis_long, Resection_location_long)

# create a field for associating cell to patient metadata if avail
patient_cell_meta = left_join(cell_meta_merged, patient_meta %>% distinct(expt_date, .keep_all = T), by = 'expt_date')
patient_cell_meta %<>% rename(cell_id = cell_inds)

# provide a simple df that has nicely structured metadata per patient and cell layer
updated_subject_metadata = patient_cell_meta %>% 
  group_by(subject_id, layer_name) %>% 
  add_count() %>% 
  ungroup() %>% 
  select(-cell_id) %>% distinct(.keep_all = T) %>% 
  arrange(expt_date) %>% 
  rename(cell_count = n) %>%
  select(subject_id, expt_date, cell_count, layer_name, recorder_name, 
         diagnosis, 
         age, sex, seizure_duration, drugs, unique_subject)

# better patient cell metadata
patient_cell_meta = left_join(cell_meta_merged, updated_subject_metadata, by = c('expt_date', 'recorder_name', 'layer_name'))
patient_cell_meta %<>% rename(cell_id = cell_inds)

# count total number of cells per patient
total_cell_count_matrix_unique = updated_subject_metadata %>% 
  filter(unique_subject, !is.na(sex), !is.na(age)) %>% 
  group_by(subject_id) %>% 
  summarize(total_cell_count = sum(cell_count)) %>% 
  ungroup()

write.csv(total_cell_count_matrix_unique, file = 'summary_tables/demographic_data_w_cell_counts.csv')



### pull ephys features (calc in matlab externally) from provided spreadsheets

ephys_sheets = c('Layer5_Homeira_With2020data_March09_20.xlsx', 'Layer5_Lihua_March06_20.xlsx', 
                 'Layer23_Homeira_March10_20.xlsx', 'Layer23_Lihua_March09_20.xlsx',
                 'Layer3_Human_Homeira_March 10_20.xlsx')

iv_ranges = c('A2:CJ10', 'A2:AY18', 'A2:AE10', 'A2:AK18', 'A2:O10')

kri_ephys_data = lapply(1:length(ephys_sheets), function(sheet_name){
  curr_path = paste0('raw-data/final_kri_data/', ephys_sheets[sheet_name])
  rmp_data = read_excel(path = curr_path, sheet = 7)
  colnames(rmp_data) = c('cell_inds', 'rmp')
  
  sag_data = read_excel(path = curr_path, sheet = 2, skip = 1)
  sag_names = paste0('sagamp', sag_data[, 1] %>% unlist) %>% make.names()
  
  sag_trans = sag_data[, -1] %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  colnames(sag_trans) = c('cell_inds', sag_names)
  
  sag_ratio_data = read_excel(path = curr_path, sheet = 14, skip = 1)
  sag_names = paste0('sag', sag_data[, 1] %>% unlist) %>% make.names()
  
  sag_ratio_trans = sag_ratio_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  colnames(sag_ratio_trans) = c('cell_inds', sag_names)
  
  ap_features = read_excel(path = curr_path, sheet = 9, skip = 1)
  ap_names = c('apamp', 'aphw', 'ap_rise', 'ap_fall')
  
  ap_features[ap_features < 0] = NA
  
  ap_trans = ap_features[, -1] %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  colnames(ap_trans) = c('cell_inds', ap_names)
  
  iv_data = read_excel(path = curr_path, sheet = 1, skip = 1, range = iv_ranges[sheet_name])
  
  rin = apply(iv_data[, -1], 2, function(vec){
   (vec/iv_data[1]) %>% unlist %>% mean(., na.rm = T) * 1000
  })
  
  if_hz_data = read_excel(path = curr_path, sheet = 3, skip = 1)
  if_hz_names = paste0('if_hz', if_hz_data[, 1] %>% unlist) %>% make.names()
  if_hz_trans = if_hz_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  colnames(if_hz_trans) = c('cell_inds', if_hz_names)
  if_hz_trans[if_hz_trans == -1] <- 0
  
  tau_data = read_excel(path = curr_path, sheet = 8, skip = 1)
  tau_names = paste0('tau', tau_data[, 1] %>% unlist) %>% make.names()
  tau_trans = tau_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  colnames(tau_trans) = c('cell_inds', tau_names)
  
  inst_freq_data = read_excel(path = curr_path, sheet = 12, skip = 1)
  inst_freq_names = paste0('inst_freq', inst_freq_data[, 1] %>% unlist) %>% make.names()
  inst_freq_trans = inst_freq_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  colnames(inst_freq_trans) = c('cell_inds', inst_freq_names)
  
  inst_freq_trans[inst_freq_trans == -1] <- 0
  
  # rebound spikes
  rebound_data = read_excel(path = curr_path, sheet = 5, skip = 1)
  rebound_names = paste0('rebound', rebound_data[, 1] %>% unlist) %>% make.names()
  rebound_trans = rebound_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  rebound_trans[-1] <- mutate_all(rebound_trans[-1], function(x) as.numeric(as.character(x)))
  colnames(rebound_trans) = c('cell_inds', rebound_names)
  #rebound_trans <- mutate_all(rebound_trans, function(x) as.numeric(as.character(x)))
  # print(rebound_trans)

  rin_df = rin %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  colnames(rin_df) = c('cell_inds', 'rin')
  # sag_names = paste0('sag', sag_data[, 1] %>% unlist) %>% make.names()
  
  rmp_data = rmp_data %>% as.data.frame()
  rmp_data = Reduce(merge, list(rmp_data, sag_trans, sag_ratio_trans, ap_trans, rin_df, if_hz_trans, 
                                tau_trans, inst_freq_trans, rebound_trans))
  return(rmp_data)
}) %>% bind_rows()

kri_ephys_data %<>% rename(cell_id = cell_inds)

kri_ephys_data$sag = kri_ephys_data$sag.400

joined_ephys_data = left_join(patient_cell_meta, kri_ephys_data %>% distinct(cell_id, .keep_all = T))

joined_ephys_data$cohort = 'KRI'

### define a final data frame that has just ephys data from cells with "complete" demographic data from patients
kri_valid_demo_cells = joined_ephys_data %>% filter(!is.na(age), 
                                                    !is.na(sex), 
                                                    unique_subject,
                                                    # diagnosis %in% c('TLE', 'Tumor')
                                                    ) %>% 
  distinct(subject_id, cell_id, .keep_all = T)


kri_valid_demo_cells %>% 
  ggplot(aes(x = age, y = sag.400)) + 
  geom_jitter(alpha = .7) + 
  geom_smooth(method = 'lm') + 
  facet_grid(~layer_name) + 
  xlab('Subject age') + ylab('Sag ratio (at -400 pA)')









