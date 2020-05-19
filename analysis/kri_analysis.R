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
demo_file_name = 'raw-data/final_kri_data/updated_ephys_datasets_Apr_2020/Demographic information Final version_March 12, 2020 (4).xlsx'
l5_sheet_name = 'Layer 5- cells'
l23_sheet_name = 'L23-cells'
l3_sheet_name = 'Layer 3'

# layer 5
ephys_meta = read_excel(path = demo_file_name, sheet = l5_sheet_name, range = 'A4:AB12', col_names = F)
ephys_meta_dates = read_excel(path =demo_file_name , sheet = l5_sheet_name, range = 'A3:AB3', col_names = F)
l5_hom = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L5', recorder_name = 'Homeira')

ephys_meta = read_excel(path = demo_file_name, sheet = l5_sheet_name, range = 'A19:Q23', col_names = F)
ephys_meta_dates = read_excel(path = demo_file_name, sheet = l5_sheet_name, range = 'A18:Q18', col_names = F)
l5_lihua = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L5', recorder_name = 'Lihua')

# layer 2-3
ephys_meta = read_excel(path = demo_file_name, sheet = l23_sheet_name, range = 'B5:I10', col_names = F)
ephys_meta_dates = read_excel(path = demo_file_name, sheet = l23_sheet_name, range = 'B4:I4', col_names = F, n_max = 1)
l23_hom = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L2.3', recorder_name = 'Homeira')

ephys_meta = read_excel(path = demo_file_name, sheet = l23_sheet_name, range = 'B19:L26', col_names = F)
ephys_meta_dates = read_excel(path = demo_file_name, sheet = l23_sheet_name, range = 'B18:L18', col_names = F)
l23_lihua = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L2.3', recorder_name = 'Lihua')

# layer 3c - just recordings made by homeira
ephys_meta = read_excel(path = demo_file_name, sheet = l3_sheet_name, range = 'A3:F10', col_names = F)
ephys_meta_dates = read_excel(path = demo_file_name, sheet = l3_sheet_name, range = 'A2:F2', col_names = F)
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
                          range = 'A1:H62',
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

ephys_sheets = c('Layer5_Homeira_With2020data_April 01_20.xlsx', 'Layer5_Lihua_April 02_20.xlsx', 
                 'Layer23_Homeira_March 28_20.xlsx', 'Layer23_Lihua_March 29_20.xlsx',
                 'Layer3_Human_Homeira_April 06_20 - Editable.xlsx', 
                 'Layer5_Homeira_April 03_20_Interneurone.xlsx', 'Layer5_Lihua_March 30_20 - Interneurone.xlsx', 
                 'Layer23_Lihua_April 04_20_Interneurone.xlsx')

cell_layer_vec = c('L5', 'L5', 'L2.3', 'L2.3', 'L3c', 'L5', 'L5', 'L2.3')
cell_type_vec = c('Pyr', 'Pyr', 'Pyr', 'Pyr', 'Pyr', 'Int', 'Int', 'Int')

iv_ranges = c('A2:BZ10', 'A2:AW18', 'A2:AE10', 'A2:AF18', 'A2:T10',
              'A2:D10', 'A2:M18', 'A2:I18')

burst_row = c('A27:BZ27', 'A29:AW29', 
              'A24:AD24', 'A33:AF33',
              'A27:T27')

library(purrr)
kri_ephys_data = lapply(1:length(ephys_sheets), function(sheet_name){
  print(ephys_sheets[sheet_name])
  curr_path = paste0('raw-data/final_kri_data/updated_ephys_datasets_Apr_2020/', ephys_sheets[sheet_name])
  rmp_data = read_excel(path = curr_path, sheet = 'Sheet7')
  colnames(rmp_data) = c('cell_inds', 'rmp')
  
  # iv_data = read_excel(path = curr_path, sheet = 'Sheet1', skip = 1, range = iv_ranges[sheet_name])
  # 
  # rin = apply(iv_data[, -1], 2, function(vec){
  #  (vec/iv_data[1]) %>% unlist %>% mean(., na.rm = T) * 1000
  # })

  if_hz_data = read_excel(path = curr_path, sheet = 'Sheet3', skip = 1)
  if_hz_names = paste0('if_hz', if_hz_data[, 1] %>% unlist) %>% make.names()
  if_hz_trans = if_hz_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  colnames(if_hz_trans) = c('cell_inds', if_hz_names)
  if_hz_trans[if_hz_trans == -1] <- 0
  
  if (cell_type_vec[sheet_name] == 'Int'){
    burst_temp = list()
    ephys_data = rmp_data
    
  } else{
   burst_temp = read_excel(path = curr_path, sheet = 'Sheet3', range = burst_row[sheet_name], col_names = F)
  colnames(burst_temp) = colnames(if_hz_data)
  burst_temp = burst_temp[, 2:ncol(burst_temp)]
  burst_temp = burst_temp  %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  
  colnames(burst_temp) = c('cell_inds', 'Burst')
  ephys_data = Reduce(merge, list(rmp_data, burst_temp))
  
  }
  
  
  # sag_data = read_excel(path = curr_path, sheet = 2, skip = 1)
  # sag_names = paste0('sagamp', sag_data[, 1] %>% unlist) %>% make.names()
  # 
  # sag_trans = sag_data[, -1] %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  # colnames(sag_trans) = c('cell_inds', sag_names)
  # 
  # sag_ratio_data = read_excel(path = curr_path, sheet = 14, skip = 1)
  # sag_names = paste0('sag', sag_data[, 1] %>% unlist) %>% make.names()
  # 
  # sag_ratio_trans = sag_ratio_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  # colnames(sag_ratio_trans) = c('cell_inds', sag_names)
  # 
  # ap_features = read_excel(path = curr_path, sheet = 9, skip = 1)
  # ap_names = c('apamp', 'aphw', 'ap_rise', 'ap_fall')
  # 
  # ap_features[ap_features < 0] = NA
  # 
  # ap_trans = ap_features[, -1] %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  # colnames(ap_trans) = c('cell_inds', ap_names)
  # 
  # iv_data = read_excel(path = curr_path, sheet = 1, skip = 1, range = iv_ranges[sheet_name])
  # 
  # rin = apply(iv_data[, -1], 2, function(vec){
  #  (vec/iv_data[1]) %>% unlist %>% mean(., na.rm = T) * 1000
  # })
  # 
  # if_hz_data = read_excel(path = curr_path, sheet = 3, skip = 1)
  # if_hz_names = paste0('if_hz', if_hz_data[, 1] %>% unlist) %>% make.names()
  # if_hz_trans = if_hz_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  # colnames(if_hz_trans) = c('cell_inds', if_hz_names)
  # if_hz_trans[if_hz_trans == -1] <- 0
  # 
  # tau_data = read_excel(path = curr_path, sheet = 8, skip = 1)
  # tau_names = paste0('tau', tau_data[, 1] %>% unlist) %>% make.names()
  # tau_trans = tau_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  # colnames(tau_trans) = c('cell_inds', tau_names)
  # 
  # inst_freq_data = read_excel(path = curr_path, sheet = 12, skip = 1)
  # inst_freq_names = paste0('inst_freq', inst_freq_data[, 1] %>% unlist) %>% make.names()
  # inst_freq_trans = inst_freq_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  # colnames(inst_freq_trans) = c('cell_inds', inst_freq_names)
  # 
  # inst_freq_trans[inst_freq_trans == -1] <- 0
  # 
  # # rebound spikes
  # rebound_data = read_excel(path = curr_path, sheet = 5, skip = 1)
  # rebound_names = paste0('rebound', rebound_data[, 1] %>% unlist) %>% make.names()
  # rebound_trans = rebound_data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  # rebound_trans[-1] <- mutate_all(rebound_trans[-1], function(x) as.numeric(as.character(x)))
  # colnames(rebound_trans) = c('cell_inds', rebound_names)
  # #rebound_trans <- mutate_all(rebound_trans, function(x) as.numeric(as.character(x)))
  # # print(rebound_trans)
  # 
  # rin_df = rin %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  # colnames(rin_df) = c('cell_inds', 'rin')
  # # sag_names = paste0('sag', sag_data[, 1] %>% unlist) %>% make.names()
  
  #ephys_data = rmp_data %>% as.data.frame()
  # rmp_data = Reduce(full_join, list(rmp_data, sag_trans, sag_ratio_trans, ap_trans, rin_df, if_hz_trans, 
  #                               tau_trans, inst_freq_trans, rebound_trans))
  #rmp_data = Reduce(merge, list(rmp_data, sag_trans, sag_ratio_trans, ap_trans, rin_df))
  ephys_data$layer_name = cell_layer_vec[sheet_name]
  ephys_data$cell_type = cell_type_vec[sheet_name]
  print(ephys_data)
  return(ephys_data)
}) %>% bind_rows()

kri_ephys_data %<>% rename(cell_id = cell_inds)

kri_ephys_data$sag = kri_ephys_data$sag.400

joined_ephys_data = left_join(patient_cell_meta, kri_ephys_data %>% select(-layer_name) %>% distinct(cell_id, .keep_all = T)) %>% 
  select(cell_id, expt_date, layer_name, cell_type, recorder_name, subject_id, everything()) %>% arrange(expt_date)

cell_meta = joined_ephys_data %>% select(cell_id, expt_date, 
                                         layer_name, cell_type, recorder_name, subject_id, rmp) %>% 
  filter(!is.na(layer_name), !is.na(cell_type)) %>% 
  arrange(expt_date)

# add a final cell as an interneuron
cell_meta[cell_meta$cell_id == '15o08022.abf', 'cell_type'] = 'Int'

# manual mapping of L5 cell to L2.3
cell_meta[cell_meta$cell_id == '2016_03_01_0047.abf', 'layer_name'] = 'L2.3'

add_abf_file_ext = function(input_vector){
  return_vector = input_vector
  matching_inds = !str_detect(input_vector, '\\.abf')
  return_vector[matching_inds] = paste0(input_vector[matching_inds], '.abf')
  return(return_vector)
}

# get gain and offset info and sometimes RMPs from metadata sheets
# read in rmp metadata from metadata sheets
base_dir = 'raw-data/final_kri_data/gain_offset_info/'

gain_files = list.files(base_dir, full.names = T)

df_list = lapply(gain_files, function(f){
  return(read_excel(f))
})


df_list = df_list %>% bind_rows() 
df_list = df_list %>% rename(cell_id = `File name`, clampfit_gain = Gain, resp_offset = `Offset voltage`, 
                             resp_channel_name = `Response channel`, command_channel_name = `Command channel`)

df_list = df_list %>% filter(!is.na(cell_id))
df_list$cell_id = add_abf_file_ext(df_list$cell_id)

df_list = right_join(df_list, joined_ephys_data)

df_list[!is.na(df_list$RMP), 'rmp'] = df_list$RMP[!is.na(df_list$RMP)]


gain_offset_df = df_list %>% select(cell_id:command_channel_name, expt_date:subject_id, rmp) %>% tbl_df() %>% distinct(cell_id, .keep_all = T)
gain_offset_df[gain_offset_df$cell_id == '15622019.abf', 'resp_offset'] = -104.1 # need to hardcode this response offset as it is wrong in originating spreadsheet

gain_offset_df[is.na(gain_offset_df$resp_offset), 'resp_offset'] = 0

write.csv(gain_offset_df, file = 'summary_tables/gain_offset_info_merged.csv')

### merge cell meta temp thing with ephys data we got from ipfx in python

kri_extracted_features = read.csv('~/Documents/GitHub/valiante_lab_abf_process/output_files/extracted_cell_features.csv') %>% distinct(cell_id, .keep_all = T)

#kri_extracted_features = left_join(joined_ephys_data, kri_extracted_features, by = 'cell_id') %>% distinct(cell_id, .keep_all = T)

kri_extracted_features = left_join(kri_extracted_features, gain_offset_df %>% select(cell_id, rmp))


library(umap)
kri_extracted_features %<>% mutate(
  rin = input_resistance, 
  apthr = threshold_v,
  apamp = peak_v - threshold_v, 
  ahpamp = threshold_v - fast_trough_v,
  aphw = width * 1000, 
  apvel = upstroke_downstroke_ratio,
  rheo = rheobase_i,
  #maxfreq = max_rate_long_square, 
  adratio = adapt, 
  sag = sag, # defined by aibs as vm change / max vm deflection
  sagamp = (v_baseline -vm_for_sag) * sag,
  fislope = fi_fit_slope, 
  avgisi = mean_isi, 
  cvisi = isi_cv, 
  latency = latency * 1000, # convert to ms
  tau = tau * 1000
)

kri_ephys_features = c('rin', 'rmp', 'apamp', 'ahpamp', 'aphw', 'apvel', 'sagamp',
                       'adratio', 'first_isi', 'avgisi', 'cvisi',
                       'sag', 'fislope',  'latency', 'avg_rate', 'tau', 'rheo',
                       'apthr',
                       'peak_t', 'fast_trough_t', 'trough_t')

new_interneuron_cell_ids = kri_extracted_features %>% 
  mutate(ahp_ratio = ahpamp / apamp) %>% 
  filter(ahp_ratio > .1, fi_fit_slope > .2) %>% 
  pull(cell_id)

# get RMP from kri sheet, not joined_ephys_data
joined_ephys_data[joined_ephys_data$cell_id %in% new_interneuron_cell_ids, 'cell_type'] = 'Int'
joined_ephys_data[is.na(joined_ephys_data$cell_type), 'cell_type'] = 'Pyr'


cell_meta = joined_ephys_data %>% mutate(cell_type = factor(cell_type),
                                         has_burst = !is.na(Burst)) %>% select(cell_id, expt_date, 
                                         layer_name, cell_type, recorder_name, subject_id, has_burst) %>% 
  filter(!is.na(layer_name), !is.na(cell_type)) %>% 
  arrange(expt_date)

# add a final cell as an interneuron
cell_meta[cell_meta$cell_id == '15o08022.abf', 'cell_type'] = 'Int'

# manual mapping of L5 cell to L2.3
cell_meta[cell_meta$cell_id == '2016_03_01_0047.abf', 'layer_name'] = 'L2.3'

cell_meta = left_join(cell_meta, gain_offset_df %>% select(cell_id, rmp))

# import abf metadata sheet 
cell_abf_metadata = read.csv('~/Documents/GitHub/valiante_lab_abf_process/output_files/cell_final_raw_meta_df.csv')

cell_meta = left_join(cell_meta, 
          cell_abf_metadata %>% select(file_time, abf_tag, cell_id, rmp_error) %>% 
            rename(acquisition_time = file_time, 
                   tag_comments = abf_tag, 
                   voltage_drift = rmp_error)) %>% 
  arrange(acquisition_time)


# go through each cell in cell_meta and ensure that there aren't abf file duplicates due to file time issues

df_list = list()
df_list = apply(cell_meta, 1, function(row){
  curr_row_id = row['cell_id'] %>% as.character()
  curr_file_time = row['acquisition_time'] 
  check_inds = which(abs(difftime(curr_file_time, cell_meta$acquisition_time, units = "mins")) < 20)
  
  check_ids = setdiff(cell_meta[check_inds, 'cell_id'] %>% unlist(), curr_row_id)
  
  if (length(check_ids) > 0){
    curr_df = cell_meta %>% filter(cell_id %in% c(curr_row_id, check_ids))
    df_list = c(df_list, curr_df)
    print(curr_df)
  } else{
    curr_df = list()
  }
  return(curr_df)

})

cells_to_drop = c('13d03029.abf', 
                  '19122022.abf', 
                  '19129072.abf', 
                  '2019_11_04_0113.abf', 
                  '2015_11_09_0085.abf')

cell_meta = cell_meta %>% filter(!cell_id %in% cells_to_drop)

cells_w_morphology = c('19122003.abf', # smaller L2 morphlogy, first on left
                       '2017_03_20_0026.abf', # synaptic blockers?
                       '19320017.abf', # small L2 cell, not currently highlighted in paper
                       '19122043.abf', # cell currently labelled as L2 350 um
                       '19129004.abf', # large L3 morphology
                       '19320022.abf', # L3c cell, synaptic blockers
                       '19219004.abf', # L3c cell, will be added to paper, synaptic blockers
                       '19219021.abf', # L3c cell, will be added to paper, synaptic blockers
                       '19128061.abf', # big L5 morphology
                       '19129015.abf', # smaller L5 morphology
                       '2020_01_06_0048.abf' # L3c cell
) 

cell_meta = cell_meta %>% mutate(has_morphology = if_else(cell_id %in% cells_w_morphology, T, F))

cell_meta = cell_meta %>% 
  select(cell_id, layer_name, cell_type, recorder_name, subject_id, acquisition_time, rmp, voltage_drift, has_morphology, tag_comments, has_burst) %>% 
  arrange(acquisition_time)

write.csv(cell_meta, file = 'summary_tables/cell_metadata.csv')

cell_meta_w_ephys = left_join(cell_meta %>% select(-rmp), 
                              kri_extracted_features %>% select(cell_id, kri_ephys_features) %>% 
                                distinct(cell_id, .keep_all = T)) %>% 
  arrange(acquisition_time)

# read zap data from summary spreadsheet
zap_data_w_intrinsic_id_final = read.csv('summary_tables/zap_data_final_summarized.csv')

#cell_meta_w_ephys[cell_meta_w_ephys$has_resonance %>% is.na, 'has_resonance'] = 'not measured'

cell_patient_ephys_combined = left_join(cell_meta %>% select(-rmp), 
                                        patient_meta %>% filter(unique_subject == T) %>%
                                          rename(resection_date = expt_date)) %>% arrange(acquisition_time)

cell_patient_ephys_combined = left_join(cell_patient_ephys_combined %>% select(-has_burst, everything()), 
                                        kri_extracted_features %>% select(cell_id, kri_ephys_features) %>% 
                                          distinct(cell_id, .keep_all = T)) %>% 
  arrange(acquisition_time)

cell_patient_ephys_combined = left_join(cell_patient_ephys_combined, zap_data_w_intrinsic_id_final %>% select(cell_id, res_center_freq, res_3dB_freq, res_sharpness, has_resonance))


write.csv(cell_patient_ephys_combined, file = 'summary_tables/cell_patient_ephys_combined.csv')

# write donor metadata sheet

write.csv(patient_meta, file = 'summary_tables/donor_metadata.csv')



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
total_cell_count_matrix_unique = updated_subject_metadata[updated_subject_metadata$subject_id == 'subject1', ] %>% 
  filter(unique_subject, !is.na(sex), !is.na(age)) %>% 
  group_by(subject_id) %>% 
  summarize(total_cell_count = sum(cell_count)) %>% 
  ungroup()

write.csv(total_cell_count_matrix_unique, file = 'summary_tables/demographic_data_w_cell_counts.csv')


cell_patient_ephys_combined = left_join(cell_meta, patient_meta %>% filter(unique_subject == T)) %>% 
  select(-rmp, everything(), rmp)  %>% 
  left_join(., kri_extracted_features %>% distinct(cell_id, .keep_all = T))

write.csv(cell_patient_ephys_combined, file = 'summary_tables/cell_patient_ephys_combined.csv')


joined_ephys_data$cohort = 'KRI'

kri_ephys_data_all_cells = joined_ephys_data
write.csv(kri_ephys_data_all_cells, file = 'summary_tables/all_cells.csv')


### define a final data frame that has just ephys data from cells with "complete" demographic data from patients
kri_valid_demo_cells = joined_ephys_data %>% filter(!is.na(age), 
                                                    !is.na(sex), 
                                                    unique_subject,
                                                    # diagnosis %in% c('TLE', 'Tumor')
                                                    ) %>% 
  distinct(subject_id, cell_id, .keep_all = T)


write.csv(kri_valid_demo_cells, file = 'summary_tables/cells_w_demographic_data.csv')



kri_valid_demo_cells %>% 
  ggplot(aes(x = age, y = sag.400)) + 
  geom_jitter(alpha = .7) + 
  geom_smooth(method = 'lm') + 
  facet_grid(~layer_name) + 
  xlab('Subject age') + ylab('Sag ratio (at -400 pA)')









