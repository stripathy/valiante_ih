# read in metadata from homeira's demographic information file

library(readxl)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(cowplot)

library(lubridate)

# custom function to import data from Homeira's data frames
parseCellDateTimes = function(cell_mappings, cell_dates, layer_name = 'L5', recorder_name = 'Homeira'){
  
  elem_df = lapply(1:length(ephys_meta_dates), function(ind){
    cell_inds = ephys_meta[!is.na(ephys_meta[, ind]), ind] %>% unlist
    date = ephys_meta_dates[[ind]] %>% as.character.Date()
    elem_df = data.frame(cell_inds = cell_inds, expt_date =  date)
    # elem_df$cell_ids = cell_inds
    # elem_df$expt_date = ephys_meta_dates[ind] %>% unlist
  }) %>% bind_rows()
  
  elem_df$layer_name = layer_name
  elem_df$recorder_name = recorder_name
  return(elem_df)
  
}

# layer 5
ephys_meta = read_excel(path = '~/valiante_ih/raw-data/Demographic information_Request_HM.xlsx', sheet = 'Layer 5- cells', skip = 3, col_names = F, n_max = 12)
ephys_meta_dates = read_excel(path = '~/valiante_ih/raw-data/Demographic information_Request_HM.xlsx', sheet = 'Layer 5- cells', skip = 2, col_names = F, n_max = 1)
l5_hom = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L5', recorder_name = 'Homeira')

ephys_meta = read_excel(path = '~/valiante_ih/raw-data/Demographic information_Request_HM.xlsx', sheet = 'Layer 5- cells', skip = 18, col_names = F, n_max = 8)
ephys_meta_dates = read_excel(path = '~/valiante_ih/raw-data/Demographic information_Request_HM.xlsx', sheet = 'Layer 5- cells', skip = 17, col_names = F, n_max = 1)
l5_lihua = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L5', recorder_name = 'Lihua')

# layer 2-3

ephys_meta = read_excel(path = '~/valiante_ih/raw-data/Demographic information_Request_HM.xlsx', sheet = 'L23-cells', range = 'B6:G12', col_names = F, )
ephys_meta_dates = read_excel(path = '~/valiante_ih/raw-data/Demographic information_Request_HM.xlsx', sheet = 'L23-cells', skip = 3, col_names = F, n_max = 1)
l23_hom = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L2.3', recorder_name = 'Homeira')

ephys_meta = read_excel(path = '~/valiante_ih/raw-data/Demographic information_Request_HM.xlsx', sheet = 'L23-cells', range = 'B20:K31', col_names = F)
ephys_meta_dates = read_excel(path = '~/valiante_ih/raw-data/Demographic information_Request_HM.xlsx', sheet = 'L23-cells', skip = 17, col_names = F, n_max = 1)
l23_lihua = parseCellDateTimes(cell_mappings = ephys_meta, cell_dates = ephys_meta_dates, layer_name = 'L2.3', recorder_name = 'Lihua')

cell_meta_merged = bind_rows(l5_hom, l5_lihua, l23_hom, l23_lihua) %>% as_tibble()
cell_meta_merged$expt_date = ymd(cell_meta_merged$expt_date)

cell_meta_merged %<>% mutate(layer_name = factor(layer_name),
                             recorder_name = factor(recorder_name))

# 
# # read in patient chart info
# 
# patient_meta = read_excel(path = '~/valiante_ih/raw-data/2018_Human_Tissue_Anon_.xlsx', 
#                           range = 'A1:L109', 
#                           col_names = T, 
#                           col_types = c('guess', 'guess', 'guess', 'skip', 'date', 'guess', 'numeric', 'guess', 'guess', 'guess', 'guess', 'guess')) %>% 
#   rename(expt_date = `Resection Date`, age = `Age At OR`)
# 
# patient_meta$expt_date = ymd(patient_meta$expt_date)
# patient_meta$subject_id = factor(make.names(patient_meta$expt_date), levels = unique(make.names(patient_meta$expt_date)))
# 
# 
# patient_cell_meta = left_join(cell_meta_merged, patient_meta) 
# patient_cell_meta %<>% rename(cell_id = cell_inds)

# read in patient chart info

patient_meta = read_excel(path = '~/valiante_ih/raw-data/Demographic information_Request_HM.xlsx', 
                          range = 'A1:H53', 
                          col_names = T, 
                          col_types = c('date', 'guess', 'guess', 'guess', 'guess', 'guess', 'guess', 'guess'), 
                          sheet = 'Demographic information-paper') %>% 
  rename(expt_date = `Date of surgery`, age = `Age (Years)`, sex = Sex, seizure_duration = "Years of seizures history", drugs = "Antiepileptic drugs")

patient_meta$expt_date = ymd(patient_meta$expt_date)
patient_meta$subject_id = factor(make.names(patient_meta$expt_date), levels = unique(make.names(patient_meta$expt_date)))

# add a field to tell if subjects are unique by their date
patient_meta$unique_subject = T
ambiguous_subject_ids = patient_meta %>% group_by(subject_id) %>% add_count() %>% filter(n > 1) %>% select(subject_id) %>% unlist()
patient_meta[patient_meta$subject_id %in% ambiguous_subject_ids, 'unique_subject'] = F

homeira_meta_structured = read.csv(file = '~/valiante_ih/raw-data/homeira_demo_structured.csv') %>% 
  rename(age = Age..years., sex = Sex, seizure_duration = Years.of.seizure.history, drugs = Antiepileptic.drugs, diagnosis = Diagnosis)

patient_meta = left_join(patient_meta %>% select(-drugs), homeira_meta_structured, by = c('age', 'sex', 'seizure_duration'))

patient_cell_meta = left_join(cell_meta_merged, patient_meta, by = 'expt_date')
patient_cell_meta %<>% rename(cell_id = cell_inds)

updated_subject_metadata = patient_cell_meta %>% 
  group_by(subject_id) %>% 
  add_count() %>% 
  ungroup() %>% 
  select(-cell_id) %>% distinct(.keep_all = T) %>% 
  arrange(`drugs`) %>% 
  rename(cell_count = n, 
         notes = X__1) %>%
  select(subject_id, expt_date, cell_count, layer_name, recorder_name, diagnosis, age, sex, seizure_duration, drugs, unique_subject)


patient_cell_meta = left_join(cell_meta_merged, updated_subject_metadata, by = c('expt_date', 'recorder_name', 'layer_name'))
patient_cell_meta %<>% rename(cell_id = cell_inds)

total_cell_count_matrix_unique = updated_subject_metadata %>% 
  filter(unique_subject, !is.na(sex), !is.na(age), !is.na(diagnosis)) %>% 
  group_by(subject_id) %>% 
  summarize(total_cell_count = sum(cell_count)) %>% 
  ungroup()

write.csv(total_cell_count_matrix_unique, file = 'demographic_data_w_cell_counts.csv')



# get ephys data - try plotting rmp initially

ephys_sheets = c('Total5Homeira-Lastversion.xlsx', 'Total5Lihua-Lastversion.xlsx', 'Total2Homeira-Lastversion.xlsx', 'Total2Lihua-Lastversion.xlsx')

iv_ranges = c('A2:CD10', 'A2:BO18', 'A2:AC10', 'A2:BC18')

rmp_data = lapply(1:length(ephys_sheets), function(sheet_name){
  curr_path = paste0('~/valiante_ih/raw-data/', ephys_sheets[sheet_name])
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

  rin_df = rin %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_inds')
  colnames(rin_df) = c('cell_inds', 'rin')
  # sag_names = paste0('sag', sag_data[, 1] %>% unlist) %>% make.names()
  
  rmp_data = rmp_data %>% as.data.frame()
  rmp_data = Reduce(merge, list(rmp_data, sag_trans, sag_ratio_trans, ap_trans, rin_df))
  return(rmp_data)
}) %>% bind_rows()

rmp_data %<>% rename(cell_id = cell_inds)

joined_ephys_data = left_join(patient_cell_meta, rmp_data)

joined_ephys_data$cohort = 'KRI'

kri_valid_demo_cells = joined_ephys_data %>% filter(!is.na(age), !is.na(sex), unique_subject, diagnosis %in% c('TLE','Tumor')) %>% 
  distinct(subject_id, cell_id, .keep_all = T)


joined_ephys_data %>% ggplot(aes(x = age, y = aphw, color = layer_name, shape = sex)) +  geom_jitter() + facet_grid(~factor(recorder_name))


joined_ephys_data %>% 
  ggplot(aes(x = layer_name, y = sagamp.400, color = layer_name)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(alpha = .5) + 
  facet_grid(~factor(recorder_name)) + 
  # ylab('Sag amp (mV, at -400 pA)') +
  xlab('cortical layer') 

kri_valid_demo_cells %>% 
  ggplot(aes(x = age, y = sag.400)) + 
  geom_jitter(alpha = .7) + 
  geom_smooth(method = 'lm') + 
  facet_grid(~layer_name) + 
  xlab('Subject age') + ylab('Sag ratio (at -400 pA)')


library(lme4)

m1 = lmer(sag.400 ~ layer_name + (1|subject_id) + age + sex, data = kri_valid_demo_cells %>% mutate(age = age/20))
m2 = lmer(sag.400 ~ layer_name + (1|subject_id) + age, data = kri_valid_demo_cells  %>% mutate(age = age/20))

m1 = lm(sagamp.400 ~ layer_name  -1 , data = kri_valid_demo_cells)


anova(m1, m2)

library(MuMIn)

r.squaredGLMM(m1)

read_excel(path = '~/valiante_ih/raw-data/Totallayer5-Homeira-Lastversion.xlsx', sheet = 7, )









