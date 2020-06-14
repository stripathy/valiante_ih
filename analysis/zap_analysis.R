temp_res_data =read.csv('~/Documents/GitHub/valiante_lab_abf_process/output_files/fred_zap_analysis.csv')

remove_zap_abf_list = c('15105068.abf', '15127064.abf', '14515402.abf') # a list of files that are suprathreshold
temp_res_data = temp_res_data %>% filter(! abf_name %in% remove_zap_abf_list)


lihua_file_path = '/Users/stripathy/Downloads/homeira_lihua_zap_files/Lihua/'

meta_path = paste0(lihua_file_path, 'L5-Resonance-Intrinsic-Lihua-April 18-2020.xlsx')

abf_sheet = read_excel(meta_path, n_max = 54)
abf_sheet$layer_name = 'L5'
abf_sheet$recorder_name = 'Lihua'

update_col_names = c('cell_id', '1', 'intrinsic_neg', '2', 'intrinsic_rmp_zd', '3', 'intrinsic_neg_zd', '11',
                     'zap_sub_rmp', '4', 'zap_sub_neg', '5', 'zap_sup_rmp', '6', 
                     'zap_sup_neg', '7', 'zap_sub_rmp_zd', '8', 'zap_sub_neg_zd', '9', 'zap_sup_rmp_zd', 
                    '10', 'tag', '13', 'blah', 'layer_name', 'recorder_name')

meta_path = paste0(lihua_file_path, 'L23-Resonance-Intrinsic-Lihua-April 18-2020 - Copy.xlsx')
abf_sheet_l23 = read_excel(meta_path, n_max = 38)
abf_sheet_l23$layer_name = 'L2.3'
abf_sheet_l23$recorder_name = 'Lihua'

# abf_sheet_l23 %<>% mutate(`ZAP-Suprathreshold...23`= as.character(`ZAP-Suprathreshold...23`))

lihua_sheet = bind_rows(abf_sheet, abf_sheet_l23)
colnames(lihua_sheet) = update_col_names

homeira_l3_path = '/Users/stripathy/Downloads/homeira_lihua_zap_files/homeira_metadata_files/L3/L3-Resonance-Intrinsic-Homeira-April 16-2020.xlsx'

meta_path = paste0(homeira_l3_path)

homeira_l3_sheet = read_excel(meta_path, n_max = 26)

update_col_names = c('cell_id', 'zap_sub_rmp', '1', 'zap_sub_neg', '5', 'zap_sup_rmp', '6', 
                      '7', 'tag')
colnames(homeira_l3_sheet) = update_col_names

homeira_l3_sheet %<>% select('cell_id', 'zap_sub_rmp', 'zap_sub_neg','zap_sup_rmp',  
             'tag')
homeira_l3_sheet$layer_name = 'L3c'
homeira_l3_sheet$recorder_name = 'Homeira'


homeira_l5_path = '/Users/stripathy/Downloads/homeira_lihua_zap_files/homeira_metadata_files/L5/L5-Resonance-Intrinsic-Homeira-April 16-2020.xlsx'

meta_path = paste0(homeira_l5_path)

homeira_l5_sheet = read_excel(meta_path, n_max = 41)

update_col_names = c('cell_id', 'zap_sub_rmp', '1', 'zap_sub_neg', '5', 'zap_sup_rmp', '6', 
                     '7', 'zap_sub_rmp_zd', '8', '9', 'tag')
colnames(homeira_l5_sheet) = update_col_names

homeira_l5_sheet %<>% select('cell_id', 'zap_sub_rmp', 'zap_sub_neg','zap_sup_rmp',  
             'zap_sub_rmp_zd', 'tag')
homeira_l5_sheet$layer_name = 'L5'
homeira_l5_sheet$recorder_name = 'Homeira'



homeira_l23_path = '/Users/stripathy/Downloads/homeira_lihua_zap_files/homeira_metadata_files/L23/L23-Resonance-Intrinsic-Homeira-April 16-2020.xlsx'

meta_path = paste0(homeira_l23_path)

homeira_l23_sheet = read_excel(meta_path, n_max = 44)

update_col_names = c('cell_id', 'zap_sub_rmp', '4', 'zap_sub_neg', '5', 'zap_sup_rmp', '6', 
                     'zap_sup_neg', '7', 'zap_sub_rmp_zd', '8', 'zap_sub_neg_zd', '9', 'zap_sup_rmp_zd', '10', '11', 'tag'
                     )
colnames(homeira_l23_sheet) = update_col_names

homeira_l23_sheet %<>% select('cell_id', 'zap_sub_rmp', 'zap_sub_neg','zap_sup_rmp',  
                              'zap_sup_neg', 'zap_sub_rmp_zd',  'zap_sub_neg_zd', 'zap_sup_rmp_zd','tag')
homeira_l23_sheet$layer_name = 'L2.3'
homeira_l23_sheet$recorder_name = 'Homeira'

data_frame_list = list(lihua_sheet, homeira_l23_sheet, homeira_l3_sheet, homeira_l5_sheet) 

abf_sheet = bind_rows(data_frame_list)

keep_columns = c('cell_id', 'intrinsic_neg',  'intrinsic_rmp_zd',  'intrinsic_neg_zd', 
                 'zap_sub_rmp', 'zap_sub_neg', 'zap_sup_rmp', 
                 'zap_sup_neg', 'zap_sub_rmp_zd',  'zap_sub_neg_zd',  "zap_sup_rmp_zd")

abf_sheet_wide = abf_sheet
colnames(abf_sheet_wide) = make.names(colnames(abf_sheet_wide), unique = T)
# now do the processing of the sheet to make it long prior to associating with zap processed data from python and fred's script

abf_sheet_temp = abf_sheet_wide %>% select(keep_columns, layer_name, recorder_name, -intrinsic_neg, -intrinsic_rmp_zd, -intrinsic_neg_zd)

abf_sheet_long = abf_sheet_temp  %>%
  pivot_longer(cols = zap_sub_rmp:zap_sup_rmp_zd,
               names_to = 'stim_type', values_to = 'abf_name') %>% filter(abf_name != '0', abf_name != 'NA') %>% 
  mutate(cell_id = paste0(trimws(cell_id), '.abf'))

temp_df = data.frame()
for (ind in 1:nrow(abf_sheet_long)){
  temp_zap_file_range = abf_sheet_long[ind, 'abf_name'] %>% unlist
  if ((temp_zap_file_range != 0) & (!is.na(temp_zap_file_range))){
    file_list = get_abf_file_list(temp_zap_file_range)
  }
  temp_df = bind_rows(temp_df, data.frame('cell_id' = abf_sheet_long[ind, 'cell_id'],  
                                          'abf_name' = file_list, 'stim_type' = abf_sheet_long[ind, 'stim_type'], 
                                          'layer_name' = abf_sheet_long[ind, 'layer_name'],
                                          'recorder_name' = abf_sheet_long[ind, 'recorder_name']))
}
abf_sheet_long_unraveled = temp_df %>% distinct(abf_name, .keep_all = T)

abf_sheet_long_first = abf_sheet_long %>% mutate(abf_name = get_first_abf_file(abf_name))


# abf_sheet_long = abf_sheet_temp %>% select(cell_id, zap_sub_rmp) %>%
#   pivot_longer(cols = zap_sub_rmp, 
#                names_to = 'stim_type', values_to = 'abf_name') %>% filter(abf_name != '0.abf', abf_name != 'NA.abf')
# 
# abf_sheet_long = abf_sheet_temp %>% select(-intrinsic_neg, -intrinsic_rmp_zd, -intrinsic_neg_zd) %>%
#   pivot_longer(cols = zap_sub_rmp:zap_sup_rmp_zd,
#                names_to = 'stim_type', values_to = 'abf_name') %>% filter(abf_name != '0.abf', abf_name != 'NA.abf')

# combine fred processed data with homeira provided zap intrinsic file mapping sheet
zap_data_full = left_join(temp_res_data %>% distinct(abf_name, .keep_all = T) , abf_sheet_long_unraveled %>% distinct(abf_name, .keep_all = T), by = 'abf_name')

zap_data = left_join(temp_res_data %>% distinct(abf_name, .keep_all = T) , abf_sheet_long_first %>% distinct(abf_name, .keep_all = T), by = 'abf_name')

zap_data = bind_rows(zap_data, zap_data_full) %>% distinct(abf_name, .keep_all = T)

# # because we use multiple protocols for zap and also to filter out files not containing valid cell ids
# zap_data_w_intrinsic_id_final = zap_data %>% 
#   filter(cell_id != '0.abf', cell_id != 'NA.abf') %>% 
#   arrange(cell_id, file_time) %>% 
#   distinct(cell_id, .keep_all = T) %>% select(-X) %>%
#   rename(res_center_freq = center_freq,
#          res_3dB_freq = X3dB_freq, 
#          res_sharpness = res_sharpness) %>%
#   mutate(has_resonance = res_center_freq > 1)

# keep only zap recordings at rmp and also to filter out files not containing valid cell ids
zap_data_w_intrinsic_id_final = zap_data %>% 
  filter(cell_id != '0.abf', cell_id != 'NA.abf', stim_type == 'zap_sub_rmp') %>% 
  arrange(cell_id, file_time) %>% 
  distinct(cell_id, .keep_all = T) %>% select(-X) %>%
  rename(res_center_freq = center_freq,
         res_3dB_freq = X3dB_freq, 
         res_sharpness = res_sharpness) %>%
  mutate(has_resonance = res_center_freq > 0) %>%
  mutate(resonance_type = case_when(res_center_freq == 0 ~ 'non-resonant',
                                   (res_center_freq > 0 & res_center_freq < 2) ~ 'weak',
                                   res_center_freq > 2 ~ 'strong'
                                   )) %>%
  mutate(resonance_type = factor(resonance_type, levels = c('non-resonant', 'weak', 'strong')))

# write csv of zap_data_w_intrinsic_id_final
# write.csv(zap_data_w_intrinsic_id_final, 'summary_tables/zap_data_final_summarized.csv')


# now merge zap data with cell meta based on my prior processing in R

zap_data_w_cell_meta = merge(cell_meta, zap_data_w_intrinsic_id_final %>% select(-layer_name), by = c('cell_id', 'recorder_name') )

zap_data_merged_final = bind_rows(zap_data_w_cell_meta %>% mutate(in_big_dataset = T), 
                                  zap_data_w_intrinsic_id_final %>% filter(! cell_id %in% zap_data_w_cell_meta$cell_id) %>% mutate(in_big_dataset = F)) %>%
  rename(zap_abf_name = abf_name) %>% arrange(file_time)

# new interneurons

new_interneuron_list = c('14d02169.abf', '14d160190.abf', '2019_11_04_0083.abf')

zap_data_merged_final[zap_data_merged_final$cell_id %in% new_interneuron_list, 'cell_type'] = 'Int'

new_pyr_list = c('19228000.abf', '14n03301.abf', '19228029.abf', '19228058.abf', 
                 '2019_11_26_0006.abf', '2020_01_27_0014.abf', '2020_01_27_0045.abf', 
                 '19129064.abf', '19129072.abf', '2019_11_26_0037.abf', '2019_11_28_0010.abf', '2019_11_28_0032.abf',
                 '2020_01_06_0101.abf', '2020_01_06_0108.abf', '2020_03_02_0019.abf')

zap_data_merged_final[zap_data_merged_final$cell_id %in% new_pyr_list, 'cell_type'] = 'Pyr'

zap_data_merged_final %>% 
  group_by(cell_type, layer_name) %>% 
  summarize(strong_res_pct = sum((resonance_type == 'strong'), na.rm = T) / n(), 
            weak_res_pct = sum((resonance_type == 'weak'), na.rm = T)/n())

# write csv of zap_data_w_intrinsic_id_final
write.csv(zap_data_merged_final, 'summary_tables/zap_data_final_summarized.csv')

library(ggbeeswarm)

p1 = zap_data_merged_final %>% filter(cell_type == 'Pyr') %>%
  filter(!res_center_freq %>% is.na) %>% 
  ggplot(aes(x = layer_name, fill = has_resonance)) + 
  geom_bar()

cell_patient_ephys_combined$cell_type = factor(cell_patient_ephys_combined$cell_type, levels = c('Pyr', 'Int'))
p2 = zap_data_merged_final %>% filter(cell_type == 'Pyr') %>% 
  ggplot(aes(x = layer_name, y = res_center_freq, color = layer_name)) +
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_quasirandom(alpha = .5, size = 2) + 
  #facet_wrap(~cell_type) + 
  ylab('Resonance frequency (Hz)') + 
  #xlab('Cortical layer') + 
  theme(legend.position = "none", axis.title.x=element_blank()) +
  theme(text=element_text(size=12))

p3 = zap_data_merged_final %>% filter(cell_type == 'Pyr') %>%
  ggplot(aes(x = layer_name, y = res_3dB_freq, color = layer_name)) +
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_quasirandom(alpha = .5, size = 2) + 
  #facet_wrap(~cell_type) + 
  ylab('3dB cutoff (Hz)') + 
  #xlab('Cortical layer') + 
  theme(legend.position = "none", axis.title.x=element_blank()) +
  theme(text=element_text(size=12))


res_freq_plot = plot_grid(p2, p3, nrow = 1)
ggsave('figures/res_freq_plot.pdf', plot = res_freq_plot, width = 4, height = 2.5, units = 'in', scale = 1, useDingbats=FALSE)

l23 = zap_data_merged_final %>% filter(cell_type == 'Pyr', layer_name == 'L2.3') %>% pull(res_3dB_freq)
l3c = zap_data_merged_final %>% filter(cell_type == 'Pyr', layer_name == 'L3c') %>% pull(res_3dB_freq)
l5 = zap_data_merged_final %>% filter(cell_type == 'Pyr', layer_name == 'L5') %>% pull(res_3dB_freq)

wilcox.test(l5, l3c)
wilcox.test(l23, l5)


cell_patient_ephys_combined$cell_type = factor(cell_patient_ephys_combined$cell_type, levels = c('Pyr', 'Int'))
p2 = zap_data_merged_final %>% filter(cell_type == 'Int') %>% 
  ggplot(aes(x = layer_name, y = res_center_freq, color = layer_name)) +
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_quasirandom(alpha = .5, size = 2) + 
  #facet_wrap(~cell_type) + 
  ylab('Resonance frequency (Hz)') + 
  #xlab('Cortical layer') + 
  theme(legend.position = "none", axis.title.x=element_blank()) +
  theme(text=element_text(size=12))

p3 = zap_data_merged_final %>% filter(cell_type == 'Int') %>%
  ggplot(aes(x = layer_name, y = res_3dB_freq, color = layer_name)) +
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_quasirandom(alpha = .5, size = 2) + 
  #facet_wrap(~cell_type) + 
  ylab('3dB cutoff (Hz)') + 
  #xlab('Cortical layer') + 
  theme(legend.position = "none", axis.title.x=element_blank()) +
  theme(text=element_text(size=12))


res_freq_plot = plot_grid(p2, p3, nrow = 1)
ggsave('figures/int_res_freq_plot.pdf', plot = res_freq_plot, width = 4, height = 3, units = 'in', scale = 1, useDingbats=FALSE)




plot_grid(p1, p2, nrow = 1)

# merge zap data with some representation of intrinsic data in R
zap_data_joined_full = zap_data_w_cell_meta %>%
  arrange(cell_id, file_time) %>%
  #filter(stim_type == 'zap_sub_rmp') %>%
  distinct(cell_id, .keep_all = T)

zap_freq_stats = zap_data_joined_full %>% filter(!is.na(res_center_freq)) %>% 
  group_by(cell_type, layer_name, recorder_name) %>% 
  summarize(resonant_cell_ratio = sum(res_center_freq > 0) / n(), cell_count = n()) %>% 
  arrange(cell_type, layer_name)

zap_data_joined_full %>% 
  filter(!is.na(res_center_freq))%>% 
  #filter(layer_name == 'L3c') %>%
  filter(res_center_freq == 0) %>%
  ggplot(aes(x = V1, y = V2, shape = cluster, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_jitter(size = 2.5, alpha = .75)
# geom_text_repel(data = complete_kri_human_ephys_umap_new %>% 
#                   filter(cell_id %in% cells_w_morphology), aes(label = cell_id)) 

complete_kri_human_ephys_umap_new_w_resonance = left_join(complete_kri_human_ephys_umap_new, zap_data_w_intrinsic_id_final %>% select(-layer_name), by = c('cell_id', 'recorder_name') )

complete_kri_human_ephys_umap_new_w_resonance = complete_kri_human_ephys_umap_new_w_resonance %>% mutate(has_resonance = case_when(
    is.na(center_freq)  ~ "not measured",
    center_freq == 0 ~ "non-resonant",
    center_freq > 0 ~ "resonant"
    #center_freq < 1 & center_freq > 0 ~ "low resonance",
    #center_freq > 1  ~ "higher resonance",
))


p1 = complete_kri_human_ephys_umap_new_w_resonance %>% filter(has_resonance != 'not measured') %>% ggplot(aes(x = V1, y = V2, color = has_resonance)) + 
  #scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_jitter(size = 2.5, alpha = .75, width = .5, height = .5)

p2 = complete_kri_human_ephys_umap_new_w_resonance  %>% ggplot(aes(x = V1, y = V2, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_jitter(size = 2.5, alpha = .75, width = .5, height = .5)

plot_grid(p2, p1, ncol = 2)



# come up with list of cells for plotting
# L2.3
'14317333.abf' # most resonant L2.3 cell
'19122003.abf' # small layer 2.3 cell
'19129004.abf' # upper layer 3 cell

# L3c
'2020_03_02_0019.abf' # most resonant L3c cell
'2020_01_06_0048.abf' # weakly resonance, L3c cell

# L5
'15420000.abf' # most resonant L5 cell
'19129015.abf' # small L5 morphology, not resonant


homeira_excel_files = glob.glob(lihua_file_path, recursive= False)
print(homeira_excel_files)

for excel_file_name in homeira_excel_files:
  excel_sheet = read_excel(excel_file_name)
abf_sheet = abf_sheet.append(excel_sheet, ignore_index=True) 
len(abf_sheet)

get_first_abf_file = function(input_string){
  return(input_string %>% str_split(., '-|to') %>% lapply(`[[`, 1) %>% str_trim() %>% paste0(., '.abf'))
}

get_abf_file_list = function(input_string){
  
  first_file = input_string %>% str_split(., '-|to') %>% unlist %>% first() %>% str_trim()
  if (first_file == input_string){
    return_string = paste0(input_string[[1]], '.abf')
    return(input_string)
  }
  last_file = input_string %>% str_split(., '-|to') %>% lapply(`[[`, -1) %>% str_trim()
  
  if (str_detect(first_file, '\\D')){ # if not entirely numeric
    first_file_num_str = first_file %>% str_split(., '\\D') %>% unlist %>% last() %>% str_trim()
    second_file_num_str = last_file %>% str_split(., '\\D') %>% unlist %>% last() %>% str_trim()
    
    num_files = as.integer(second_file_num_str) - as.integer(first_file_num_str)
    
    if ((num_files > 30 ) || (num_files < 0)){
      return(first_file)
    }
    add_ints = seq(first_file_num_str, second_file_num_str, by = 1)
    if (nchar(add_ints[1]) < nchar(first_file_num_str)){
      add_ints = str_pad(add_ints, width = nchar(first_file_num_str), "0", side = "left")
    }
    first_file_base = substr(first_file,1, (nchar(first_file) - nchar(first_file_num_str)))
    file_list = paste0(first_file_base, add_ints, '.abf')
  }else{
    first_file_int = as.integer(first_file)
    last_file_int = as.integer(last_file)
    
    num_files = last_file_int - first_file_int
    
    if ((num_files > 30 ) || (num_files < 0)){
      return(first_file)
    }
    
    file_list = seq(first_file_int, last_file_int, by = 1) %>% paste0(., '.abf')
  }
  return(file_list)
}

missing_zap_file_ids