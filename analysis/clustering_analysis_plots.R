library(ggbeeswarm)
library(ggrepel)

cell_patient_ephys_combined = read.csv('summary_tables/cell_patient_ephys_combined.csv')
cell_patient_ephys_combined$layer_name = plyr::mapvalues(cell_patient_ephys_combined$layer_name, from = levels(cell_patient_ephys_combined$layer_name), to = c("Layer 2/3", "Layer 3c", "Layer 5"))

cell_patient_ephys_combined %>% ggplot(aes(x = layer_name, fill = cell_type)) + 
  geom_bar() + facet_wrap(~recorder_name) + 
  xlab('Neocortex layer') + ylab('Cell count')

p1 = cell_patient_ephys_combined %>% filter(cell_type == 'Pyr') %>%
  filter(!res_center_freq %>% is.na) %>% 
  ggplot(aes(x = layer_name, fill = has_resonance)) + 
  geom_bar() + facet_wrap(~cell_type)

cell_patient_ephys_combined$cell_type = factor(cell_patient_ephys_combined$cell_type, levels = c('Pyr', 'Int'))
p2 = cell_patient_ephys_combined %>%
  filter(res_center_freq > 1) %>% filter(cell_type == 'Pyr') %>%
  ggplot(aes(x = layer_name, y = res_center_freq, color = layer_name)) +
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_quasirandom(alpha = .75, size = 2) + 
  facet_wrap(~cell_type) + 
  ylab('Resonance frequency (Hz)') + 
  xlab('Cortical layer')

plot_grid(p1, p2, nrow = 1)


library(umap)
# kri_extracted_features

kri_ephys_features = c('rin', 'rmp', 'apamp', 'ahpamp', 'aphw', 'apvel', 'sagamp',
                       #'adratio', 'first_isi', 'avgisi', 'cvisi',
                       'sag', 'fislope',  'latency', 'avg_rate', 'tau', 'rheo',
                       'apthr')
use_cell_types = c('Pyr')

#kri_ephys_features = c('aphw', 'apvel', 'adratio', 'fislope', 'ahpamp', 'avgisi', 'avg_rate')
#use_cell_types = c('Pyr', 'Int', NA)

# reassign big Layer3 cell to Layer 2/3 for consistency
cell_patient_ephys_combined[cell_patient_ephys_combined$cell_id == '19129004.abf', 'layer_name'] = 'Layer 2/3'

# add new field to data frame to indicate naming of cells in Fig 1 



complete_kri_human_ephys = cell_patient_ephys_combined %>% 
  filter(cell_type %in% use_cell_types) %>% 
  dplyr::select(cell_id, kri_ephys_features, ) %>% drop_na %>% as.data.frame() %>% select(-cell_id)
rownames(complete_kri_human_ephys) = cell_patient_ephys_combined %>% 
  filter(cell_type %in% use_cell_types) %>% 
  dplyr::select(c(kri_ephys_features, 'cell_id'), ) %>% drop_na %>% pull(cell_id)


complete_kri_human_ephys_umap = umap(complete_kri_human_ephys)

complete_kri_human_ephys_umap_new = merge(complete_kri_human_ephys_umap$layout %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_id'), 
                                          cell_patient_ephys_combined, by = 'cell_id')
# k2 <- kmeans(complete_aibs_human_ephys, centers = 3, nstart = 25)
# cluster_id_df = k2$cluster %>% as.data.frame() %>% tibble::rownames_to_column(var = 'name')
# colnames(cluster_id_df) = c('name', 'cluster')
# aibs_human_ephys_umap_new = merge(aibs_human_ephys_umap_new, cluster_id_df, by = 'name')
# aibs_human_ephys_umap_new$cluster = aibs_human_ephys_umap_new$cluster %>% factor()



k2 <- kmeans(complete_kri_human_ephys, centers = 4, nstart = 25)
cluster_id_df = k2$cluster %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_id')
colnames(cluster_id_df) = c('cell_id', 'cluster')
complete_kri_human_ephys_umap_new = merge(complete_kri_human_ephys_umap_new, cluster_id_df, by = 'cell_id')
complete_kri_human_ephys_umap_new$cluster = complete_kri_human_ephys_umap_new$cluster %>% factor()

complete_kri_human_ephys_umap_new$layer_name = factor(complete_kri_human_ephys_umap_new$layer_name, 
                                                      levels = c('Layer 5', 'Layer 2/3', 'Layer 3c'))
res_types = c('non-resonant (fR = 0 Hz)', 'weak (0 < fR < 2 Hz)', 'strong (fR > 2 Hz)')

complete_kri_human_ephys_umap_new$resonance_type = factor(complete_kri_human_ephys_umap_new$resonance_type, 
                                                      levels = c('non-resonant', 'weak', 'strong'))
# complete_kri_human_ephys_umap_new$resonance_type  = plyr::mapvalues(complete_kri_human_ephys_umap_new$resonance_type , 
#                                                                     from = levels(cell_patient_ephys_combined$resonance_type), to = res_types)

                                                      

p0 = complete_kri_human_ephys_umap_new  %>% arrange(layer_name) %>%
  ggplot(aes(x = V1, y = V2, color = layer_name)) + 
  scale_color_manual(values = c('red', 'blue', 'turquoise4'), name = 'Layer') + 
  geom_point(size = 2, alpha = .75) + 
  geom_text_repel(data = complete_kri_human_ephys_umap_new %>%
                    filter(cell_id %in% cells_w_morphology), aes(label = cell_id), 
                  nudge_y = -1, nudge_x = 1) +
  ylab('UMAP 2') + xlab('UMAP 1') + 
  theme(legend.position="top") 

p1 = complete_kri_human_ephys_umap_new %>% filter(has_resonance %in% c(T, F)) %>% arrange(resonance_type) %>%
  ggplot(aes(x = V1, y = V2, color = resonance_type)) + 
  scale_color_manual(values = c('lightgrey', 'black', 'red'), name = 'Resonance') + 
  #scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_jitter(size = 1.5, width = .25, height = .25, alpha = .5) + 
  #scale_color_discrete(name = "Resonant") + 
  #scale_alpha_manual(guide='none', values = list(a = 0.2, point = 1))
  # geom_text_repel(data = complete_kri_human_ephys_umap_new %>% 
  #                   filter(cell_id %in% cells_w_morphology), aes(label = cell_id)) +
  ylab('UMAP 2') + xlab('UMAP 1') + 
  theme(legend.position="right")

complete_kri_human_ephys_umap_new = complete_kri_human_ephys_umap_new %>% mutate(has_burst_algo_hero = ((1/first_isi) > 75 & !(is.na(first_isi))))
complete_kri_human_ephys_umap_new = complete_kri_human_ephys_umap_new %>% mutate(has_burst_def = case_when(
  (has_burst_algo == T)~ 'strong',
  (has_burst_algo == F & has_burst_algo_hero == T) ~ 'weak',
  has_burst_algo_hero == F ~ 'none'
))
complete_kri_human_ephys_umap_new$has_burst_def = factor(complete_kri_human_ephys_umap_new$has_burst_def, levels = c('none', 'weak', 'strong'))

p2 = complete_kri_human_ephys_umap_new %>% 
  arrange(has_burst_def) %>%
  ggplot(aes(x = V1, y = V2, color = has_burst_def)) + 
  scale_color_manual(values = c('lightgrey', 'black', 'red'), name = "Bursting") + 
  geom_jitter(size = 1.5, width = .25, height = .25, alpha = .5) + 
  # geom_text_repel(data = complete_kri_human_ephys_umap_new %>% 
  #                   filter(cell_id %in% cells_w_morphology), aes(label = cell_id)) +
  ylab('UMAP 2') + xlab('UMAP 1') + 
  theme(legend.position="right") #+ scale_color_discrete(name = "Bursting")
# 
# p2 = complete_kri_human_ephys_umap_new %>% filter(has_burst_algo_hero %in% c(T, F)) %>%
#   arrange(has_burst_algo) %>%
#   ggplot(aes(x = V1, y = V2, color = has_burst_algo_hero)) + 
#   scale_color_manual(values = c('lightgrey', 'red'), name = "Bursting") + 
#   geom_jitter(size = 1.5, width = .25, height = .25, alpha = .5) + 
#   # geom_text_repel(data = complete_kri_human_ephys_umap_new %>% 
#   #                   filter(cell_id %in% cells_w_morphology), aes(label = cell_id)) +
#   ylab('UMAP 2') + xlab('UMAP 1') + 
#   theme(legend.position="right") #+ scale_color_discrete(name = "Bursting")

p3 = complete_kri_human_ephys_umap_new %>% 
  #arrange(has_burst_algo) %>%
  ggplot(aes(x = V1, y = V2, color = rin)) + 
  scale_color_steps(name = 'Rin (MOhm)', breaks = c(25, 50, 75, 100, 125)) + 
  geom_jitter(size = 1.5, alpha = .5, width = .25, height = .25) + 
  # geom_text_repel(data = complete_kri_human_ephys_umap_new %>% 
  #                   filter(cell_id %in% cells_w_morphology), aes(label = cell_id)) +
  ylab('UMAP 2') + xlab('UMAP 1') + 
  theme(legend.position="right") #+ scale_color_discrete(name = "Bursting")

right_plot = plot_grid(p3, p2, p1, ncol = 1, align = 'hv')
full_plot = plot_grid(p0, right_plot, nrow = 1, rel_widths = c(1, .75)) 
full_plot

ggsave(file = 'figures/Fig7_clustering_r_gen.pdf', plot = full_plot, width = 8, height = 5, device = "pdf")

#plot_grid(p0, p1, p2, nrow = 1)



