# try to do some basic clustering with Homeira and lihua's data
library(tidyr)
library(tidyverse)
df %>% drop_na

aibs_human_ephys %>% head

ephys_features = c('rmp', 'rin', 'apthr', 'apamp', 'aphw', 'apvel', 'rheo', 'maxfreq', 'adratio', 'fislope', 'avgisi', 'sag')


df =aibs_human_ephys %>% filter(dendrite_type %in% c('spiny')) 


df = df %>% select(name, ephys_features)

complete_aibs_human_ephys = df %>% drop_na
rownames(complete_aibs_human_ephys) = df %>% drop_na %>% pull(name)

complete_aibs_human_ephys = complete_aibs_human_ephys %>% select(-name)

library(umap)
aibs_human_ephys_umap = umap(complete_aibs_human_ephys)

aibs_human_ephys_umap_new = merge(aibs_human_ephys_umap$layout %>% as.data.frame() %>% tibble::rownames_to_column(var = 'name'), aibs_human_ephys, by = 'name')
k2 <- kmeans(complete_aibs_human_ephys, centers = 3, nstart = 25)
cluster_id_df = k2$cluster %>% as.data.frame() %>% tibble::rownames_to_column(var = 'name')
colnames(cluster_id_df) = c('name', 'cluster')
aibs_human_ephys_umap_new = merge(aibs_human_ephys_umap_new, cluster_id_df, by = 'name')
aibs_human_ephys_umap_new$cluster = aibs_human_ephys_umap_new$cluster %>% factor()


aibs_human_ephys_umap_new %>% filter(layer_name %in% c('L2', 'L3', 'L5')) %>% 
  ggplot(aes(x = V1, y = V2, color = layer_name, shape = cluster)) + 
  geom_point(size = 3, alpha = .7) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) 

aibs_human_ephys_umap_new %>% filter(V1 > 2, layer_name == 'L5', dendrite_type == 'spiny') %>% head

library(cluster)    # clustering algorithms

str(k2)

## try clustering kri data

joined_ephys_data


kri_ephys_features = c('rin', 'rmp', 'apamp', 'aphw', 'sagamp.400', 'inst_freq300', 'if_hz300', 'if_hz200', 'if_hz100', 'tau.300', 'sag.300', 'rebound.300')

complete_kri_human_ephys = joined_ephys_data %>% 
  filter(layer_name %in% c('L2.3', 'L3c')) %>% 
  dplyr::select(cell_id, kri_ephys_features, ) %>% drop_na %>% as.data.frame() %>% select(-cell_id)
rownames(complete_kri_human_ephys) = joined_ephys_data %>% 
  filter(layer_name %in% c('L2.3', 'L3c')) %>% 
  dplyr::select(c(kri_ephys_features, 'cell_id'), ) %>% drop_na %>% pull(cell_id)


complete_kri_human_ephys_umap = umap(complete_kri_human_ephys)

complete_kri_human_ephys_umap_new = merge(complete_kri_human_ephys_umap$layout %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_id'), joined_ephys_data, by = 'cell_id')
# k2 <- kmeans(complete_aibs_human_ephys, centers = 3, nstart = 25)
# cluster_id_df = k2$cluster %>% as.data.frame() %>% tibble::rownames_to_column(var = 'name')
# colnames(cluster_id_df) = c('name', 'cluster')
# aibs_human_ephys_umap_new = merge(aibs_human_ephys_umap_new, cluster_id_df, by = 'name')
# aibs_human_ephys_umap_new$cluster = aibs_human_ephys_umap_new$cluster %>% factor()



k2 <- kmeans(complete_kri_human_ephys, centers = 2, nstart = 25)
cluster_id_df = k2$cluster %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cell_id')
colnames(cluster_id_df) = c('cell_id', 'cluster')
complete_kri_human_ephys_umap_new = merge(complete_kri_human_ephys_umap_new, cluster_id_df, by = 'cell_id')
complete_kri_human_ephys_umap_new$cluster = complete_kri_human_ephys_umap_new$cluster %>% factor()

complete_kri_human_ephys_umap_new %>% ggplot(aes(x = V1, y = V2, color = layer_name, shape = cluster)) + geom_point(size = 2)

