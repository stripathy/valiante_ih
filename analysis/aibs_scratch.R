library(stringr)
library(lme4)
library(magrittr)
library(dplyr)
library(ggplot2)
library(stats)
# get ephys data from aibs

# download summary cell types data from aibs cell types portal
aibs_url = 'http://celltypes.brain-map.org/cell_types_specimen_details.csv'

aibs_ephys_meta = read.csv(url(aibs_url))

# get computed ephys features based on Shreejoy's custom scripts - these are mostly the precomputed features plus a couple extras
# based on the 2018 Fall release of data
aibs_ephys_orig =read.csv(file = 'raw-data/aibs_ephys/aibs_aggregated_ephys_plus_morpho_v10.csv')

aibs_ephys = aibs_ephys_orig %>% filter(species == 'Homo Sapiens')

aibs_ephys %<>% mutate(
                        rmp = vrest, 
                        rin = ri, 
                        apthr = threshold_v_long_square,
                        apamp = peak_v_long_square - threshold_v_long_square, 
                        ahpamp = threshold_v_long_square - fast_trough_v_long_square,
                        aphw = rheo_first_spike_hw * 1000, 
                        apvel = upstroke_downstroke_ratio_long_square,
                        rheo = threshold_i_long_square,
                        maxfreq = max_rate_long_square, 
                        adratio = hero_first_mean_adratio, 
                        # sag = 1 - sag, # defined by aibs as vm change / max vm deflection
                        sagamp.400 = (rmp -vm_for_sag) * sag,
                        fislope = f_i_curve_slope, 
                        avgisi = avg_isi, 
                        cvisi = hero_isi_cv, 
                        latency = latency * 1000 # convert to ms
)


aibs_ephys_meta_small = aibs_ephys_meta %>% select(contains('donor'), 'specimen__id', 'csl__normalized_depth') %>%
  rename_at(vars(contains('donor')), funs(sub('donor__', '', .))) %>%
  rename(specimen_id = specimen__id, 
         depth_norm = csl__normalized_depth, 
         donor_id = id, donor_name = name)

aibs_ephys_meta_small %<>% filter(species == 'Homo Sapiens')

aibs_ephys_meta_small$age = lapply(aibs_ephys_meta_small$age, function(age_str) strsplit(age_str[[1]] %>% as.character, split = ' ')[[1]][1] %>% as.numeric) %>% unlist


aibs_human_ephys = merge(aibs_ephys, aibs_ephys_meta_small, by = 'specimen_id')

aibs_human_ephys %<>% mutate(layer_name = paste0('L', structure_layer_name), 
                       recorder_name = 'AIBS',
                       cell_id = factor(specimen_id),
                       subject_id = donor_name, 
                       diagnosis = disease_state.x, 
                       cohort = 'AIBS')

aibs_human_ephys$diagnosis = plyr::mapvalues(aibs_human_ephys$diagnosis, from = c('epilepsy', 'tumor'), to = c('TLE', 'Tumor'))
# aibs_human_ephys$layer_name = plyr::mapvalues(aibs_human_ephys$layer_name, from = c('L2', 'L3'), to = c('L2.3', 'L2.3'))
aibs_human_ephys$sex = plyr::mapvalues(aibs_human_ephys$sex, from = c('Male', 'Female'), to = c('M', 'F'))


aibs_human_ephys %>% filter(dendrite_type == 'spiny') %>% 
  ggplot(aes(x = structure_layer_name, y = sagamp.400, color = structure_layer_name, shape = disease_state.x)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(size = 2, alpha = .5) +
  xlab('cortical layer') + 
  # ylab('Sag amp (mV, at max hyperpol step)') + 
  facet_wrap(~disease_state.x)

aibs_human_ephys %>% filter(dendrite_type == 'spiny') %>% 
  ggplot(aes(x = structure_layer_name, y = sagamp.400, color = structure_layer_name)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(size = 2, alpha = .5) +
  xlab('cortical layer') + ylab('Sag amp (mV, at max hyperpol step)') + facet_wrap(~structure_area_abbrev)

aibs_human_ephys %>% filter(dendrite_type == 'spiny') %>% 
  ggplot(aes(x = age, y = sag,  color = disease_state.x, group = structure_layer_name)) + 
  stat_smooth(method = 'lm', color = 'black') + 
  geom_point(size = 2, alpha = .6) + facet_wrap(~structure_layer_name) + 
  xlab('Patient age (years)') +  ylab('Sag amp (mV, at max hyperpol step)') 

      
regress_dataset = aibs_human_ephys %>% filter(dendrite_type == 'spiny')
m1 = lmer(sagamp.400 ~ structure_layer_name  + (1|subject_id) + age + sex + diagnosis - 1, data = regress_dataset %>% mutate(age = age/15))
m2 = lmer(sagamp.400 ~ structure_layer_name + (1|subject_id) + sex  - 1, data = regress_dataset)


m1 = lmer(sagamp.400 ~ age + (1|structure_layer_name) + (1|disease_state.x) + (1|subject_id) + (1|sex), data = aibs_human_ephys %>% filter(dendrite_type == 'spiny'), REML = F)
m2 = lmer(sagamp.400 ~ (1|structure_layer_name) + (1|disease_state.x) + (1|subject_id) + (1|sex) , data = aibs_human_ephys %>% filter(dendrite_type == 'spiny'), REML = F)

anova(m1, m2)

m1 = lmer(sagamp.400 ~ sex + (1|structure_layer_name) + (1|disease_state.x) + (1|subject_id) + (1|age), data = aibs_ephys %>% filter(dendrite_type == 'spiny'), REML = F)
m2 = lmer(sagamp.400 ~ (1|structure_layer_name) + (1|disease_state.x) + (1|subject_id) + (1|age) , data = aibs_ephys %>% filter(dendrite_type == 'spiny'), REML = F)

anova(m1, m2)

m1 = lmer(sagamp.400 ~ depth_norm + (1|subject_id), data = aibs_ephys %>% filter(dendrite_type == 'spiny'), REML = F)
m2 = lmer(sagamp.400 ~ (1|structure_layer_name) + (1|disease_state.x) + (1|subject_id) + (1|age) , data = aibs_ephys %>% filter(dendrite_type == 'spiny'), REML = F)



aibs_ephys_orig =read.csv(file = 'raw-data/aibs_ephys/aibs_aggregated_ephys_plus_morpho_v10.csv')

aibs_ephys_mouse = aibs_ephys_orig %>% filter(species == 'Mus musculus')

aibs_ephys_mouse %<>% mutate(
  rmp = vrest, 
  rin = ri, 
  apthr = threshold_v_long_square,
  apamp = peak_v_long_square - threshold_v_long_square, 
  ahpamp = threshold_v_long_square - fast_trough_v_long_square,
  aphw = rheo_first_spike_hw * 1000, 
  apvel = upstroke_downstroke_ratio_long_square,
  rheo = threshold_i_long_square,
  maxfreq = max_rate_long_square, 
  adratio = hero_first_mean_adratio, 
  # sag = 1 - sag, # defined by aibs as vm change / max vm deflection
  sagamp.400 = (rmp-vm_for_sag) * sag,
  fislope = f_i_curve_slope, 
  avgisi = avg_isi, 
  cvisi = hero_isi_cv, 
  latency = latency * 1000 # convert to ms
)


aibs_ephys_mouse %>% filter(dendrite_type == 'spiny') %>% 
  ggplot(aes(x = structure_layer_name, y = rin, color = structure_layer_name)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(size = 2, alpha = .5) +
  xlab('cortical layer') + ylab('Sag amp (mV, at max hyperpol step)')

aibs_ephys_meta_small$age = lapply(aibs_ephys_meta_small$age, function(age_str) strsplit(age_str[[1]] %>% as.character, split = ' ')[[1]][1] %>% as.numeric) %>% unlist


aibs_ephys = merge(aibs_ephys, aibs_ephys_meta_small, by = 'specimen_id')

aibs_ephys %>% filter(dendrite_type == 'spiny') %>% 
  ggplot(aes(x = age, y = rmp, color = structure_layer_name, shape = disease_state.x)) + 
  geom_jitter(size = 2, alpha = .5) + facet_wrap(~structure_layer_name)

