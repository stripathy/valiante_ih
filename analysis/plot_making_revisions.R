library(lme4)
library(broom)
library(MuMIn)
library(cowplot)

theme_set(theme_cowplot(14) )

# read summarized aibs's data


aibs_human_pyr_ephys = aibs_human_ephys %>% filter(dendrite_type == 'spiny', 
                                                   structure_layer_name %in% c(2, 3, 5)) %>%
  rename(layer = structure_layer_name)
aibs_human_pyr_ephys$diagnosis  = plyr::mapvalues(aibs_human_pyr_ephys$diagnosis, c('TLE', 'Tumor'), c('Epilepsy', 'Tumor'))
aibs_human_pyr_ephys$cell_type  = aibs_human_pyr_ephys$dendrite_type



### define a final data frame that has just ephys data from cells with "complete" demographic data from patients
cell_patient_ephys_combined = read.csv('summary_tables/cell_patient_ephys_combined.csv')

kri_valid_demo_cells = cell_patient_ephys_combined %>% filter(!is.na(age), 
                                                              !is.na(sex), 
                                                              unique_subject,
                                                              # diagnosis %in% c('TLE', 'Tumor')
) %>% 
  distinct(subject_id, cell_id, .keep_all = T) %>% select(-has_burst)



kri_valid_demo_cells$layer = plyr::mapvalues(kri_valid_demo_cells$layer_name, from = c('L2.3', 'L5'), to = c('Layer 2/3', 'Layer 5'))
kri_valid_demo_cells$source = 'Krembil (Toronto)'

aibs_human_pyr_ephys$layer = plyr::mapvalues(aibs_human_pyr_ephys$layer, from = c('2', '3', '5'), to = c('Layer 2', 'Layer 3', 'Layer 5'))
aibs_human_pyr_ephys$source = 'Allen (Seattle)'

comb_plot_dataset = bind_rows(kri_valid_demo_cells %>% filter(cell_type == 'Pyr'), aibs_human_pyr_ephys %>% filter(layer %in% c('Layer 2', 'Layer 3', 'Layer 5')))
comb_plot_dataset$source = factor(comb_plot_dataset$source, levels = c('Krembil (Toronto)', 'Allen (Seattle)'))
comb_plot_dataset$diagnosis = factor(comb_plot_dataset$diagnosis, levels = c('Epilepsy', 'Tumor', 'Other'))


### supplemental figure showing that sag amplitudes are correlated with sag ratio

sag_ratio_vs_amp = comb_plot_dataset %>% filter(source == 'Krembil (Toronto)') %>% 
  ggplot(aes(x = sagamp.400, y = sag.400, group = layer)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  # scale_color_manual(values = c('black', 'red', 'purple')) + 
  facet_grid(~layer, scales = "free_x") +  
  xlab('Sag amplitude (mV)') + ylab('Sag ratio')

ggsave('figures/sag_ratio_vs_amp.pdf', plot = sag_ratio_vs_amp, width = 6, height = 3, units = 'in', scale = 1, useDingbats=FALSE)

sources = c('Krembil (Toronto)', 'Allen (Seattle)')
layers = c('Layer 2/3', 'Layer 5')

sag_formula = 'sag ~ (1|subject_id) + age'
sag_formula_no_demo = 'sag ~ (1|subject_id)'
cnt = 1
ret_df = data.frame(source=character(),
                    layer=character(), 
                    beta_age=double(), 
                    std.err_age=double(),
                    age_corr_val = double(),
                    anova_p = double()) 

for (i in 1:length(sources)){
  for (j in 1:length(layers)){
    curr_cohort = sources[i]
    curr_layer = layers[j]
    
    if ((curr_cohort == "Allen (Seattle)") & (curr_layer  == "Layer 2/3")){
      use_layers = c('Layer 2', 'Layer 3')
    } else{
      use_layers = curr_layer
    }
    
    lmer_mod = lmer(sag_formula, 
                    data = comb_plot_dataset %>% filter(source == curr_cohort, layer %in% use_layers))
    lmer_mod_no_demo = lmer(sag_formula_no_demo, 
                            data = comb_plot_dataset %>% filter(source == curr_cohort, layer %in% use_layers))
    
    anova_pvalue_df = anova(lmer_mod,lmer_mod_no_demo) %>% tidy()
    anova_pvalue = anova_pvalue_df[2, 'p.value']
    
    lmer_df = lmer_mod %>% tidy
    beta_age = lmer_df[2, 'estimate']
    sd_age = lmer_df[2, 'std.error']
    
    # uses the MuMIn toolbox and method described here: https://stats.stackexchange.com/questions/327179/find-repeated-measure-correlation-coefficient-using-linear-mixed-model
    fixed_effects_corr_value = sqrt(r.squaredGLMM(lmer_mod)[1]) 
    
    
    ret_df = rbind(ret_df, data.frame(source = curr_cohort %>% as.character(), 
                                      layer = curr_layer %>% as.character(), 
                                      beta_age = beta_age %>% as.double(), 
                                      std.err_age = sd_age %>% as.double(), 
                                      age_corr_val = fixed_effects_corr_value %>% as.double(),
                                      anova_p = anova_pvalue %>% as.double()
    )
    )
    
    cnt = cnt + 1
    
  }
}



new_df = comb_plot_dataset 
new_df$colors = 'grey35'
new_df$sizes = 2
new_df[new_df$cell_id == '15o08038.abf', 'colors'] = 'blue'
#new_df[new_df$cell_id == '13d03008.abf', 'colors'] = 'red'
new_df[new_df$cell_id == '19129024.abf', 'colors'] = 'red'
#new_df[new_df$cell_id %in% c('13d03008.abf','15o08038.abf') , 'sizes'] = 6
new_df[new_df$cell_id %in% c('15o08038.abf','19129024.abf') , 'sizes'] = 6

krembil_xlim = new_df %>% filter(source == 'Krembil (Toronto)', layer == 'Layer 5') %>% pull(age) %>% range
krembil_xlim[1] = krembil_xlim[1] - 1
krembil_xlim[2] = krembil_xlim[2] + 1

allen_xlim = new_df %>% filter(source == 'Allen (Seattle)') %>% pull(age) %>% range
allen_xlim[1] = allen_xlim[1] - 1
allen_xlim[2] = allen_xlim[2] + 1

sag_ratio_vs_age_l5_krembil = new_df %>% filter(source == 'Krembil (Toronto)', layer == 'Layer 5') %>% 
  ggplot(aes(x = age, y = sag, group = layer, color = colors, size = sizes)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_identity() + 
  scale_size_identity() + 
  xlim(krembil_xlim) + 
  ylim(c(0, .325)) + 
  xlab('Age (years)') + ylab('Layer 5, Sag ratio')

sag_ratio_vs_age_l23_krembil = new_df %>% filter(source == 'Krembil (Toronto)', layer == 'Layer 2/3') %>% 
  ggplot(aes(x = age, y = sag, group = layer, color = colors, size = sizes)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_identity() + 
  scale_size_identity() + 
  xlim(krembil_xlim) + 
  ylim(c(0, .325)) + 
  xlab('Age (years)') + ylab('Layer 2/3, Sag ratio')

sag_ratio_vs_age_l5_allen = new_df %>% filter(source == 'Allen (Seattle)', layer == 'Layer 5') %>% 
  ggplot(aes(x = age, y = sag, group = layer, color = colors, size = sizes)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_identity() +
  scale_size_identity() + 
  ylim(c(0, .325)) + 
  # xlim(allen_xlim) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  xlab('Age (years)') + ylab('Layer 5, Sag ratio')

sag_ratio_vs_age_l23_allen = new_df %>% filter(source == 'Allen (Seattle)', layer %in% c('Layer 2', 'Layer 3')) %>% 
  ggplot(aes(x = age, y = sag, color = colors, size = sizes)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_identity() + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_size_identity() + 
  ylim(c(0, .325)) + 
  # xlim(allen_xlim) +
  xlab('Age (years)') + ylab('Layer 2/3, Sag ratio')

sag_ratio_all_panels = plot_grid(sag_ratio_vs_age_l5_krembil, sag_ratio_vs_age_l5_allen, sag_ratio_vs_age_l23_krembil, sag_ratio_vs_age_l23_allen, ncol = 2, align = 'v')

ggsave('figures/sag_ratio_all_panels.pdf', plot = sag_ratio_all_panels, width = 5, height = 6, units = 'in', scale = 1, useDingbats=FALSE)
