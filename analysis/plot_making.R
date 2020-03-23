
aibs_human_pyr_ephys = aibs_human_ephys %>% filter(dendrite_type == 'spiny', 
                                                   structure_layer_name %in% c(2, 3, 5)) %>%
                                            rename(layer = structure_layer_name)
aibs_human_pyr_ephys$diagnosis  = plyr::mapvalues(aibs_human_pyr_ephys$diagnosis, c('TLE', 'Tumor'), c('Epilepsy', 'Tumor'))


# aibs_human_pyr_ephys = aibs_human_ephys %>% filter(dendrite_type == 'spiny') %>% rename(layer = structure_layer_name)
aibs_age_sd = aibs_human_pyr_ephys %>% distinct(subject_id, .keep_all = T) %>% pull(age) %>% sd

aibs_lmer_model = lmer(sag ~ layer + (1|subject_id) + age + sex - 1, 
                       data = aibs_human_pyr_ephys %>% mutate(age = age/25) )


library(broom)
aibs_model_df = tidy(aibs_lmer_model) %>% filter(group == 'fixed') 

aibs_model_df$term = factor(aibs_model_df$term, levels = c('age', 'sexM', 'layer2',  'layer3', 'layer5'))
aibs_model_df$var_type = factor(c('layer', 'layer', 'layer', 'demo', 'demo'), levels = c('layer', 'demo'))
aibs_model_df$source = 'Allen'

kri_valid_demo_cells = joined_ephys_data %>% filter(!is.na(age), !is.na(sex), unique_subject) %>% rename(layer = layer_name)
kri_valid_demo_cells$diagnosis  = plyr::mapvalues(kri_valid_demo_cells$diagnosis, c('TLE', 'Tumor', NA, ''), c('Epilepsy', 'Tumor', 'Other', 'Other'))
kri_valid_demo_cells$diagnosis[is.na(kri_valid_demo_cells$diagnosis)] = 'Other'

kri_valid_demo_cells$diagnosis = factor(kri_valid_demo_cells$diagnosis, levels = c('Epilepsy', 'Tumor', 'Other'))

kri_age_sd = kri_valid_demo_cells %>% distinct(subject_id, .keep_all = T) %>% pull(age) %>% sd

kri_lmer_model = lmer(sag.400  ~ layer + (1|subject_id) + age + sex - 1 , data = kri_valid_demo_cells %>% mutate(age = age/25))


library(broom)
kri_model_df = tidy(kri_lmer_model) %>% filter(group == 'fixed') 
kri_model_df$term = plyr::mapvalues(kri_model_df$term, from = c('layerL2.3', 'layerL5'), to = c('layer2/3', 'layer5'))
kri_model_df$term = factor(kri_model_df$term, levels = c('age', 'sexM', 'layer2/3', 'layer5'))

kri_model_df$var_type = factor(c('layer', 'layer', 'demo', 'demo'), levels = c('layer', 'demo'))
kri_model_df$source = 'Krembil'


kri_model_df %>% ggplot(aes(x = term, y = estimate, color = var_type)) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  geom_pointrange(aes(x = term, y = estimate, 
                      ymin = estimate - std.error, ymax = estimate + std.error)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab('') + ylab('Sag ratio coeff')

model_df_merged = bind_rows(kri_model_df, aibs_model_df)
model_df_merged$term = factor(model_df_merged$term, levels = c('layer2', 'layer2/3', 'layer3', 'layer4', 'layer5', 'layer6', 'diagnosisTumor', 'age', 'sexM'))
model_df_merged$source = factor(model_df_merged$source, levels = c('Krembil', 'Allen'))
model_df_merged$var_type = factor(model_df_merged$var_type, levels = c('layer', 'demo'))


library(MuMIn)

r.squaredGLMM(kri_lmer_model)

model_coeff_comb_plot = model_df_merged %>% ggplot(aes(x = term, y = estimate, color = var_type)) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  geom_pointrange(aes(x = term, y = estimate, 
                      ymin = estimate - std.error, ymax = estimate + std.error)) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  xlab('') + ylab('Sag ratio coeff') + 
  facet_wrap(~source, scales = 'free')


kri_valid_demo_cells$layer = plyr::mapvalues(kri_valid_demo_cells$layer_name, from = c('L2.3', 'L5'), to = c('layer2/3', 'layer5'))
kri_valid_demo_cells$source = 'Krembil (Toronto)'

aibs_human_pyr_ephys$layer = plyr::mapvalues(aibs_human_pyr_ephys$layer_name, from = c('2', '3', '5'), to = c('layer2', 'layer3', 'layer5'))
aibs_human_pyr_ephys$source = 'Allen (Seattle)'

comb_plot_dataset = bind_rows(kri_valid_demo_cells, aibs_human_pyr_ephys %>% filter(layer %in% c('layer2', 'layer3', 'layer5')))
comb_plot_dataset$source = factor(comb_plot_dataset$source, levels = c('Krembil (Toronto)', 'Allen (Seattle)'))
comb_plot_dataset$diagnosis = factor(comb_plot_dataset$diagnosis, levels = c('Epilepsy', 'Tumor', 'Other'))


kri_v_age = comb_plot_dataset %>% filter(source == 'Krembil (Toronto)', layer == 'layer5') %>% 
  ggplot(aes(x = age, y = sag.400)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  # scale_color_manual(values = c('black', 'red', 'purple')) + 
  xlab('Patient age (years)') + ylab('Sag ratio')

aibs_v_age = comb_plot_dataset %>% filter(source == 'Allen (Seattle)') %>% 
  ggplot(aes(x = age, y = sag, group = layer)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  # scale_color_manual(values = c('black', 'red')) +
  facet_grid(~layer) + theme(legend.position="none") +
  xlab('Patient age (years)') + ylab('Sag ratio')

sag_by_age_comb_plot = plot_grid(kri_v_age, aibs_v_age, nrow = 1, rel_widths = c(1.3, 1.2))


kri_lmer_model = lmer(sag.400 ~ (1|subject_id) + age, 
                      data = kri_valid_demo_cells %>% mutate(age = age/20) %>% filter(layer %in% c('L5')) )

kri_lmer_model_no_demo = lmer(sag.400 ~ (1|subject_id)  , 
                      data = kri_valid_demo_cells %>% mutate(age = age/20) %>% filter(layer %in% c('L5')) )

anova(kri_lmer_model, kri_lmer_model_no_demo)

aibs_lmer_model = lmer(sag ~ (1|subject_id) + age, 
                       data = aibs_human_pyr_ephys %>% mutate(age = age/20) %>% filter(layer %in% c('3')) )

aibs_lmer_model_no_demo = lmer(sag ~ (1|subject_id), 
                       data = aibs_human_pyr_ephys %>% mutate(age = age/20) %>% filter(layer %in% c('3')) )

anova(aibs_lmer_model, aibs_lmer_model_no_demo)


model_list = list(kri_lmer_model_no_demo, kri_lmer_model, aibs_lmer_model_no_demo, aibs_lmer_model)

model_rsq = lapply(model_list, function(l) r.squaredGLMM(l)[1][1]) %>% unlist

model_rsq_df = data.frame(models= factor(c('layer only', 'layer + demographics', 'layer only', 'layer + demographics'), 
                                         levels = c('layer only', 'layer + demographics') ),
                          source = factor(c('Krembil', 'Krembil', 'Allen', 'Allen'), levels = c('Krembil', 'Allen')),
                          rsq = model_rsq*100)

model_rsq_df = unite(model_rsq_df, col = 'source_by_model', source, models, remove = F)

model_rsq_df$source_by_model = factor(model_rsq_df$source_by_model, levels = 
                                        c('Krembil_layer_only', 'Krembil_layer+demo',
                                          'Allen_layer_only', 'Allen_layer+demo'))

model_rsq_plot = model_rsq_df %>% ggplot(aes(x = source, y = rsq, fill = models)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  ylab('Model variance exp. (%)') + xlab('') + theme(legend.position="bottom")


top_plot = plot_grid(model_coeff_comb_plot, model_rsq_plot, nrow = 1, rel_widths = c(2, 1))
final_demo_plot = plot_grid(sag_by_age_comb_plot, top_plot, nrow = 2)

ggsave('figures/final_demo_plot.pdf', plot = final_demo_plot, width = 12, height = 7, units = 'in', scale = 1, useDingbats=FALSE)




### supplemental fig for Layer 3 Allen data

aibs_human_pyr_ephys_l3 = aibs_human_ephys %>% filter(dendrite_type == 'spiny', !is.na(depth_norm)) %>%
  rename(layer = structure_layer_name)
aibs_human_pyr_ephys_l3$diagnosis  = plyr::mapvalues(aibs_human_pyr_ephys_l3$diagnosis, c('TLE', 'Tumor'), c('Epilepsy', 'Tumor'))


m1 = lmer(sag ~ depth_norm + age + (1|subject_id), data = aibs_human_pyr_ephys_l3 %>% filter(layer %in% c(2, 3)))
m2 = lmer(sag ~ depth_norm + (1|subject_id), data = aibs_human_pyr_ephys_l3 %>% filter(layer %in% c(2, 3)))

anova(m1, m2)


sag_lm = lm(sag ~ depth_norm, data = aibs_human_pyr_ephys_l3)
aibs_human_pyr_ephys_l3$sag.norm = aibs_human_pyr_ephys_l3$sag - predict(sag_lm, data = aibs_human_pyr_ephys_l3)

aibs_human_pyr_ephys_l3 %>% filter(layer == 3) %>% ggplot(aes(x = age, y = sag)) + geom_point() + stat_smooth(method = "lm")

aibs_human_pyr_ephys_l3 %>% filter(layer == 3) %>% ggplot(aes(x = age, y = sag.norm)) + geom_point() + stat_smooth(method = "lm")
# aibs_human_pyr_ephys = aibs_human_ephys %>% filter(dendrite_type == 'spiny') %>% rename(layer = structure_layer_name)


library(broom)
aibs_model_df = tidy(aibs_lmer_model) %>% filter(group == 'fixed') 

aibs_model_df$term = factor(aibs_model_df$term, levels = c('age', 'sexM', 'layer2',  'layer5'))
aibs_model_df$var_type = factor(c('layer', 'layer', 'demo', 'demo'), levels = c('layer', 'demo'))
aibs_model_df$source = 'Allen'



kri_valid_demo_cells = joined_ephys_data %>% filter(!is.na(age), !is.na(sex), unique_subject) %>% rename(layer = layer_name)
kri_valid_demo_cells$diagnosis  = plyr::mapvalues(kri_valid_demo_cells$diagnosis, c('TLE', 'Tumor', NA, ''), c('Epilepsy', 'Tumor', 'Other', 'Other'))
kri_valid_demo_cells$diagnosis[is.na(kri_valid_demo_cells$diagnosis)] = 'Other'

kri_valid_demo_cells$diagnosis = factor(kri_valid_demo_cells$diagnosis, levels = c('Epilepsy', 'Tumor', 'Other'))

kri_age_sd = kri_valid_demo_cells %>% distinct(subject_id, .keep_all = T) %>% pull(age) %>% sd

kri_lmer_model = lmer(sag.400  ~ layer + (1|subject_id) + age + sex - 1 , data = kri_valid_demo_cells %>% mutate(age = age/25))


library(broom)
kri_model_df = tidy(kri_lmer_model) %>% filter(group == 'fixed') 
kri_model_df$term = plyr::mapvalues(kri_model_df$term, from = c('layerL2.3', 'layerL5'), to = c('layer2/3', 'layer5'))
kri_model_df$term = factor(kri_model_df$term, levels = c('age', 'sexM', 'layer2/3', 'layer5'))

kri_model_df$var_type = factor(c('layer', 'layer', 'demo', 'demo'), levels = c('layer', 'demo'))
kri_model_df$source = 'Krembil'


kri_model_df %>% ggplot(aes(x = term, y = estimate, color = var_type)) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  geom_pointrange(aes(x = term, y = estimate, 
                      ymin = estimate - std.error, ymax = estimate + std.error)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab('') + ylab('Sag ratio coeff')

model_df_merged = bind_rows(kri_model_df, aibs_model_df)
model_df_merged$term = factor(model_df_merged$term, levels = c('layer2', 'layer2/3', 'layer3', 'layer4', 'layer5', 'layer6', 'diagnosisTumor', 'age', 'sexM'))
model_df_merged$source = factor(model_df_merged$source, levels = c('Krembil', 'Allen'))
model_df_merged$var_type = factor(model_df_merged$var_type, levels = c('layer', 'demo'))


library(MuMIn)

r.squaredGLMM(kri_lmer_model)

model_coeff_comb_plot = model_df_merged %>% ggplot(aes(x = term, y = estimate, color = var_type)) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  geom_pointrange(aes(x = term, y = estimate, 
                      ymin = estimate - std.error, ymax = estimate + std.error)) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  xlab('') + ylab('Sag ratio coeff') + 
  facet_wrap(~source, scales = 'free')


kri_valid_demo_cells$layer = plyr::mapvalues(kri_valid_demo_cells$layer, from = c('L2.3', 'L5'), to = c('layer2/3', 'layer5'))
kri_valid_demo_cells$source = 'Krembil (Toronto)'

aibs_human_pyr_ephys$layer = plyr::mapvalues(aibs_human_pyr_ephys$layer, from = c('2', '3', '5'), to = c('layer2', 'layer3', 'layer5'))
aibs_human_pyr_ephys$source = 'Allen (Seattle)'

comb_plot_dataset = bind_rows(kri_valid_demo_cells, aibs_human_pyr_ephys %>% filter(layer %in% c('layer2', 'layer3', 'layer5')))
comb_plot_dataset$source = factor(comb_plot_dataset$source, levels = c('Krembil (Toronto)', 'Allen (Seattle)'))
comb_plot_dataset$diagnosis = factor(comb_plot_dataset$diagnosis, levels = c('Epilepsy', 'Tumor', 'Other'))


kri_v_age = comb_plot_dataset %>% filter(source == 'Krembil (Toronto)') %>% 
  ggplot(aes(x = age, y = sag.400, color = diagnosis, group = layer)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_manual(values = c('black', 'red', 'purple')) + 
  facet_grid(~layer) +  
  xlab('Patient age (years)') + ylab('Sag ratio')


### supplemental figure showing that sag amplitudes are correlated with sag ratio

sag_ratio_vs_amp = comb_plot_dataset %>% filter(source == 'Krembil (Toronto)') %>% 
  ggplot(aes(x = sagamp.400, y = sag.400, group = layer)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  # scale_color_manual(values = c('black', 'red', 'purple')) + 
  facet_grid(~layer, scales = "free_x") +  
  xlab('Sag amplitude (mV)') + ylab('Sag ratio')

ggsave('figures/sag_ratio_vs_amp.pdf', plot = sag_ratio_vs_amp, width = 6, height = 3, units = 'in', scale = 1, useDingbats=FALSE)

# 13d03008; L5; 58 years old
# 15o08038; L5; 21 years old

new_df = comb_plot_dataset 
new_df$colors = 'grey35'
new_df$sizes = 2
new_df[new_df$cell_id == '15o08038.abf', 'colors'] = 'blue'
new_df[new_df$cell_id == '13d03008.abf', 'colors'] = 'red'
new_df[new_df$cell_id %in% c('13d03008.abf','15o08038.abf') , 'sizes'] = 6
new_df[new_df$cohort == 'KRI', 'sag'] = new_df[new_df$cohort == 'KRI', 'sag.400']



sag_ratio_vs_amp = new_df %>% filter(source == 'Krembil (Toronto)', layer == 'layer5') %>% 
  ggplot(aes(x = age, y = sag.400, group = layer, color = colors, size = sizes)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_identity() + 
  scale_size_identity() + 
  xlab('Age (years)') + ylab('Sag ratio')

krembil_lmer_model = lmer(sag.400 ~ (1|subject_id) + age, 
                       data = new_df %>% filter(source == 'Krembil (Toronto)', layer == 'layer2/3'))

krembil_lmer_model_no_demo = lmer(sag.400 ~ (1|subject_id), 
                       data = new_df %>% filter(source == 'Krembil (Toronto)', layer == 'layer5'))

anova(krembil_lmer_model,krembil_lmer_model_no_demo)

library(broom)

kri_fixed_effects_rsq = r.squaredGLMM(krembil_lmer_model)
kri_fixed_effects_corr = sqrt(kri_fixed_effects_rsq[1])

sources = c('Krembil (Toronto)', 'Allen (Seattle)')
layers = c('layer5', 'layer2/3')

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
    
    if ((curr_cohort == "Allen (Seattle)") & (curr_layer  == "layer2/3")){
      use_layers = c('layer2','layer3')
    } else{
      use_layers = curr_layer
    }
    
    lmer_mod = lmer(sag_formula, 
                    data = new_df %>% filter(source == curr_cohort, layer %in% use_layers))
    lmer_mod_no_demo = lmer(sag_formula_no_demo, 
                    data = new_df %>% filter(source == curr_cohort, layer %in% use_layers))
    
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





allen_lmer_model = lmer(sag ~ (1|subject_id) + age, 
                          data = new_df %>% filter(source == 'Allen (Seattle)', layer %in% c('layer2', 'layer3')))

allen_lmer_model_no_demo = lmer(sag ~ (1|subject_id) , 
                                  data = new_df %>% filter(source == 'Allen (Seattle)', layer %in% c('layer5')))

anova(allen_lmer_model,allen_lmer_model_no_demo)

allen_fixed_effects_rsq = r.squaredGLMM(allen_lmer_model)
allen_fixed_effects_corr = sqrt(allen_fixed_effects_rsq[1])

new_df = comb_plot_dataset 
new_df$colors = 'grey35'
new_df$sizes = 2
new_df[new_df$cell_id == '15o08038.abf', 'colors'] = 'blue'
new_df[new_df$cell_id == '13d03008.abf', 'colors'] = 'red'
new_df[new_df$cell_id %in% c('13d03008.abf','15o08038.abf') , 'sizes'] = 6

krembil_xlim = new_df %>% filter(source == 'Krembil (Toronto)', layer == 'layer5') %>% pull(age) %>% range
krembil_xlim[1] = krembil_xlim[1] - 1
krembil_xlim[2] = krembil_xlim[2] + 1

allen_xlim = new_df %>% filter(source == 'Allen (Seattle)') %>% pull(age) %>% range
allen_xlim[1] = allen_xlim[1] - 1
allen_xlim[2] = allen_xlim[2] + 1

sag_ratio_vs_age_l5_krembil = new_df %>% filter(source == 'Krembil (Toronto)', layer == 'layer5') %>% 
  ggplot(aes(x = age, y = sag.400, group = layer, color = colors, size = sizes)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_identity() + 
  scale_size_identity() + 
  xlim(krembil_xlim) + 
  xlab('Age (years)') + ylab('Sag ratio')

sag_ratio_vs_age_l23_krembil = new_df %>% filter(source == 'Krembil (Toronto)', layer == 'layer2/3') %>% 
  ggplot(aes(x = age, y = sag.400, group = layer, color = colors, size = sizes)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_identity() + 
  scale_size_identity() + 
  xlim(krembil_xlim) + 
  xlab('Age (years)') + ylab('Sag ratio')

sag_ratio_vs_age_l5_allen = new_df %>% filter(source == 'Allen (Seattle)', layer == 'layer5') %>% 
  ggplot(aes(x = age, y = sag, group = layer, color = colors, size = sizes)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_identity() +
  scale_size_identity() + 
  # xlim(allen_xlim) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  xlab('Age (years)') + ylab('Sag ratio')

sag_ratio_vs_age_l23_allen = new_df %>% filter(source == 'Allen (Seattle)', layer %in% c('layer2', 'layer3')) %>% 
  ggplot(aes(x = age, y = sag, color = colors, size = sizes)) + 
  geom_smooth(method = 'lm', color = 'grey') + 
  geom_jitter(alpha = .5) + 
  scale_color_identity() + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_size_identity() + 
  # xlim(allen_xlim) +
  xlab('Age (years)') + ylab('Sag ratio')

sag_ratio_all_panels = plot_grid(sag_ratio_vs_age_l5_krembil, sag_ratio_vs_age_l5_allen, sag_ratio_vs_age_l23_krembil, sag_ratio_vs_age_l23_allen, ncol = 2, align = 'v')

ggsave('figures/sag_ratio_all_panels.pdf', plot = sag_ratio_all_panels, width = 5, height = 6, units = 'in', scale = 1, useDingbats=FALSE)

