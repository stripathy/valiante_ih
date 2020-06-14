library(ggbeeswarm)
kri_fi_curves = read.csv('~/Documents/GitHub/valiante_lab_abf_process/output_files/fi_curves.csv')

# read cell meta
cell_meta = read.csv('summary_tables/cell_metadata.csv')

# restrict plot to just pyramidal cells
#cell_meta = cell_meta %>% filter(cell_type == 'Pyr')


cell_meta$layer_name = plyr::mapvalues(cell_meta$layer_name, from = levels(cell_meta$layer_name), to = c("Layer 2/3", "Layer 3C", "Layer 5"))

# read cell ephys data
cell_patient_ephys_combined = read.csv('summary_tables/cell_patient_ephys_combined.csv')
cell_patient_ephys_combined$layer_name = plyr::mapvalues(cell_patient_ephys_combined$layer_name, from = levels(cell_patient_ephys_combined$layer_name), to = c("Layer 2/3", "Layer 3C", "Layer 5"))


kri_fi_curves = merge(kri_fi_curves, cell_meta %>% distinct(cell_id, .keep_all = T), by = 'cell_id')



#

kri_ephys_long = cell_patient_ephys_combined %>% select(kri_ephys_features, cell_id, 
                                                              recorder_name, layer_name, cell_type) %>% 
  unite(col = 'layer_cell_type', cell_type, layer_name, remove = F) %>%
gather(key = 'ephys_name', value = 'ephys_val', -cell_id, -recorder_name, -layer_cell_type, -layer_name, -cell_type)

rheo_fig = kri_ephys_long %>% filter(ephys_name == 'adapt') %>%
ggplot(aes(x = layer_cell_type, y = ephys_val, color = layer_name)) +
scale_color_manual(values = c('blue', 'turquoise4', 'red')) +
geom_boxplot(outlier.alpha = 0)  + geom_jitter(alpha = .5, width = .25) + ylab('AP amplitude (mV)') + xlab('')


kri_fi_curves %>% filter(avg_rate > 0, cell_type == 'Pyr')

kri_fi_curves %>% filter(avg_rate > 0) %>% 
  ggplot(aes(x = norm_stim_amp, y = avg_rate, group = cell_id, color = cell_type)) + 
  geom_line(alpha = .3)  +
  xlim(c(0, 5))

kri_fi_curves %>% filter(avg_rate > 0, cell_type == 'Pyr', layer_name == 'L2.3') %>% 
  ggplot(aes(x = norm_stim_amp, y = avg_rate, group = cell_id, color = layer_name)) + geom_line(alpha = .7)  +xlim(c(0, 3)) + geom_smooth()

kri_fi_curves %>% filter(avg_rate > 0, cell_type == 'Pyr') %>% 
  ggplot(aes(x = stim_amp, y = avg_rate, group = cell_id, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_line(alpha = .1)  + xlim(0, 400) + facet_wrap(~layer_name) + ylim(0, 50)

abf_files_with_fi_curves = unique(kri_fi_curves$cell_id)

# max_rate = 0
# tol_ratio = .2
# inverting_cells = c()
# for (ind in 1:length(abf_files_with_fi_curves)){
#   curr_cell_id = abf_files_with_fi_curves[ind] %>% as.character()
#   max_freq = kri_fi_curves %>% filter(cell_id == curr_cell_id) %>% select(avg_rate) %>% max()
#   cell_fi_curves = kri_fi_curves %>% filter(cell_id == curr_cell_id)
#   last_freq = cell_fi_curves[nrow(cell_fi_curves), 'avg_rate']
#   if (max_freq > last_freq * (1 + tol_ratio)){
#     inverting_cells = c(inverting_cells, curr_cell_id)
#   }
# }
# 
# inverting_cell_df = cell_meta %>% filter(cell_id %in% inverting_cells)
# write.csv(inverting_cell_df, 'summary_tables/inverting_cells.csv')


kri_ap_curves = read.csv('~/Documents/GitHub/valiante_lab_abf_process/output_files/ap_features.csv')
kri_ap_curves = merge(kri_ap_curves %>% arrange(cell_id, X), cell_meta %>% distinct(cell_id, .keep_all = T), by = 'cell_id') 
    
kri_ap_curves %>% filter(cell_type == 'Pyr') %>% 
  ggplot(aes(x = t*1000, y = v, group = cell_id, color = layer_name)) + 
  geom_line(alpha = .1) + 
  facet_wrap(~layer_name) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  xlab('time (ms)') + ylab('voltage (mV)') + xlim(c(0, 10))

kri_ap_curves  %>% filter(cell_type == 'Pyr') %>% 
  ggplot(aes(x = v, y = dvdt, color = layer_name)) + 
  geom_path(alpha = .1) + 
  facet_wrap(~layer_name) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  xlab('voltage (mV)') + ylab('dV dt (mV)')

ap_df = kri_ap_curves %>% filter(cell_id == '13d03005.abf')

library(zoo)

resample_ap_timeseries = function(ap_df){
  zs <- zoo(kri_ap_curves %>% filter(cell_id == '2020_03_02_0023.abf') %>% pull(v), kri_ap_curves %>% filter(cell_id == '2020_03_02_0023.abf') %>% pull(t))    # high freq
  #print(head(ap_df))
  sampling_rate = 1 / (ap_df[2, 't'] - ap_df[1, 't'])
  print(sampling_rate)
  if (sampling_rate <= 25000){
    
    zv <- zoo(ap_df$v, ap_df$t)    # low freq
    #print(zv)
    z <- merge(zs,zv)
    #print(z)
    zdvdt <- zoo(ap_df %>% pull(dvdt), ap_df %>% pull(t))    # low freq
    z <- merge(z, zdvdt)
    
    #print(z)
    
    # Interpolate calibration data (na.spline could also be used)
    z$zv <- na.spline(z$zv)
    z$zdvdt <- na.spline(z$zdvdt)
    z = z[!is.na(z$zs), ]
    ret_df = data.frame('t' = index(z), 'v' = z$zv, 'dvdt' = z$zdvdt)
    return(ret_df)
  }else{
    return(NULL)
    return(ap_df %>% as.data.frame %>% select(t, v, dvdt))
  }
}


unique_cell_ids = kri_fi_curves %$% cell_id %>% unique() %>% as.character()

ret_list = lapply(unique_cell_ids, function(curr_cell_id){
  
  ap_df = kri_ap_curves %>% filter(cell_id == curr_cell_id)
  print(curr_cell_id)
  #print(ap_df)
  new_ap_df = resample_ap_timeseries(ap_df)
  new_ap_df$cell_id = curr_cell_id
  new_ap_df = data.frame(new_ap_df)
})

ret_list_df = ret_list %>% bind_rows()

mean_ap_wave_features = left_join(ret_list_df, cell_meta) 

mean_ap_wave_features = mean_ap_wave_features %>% 
  group_by(layer_name, cell_type, t) %>% 
  summarize(sd_v = sd(v, na.rm = T),
            sem_v = sd_v / sqrt(n() - 1),
            v = mean(v, na.rm = T), 
            sd_dvdt = sd(dvdt, na.rm = T),
            sem_dvdt = sd_dvdt / sqrt(n() - 1),
            dvdt = mean(dvdt), 
            ) %>% arrange(layer_name, cell_type, t)

kri_ap_curves = bind_rows(kri_ap_curves %>% mutate(type = 'single_cell'), mean_ap_wave_features %>% mutate(type = 'group_average'))

ap_waveform_plot = kri_ap_curves %>% filter(cell_type == 'Pyr') %>% 
  ggplot(aes(x = t*1000, y = v, group = cell_id, color = layer_name, alpha = type)) + 
  geom_line(aes(alpha = type)) + scale_alpha_discrete(range = c(1, .1)) + 
  facet_wrap(~layer_name) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  xlab('time (ms)') + ylab('voltage (mV)') + xlim(c(0, 10)) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ap_phase_plot = kri_ap_curves %>% filter(cell_type == 'Pyr') %>% 
  ggplot(aes(x = v, y = dvdt, color = layer_name, alpha = type)) + 
  geom_path(aes(alpha = type)) + scale_alpha_discrete(range = c(1, .1)) + 
  facet_wrap(~layer_name) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  xlab('voltage (mV)') + ylab('dV/dt (mV/ms)') + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

rheo_fig = kri_ephys_long  %>% filter(ephys_name == 'adapt') %>% 
  ggplot(aes(x = layer_cell_type, y = ephys_val, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_boxplot(outlier.alpha = 0)  + geom_jitter(alpha = .5, width = .25) + ylab('AP amplitude (mV)') + xlab('')

a = kri_ephys_long %>% filter(ephys_name %in% kri_ephys_features) %>% filter(ephys_name == 'rin', layer_name %in% c('L5', 'L2.3', 'L3c')) %>% select(ephys_val, layer_name)
t.test(a %>% filter(layer_name == 'L2.3') %$% ephys_val, a %>% filter(layer_name == 'L5') %$% ephys_val)

figure <- ggarrange(ap_waveform_plot, ap_phase_plot, 
                    labels = c("A", "B"), common.legend = TRUE, legend = 'none',
                    ncol = 1, nrow = 2)

mean_ap_wave_features %>% arrange(layer_name, t) %>%ggplot(aes(x = t*1000, y = v, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_line() + xlim(c(0, 10))

mean_ap_wave_features %>% arrange(layer_name, t) %>%ggplot(aes(x = v, y = dvdt, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_path() 


k = kri_fi_curves %>% filter(avg_rate > 0, cell_type == 'Pyr')  %>% filter(cell_id == '13d03005.abf') %>% arrange(norm_stim_amp)

resample_fi_curves = function(fi_df){
  target_fi_norm_amps = seq(1, 5, by = .1)
  resampled_fi_curves = approx(fi_df$norm_stim_amp, y = fi_df$avg_rate, xout = target_fi_norm_amps) %$% y
  new_df = data.frame('norm_stim_amp' = target_fi_norm_amps, 'avg_rate' = resampled_fi_curves)
  return(new_df)
}

resample_unnorm_fi_curves = function(fi_df){
  target_fi_norm_amps = seq(5, 2000, by = 5)
  resampled_fi_curves = approx(fi_df$stim_amp, y = fi_df$avg_rate, xout = target_fi_norm_amps) %$% y
  new_df = data.frame('stim_amp' = target_fi_norm_amps, 'avg_rate' = resampled_fi_curves)
  return(new_df)
}
# attempt fi curve interpolation

fi_resampled_df = lapply(unique_cell_ids, function(curr_cell_id){
  
  ap_df = kri_fi_curves %>% filter(cell_id == curr_cell_id)
  new_ap_df = resample_unnorm_fi_curves(ap_df)
  new_ap_df$cell_id = curr_cell_id
  return(new_ap_df)
}) %>% bind_rows()

resampled_fi_features = left_join(fi_resampled_df, cell_meta) 

mean_fi_features = resampled_fi_features %>% 
  group_by(cell_type, layer_name, stim_amp) %>% 
  summarize(sd_rate = sd(avg_rate, na.rm = T), 
            sem_rate = sd_rate / sqrt(n() -1 ),
            avg_rate = mean(avg_rate, na.rm = T))

group_avg_fi_features = resampled_fi_features %>% filter(stim_amp == 300) %>%
  group_by(cell_type, layer_name) %>% 
  summarize(mean_rate = mean(avg_rate, na.rm = T), 
            sd_rate = sd(avg_rate, na.rm = T),
            n_rate = n(),
            sem_rate = sd_rate / sqrt(n_rate - 1))

kri_fi_curves %>% filter(avg_rate > 0) %>% 
  ggplot(aes(x = stim_amp, y = avg_rate, group = cell_id, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_line(alpha = .1)  + xlim(0, 400) + facet_wrap(~layer_name*cell_type) + ylim(0, 50)


library("ggpubr")
theme_set(
  theme_cowplot(14) +
    theme(legend.position = "bottom")
)

ap_waveform_fig = mean_ap_wave_features %>% filter(cell_type == 'Pyr') %>% arrange(layer_name, t) %>%
  ggplot(aes(x = t*1000, y = v, color = layer_name)) + 
  #geom_ribbon(aes(ymin = v - sem_v, ymax = v + sem_v, fill = layer_name), alpha = .25, color = NA) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  #scale_fill_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_line() + xlim(c(0, 7.5)) + ylim(c(-50, 40)) + ylab('Voltage (mV)') + xlab('Time (ms)')

ap_dvdt_fig = mean_ap_wave_features %>% filter(cell_type == 'Pyr') %>% arrange(layer_name, t) %>% 
  ggplot(aes(x = v, y = dvdt, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_path() + xlab('Voltage (mV)') + ylab('Voltage deriv. (mV/ms)')

rheo_fig = kri_ephys_long %>% filter(ephys_name %in% kri_ephys_features) %>% filter(ephys_name == 'rheo') %>% 
  ggplot(aes(x = layer_name, y = ephys_val, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_boxplot(outlier.alpha = 0)  + geom_quasirandom(alpha = .5, width = .25) + ylab('Rheobase (pA)') + xlab('')

aphw_fig = kri_ephys_long %>% filter(ephys_name %in% kri_ephys_features) %>% filter(ephys_name == 'aphw') %>% 
  ggplot(aes(x = layer_name, y = ephys_val, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_boxplot(outlier.alpha = 0)  + geom_quasirandom(alpha = .5, width = .25) + ylab('AP half-width (mV)') + xlab('')

adratio_fig = kri_ephys_long %>% filter(ephys_name %in% kri_ephys_features) %>% filter(ephys_name == 'adratio') %>% 
  ggplot(aes(x = layer_name, y = ephys_val, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_boxplot(outlier.alpha = 0)  + geom_quasirandom(alpha = .5, width = .25) + ylab('Adaptation index') + xlab('')

fi_slope = kri_ephys_long %>% filter(ephys_name %in% kri_ephys_features) %>% filter(ephys_name == 'rheo') %>% 
  ggplot(aes(x = layer_name, y = ephys_val, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_boxplot(outlier.alpha = 0)  + geom_quasirandom(alpha = .5, width = .25) + ylab('Rheobase (pA)') + xlab('')

fi_fig = mean_fi_features %>% filter(cell_type == 'Pyr') %>% ggplot(aes(x = stim_amp, y = avg_rate , color = layer_name)) + 
  geom_ribbon(aes(ymin = avg_rate - sem_rate, ymax = avg_rate + sem_rate, fill = layer_name), alpha = .25, color = NA) + 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  scale_fill_manual(values = c('blue', 'turquoise4', 'red')) + 
  geom_line(aes(y = avg_rate)) + xlim(c(0, 300)) + ylim(c(0, 25)) + ylab('Firing rate (Hz)') + xlab('Current (pA)')

kri_fi_curves = bind_rows(kri_fi_curves %>% mutate(type = 'single_cell'), mean_fi_features %>% mutate(type = 'group_average'))

fi_plot = kri_fi_curves  %>% filter(cell_type == 'Pyr' | type == 'group_average', avg_rate > 0) %>% 
  ggplot(aes(x = stim_amp, y = avg_rate, color = layer_name, alpha = type, group = cell_id)) + 
  geom_line(aes(alpha = type)) + scale_alpha_discrete(range = c(1, .1)) + 
  facet_wrap(~layer_name) +  xlim(c(0, 300)) + ylim(0, 50)+ 
  scale_color_manual(values = c('blue', 'turquoise4', 'red')) + 
  ylab('Firing rate (Hz)') + xlab('Injected current (pA)') + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# figure <- ggarrange(ap_waveform_plot, ap_phase_plot, fi_plot, 
#                     labels = c("A", "B", "C"), common.legend = TRUE, legend = 'none',
#                     ncol = 1, nrow = 3, align = 'hv')

figure <- ggarrange(ap_waveform_fig, ap_dvdt_fig, aphw_fig, fi_fig, rheo_fig, adratio_fig,
                    labels = c("A", "B", "C", "D", "E", "F"), common.legend = TRUE, legend = 'none',
                    ncol = 3, nrow = 2, align = 'hv')
figure

ggsave(file = 'figures/Fig3_spike_waveform_train.pdf', plot = figure, width = 12, height = 8, device = "pdf")




ap_waveform_fig = mean_ap_wave_features %>% filter(cell_type == 'Int') %>% arrange(layer_name, t) %>%
  ggplot(aes(x = t*1000, y = v, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'red', 'turquoise4')) + 
  geom_line() + xlim(c(0, 2)) + ylim(c(-50, 15)) + ylab('Voltage (mV)') + xlab('Time (ms)')

ap_dvdt_fig = mean_ap_wave_features %>% filter(cell_type == 'Int') %>% arrange(layer_name, t) %>% 
  ggplot(aes(x = v, y = dvdt, color = layer_name)) + 
  scale_color_manual(values = c('blue', 'red', 'turquoise4')) + 
  geom_path() + xlab('Voltage (mV)') + ylab('Voltage deriv. (mV/ms)')
  
  figure <- ggarrange(ap_waveform_fig, ap_dvdt_fig, common.legend = TRUE, legend = 'none',
                      ncol = 2, nrow = 1, align = 'hv')
  ggsave(file = 'figures/interneuron_spike_waveforms.pdf', plot = figure, width = 5, height = 2.5, device = "pdf")
  
  geom_path() + xlab('Voltage (mV)') + ylab('Voltage deriv. (mV/ms)')