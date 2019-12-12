neuroelectro_df = read.csv('~/ephys_analysis/data/neuroelectro_csv_filtered.csv')

neuroelectro_df$sagratio[!is.na(neuroelectro_df$sagratio) & (neuroelectro_df$sagratio < .5)] = 1 - neuroelectro_df$sagratio[!is.na(neuroelectro_df$sagratio) & (neuroelectro_df$sagratio < .5)]

neuroelectro_df %>% 
  filter(NeuronName %in% c('Neocortex pyramidal cell layer 2-3', 'Neocortex pyramidal cell layer 5-6', 
                           'Hippocampus CA1 pyramidal cell'), !is.na(sagratio) | !is.na(sagamp)) %>% 
  group_by(Pmid, NeuronName) %>% 
  mutate(Age_count = length(unique(AnimalAge))) %>% 
  mutate(sagratio = 1 - sagratio) %>% 
  # filter(Age_count > 1) %>% 
  select(Age_count, NeuronName, sagratio, sagamp, AnimalAge, Species) %>% 
  ungroup() %>% 
  ggplot(aes(x = AnimalAge, y = sagamp, color = NeuronName)) + geom_point() + 
  facet_grid(~NeuronName*Species)
