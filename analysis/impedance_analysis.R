# do initial analysis of Fred's resonance analysis

library(readxl)
library(tidyverse)

lihua_l23 = read_excel(path = 'raw-data/impedance_calc_temp/homeira_L5.impedance_info.xlsx', sheet = 2)

path = 'raw-data/impedance_calc_temp/homeira_L5.impedance_info.xlsx'
sheets_list = path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = path)

lapply(sheets_list, function(sheet){
  
}
