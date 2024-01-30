# Testing filtering/decontamination of samples from extraction blanks in specific extraction plate ## 
### Sept. 23, 2020
## Code by Mareli Sánchez Juliá and Bolívar Aponte Rolón

library(readxl)
library(tidyverse)

data <- read_excel( "Pool1R1_95zotu_testFilter.xlsx")
View(data)

filter_data <- read_excel("ExtractionBlankTables_ELT.xlsx")
View(filter_data)


filter_data <-filter_data %>% 
  pivot_longer(cols = 1:12, names_to = "samples", values_to = "values") %>%
  select(-samples)

new_data <- data %>% t() %>%  as.data.frame()
names(new_data) <- new_data %>% slice(1) %>% unlist() #%>%
new_data <- rownames_to_column(new_data, "sample_names")
new_data <- new_data %>% slice(-1)

new_data_filtered <- semi_join(new_data, filter_data, by= c("sample_names" = "values"))
filter_filter <- anti_join( filter_data, new_data_filtered, by = c("values"="sample_names"))


#new_data <- new_data %>% slice(2:323)
View(new_data_filtered)

