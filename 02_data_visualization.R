# packages ----------------------------------------------------------------

library(tidyverse)
library(here)


# import clonesets --------------------------------------------------------

filenames <- list.files(here("outs", "clonesets"), full.names = T)
filenames_short <- list.files(here("outs", "clonesets"), full.names = F) %>% str_remove(pattern = ".txt")

clonesets <- map(filenames, read_tsv)
names(clonesets) <- filenames_short

clonesets2 <- bind_rows(clonesets, .id = "sample_id")


# add metadata ------------------------------------------------------------

clonesets_w_meta <- clonesets2 %>%
  left_join(read_csv(here("data", "sra_meta.csv")), by = c("sample_id" = "Run"))

clonesets_w_meta %>%
  mutate(TREATMENT = case_when(str_detect(.$TREATMENT, pattern = "^[0-9] ") ~ paste0(0, .$TREATMENT),
                               TRUE ~ .$TREATMENT)) %>% #just data cleaning to make the time points look good and in the correct order on plots 
  mutate(sample_treatment = paste(sample_id, TREATMENT, sep = "_")) %>% #will be colnames for wide format data
  select(sample_treatment, cdr3aa, cdr3nt, v, j, freq) %>%
  pivot_wider(names_from = sample_treatment, values_from = freq, values_fill = 10^-6, id_cols = c("cdr3aa", "cdr3nt", "v", "j")) %>% #transform data to wide format (to make all clonotypes present in all samples in pseudofrequency = 10^-6)
  pivot_longer(names_to = "sample_treatment", values_to = "freq", cols = starts_with("SRR")) %>% # return back to long format
  separate(sample_treatment, into = c("sample_id", "treatment"), sep = "_") %>%
  mutate(line_group = paste(cdr3aa, cdr3nt, v, j, sep = "_")) %>%
  group_by(cdr3aa, cdr3nt, v, j, treatment) %>%
  mutate(mean_freq = mean(freq)) %>%
  ungroup() %>%
  mutate(treatment_sample_id = paste(treatment, sample_id, sep = "_")) -> plot_data

plot_data %>%
  distinct(cdr3aa, mean_freq) %>%
  filter(mean_freq > 0.5) %>%
  select(cdr3aa) %>% pull %>% unique -> expanded_clones


  
ggplot(data = plot_data, aes(x = treatment)) +
  geom_line(aes(group = line_group, y = mean_freq), alpha = 0.5, size = 0.2) +
  geom_line(data = plot_data %>% filter(cdr3aa %in% expanded_clones), aes(group = line_group, y = mean_freq, color = cdr3aa), alpha = 0.8, size = 0.6) +
  geom_point(aes(y = freq), size = 0.2, alpha = 0.5) +
  theme_minimal() + 
  theme(legend.position = "none", axis.text.x=element_text(angle = 90),
        panel.grid.major.y = element_line(size = 0.05, linetype = 'solid', colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("number of immunizations") +
  ylab("clone frequency") +
  ggtitle("tracking of T cell clones through generations of prime-boost-boost challenges")


ggsave(plot = last_plot(), filename = "fig1.pdf", path = here("outs"), height = 6, width = 8)





