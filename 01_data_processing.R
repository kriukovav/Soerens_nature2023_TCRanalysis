
# packages ----------------------------------------------------------------

library(tidyverse)
library(here)

sra_toolkit_directory <- "/home/vkriukova/software/sra-toolkit/sratoolkit.3.0.1-ubuntu64/bin"
slurm_submit_prefetch.sh_path <- here("slurm_submit_prefetch.sh")
slurm_submit_fastqdump.sh_path <- here("slurm_submit_fastqdump.sh")
mixcr_directory <- "/home/vkriukova/software/mixcr/mixcr" #v4.1.0
slurm_submit_mixcr.sh_path <- here("slurm_submit_mixcr.sh")

# prepare for prefetch ----------------------------------------------------

read_csv(here("data", "sra_meta.csv")) %>% 
  select(Run) %>%
  write_tsv("SRRxxx.tsv", col_names = F)

# prefetch ----------------------------------------------------------

paste(slurm_submit_prefetch.sh_path, sra_toolkit_directory) -> slurm_submit_prefetch.sh #create a variable - shell command, which runs the slurm_submit_prefetch.sh script. This script further executes run_prefetch.sh script in sbatches
system(slurm_submit_prefetch.sh) #run slurm_submit_prefetch.sh script

# fasterq-dumb ----------------------------------------------------------

paste(slurm_submit_fastqdump.sh_path, sra_toolkit_directory) -> slurm_submit_fastqdump.sh #create a variable - shell command, which runs the slurm_submit_fastqdump.sh script. This script further executes run_fastqdumb.sh script in sbatches
system(slurm_submit_fastqdump.sh) #run slurm_submit_fastqdump.sh script

# prepare for mixcr -------------------------------------------------------

ngs_directory <- "/projects/cdr3_common/user/vkriukova/soerens_nature2023" #locate the folder with ngs data

filenames <- list.files(ngs_directory, recursive = T, full.names = T) #full paths to fastq files
filenames_short <- str_remove_all(list.files(ngs_directory, recursive = T, full.names = F), pattern = "fastq_|\\/.*") #sample names

data.frame(filenames, "sample_id" = filenames_short) %>%
  mutate(R = case_when(str_detect(.$filenames, pattern = "_1") ~ "R1",
                       str_detect(.$filenames, pattern = "_2") ~ "R2")) %>%
  pivot_wider(names_from = "R", values_from = "filenames") %>%
  relocate(sample_id, .after = R2) %>% 
  write_tsv(here("filenames.txt"), col_names = F) #generate the filenames.txt file containing full paths to the ngs files and sample_ids


# mixcr -------------------------------------------------------------------

mixcr_directory <- "/home/vkriukova/software/mixcrV4.1.2/mixcr" #v4.1.0

paste(slurm_submit_mixcr.sh_path, mixcr_directory) -> slurm_submit_mixcr.sh #create a variable - shell command, which runs the slurm_submit_mixcr.sh script. This script further executes run_mixcr.sh script in slurm batches
system(slurm_submit_mixcr.sh) #run slurm_submit_mixcr.sh script


# mixcr alignment report --------------------------------------------------

filenames <- list.files(here("outs", "mixcr", "stdout_files"), full.names = T) #full paths to mixcr stdout files
filenames_short <- str_remove(list.files(here("outs", "mixcr", "stdout_files"), full.names = F), pattern = ".txt") #sample names

t <- map(filenames, function(x) read_lines(x))
names(t) <- filenames_short

all_data <- imap(t, function(x, y) data.frame("sample_id" = y, "text" = x)) %>%
  bind_rows

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "Total sequencing reads:")) %>%
  mutate(text = str_remove_all(text, "Total sequencing reads: "),
         text = as.double(text)) %>%
  rename(total_sequencing_reads = text) ->
  total_sequencing_reads

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "Successfully aligned reads:")) %>%
  mutate(text = str_remove_all(text, "Successfully aligned reads: | \\(.*"),
         text = as.double(text)) %>%
  rename(total_aligned_reads = text) ->
  total_aligned_reads


all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "TRA chains:")) %>%
  slice_head() %>%
  mutate(text = str_remove_all(text, "TRA chains: ")) %>%
  separate(col = text, into = c("tra_reads_count", "tra_reads_fraction"), " ") %>%
  mutate(tra_reads_fraction = str_remove_all(tra_reads_fraction, "[\\(\\)%]"),
         tra_reads_fraction = as.double(tra_reads_fraction)) ->
  tra_reads

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "TRB chains:")) %>%
  slice_head() %>%
  mutate(text = str_remove_all(text, "TRB chains: ")) %>%
  separate(col = text, into = c("trb_reads_count", "trb_reads_fraction"), " ") %>%
  mutate(trb_reads_fraction = str_remove_all(trb_reads_fraction, "[\\(\\)%]"),
         trb_reads_fraction = as.double(trb_reads_fraction)) ->
  trb_reads


all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "Final clonotype count: ")) %>%
  mutate(text = str_remove_all(text, "Final clonotype count: "),
         text = as.double(text)) %>%
  rename(final_clonotype_count = text) ->
  final_clonotype_count

all_data %>%
  group_by(sample_id) %>%
  filter(str_detect(text, "Average number of reads per clonotype: ")) %>%
  slice_tail(n = 1) %>%
  mutate(text = str_remove_all(text, "Average number of reads per clonotype: "),
         text = as.double(text)) %>%
  rename(average_number_reads_per_clonotype = text) ->
  average_number_reads_per_clonotype

alignment_report <- total_sequencing_reads %>%
  left_join(total_aligned_reads) %>%
  left_join(tra_reads) %>%
  left_join(trb_reads) %>%
  left_join(final_clonotype_count) %>%
  left_join(average_number_reads_per_clonotype) %>%
  left_join(read_csv(here("data", "sra_meta.csv")), by = c("sample_id" = "Run"))

write_tsv(alignment_report, here("outs", "mixcr", "alignment_report.txt"))

# collect final clonesets -------------------------------------------------

filenames <- list.files(here("outs", "mixcr"), pattern = "clones\\_TRB\\.tsv", full.names = T)
filenames_short <- str_remove(list.files(here("outs", "mixcr"), pattern = "clones\\_TRB\\.tsv", full.names = F), pattern = "\\.clones.*")

clonesets <- map(filenames, function(x) read_tsv(x))
names(clonesets) <- filenames_short


clonesets %>%
  map(function(x){
    select(x, readCount, nSeqCDR3, aaSeqCDR3, allVHitsWithScore, allDHitsWithScore, allJHitsWithScore) %>%
      mutate(allVHitsWithScore = str_remove(.$allVHitsWithScore, pattern = "\\*.*")) %>%
      mutate(allJHitsWithScore = str_remove(.$allJHitsWithScore, pattern = "\\*.*")) %>%
      mutate(allDHitsWithScore = str_remove(.$allDHitsWithScore, pattern = "\\*.*")) -> temp_table
    names(temp_table) <- c("count", "cdr3nt", "cdr3aa", "v", "d", "j")
    group_by(temp_table, cdr3nt, cdr3aa, v, d, j) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      filter(str_detect(.$cdr3aa, pattern = "^C.*F$")) %>% #remove CDR3 clonotypes not beginning with conserved C or ending with conserved F
      filter(!str_detect(.$cdr3aa, pattern = "\\_|\\*")) %>%  #remove CDR3 clonotypes with frameshifts or stopcodons
      mutate(freq = count/sum(count)) %>%
      relocate(count, freq, .before = "cdr3nt") %>%
      arrange(desc(count))
  }) -> clonesets_clear

dir.create(path = here("outs", "clonesets"))
iwalk(clonesets_clear, ~ write_tsv(.x, here("outs", "clonesets", paste0(.y, ".txt"))))










