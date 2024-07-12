#### PSVs & SNVs ####

#load libraries
library(tidyverse)


args = commandArgs(trailingOnly=TRUE)
#load in all variant calls
SNV_table = read_delim(args[1])
#load in PSV positions
PSVs = read_delim(args[2], delim = "\t")

#load in bed file containing SMN1/2-specific positions
SMN2_specific_positions_bed = read_delim(args[3], delim = "\t", col_names = FALSE)
colnames(SMN2_specific_positions_bed) <- c("chr", "start", "POS", "annotation")

prefix = args[4]

###PROCESSING SNVs###
#add column for genotype_checked
##if cov <3 --> NA
##if cov >=3 & !is.na(ALT_value) --> ALT
##if cov >=3 & ALT_value == NA --> REF
SNV_table <- SNV_table %>%
  mutate(genotype_checked = case_when(depth >= 3 & AF >= 0.5 ~ "ALT",
                                      depth >= 3 & AF <= 0.5 ~ "REF",
                                      depth >= 3 & is.na(AF) ~ "REF",
  )) %>%
  left_join(select(PSVs, POS, PSV), by = "POS") #add PSV column

#make table with SMN1/2 copy type
#filter for PSV13 and add SMN copy type
PSV13_copy_type <- SNV_table %>% filter(PSV == "PSV13") %>%
  mutate(SMN_copy_type = case_when(genotype_checked == "REF" ~ "SMN1",
                                   genotype_checked == "ALT" ~ "SMN2"))

#add SMN_copy_type to SNV table
SNV_table <- left_join(SNV_table, select(PSV13_copy_type, sample_name, hap, SMN_copy_type), by = c("sample_name", "hap")) %>%
  mutate(sample_hap = paste(sample_name, "_", hap, sep = ""))


###SMN1/2 COPY-SPECIFIC VARIANTS###
#Look at SNVs at SMN2-specific positions, but for SMN2_specific_positions_bed
SNV_at_SMN2_specific_positions <- SNV_table %>%
  filter(POS %in% SMN2_specific_positions_bed$POS) %>% #filter for SMN2-specific positions
  mutate(environment = case_when(POS != 71407128 & genotype_checked == "ALT" ~ 2,
                                 POS != 71407128 & genotype_checked == "REF" ~ 1,
                                 POS == 71407128 & genotype_checked == "ALT" ~ 1,
                                 POS == 71407128 & genotype_checked == "REF" ~ 2)) #PSV5 behaves differently because it is a hybrid PSV in reference genome T2T-CHM13

#pivot the table to show all SNVs next to each other at SMN2-specific positions
SNV_at_SMN2_specific_positions_pivoted <- SNV_at_SMN2_specific_positions %>%
  select(sample_name, hap, POS, SMN_copy_type, environment) %>%
  arrange(POS) %>%
  pivot_wider(names_from = POS, values_from = environment) %>%
  arrange(sample_name, hap)
write_tsv(SNV_at_SMN2_specific_positions_pivoted, paste(prefix, "_SNVs_at_SMN2_specific_positions.tsv", sep = ""), quote="none", na="")

