#### PSVs & SNVs ####
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
#load in all variant calls
SNV_table = read_delim(args[1])
#load in PSV positions
PSVs = read_delim(args[2], delim = "\t")

prefix = args[3]

###LOADING IN DATA###
#load in all variant calls
#SNV_table <- read_delim("vcf_depth_merged_all_haps.tsv", 
#load in PSV positions
#PSVs <- read_delim("PSV_liftover_hg19_to_T2T_CHM13.txt", delim = "\t")

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

#add copy_type to SNV table
SNV_table <- left_join(SNV_table, select(PSV13_copy_type, sample_name, hap, SMN_copy_type), by = c("sample_name", "hap"))

#add column with position_ref_alt like in Chen 2023 supplementary table
SNV_table <- SNV_table %>%
  mutate(sample_hap = paste(sample_name, "_", hap, sep = "")) %>%
  mutate(pos_ref_alt = case_when(genotype_checked == "ALT" ~ paste(POS, "_", REF, "_", ALT_value, sep = ""))) %>%
  mutate(pos_ref_alt_present = case_when(!is.na(pos_ref_alt) ~ "x"))
#check how many unique variants there are; there are 2497
length(unique(SNV_table$pos_ref_alt))

#make a big pivoted table with all variants REF/ALT
SNVs_pivoted <- SNV_table %>%
  select(POS, sample_name, hap, genotype_checked) %>%
  pivot_wider(names_from = POS, values_from = genotype_checked)

#make table like in chen 2023
SNVs_pivoted_paraphase_suppl <- SNV_table %>%
  select(sample_hap,pos_ref_alt, pos_ref_alt_present) %>%
  pivot_wider(names_from = sample_hap, values_from = pos_ref_alt_present, values_fill = NA, values_fn = length) %>%
  filter(!is.na(pos_ref_alt))
#write_tsv(SNVs_pivoted_paraphase_suppl, "SNVs_pivoted_paraphase_suppl_made_in_R.tsv")
write_tsv(SNVs_pivoted_paraphase_suppl, paste(prefix, "_SNVs_pivoted_paraphase_suppl_made_in_R.tsv", sep = ""), quote = "none")


#write a tsv file with sample, hap and copy type
list_haplotypes_copy_type <- SNV_table %>%
  filter(PSV == "PSV13") %>%
  select(sample_hap, SMN_copy_type)
#write_tsv(list_haplotypes_copy_type, "list_haplotypes_copy_type.tsv")
write_tsv(list_haplotypes_copy_type, paste(prefix, "_list_haplotypes_copy_type.tsv", sep = ""), quote = "none")

