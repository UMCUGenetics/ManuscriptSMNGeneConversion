#### PSVs & SNVs ####

# Load library
library(tidyverse)

###LOADING IN DATA###
# Retrieve command-line arguments passed to the script when it is executed
args = commandArgs(trailingOnly=TRUE)

# Load in all variant calls
SNV_table = read_delim(args[1])

# Load in PSV positions
PSVs = read_delim(args[2], delim = "\t")

# Load in bed file containing SMN1/2-specific positions
SMN_specific_positions_bed = read_delim(args[3], delim = "\t", col_names = FALSE)
colnames(SMN_specific_positions_bed) <- c("chr", "start", "POS", "annotation")

# Define prefix, either "SMA" or "1000G"
prefix = args[4]

###PROCESSING SNVs###
# Add column for genotype_checked
# If read depth is at least 3 and allele frequency at least 0.5 --> genotype_checked = "ALT"
# If read depth is at least 3 and allele frequency lower than 0.5 --> genotype_checked = "REF"
# If read depth is at least 3 and no variant call has been made --> genotype_checked = "REF"
SNV_table <- SNV_table %>%
  mutate(genotype_checked = case_when(depth >= 3 & AF >= 0.5 ~ "ALT",
                                      depth >= 3 & AF < 0.5 ~ "REF",
                                      depth >= 3 & is.na(AF) ~ "REF",
  )) %>%
  left_join(select(PSVs, POS, PSV), by = "POS") #add PSV column

# Make table with SMN1/2 copy type
# Filter for PSV13 and add SMN copy type
PSV13_copy_type <- SNV_table %>% filter(PSV == "PSV13") %>%
  mutate(SMN_copy_type = case_when(genotype_checked == "REF" ~ "SMN1",
                                   genotype_checked == "ALT" ~ "SMN2"))

# Add SMN_copy_type to SNV table
SNV_table <- left_join(SNV_table, select(PSV13_copy_type, sample_name, hap, SMN_copy_type), by = c("sample_name", "hap")) %>%
  mutate(sample_hap = paste(sample_name, "_", hap, sep = ""))


###SMN1/2 COPY-SPECIFIC VARIANTS###
# Look at SNVs at SMN1/2-specific positions, for SMN_specific_positions_bed
SNV_at_SMN_specific_positions <- SNV_table %>%
  left_join(SMN_specific_positions_bed) %>%
  filter(POS %in% SMN_specific_positions_bed$POS) %>% #filter for PSVs and SMN1/2-specific positions
  mutate(position = case_when(!is.na(PSV) ~ PSV,
                              is.na(PSV) ~ as.character(POS))) %>%
  mutate(environment = case_when(annotation == "SMN1-spec" & genotype_checked == "REF" ~ 2,
                                 annotation == "SMN1-spec" & genotype_checked == "ALT" ~ 1,
                                 annotation == "SMN2-spec" & genotype_checked == "REF" ~ 1,
                                 annotation == "SMN2-spec" & genotype_checked == "ALT" ~ 2,
                                 !is.na(PSV) & PSV != "PSV5" & genotype_checked == "REF" ~ 1,
                                 !is.na(PSV) & PSV != "PSV5" & genotype_checked == "ALT" ~ 2,
                                 !is.na(PSV) & PSV == "PSV5" & genotype_checked == "REF" ~ 2,
                                 !is.na(PSV) & PSV == "PSV5" & genotype_checked == "ALT" ~ 1)) #PSV5 behaves differently because it is a hybrid PSV in reference genome T2T-CHM13

# Pivot the table to show all SNVs next to each other at SMN1/2-specific positions
SNV_at_SMN_specific_positions_pivoted <- SNV_at_SMN_specific_positions %>%
  arrange(POS) %>%
  select(sample_name, hap, position, SMN_copy_type, environment) %>%
  pivot_wider(names_from = position, values_from = environment) %>%
  arrange(sample_name, hap)
SNV_at_SMN_specific_positions_pivoted_incl_empty_haplotypes <- SNV_table %>%
  select(sample_name, hap) %>%
  distinct() %>%
  arrange(sample_name, hap) %>% 
  left_join(SNV_at_SMN_specific_positions_pivoted)
write_tsv(SNV_at_SMN_specific_positions_pivoted_incl_empty_haplotypes, paste(prefix, "_SNVs_at_SMN1-2_specific_positions.tsv", sep = ""), quote="none", na="")
