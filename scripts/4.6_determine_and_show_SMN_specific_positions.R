#### PSVs & SNVs ####

# Load library
library(tidyverse)

###LOADING IN DATA###
# Retrieve command-line arguments passed to the script when it is executed
args = commandArgs(trailingOnly=TRUE)

# Load in all variant calls
SNV_table = read_delim(args[1], delim = "\t")

# Load in PSV positions
PSVs = read_delim(args[2], delim = "\t")

# Define prefix, either "SMA" or "1000G"
prefix = args[3]


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

###DETERMINE SMN1/2 COPY-SPECIFIC VARIANTS###

# Make table with variants divided by SMN1 and SMN2
SNV_tally_SMN <- SNV_table %>%
  select(POS, sample_hap, genotype_checked, SMN_copy_type) %>%
  group_by(POS,SMN_copy_type,genotype_checked) %>%
  tally() %>%
  mutate(SMN_copy_genotype = paste(as.character(SMN_copy_type), "_",
                                   as.character(genotype_checked), sep="")) %>%
  ungroup() %>%
  group_by(POS, SMN_copy_genotype) %>%
  summarise(N = sum(n)) %>%
  pivot_wider(names_from = SMN_copy_genotype, values_from = N) %>%
  select(POS, SMN1_REF, SMN1_ALT, SMN1_NA, SMN2_REF, SMN2_ALT, SMN2_NA) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(SMN1_tot_haps_non_NA = SMN1_REF+SMN1_ALT) %>%
  mutate(SMN2_tot_haps_non_NA = SMN2_REF+SMN2_ALT) %>%
  mutate(ratio_SMN1_haps = SMN1_ALT / SMN1_tot_haps_non_NA) %>%
  mutate(ratio_SMN2_haps = SMN2_ALT / SMN2_tot_haps_non_NA) %>%
  left_join(PSVs, by = "POS") %>%
  mutate(ratio_SMN1_haps = case_when(PSV != "PSV5" | is.na(PSV) ~ ratio_SMN1_haps,
                                     PSV == "PSV5" ~ 1 - ratio_SMN1_haps),
         ratio_SMN2_haps = case_when(PSV != "PSV5" | is.na(PSV) ~ ratio_SMN2_haps,
                                     PSV == "PSV5" ~ 1 - ratio_SMN2_haps)) %>% #PSV5 behaves differently because it is a hybrid PSV in reference genome T2T-CHM13
  mutate(copy_specific = case_when(ratio_SMN1_haps >=0.9 & ratio_SMN2_haps <=0.1 ~ "SMN1-spec",
                                   ratio_SMN2_haps >=0.9 & ratio_SMN1_haps <=0.1 ~ "SMN2-spec"))

# Output variants that are SMN1- or SMN2-specific with ratios
SNV_tally_SMN_filtered <- SNV_tally_SMN %>%
  filter(!is.na(copy_specific) | !is.na(PSV)) %>%
  select(-start, -col1, -chrom) %>%
  filter(SMN1_tot_haps_non_NA >= 20 & SMN2_tot_haps_non_NA >= 20)
write_tsv(SNV_tally_SMN_filtered, paste(prefix, "_SNV_tally_at_PSVs_and_SMN1-2_specific_positions.tsv", sep = ""), quote="none", na="")

# Make a bed file with the PSVs and copy-specific positions
SMN_specific_positions_bed <- SNV_tally_SMN_filtered %>%
  ungroup() %>%
  mutate(annotation = case_when(!is.na(PSV) ~ PSV,
                                is.na(PSV) ~ copy_specific)) %>%
  mutate(chr = "chr5") %>%
  mutate(start = POS - 1) %>%
  mutate(stop = POS) %>%
  select(chr, start, stop, annotation)

# Write bed files with PSVs and SMN1/2-specific positions
write.table(SMN_specific_positions_bed, paste(prefix, "_PSVs_and_SMN1-2_specific_positions.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep="\t", quote=FALSE)

# Look at SNVs at SMN1/2-specific positions, for SMN_specific_positions_bed
colnames(SMN_specific_positions_bed) <- c("chr", "start", "POS", "annotation")
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

