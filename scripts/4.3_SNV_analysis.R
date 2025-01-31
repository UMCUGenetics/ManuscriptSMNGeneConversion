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
SNV_table <- left_join(SNV_table, select(PSV13_copy_type, sample_name, hap, SMN_copy_type), by = c("sample_name", "hap"))

# Add column with position_ref_alt
SNV_table <- SNV_table %>%
  mutate(sample_hap = paste(sample_name, "_", hap, sep = "")) %>%
  mutate(pos_ref_alt = case_when(genotype_checked == "ALT" ~ paste(POS, "_", REF, "_", ALT_value, sep = ""))) %>%
  mutate(pos_ref_alt_present = case_when(!is.na(pos_ref_alt) ~ "x"))

# Report the total number of variants
print("The total number of unique variants is:")
length(unique(filter(SNV_table, !is.na(pos_ref_alt))$pos_ref_alt))
print("The total number of unique SNVs is:")
length(unique(filter(filter(SNV_table, !is.na(pos_ref_alt)), variant_type == "SNV")$pos_ref_alt))
print("The total number of unique INDELs is:")
length(unique(filter(filter(SNV_table, !is.na(pos_ref_alt)), variant_type == "INS" | variant_type == "DEL")$pos_ref_alt))

# Count number of unique variants per sample
unique_variants_per_sample_count <- filter(SNV_table, !is.na(pos_ref_alt)) %>%
  group_by(sample_name) %>%
  summarise(unique_variants_per_sample = n_distinct(pos_ref_alt))

# Report the number of variants per sample
print("The median number of unique variants per sample is:")
median(unique_variants_per_sample_count$unique_variants_per_sample)
print("The minimum number of unique variants per sample is:")
min(unique_variants_per_sample_count$unique_variants_per_sample)
print("The maximum number of unique variants per sample is:")
max(unique_variants_per_sample_count$unique_variants_per_sample)

# Make a pivoted table with all variants REF/ALT
SNVs_pivoted <- SNV_table %>%
  select(POS, sample_name, hap, genotype_checked) %>%
  pivot_wider(names_from = POS, values_from = genotype_checked)

# Make table to report variants in supplementary file
# Include empty haplotypes
SNVs_all_pivoted_supplementary <- SNV_table %>%
  select(sample_hap,pos_ref_alt, pos_ref_alt_present) %>%
  pivot_wider(names_from = sample_hap, values_from = pos_ref_alt_present, values_fill = NA, values_fn = length) %>%
  filter(!is.na(pos_ref_alt))
write_tsv(SNVs_all_pivoted_supplementary, paste(prefix, "_variants_pivoted_supplementary.tsv", sep = ""), quote = "none", na="")

# Write a tsv file with sample, hap and copy type
# Include empty haplotypes
list_haplotypes_copy_type <- SNV_table %>%
  filter(PSV == "PSV13") %>%
  select(sample_hap, SMN_copy_type)
list_haplotypes_copy_type_incl_empty_haplotypes <- SNV_table %>%
  select(sample_hap) %>%
  distinct() %>%
  arrange(sample_hap) %>% 
  left_join(list_haplotypes_copy_type)
write_tsv(list_haplotypes_copy_type_incl_empty_haplotypes, paste(prefix, "_list_haplotypes_copy_type.tsv", sep = ""), quote = "none", na="")
