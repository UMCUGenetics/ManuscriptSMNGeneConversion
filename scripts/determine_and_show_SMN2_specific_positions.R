#### PSVs & SNVs ####

#load libraries
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
SNV_table = read_delim(args[1], delim = "\t")
PSVs = read_delim(args[2], delim = "\t")
prefix = args[3]


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


###DETERMINE SMN1/2 COPY-SPECIFIC VARIANTS###

#make table with variants divided by SMN1 and SMN2
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
  mutate(copy_specific = case_when(ratio_SMN1_haps >=0.9 & ratio_SMN2_haps <=0.2 ~ "SMN1-env",
                                   ratio_SMN2_haps >=0.9 & ratio_SMN1_haps <=0.2 ~ "SMN2-env"))

SNV_tally_SMN_filtered <- SNV_tally_SMN %>%
  filter(!is.na(copy_specific)) %>%
  filter(SMN1_tot_haps_non_NA >= 20 & SMN2_tot_haps_non_NA >= 20)

SMN2_specific_positions_bed <- SNV_tally_SMN_filtered %>%
  ungroup() %>%
  filter(copy_specific == "SMN2-env") %>% #filter for only SMN2-specific positions
  mutate(annotation = case_when(!is.na(PSV) ~ PSV,
                                is.na(PSV) ~ copy_specific)) %>%
  mutate(chr = "chr5") %>%
  mutate(start = POS - 1) %>%
  mutate(stop = POS) %>%
  select(chr, start, stop, annotation)

#write bed files with SMN2 specific positions
write.table(SMN2_specific_positions_bed, paste(prefix, "_SMN2_specific_positions_SMN2_0.9_SMN1_0.2.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep="\t", quote=FALSE)

#Look at SNVs at SMN2-specific positions, but for SMN2_specific_positions_bed
colnames(SMN2_specific_positions_bed) <- c("chr", "start", "POS", "annotation")
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

