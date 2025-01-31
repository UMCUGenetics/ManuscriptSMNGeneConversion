# ManuscriptSMNGeneConversion
Scripts and files related to the SMN gene conversion manuscript. The preprint of this manuscript can currently be found at: https://doi.org/10.1101/2024.07.16.24310417

# Download repository
To download the main branch of the repository:
```bash
git clone https://github.com/UMCUGenetics/ManuscriptSMNGeneConversion.git
```
To download a specific branch:
```bash
git clone -b release_v1.1.0 https://github.com/UMCUGenetics/ManuscriptSMNGeneConversion.git
```

# Make virtual env for python script
Note: this was developed and tested with Python 3.6.8. Using other versions of Python might give errors.
```bash
cd ManuscriptSMNGeneConversion/scripts
python3.6 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

# Methods

## 1) How to make the homopolymer BED file.
get_homopolymers_from_fasta.py 	create BED file with homopolymer regions in the reference genome for length x-bp
```bash
source <workflow_folder>/scripts/venv/bin/activate
python get_homopolymers_from_fasta.py <path_to_reference_fasta> <output_file> --homopolymer_len <int>
```
* --homopolymer_len = minimum length of homopolymer to consider in output [default 3]

## 2) How to make a masked reference:

### 2.1 Determining masking coordinates with segmental duplication analysis
Segmental duplication analysis was performed by running script: 2.1_numcer_analysis.sh. 
This script requires installing the [MUMmer (v4.0.0)](https://mummer4.github.io/manual/manual.html) system. Adjust path of mummer in the script for personal use.

The script runs `nucmer` for genome alignment, `delta-filter` for filtering the delta file for 95% identity and alignment length of 10kb, and `show-coords` for converting the delta file into a coordinates file. 
Columns 1, 2 and 8 are extracted to create a BED file, where bed segment is named by reference-query alignment (a segments) and by query-reference alignment (b segments). 

```bash
scripts/./2.1_nucmer_analysis.sh -o <output_name> -r <reference> -q <query>
```
* output_name = output name for the files this script outputs
* reference = reference fasta (we used chr5.fa of T2T-CHM13)
* query = query fasta (we used chr5.fa of T2T-CHM13)

Coordinate for masking can be extracted from the BED file into a new BED file and used in step 2.2 Masked reference genome.

### 2.2 Masked reference genome
Tested with BEDtools v2.25.0
```bash
bedtools maskfasta -fi <reference_genome> -bed <mask_bed_file> -fo <output_file>
```
* reference_genome = full path to reference genome .fa(sta)
* mask_bed_file = BED file with-to-mask regions (.e.g. <repo_folder>/datafiles/masking_approach_mask_coords_20a_22a.bed as used in the manuscript)
* output = output name for masker reference genome (e.g. T2T-CHM13 v2.0/hs1 as used in the manuscript)

## 3) Calculate statistics between ONT and Illumina data with loop_illumina_stats.py
This script calculates sensitivity and precision for SNV (only) between Illumina (GATK ploidy VCF) and ONT data (Clair3 haplotype variant calling).

```bash
source <workflow_folder>/scripts/venv/bin/activate
python loop_illumina_stats.py <input_folder_ont> <input_folder_illumina> <roi> <sample_file> <translation_file> <analysisID>
```
* input_folder_ont = full path to folder containing all Clair3 VCF files from ONT analysis
* input_folder_illumina = full path to folder containing all GATK ploidy VCF files from Illumina analysis
* roi =  region of interest (chr:start-stop) i.e. the long-range PCR coordinates on the reference genome. Statistics will be calculated within this region.
* sample_file = tsv file for renaming sampleIDs: sampleID, Illumina sampleID, SMN1/2 ploidy
* translation_file = tsv file with sampleID, ONT sampleID, SMN1/2 ploidy
* analysisID =  name for output folder

## 4) Postprocessing data for manuscript results and figures

To postprocess the SMA and 1000G data for the manuscript, seven main scripts are used.

1) 4.1_select_bed_or_region_phasing.sh
2) 4.2_vcf_parse_merge_depth.sh was used to create a tsv file containing the Clair3 variants and read depth on eacht variant position in each haplotype for each sample
3) 4.3_SNV_analysis.R was used to to determine SMN_copy_type (SMN1 or SMN2) based on PSV13 and output TSV file with all called variants
4) 4.4_split_reference_genome.py was used to slice the reference genome for the contig of interest (e.g. chr5 for SMA)
5) 4.5_create_fasta_roi.sh was used to create fasta sequences of each haplotype based on original reference contig and detected variants.
6) 4.6_determine_and_show_SMN_specific_positions.R was used to determine SMN1/SMN2 specific positions.
7) 4.7_load_bed_and_show_SMN_specific_positions.R  was used to determine SMN1/SMN2 specific positions based on a BED input file.

### 4.1) 4.1_select_bed_or_region_phasing.sh
This script selects the best phasing method ('bed' or 'region') based on read depth on PSV positions, and moves the unselected files into a separate folder.

```bash
sh 4.1_select_bed_or_region_phasing.sh -i <path_to_input_folder> -p <path_to_PSV_bed_file>
```

* path_to_input_folder = path to input folder. Within this folder, specific samplefolders should exist (sample names starting with 'SMA' or 'HG'). Within each samplefolder, bam_files_haplotagged/, bam_files_haplotagged_split/, clair3/, sniffles2/ and vcf/ folders should be present (i.e. the output of HapSMA).
* path_to_PSV_bed_file = path to a bed file containing PSV positions. See 'PSV_SMN1_minus_PSV8_liftover_hg19_to_T2T_CHM13.bed' in the <repo_folder>/datafiles folder.

The unselected files will be moved into: {input_dir}/phasing_not_selected/ 
The user can keep this folder or remove it if desired.

### 4.2) 4.2_vcf_parse_merge_depth.sh
Make TSV file for read depth of each variant position in each haplotype for each sample.

1) Loop over clair3 VCF files and make .tsv file
2) Make BED file for all variant positions
3) Calculate depth for all variant positions in BED file
4) Merge depth files into analysis specific .tsv files

Note:
* change /path/to/samtools in vcf_parse_merge_depth.sh to excecutable/binary or docker/singularity command of samtools (tested with samtools 1.17).
* make sure the python virtual environment has been made in the repo folder (see above: Make virtual env for python script).
* the script contains specific regex for SMA/1000G sampleIDs used in the manuscript. These need to be changed if other sampleID are used.
* 4.2_vcf_parse_merge_depth.sh will run the scripts 4.2a_vcf_parser.py and 4.2b_merging_variant_depth_files.py from the repo folder.

```bash
sh vcf_parse_merge_depth.sh -o <path_to_output_folder> -i <path_to_input_folder> -s <path_to_repo_folder>
```
* path_to_output_folder = path to output folder
* path_to_input_folder =  path to input folder. Input folder should contain SMA/ and 1000G/ folder. Within SMA/ and 1000/ specific samplefolder should exist. Within each samplefolder clair3/ and bam_files_haplotagged_split/ should be present. clair3/ folder includes the clair3 VCF and index, bam_files_haplotagged_split/ contains the haplotype specific BAM files + index.
* path_to_repo_folder = path to repo folder containing the scripts, e.g. ManuscriptSMNGeneConversion/scripts.

4.2_vcf_parse_merge_depth.sh will produce SMA/vcf_depth_merged_all_haps.tsv and 1000G/vcf_depth_merged_all_haps.tsv TSV files that will be the input file of step 4.3, 4.6, and 4.7.


### 4.3) 4.3_SNV_analysis.R

Determine SMN_copy_type based on PSV13 and output a file with all variants.

Note: script was runned and tested using rocker tidyverse v4.4 image.

```bash
Rscript 4.3_SNV_analysis.R <input_SNV_table> <PSV_file> <prefix>
```
* input_SNV_table = .tsv file (vcf_depth_merged_all_haps.tsv) produced in step 4.2 (4.2_vcf_parse_merge_depth.sh)
* PSV_file = PSV positions for the use reference genome. See repo/datafiles/PSV_liftover_hg19_to_T2T_CHM13.txt for T2T-CHM13 positions.
* prefix (e.g. SMA/1000G)

output files will be stored in working directory.

The output of this script will result in two .tsv files:
* {prefix}_list_haplotypes_copy_type.tsv     outputs SMN1, SMN2, of NA of SMN_copy_type for each sample based on PSV13
* {prefix}_variants_pivoted_supplementary.tsv   table of variants for each position for each sample
These output TSV files will be used in step 4.5.

### 4.4) 4.4_split_reference_genome.py
Slice reference genome for contig of interest.

```bash
source <workflow_folder>/scripts/venv/bin/activate
python 4.4_split_reference_genome.py <path_to_fasta> <contig>
```
* path_to_fasta = full path the reference genome fa/fasta file
* contig = ID of contig that needs to be in the output (e.g. chr5).
* output will be written to {contig}.fa


### 4.5) 4.5_create_fasta_roi.sh

Create fasta sequence of haplotype based on original contig and detected genetic variants within the region-of-interest.

Note:
* change /path/to/samtools in vcf_parse_merge_depth.sh to excecutable/binary or docker/singularity command of samtools (tested with samtools 1.17).
* make sure the python virtual environment has been made in the repo folder (see above: Make virtual env for python script)
* the script contains specific regex for SMA/1000G sampleIDs used in the manuscript. These need to be changed if other sampleID are used.
* 4.5_create_fasta_roi.sh will run the scripts 4.5a_prep_variant_input_fasta_maker.py and 4.5b_create_new_fasta_single.py from the repo folder.

1) Make depth files for region of interest for each haplotype specific BAM using samtools
2) Convert output of step 4.3_SNV_analysis.R using 4.5a_prep_variant_input_fasta_maker.py and output one file per haplotype
3) Make variant correct reference sequences (FASTA) for each haplotype for each sample using 4.5b_create_new_fasta_single.py
4) Concatenate all FASTA files into SMN1 or SMN2 haplotype output files.

```bash
sh 4.5_create_fasta_roi.sh -o <output_dir> -i <input_dir> -f <contig_fasta> -r <ROI> -c <copy_type_file> -p <pivot_file> -a <analysis> -s <path_to_repo_folder>
```
* output_dir =  full path to output folder.
* input_dir = full path to input folder. This should be same input folder as used in step 4.2.
* contig_fasta = full path to fa(sta) file from the selected contig (e.g. chr5). This is the output file from step 4.4 (e.g. chr5.fa).
* ROI = region of interest. e.g. 71274893-71447410.
* copy_type_file = full path to copy_type file as generated in step 4.3 (e.g. SMA_list_haplotypes_copy_type.tsv).
* pivot_file = full path to pivoted file as generated in step 4.3  (e.g. SMA_variants_pivoted_supplementary.tsv).
* analysis = type of analysis, e.g. SMA or 1000G.
* path_to_repo_folder = path to repo folder containing the scripts, e.g. ManuscriptSMNGeneConversion/scripts.


### 4.6) 4.6_determine_and_show_SMN_specific_positions.R

Determine SMN1/SMN2 specific positions. This script should be run on control samples containing SMN1 and SMN2 (ideally 2xSMN1 and 2xSMN2 per sample). Criteria for calling a position as SMN2-specific are that the variant is present in at least 90% of the called SMN2 haplotypes and at maximum in 10% of the called SMN1 haplotypes, based on at least 20 called haplotypes for both SMN1 and SMN2. These parameters can be modified, e.g. using less strict cutoffs for discovery purposes.

Note: script was runned and tested using rocker tidyverse v4.4 image.

```bash
Rscript 4.6_determine_and_show_SMN_specific_positions.R <input_SNV_table> <PSV_file> <prefix>
```

* input_SNV_table = .tsv file (vcf_depth_merged_all_haps.tsv) produced in step 4.2 (4.2_vcf_parse_merge_depth.sh)
* PSV_file = PSV positions for the used reference genome. See <repo_folder>/datafiles/PSV_liftover_hg19_to_T2T_CHM13.txt for CHM13 positions.
* prefix (e.g. SMA/1000G)

The output of this script will result in three output files:
* {prefix}_SNV_tally_at_PSVs_and_SMN1-2_specific_positions.tsv   TSV file containing all SMN1- or SMN2-specific positions including variant counts and ratios in SMN1 and SMN2 haplotypes
* {prefix}_PSVs_and_SMN1-2_specific_positions.bed   BED file containing all SMN1- or SMN2-specific positions
* {prefix}_SNVs_at_SMN1-2_specific_positions.tsv   TSV file containing variant calls (1 for SMN1 environment SNV, 2 for SMN2 environment SNV) at SMN1/2-specific variant positions per haplotype.

### 4.7) 4.7_load_bed_and_show_SMN_specific_positions.R

Determine variants at SMN1/SMN2 specific positions based on a specific BED input file (generated by 4.6_determine_and_show_SMN_specific_positions.R).

Note: script was runned and tested using rocker tidyverse v4.4 image.

```bash
Rscript 4.7_load_bed_and_show_SMN_specific_positions.R <input_SNV_table> <PSV_file> <SMN_specific_positions_bed> <prefix>
```

* input_SNV_table = .tsv file (vcf_depth_merged_all_haps.tsv) produced in step 4.2 (4.2_vcf_parse_merge_depth.sh)
* PSV_file = PSV positions for the used reference genome. See <repo_folder>/datafiles/PSV_liftover_hg19_to_T2T_CHM13.txt for CHM13 positions.
* SMN_specific_positions_bed = {prefix}_PSVs_and_SMN1-2_specific_positions.bed file created in step 4.6 (e.g <repo_folder>/datafiles/1000G_PSVs_and_SMN1-2_specific_positions.bed
* prefix = prefix (e.g. SMA/1000G)

The output of this script will result in a .tsv output file:
* {prefix}_SNVs_at_SMN1-2_specific_positions.tsv    tsv table containing variant calls (1 for SMN1 environment SNV, 2 for SMN2 environment SNV) at SMN1/2-specific variant positions per haplotype.
