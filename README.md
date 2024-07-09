# ManuscriptSMNGeneConversion
Scripts and files related to the SMN gene conversion manuscript

# Make virtual env for python script
```bash
cd scripts
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

# Methods

## 1) How to make the homopolymer BED file.
get_homopolymers_from_fasta.py 	create BED file with homopolymer region in the reference genome for length x-bp
```bash
source <workflow_folder>/scripts/venv/bin/activate
python get_homopolymers_from_fasta.py <path_to_reference_fasta> <output_file> --homopolymer_len <int>
```
* --homopolymer_len = minimum length of homopolymer to consider in output [default 3]

## 2) How to make a masked reference:
Tested with BEDtools v2.25.0
```bash
bedtools maskfasta -fi <reference_genome> -bed <mask_bed_file> -fo <output_file>
```
* reference_genome = full path to reference genome .fa(sta)
* mask_bed_file = BED file with-to-mask regions (.e.g. ../datafiles/masking_approach_mask_coords_20a_22a.bed as used in the manuscript)
* output = output name for masker reference genome (e.g. T2T-CHM13 v2.0/hs1 as used in the manuscript)

## 3) Calculate statistics between ONT and Illumina data with loop_illumina_stats.py
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

To postprocess the SMA and 1000G data for the manuscript 4 main scripts are used.
1) vcf_parse_merge_depth.sh was used to create a tsv file based on the readdepth of each variant position in each haplotype for each sample
2) SNV_analysis_paraphase.R was used to to determine SMN_copy_type based on PSV13 and output TSV file with variants
3) split_reference_genome.py was used to slice the reference genome for the contig of interest (e.g. chr5 for SMA)
4) create_fasta_roi_SMA.sh was used to create fasta sequences of each haplotype based on original reference contig and detected SNVs.

### 4.1) vcf_parse_merge_depth.sh, to make TSV file for readdepth of each variant position in each haplotype for each sample

1) Loop over clair3 VCF files and make .tsv file
2) Make BED file for all variant positions
3) Calculate depth for all variant positions in BED file
4) Merge depth files into analysis specific .tsv files

Note: change path to binairy for samtools in the script: samtools=/path/to/samtoolsript assumes binairy of samtools in path (tested with samtools 1.17)
Alternatively change this to docker/singularity command.

```bash
vcf_parse_merge_depth.sh -o <path_to_output_folder> -i <path_to_input_folder>
```
* path_to_output_folder = path to output folder
* path_to_input_folder =  path to input folder. Input folder should contain SMA/ and 1000G/ folder. Within SMA/ and 1000/ specific samplefolder should excist. Within each samplefolder clair3/ and bam_files_haplotagged_split/ should be present. clair3/ folder includes the clair3 VCF and index, bam_files_haplotagged_split/ contains the haplotype specific BAM files + index.

Note:
* the script contains specific regex for SMA/1000G sampleIDs used in the manuscript. These need to be changed if other sampleID are used.
* vcf_parse_merge_depth.sh will run vcf_parser.py and merging_variant_depth_files.py

vcf_parse_merge_depth.sh will produce SMA/vcf_depth_merged_all_haps.tsv TSV file that will be the input file of step 2.


### 4.2) SNV_analysis_paraphase.R to determine SMN_copy_type based on PSV13 and output TSV file with variants

Note: script was runned and tested using rocker tidyverse v4.4 image.

```bash
Rscript SNV_analysis_paraphase.R <input_SNV_table> <PSV_file>
```
* input_SNV_table = .tsv file (vcf_depth_merged_all_haps.tsv) produced in step 1 (vcf_parse_merge_depth.sh)
* PSV_file = PSV positions for the use reference genome. See repo/datafiles/PSV_liftover_hg19_to_T2T_CHM13.txt for CHM13 positions.


### 4.3) split_reference_genome.py, slice reference genome for contig of interest.
```bash
source <workflow_folder>/scripts//venv/bin/activate
python split_reference_genome.py <path_to_fasta> <contig>
```
* path_to_fasta = full path the reference genome fa/fasta file
* contig = ID of contig that needs to be in the output (e.g. chr5).
* output will be written to {contig}.fa


### 4.4)create_fasta_roi_SMA.sh, create fasta sequence of haplotype based on original contig and detected SNVs within the region-of-interest.

Note: change path to binairy for samtools in the script: samtools=/path/to/samtoolsript assumes binairy of samtools in path (tested with samtools 1.17)\
Alternatively change this to docker/singularity command.


1) Make depth files for region of interest for each haplotype specific BAM using samtools
2) Convert output of step 4.1 using prep_variant_input_fasta_maker.py and output files for step 4.3
3) Make variant correct reference sequences (FASTA) for each haplotype for each sample using creating_new_fasta_single.py
4) Concatenate all FASTA into SMN1 or SMN2 haplotype output files.

```bash
create_fasta_roi_SMA.sh -o <output_dir> -i <input_dir> -f <contig_fasta> -r <ROI> -c <copy_type_file> -p <pivot_file> -a <analysis>
```
* output_dir =  full path to output folder.
* input_dir = full path to input folder. This should be the output folder used in step 1.
* contig_fasta = full path to fa(sta) file from the selected contig (e.g. chr5). This is the output file from step 3 (e.g. chr5.fa).
* ROI = region of interest. e.g. 71274893-71447410.
* copy_type_file = full path to copy_type file as generated in step 2 (e.g. SMA_list_haplotypes_copy_type.tsv).
* pivot_file = full path to pivoted file as generated in step 2  (e.g. SMA_SNVs_pivoted_paraphase_suppl_made_in_R.tsv).
* analysis = type of analysis, e.g. SMA or 1000G.

