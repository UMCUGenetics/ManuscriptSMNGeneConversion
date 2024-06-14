# ManuscriptSMNGeneConversion
Scripts and files related to the SMN gene conversion manuscript

# Make virtual env for python script
```bash
cd Scripts
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```



> Rscript for segmental duplications?

> how to make a masked reference genome

> how to plot results (R scripts?)

> script for annotation?

# Methods

## how to make the homopolymer BED file.
get_homopolymers_from_fasta.py 	create BED file with homopolymer region in the reference genome for length x-bp
```bash
get_homopolymers_from_fasta.py <path_to_reference_fasta> <output_file> --homopolymer_len <int>
```
--homopolymer_len = minimum length of homopolymer to consider in output [default 3]

## how to index the reference genome.
### .fai with samtools (tested version 1.15)
```bash
samtools faidx <path_to_reference_fasta>
```

### .dict with Picard (tested version 3.0.0)
```bash
picard CreateSequenceDictionary -R <path_to_reference_fasta> -O <path_to_reference_fasta>.dict
```

### .mmi with minimap2 (tested version 2.26)
(for rebasecalling only) .mmi using minimap2 -d (tested version 2.26)
```bash
minimap2 -d <path_to_reference_fasta>.mmi <path_to_reference_fasta>
```

## how to make the masked reference genome

TODO Demi

## Calculate statistics between ONT and Illumina data with loop_illumina_stats.py
```bash
loop_illumina_stats.py <input_folder_ont> <input_folder_illumina> <roi> <sample_file> <translation_file> <analysisID>
```

input_folder_ont = full path to folder containing all Clair3 VCF files from ONT analysis

input_folder_illumina = full path to folder containing all GATK ploidy VCF files from Illumina analysis

roi =  region of interest (chr:start-stop) i.e. the long-range PCR coordinates on the reference genome. Statistics will be calculated within this region.

sample_file = tsv file for renaming sampleIDs: Illumina sampleID, ONT sampleID, SMN1/2 ploidy

translation_file = tsv file with Illumina sampleID, ONT experiment sampleID, SMN1/2 ploidy

analysisID =  name for output folder


#### TODO > change default params to paper?
Optional
  --min_gq_ont MIN_GQ_ONT	minimum genotype quality variant [default 10]
  --min_dp_ont MIN_DP_ONT	minimum depth variant [default 4]
  --min_af_ont MIN_AF_ONT	minimum allele frequency variant [default list [0,0.5, 0.85]]
  --min_dp_ill MIN_DP_ILL	minimum depth variant [default 4]
  --min_af_ill MIN_AF_ILL       minimum allele frequency variant [default 0.85]

## Postprocessing data for manuscript results and figures
### Step 1 postprocessing:  vcf_parse_merge_depth.sh, to make TSV file for readdepth of each variant position in each haplotype for each sample

1) Loop over clair3 VCF files and make .tsv file
2) Make BED file for all variant positions
3) Calculate depth for all variant positions in BED file
4) Merge depth files into analysis specific .tsv files

Note: change path to binairy for samtools in the script: samtools=/path/to/samtoolsript assumes binairy of samtools in path (tested with samtools 1.17)
Alternatively change this to docker/singularity command.

```bash
vcf_parse_merge_depth.sh -o <path_to_output_folder> -i <path_to_input_folder>
```

path_to_output_folder = path to output folder

path_to_input_folder =  path to input folder. Input folder should contain SMA/ and 1000G/ folder. Within SMA/ and 1000/ specific samplefolder should excist. Within each samplefolder clair3/ and bam_files_haplotagged_split/ should be present. clair3/ folder includes the clair3 VCF and index, bam_files_haplotagged_split/ contains the haplotype specific BAM files + index.

Note:
* the script contains specific regex for SMA/1000G sampleID used in the manuscript.
* vcf_parse_merge_depth.sh will run vcf_parser.py and merging_variant_depth_files.py

### TODO path_to_output_folder/SMA/vcf_depth_merged_all_haps.tsv and path_to_output_folder/1000G/vcf_depth_merged_all_haps.tsv will be used in step 2?


### Step 2 postprocessing: "R-script om de excel met alle vairanten te maken"


2) Dan komt er een R script om de excel met alle vairanten te maken. Deze zal ik zelf aanpassen, want is nog redelijk rommelig en heb ik zelf geschreven (en er zit nog meer in dan alleen het maken van onderstaande file).
Script locatie: /hpc/dlab_haaften_mzwartkruis/Martin_C9/ONT_SMA_1000G_DATAFREEZE_20240223/analysis/SMA/SNV_analysis_copy_2.R
Excel met alle varianten:  /hpc/dlab_haaften_mzwartkruis/Martin_C9/ONT_SMA_1000G_DATAFREEZE_20240223/analysis/SMA/SNVs_pivoted_paraphase_suppl_made_in_R.xlsx
(vergelijkbaar voor 1000G)


### Step 3 postprocessing: create fasta sequence based on region-of-interest.

3) create_fasta_roi_SMA.sh -o <output_dir> -i <input_dir> -r <ROI>
*Part 1: make depth files for region of interest (contains all positions)
 i.      Iterate over bam files
 ii.      Run samtools depth over every bam file
 iii.      Add correct header and remove files without header
Part 2: prepare variant input file (excel) for the fasta maker script. Input file is now hardcoded, it’s better to make it into an argument
 i.      Script: prep_variant_input_fasta_maker.py
Part 3: make fasta files of all haplotypes
 i.      Iterate over variant files (one per haplotype)
 ii.      Per haplotype, create fasta file with creating_new_fasta_single.py
Part 4: concatenate all SMN1 fasta files and all SMN2 fasta files
  i.      This is based on the following file that is also made by the R script from step 2: list_haplotypes_copy_type.tsv (now hardcoded, best to make it into argument as well)
  ii.      Dit stuk heeft Demi geloof ik ook vooral geschreven dus kent zij beter dan ik ;)


