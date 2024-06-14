#!/bin/bash
#SBATCH --job-name=create_fasta_roi
#SBATCH --time=01:00:00
#SBATCH --gres=tmpspace:50G
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --ntasks=6
#SBATCH --mail-user=mzwartk2@umcutrecht.nl
#SBATCH --output=/hpc/dlab_haaften/mzwartkruis/Martin_C9/ONT_SMA_1000G_DATAFREEZE_20240223/slurm_output/create_fasta_roi_%j.out


# Stop execution of script whenever a non-zero exit status is returned
set -e

# Help function so other users can see what input is required
Help()
{
  # Display Help
  echo "Executing script to make fasta files for all haplotypes"
  echo
  echo
  echo "Syntax: Input requirements [-h:o:i:r]"
  echo "options:"
  echo "h     Print help"
  echo "o     Output directory, e.g. /hpc/dlab_haaften/mzwartkruis/Martin_C9/ONT_SMA_1000G_DATAFREEZE_20240223/analysis"
  echo "i     Input directory, requires both SMA and 1000G subdirectory, e.g. /hpc/dlab_haaften/mzwartkruis/Martin_C9/ONT_SMA_1000G_DATAFREEZE_20240223/data"
  echo "r     region of interest on chr5, e.g. 71274893-71447410"
}

### Files required ###
while getopts h:o:i:r: flag
do
  case "${flag}" in
    o) output_dir=${OPTARG};;
    i) input_dir=${OPTARG};;
    r) ROI=${OPTARG};;
    h) # display Help
        Help
        exit 1;;
    *) Help
        echo "Command not found"
        exit 1;;
  esac
done

#define tools
samtools="singularity exec -B /hpc/:/hpc/ -B $TMPDIR:$TMPDIR /hpc/dlab_haaften/dgommers/tools/samtools:1.17--hd87286a_1.img samtools"


#Part 1
#Creating depth files of region of interest
#SMA
mkdir -p "${output_dir}"
mkdir -p "${output_dir}/SMA/"
mkdir -p "${output_dir}/SMA/fasta_haplotypes"
mkdir -p "${output_dir}/SMA/fasta_haplotypes/depth_files"

echo "Creating depth files per sample haplotype".
echo ${input_dir}/SMA/{SMA,HG}*/bam_files_haplotagged_split
# Iterate over directories
for dir in ${input_dir}/SMA/{SMA,HG}*/bam_files_haplotagged_split; do
  if [ -d "${dir}" ]; then
    echo "Processing files in directory: ${dir}"

    # Iterate over bam files of specified phasing method in the directory
    for bam_file_path in ${dir}/*.bam; do
      if [ -f "${bam_file_path}" ]; then
      echo "Processing file: ${bam_file_path}"
      bam_file=$(basename $bam_file_path)
      sample_name=$(echo $bam_file | grep -Eo "HG[0-9]+|(SMA[0-9]+)((_blood)?(_fib_P[0-9]+)*(_P[0-9]+)?)")
      echo "${sample_name}"
      hap=$(echo $bam_file | grep -Eo "hap[0-9]+")
      echo "${hap}"


      #Run samtools depth for every bam file
      $samtools depth -J -H --min-MQ 5 -aa -r "chr5:${ROI}" -o "${output_dir}/SMA/fasta_haplotypes/depth_files/${sample_name}_${hap}_depth_all_pos_roi.bed" ${bam_file_path}
      echo "Samtools depth for ${sample_name}_${hap} completed on ROI"
      #add correct header
      awk 'BEGIN {OFS="\t"} NR==1 {print "CHROM", "POS", "depth"; next} {print}' "${output_dir}/SMA/fasta_haplotypes/depth_files/${sample_name}_${hap}_depth_all_pos_roi.bed" > "${output_dir}/SMA/fasta_haplotypes/depth_files/${sample_name}_${hap}_depth_all_pos_roi_header.bed"
      #remove files without header
      rm "${output_dir}/SMA/fasta_haplotypes/depth_files/${sample_name}_${hap}_depth_all_pos_roi.bed"

      fi
    done
  else
  echo "not right directory given"
  fi
done


#Part 2
#Preparing variant input files for fasta maker script
mkdir -p "${output_dir}/SMA/fasta_haplotypes/pos_ref_alt_files"
cd "${output_dir}/SMA/fasta_haplotypes/pos_ref_alt_files"

/hpc/dlab_haaften/mzwartkruis/Martin_C9/ONT_SMA_1000G_DATAFREEZE_20240223/scripts/./prep_variant_input_fasta_maker.py ${output_dir}/SMA/SNVs_pivoted_paraphase_suppl_made_in_R.xlsx pos_ref_alt_

#Part 3
#making fasta files of all haplotypes
chr5_fa="/hpc/dlab_haaften/mzwartkruis/Martin_C9/ONT_SMA_1000G_DATAFREEZE_20240223/scripts/chr5.fa"

mkdir -p "${output_dir}/SMA/fasta_haplotypes/fasta"
mkdir -p "${output_dir}/SMA/fasta_haplotypes/log"

# Iterate over haplotypes
for alt_file_path in ${output_dir}/SMA/fasta_haplotypes/pos_ref_alt_files/*.tsv; do
  if [ -f "${alt_file_path}" ]; then
    echo "Processing file: ${alt_file_path}"
    alt_file=$(basename $alt_file_path)
    sample_name=$(echo $alt_file | grep -Eo "HG[0-9]+|(SMA[0-9]+)((_blood)?(_fib_P[0-9]+)*(_P[0-9]+)?)")
    echo "${sample_name}"
    hap=$(echo $alt_file | grep -Eo "hap[0-9]+")
    echo "${hap}"
    echo "Creating a fasta file for ${sample_name} ${hap}"
    # Execute python script
    python3 /hpc/dlab_haaften/mzwartkruis/Martin_C9/ONT_SMA_1000G_DATAFREEZE_20240223/scripts/./creating_new_fasta_single.py ${chr5_fa} ${alt_file_path} ${ROI} ${output_dir}/SMA/fasta_haplotypes/fasta/${sample_name}_${hap}.fa ${output_dir}/SMA/fasta_haplotypes/log/${sample_name}_${hap}.txt ${output_dir}/SMA/fasta_haplotypes/depth_files/*${sample_name}_${hap}*.bed ${sample_name}_${hap}
    #usage: /hpc/dlab_haaften/dgommers/scripts/./creating_new_fasta_demi.py <reference_fasta> <pos_ref_alt_file.tsv> <ROI, eg. 71274893-71447410> <output_name.fa> <logfile_name.txt> <depth_file.bed> <sample_hap>
  fi
done


#Part 4
#Concatenate fasta files per copy type (SMN1 or SMN2)
copy_type_file=/hpc/dlab_haaften/mzwartkruis/Martin_C9/ONT_SMA_1000G_DATAFREEZE_20240223/analysis/SMA/list_haplotypes_copy_type.tsv

#exclude header
tail -n +2 ${copy_type_file} > ${output_dir}/SMA/fasta_haplotypes/fasta/sample_names_file.txt

#Making sure that there are no empty lines in the file
tr -d '\r' < ${output_dir}/SMA/fasta_haplotypes/fasta/sample_names_file.txt > ${output_dir}/SMA/fasta_haplotypes/fasta/sample_names_file_unix.txt

# Make directories for SMN1 and SMN2
mkdir -p ${output_dir}/SMA/fasta_haplotypes/fasta/SMN1
mkdir -p ${output_dir}/SMA/fasta_haplotypes/fasta/SMN2
mkdir -p ${output_dir}/SMA/fasta_haplotypes/fasta/merged

while IFS= read -r line; do
  echo $line

  # Extract sample name
  sample=$(echo "$line" | awk '{print $1}')
  echo ${sample}

  # Extract smn copy
	smn_copy=$(echo "$line" | awk '{print $2}')
  echo ${smn_copy}


  if [ -n "$sample" ] && [ -n "$smn_copy" ] && [ "$smn_copy" != "NA" ]; then
    cp "${output_dir}/SMA/fasta_haplotypes/fasta/${sample}.fa" "${output_dir}/SMA/fasta_haplotypes/fasta/${smn_copy}"
  else
    echo "Invalid sample or SMN copy."
  fi


done < ${output_dir}/SMA/fasta_haplotypes/fasta/sample_names_file_unix.txt

# Removing merged fastas if they already exist
#if [ -d ${output_dir}/SMA/fasta_haplotypes/fasta/merged/SMN1_haplotypes.fa]; then
#  rm ${output_dir}/SMA/fasta_haplotypes/fasta/merged/SMN2_haplotypes.fa
#else
#  continue
#fi

# Concatenating into haplotype files
cat ${output_dir}/SMA/fasta_haplotypes/fasta/SMN1/* >> ${output_dir}/SMA/fasta_haplotypes/fasta/merged/SMN1_haplotypes.fa
cat ${output_dir}/SMA/fasta_haplotypes/fasta/SMN2/* >> ${output_dir}/SMA/fasta_haplotypes/fasta/merged/SMN2_haplotypes.fa

echo "End of script"
