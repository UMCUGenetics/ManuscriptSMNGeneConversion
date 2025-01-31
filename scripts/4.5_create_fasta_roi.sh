#!/bin/bash
#SBATCH --job-name=create_fasta_roi
#SBATCH --time=01:00:00
#SBATCH --gres=tmpspace:50G
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=6
#SBATCH --export=NONE
#SBATCH --output=create_fasta_roi_%j.out


# Stop execution of script whenever a non-zero exit status is returned
set -e

# Help function so other users can see what input is required
Help()
{
  # Display Help
  echo "Executing script to make fasta files for all haplotypes"
  echo
  echo
  echo "Syntax: Input requirements [-h:o:i:f:r:c:p:a:s]"
  echo "options:"
  echo "h     print help"
  echo "o     output directory"
  echo "i     input directory, requires both SMA and/or 1000G subdirectory"
  echo "f     full path to FASTA of reference contig, e.g. chr5.fa"
  echo "r     region of interest on contig, e.g. 71274893-71447410"
  echo "c     full path to copy_type_file"
  echo "p     full path to SNVs_pivoted_paraphase file"
  echo "a     dataset analysis type (e.g. SMA or 1000G)"
  echo "s     Repo folder, containing the scripts"
}

### Files required ###
while getopts h:o:i:f:r:c:p:a:s: flag
do
  case "${flag}" in
    o) output_dir=${OPTARG};;
    i) input_dir=${OPTARG};;
    f) contig_fasta=${OPTARG};;
    r) ROI=${OPTARG};;
    c) copy_type_file=${OPTARG};;
    p) pivot_file=${OPTARG};;
    a) analysis=${OPTARG};;
    s) repo_folder=${OPTARG};;
    h) # display Help
        Help
        exit 1;;
    *) Help
        echo "Command not found"
        exit 1;;
  esac
done

#define tools
samtools=/path/to/samtools


## Load virtual venv
source ${repo_folder}/venv/bin/activate

#Part 1
#Creating depth files of region of interest
mkdir -p "${output_dir}/${analysis}/fasta_haplotypes/depth_files"

echo "Creating depth files per sample haplotype"
echo ${input_dir}/${analysis}/{SMA,HG}*/bam_files_haplotagged_split


# Iterate over directories
for dir in ${input_dir}/${analysis}/{SMA,HG}*/bam_files_haplotagged_split; do
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
      $samtools depth -J -H --min-MQ 5 -aa -r "chr5:${ROI}" -o "${output_dir}/${analysis}/fasta_haplotypes/depth_files/${sample_name}_${hap}_depth_all_pos_roi.bed" ${bam_file_path}
      echo "samtools depth for ${sample_name}_${hap} completed on ROI"
      #add correct header
      awk 'BEGIN {OFS="\t"} NR==1 {print "CHROM", "POS", "depth"; next} {print}' "${output_dir}/${analysis}/fasta_haplotypes/depth_files/${sample_name}_${hap}_depth_all_pos_roi.bed" > "${output_dir}/${analysis}/fasta_haplotypes/depth_files/${sample_name}_${hap}_depth_all_pos_roi_header.bed"
      #remove files without header
      rm "${output_dir}/${analysis}/fasta_haplotypes/depth_files/${sample_name}_${hap}_depth_all_pos_roi.bed"

      fi
    done
  else
  echo "wrong directory provided as input"
  fi
done


#Part 2
mkdir -p "${output_dir}/${analysis}/fasta_haplotypes/pos_ref_alt_files"
cd "${output_dir}/${analysis}/fasta_haplotypes/pos_ref_alt_files"


$repo_folder/4.5a_prep_variant_input_fasta_maker.py $pivot_file pos_ref_alt_

#Part 3
#making fasta files of all haplotypes

mkdir -p "${output_dir}/${analysis}/fasta_haplotypes/fasta"
mkdir -p "${output_dir}/${analysis}/fasta_haplotypes/log"

# Iterate over haplotypes
for alt_file_path in ${output_dir}/${analysis}/fasta_haplotypes/pos_ref_alt_files/*.tsv; do
  if [ -f "${alt_file_path}" ]; then
    echo "Processing file: ${alt_file_path}"
    alt_file=$(basename $alt_file_path)
    sample_name=$(echo $alt_file | grep -Eo "HG[0-9]+|(SMA[0-9]+)((_blood)?(_fib_P[0-9]+)*(_P[0-9]+)?)")
    echo "${sample_name}"
    hap=$(echo $alt_file | grep -Eo "hap[0-9]+")
    echo "${hap}"
    echo "Creating a fasta file for ${sample_name} ${hap}"
    # Execute python script
    python3 $repo_folder/4.5b_create_new_fasta_single.py ${contig_fasta} ${alt_file_path} ${ROI} ${output_dir}/${analysis}/fasta_haplotypes/fasta/${sample_name}_${hap}.fa ${output_dir}/${analysis}/fasta_haplotypes/log/${sample_name}_${hap}.txt ${output_dir}/${analysis}/fasta_haplotypes/depth_files/*${sample_name}_${hap}*.bed ${sample_name}_${hap}
  fi
done


#Part 4
#Concatenate fasta files per copy type (SMN1 or SMN2)

#exclude header
tail -n +2 ${copy_type_file} > ${output_dir}/${analysis}/fasta_haplotypes/fasta/sample_names_file.txt

#Making sure that there are no empty lines in the file
tr -d '\r' < ${output_dir}/${analysis}/fasta_haplotypes/fasta/sample_names_file.txt > ${output_dir}/${analysis}/fasta_haplotypes/fasta/sample_names_file_unix.txt

# Make directories for SMN1 and SMN2
mkdir -p ${output_dir}/${analysis}/fasta_haplotypes/fasta/SMN1
mkdir -p ${output_dir}/${analysis}/fasta_haplotypes/fasta/SMN2
mkdir -p ${output_dir}/${analysis}/fasta_haplotypes/fasta/merged

while IFS= read -r line; do
  echo $line

  # Extract sample name
  sample=$(echo "$line" | awk '{print $1}')
  echo ${sample}

  # Extract smn copy
	smn_copy=$(echo "$line" | awk '{print $2}')
  echo ${smn_copy}


  if [ -n "$sample" ] && [ -n "$smn_copy" ] && [ "$smn_copy" != "NA" ]; then
    cp "${output_dir}/${analysis}/fasta_haplotypes/fasta/${sample}.fa" "${output_dir}/${analysis}/fasta_haplotypes/fasta/${smn_copy}"
  else
    echo "Invalid sample or SMN copy."
  fi


done < ${output_dir}/${analysis}/fasta_haplotypes/fasta/sample_names_file_unix.txt

# Concatenating into haplotype files
cat ${output_dir}/${analysis}/fasta_haplotypes/fasta/SMN1/* >> ${output_dir}/${analysis}/fasta_haplotypes/fasta/merged/SMN1_haplotypes.fa
cat ${output_dir}/${analysis}/fasta_haplotypes/fasta/SMN2/* >> ${output_dir}/${analysis}/fasta_haplotypes/fasta/merged/SMN2_haplotypes.fa

echo "End of script"
