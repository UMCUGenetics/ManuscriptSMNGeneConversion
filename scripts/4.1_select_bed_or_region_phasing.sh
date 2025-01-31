#!/bin/bash
#SBATCH --job-name=select_bed_or_region
#SBATCH --time=03:00:00
#SBATCH --gres=tmpspace:40G
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=8
#SBATCH --output=select_bed_or_region_%j.out
# ---------------------------------------------------------------------------
# [Maria Zwartkruis]
# select_bed_or_region
# ---------------------------------------------------------------------------


# Settings
samtools=/path/to/samtools

# Stop execution of script whenever a non-zero exit status is returned
set -e

# Help function so other users can see what input is required
Help()
{
  # Display Help
  echo "This script selects the best phasing method ('bed' or 'region') and moves the unselected files into a separate folder"
  echo
  echo "Syntax: Input requirements [-h:i:p]"
  echo "options:"
  echo "h     Print help"
  echo "i     Input directory, containing the HapSMA results in one directory per sample (sample names have to start with SMA or HG)"
  echo "p     Path to PSV bed file"
}

### Files required ###
while getopts h:i:p: flag
do
  case "${flag}" in
    i) input_dir=${OPTARG};;
    p) PSV_bed_file=${OPTARG};;
    h) # display Help
        Help
        exit 1;;
    *) Help
        echo "Command not found"
        exit 1;;
  esac
done

## Echo the input parameters
echo "Input directory is ${input_dir}"
echo "PSV bed file is ${PSV_bed_file}"

# Make a directory to move the old phasing files to
mkdir -p ${input_dir}/select_bed_or_region

# Iterate over directories
for dir in ${input_dir}/{SMA,HG}*/bam_files_haplotagged_split; do
  if [ -d "${dir}" ]; then
    echo "Processing files in directory: ${dir}"
    sample_name=$(echo ${dir} | grep -Eo "HG[0-9]+|(SMA[0-9]+)((_blood)?(_fib_P[0-9]+)*(_P[0-9]+)?)" | uniq)
    echo "${sample_name}"

    # Iterate over bed phasing bam files in the directory
    for bam_file_path_bed in ${dir}/*bed*.bam; do
      if [ -f "${bam_file_path_bed}" ]; then
      echo "Processing file: ${bam_file_path_bed}"
      echo "${sample_name}"
      hap=$(echo $bam_file_path_bed | grep -Eo "hap[0-9]+")
      echo "${hap}"

      # Create zero-depth entries from the BED file
      awk '{for (i=$2+1; i<=$3; i++) print $1 "\t" i "\t" 0}' "${PSV_bed_file}" > "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_prefilled_depth_0.bed"

      # Run samtools depth for every bam file
      $samtools depth -@ 8 -a -b "${PSV_bed_file}" -o "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_actual_depth.bed" ${bam_file_path_bed}
      echo "Samtools depth for ${sample_name}_${hap}_bed completed"

      # Merge prefilled and actual depth file and remove duplicates, prioritizing actual depth values.
      #cat "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_actual_depth.bed" "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_prefilled_depth_0.bed" | sort -k1,1 -k2,2n | awk '!seen[$1,$2]++' > "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_depth.bed"

      cat "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_actual_depth.bed" "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_prefilled_depth_0.bed" \
    | sort -k1,1 -k2,2n \
    | awk '{
        key = $1 FS $2;
        if (!(key in seen)) {
            seen[key] = $0;
        } else if ($3 > 0) {
            seen[key] = $0;
        }
    } END {
        for (key in seen) {
            print seen[key];
        }
    }' > "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_depth.bed"

      # Clean old files
      rm "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_prefilled_depth_0.bed" "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_actual_depth.bed"

      # Merge depth files of all haplotypes
      cat "${input_dir}/select_bed_or_region/${sample_name}_${hap}_bed_depth.bed" >> "${input_dir}/select_bed_or_region/${sample_name}_merged_bed_depth.bed"
      echo "Depth files for ${sample_name} have been merged."

      fi
    done

    # Iterate over region phasing bam files in the directory
    for bam_file_path_region in ${dir}/*region*.bam; do
      if [ -f "${bam_file_path_region}" ]; then
      echo "Processing file: ${bam_file_path_region}"
      echo "${sample_name}"
      hap=$(echo $bam_file_path_region | grep -Eo "hap[0-9]+")
      echo "${hap}"

      # Create zero-depth entries from the BED file
      awk '{for (i=$2+1; i<=$3; i++) print $1 "\t" i "\t" 0}' "${PSV_bed_file}" > "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_prefilled_depth_0.bed"

      #Run samtools depth for every bam file
      $samtools depth -@ 8 -a -b "${PSV_bed_file}" -o "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_actual_depth.bed" ${bam_file_path_region}
      echo "Samtools depth for ${sample_name}_${hap}_region completed"

      # Merge prefilled and actual depth file and remove duplicates, prioritizing actual depth values.
      #cat "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_actual_depth.bed" "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_prefilled_depth_0.bed" | sort -k1,1 -k2,2n | awk '!seen[$1,$2]++' > "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_depth.bed"

      cat "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_actual_depth.bed" "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_prefilled_depth_0.bed" \
    | sort -k1,1 -k2,2n \
    | awk '{
        key = $1 FS $2;
        if (!(key in seen)) {
            seen[key] = $0;
        } else if ($3 > 0) {
            seen[key] = $0;
        }
    } END {
        for (key in seen) {
            print seen[key];
        }
    }' > "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_depth.bed"

      # Clean old files
      rm "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_prefilled_depth_0.bed" "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_actual_depth.bed"

      # Merge depth files of all haplotypes
      cat "${input_dir}/select_bed_or_region/${sample_name}_${hap}_region_depth.bed" >> "${input_dir}/select_bed_or_region/${sample_name}_merged_region_depth.bed"
      echo "Depth files for ${sample_name} have been merged."

      fi
    done

    count_bed=$(awk '$3 < 3 {count++} END {print count+0}' "${input_dir}/select_bed_or_region/${sample_name}_merged_bed_depth.bed")
    echo "Number of PSVs with read depth below 3 in bed phasing of ${sample_name} is ${count_bed}."
    count_region=$(awk '$3 < 3 {count++} END {print count+0}' "${input_dir}/select_bed_or_region/${sample_name}_merged_region_depth.bed")
    echo "Number of PSVs with read depth below 3 in region phasing of ${sample_name} is ${count_region}."

    if [ "${count_bed}" -lt "${count_region}" ]; then
      echo "Bed is best phasing method for ${sample_name}."
      echo -e "${sample_name}\tbed" >> "${input_dir}/select_bed_or_region/bed_or_region_overview.txt"
    elif [ "${count_bed}" -gt "${count_region}" ]; then
      echo "Region is best phasing method for ${sample_name}."
      echo -e "${sample_name}\tregion" >> "${input_dir}/select_bed_or_region/bed_or_region_overview.txt"
    else
      echo "Bed and region phasing are equal, region phasing is chosen."
      echo -e "${sample_name}\tregion" >> "${input_dir}/select_bed_or_region/bed_or_region_overview.txt"
    fi

    echo "
    "

  fi
done

# Make a variable of the bed/region selection overview
file_overview="${input_dir}/select_bed_or_region/bed_or_region_overview.txt"

# Make a directory to move the old phasing files to
mkdir -p ${input_dir}/phasing_not_selected

for sample_name in $(awk '{print $1}'  ${file_overview})
do
  sample=$(awk '/'${sample_name}'/ {print $1}' ${file_overview})
  best_phasing=$(awk '/'${sample_name}'/ {print $2}' ${file_overview})
  echo "Sample is ${sample}."
  mkdir -p ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged
  mkdir -p ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged_split
  mkdir -p ${input_dir}/phasing_not_selected/${sample}/clair3
  mkdir -p ${input_dir}/phasing_not_selected/${sample}/sniffles2
  mkdir -p ${input_dir}/phasing_not_selected/${sample}/vcf
  echo "Best phasing is ${best_phasing}."
  if [ "$best_phasing" == "region" ]; then
    mv ${input_dir}/${sample}*/bam_files_haplotagged/*_bed.bam ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged
    mv ${input_dir}/${sample}*/bam_files_haplotagged/*_bed.bam.bai ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged
    mv ${input_dir}/${sample}*/bam_files_haplotagged_split/*bed*.bam ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged_split
    mv ${input_dir}/${sample}*/bam_files_haplotagged_split/*bed*.bam.bai ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged_split
    mv ${input_dir}/${sample}*/clair3/*bed*.vcf.gz ${input_dir}/phasing_not_selected/${sample}/clair3
    mv ${input_dir}/${sample}*/clair3/*bed*.vcf.gz.tbi ${input_dir}/phasing_not_selected/${sample}/clair3
    mv ${input_dir}/${sample}*/sniffles2/*bed*.vcf ${input_dir}/phasing_not_selected/${sample}/sniffles2
    mv ${input_dir}/${sample}*/vcf/*bed*.vcf.gz ${input_dir}/phasing_not_selected/${sample}/vcf
    mv ${input_dir}/${sample}*/vcf/*bed*.vcf.gz.tbi ${input_dir}/phasing_not_selected/${sample}/vcf
    echo "Moved bed phasing files for ${sample}."
  elif [ "$best_phasing" == "bed" ]; then
    mv ${input_dir}/${sample}*/bam_files_haplotagged/*_region.bam ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged
    mv ${input_dir}/${sample}*/bam_files_haplotagged/*_region.bam.bai ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged
    mv ${input_dir}/${sample}*/bam_files_haplotagged_split/*region*.bam ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged_split
    mv ${input_dir}/${sample}*/bam_files_haplotagged_split/*region*.bam.bai ${input_dir}/phasing_not_selected/${sample}/bam_files_haplotagged_split
    mv ${input_dir}/${sample}*/clair3/*region*.vcf.gz ${input_dir}/phasing_not_selected/${sample}/clair3
    mv ${input_dir}/${sample}*/clair3/*region*.vcf.gz.tbi ${input_dir}/phasing_not_selected/${sample}/clair3
    mv ${input_dir}/${sample}*/sniffles2/*region*.vcf ${input_dir}/phasing_not_selected/${sample}/sniffles2
    mv ${input_dir}/${sample}*/vcf/*region*.vcf.gz ${input_dir}/phasing_not_selected/${sample}/vcf
    mv ${input_dir}/${sample}*/vcf/*region*.vcf.gz.tbi ${input_dir}/phasing_not_selected/${sample}/vcf
    echo "Moved region phasing files for ${sample}."
  else
    echo "Overview file format is not valid."
    exit 1
  fi
done

echo "End of script."
