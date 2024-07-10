#!/bin/bash
#SBATCH --job-name=vcf_parse_merge_depth
#SBATCH --time=01:00:00
#SBATCH --gres=tmpspace:50G
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --ntasks=6
#SBATCH --output=vcf_parse_merge_depth_%j.out

# Settings
samtools=/path/to/samtools

# Stop execution of script whenever a non-zero exit status is returned
set -e

# Help function so other users can see what input is required
Help()
{
  # Display Help
  echo "Executing for multiple vcf files"
  echo
  echo
  echo "Syntax: Input requirements [-d:h:o:i:f]"
  echo "options:"
  echo "h     Print help"
  echo "o     Output directory"
  echo "i     Input directory, requires both SMA and 1000G subdirectory"
}

### Files required ###
while getopts h:o:i: flag
do
  case "${flag}" in
    o) output_dir=${OPTARG};;
    i) input_dir=${OPTARG};;
    h) # display Help
        Help
        exit 1;;
    *) Help
        echo "Command not found"
        exit 1;;
  esac
done

## Repository wkdir
repo_folder="$(dirname "$(readlink -f "$0")")"

## Load virtual venv
source ${repo_folder}/venv/bin/activate

#Part 1 of script
#parsing VCF Files

# Directory containing subdirectories with VCF files
mkdir -p "${output_dir}/SMA/vcfs_parsed/"
mkdir -p "${output_dir}/1000G/vcfs_parsed/"

# Iterate over VCF files in the directory
# SMA
for vcf_file in ${input_dir}/SMA/{SMA,HG}*/clair3/*.vcf.gz; do
  if [ -f "${vcf_file}" ]; then
    echo "Processing file: ${vcf_file}"
    file_name=$(basename "${vcf_file}" .vcf.gz)
    # Execute vcf parser python script
    ${repo_folder}/vcf_parser.py "${vcf_file}" "${output_dir}/SMA/vcfs_parsed/${file_name}.tsv"
  fi
done

echo "Concatenating all SMA .tsv files into one file. Contains one header"
awk 'FNR == 1 && NR!=1 { while (/^CHROM/) getline; }
	1 {print}
' ${output_dir}/SMA/vcfs_parsed/*filter.tsv > "${output_dir}/SMA/vcfs_parsed/all_variants_all_haps_sample.tsv"

# Iterate over VCF files in the directory
# 1000G
for vcf_file in ${input_dir}/1000G/{SMA,HG}*/clair3/*.vcf.gz; do
  if [ -f "${vcf_file}" ]; then
    echo "Processing file: ${vcf_file}"
    file_name=$(basename "${vcf_file}" .vcf.gz)
    # Execute vcf parser python script
    ${repo_folder}/vcf_parser.py "${vcf_file}" "${output_dir}/1000G/vcfs_parsed/${file_name}.tsv"
  fi
done

echo "Concatenating all 1000G .tsv files into one file. Contains one header"
awk 'FNR == 1 && NR!=1 { while (/^CHROM/) getline; }
	1 {print}
' ${output_dir}/1000G/vcfs_parsed/*filter.tsv > "${output_dir}/1000G/vcfs_parsed/all_variants_all_haps_sample.tsv"


#Part 2 of script
# Creating input bed file for samtools depth based on provided vcf parsed files
# Preparing SMA bed file
mkdir -p "${output_dir}/SMA/bed_files_variants/"
mkdir -p "${output_dir}/1000G/bed_files_variants/"

input_SMA=${output_dir}/SMA/vcfs_parsed/all_variants_all_haps_sample.tsv
input_1000G=${output_dir}/1000G/vcfs_parsed/all_variants_all_haps_sample.tsv

echo "Cutting the first two columns"
cat ${input_SMA} | cut -f1,2 > ${output_dir}/SMA/bed_files_variants/all_variants_SMA.bed

echo "-------SMA variant file-------"
echo "Creating new column to print out BED positions"
awk '{print $0, $2-1}' ${output_dir}/SMA/bed_files_variants/all_variants_SMA.bed > tmp_file && mv tmp_file ${output_dir}/SMA/bed_files_variants/all_variants_SMA_pos.bed
awk 'BEGIN{OFS = "\t"} {print $1, $3, $2}' ${output_dir}/SMA/bed_files_variants/all_variants_SMA_pos.bed > ${output_dir}/SMA/bed_files_variants/all_variants_SMA_pos_complete.bed

echo "Checking existence of positions, only include unique positions"
awk '!seen[$2,$3]++' ${output_dir}/SMA/bed_files_variants/all_variants_SMA_pos_complete.bed > ${output_dir}/SMA/bed_files_variants/all_variants_SMA_uniq.bed


# Preparing 1000G bed file
cat ${input_1000G} | cut -f1,2 > ${output_dir}/1000G/bed_files_variants/all_variants_1000G.bed

echo "-------1000G variant file-------"
echo "Creating new column to print out BED positions"

awk '{print $0, $2-1}' ${output_dir}/1000G/bed_files_variants/all_variants_1000G.bed > tmp_file && mv tmp_file ${output_dir}/1000G/bed_files_variants/all_variants_1000G_pos.bed
awk 'BEGIN{OFS = "\t"} {print $1, $3, $2}' ${output_dir}/1000G/bed_files_variants/all_variants_1000G_pos.bed > ${output_dir}/1000G/bed_files_variants/all_variants_1000G_pos_complete.bed

echo "Checking existence of positions, only include unique positions"
awk '!seen[$2,$3]++' ${output_dir}/1000G/bed_files_variants/all_variants_1000G_pos_complete.bed > ${output_dir}/1000G/bed_files_variants/all_variants_1000G_uniq.bed



# Concatenating both bed files together:
echo "excluding headers of both files"
tail -n +2 "${output_dir}/1000G/bed_files_variants/all_variants_1000G_uniq.bed" > "${output_dir}/1000G/bed_files_variants/all_variants_1000G_uniq.tmp" && mv "${output_dir}/1000G/bed_files_variants/all_variants_1000G_uniq.tmp" "${output_dir}/1000G/bed_files_variants/all_variants_1000G_uniq_noheader.bed"

tail -n +2 "${output_dir}/SMA/bed_files_variants/all_variants_SMA_uniq.bed" > "${output_dir}/SMA/bed_files_variants/all_variants_SMA_uniq.tmp" && mv "${output_dir}/SMA/bed_files_variants/all_variants_SMA_uniq.tmp" "${output_dir}/SMA/bed_files_variants/all_variants_SMA_uniq_noheader.bed"

echo "concatenating both bed files together and excluding duplicate positions"
cat "${output_dir}/SMA/bed_files_variants/all_variants_SMA_uniq_noheader.bed" "${output_dir}/1000G/bed_files_variants/all_variants_1000G_uniq_noheader.bed" | awk '!seen[$2,$3]++' > ${output_dir}/all_variants_positions_uniq.bed

#part 3
#make depth file for all variant positions

#define bed file
bed_file=${output_dir}/all_variants_positions_uniq.bed

mkdir -p "${output_dir}/SMA/depth_variant_positions"
mkdir -p "${output_dir}/1000G/depth_variant_positions"


# Iterate over directories
# SMA
for dir in ${input_dir}/SMA/{SMA,HG}*/bam_files_haplotagged_split; do
  if [ -d "${dir}" ]; then
    echo "Processing files in directory: ${dir}"

    # Iterate over bam files in the directory
    for bam_file_path in ${dir}/*.bam; do
      if [ -f "${bam_file_path}" ]; then
      echo "Processing file: ${bam_file_path}"
      bam_file=$(basename $bam_file_path)
      sample_name=$(echo $bam_file | grep -Eo "HG[0-9]+|(SMA[0-9]+)((_blood)?(_fib_P[0-9]+)*(_P[0-9]+)?)")
      echo "${sample_name}"
      hap=$(echo $bam_file | grep -Eo "hap[0-9]+")
      echo "${hap}"

      #Run samtools depth for every bam file
      $samtools depth -@ 2 -H -b "${bed_file}" -o "${output_dir}/SMA/depth_variant_positions/${sample_name}_${hap}_depth.bed" ${bam_file_path}
      echo "Samtools depth for ${sample_name}_${hap} completed"

      #add correct header per bed file
      awk -v sample_name="${sample_name}" -v hap="${hap}" 'BEGIN{OFS="\t"} NR==1 {print "CHROM", $2, "depth", "sample_name", "hap"; next} {print $0, sample_name, hap}' "${output_dir}/SMA/depth_variant_positions/${sample_name}_${hap}_depth.bed" > "${output_dir}/SMA/depth_variant_positions/${sample_name}_${hap}_depth_header.bed"
      fi
    done
  fi
done

#once all files have proper headers you can do this:
#merge all bed files, but take the header of only the first file
awk 'FNR == 1 && NR!=1 { while (/^CHROM/) getline; }
1 {print}
' ${output_dir}/SMA/depth_variant_positions/*header.bed > ${output_dir}/SMA/depth_variant_positions/depth_all_variant_positions_all_haps.bed


# Iterate over directories
# 1000G
for dir in ${input_dir}/1000G/{SMA,HG}*/bam_files_haplotagged_split; do
  if [ -d "${dir}" ]; then
    echo "Processing files in directory: ${dir}"

    # Iterate over bam files in the directory
    for bam_file_path in ${dir}/*.bam; do
      if [ -f "${bam_file_path}" ]; then
      echo "Processing file: ${bam_file_path}"
      bam_file=$(basename $bam_file_path)
      sample_name=$(echo $bam_file | grep -Eo "HG[0-9]+|(SMA[0-9]+)((_blood)?(_fib_P[0-9]+)*(_P[0-9]+)?)")
      echo "${sample_name}"
      hap=$(echo $bam_file | grep -Eo "hap[0-9]+")
      echo "${hap}"

      #Run samtools depth for every bam file
      $samtools depth -@ 2 -H -b "${bed_file}" -o "${output_dir}/1000G/depth_variant_positions/${sample_name}_${hap}_depth.bed" ${bam_file_path}
      echo "Samtool depth for ${sample_name}_${hap} completed"

      #add correct header per bed file
      awk -v sample_name="${sample_name}" -v hap="${hap}" 'BEGIN{OFS="\t"} NR==1 {print "CHROM", $2, "depth", "sample_name", "hap"; next} {print $0, sample_name, hap}' "${output_dir}/1000G/depth_variant_positions/${sample_name}_${hap}_depth.bed" > "${output_dir}/1000G/depth_variant_positions/${sample_name}_${hap}_depth_header.bed"
      fi
    done
  fi
done

#once all files have proper headers you can do this:
#merge all bed files, but take the header of only the first file
awk 'FNR == 1 && NR!=1 { while (/^CHROM/) getline; }
1 {print}
' ${output_dir}/1000G/depth_variant_positions/*header.bed > ${output_dir}/1000G/depth_variant_positions/depth_all_variant_positions_all_haps.bed


#part 4: merge depth column with merged vcf tsv
#SMA
${repo_folder}/merging_variant_depth_files.py "${output_dir}/SMA/vcfs_parsed/all_variants_all_haps_sample.tsv" "${output_dir}/SMA/depth_variant_positions/depth_all_variant_positions_all_haps.bed" "${output_dir}/SMA/vcf_depth_merged_all_haps.tsv"

#1000G
${repo_folder}/merging_variant_depth_files.py "${output_dir}/1000G/vcfs_parsed/all_variants_all_haps_sample.tsv" "${output_dir}/1000G/depth_variant_positions/depth_all_variant_positions_all_haps.bed" "${output_dir}/1000G/vcf_depth_merged_all_haps.tsv"


echo "End of script"
