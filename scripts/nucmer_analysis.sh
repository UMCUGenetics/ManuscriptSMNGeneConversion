#!/bin/bash
#SBATCH --job-name=nucmer_analysis
#SBATCH --time=00:10:00
#SBATCH --gres=tmpspace:1G
#SBATCH --mem=2G
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=10
#SBATCH --output=~/nucmer_analysis_%j.out


# Stop execution of script whenever a non-zero exit status is returned
set -e

# Help function so other users can see what input is required
Help()
{
  # Display Help
  echo "This script aligns two fasta files with nucmer and generates a filtered delta file and a coordinate file."
  echo
  echo
  echo "Syntax: Input requirements [-h:r:q:o]"
  echo "options:"
  echo "h     Print help"
  echo "r     Reference query (.fasta or .fa)"
  echo "q     Query (.fasta or .fa)"
  echo "o     Output name"
}

### Files required ###
while getopts h:o:r:q: flag
do
  case "${flag}" in
    o) output_name=${OPTARG};;
    r) reference_query=${OPTARG};;
    q) query=${OPTARG};;
    h) # display Help
        Help
        exit 1;;
    *) Help
        echo "Command not found"
        exit 1;;
  esac
done

# Load tools
nucmer="~/mummer.4.0.0/bin/./nucmer"
deltafilter="~/mummer.4.0.0/bin/./delta-filter"
showcoords="~/mummer.4.0.0/bin/./show-coords"

# Rename contigs based on file name
# Get the basename of the file (without extension
basename_reference_query=$(basename "$reference_query" .fasta)
basename_query=$(basename "$query" .fasta)

echo "Running nucmer"
# Run nucmer
${nucmer} -t 10 ${reference_query} ${query} -p ${output_name}

# Filter delta file
${deltafilter} -i 95 -l 10000 ${output_name}.delta > ${output_name}_95i_10k.delta

# Show coordinates from delta file
${showcoords} -clrTdH ${output_name}_95i_10k.delta > ${output_name}_95i_10k.coords.txt

# Convert coordinate file to bed file by selecting column 1, 2, and 8
${showcoords} -r ${output_name}_95i_10k.delta -T -H | cut -f1,2,8 | awk '{if($1 < $2){print $3,$1-1,$2} else {print $3,$2-1,$1}}' | tr " " "\t" > ${output_name}_95i_10k.delta.coords.bed

echo "Done"
