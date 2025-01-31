# -*- coding: utf-8 -*-

# Author: Demi Gommers

# Required modules:
import vcfpy
import pandas as pd
import sys
import os
import re


def extract_records_and_format(reader):
    records = []
    df_FORMAT_list = []

    for record in reader:
        try:
            line = [record.CHROM, record.POS, record.REF, record.ALT, record.QUAL, record.FILTER]
            # GT:GQ:DP:AD:AF
            for call in record.calls:
                gt = call.data.get('GT')
                gq = call.data.get('GQ')
                dp = call.data.get('DP')
                ad = call.data.get('AD')
                af = call.data.get('AF')
                df_FORMAT_list.append([gt, gq, dp, ad, af])

            records.append(line)
        except vcfpy.exceptions.CannotConvertValue:
            pass

    return records, df_FORMAT_list


def create_dataframe(records, df_FORMAT_list):
    records_df = pd.DataFrame(records, columns=['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER'])
    df_FORMAT = pd.DataFrame(df_FORMAT_list, columns=['GT', 'GQ', 'DP', 'AD', 'AF'])
    df_records_FORMAT = pd.concat([records_df, df_FORMAT], axis=1)

    return df_records_FORMAT


def save_dataframe_to_csv(df, output_file):
    df.to_csv(output_file, sep='\t', index=False)
    print("File is saved in the provided directory")


def extract_sample_and_haplotype(input_file_path):
    # Extract the file name from the path
    file_name = os.path.basename(input_file_path)

    # Regular expressions to extract sample name and haplotype
    sample_pattern = r'HG\d+|(SMA\d+)((?:_blood)?(?:_fib_P\d+)*(?:_P\d+)?)'
    haplotype_pattern = r'hap[0-9]+'

    sample_match = re.search(sample_pattern, file_name)
    haplotype_match = re.search(haplotype_pattern, file_name)

    sample_name = sample_match.group(0) if sample_match else None
    haplotype = haplotype_match.group(0) if haplotype_match else None

    return sample_name, haplotype


def extract_value(df, column='ALT'):
    # Converting value of ALT column into a string
    df['ALT'] = df[column].astype(str)

    # Extracting the ALT value of the variant into new column
    df['ALT_value'] = df[column].str.extract(r"value='(.*?)'")

    # Extracting the variant type into new column
    df['variant_type'] = df[column].str.extract(r"(SNV|DEL|INS)")

    return df


def add_sample_name_and_haplotype(df, input_file_path):
    # Extract sample name and haplotype from the input file path
    sample_name, haplotype = extract_sample_and_haplotype(input_file_path)
    df['sample_name'] = sample_name
    df['hap'] = haplotype

    # Print sample name and haplotype:
    print("Sample Name:", sample_name)
    print("Haplotype:", haplotype)


def main():
    vcf_file = sys.argv[1]
    output_file = sys.argv[2]
    reader = vcfpy.Reader.from_path(vcf_file)
    records, df_FORMAT_list = extract_records_and_format(reader)
    df_records_FORMAT = create_dataframe(records, df_FORMAT_list)
    add_sample_name_and_haplotype(df_records_FORMAT, vcf_file)
    df_records_FORMAT = extract_value(df_records_FORMAT)

    # Excluding the ALT column
    df_records_FORMAT = df_records_FORMAT.drop("ALT", axis=1)
    print(df_records_FORMAT)
    save_dataframe_to_csv(df_records_FORMAT, output_file)
    print("End of script")


if __name__ == "__main__":
    main()
