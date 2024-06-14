#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Demi Gommers
import pandas as pd
import sys


def load_file(file, header):
    print("Loading file as dataframe")
    df = pd.read_csv(file, sep='\t', header=header)
    return df


def save_dataframe_to_csv(df, output_file):
    df.to_csv(output_file, sep='\t', index=False)
    print("File is saved in the provided directory")


def merging_dataframes_on_columns(df1, df2, columns=["POS", "sample_name", "hap", "CHROM"], how="left"):
    print("Merged dataframes on provided columns")
    merged_df = pd.merge(df1, df2, on=columns, how=how)
    print("NAs are now empty")
    merged_df.fillna("", inplace=True)
    return merged_df


def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <variant_file> <depth_file> <output_file>")
        sys.exit(1)

    variant_file = sys.argv[1]
    depth_file = sys.argv[2]
    output_file = sys.argv[3]

    # Loading variant file into df
    variant_df = load_file(variant_file, header=0)

    # Loading depth file into df
    depth_df = load_file(depth_file, header=0)

    # Merging the depth df and the variants df
    merged_df = merging_dataframes_on_columns(depth_df, variant_df)
    save_dataframe_to_csv(merged_df, output_file)
    print("End of script")


if __name__ == "__main__":
    main()
