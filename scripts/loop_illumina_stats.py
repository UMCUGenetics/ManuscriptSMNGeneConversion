#! /usr/bin/env python3

import argparse
import vcf
import glob
import os


def DetermineOverlap_ONT(input_folder, roi, min_gq, min_af, min_dp, position_list):
    vcf_files = glob.glob("{}/*.vcf.gz".format(input_folder), recursive=True)
    chrom_roi = roi.split(":")[0]
    start_roi, stop_roi = roi.split(":")[1].split("-")
    variant_list = []
    for vcf_file in vcf_files:
        vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
        sampleid = vcf_reader.samples[0]  # Assume single sample VCF here!
        for record in vcf_reader:
            if not record.is_snp:
                continue
            chrom = record.CHROM
            start = int(record.POS) - 1
            stop = record.POS
            ref = record.REF
            if len(record.ALT) > 1:
                alt = "_".join(record.ALT)
            else:
                alt = record.ALT[0]
            pass_filter = "_".join(record.FILTER)
            if chrom == chrom_roi and start >= int(start_roi) and stop <= int(stop_roi):
                gq = float(record.genotype(sampleid)['GQ'])
                dp = float(record.genotype(sampleid)['DP'])
                af = float(record.genotype(sampleid)['AF'])
                qc_pass = "failed"
                if gq >= min_gq and af >= min_af and dp >= min_dp and record.is_snp:
                    qc_pass = "passed"
                variant_list.append(
                    f"{sampleid}\t{chrom}\t{start}\t{stop}\t{ref}\t{alt}\t{gq}\t{dp}\t{af}\t{qc_pass}\t{pass_filter}"
                )
                position = f"{chrom}_{start}_{stop}"
                if position not in position_list:
                    position_list.append(position)

    return variant_list, position_list


def DetermineOverlap_ILL(input_folder, roi, sample_dict, min_dp, min_af):
    vcf_files = glob.glob("{}/*.vcf.gz".format(input_folder), recursive=True)
    chrom_roi = roi.split(":")[0]
    start_roi, stop_roi = roi.split(":")[1].split("-")
    variant_list = []
    position_list = []
    for vcf_file in vcf_files:
        vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
        sampleid = vcf_reader.samples[0]  # Assume single sample VCF here!
        for record in vcf_reader:
            chrom = record.CHROM
            start = int(record.POS) - 1
            stop = record.POS
            ref = record.REF
            if record.genotype(sampleid).is_variant and record.is_snp:
                if chrom == chrom_roi and start >= int(start_roi) and stop <= int(stop_roi) and len(record.ALT) == 1:
                    if record.ALT[0] == "*":
                        continue
                    alt = record.ALT[0]
                    gq = float(record.genotype(sampleid)['GQ'])
                    dp = float(record.genotype(sampleid)['DP'])
                    gt = str(record.genotype(sampleid)['GT'])
                    ad = record.genotype(sampleid)['AD']
                    af = (ad[1]/dp)*100
                    sample_id_ont = sample_dict[sampleid]
                    if dp >= min_dp and af >= min_af:
                        variant_list.append(
                            f"{sample_id_ont}\t{sampleid}\t{chrom}\t{start}\t{stop}\t{ref}\t{alt}\t{gq}\t{dp}\t{gt}\t{af:.2f}\tpassed"
                        )
                        position = f"{chrom}_{start}_{stop}"
                        if position not in position_list:
                            position_list.append(position)

    return variant_list, position_list


def sample_translation(sample_file):
    sample_dict = {}
    for line in sample_file:
        splitline = line.split()
        if splitline[0] not in sample_dict:
            sample_dict[splitline[0]] = splitline[1]
    return sample_dict


def DetermineOverlap(position_dic, roi, sample_id, analysis, qc):
    roi_chrom = roi.split(":")[0]
    roi_start, roi_stop = roi.split(":")[1].split("-")

    positions_called_ill_ont_passed = 0
    positions_called_ill_ont_failed = 0
    positions_called_ill_only = 0
    positions_called_ont_onlyP = 0
    positions_called_ont_onlyF = 0
    positions_no_call = 0
    positions_total = 0

    with open(f"{analysis}/{sample_id}_{qc}_stats.bed", 'w') as write_file:
        for pos in position_dic:
            chrom, start, stop = pos.split("_")
            if chrom == roi_chrom and start <= roi_stop and stop >= roi_start:
                illumina_filter = None

                if sample_id in position_dic[pos]["illumina"]:  # Assumes all Illumina input in hardfiltered passed
                    illumina_filter = "passed"

                # check if at least one passed in ONT haplotypes
                ont_filter = None
                if sample_id in position_dic[pos]["ont"]:
                    ont_filter = "failed"
                    ont_geno = position_dic[pos]["ont"][sample_id][0][4]
                    for haplotype in position_dic[pos]["ont"][sample_id]:
                        if haplotype[8] == "passed":  # if one haplo is passed, assume whole position passed for ont
                            ont_filter = "passed"
                        if haplotype[4] is not ont_geno:
                            write_file.write(
                               f"{chrom}\t{start}\t{stop}\tWARNING: sample {sample_id} has different alleles within haplotypes {pos} has ALT {haplotype[4]} and ALT {ont_geno} continued with first allele\n"
                            )

                ill_line = "n/a"
                ont_line = "n/a"

                if sample_id in position_dic[pos]["illumina"]:
                    ill_line = "_".join(position_dic[pos]["illumina"][sample_id])
                if sample_id in position_dic[pos]["ont"]:
                    ont_line = position_dic[pos]["ont"][sample_id]
                    write_line = []
                    for haplo in ont_line:
                        write_line.append("_".join(haplo))
                    ont_line = "/".join(write_line)

                if illumina_filter == "passed" and ont_filter == "passed":
                    positions_called_ill_ont_passed += 1
                    write_file.write(f"{chrom}\t{start}\t{stop}\tdetected_ill_ont\t{sample_id}\t{ill_line}\t{ont_line}\n")
                elif illumina_filter == "passed" and ont_filter == "failed":  # Although genotype in ONT,all are not passed and considered fn
                    positions_called_ill_ont_failed += 1
                    positions_called_ill_only += 1
                    write_file.write(f"{chrom}\t{start}\t{stop}\tdetected_ill_only_ont_failed\t{sample_id}\t{ill_line}\t{ont_line}\n")
                elif illumina_filter == "passed" and not ont_filter:
                    positions_called_ill_only += 1
                    write_file.write(f"{chrom}\t{start}\t{stop}\tdetected_ill_only\t{sample_id}\t{ill_line}\t{ont_line}\n")
                elif not illumina_filter and ont_filter == "failed":
                    positions_called_ont_onlyF += 1
                    write_file.write(f"{chrom}\t{start}\t{stop}\tno_call_ill_ont_failed\t{sample_id}\t{ill_line}\t{ont_line}\n")
                elif not illumina_filter and ont_filter == "passed":
                    positions_called_ont_onlyP += 1
                    write_file.write(f"{chrom}\t{start}\t{stop}\tdetected_ont_only\t{sample_id}\t{ill_line}\t{ont_line}\n")
                elif not illumina_filter and not ont_filter:
                    write_file.write(f"{chrom}\t{start}\t{stop}\tno_call_ill_ont\t{sample_id}\t{ill_line}\t{ont_line}\n")
                    positions_no_call += 1
                positions_total += 1

    try:
        sensitivity = (positions_called_ill_ont_passed/(positions_called_ill_ont_passed + positions_called_ill_only)) * 100
        sensitivity = f"{sensitivity:.2f}"
    except:
        sensitivity = "n/a"
    try:
        sensitivity_failed = ((positions_called_ill_ont_passed + positions_called_ill_ont_failed) / (positions_called_ill_ont_passed + positions_called_ill_only)) * 100
        sensitivity_failed = f"{sensitivity_failed:.2f}"
    except:
        sensitivity_failed = "n/a"

    try:
        precision = (positions_called_ill_ont_passed/(positions_called_ill_ont_passed + positions_called_ont_onlyP)) * 100
        precision = f"{precision:.2f}"
    except:
        precision = "n/a"

    try:
        precision_failed = (positions_called_ill_ont_passed/(positions_called_ill_ont_passed + positions_called_ont_onlyP + positions_called_ont_onlyF)) * 100
        precision_failed = f"{precision_failed:.2f}"
    except:
        precision_failed = "n/a"

    print(f"{sample_id}\t{positions_total}\t{positions_called_ill_only}\t{positions_called_ont_onlyP}\t{positions_called_ill_ont_passed}\t{positions_no_call}\t{sensitivity}\t{precision}\t{positions_called_ill_ont_failed}\t{sensitivity_failed}\t{positions_called_ont_onlyF}\t{precision_failed}")
    return positions_total, positions_called_ill_only, positions_called_ont_onlyP, positions_called_ill_ont_passed, positions_no_call, positions_called_ill_ont_failed, positions_called_ont_onlyF


def AddToDicIllumina(data, position_dic, method, sample_id_ill, sample_id_ont):
    sample_dic = {}
    for line in data:
        if sample_id_ill not in line:
            continue
        sample_id = sample_id_ont
        splitline = line.split()
        if sample_id_ill not in sample_dic:
            sample_dic[sample_id] = {"illumina": True, "ont": False}
        position = f"{splitline[2]}_{splitline[3]}_{splitline[4]}"
        if position in position_dic:
            if splitline[11] == "passed":  # include only if passed
                if sample_id not in position_dic[position]["illumina"]:
                    position_dic[position]["illumina"][sample_id] = []
                position_dic[position]["illumina"][sample_id] = splitline[2:]
    return position_dic, sample_dic


def AddToDicONT(data, position_dic, method, sample_dic, sample_id):
    for line in data:
        if sample_id not in line:
            continue
        splitline = line.split()
        if sample_id not in sample_dic:
            sample_dic[sample_id] = {"illumina": False, "ont": True}
        elif sample_id in sample_dic:
            sample_dic[sample_id]["ont"] = True

        position = f"{splitline[1]}_{splitline[2]}_{splitline[3]}"
        if position in position_dic:
            if sample_id not in position_dic[position]["ont"]:
                position_dic[position]["ont"][sample_id] = []
            # Added all haplotype variants if multiple are detected in ont
            position_dic[position]["ont"][sample_id].append(splitline[1:])
    return position_dic, sample_dic


def BEDdict(position_list):
    position_dic = {}
    for line in position_list:
        position = line.rstrip()
        if position not in position_dic:
            position_dic[position] = {"illumina": {}, "ont": {}}
    return position_dic


def TranslationDict(translation_file):
    translation = {}
    for line in translation_file:
        splitline = line.split()
        if splitline[0] not in translation:
            translation[splitline[0]] = splitline[1]
    return translation


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_folder_ont', help='path to folder containing all VCF files from ONT analysis')
    parser.add_argument('input_folder_ill', help='path to folder containing all VCF files from Illumina analysis')
    parser.add_argument('roi', help='region of interest chr:start-stop')
    parser.add_argument('sample_file_illumina', type=argparse.FileType('r'), help='File file sample_ids (for renaming)')
    parser.add_argument('translation_file', type=argparse.FileType('r'), help='TSV file with illumina, ont, ploidy')
    parser.add_argument('analysisID', help='ID of analysis (output folder)')
    parser.add_argument('--min_gq_ont', default=10, type=float, help='minimum genotype quality variant [default 10]')
    parser.add_argument('--min_dp_ont', default=4, type=float, help='minimum depth variant [default 4]')
    parser.add_argument(
        '--min_af_ont',
        default=[0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00],
        type=list,
        help='minimum allele frequency variant [default list [0, 0.5, 0.85]]'
    )
    parser.add_argument('--min_dp_ill', default=10, type=float, help='minimum depth variant [default 4]')
    parser.add_argument('--min_af_ill', default=15, type=float, help='minimum allele frequency variant [default 0.85]')
    args = parser.parse_args()

    translation = TranslationDict(args.translation_file)
    sample_dict_illumina = sample_translation(args.sample_file_illumina)
    ill_list, position_list_ill = DetermineOverlap_ILL(args.input_folder_ill, args.roi, sample_dict_illumina, args.min_dp_ill, args.min_af_ill)

    os.system(f"rm -rf {args.analysisID}")
    os.system(f"mkdir {args.analysisID}")

    for min_af in args.min_af_ont:
        print("analysis\tsample_id\ttotal_positions\tcalled_illumina_only\tcalled_ont_onlyP\tcalled_ont_illumina\tnot_called_ont_illumina\tsensitivity\tprecision\tcalled_ont_failed_illumina\tsensitivity_ont_failed\tcalled_ont_onlyF\tprecision_ont_failed")
        ont_list, position_list = DetermineOverlap_ONT(args.input_folder_ont, args.roi, args.min_gq_ont, min_af, args.min_dp_ont, position_list_ill)

        final_stats = []
        for sample in translation:
            position_dic = BEDdict(position_list)
            position_dic, sample_dic = AddToDicIllumina(ill_list, position_dic, "illumina", sample, translation[sample])
            position_dic, sample_dic = AddToDicONT(ont_list, position_dic, "ont", sample_dic, translation[sample])
            stats = DetermineOverlap(position_dic, args.roi, translation[sample], args.analysisID, min_af)
            final_stats.append(stats)

        total_positions = 0
        called_ill_only = 0
        called_ont_onlyP = 0
        called_ont_onlyF = 0
        called_ill_ont = 0
        not_called_ont_illumina = 0
        called_ill_ont_failed = 0
        for item in final_stats:
            total_positions += item[0]
            called_ill_only += item[1]
            called_ont_onlyP += item[2]
            called_ill_ont += item[3]
            not_called_ont_illumina += item[4]
            called_ill_ont_failed += item[5]
            called_ont_onlyF += item[6]
        try:
            sensitivity = (called_ill_ont/(called_ill_ont + called_ill_only)) * 100
            sensitivity = f"{sensitivity:.2f}"
        except:
            sensitivity = "n/a"
        try:
            sensitivity_failed = ((called_ill_ont + called_ill_ont_failed) / (called_ill_ont + called_ill_only)) * 100
            sensitivity_failed = f"{sensitivity_failed:.2f}"
        except:
            sensitivity = "n/a"
        try:
            precision = (called_ill_ont/(called_ill_ont + called_ont_onlyP)) * 100
            precision = f"{precision:.2f}"
        except:
            precision = "n/a"

        try:
            precision_failed = (called_ill_ont/(called_ill_ont + called_ont_onlyP + called_ont_onlyF)) * 100
            precision_failed = f"{precision_failed:.2f}"
        except:
            precision_failed = "n/a"
        print(f"Total\t{min_af}\t{total_positions}\t{called_ill_only}\t{called_ont_onlyP}\t{called_ill_ont}\t{not_called_ont_illumina}\t{sensitivity}\t{precision}\t{called_ill_ont_failed}\t{sensitivity_failed}\t{called_ont_onlyF}\t{precision_failed}")
