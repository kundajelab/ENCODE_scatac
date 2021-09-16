"""
Adapted from https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_lib_log_parser.py
"""

from collections import OrderedDict
import json
import os

def to_int(var):
    try:
        return int(var)
    except ValueError:
        return None


def to_float(var):
    try:
        return float(var)
    except ValueError:
        return None


def to_bool(var):
    return var.lower() in set(['true', 't', 'ok', 'yes', '1'])

def parse_frac_mito_qc(txt):
    result = OrderedDict()
    with open(txt, 'r') as fp:
        for line in fp.read().strip('\n').split('\n'):
            k, v = line.split('\t')
            if k.startswith('frac_'):
                result[k] = float(v)
            else:
                result[k] = int(v)
    return result


def parse_flagstat_qc(txt):
    result = OrderedDict()
    if not txt:
        return result
    total = ''
    total_qc_failed = ''
    duplicates = ''
    duplicates_qc_failed = ''
    mapped = ''
    mapped_qc_failed = ''
    mapped_pct = ''
    paired = ''
    paired_qc_failed = ''
    read1 = ''
    read1_qc_failed = ''
    read2 = ''
    read2_qc_failed = ''
    paired_properly = ''
    paired_properly_qc_failed = ''
    paired_properly_pct = ''
    with_itself = ''
    with_itself_qc_failed = ''
    singletons = ''
    singletons_qc_failed = ''
    singletons_pct = ''
    diff_chroms = ''
    diff_chroms_qc_failed = ''

    delimiter_pass_fail = ' + '
    with open(txt, 'r') as f:
        for line in f:
            if ' total ' in line:
                if ' in total ' in line:
                    tmp1 = line.split(' in total ')
                else:
                    tmp1 = line.split(' total ')
                line1 = tmp1[0]
                tmp1 = line1.split(delimiter_pass_fail)
                total = tmp1[0]
                total_qc_failed = tmp1[1]
            if ' duplicates' in line:
                tmp2 = line.split(' duplicates')
                line2 = tmp2[0]
                tmp2 = line2.split(delimiter_pass_fail)
                duplicates = tmp2[0]
                duplicates_qc_failed = tmp2[1]
            if ' mapped (' in line:
                tmp3 = line.split(' mapped (')
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(delimiter_pass_fail)
                mapped = tmp3_1[0]
                mapped_qc_failed = tmp3_1[1]
                line3_2 = tmp3[1]
                tmp3_2 = line3_2.split(':')
                mapped_pct = tmp3_2[0]  # .replace('%','')
            if ' paired in sequencing' in line:
                tmp2 = line.split(' paired in sequencing')
                line2 = tmp2[0]
                tmp2 = line2.split(delimiter_pass_fail)
                paired = tmp2[0]
                paired_qc_failed = tmp2[1]
            if ' read1' in line:
                tmp2 = line.split(' read1')
                line2 = tmp2[0]
                tmp2 = line2.split(delimiter_pass_fail)
                read1 = tmp2[0]
                read1_qc_failed = tmp2[1]
            if ' read2' in line:
                tmp2 = line.split(' read2')
                line2 = tmp2[0]
                tmp2 = line2.split(delimiter_pass_fail)
                read2 = tmp2[0]
                read2_qc_failed = tmp2[1]
            if ' properly paired (' in line:
                tmp3 = line.split(' properly paired (')
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(delimiter_pass_fail)
                paired_properly = tmp3_1[0]
                paired_properly_qc_failed = tmp3_1[1]
                line3_2 = tmp3[1]
                tmp3_2 = line3_2.split(':')
                paired_properly_pct = tmp3_2[0]  # .replace('%','')
            if ' with itself and mate mapped' in line:
                tmp3 = line.split(' with itself and mate mapped')
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(delimiter_pass_fail)
                with_itself = tmp3_1[0]
                with_itself_qc_failed = tmp3_1[1]
            if ' singletons (' in line:
                tmp3 = line.split(' singletons (')
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(delimiter_pass_fail)
                singletons = tmp3_1[0]
                singletons_qc_failed = tmp3_1[1]
                line3_2 = tmp3[1]
                tmp3_2 = line3_2.split(':')
                singletons_pct = tmp3_2[0]  # .replace('%','')
            if ' with mate mapped to a different chr' in line:
                tmp3 = line.split(' with mate mapped to a different chr')
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(delimiter_pass_fail)
                diff_chroms = tmp3_1[0]
                diff_chroms_qc_failed = tmp3_1[1]
    if total:
        result['total_reads'] = int(total)
    if total_qc_failed:
        result['total_reads_qc_failed'] = int(total_qc_failed)
    if duplicates:
        result['duplicate_reads'] = int(duplicates)
    if duplicates_qc_failed:
        result['duplicate_reads_qc_failed'] = int(duplicates_qc_failed)
    if mapped:
        result['mapped_reads'] = int(mapped)
    if mapped_qc_failed:
        result['mapped_reads_qc_failed'] = int(mapped_qc_failed)
    if mapped_pct:
        if 'nan' not in mapped_pct and 'N/A' not in mapped_pct \
                and 'NA' not in mapped_pct:
            if '%' in mapped_pct:
                mapped_pct = mapped_pct.replace('%', '')
                result['pct_mapped_reads'] = float(mapped_pct)
            else:
                result['pct_mapped_reads'] = 100.0 * float(mapped_pct)
        else:
            result['pct_mapped_reads'] = 0.0
    if paired:
        result['paired_reads'] = int(paired)
    if paired_qc_failed:
        result['paired_reads_qc_failed'] = int(paired_qc_failed)
    if read1:
        result['read1'] = int(read1)
    if read1_qc_failed:
        result['read1_qc_failed'] = int(read1_qc_failed)
    if read2:
        result['read2'] = int(read2)
    if read2_qc_failed:
        result['read2_qc_failed'] = int(read2_qc_failed)
    if paired_properly:
        result['properly_paired_reads'] = int(paired_properly)
    if paired_properly_qc_failed:
        result['properly_paired_reads_qc_failed'] = int(
            paired_properly_qc_failed)
    if paired_properly_pct:
        if 'nan' not in paired_properly_pct and \
                'N/A' not in paired_properly_pct \
                and 'NA' not in paired_properly_pct:
            if '%' in paired_properly_pct:
                paired_properly_pct = paired_properly_pct.replace('%', '')
                result['pct_properly_paired_reads'] = float(
                    paired_properly_pct)
            else:
                result['pct_properly_paired_reads'] = 100.0 * \
                    float(paired_properly_pct)
        else:
            result['pct_properly_paired_reads'] = 0.0
    if with_itself:
        result['with_itself'] = int(with_itself)
    if with_itself_qc_failed:
        result['with_itself_qc_failed'] = int(with_itself_qc_failed)
    if singletons:
        result['singletons'] = int(singletons)
    if singletons_qc_failed:
        result['singletons_qc_failed'] = int(singletons_qc_failed)
    if singletons_pct:
        if 'nan' not in singletons_pct and 'N/A' not in singletons_pct \
                and 'NA' not in singletons_pct:
            if '%' in singletons_pct:
                singletons_pct = singletons_pct.replace('%', '')
                result['pct_singletons'] = float(singletons_pct)
            else:
                result['pct_singletons'] = 100.0 * float(singletons_pct)
        else:
            result['pct_singletons'] = 0.0
    if diff_chroms:
        result['diff_chroms'] = int(diff_chroms)
    if diff_chroms_qc_failed:
        result['diff_chroms_qc_failed'] = int(diff_chroms_qc_failed)
    return result


def parse_dup_qc(txt):
    result = OrderedDict()
    if not txt:
        return result
    paired_reads = ''
    unpaired_reads = ''
    unmapped_reads = ''
    unpaired_dupes = ''
    paired_dupes = ''
    paired_opt_dupes = ''
    dupes_pct = ''

    picard_log_found = False
    # picard markdup
    with open(txt, 'r') as f:
        header = ''  # if 'UNPAIRED_READS_EXAMINED' in header
        content = ''
        for line in f:
            if header:
                content = line.replace(',', '.')
                picard_log_found = True
                break
            if 'UNPAIRED_READS_EXAMINED' in line:
                header = line
    if picard_log_found:
        header_items = header.split('\t')
        content_items = content.split('\t')
        m = dict(zip(header_items, content_items))
        unpaired_reads = m['UNPAIRED_READS_EXAMINED']
        paired_reads = m['READ_PAIRS_EXAMINED']
        unmapped_reads = m['UNMAPPED_READS']
        unpaired_dupes = m['UNPAIRED_READ_DUPLICATES']
        paired_dupes = m['READ_PAIR_DUPLICATES']
        paired_opt_dupes = m['READ_PAIR_OPTICAL_DUPLICATES']
        if 'PERCENT_DUPLICATION' in m:
            dupes_pct = m['PERCENT_DUPLICATION']
        else:
            dupes_pct = '0'
    else:
        # sambamba markdup
        with open(txt, 'r') as f:
            for line in f:
                if ' end pairs' in line:
                    tmp1 = line.strip().split(' ')
                    paired_reads = tmp1[1]
                if ' single ends ' in line:
                    tmp1 = line.strip().split(' ')
                    unpaired_reads = tmp1[1]
                    unmapped_reads = tmp1[6]
                if 'found ' in line:
                    tmp1 = line.strip().split(' ')
                    if paired_reads == '0':
                        unpaired_dupes = tmp1[1]  # SE
                        paired_dupes = 0
                    else:
                        unpaired_dupes = 0
                        paired_dupes = str(int(tmp1[1])/2)  # PE
                if paired_reads == '0':  # SE
                    dupes_pct = '{0:.2f}'.format(
                                float(unpaired_dupes)/float(unpaired_reads))
                elif paired_reads:
                    dupes_pct = '{0:.2f}'.format(
                                float(paired_dupes)/float(paired_reads))
    if unpaired_reads:
        result['unpaired_reads'] = int(unpaired_reads)
    if paired_reads:
        result['paired_reads'] = int(paired_reads)
    if unmapped_reads:
        result['unmapped_reads'] = int(unmapped_reads)
    if unpaired_dupes:
        result['unpaired_duplicate_reads'] = int(unpaired_dupes)
    if paired_dupes:
        result['paired_duplicate_reads'] = int(paired_dupes)
    if paired_opt_dupes:
        result['paired_optical_duplicate_reads'] = int(paired_opt_dupes)
    if dupes_pct:
        result['pct_duplicate_reads'] = float(dupes_pct)*100.0
    return result


def parse_lib_complexity_qc(txt):
    result = OrderedDict()
    if not txt:
        return result
    with open(txt, 'r') as f:
        for line in f:
            arr = line.strip().split('\t')
            break
    result['total_fragments'] = to_int(arr[0])
    result['distinct_fragments'] = to_int(arr[1])
    result['positions_with_one_read'] = to_int(arr[2])
    result['positions_with_one_read'] = to_int(arr[3])
    result['NRF'] = to_float(arr[4])
    result['PBC1'] = to_float(arr[5])
    result['PBC2'] = to_float(arr[6])
    return result


def parse_picard_est_lib_size_qc(txt):
    result = OrderedDict()
    if not txt:
        return result
    with open(txt, 'r') as f:
        val = f.readlines()[0].strip()
    result['picard_est_lib_size'] = float(val)
    return result


def build_quality_metric_header(sample, config, data_path, out_path):
    lab = config["dcc_lab"]
    data_alias = f"{lab}:{sample}${os.path.basename(data_path)}"
    alias = f"{lab}:{sample}${os.path.basename(out_path)}"
    h = OrderedDict({
        "lab": lab,
        "award": config["dcc_award"],
        "quality_metric_of": data_alias,
        "aliases": [alias],
    })
    return h


def write_json(data, out_path):
    with open(out_path, "w") as f:
        json.dump(data, f, indent=4)


try:
    out_group = snakemake.params['out_group']
    sample = snakemake.params['sample']
    data_path = snakemake.input['data_file']
    config = snakemake.config

    if out_group == "fastqs":
        pass

    elif out_group == "mapping":
        alignment_stats_out = snakemake.output['alignment_stats']
        samstats_raw = snakemake.input['samstats_raw']

        a = parse_flagstat_qc(samstats_raw)
        h = build_quality_metric_header(sample, config, data_path)
        alignment_stats = h | a

        write_json(alignment_stats, alignment_stats_out)

    elif out_group == "filtering":
        alignment_stats_out = snakemake.output['alignment_stats']
        lib_comp_stats_out = snakemake.output['lib_comp_stats']
        samstats_filtered = snakemake.input['samstats_filtered']
        picard_markdup = snakemake.input['picard_markdup']
        pbc_stats = snakemake.input['pbc_stats']
        frac_mito = snakemake.input['frac_mito']

        s = parse_flagstat_qc(samstats_filtered)
        p = parse_picard_est_lib_size_qc(picard_markdup)
        l = parse_lib_complexity_qc(pbc_stats)
        m = parse_frac_mito_qc(frac_mito)
        h = build_quality_metric_header(sample, config, data_path)
        alignment_stats = h | s | m
        lib_comp_stats = h | p | l

        write_json(alignment_stats, alignment_stats_out)
        write_json(lib_comp_stats, lib_comp_stats_out)

    elif out_group == "fragments":
        pass

    elif out_group == "archr":
        pass


except NameError:
    pass