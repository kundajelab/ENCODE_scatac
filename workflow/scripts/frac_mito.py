import os
import re
import subprocess
import logging
import signal
import time
import sys
from collections import OrderedDict

"""
Derived from: https://github.com/ENCODE-DCC/atac-seq-pipeline
"""

mito_bam_path = snakemake.input['bam_mito']
no_mito_bam_path = snakemake.input['bam_no_mito']

mito_tmp_path = snakemake.output['tmp_sort_mito']
no_mito_tmp_path = snakemake.output['tmp_sort_no_mito']

mito_samstats_qc_path = snakemake.output['qc_samstats_mito']
no_mito_samstats_qc_path = snakemake.output['qc_samstats_no_mito']
frac_mito_qc_path = snakemake.output['qc_frac_mito']

threads = snakemake.threads
log_path, = snakemake.log

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    filename=log_path)
log = logging.getLogger(__name__)

def get_ticks():
    """Returns ticks.
        - Python3: Use time.perf_counter().
        - Python2: Use time.time().
    """
    return getattr(time, 'perf_counter', getattr(time, 'time'))()

def run_shell_cmd(cmd):
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],  # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid)  # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    t0 = get_ticks()
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    t1 = get_ticks()
    err_str = (
        'PID={pid}, PGID={pgid}, RC={rc}, DURATION_SEC={dur:.1f}\n'
        'STDERR={stde}\nSTDOUT={stdo}'
    ).format(
        pid=pid, pgid=pgid, rc=rc, dur=t1 - t0, stde=stderr.strip(), stdo=stdout.strip()
    )
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')

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

def frac_mito(non_mito_samstat, mito_samstat, qc_path):
    non_mito_samstat_dict = parse_flagstat_qc(non_mito_samstat)
    mito_samstat_dict = parse_flagstat_qc(mito_samstat)

    if 'mapped' in non_mito_samstat_dict:
        # backward compatibility (old key name was 'total')
        key_mapped = 'mapped'
    elif 'mapped_reads' in non_mito_samstat_dict:
        key_mapped = 'mapped_reads'
    Rn = non_mito_samstat_dict[key_mapped]

    if 'mapped' in mito_samstat_dict:
        # backward compatibility (old key name was 'total')
        key_mapped = 'mapped'
    elif 'mapped_reads' in mito_samstat_dict:
        key_mapped = 'mapped_reads'
    Rm = mito_samstat_dict[key_mapped]

    frac = float(Rm)/float(Rn + Rm)
    with open(qc_path, 'w') as fp:
        fp.write('non_mito_reads\t{}\n'.format(Rn))
        fp.write('mito_reads\t{}\n'.format(Rm))
        fp.write('frac_mito_reads\t{}\n'.format(frac))

def samstat(bam_path, qc_path, threads):
    run_shell_cmd(
        'samtools sort -n {bam} -@ {threads} -O sam | '
        'SAMstats --sorted_sam_file - --outf {samstat_qc}'.format(
            bam=bam_path,
            threads=threads,
            samstat_qc=qc_path,
        )
    )


samstat(mito_bam_path, mito_samstats_qc_path, threads)
samstat(no_mito_bam_path, no_mito_samstats_qc_path, threads)
frac_mito(no_mito_samstats_qc_path, mito_samstats_qc_path, frac_mito_qc_path)