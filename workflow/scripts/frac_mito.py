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
    mapped = ''
    delimiter_pass_fail = ' + '
    with open(txt, 'r') as f:
        for line in f:
            if ' mapped (' in line:
                tmp3 = line.split(' mapped (')
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(delimiter_pass_fail)
                mapped = tmp3_1[0]

    return int(mapped)

def frac_mito(non_mito_samstat, mito_samstat, qc_path):
    Rn = parse_flagstat_qc(non_mito_samstat)
    Rm = parse_flagstat_qc(mito_samstat)

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