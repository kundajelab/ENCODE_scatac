import sys
import gzip
import pysam

def bam_to_frag(in_path, out_path, shift_plus=4, shift_minus=-4):
    """
    Convert BAM file to a fragment file format, while adding +4/-4 coordinate adjustment
    BAM should be pre-filtered for PCR duplicates, secondary alignments, and unpaired reads
    """

    input = pysam.AlignmentFile(in_path, "rb")
    with open(out_path, "w") as out_file:
        buf = []
        curr_pos = None
        for read in input:
            # if not ((read.flag & 80 == 80) or (read.flag & 160 == 160)): 
            if not ((read.flag & 64 == 64)): 
                continue # ignore coordinate-wise second read in pair
            
            chromosome = read.reference_name
            start = read.reference_start + shift_plus
            end = start + read.template_length + shift_minus
            cell_barcode = read.get_tag("CB")
            data = (chromosome, start, end, cell_barcode, 1, read.template_length, read.next_reference_start) ####
            pos = (chromosome, start)

            if pos == curr_pos:
                buf.append(data)
            else:
                buf.sort()
                for i in buf:
                    print(*i, sep="\t", file=out_file)
                buf.clear()
                buf.append(data)
                curr_pos = pos

if __name__ == '__main__':
    try:
        in_path, = snakemake.input
        out_path, = snakemake.output

        shift_plus = snakemake.params['shift_plus']
        shift_minus = snakemake.params['shift_minus']

        bam_to_frag(in_path, out_path, shift_plus=shift_plus, shift_minus=shift_minus)

    except NameError:
        bam_to_frag('/dev/stdin', '/dev/stdout')