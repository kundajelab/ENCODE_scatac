import pysam

def filter_mito(in_path, out_path, qc_path):
    """
    Removes mitochondrial alignments from BAM
    Calculates number of mapped mitochondrial and non-mitochondrial reads (not alignments)
    Assumes mitochondrial chromosome is "chrM"
    """

    infile = pysam.AlignmentFile(in_path, "rb")
    outfile = pysam.AlignmentFile(out_path, "wb", template=infile)

    num_mito = 0
    num_non_mito = 0
    for a in infile:
        if a.reference_name == "chrM":
            if a.flag & 260 == 0: # Alignment is mapped and is primary
                num_mito += 1
        else:
            if a.flag & 260 == 0:
                num_non_mito += 1
            outfile.write(a)

    with open(qc_path, "w") as qc_file:
        print("Non-Mitochondrial\tMitochondrial", file=qc_file)
        print(f"{num_non_mito}\t{num_mito}", file=qc_file)

if __name__ == '__main__':
    try:
        in_path, = snakemake.input
        out_path = snakemake.output["bam"]
        qc_path = snakemake.output["qc"]

        filter_mito(in_path, out_path, qc_path)

    except NameError:
        filter_mito('/dev/stdin', '/dev/stdout', '/dev/stderr')