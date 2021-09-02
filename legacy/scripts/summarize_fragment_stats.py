import re

input = snakemake.input

with open(snakemake.output[0], "w") as out:
    out.write("Barcode matching: ")
    out.write(open(input.barcode_matching, "r").readline())

    adapter_stats = open(input.adapter_stats).read()
    adapter_pass = int(re.search("reads passed filter: (\d+)", adapter_stats).group(1)) // 2
    adapter_trimmed = int(re.search("reads with adapter trimmed: (\d+)", adapter_stats).group(1)) // 2
    adapter_percent = adapter_trimmed / adapter_pass * 100
    out.write(f"Adapter trimming: {adapter_trimmed}/{adapter_pass} reads trimmed ({adapter_percent:.2f}%)\n")


    total_align = int(open(input.total_align).read())
    align_percent = total_align/adapter_pass * 100
    out.write(f"Total high-quality alignments: {total_align}/{adapter_pass} ({align_percent:.2f}%)\n")

    unique_align = int(open(input.unique_align).read())
    unique_percent = unique_align / total_align * 100
    out.write(f"Total unique fragments: {unique_align}/{total_align} ({unique_percent:.2f}%)\n")

    chrM_count = unique_align - int(open(input.no_mito).read())
    chrM_percent = chrM_count / unique_align * 100
    out.write(f"Mitochondrial fragments: {chrM_count}/{unique_align} ({chrM_percent:.2f}%)\n")



    #TotalReadPairs\\tDistinctReadPairs\\tOneReadPair\\tTwoReadPairs\\tNRF=Distinct/Total\\tPBC1=OnePair/Distinct\\tPBC2=OnePair/TwoPair
    total_reads, unique_reads, _, _, NRF, PCB1, PCB2 = open(input.duplicate_stats).readlines()[1].strip().split("\t")
    total_reads = int(total_reads)
    unique_reads = int(unique_reads)

    blacklist_reads = unique_align - chrM_count - unique_reads 
    blacklist_percent = blacklist_reads / (unique_align - chrM_count)
    out.write(f"Blacklist fragments: {blacklist_reads}/{unique_align - chrM_count} ({blacklist_percent:.2f}%)\n")
    out.write(f"PCR bottleneck coefficients: PCB1 {float(PCB1):.2f}, PCB2 {float(PCB2):.2f}\n")