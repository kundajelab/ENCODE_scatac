import sys
import pysam

# Convert bam file to a fragment file format, while adding +4/-4 coordinate adjustment

input = pysam.AlignmentFile(sys.argv[1], "rb")

prev_read = next(input)
while True:
    read = next(input, None)
    if read is None or prev_read is None:
        break
    if read.query_name != prev_read.query_name:
        print("Read with missing mate:", prev_read.query_name, file=sys.stderr)
        prev_read = read
        continue

    if prev_read.is_reverse:
        prev_read, read = read, prev_read

    chromosome = read.reference_name
    start = prev_read.reference_start + 4
    end = start + prev_read.template_length - 4
    cell_barcode = read.get_tag("CB")
    print(chromosome, start, end, cell_barcode, sep="\t")
    prev_read = next(input, None)

