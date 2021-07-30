#!/usr/bin/env python

# Caleb Lareau, Surag Nair; Stanford University
# Implemented: 1 August 2021
# This program will error correct barcodes
# From 10x sequencing data from scATAC

##### IMPORT MODULES #####
import os
import re
import regex
import sys
import gzip
import itertools

from optparse import OptionParser
from multiprocessing import Pool, freeze_support
from itertools import repeat

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from fuzzysearch import find_near_matches

#### OPTIONS ####
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process raw .fastq reads and make data suitable for downstream processes"

opts.add_option("-a", "--fastq1", help="<Read1> Accepts fastq.gz")
opts.add_option("-b", "--fastq2", help="<Read2> Accepts fastq.gz")
opts.add_option("-c", "--fastq3", help="<Read3> Accepts fastq.gz")
opts.add_option("-f", "--barcodesFile", help="<gzip of the universe of valid barcodes (whitelist) to check")

opts.add_option("-m", "--max-barcode-dist", default=2, type=int, help="Maximum edit distance allowed between a barcode and an entry in whitelist")
opts.add_option("-n", "--nreads", default = 10000000, help="Number of reads to process in a given chunk")
opts.add_option("-r", "--ncores", default = 8, help="Number of cores for parallel processing.")

opts.add_option("-o", "--output", help="Output sample convention")

options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

# Define barcodes
barcodesfilepath = options.barcodesFile
with gzip.open(barcodesfilepath, "rt") as my_file:
    barcodesR = my_file.readlines()
barcodes = [barcode.rstrip() for barcode in barcodesR]
print("Found and imported " + str(len(barcodes)) + " barcodes")
global barcodes_set
barcodes_set = set(barcodes)


def batch_iterator(iterator, batch_size):
    """
    Returns lists of tuples of length batch_size.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def chunk_writer_gzip(filename, what):
    '''
    Basic function to write a chunk of a fastq file
    to a gzipped file
    '''
    with gzip.open(filename, 'wt') as out_write:
                out_write.writelines(what)
    return(filename)            

            
def formatRead(raw_barcode, corrected_barcode, title, sequence, quality):
    """
    Takes three components of fastq file + raw and corrected barcodes and stich them together in a string
    """
        
    # for bowtie, anything after space in name will go to SAM
    # remove existing comments as they may not be properly formatted
    mod_title = title.split(" ")[0]
    
    mod_title +=  " CB:Z:" + corrected_barcode + "\tCR:Z:" + raw_barcode

    return("@%s\n%s\n+\n%s\n" % (mod_title, sequence, quality))


#-----
# This code is modified from CellRanger-ATAC, base qualities are not factored in however
# https://github.com/10XGenomics/cellranger-atac/blob/main/mro/atac/stages/processing/attach_bcs/__init__.py
#-----
DNA_ALPHABET = 'AGCT'
ALPHABET_MINUS = {char: {c for c in DNA_ALPHABET if c != char} for char in DNA_ALPHABET}
ALPHABET_MINUS['N'] = set(DNA_ALPHABET)

def gen_nearby_seqs(seq, maxdist=3):
    
    allowed_indices = [i for i in range(len(seq)) if seq[i] != 'N']
    required_indices = tuple([i for i in range(len(seq)) if seq[i] == 'N'])
    mindist = len(required_indices)
    if mindist > maxdist:
        return

    for dist in range(mindist + 1, maxdist + 1):
        for modified_indices in itertools.combinations(allowed_indices, dist - mindist):
            indices = set(modified_indices + required_indices)

            for substitutions in itertools.product(
                    *[ALPHABET_MINUS[base] if i in indices else base for i, base in enumerate(seq)]):
                new_seq = ''.join(substitutions)
                if new_seq in barcodes_set:
                    yield new_seq
                    
#------ 

def correct_barcode(seq, maxdist=2):
    # returns (corrected sequence, edit distance) if found
    # else returns None
    if seq in barcodes_set:
        return (seq, 0)

    for test_str in gen_nearby_seqs(seq, maxdist):
        edit_dist = len([i for i in range(len(test_str)) if test_str[i]!=seq[i]])
        return(test_str, edit_dist)

    return (None, None)

    
def debarcode_trio(trio, max_barcode_dist):
    """
    Function that is called in parallel
    """
    # Parse out inputs
    listRead1 = trio[0]; listRead2 = trio[1]; listRead3 = trio[2]
    
    # parameters to return
    fq1 = ""
    fq2 = ""
    
    # Grab attributes
    title1 = listRead1[0]; sequence1 = listRead1[1]; quality1 = listRead1[2]
    title2 = listRead2[0]; sequence2 = listRead2[1]; quality2 = listRead2[2]
    title3 = listRead3[0]; sequence3 = listRead3[1]; quality3 = listRead3[2]

    corrected_barcode, edit_dist = correct_barcode(sequence2, maxdist=max_barcode_dist)
    #if(corrected_barcode != sequence2):
    #    print("was " + sequence2 + "   now: " + corrected_barcode)
    
    # Return the barcode with underscores + the biological sequence learned
    if corrected_barcode is not None:
        ofq1 = formatRead(sequence2, corrected_barcode, title1, sequence1, quality1)
        ofq2 = formatRead(sequence2, corrected_barcode, title3, sequence3, quality3)
        return(ofq1, ofq2, edit_dist)
    else:
        return None


def write_qc(edit_dists, outfname):
    with open(outfname, 'w') as f:
        f.write("Total barcodes processed\t{}\n".format(len(edit_dists)))
        f.write("Barcodes accepted\t{}\t{:.2f}\n".format(len(edit_dists)-edit_dists.count(-1), 100*(len(edit_dists)-edit_dists.count(-1))/len(edit_dists)))
        for x in range(max(edit_dists)+1):
            f.write("Barcodes with edit distance {}\t{}\t{:.2f}\n".format(x, edit_dists.count(x), 100*edit_dists.count(x)/len(edit_dists)))


if __name__ == "__main__": 
    ##### INPUTS #####
    #import pdb;pdb.set_trace()
    af = options.fastq1
    bf = options.fastq2
    cf = options.fastq3

    outname = options.output
    o = options.output
    cpu = int(options.ncores)
    n = int(options.nreads)

    # Parse input files
    extension = af.split('.')[-1]
    if extension == "fastq" or extension == "fq":
        sys.exist("Quitting... GZIP your .fastq files!")
    elif extension == "gz":
        print("Found supplied .fastq.gz files")
    else:
        sys.exit("ERROR! The input files (-a , -b, -c) a *.fastq.gz")
    print(options)
    edit_dists = []

    with gzip.open(af, "rt") as f1:
        with gzip.open(bf, "rt") as f2:
                with gzip.open(cf, "rt") as f3:

                    # Establish iterators
                    it1 = batch_iterator(FastqGeneralIterator(f1), n)
                    it2 = batch_iterator(FastqGeneralIterator(f2), n)
                    it3 = batch_iterator(FastqGeneralIterator(f3), n)

                    # iterate over batches of length n
                
                    for i, batch1 in enumerate(it1):
                        batch2 = it2.__next__()
                        batch3 = it3.__next__()
                        output = o 
            
                        # parallel process the barcode processing and accounting of failures.
                        pool = Pool(processes=cpu)

                        pm = pool.starmap(debarcode_trio, zip(zip(batch1, batch2, batch3), [options.max_barcode_dist]*len(batch1)))
                        pool.close()
            
                        # edit distances
                        edit_dists += [item[2] if item is not None else -1 for item in pm]

                        # Aggregate output
                        fastq1 = [item[0] for item in pm if item is not None]
                        fastq2 = [item[1] for item in pm if item is not None]

                        # Export one chunk in parallel
                        filename1 = output +'_R1.fastq.gz'
                        filename2 = output +'_R2.fastq.gz'
            
                        pool = Pool(processes=2)
                        toke = pool.starmap(chunk_writer_gzip, [(filename1, fastq1), (filename2, fastq2)])
                        pool.close()
    
    write_qc(edit_dists, output+"_barcode.qc.tsv")
    
