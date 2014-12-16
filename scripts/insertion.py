"""
pull out insertion sites from an alignment file aligned with bwa mem
BWA-mem sets the supplementary alignment flag (2048) and sets a SA tag for
reads with supplementary alignments
"""

import pysam
from argparse import ArgumentParser

SUPPLEMENTARY_FLAG = 2048

def get_SA_tag(read):
    SA = [x for x in read.tags if x[0] == 'SA']
    SA = SA[0] if SA else None
    return SA

def is_supplementary(read):
    return get_SA_tag(read)

def supplementary_contig(read):
    SA = get_SA_tag(read)
    contig = SA[1].split(',')[0]
    return contig

def is_chimera(read, contig, contig_name):
    scontig = supplementary_contig(read)
    if (read.rname == contig) & (scontig != contig_name):
        return True
    elif (read.rname != contig) & (scontig == contig_name):
        return True
    else:
        return False


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("bamfile", help="BAM file to call insertion sites from")
    parser.add_argument("virus_contig", help="Name of virus contig")
    args = parser.parse_args()

    with pysam.Samfile(args.bamfile, "rb") as in_handle:
        contig = in_handle.gettid(args.virus_contig)
        for read in in_handle:
            if read.is_secondary:
                continue
            if not is_supplementary(read):
                continue
            if not is_chimera(read, contig, args.virus_contig):
                continue
            print read
