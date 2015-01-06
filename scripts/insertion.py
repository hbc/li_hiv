"""
pull out insertion sites from an alignment file aligned with bwa mem
BWA-mem sets the supplementary alignment flag (2048) and sets a SA tag for
reads with supplementary alignments
"""

import pysam
import re
from argparse import ArgumentParser

SUPPLEMENTARY_FLAG = 2048

def is_supplementary(read):
    return read.flag > SUPPLEMENTARY_FLAG

def get_tag(read, tag):
    tag = [x for x in read.tags if x[0] == tag]
    tag = tag[0] if tag else None
    tag = tag[1] if tag else None
    return tag

def get_SA_items(read):
    SA = [x for x in read.tags if x[0] == 'SA']
    SA = SA[0] if SA else None
    items = SA[1].split(",")
    chrom = items[0]
    pos = items[1]
    strand = items[2]
    cigar = items[3]
    bases = [int(x) for x in re.compile("[A-Z]").split(cigar) if x]
    return chrom, pos, int(pos) + sum(bases), strand

def get_SA_tag(read):
    SA = [x for x in read.tags if x[0] == 'SA']
    SA = SA[0] if SA else None
    return SA

def has_supplementary(read):
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

def parse_cigar_tuples(tuples):
    d = {"insertions": 0,
         "deletions": 0,
         "clipped": 0,
         "matched": 0,
         "other": 0}
    for t in tuples:
        if t[0] == 4 or t[0] == 5:
            d["clipped"] += t[1]
        elif t[0] == 1:
            d["insertions"] += t[1]
        elif t[0] == 2:
            d["deletions"] += t[1]
        elif t[0] == 0:
            d["matched"] += t[1]
        else:
            d["other"] += t[1]
    return d


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("bamfile", help="BAM file to call insertion sites from")
    parser.add_argument("virus_contig", help="Name of virus contig")
    args = parser.parse_args()

    HEADER = "\t".join(["rid", "first", "chrom", "pos", "end", "strand", "mapq",
                        "SA_chrom", "SA_pos", "SA_end", "SA_strand",
                        "as", "xs", "clipped", "insertions", "deletions",
                        "matched", "other"])
    FORMAT = ("{rid}\t{first}\t{chrom}\t{pos}\t{end}\t{strand}\t{mapq}\t{SA_chrom}\t"
              "{SA_pos}\t{SA_end}\t{SA_strand}\t{AS}\t{XS}\t{clipped}\t"
              "{insertions}\t{deletions}\t{matched}\t{other}")

    with pysam.Samfile(args.bamfile, "rb") as in_handle:
        contig = in_handle.gettid(args.virus_contig)
        print HEADER
        for read in in_handle:
            if read.is_secondary:
                continue
            if is_supplementary(read):
                continue
            if not has_supplementary(read):
                continue
            if not is_chimera(read, contig, args.virus_contig):
                continue
            if read.is_duplicate:
                continue
            rid = read.qname
            first = read.is_read1
            chrom = in_handle.getrname(read.tid)
            if chrom == args.virus_contig:
                pass
            strand = "-" if read.is_reverse else "+"
            pos = read.pos
            end = read.pos + read.reference_length
            AS = get_tag(read, "AS")
            XS = get_tag(read, "XS")
            mapq = read.mapq
            SA_chrom, SA_pos, SA_end, SA_strand = get_SA_items(read)
            cigar = parse_cigar_tuples(read.cigar)
            clipped = cigar["clipped"]
            insertions = cigar["insertions"]
            deletions = cigar["deletions"]
            matched = cigar["matched"]
            other = cigar["other"]
            print FORMAT.format(**locals())
