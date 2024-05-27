import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--diag", dest="DI", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


DIAG = d(lambda: d(str))

for l in load_data(options.DI):
    a = l.rstrip().split()
    if l.startswith("Chrom"):
        continue
    DIAG[a[0]][a[1]] = a[2]

for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    a = l.rstrip().split()

    # get names from headers
    if l.startswith("#"):
        header = a[9:]
        print("ID\tChrom\tPos\tAllele\tFreq")
        continue

    if a[1] not in DIAG[a[0]]:
        continue

    # obtain alleles
    REF = a[3]
    ALT = a[4]
    ALLELE = [REF, ALT]

    # ignore tri- and tetra-allelic SNPs
    if len(ALT) > 1 or len(REF) > 1:
        continue

    pops = a[9:]
    for i in range(len(pops)):
        if pops[i].split(":")[-1] == ".":
            continue
        if DIAG[a[0]][a[1]] == REF:
            print(header[i]+"\t"+"\t".join(a[:2])+"\t"+DIAG[a[0]]
                  [a[1]]+"\t"+str(1-float(pops[i].split(":")[-1])))
        elif DIAG[a[0]][a[1]] == ALT:
            print(header[i]+"\t"+"\t".join(a[:2])+"\t"+DIAG[a[0]]
                  [a[1]]+"\t"+str(float(pops[i].split(":")[-1])))
