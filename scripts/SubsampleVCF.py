import sys
from collections import defaultdict as d
import random
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--snps", dest="SNPs", help="Number of SNPs")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    """ import data either from a gzipped or or uncrompessed file or from STDIN"""
    import gzip
    
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


SNPs = d(str)

for l in load_data(options.IN):
    if l.startswith("#"):
        print(l.rstrip())
        continue
    a = l.rstrip().split()
    SNPs[a[0] + "_" + a[1]] = l.rstrip()

KEYS = random.sample(SNPs.keys(), int(options.SNPs))

KEYS2 = d(list)

for k in KEYS:
    C, P = k.split("_")
    KEYS2[C].append(int(P))

for k, v in sorted(KEYS2.items()):
    for p in sorted(v):
        print(SNPs[k + "_" + str(p)])
