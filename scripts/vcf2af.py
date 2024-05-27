import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import gzip

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Input file")

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


out2 = gzip.open(options.OUT+"_weight.csv.gz", "wt")
out1 = gzip.open(options.OUT+"_freq.csv.gz", "wt")
for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    a = l.rstrip().split()
    if l.startswith("#"):
        header = a[:2]
        header.extend(a[9:])
        out1.write("\t".join(header)+"\n")
        out2.write("\t".join(header)+"\n")
        continue
    LIST1, LIST2 = a[:2], a[:2]
    LIST1.extend([x.split(":")[-1] for x in a[9:]])
    LIST2.extend([x.split(":")[3] for x in a[9:]])
    out1.write("\t".join(LIST1)+"\n")
    out2.write("\t".join(LIST2)+"\n")
