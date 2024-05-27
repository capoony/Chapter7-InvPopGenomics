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
parser.add_option("--breakpoints", dest="BP", help="Input file")
parser.add_option("--range", dest="RA", help="Input file")
parser.add_option("--output", dest="out", help="Input file")
parser.add_option("--MinCov", dest="MC",
                  help="numerical parameter", default=4)
parser.add_option("--Variant", dest="VA",
                  help="numerical parameter", default=1)

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


DiaOut = open(options.out+"_diag.txt", "wt")
VaOut = open(options.out+"_variants.txt", "wt")
DiaOut.write("Chrom\tPos\tINV\tST\n")
VaOut.write("Sample\tChrom\tPos\tVariant\n")
Chrom, Start, End = options.BP.split(",")
print(Chrom, Start, End)
Range = int(options.RA)
MINSt = int(Start)-Range
MAXSt = int(Start)+Range
MINEn = int(End)-Range
MAXEn = int(End)+Range
VA = d(str)
for l in load_data(options.VA):
    if l.startswith("Stock"):
        continue
    a = l.rstrip().split(",")
    VA[a[0]] = a[2]

Positions = d(str)
for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    a = l.rstrip().split()

    # get names from headers
    if l.startswith("#"):
        header = [x.split("/")[-1].split(".bam")[0] for x in a[9:]]
        continue

    if a[0] != Chrom:
        continue
    if (int(a[1]) >= MINSt and int(a[1]) <= MAXSt) or (int(a[1]) >= MINEn and int(a[1]) <= MAXEn):

        # obtain alleles
        REF = a[3]
        ALT = a[4]
        ALLELE = [REF, ALT]

        # ignore tri- and tetra-allelic SNPs
        if len(ALT) > 1 or len(REF) > 1:
            continue

        pops = a[9:]
        KEY = a[8].split(":")
        TestGT = []
        VaGT = d(list)
        TYPE = d(str)

        for i in range(len(header)):
            if header[i] not in VA:
                continue
            VAL = pops[i].split(":")
            DICT = dict(zip(*(KEY, VAL)))

            if DICT["GT"] == ".":
                continue

            if DICT["GT"] == "0":
                A = 0
                B = 1
            else:
                A = 1
                B = 0

            # only consider DICT["GT"] if (1) the read-depth > than MC, (2) the Posterior Likelihood of the DICT["GT"] is > 50 and the PL of the other (non-called) DICT["GT"] is < 30 otherwise mark as ambiguous
            # and int(PL[A]) > 50 and int(PL[B]) < 30:
            if int(DICT["DP"]) >= int(options.MC):
                VaGT[VA[header[i]]].append(DICT["GT"])
                TestGT.append(DICT["GT"])
        # print(a[:2])
        if len(TestGT) != len(VA):
            # print(VA, TestGT)
            continue

        if len(list(set(VaGT["ST"]))) == 1 and len(list(set(VaGT["INV"]))) == 1 and VaGT["INV"][0] != VaGT["ST"][0]:
            TYPE[VaGT["INV"][0]] = "INV"
            TYPE[VaGT["ST"][0]] = "ST"
            DiaOut.write(a[0]+"\t"+a[1]+"\t"+ALLELE[int(VaGT["INV"][0])] +
                         "\t"+ALLELE[int(VaGT["ST"][0])]+"\n")
        else:
            # print(VaGT["ST"], VaGT["INV"])
            continue

        # loop through all samples
        for i in range(len(header)):
            if header[i] in VA:
                continue
            VAL = pops[i].split(":")
            DICT = dict(zip(*(KEY, VAL)))

            if DICT["GT"] == ".":
                continue

            if DICT["GT"] == "0":
                A = 0
                B = 1
            else:
                A = 1
                B = 0

            # only consider DICT["GT"] if (1) the read-depth > than MC, (2) the Posterior Likelihood of the DICT["GT"] is > 50 and the PL of the other (non-called) DICT["GT"] is < 30 otherwise mark as ambiguous
            # and int(PL[A]) > 50 and int(PL[B]) < 30:
            if int(DICT["DP"]) >= int(options.MC):
                VaOut.write(header[i]+"\t"+a[0]+"\t" +
                            a[1]+"\t"+TYPE[DICT["GT"]]+"\n")
