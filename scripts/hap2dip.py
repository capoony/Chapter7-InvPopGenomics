from optparse import OptionParser, OptionGroup
import gzip

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--logical", dest="log",
                  help="logical parameter", action="store_true")
parser.add_option("--param", dest="param",
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


hap_to_dip = {'0': '0/0', '1': '1/1', '.': './.'}

with load_data(options.IN) as file, gzip.open(options.OUT, "wt") as s_file:
    for line in file:
        if line[0] == '#':
            s_file.write(line)
        else:
            rowstr = line.split('\t')[9:]
            for i in range(len(rowstr)):
                rowstr[i] = hap_to_dip[rowstr[i][0]] + rowstr[i][1:]
            row = '\t'.join(str(line).split('\t')[
                            :9]) + '\t' + '\t'.join(rowstr)
            s_file.write(row)
