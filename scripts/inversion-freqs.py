import sys
from collections import defaultdict as d


def sync2freqh(x):
    ''' convert string in SYNC format to dictionary of freqencies where x is a string in sync format'''
    from collections import defaultdict as d
    if x == ".:.:.:.:.:." or x == "0:0:0:0:0:0":
        return ({"A": "na", "T": "na", "C": "na", "G": "na"}, 0)
        return "na", "na"
    nuc = ["A", "T", "C", "G"]
    counts = [int(X) for X in x.split(":")[:4]]
    if sum(counts) == 0:
        return ({"A": 0.0, "T": 0.0, "C": 0.0, "G": 0.0}, 0)
    CO = {X: Y for X, Y in zip(*[nuc, counts])}
    # print(CO, list(zip(*[nuc,counts])))
    h = d(float)
    for k, v in CO.items():
        h[k] = v / float(sum(CO.values()))
    return h, sum(CO.values())


def median(x):
    ''' calculate median '''
    mid = int(len(x)/2)
    sort = sorted(x)
    if len(x) == 0:
        return "NA"
    if len(x) % 2 == 0:
        lower = sort[mid-1]
        upper = sort[mid]
        return (float(lower)+float(upper))/2.0
    else:
        return sort[mid]


invmarker = open(sys.argv[1], "r")
data = open(sys.argv[2], "r")
names = sys.argv[3].split(",")
invh = d(list)

for l in invmarker:
    if l.startswith("inversion") or l.startswith("#"):
        continue
    a = l.split()
    invh[a[1] + ":" + a[2]] = [a[0], a[3]]

invdata = d(lambda: d(list))
Inv = []
for l in data:
    a = l.rstrip().split()
    pops = a[3:]
    if a[0] + ":" + a[1] not in invh:
        continue
    # print(len(pops),len(names))
    I, A = invh[a[0] + ":" + a[1]]
    Inv.append(I)
    for i in range(len(names)):
        # print(pops[i],A,sync2freqh(pops[i]))
        invdata[names[i]][I].append(sync2freqh(pops[i])[0][A])

Inversions = sorted(list(set(Inv)))

print("Sample\t" + "\t".join(Inversions))
for I, v in sorted(invdata.items()):
    AF = []
    for P, V in sorted(v.items()):
        # print(V)
        AF.append(
            str(median([x for x in V if x != "na"])))
    print(I + "\t" + "\t".join(AF))
