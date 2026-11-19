import re, sys, collections
tax_file   = sys.argv[1]   # gtdb_classification.tsv
split_file = sys.argv[2]   # fam00030_split_c0.3.tsv (family <tab> genome)

def key(s):
    m = re.search(r'(\d{9}\.\d)', s)
    return m.group(1) if m else s

species = {}
with open(tax_file) as f:
    next(f)  # header
    for line in f:
        p = line.rstrip("\n").split("\t")
        sp = [x for x in p[1].split(";") if x.startswith("s__")]
        species[key(p[0])] = sp[0][3:] if sp and sp[0] != "s__" else "unassigned"

tab = collections.defaultdict(collections.Counter)
missing = 0
with open(split_file) as f:
    for line in f:
        fam, g = line.rstrip("\n").split("\t")[:2]
        sp = species.get(key(g))
        if sp is None:
            missing += 1; sp = "NOT_IN_GTDB_TABLE"
        tab[fam][sp] += 1

for fam in sorted(tab):
    print(f"\n{fam}  ({sum(tab[fam].values())} genomes)")
    for sp, n in tab[fam].most_common():
        print(f"   {n:3d}  {sp}")
print(f"\nGenomes without a match: {missing}")
