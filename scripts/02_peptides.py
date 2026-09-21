#!/usr/bin/env python3
"""
Steps 2 and 4: translate a CDS pair, enumerate the k-mers spanning the substitution,
and emit the paired input the Peptide Variant Comparison server expects.

Two mistakes this exists to prevent, both of which happened on the KRAS run:

  - Submitting the CDS instead of the protein. A, T, C and G are all valid amino acid
    letters, so the predictor accepts nucleotides and silently scans the wrong sequence.
    --check reports the residue numbering so the substitution can be confirmed.
  - Trusting a file name. A file labelled G12D contained G13D; KRAS has glycines at both
    positions. This reports the substitution it actually finds rather than the one claimed.

  ./02_peptides.py --wt wt.fa --mut mut.fa --length 9 --pairs pairs.tsv
"""
import argparse

CODONS = {}
_b = "TCAG"
_aa = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
_i = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODONS[_x + _y + _z] = _aa[_i]
            _i += 1

# The pMHC immunogenicity model masks these by default: they anchor into the MHC groove
# and are not seen by the TCR, so a substitution there cannot move that model's score.
MASKED = {1, 2}


def read_fasta(path):
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper().replace(" ", "")


def translate(dna):
    dna = "".join(c for c in dna if c in "ATCG")
    aa = "".join(CODONS.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))
    return aa.split("*")[0] if "*" in aa else aa


def as_protein(seq):
    """Accept either a CDS or a protein sequence; translate only if it looks like DNA."""
    return translate(seq) if set(seq) <= set("ATCGN") else seq


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--wt", required=True)
    p.add_argument("--mut", required=True)
    p.add_argument("--length", type=int, default=9)
    p.add_argument("--pairs", default=None)
    a = p.parse_args()

    wt = as_protein(read_fasta(a.wt))
    mu = as_protein(read_fasta(a.mut))
    if len(wt) != len(mu):
        raise SystemExit("lengths differ (%d vs %d); this script handles substitutions only"
                         % (len(wt), len(mu)))

    diffs = [(i + 1, x, y) for i, (x, y) in enumerate(zip(wt, mu)) if x != y]
    if len(diffs) != 1:
        raise SystemExit("expected exactly one substitution, found %d: %s" % (len(diffs), diffs))
    pos, ref, alt = diffs[0]
    print("protein length %d, substitution %s%d%s" % (len(wt), ref, pos, alt))
    print("context  WT %s" % wt[max(0, pos - 6):pos + 5])
    print("        MUT %s" % mu[max(0, pos - 6):pos + 5])

    k = a.length
    rows = []
    print("\n%d-mers spanning residue %d:" % (k, pos))
    print("  %-7s %-12s %-12s %-5s %s" % ("window", "wild type", "mutant", "pos", "seen by pMHC model"))
    for s in range(max(0, pos - k), pos):
        if s + k > len(wt):
            break
        A, B = wt[s:s + k], mu[s:s + k]
        if A == B:
            continue
        pp = pos - s
        masked = pp in MASKED or pp == k
        print("  %-7s %-12s %-12s P%-4d %s"
              % ("%d-%d" % (s + 1, s + k), A, B, pp, "no, masked" if masked else "yes"))
        rows.append((A, B))

    if a.pairs:
        with open(a.pairs, "w") as f:
            for A, B in rows:
                f.write("%s\t%s\n" % (A, B))
        print("\nwrote %s (%d pairs, wild type in column 1)" % (a.pairs, len(rows)))


if __name__ == "__main__":
    main()
