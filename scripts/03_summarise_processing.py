#!/usr/bin/env python3
"""
Step 2: rank an IEDB MHC-I Processing table and compare the two sequences.

The tool is given both the wild-type and the mutant protein in one submission, so the
output carries a `seq #` column: 1 is whichever FASTA record came first. This reads the
table, prints the top N by processing total score for each sequence, and separately
reports every peptide spanning the substituted residue, which is usually not the same
set.

  ./03_summarise_processing.py results/iedb_9mer/peptide_table.csv --position 598
"""
import argparse, csv


def num(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("table")
    p.add_argument("--position", type=int, default=None,
                   help="residue number of the substitution, to list the peptides covering it")
    p.add_argument("--top", type=int, default=10)
    a = p.parse_args()

    rows = list(csv.DictReader(open(a.table)))
    seqs = sorted({r["seq #"] for r in rows})
    ranked = {}
    for s in seqs:
        sub = [r for r in rows if r["seq #"] == s]
        sub.sort(key=lambda r: -num(r, "processing total score"))
        for i, r in enumerate(sub, 1):
            r["_rank"] = i
        ranked[s] = sub

    length = rows[0]["peptide length"] if rows else "?"
    print("%s rows, %s-mers, %d sequence(s)" % (len(rows), length, len(seqs)))

    for s in seqs:
        print("\nseq %s, top %d by processing total score" % (s, a.top))
        print("  %-4s %-16s %-7s %8s %8s %8s %9s"
              % ("#", "peptide", "start", "proteas", "tap", "mhc", "total"))
        for r in ranked[s][:a.top]:
            print("  %-4d %-16s %-7s %8.3f %8.3f %8.3f %9.3f"
                  % (r["_rank"], r["peptide"], r["start"], num(r, "proteasome score"),
                     num(r, "tap score"), num(r, "mhc score"),
                     num(r, "processing total score")))

    # Identical top-10 lists across the two sequences are the expected outcome when the
    # substitution falls outside every strong-binding region, which is worth stating rather
    # than leaving the reader to compare two tables by eye.
    if len(seqs) == 2:
        a_top = [r["peptide"] for r in ranked[seqs[0]][:a.top]]
        b_top = [r["peptide"] for r in ranked[seqs[1]][:a.top]]
        print("\ntop %d identical between the two sequences: %s"
              % (a.top, "yes" if a_top == b_top else "no"))

    if a.position and len(seqs) == 2:
        print("\npeptides covering residue %d" % a.position)
        print("  %-9s %-16s %-6s %9s  %-16s %-6s %9s %9s"
              % ("window", "seq1", "rank", "total", "seq2", "rank", "total", "delta"))
        cov = {}
        for s in seqs:
            cov[s] = sorted((r for r in ranked[s]
                             if int(r["start"]) <= a.position <= int(r["end"])),
                            key=lambda r: int(r["start"]))
        for x, y in zip(cov[seqs[0]], cov[seqs[1]]):
            tx, ty = num(x, "processing total score"), num(y, "processing total score")
            print("  %-9s %-16s %-6d %9.3f  %-16s %-6d %9.3f %+9.3f"
                  % ("%s-%s" % (x["start"], x["end"]), x["peptide"], x["_rank"], tx,
                     y["peptide"], y["_rank"], ty, ty - tx))


if __name__ == "__main__":
    main()
