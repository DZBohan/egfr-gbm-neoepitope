#!/usr/bin/env python3
"""
Step 3 (Expression Integration): pull mRNA and protein for one gene from cBioPortal
and compare mutated against non-mutated samples.

Nothing here is specific to EGFR; --gene and --study select the combination. Defaults
are EGFR in the CPTAC glioblastoma cohort, whose protein data is mass spectrometry
rather than RPPA. That matters: the assignment's final question is about mass
spectrometry, and an antibody array cannot speak to it.

  ./01_expression.py --study gbm_cptac_2021 --gene EGFR --entrez 1956 --change G598V
"""
import argparse, json, math, statistics as st, urllib.request

API = "https://www.cbioportal.org/api"


def post(path, payload):
    req = urllib.request.Request(
        API + path, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def values(study, profile, entrez):
    """Return {sample: float} for one molecular profile, dropping missing entries."""
    d = post("/molecular-profiles/%s_%s/molecular-data/fetch" % (study, profile),
             {"entrezGeneIds": [entrez], "sampleListId": study + "_all"})
    return {x["sampleId"]: float(x["value"]) for x in d
            if x.get("value") not in (None, "NA")}


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    return num / den if den else float("nan")


def summarise(name, group, data):
    a = [v for s, v in data.items() if s in group]
    b = [v for s, v in data.items() if s not in group]
    if not a:
        return "%s: no data in the mutated group" % name
    return ("%-9s mutated n=%-3d median %8.3f (%.3f to %.3f) | rest n=%-3d median %8.3f"
            % (name, len(a), st.median(a), min(a), max(a), len(b), st.median(b)))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--study", default="gbm_cptac_2021")
    p.add_argument("--gene", default="EGFR")
    p.add_argument("--entrez", type=int, default=1956)
    p.add_argument("--protein-profile", default="protein_quantification")
    p.add_argument("--mrna-profile", default="mrna")
    p.add_argument("--cna-profile", default="gistic")
    p.add_argument("--change", default=None,
                   help="amino acid change to isolate, e.g. G598V")
    p.add_argument("--out", default=None)
    a = p.parse_args()

    prot = values(a.study, a.protein_profile, a.entrez)
    mrna = values(a.study, a.mrna_profile, a.entrez)
    cna = {k: int(v) for k, v in values(a.study, a.cna_profile, a.entrez).items()}
    mut = post("/molecular-profiles/%s_mutations/mutations/fetch" % a.study,
               {"entrezGeneIds": [a.entrez], "sampleListId": a.study + "_all"})

    any_mut = {x["sampleId"] for x in mut}
    picked = ({x["sampleId"] for x in mut if x.get("proteinChange") == a.change}
              if a.change else any_mut)

    print("%s in %s" % (a.gene, a.study))
    print("  protein n=%d  mRNA n=%d  CNA n=%d  mutated samples n=%d"
          % (len(prot), len(mrna), len(cna), len(any_mut)))

    changes = {}
    for x in mut:
        changes[x.get("proteinChange")] = changes.get(x.get("proteinChange"), 0) + 1
    top = sorted(changes.items(), key=lambda kv: -kv[1])[:8]
    print("  mutation spectrum: " + ", ".join("%s x%d" % (k, v) for k, v in top))

    print("\nany mutation vs none")
    print("  " + summarise("protein", any_mut, prot))
    print("  " + summarise("mRNA", any_mut, mrna))

    if a.change:
        print("\n%s vs rest" % a.change)
        print("  " + summarise("protein", picked, prot))
        print("  " + summarise("mRNA", picked, mrna))

    # Amplification is the confounder that matters for a dosage-driven oncogene: if every
    # mutated sample is also amplified, a difference in expression says nothing about the
    # mutation. Stratify rather than report the raw contrast alone.
    amp = {s for s, v in cna.items() if v >= 2}
    print("\ncopy number")
    print("  amplified %d/%d samples; of the mutated, %d/%d are amplified"
          % (len(amp), len(cna), len(any_mut & amp), len(any_mut & set(cna))))
    for label, grp in (("amplified + mutated", any_mut & amp),
                       ("amplified, no mutation", amp - any_mut),
                       ("neither", set(cna) - amp - any_mut)):
        v = [prot[s] for s in grp if s in prot]
        if v:
            print("    %-24s n=%-3d protein median %+.3f" % (label, len(v), st.median(v)))

    both = sorted(set(prot) & set(mrna))
    if len(both) > 2:
        r = pearson([mrna[s] for s in both], [prot[s] for s in both])
        print("\nprotein vs mRNA: n=%d, Pearson r=%.3f" % (len(both), r))

    if a.out:
        rows = [{"sample": s, "mrna": mrna.get(s), "protein": prot.get(s),
                 "cna": cna.get(s), "mutated": s in any_mut, "picked": s in picked}
                for s in sorted(set(prot) | set(mrna))]
        json.dump(rows, open(a.out, "w"), indent=1)
        print("\nwrote %s (%d samples)" % (a.out, len(rows)))


if __name__ == "__main__":
    main()
