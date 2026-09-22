#!/usr/bin/env python3
"""Build the assignment from archived result tables; no network or third-party packages."""
import base64, csv, json, math, statistics as st
from pathlib import Path
from html import escape as e

ROOT = Path(__file__).resolve().parents[1]
def readcsv(path, delimiter=','):
    with (ROOT/path).open(newline='') as f: return list(csv.DictReader(f, delimiter=delimiter))
def num(r,k): return float(r[k])
def f(v,n=3): return f'{float(v):.{n}f}'
def table(caption,headers,rows):
    return '<div class="table-wrap"><table><caption>'+caption+'</caption><thead><tr>'+''.join('<th scope="col">'+e(str(x))+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+e(str(x))+'</td>' for x in r)+'</tr>' for r in rows)+'</tbody></table></div>\n'
def fig(name,caption):
    image_data = base64.b64encode((ROOT/'figures'/name).read_bytes()).decode('ascii')
    return f'<figure><img src="data:image/jpeg;base64,{image_data}" alt="{e(caption)}"><figcaption>{caption}</figcaption></figure>\n'
def source(path): return f'<p class="source">Data: <a href="../{path}">{path}</a>.</p>'
def cite(url,label): return f'<a href="{url}">{label}</a>'
def fasta(path):
    lines=(ROOT/path).read_text().splitlines();seq=''.join(lines[1:]);parts=[]
    for start in range(0,len(seq),60):parts.append(''.join(f'<span class="variant">{aa}</span>' if i==597 else aa for i,aa in enumerate(seq[start:start+60],start)))
    return '<pre class="fasta">'+e(lines[0])+'\n'+'\n'.join(parts)+'</pre>',seq
fw,wt=fasta('data/egfr_wt.fa');fm,mt=fasta('data/egfr_g598v.fa')
assert len(wt)==len(mt)==1210 and [(i+1,a,b) for i,(a,b) in enumerate(zip(wt,mt)) if a!=b]==[(598,'G','V')]
ranked={}
for k in [9,13]:
    rows=readcsv(f'results/iedb_{k}mer/peptide_table.csv')
    for s,seq in [('1',wt),('2',mt)]:
        rr=sorted((r for r in rows if r['seq #']==s),key=lambda r:-num(r,'processing total score'))
        assert len(rr)==1210-k+1
        for rank,r in enumerate(rr,1):
            r['rank']=rank
            assert r['peptide']==seq[int(r['start'])-1:int(r['end'])]
            assert math.isclose(num(r,'processing total score'),sum(num(r,c) for c in ['proteasome score','tap score','mhc score']),abs_tol=1e-10)
            assert math.isclose(num(r,'mhc score'),-math.log10(num(r,'netmhcpan_ba ic50')),abs_tol=1e-10)
        ranked[k,s]=rr
vc=readcsv('results/iedb_variant_comparison/peptide_table.csv')
tcia=readcsv('results/tcia_gbm_egfr/neoantigens_gbm_egfr.tsv','\t');tpeps={r['peptide'] for r in tcia}
expr={k:{x['sampleId']:float(x['value']) for x in json.loads((ROOT/f'results/cbioportal_api/{k}.json').read_text())} for k in ['mrna','protein_quantification','gistic','mrna_zscores']}
mutations=json.loads((ROOT/'results/cbioportal_api/mutations.json').read_text());anymut={r['sampleId'] for r in mutations};picked={r['sampleId'] for r in mutations if r['proteinChange']=='G598V'}
common=set(expr['mrna'])&set(expr['protein_quantification']);amp={s for s,v in expr['gistic'].items() if v>=2}
with (ROOT/'results/cbioportal_api/sample_table.tsv').open('w',newline='') as out:
    writer=csv.writer(out,delimiter='\t');writer.writerow(['sample','mrna_profile_value','protein_profile_value','mrna_zscore','gistic','any_EGFR_mutation','G598V'])
    for s in sorted(common):writer.writerow([s,expr['mrna'][s],expr['protein_quantification'][s],expr['mrna_zscores'][s],expr['gistic'].get(s,'NA'),int(s in anymut),int(s in picked)])
processing='https://tools.iedb.org/processing/help/'
peters='https://tools.iedb.org/static/pdf/peters_2003_joi.pdf'
calis='https://tools.iedb.org/immunogenicity/help/'
icer='https://academic.oup.com/narcancer/article/6/1/zcae002/7591107'
dtu='https://services.healthtech.dtu.dk/services/ICERFIRE-1.0/'
html=['''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>EGFR G598V in Glioblastoma: An Immunology Assignment</title>
<style>
:root{color-scheme:light;--ink:#192d3c;--accent:#24566b;--line:#ccd8df}
*{box-sizing:border-box}body{margin:0;background:#edf2f4;color:var(--ink);font:16px/1.65 Georgia,"Times New Roman",serif}main{max-width:1160px;margin:32px auto;background:#fff;padding:48px 54px;box-shadow:0 3px 22px #182d3c12}h1,h2,h3,nav,th,caption,.eyebrow{font-family:Arial,Helvetica,sans-serif}h1{font-size:2.25rem;line-height:1.2;margin:8px 0 20px}h2{font-size:1.55rem;border-top:3px solid var(--accent);padding-top:20px;margin-top:50px}h3{font-size:1.1rem;margin-top:28px}.eyebrow{color:var(--accent);letter-spacing:.12em;font-size:.8rem;text-transform:uppercase}.byline{font-family:Arial,Helvetica,sans-serif;font-size:1.05rem;color:var(--accent);margin:-8px 0 14px}a{color:#165c7a}nav{display:flex;flex-wrap:wrap;gap:12px;margin:24px 0}nav a{background:#edf4f7;padding:7px 12px;text-decoration:none}p{margin:12px 0}.summary{border-left:4px solid var(--accent);background:#eff6f8;padding:14px 20px;margin:20px 0}.limitation{border-left:4px solid #997343;background:#faf5ed;padding:14px 20px}.table-wrap{overflow-x:auto;margin:20px 0}table{border-collapse:collapse;width:100%;font:13px/1.45 Arial,Helvetica,sans-serif;font-variant-numeric:tabular-nums}caption{text-align:left;font-weight:bold;font-size:14px;margin-bottom:10px}th,td{padding:8px 9px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}th{background:#eaf1f5}tbody tr:nth-child(even){background:#f7f9fa}td{white-space:nowrap}pre,code{font-family:"DejaVu Sans Mono",Consolas,monospace}code{font-size:.9em}pre{background:#f4f6f8;padding:18px;overflow-x:auto;font-size:14px;line-height:1.65}.variant{color:#c51627;font-weight:bold;text-decoration:underline}figure{margin:28px 0}img{display:block;max-width:100%;height:auto;margin:auto;max-height:850px}figcaption,.source{font-size:.88rem;color:#435866}figcaption{margin-top:12px}.source{overflow-wrap:anywhere}li{margin:9px 0}footer{border-top:1px solid var(--line);margin-top:40px;padding-top:16px;color:#435866;font-size:.85rem}@media(max-width:700px){main{margin:0;padding:22px 18px}h1{font-size:1.8rem}pre{font-size:11px}table{font-size:12px}}@page {
 size: A4;
 margin: 17mm 15mm 19mm;
 @bottom-center { content: "Page " counter(page) " of " counter(pages); font: 9pt Arial, sans-serif; color: #435866; }
}
@media print {
 body { background: white; font-size: 10pt; line-height: 1.45; }
 main { max-width: none; margin: 0; padding: 0; box-shadow: none; }
 h1 { font-size: 23pt; }
 h2 { font-size: 15pt; margin-top: 20pt; padding-top: 10pt; break-after: avoid; }
 h3 { font-size: 11pt; margin-top: 14pt; break-after: avoid; }
 p { margin: 7pt 0; orphans: 3; widows: 3; }
 table, figure, tr, .summary, .limitation, pre { break-inside: avoid; }
 .table-wrap { overflow: visible; margin: 12pt 0; break-inside: avoid; }
 table { width: 100%; font-size: 8pt; line-height: 1.3; }
 thead { display: table-header-group; }
 th, td { padding: 4pt 3pt; font-size: inherit; }
 th { white-space: normal; }
 td { white-space: normal; overflow-wrap: anywhere; }
 caption { font-size: 9pt; margin-bottom: 6pt; break-after: avoid; }
 figure { margin: 8pt 0; }
 img { max-width: 100%; max-height: 105mm; width: auto; height: auto; }
 figcaption { margin-top: 6pt; font-size: 8.5pt; line-height: 1.35; }
 a { color: inherit; text-decoration: none; }
 pre { font-size: 8.2pt; line-height: 1.35; padding: 10pt; white-space: pre; overflow: visible; }
 .variant { color: #c51627 !important; print-color-adjust: exact; }
 .summary, .limitation { padding: 9pt 12pt; margin: 12pt 0; }
 .source { display: none; }
 footer { display: none; }
 nav { display: none; }
}

</style></head><body><main>
<header><p class="eyebrow">Immunology course assignment</p><h1>EGFR G598V in Glioblastoma</h1><p class="byline">Bohan Zhang</p><p>From protein sequence and antigen processing to expression, immunogenicity, and the limits of database validation.</p></header>
<nav aria-label="Assignment steps"><a href="#step1">1. Sequences</a><a href="#step2">2. Processing</a><a href="#step3">3. Expression</a><a href="#step4">4. Immunogenicity</a><a href="#step5">5. Validation</a></nav>
<p>This analysis evaluates a single EGFR substitution using HLA-A*02:01 as a representative common class I allele. Scores describe computational predictions; expression measurements are gene-level observations. Neither establishes that a mutant peptide is naturally presented or recognized by T cells.</p>
<section id="step1"><h2>Step 1. Generate Neoantigen and Native Antigen Sequences</h2>
<p>The selected alteration is EGFR p.Gly598Val (G598V) in glioblastoma. The protein FASTA records identify transcript <code>ENST00000275493</code> and contain 1,210 amino acids each. Direct comparison finds exactly one difference: glycine in the wild type is replaced by valine at residue 598, using one-based numbering of the full precursor protein.</p>

''']
html.append('<p>Sequence context, residues 586 to 610:</p><pre>WT     '+wt[585:597]+'<span class="variant">G</span>'+wt[598:610]+'\nG598V  '+mt[585:597]+'<span class="variant">V</span>'+mt[598:610]+'</pre>')
html.append('<h3>Wild-type protein</h3>'+fw+'<h3>G598V protein</h3>'+fm+'<p>The homologous residue at position 598 is red and underlined in both FASTA records.</p></section>')
html.append('''<section id="step2"><h2>Step 2. Intracellular Processing and MHC Class I Binding</h2>
<p>IEDB T Cell Prediction Class I was run with MHC-I Processing Basic Processing Predictions using the immunoproteasome option, with NetMHCpan 4.1 BA running simultaneously for HLA-A*02:01. Sequence 1 is wild type and sequence 2 is G598V. Each full-length protein yields 1,202 nine-residue windows or 1,198 thirteen-residue windows.</p>''')
html.append(f'<p>Tables are sorted by decreasing <code>processing total score</code>. In these exports, total = proteasome + TAP + MHC, whereas <code>processing score</code> alone = proteasome + TAP. MHC score = −log<sub>10</sub>(predicted IC50 in nM). Higher scores and lower IC50 values are favorable, but these are not calibrated probabilities of cleavage, transport, or presentation. See {cite(processing,"IEDB score definitions")}.</p>')
for k in [9,13]:
    html.append(f'<h3>Top ten {k}-mers from each protein</h3>')
    for s,label in [('1','Wild type'),('2','G598V')]:
        rr=ranked[k,s]
        html.append(table(f'{label}: {k}-mers ranked by processing total score',['Rank','Peptide','Start','Proteasome','TAP','MHC','Total','IC50 (nM)'],[[r['rank'],r['peptide'],r['start']]+[f(num(r,c)) for c in ['proteasome score','tap score','mhc score','processing total score']]+[f(num(r,'netmhcpan_ba ic50'),2)] for r in rr[:10]]))
    html.append(source(f'results/iedb_{k}mer/peptide_table.csv'))
html.append('<div class="summary">The top ten sequences and their scores are identical between wild type and G598V at each length. None spans residue 598. They are high-scoring EGFR peptides, not mutation-specific neoantigen candidates.</div><h3>Windows spanning the mutation</h3><p>Ranks below are ordinal positions within each protein and peptide length, calculated from unrounded totals. They are not percentile ranks. Displayed deltas are calculated before rounding.</p>')
for k in [9,13]:
    aa={r['start']:r for r in ranked[k,'1']};bb={r['start']:r for r in ranked[k,'2']};rr=[]
    for start in range(598-k+1,599):
        a,b=aa[str(start)],bb[str(start)]
        rr.append([f'{start}–{start+k-1}',a['peptide'],a['rank'],f(num(a,'processing total score')),b['peptide'],b['rank'],f(num(b,'processing total score')),f"{num(b,'processing total score')-num(a,'processing total score'):+.3f}"])
    html.append(table(f'Mutation-containing {k}-mers',['Window','WT peptide','WT rank','WT total','Mutant peptide','Mutant rank','Mutant total','Δ total'],rr))
html.append('''<p>The largest total-score increase occurs when G598V is the C-terminal residue: HCVKTCPAG → HCVKTCPAV at nine residues, and IDGPHCVKTCPAG → IDGPHCVKTCPAV at thirteen. A valine side chain is more compatible with a hydrophobic HLA-A*02:01 C-terminal anchor than glycine, so improved binding is biologically plausible. This is a mechanistic interpretation, not a measured structure or affinity. The total-score improvement also contains processing effects:</p>''')
rr=[]
for k in [9,13]:
 a,b=[next(r for r in ranked[k,s] if r['end']=='598') for s in ['1','2']]
 rr.append([k]+[f"{num(b,c)-num(a,c):+.6f}" for c in ['proteasome score','tap score','mhc score','processing total score']]+[f(num(a,'netmhcpan_ba ic50'),2),f(num(b,'netmhcpan_ba ic50'),2)])
html.append(table('Decomposition of the C-terminal G598V effect',['Length','Δ proteasome','Δ TAP','Δ MHC','Δ total','WT IC50','Mutant IC50'],rr))
html.append('<p>Even after improvement, these predicted affinities remain above 500 nM. None of the nine mutant 9-mers meets that affinity cutoff. VVMGENNTL has the highest processing total among the nine mutant windows, but its BA IC50 is 1,808.14 nM. Relative rank alone should not be called strong binding.</p><h3>Why do the 9-mer and 13-mer results differ?</h3>')
rr=[]
for c in ['proteasome score','tap score','mhc score','processing total score']:
 rr.append([c]+[f(st.mean(num(r,c) for r in ranked[k,'1']),6) for k in [9,13]]+[f(max(num(r,c) for r in ranked[k,'1']),6) for k in [9,13]])
html.append(table('Wild-type score distributions',['Column','9-mer mean','13-mer mean','9-mer maximum','13-mer maximum'],rr))
html.append(f'''<p>The numbers of WT peptides with predicted IC50 &lt; 500 nM are 31/1,202 and 11/1,198, respectively. The mean total decreases by 0.040401 for 13-mers; the MHC term accounts for 0.038719 of this difference. This arithmetic describes these two window sets, not a controlled experiment isolating length.</p>
<p>The proteasome score concerns the C-terminal cleavage site and its local protein context. Identical cleavage sites have identical proteasome scores in these runs. TAP scoring emphasizes terminal sequence and potential N-terminal precursors, so changing the window can change its score even when the C terminus is fixed. Similar whole-protein means do not demonstrate biological independence from length. For example, the WT windows ending at 598 have TAP scores of −0.665889 at nine residues and −0.822274 at thirteen. The terminal-sequence approximation is described by {cite(peters,'Peters et al. (2003)')}.</p>
<p>Class I grooves commonly accommodate short peptides; longer ligands can use altered conformations, including central bulging. This gives a plausible context for the weaker binding-score tail here, but neither a universally flat 9-mer nor a fixed energetic penalty for every 13-mer can be inferred. {cite('https://pmc.ncbi.nlm.nih.gov/articles/PMC3653566/','Structural studies of peptide length')} demonstrate the importance of peptide conformation.</p>
<p class="limitation">Only one HLA allele and two lengths were examined. Scores omit patient-specific antigen abundance and much of cellular processing, including variable trimming and degradation. Top-ranked full-protein windows may also arise from regions handled differently during protein maturation. Absolute scores share a computational definition across runs but are not calibrated cross-length probabilities; ranks refer to different candidate sets.</p></section>''')
html.append('''<section id="step3"><h2>Step 3. Expression Integration</h2>
<p>The representative cohort is <a href="https://www.cbioportal.org/study/summary?id=gbm_cptac_2021">cBioPortal <code>gbm_cptac_2021</code></a>, Glioblastoma (CPTAC, Cell 2021), queried for EGFR (Entrez 1956). Protein abundance was measured by mass spectrometry. The screenshots show the Plots tab with Protein vs mRNA and G598V entered in Search Mutation(s). They contain 99 samples with both measurements.</p>''')
html.append(fig('step3_protein_vs_mrna_g598v.jpg','Figure 1. EGFR protein abundance versus log2-transformed mRNA z-score. Enlarged points identify G598V samples. The screenshot reports Pearson r = 0.92 and Spearman ρ = 0.93. These measurements represent total EGFR, not separate wild-type and mutant molecules.'))
html.append(fig('step3_protein_vs_mrna_g598v_cna.jpg','Figure 2. The same expression comparison with copy-number annotation. Red outlines denote amplification; black outlines indicate missing CNA profiling. Five of the six G598V samples have amplification calls and one lacks a CNA call.'))
html.append('<h3>Expression in mutant-bearing versus non-mutated samples</h3><p>The following statistics come from cBioPortal API queries for this cohort. Group membership uses recorded EGFR mutations; “no mutation” means no mutation recorded in this profile, not proof of a completely wild-type tumor.</p>')
rr=[]
for label,group in [('Any EGFR mutation',anymut),('No recorded EGFR mutation',common-anymut),('G598V',picked),('All samples without G598V',common-picked)]:
    ss=common&group;rr.append([label,len(ss),f(st.median(expr['protein_quantification'][s] for s in ss)),f'{st.median(expr["mrna"][s] for s in ss):,.1f}'])
html.append(table('EGFR expression by mutation group',['Group','n','Protein median','mRNA median (profile units)'],rr))
html.append(source('results/cbioportal_api/sample_table.tsv'))
html.append('''<p>The six G598V tumors have higher gene-level mRNA and protein values than the 82 samples with no recorded EGFR mutation. The 93-sample “without G598V” group also includes other EGFR mutations and should not be equated with the non-mutated group. Protein values are the portal abundance-ratio profile, not absolute concentrations. The raw mRNA profile is described by cBioPortal as UQ-normalized FPKM, median-centered by gene; its numeric scale must not be interpreted as TPM or compared directly with ICERFIRE expression inputs.</p>''')
rawr=st.correlation([expr['mrna'][s] for s in sorted(common)],[expr['protein_quantification'][s] for s in sorted(common)])
zr=st.correlation([expr['mrna_zscores'][s] for s in sorted(common)],[expr['protein_quantification'][s] for s in sorted(common)])
html.append(f'<p>Pearson correlation is {rawr:.3f} for the raw mRNA profile and {zr:.3f} for the log2-based mRNA z-score profile. The latter reproduces the displayed 0.92. Agreement across molecular levels supports a coherent association, but does not exclude shared technical or biological confounding.</p>')
rr=[]
for label,grp in [('Amplified, EGFR mutated',amp&anymut),('Amplified, no recorded EGFR mutation',amp-anymut),('Not amplified, no recorded EGFR mutation',set(expr['gistic'])-amp-anymut)]:
    ss=grp&common;rr.append([label,len(ss),f(st.median(expr['protein_quantification'][s] for s in ss))])
html.append(table('Copy-number stratification (GISTIC ≥ 2 defines amplification)',['Group','n','Protein median'],rr))
html.append('''<div class="summary">CNA is available for 96 samples; 47 are amplified. All 15 EGFR-mutated samples with CNA calls are amplified. The substantial expression contrast is confounded by amplification. These observational data do not establish that G598V increases expression or that amplification is the sole cause of the difference.</div>
<h3>Consequences for mass-spectrometry identification</h3>
<p>High parent-protein abundance can improve the opportunity to detect EGFR peptides, but this gene-level measurement does not quantify mutant protein, mutant allele dosage, or HLA-bound antigen. Tumor purity, subclonality, allele-specific expression, peptide ionization, and sampling depth all affect recovery.</p>
<p>For a complete tryptic digest with no missed cleavages, the site-containing peptide is <code>TCPAGVMGENNTLVWK</code> in WT and <code>TCPAVVMGENNTLVWK</code> in G598V, residues 594 to 609. Missed cleavages or other proteases can yield additional site-containing peptides. Identification requires a variant-aware sequence search and supporting fragment ions. A reference-only search cannot directly assign the mutant sequence.</p>
<p>The G→V substitution adds C<sub>3</sub>H<sub>6</sub>, approximately 42.047 Da, the same elemental mass increment as trimethylation. Acetylation is approximately 42.011 Da and can be distinguishable by sufficiently accurate precursor mass. Modification site chemistry and site-localizing fragment ions remain essential. The mass increments are tabulated in <a href="https://www.unimod.org/modifications_view.php?editid1=37">Unimod: trimethylation</a> and <a href="https://www.unimod.org/modifications_view.php?editid1=1">Unimod: acetylation</a>. Coexisting WT protein is possible, but its presence and abundance relative to mutant protein cannot be established from these gene-level profiles.</p>
<p class="limitation">Bulk tryptic proteomics and HLA immunopeptidomics answer different questions. Detecting a mutant tryptic EGFR peptide would support translation of the variant; demonstrating presentation requires HLA isolation and identification of the naturally bound peptide, usually without tryptic digestion of that ligand pool.</p></section>''')
html.append('''<section id="step4"><h2>Step 4. Structural Immunogenicity</h2>
<p>The input pairs are the nine mutation-spanning 9-mers from Step 2, with WT as peptide A and G598V as peptide B on the same input line. The global top-ten peptides are identical between proteins and therefore are not mutation-specific candidates. The Peptide Variant Comparison results include Class I pMHC Immunogenicity, Neo-Epitope Immunogenicity (ICERFIRE 1.0), and NetMHCpan 4.1 EL context for HLA-A*02:01.</p>''')
html.append(table('Paired peptide immunogenicity and ICERFIRE outputs',['Pair','Peptide A (WT)','Peptide B (G598V)','pMHC A','pMHC B','Δ pMHC','ICERFIRE prediction','ICERFIRE percentile'],[[r['seq #'],r['peptideA'],r['peptideB'],f(num(r,'peptideA immunogenicity score'),5),f(num(r,'peptideB immunogenicity score'),5),f"{num(r,'difference (peptideB - peptideA) immunogenicity score'):+.5f}",f(num(r,'icerfire prediction'),6),f(num(r,'icerfire percentile rank'),2)] for r in vc]))
html.append(source('results/iedb_variant_comparison/peptide_table.csv'))
html.append('<p>Higher pMHC scores and ICERFIRE predictions are favorable; lower percentile ranks are favorable. The export contains one ICERFIRE immunogenicity prediction per WT/mutant pair, not separate WT and mutant ICERFIRE immunogenicity scores. It therefore cannot directly quantify a WT-to-mutant immunogenicity increase for that model.</p>')
html.append(table('Presentation ranks must be distinguished from final ICERFIRE ranks',['Pair','NetMHCpan 4.1 EL A','NetMHCpan 4.1 EL B','ICERFIRE internal EL A','ICERFIRE internal EL B','Mutant ICORE'],[[r['seq #'],r['peptideA netmhcpan_el percentile'],r['peptideB netmhcpan_el percentile'],r['icerfire peptide_a_el_rank'],r['icerfire peptide_b_el_rank'],r['icerfire peptide_b_icore']] for r in vc]))
html.append(fig('step4_immunogenicity_score_scatter.jpg','Figure 3. Peptide B versus peptide A pMHC immunogenicity score. Three pairs lie exactly on the identity line; six lie slightly above it. The small differences are visible in the numeric table even where they are difficult to resolve in the screenshot.'))
html.append(fig('step4_icerfire_el_rank_scatter.jpg','Figure 4. Peptide B versus peptide A ICERFIRE internal EL rank, as identified by the screenshot metric. This is not the final ICERFIRE immunogenicity percentile. Six pairs improve, two worsen, and one is unchanged; pairs 3 and 4 overlap. The largest improvement is 80 to 10.1.'))
html.append(f'''<h3>Are the mutant peptides more immunogenic?</h3>
<p>The pMHC model predicts six small positive differences (0.00240 to 0.00744) and three exact ties. These are real model differences, not rounding artifacts, but there is no uncertainty estimate or validated minimum difference here that establishes biological importance. The ties occur at mutation positions P9, P2, and P1, which are excluded by the default mask. A masked residue contributes nothing to this calculation; this does not prove that it is inaccessible to every TCR or cannot indirectly change peptide conformation. See {cite(calis,'IEDB immunogenicity documentation')}.</p>
<p>ICERFIRE ranks TCPAVVMGE highest in this set (prediction 0.226528; percentile 15.86), followed by AVVMGENNT and CPAVVMGEN. No experimental response threshold was specified, so these ranks do not establish that all candidates are non-immunogenic. VVMGENNTL has a NetMHCpan 4.1 EL rank of 0.93 but an ICERFIRE percentile of 42.28. It has relatively favorable presentation prediction, not demonstrated mutant-specific T-cell recognition.</p>
<h3>Why can the models disagree?</h3>
<p>The pMHC method uses positional amino-acid features. ICERFIRE instead uses an ensemble model with a selected ICORE and presentation information. Its published consensus model uses anchor masking, BLOSUM mutation features, and expression. Thus, the difference cannot be explained by claiming that ICERFIRE has no mask. Pairs 3 and 4 select the same mutant ICORE, KTCPAVVM, and have identical exported ICERFIRE outputs. See {cite(icer,'Wan et al. (2024)')}.</p>
<p>The DTU service describes self-similarity among its features, but the paper explicitly says the final consensus model does not use the self-similarity feature. The exact deployed wrapper configuration cannot be determined from the exported results. A reported similarity column therefore does not demonstrate a fixed self-similarity penalty or explain a particular score causally. HCVKTCPAV has similarity 0.967956; VVMGENNTL is higher at 0.971274. Neither is evidence of measured tolerance. See {cite(dtu,'ICERFIRE service documentation')} and {cite(icer,'the consensus-model description')}.</p>
<p>All pairs report <code>total_gene_tpm = 6.071</code>. It is not a tumor-specific measurement in this analysis. The DTU service can obtain reference expression or accept user TPM, but this export does not document the origin of 6.071 in the IEDB wrapper. CPTAC FPKM-profile values cannot replace TPM directly; changing expression requires a rerun and need not improve every random-forest prediction monotonically. These scores are not established lower bounds.</p>
<p class="limitation">Neither tool reconstructs this peptide-HLA-TCR structure or measures a patient's T-cell repertoire. Presentation rank, affinity, immunogenicity score, and ICERFIRE percentile are distinct outputs. Differences in training data, features, and ICORE selection can change their candidate ordering.</p></section>''')
html.append('''<section id="step5"><h2>Step 5. LC-MS/MS Benchmark Validation and Future Experiments</h2>
<h3>TCIA comparison and precision calculation</h3>
<p>The Cancer Immunome Atlas export was filtered to disease GBM and gene EGFR through the Neoantigens table. It contains 111 rows, 78 distinct peptide sequences, and 35 patient identifiers, with peptide lengths from 8 to 11. Every HLA-alleles field is NA.</p>''')
html.append(f'<p class="limitation"><strong>This export does not provide an LC-MS/MS benchmark.</strong> TCIA neoantigen entries are computational candidates derived from genomic/transcriptomic data, not peptide-spectrum observations from an HLA pull-down. Consequently, an entry does not establish MS detection and absence does not establish a false positive. The distinction is supported by the {cite("https://pubmed.ncbi.nlm.nih.gov/28052254/","TCIA primary publication")} and the {cite("https://pmc.ncbi.nlm.nih.gov/articles/PMC11020248/","TSAFinder study’s description of the downloaded TCIA predictions")}.</p>')
rr=[]
for r in vc:
 p=r['peptideB'];patients=sorted({x['patientBarcode'] for x in tcia if x['peptide']==p});longer=sorted(q for q in tpeps if p in q and len(q)>len(p));rr.append([p,'Yes' if patients else 'No',len(patients),', '.join(longer) or 'None'])
html.append(table('Exact sequence and longer-peptide overlap with TCIA',['Step 4 mutant peptide','Exact match','Exact-match patients','Longer TCIA entries containing the 9-mer'],rr));html.append(source('results/tcia_gbm_egfr/neoantigens_gbm_egfr.tsv'))
html.append('''<p>If TCIA membership is provisionally treated as a reference label, exact matching gives nominal TP = 1 and nominal FP = 8, with TP/(TP + FP) = 1/9 = 11.1%. These should be reported as <strong>one matched and eight unmatched predictions</strong>, not experimentally established TP and FP. Allowing a longer peptide to contain the candidate gives 3/9 = 33.3% overlap, with nominal counts 3 and 6. This second metric is substring recovery, not exact peptide validation or mutation-site precision.</p>
<div class="summary">Observed MS TP and FP counts are unknown, so experimental positive predictive value cannot be estimated from these files. A 0% MS precision is also unjustified: unavailable validation labels are not negative results.</div>
<p>VVMGENNTL occurs in four patients: TCGA-06-0174, TCGA-12-0616, TCGA-19-2620, and TCGA-28-5213. Six distinct TCIA peptide sequences are compatible with the local G598V sequence, including GPHCVKTCPAV. These entries span seven patients in the filtered EGFR table. The export has no genomic variant column or HLA assignment, so sequence compatibility alone does not verify each patient's mutation or HLA restriction and is not a population prevalence estimate.</p>
<h3>Does the result justify experimental testing?</h3>
<p>The overlap fraction cannot estimate the chance of experimental success or justify rejecting the eight unmatched peptides. VVMGENNTL is a reasonable initial presentation candidate because it leads the mutation-containing 9-mers by processing total, has a relatively favorable EL rank, and has an exact TCIA sequence match. However, its WT counterpart has a slightly better NetMHCpan EL rank (0.78 versus 0.93), so the mutation did not create the favorable EL prediction. Both sequences need experimental comparison. HCVKTCPAV provides a useful contrasting anchor-change candidate, while TCPAVVMGE illustrates the different ICERFIRE ordering. None is established as an immunogenic neoantigen.</p>
<p>All nine candidates overlap one substitution. Their outcomes are dependent, and TCIA uses a different cohort and selection pipeline. Sequence-only agreement does not constitute an independent experimental replication. Choice of a denominator, lengths, HLA alleles, and positivity threshold must be fixed before a meaningful benchmark is computed.</p>
<h3>Proposed experimental sequence</h3>
<ol>
<li><strong>Establish endogenous presentation.</strong> Use an HLA-A*02:01-positive, G598V-positive glioma model with sequence-confirmed controls. Isolate HLA-bound ligands and perform LC-MS/MS against a variant-aware database with target-decoy error control. Compare candidate fragmentation and retention time with synthetic isotope-labeled standards and use targeted acquisition where appropriate. Include matched WT, mutant-negative, and HLA-negative or HLA-blocked controls. Confirming bulk EGFR expression is insufficient.</li>
<li><strong>Measure mutant-specific T-cell recognition.</strong> Compare mutant and WT peptide titrations using HLA-matched donor or patient T cells, peptide-HLA multimers, IFN-γ ELISpot, and intracellular cytokine assays. Require recognition and killing of cells expressing the full-length mutant antigen without exogenous peptide loading. Test WT cross-reactivity and HLA dependence; peptide-pulsed target recognition alone bypasses natural processing.</li>
<li><strong>Test immunogenicity and tumor control in vivo.</strong> An HLA-A2 transgenic system such as HHD can test induction of relevant CD8 responses. Use a compatible syngeneic orthotopic glioma engineered with the matching HLA construct and full-length antigen, with isogenic WT and mutant tumor controls. Compare mutant vaccination, WT vaccination, and adjuvant-only groups; predefine randomization, blinded tumor assessment, tumor burden, survival, and immune readouts. T-cell depletion or HLA-loss controls can test response dependence. This is a proposed model, not an experiment already performed.</li>
</ol>''')
html.append(f'''<p>The original HHD system uses an HLA-A2.1/H-2D<sup>b</sup> monochain on an H-2D<sup>b</sup>/β2-microglobulin double-knockout background ({cite('https://pmc.ncbi.nlm.nih.gov/articles/PMC2196346/','Pascolo et al., 1997')}). A human glioma xenograft is not interchangeable with a syngeneic immunocompetent model. Human EGFR and murine tolerance differences can exaggerate apparent specificity, so human HLA-matched assays remain necessary.</p>
<h3>Glioblastoma literature and translational limits</h3>
<p>The phase 3 ACT IV trial of EGFRvIII-targeted rindopepimut did not improve survival. This supports caution about extrapolating antigen-specific responses to clinical benefit; it does not prove that every single-antigen vaccine must fail or that antigen loss alone caused the result ({cite('https://doi.org/10.1016/S1470-2045(17)30517-X','Weller et al., 2017')}).</p>
<p>Recent GBM studies support investigating multiple personalized targets: a phase 1 DNA vaccine study included up to 40 neoantigens and reported vaccine-associated T-cell responses ({cite('https://www.nature.com/articles/s43018-026-01163-w','Garfinkle et al., 2026')}). A single-arm phase Ib neoantigen-pulsed dendritic-cell trial treated 11 patients and reported median PFS of 16.2 months from surgery ({cite('https://www.nature.com/articles/s41467-026-75066-w','Zhang et al., 2026')}). These early studies do not establish randomized survival benefit or validate G598V. A confirmed G598V epitope might be considered in a broader antigen strategy, subject to HLA matching, tumor heterogeneity, natural presentation, and mutant-specific recognition.</p></section>
<footer>Data provenance: supplied FASTA files, archived IEDB and TCIA result tables, retained screenshots, and an additional cBioPortal API snapshot in <code>results/cbioportal_api/</code>. Tables are generated directly from these files by <code>scripts/04_build_assignment.py</code>. Screenshots are embedded in this document; their source files are in figures/.</footer></main></body></html>''')
out='\n'.join(html)
assert '\u2014' not in out
(ROOT/'docs/assignment.html').write_text(out)
print('Wrote docs/assignment.html',len(out),'characters')
