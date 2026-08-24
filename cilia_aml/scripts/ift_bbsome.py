#!/usr/bin/env python3
"""Report the IFT / BBSome machinery in the AML miR-551b analysis, gene by gene."""
import pandas as pd

FAM = {
    "IFT-A (retrograde complex)": ["IFT122", "IFT140", "IFT43", "WDR19", "WDR35", "TTC21B"],
    "IFT-B (anterograde complex)": ["IFT20", "IFT22", "IFT27", "IFT46", "IFT52", "IFT57",
                                    "IFT74", "IFT80", "IFT81", "IFT88", "IFT172", "HSPB11",
                                    "RABL4", "TRAF3IP1", "TTC26", "TTC30A", "TTC30B", "CLUAP1"],
    "BBSome + chaperonins": ["BBS1", "BBS2", "BBS4", "BBS5", "BBS7", "BBS9", "BBS10", "BBS12",
                             "BBIP1", "ARL6", "MKKS", "TTC8", "CCDC28B"],
    "Dynein-2 motor (retrograde)": ["DYNC2H1", "DYNC2LI1", "WDR34", "WDR60", "TCTEX1D2"],
    "Kinesin-2 motor (anterograde)": ["KIF3A", "KIF3B", "KIF3C", "KIFAP3"],
}

cil = {l.strip().upper() for l in open("cilia_aml/data/cilia_genes.txt") if l.strip()}
de = pd.read_csv("cilia_aml/de_paired/de_results.csv", index_col=0)
de["sym"] = de["symbol"].astype(str).str.upper()
d = de.drop_duplicates("sym").set_index("sym")

rows = []
print("=" * 80)
print("IFT / BBSome MACHINERY IN AML   (anti-miR-551b vs Control, paired model, n=4 vs 4)")
print("=" * 80)
for fam, genes in FAM.items():
    print("\n### " + fam)
    print("    {:<10} {:<6} {:>9} {:>8} {:>7} {:>6}   {}".format(
        "gene", "inSet", "baseMean", "log2FC", "nom p", "padj", "note"))
    for g in genes:
        inset = "yes" if g in cil else "-"
        if g in d.index:
            r = d.loc[g]
            bm, lfc = float(r.baseMean), float(r.log2FoldChange)
            p = float(r.pvalue) if pd.notna(r.pvalue) else float("nan")
            q = float(r.padj) if pd.notna(r.padj) else float("nan")
            note = "very low expr" if bm < 10 else ("low expr" if bm < 50 else "")
            if pd.notna(q) and q < 0.05:
                note = (note + " SIGNIFICANT").strip()
            print("    {:<10} {:<6} {:>9.0f} {:>+8.2f} {:>7.3f} {:>6.2f}   {}".format(
                g, inset, bm, lfc, p, q, note))
            rows.append(dict(family=fam, gene=g, in_683_set=g in cil, baseMean=round(bm, 1),
                             log2FC=round(lfc, 3), pvalue=p, padj=q))
        else:
            print("    {:<10} {:<6} {:>9} {:>8} {:>7} {:>6}   NOT DETECTED".format(
                g, inset, "--", "--", "--", "--"))
            rows.append(dict(family=fam, gene=g, in_683_set=g in cil, baseMean=None,
                             log2FC=None, pvalue=None, padj=None))

t = pd.DataFrame(rows)
det = t[t.baseMean.notna()]
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print("  genes queried .......................... {}".format(len(t)))
print("  detected / expressed ................... {}".format(len(det)))
print("  NOT detected ........................... {}".format(int(t.baseMean.isna().sum())))
print("  significant (FDR<0.05) ................. {}".format(int((det.padj < 0.05).sum())))
print("  nominally p<0.05 ....................... {}".format(int((det.pvalue < 0.05).sum())))
print("  median baseMean ........................ {:.0f}".format(det.baseMean.median()))
print("  mean |log2FC|, IFT/BBSome .............. {:.3f}".format(det.log2FC.abs().mean()))
print("  mean |log2FC|, whole transcriptome ..... {:.3f}".format(de.log2FoldChange.abs().mean()))
imax = det.log2FC.abs().idxmax()
print("  largest shift .......................... {} ({:+.2f}, padj={:.2f})".format(
    det.loc[imax, "gene"], det.loc[imax, "log2FC"], det.loc[imax, "padj"]))
print("  members of the 683-gene set ............ {}/{}".format(int(t.in_683_set.sum()), len(t)))

# expression strata
print("\n  expression distribution of detected members:")
for lo, hi, lab in [(1000, 1e12, "high   (>1000)"), (100, 1000, "medium (100-1000)"),
                    (10, 100, "low    (10-100)"), (0, 10, "v.low  (<10)")]:
    n = int(((det.baseMean >= lo) & (det.baseMean < hi)).sum())
    print("    {:<20} {:>3} genes".format(lab, n))

t.to_csv("cilia_aml/de_paired/IFT_BBSome_in_AML.csv", index=False)
print("\n  wrote cilia_aml/de_paired/IFT_BBSome_in_AML.csv")
try:
    t.to_excel("cilia_aml/de_paired/IFT_BBSome_in_AML.xlsx", index=False)
    print("  wrote cilia_aml/de_paired/IFT_BBSome_in_AML.xlsx")
except Exception as exc:
    print("  xlsx skipped: {!r}".format(exc))
