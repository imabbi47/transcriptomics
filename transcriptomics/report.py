#!/usr/bin/env python3
"""Bundle DE + enrichment outputs into one self-contained HTML report.

Reads an output directory produced by de.py / enrich.py (de_results.csv, the
PNG plots, enrichr_*.csv) and writes a single portable report.html with every
image embedded as a data URI. Stdlib only.

    python report.py --outdir results --title "GSE… — A vs B" --out report.html
"""
from __future__ import annotations

import argparse
import base64
import csv
import html
import os


def data_uri(path):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode("ascii")


def read_rows(path):
    if not os.path.exists(path):
        return [], []
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    return (rows[0], rows[1:]) if rows else ([], [])


def col(header, name):
    try:
        return header.index(name)
    except ValueError:
        return None


def fmt_num(v, kind):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return html.escape(str(v))
    if kind == "int":
        return f"{x:,.0f}"
    if kind == "lfc":
        return f"{x:+.2f}"
    if kind == "sci":
        return f"{x:.1e}"
    return f"{x:g}"


def de_summary(header, rows, alpha):
    ip, il = col(header, "padj"), col(header, "log2FoldChange")
    tested = up = down = 0
    for r in rows:
        try:
            padj = float(r[ip])
        except (ValueError, IndexError, TypeError):
            continue
        tested += 1
        if padj < alpha:
            lfc = float(r[il])
            up += lfc > 0
            down += lfc < 0
    return tested, up, down


def de_table(header, rows, n=20):
    isym, ibm = col(header, "symbol"), col(header, "baseMean")
    ilfc, ipadj = col(header, "log2FoldChange"), col(header, "padj")
    out = []
    for r in rows[:n]:
        out.append((
            html.escape(r[isym]) if isym is not None else "",
            fmt_num(r[ibm], "int"),
            fmt_num(r[ilfc], "lfc"),
            fmt_num(r[ipadj], "sci"),
        ))
    return out


def enr_table(path, n=10):
    header, rows = read_rows(path)
    if not header:
        return []
    igs, it = col(header, "Gene_set"), col(header, "Term")
    iov, iap = col(header, "Overlap"), col(header, "Adjusted P-value")
    out = []
    for r in rows[:n]:
        out.append((
            html.escape(r[igs]) if igs is not None else "",
            html.escape(r[it]) if it is not None else "",
            fmt_num(r[iap], "sci") if iap is not None else "",
            html.escape(r[iov]) if iov is not None else "",
        ))
    return out


CSS = """
*{box-sizing:border-box}
body{margin:0}
.rep{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
  color:#1b2430;background:#f7f8f6;line-height:1.55;
  max-width:960px;margin:0 auto;padding:0 22px 60px}
.rep h1{font-size:1.7rem;margin:.2em 0 .1em;letter-spacing:-.01em}
.rep .sub{color:#5b6b7a;margin:0 0 4px}
.rep .acc{font-family:ui-monospace,Menlo,Consolas,monospace;color:#0e7c73;font-weight:600}
.rep h2{font-size:1.15rem;margin:34px 0 12px;padding-bottom:6px;border-bottom:2px solid #0e7c73;color:#0e5a54}
.rep header{padding:34px 0 6px}
.cards{display:flex;flex-wrap:wrap;gap:12px;margin:14px 0 6px}
.card{background:#fff;border:1px solid #e2e6e1;border-radius:10px;padding:14px 18px;min-width:130px;flex:1}
.card .n{font-size:1.6rem;font-weight:700;font-variant-numeric:tabular-nums}
.card .l{font-size:.8rem;color:#5b6b7a;text-transform:uppercase;letter-spacing:.04em}
.card.up .n{color:#c0392b}.card.down .n{color:#2c6fb0}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px}
@media(max-width:680px){.grid2{grid-template-columns:1fr}}
.rep img{max-width:100%;border:1px solid #e2e6e1;border-radius:10px;background:#fff}
figure{margin:0}
figcaption{font-size:.82rem;color:#5b6b7a;margin-top:5px;text-align:center}
table{border-collapse:collapse;width:100%;font-size:.85rem;background:#fff;
  border:1px solid #e2e6e1;border-radius:10px;overflow:hidden}
th,td{padding:7px 11px;text-align:left;border-bottom:1px solid #eef1ee}
th{background:#eef4f3;color:#0e5a54;font-weight:600}
td.num{text-align:right;font-variant-numeric:tabular-nums;font-family:ui-monospace,Menlo,Consolas,monospace}
tr:last-child td{border-bottom:none}
.wrap-x{overflow-x:auto}
.rep footer{margin-top:40px;padding-top:14px;border-top:1px solid #e2e6e1;color:#5b6b7a;font-size:.8rem}
.up-h{color:#c0392b}.down-h{color:#2c6fb0}
"""


def img_fig(uri, caption):
    if not uri:
        return f'<figcaption>({caption} not available)</figcaption>'
    return f'<figure><img src="{uri}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>'


def table_html(rows, headers, numcols):
    thead = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
    body = []
    for r in rows:
        cells = "".join(
            f'<td class="num">{c}</td>' if i in numcols else f"<td>{c}</td>"
            for i, c in enumerate(r)
        )
        body.append(f"<tr>{cells}</tr>")
    return f'<div class="wrap-x"><table><thead><tr>{thead}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--accession", default="")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    header, rows = read_rows(os.path.join(args.outdir, "de_results.csv"))
    tested, up, down = de_summary(header, rows, args.alpha)
    top = de_table(header, rows)
    up_paths = enr_table(os.path.join(args.outdir, "enrichr_up.csv"))
    down_paths = enr_table(os.path.join(args.outdir, "enrichr_down.csv"))

    pca = data_uri(os.path.join(args.outdir, "pca.png"))
    volcano = data_uri(os.path.join(args.outdir, "volcano.png"))
    ma = data_uri(os.path.join(args.outdir, "ma.png"))
    enr_up_img = data_uri(os.path.join(args.outdir, "enrichment_up.png"))
    enr_down_img = data_uri(os.path.join(args.outdir, "enrichment_down.png"))

    acc = f'<span class="acc">{html.escape(args.accession)}</span> · ' if args.accession else ""
    parts = [
        f"<title>{html.escape(args.title)}</title>",
        f"<style>{CSS}</style>",
        '<div class="rep">',
        "<header>",
        f"<h1>{html.escape(args.title)}</h1>",
        f'<p class="sub">{acc}{html.escape(args.subtitle)}</p>',
        "</header>",
        '<div class="cards">',
        f'<div class="card"><div class="n">{up + down:,}</div><div class="l">Significant (padj&lt;{args.alpha})</div></div>',
        f'<div class="card up"><div class="n">{up:,}</div><div class="l">Up-regulated</div></div>',
        f'<div class="card down"><div class="n">{down:,}</div><div class="l">Down-regulated</div></div>',
        f'<div class="card"><div class="n">{tested:,}</div><div class="l">Genes tested</div></div>',
        "</div>",
        "<h2>Sample overview</h2>",
        img_fig(pca, "PCA of top variable genes"),
        "<h2>Differential expression</h2>",
        f'<div class="grid2">{img_fig(volcano, "Volcano")}{img_fig(ma, "MA plot")}</div>',
        "<h3>Top 20 genes by adjusted p-value</h3>",
        table_html(top, ["Gene", "Base mean", "log2 FC", "Adj. p"], {1, 2, 3}),
        "<h2>Pathway enrichment</h2>",
        '<h3 class="up-h">Up-regulated pathways</h3>',
        table_html(up_paths, ["Library", "Term", "Adj. p", "Overlap"], {2}) if up_paths else "<p>(none)</p>",
        img_fig(enr_up_img, "Enriched pathways — up-regulated"),
        '<h3 class="down-h">Down-regulated pathways</h3>',
        table_html(down_paths, ["Library", "Term", "Adj. p", "Overlap"], {2}) if down_paths else "<p>(none)</p>",
        img_fig(enr_down_img, "Enriched pathways — down-regulated"),
        '<footer>Generated by the <b>transcriptomics</b> pipeline · DESeq2 (pyDESeq2) + Enrichr (gseapy). '
        'Counts → differential expression → pathway enrichment.</footer>',
        "</div>",
    ]
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(parts))
    print(f"[report] wrote {args.out}  ({os.path.getsize(args.out)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
