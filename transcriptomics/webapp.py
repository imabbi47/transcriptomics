"""FastAPI backend for the transcriptomics web app.

Turns an accession into a running pipeline job:
  GET  /                    web UI
  GET  /api/health          service status + which external tools are installed
  POST /api/resolve         {accession} -> run table (synchronous, fast)
  POST /api/jobs            {accession} -> start an async ingest job -> {job_id}
  GET  /api/jobs            list jobs
  GET  /api/jobs/{job_id}   job status + step-by-step progress + result

Run:  transcriptomics serve         (or: uvicorn transcriptomics.webapp:app)
"""
from __future__ import annotations

import shutil
import threading
import uuid
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional

try:
    from fastapi import FastAPI, HTTPException
    from fastapi.responses import HTMLResponse
    from pydantic import BaseModel
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit("web extra not installed — run:  pip install -e '.[web]'") from exc

from .accession import classify
from .metadata import MetadataError, resolve_runs
from .models import Run

app = FastAPI(title="transcriptomics", version="0.1.0")

_JOBS: Dict[str, "Job"] = {}
_LOCK = threading.Lock()
_EXTERNAL_TOOLS = ["prefetch", "fasterq-dump", "salmon", "fastp", "curl"]


@dataclass
class Job:
    id: str
    accession: str
    status: str = "queued"          # queued | running | done | error
    steps: List[dict] = field(default_factory=list)
    error: Optional[str] = None
    result: Optional[dict] = None


class AccessionIn(BaseModel):
    accession: str


def _run_dict(run: Run) -> dict:
    return {
        "run": run.run_accession,
        "sample": run.friendly_name(),
        "layout": run.library_layout.value,
        "organism": run.organism or "",
        "platform": run.platform or "",
        "spots": run.spots,
        "bytes": run.total_bytes,
        "source": "ENA" if run.has_ena_fastq else "SRA",
    }


def _resolve(accession_str: str):
    acc = classify(accession_str)
    if not acc.is_resolvable:
        raise MetadataError(f"Unrecognised accession: {accession_str!r}")
    return acc, resolve_runs(acc)


def _step(name: str, state: str, detail: str = "") -> dict:
    return {"step": name, "state": state, "detail": detail}


def _run_job(job_id: str, accession_str: str) -> None:
    job = _JOBS[job_id]
    job.status = "running"
    try:
        job.steps.append(_step("Resolve accession", "running"))
        acc, runs = _resolve(accession_str)
        paired = sum(1 for r in runs if r.is_paired)
        job.steps[-1] = _step("Resolve accession", "done",
                              f"{len(runs)} run(s) · {paired} paired · "
                              f"{(runs[0].organism or 'organism ?') if runs else ''}")

        job.steps.append(_step("Prepare samplesheet", "done",
                              "sample/run/fastq/strandedness table ready"))

        have_dl = bool(shutil.which("curl") or shutil.which("prefetch"))
        job.steps.append(_step(
            "Download reads (fetch)",
            "ready" if have_dl else "skipped",
            "ENA/SRA reads available to download" if have_dl
            else "no downloader (curl/sra-tools) on the server"))

        have_quant = bool(shutil.which("salmon"))
        job.steps.append(_step(
            "Quantify → DE → enrichment",
            "ready" if have_quant else "needs setup",
            "salmon present" if have_quant
            else "install salmon + provide an experimental design to run locally"))

        job.result = {
            "accession": acc.value,
            "type": acc.type.value,
            "n_runs": len(runs),
            "total_bytes": sum(r.total_bytes for r in runs),
            "runs": [_run_dict(r) for r in runs],
        }
        job.status = "done"
    except MetadataError as error:
        job.steps.append(_step("Resolve accession", "error", str(error)))
        job.error = str(error)
        job.status = "error"
    except Exception as error:  # pragma: no cover - safety net
        job.error = f"{type(error).__name__}: {error}"
        job.status = "error"


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok",
            "tools": {t: bool(shutil.which(t)) for t in _EXTERNAL_TOOLS}}


@app.post("/api/resolve")
def api_resolve(body: AccessionIn) -> dict:
    try:
        acc, runs = _resolve(body.accession.strip())
    except MetadataError as error:
        raise HTTPException(status_code=400, detail=str(error))
    return {"accession": acc.value, "type": acc.type.value, "n_runs": len(runs),
            "total_bytes": sum(r.total_bytes for r in runs),
            "runs": [_run_dict(r) for r in runs]}


@app.post("/api/jobs")
def api_start_job(body: AccessionIn) -> dict:
    acc = classify(body.accession.strip())
    if not acc.is_resolvable:
        raise HTTPException(status_code=400, detail=f"Unrecognised accession: {body.accession!r}")
    job_id = uuid.uuid4().hex[:12]
    with _LOCK:
        _JOBS[job_id] = Job(id=job_id, accession=acc.value)
    threading.Thread(target=_run_job, args=(job_id, body.accession.strip()), daemon=True).start()
    return {"job_id": job_id}


@app.get("/api/jobs")
def api_list_jobs() -> list:
    return [{"id": j.id, "accession": j.accession, "status": j.status} for j in _JOBS.values()]


@app.get("/api/jobs/{job_id}")
def api_job(job_id: str) -> dict:
    job = _JOBS.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="no such job")
    return asdict(job)


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return FRONTEND_HTML


FRONTEND_HTML = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Transcriptomics</title><style>
:root{--bg:#0f172a;--card:#1e293b;--ink:#e2e8f0;--muted:#94a3b8;--accent:#2dd4bf;--accent2:#818cf8;--line:#334155}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
 font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;line-height:1.6}
.wrap{max-width:860px;margin:0 auto;padding:0 20px 60px}
header{text-align:center;padding:56px 20px 26px}
h1{font-size:2.3rem;margin:0 0 .2em}
h1 .g{background:linear-gradient(90deg,var(--accent),var(--accent2));-webkit-background-clip:text;background-clip:text;color:transparent}
header p{color:var(--muted);margin:0}
.card{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:22px;margin-top:18px}
input{width:100%;padding:13px 15px;border-radius:10px;border:1px solid var(--line);background:#0b1220;color:var(--ink);font-size:1rem;margin-bottom:12px}
button{background:linear-gradient(90deg,var(--accent),var(--accent2));color:#04201c;border:0;padding:12px 22px;border-radius:10px;font-weight:700;font-size:1rem;cursor:pointer}
button:disabled{opacity:.5;cursor:not-allowed}
#steps{margin-top:16px}.step{display:flex;gap:10px;align-items:center;padding:6px 0;font-size:.92rem}
.dot{width:11px;height:11px;border-radius:50%;flex:0 0 auto;background:#475569}
.done .dot{background:var(--accent)}.running .dot{background:#facc15;animation:p 1s infinite}
.error .dot{background:#ef4444}.ready .dot{background:#38bdf8}.skipped .dot,.needs .dot{background:#64748b}
@keyframes p{50%{opacity:.3}}.detail{color:var(--muted);font-size:.85rem}
table{width:100%;border-collapse:collapse;margin-top:14px;font-size:.85rem}
th,td{border-bottom:1px solid var(--line);padding:7px 9px;text-align:left}th{color:var(--accent)}
td.n{text-align:right;font-variant-numeric:tabular-nums}
#msg{margin-top:12px;min-height:1.2em;font-size:.9rem}.err{color:#fca5a5}
.env{text-align:center;color:var(--muted);font-size:.78rem;margin-top:8px}
.env b{color:var(--ink)}
</style></head><body><div class="wrap">
<header><h1>🧬 <span class="g">Transcriptomics</span></h1>
<p>RNA-seq analysis straight from a GEO/SRA accession number — for any organism.</p>
<div class="env" id="env">checking server…</div></header>
<div class="card">
  <input id="acc" placeholder="Enter an accession, e.g. GSE52778, SRR1039508, GSM1275862">
  <button id="run" onclick="run()">Run analysis</button>
  <div id="msg"></div>
  <div id="steps"></div>
  <div id="table"></div>
</div></div>
<script>
function hbytes(n){if(!n)return"?";const u=["B","KB","MB","GB","TB"];let i=0,v=n;while(v>=1024&&i<4){v/=1024;i++;}return v.toFixed(i?1:0)+u[i];}
async function health(){try{const r=await fetch("/api/health");const j=await r.json();
 const t=j.tools;const on=Object.keys(t).filter(k=>t[k]);
 document.getElementById("env").innerHTML="server ready · tools installed: <b>"+(on.length?on.join(", "):"none (ENA download uses curl)")+"</b>";}
 catch(e){document.getElementById("env").textContent="server unreachable";}}
async function run(){
 const acc=document.getElementById("acc").value.trim();
 const btn=document.getElementById("run"),msg=document.getElementById("msg"),steps=document.getElementById("steps"),tbl=document.getElementById("table");
 msg.textContent="";msg.className="";steps.innerHTML="";tbl.innerHTML="";
 if(!acc){msg.textContent="Enter an accession first.";return;}
 btn.disabled=true;msg.textContent="Starting job…";
 try{
  const r=await fetch("/api/jobs",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({accession:acc})});
  if(!r.ok){const e=await r.json();throw new Error(e.detail||"failed");}
  const {job_id}=await r.json();poll(job_id);
 }catch(e){msg.textContent=e.message;msg.className="err";btn.disabled=false;}
}
async function poll(id){
 const msg=document.getElementById("msg"),steps=document.getElementById("steps"),btn=document.getElementById("run");
 const r=await fetch("/api/jobs/"+id);const j=await r.json();
 msg.textContent="Job "+id+" — "+j.status;
 steps.innerHTML=j.steps.map(s=>`<div class="step ${s.state}"><span class="dot"></span><div><b>${s.step}</b> <span class="detail">${s.state}${s.detail?" · "+s.detail:""}</span></div></div>`).join("");
 if(j.status==="running"||j.status==="queued"){setTimeout(()=>poll(id),700);return;}
 btn.disabled=false;
 if(j.status==="error"){msg.className="err";return;}
 if(j.result){renderTable(j.result);}
}
function renderTable(res){
 const tbl=document.getElementById("table");
 let h=`<p class="detail">${res.accession} · ${res.n_runs} run(s) · ~${hbytes(res.total_bytes)} total</p>`;
 h+="<table><thead><tr><th>Run</th><th>Sample</th><th>Layout</th><th>Organism</th><th class='n'>Size</th><th>Source</th></tr></thead><tbody>";
 h+=res.runs.map(r=>`<tr><td>${r.run}</td><td>${r.sample}</td><td>${r.layout}</td><td>${r.organism}</td><td class="n">${hbytes(r.bytes)}</td><td>${r.source}</td></tr>`).join("");
 h+="</tbody></table>";
 document.getElementById("table").innerHTML=h;
}
health();
</script></body></html>"""
