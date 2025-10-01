#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  combine_reports.py — TrackTx Cohort SPA (polished, single-file, offline)║
# ║                                                                          ║
# ║  Inputs  : per-sample *.summary.json/*.report.json (also *.json.gz)      ║
# ║            paths may be files or directories (expanded upstream & here)  ║
# ║  Outputs : global_summary.{html,tsv,json} + optional region totals TSV   ║
# ║  Design  : No external deps/CDNs; embedded CSS/JS; exportable charts     ║
# ║  Science : Surfaces unlocalized fraction; PI & density medians; QC meds  ║
# ╚══════════════════════════════════════════════════════════════════════════╝

from __future__ import annotations
import argparse, json, math, os, sys, gzip, io, re, datetime
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd

# ─────────────────────────── CLI ───────────────────────────
ap = argparse.ArgumentParser(description="Combine TrackTx per-sample summaries into a cohort report.")
ap.add_argument("--inputs", nargs="+", required=True, help="Per-sample JSONs (files or directories). Supports *.json, *.json.gz.")
ap.add_argument("--out-tsv", required=True)
ap.add_argument("--out-json", required=True)
ap.add_argument("--out-html", required=True)
ap.add_argument("--out-regions", default=None, help="Optional TSV of summed functional-region totals.")
ap.add_argument("--pipeline-version", default="unknown", help="Pipeline version for metadata.")
ap.add_argument("--run-name", default="unnamed", help="Run name for metadata.")
ap.add_argument("--duration", default="unknown", help="Run duration for metadata.")
ap.add_argument("--profile", default="unknown", help="Execution profile for metadata.")
args = ap.parse_args()
print(f"[combine_py] start ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)

# ───────────────────────── Utilities ───────────────────────
def _walk_inputs(paths: List[str]) -> List[str]:
    out: List[str] = []
    for p in paths:
        if not p: continue
        p = str(p)
        if os.path.isdir(p):
            for root, _, files in os.walk(p):
                for fn in files:
                    lf = fn.lower()
                    if lf.endswith((".summary.json",".report.json",".json",".json.gz")):
                        out.append(os.path.join(root, fn))
        else:
            lf = p.lower()
            if lf.endswith((".summary.json",".report.json",".json",".json.gz")):
                out.append(p)
    # stable order
    return sorted(set(out))

def _read_json_any(path: str) -> Optional[Dict[str, Any]]:
    try:
        if not os.path.exists(path): return None
        if path.lower().endswith(".gz"):
            with gzip.open(path, "rb") as fh:
                return json.loads(fh.read().decode("utf-8", "replace"))
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception as e:
        sys.stderr.write(f"WARN: failed to read {path}: {e}\n")
        return None

def _is_num(x: Any) -> bool:
    return isinstance(x, (int, float)) and not (isinstance(x, float) and (math.isnan(x) or math.isinf(x)))

def _to_float(x: Any) -> Optional[float]:
    try:
        v = float(x)
        return None if math.isnan(v) or math.isinf(v) else v
    except Exception:
        return None

def _sanitize(obj: Any) -> Any:
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)): return None
    if isinstance(obj, dict):  return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, list):  return [_sanitize(v) for v in obj]
    return obj

# ───────────── Schema normalization & science helpers ─────────────
_UNLOC_RE = re.compile(r"non[-\s_]?localized|unlocalized|unlocalised", re.I)

def _normalize_sample(js: Dict[str, Any], fallback_name: str) -> Optional[Dict[str, Any]]:
    # Accept both new/legacy keys; guard missing content
    sample_id = str(js.get("sample") or js.get("sample_id") or fallback_name or "NA")
    cond      = js.get("condition", "NA")
    tp        = js.get("timepoint", "NA")
    rep       = js.get("replicate", "NA")

    # Regions list: [{"region": name, "reads": float, "region_count": int}, ...]
    regions = js.get("regions") or []
    func_totals: Dict[str, float] = {}
    region_counts: Dict[str, int] = {}
    for it in regions:
        try:
            k = str(it.get("region"))
            v = float(it.get("reads", 0) or 0.0)
            c = int(it.get("region_count", 0) or 0)
            if k:
                func_totals[k] = func_totals.get(k, 0.0) + v
                region_counts[k] = region_counts.get(k, 0) + c
        except Exception:
            pass
    reads_total_functional = sum(func_totals.values()) if func_totals else _to_float(js.get("reads_total_functional")) or 0.0

    # Unlocalized fraction (computed if not present): sum of “non-localized” labels / total
    unloc = js.get("unlocalized_fraction")
    if unloc is None and func_totals and reads_total_functional:
        unl_reads = sum(v for k, v in func_totals.items() if _UNLOC_RE.search(str(k) or ""))
        unloc = unl_reads / reads_total_functional if reads_total_functional > 0 else None
    unloc = _to_float(unloc)

    # Extract nested fields from metrics and qc objects
    metrics = js.get("metrics") or {}
    qc = js.get("qc") or {}
    
    # Determine duplicate percentage - prefer UMI deduplication if available
    dup_percent = None
    if qc.get("umi_deduplication_enabled", False) and qc.get("umi_deduplication_percent") is not None:
        dup_percent = qc.get("umi_deduplication_percent")
    else:
        dup_percent = qc.get("duplicate_percent") or js.get("duplicate_percent")
    
    rec = dict(
        sample_id=sample_id,
        condition=cond, timepoint=tp, replicate=rep,
        divergent_regions=metrics.get("divergent_regions"),
        total_regions=metrics.get("total_functional_regions") or metrics.get("total_regions"),
        reads_total_functional=reads_total_functional,
        median_pausing_index=metrics.get("median_pausing_index") or metrics.get("median_pi"),
        median_density=metrics.get("median_functional_cpm") or metrics.get("median_density"),
        cpm_factor=metrics.get("cpm_factor"),
        crpmsi_factor=metrics.get("sicpm_factor") or metrics.get("crpmsi_factor"),
        input_reads=qc.get("total_reads_raw") or js.get("input_reads"),
        dedup_reads=qc.get("dedup_reads_mapq_ge") or js.get("dedup_reads"),
        duplicate_percent=dup_percent,
        umi_deduplication_enabled=qc.get("umi_deduplication_enabled", False),
        umi_deduplication_percent=qc.get("umi_deduplication_percent"),
        unlocalized_fraction=unloc,
        func_totals=func_totals,
        region_counts=region_counts,
        regions=regions
    )
    return rec

def _collect_region_keys(samples: List[Dict[str, Any]]) -> List[str]:
    keys = {k for s in samples for k in (s.get("func_totals") or {}).keys()}
    return sorted(keys)

def _sum_region_totals(samples: List[Dict[str, Any]]) -> Dict[str, float]:
    agg: Dict[str, float] = {}
    for s in samples:
        for k, v in (s.get("func_totals") or {}).items():
            try: agg[k] = agg.get(k, 0.0) + float(v or 0.0)
            except Exception: pass
    return agg

# ─────────────────────── Load & normalize ───────────────────────
in_paths = _walk_inputs(args.inputs)
if not in_paths:
    sys.stderr.write("ERROR: no usable input paths.\n")
    sys.exit(2)

samples: List[Dict[str, Any]] = []
skipped: List[str] = []
for p in in_paths:
    js = _read_json_any(p)
    if not js:
        skipped.append(p); continue
    rec = _normalize_sample(js, os.path.basename(p).split(".")[0])
    if rec:
        samples.append(rec)
    else:
        skipped.append(p)

if not samples:
    sys.stderr.write("ERROR: no valid sample JSONs loaded.\n")
    sys.exit(3)

region_keys = _collect_region_keys(samples)
region_totals = _sum_region_totals(samples)

# ─────────────────────── Build TSV rows ───────────────────────
core_cols = [
    "sample_id","condition","timepoint","replicate",
    "input_reads","dedup_reads","duplicate_percent",
    "divergent_regions","total_regions","reads_total_functional",
    "median_pausing_index","median_density",
    "cpm_factor","crpmsi_factor","unlocalized_fraction"
]
rows: List[Dict[str, Any]] = []
for s in samples:
    r = {k: s.get(k) for k in core_cols}
    for k in region_keys:
        r[f"func_{k}"] = (s.get("func_totals") or {}).get(k, 0)
        r[f"count_{k}"] = (s.get("region_counts") or {}).get(k, 0)
    rows.append(r)

df = pd.DataFrame(rows)
df.to_csv(args.out_tsv, sep="\t", index=False)

# ─────────────────────── Write cohort JSON ───────────────────────
payload = _sanitize(dict(
    n_samples=len(samples),
    rows=rows,
    samples=samples,
    region_totals=region_totals,
    region_keys=region_keys,
    columns=core_cols + [f"func_{k}" for k in region_keys] + [f"count_{k}" for k in region_keys],
    skipped_inputs=skipped
))
with open(args.out_json, "w", encoding="utf-8") as fh:
    json.dump(payload, fh, indent=2)

# Optional global region totals TSV
if args.out_regions and region_totals:
    with open(args.out_regions, "w", encoding="utf-8") as fh:
        fh.write("region\treads\n")
        for k in region_keys:
            fh.write(f"{k}\t{region_totals.get(k,0)}\n")

# ─────────────────────── HTML (single file) ───────────────────────
DATA_JSON = json.dumps(payload, ensure_ascii=False)

CSS = r"""
<style>
:root{
  --bg:#f8fafc; --fg:#0b1220; --muted:#667085; --line:#e5e7eb; --card:#ffffff;
  --pri:#2563eb; --ok:#10b981; --warn:#f59e0b; --fail:#ef4444; --chip:#e5e7eb;
}
@media (prefers-color-scheme: dark){
  :root{ --bg:#0b0d10; --fg:#e5e7eb; --muted:#9aa4b2; --line:#1f2937; --card:#11161c; --chip:#374151; --pri:#60a5fa; }
}
* { box-sizing: border-box; }
html,body{margin:0;background:var(--bg);color:var(--fg);font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;line-height:1.5}
.app{display:flex;flex-direction:column;min-height:100vh}
.header{background:var(--card);border-bottom:1px solid var(--line);padding:1rem 2rem;position:sticky;top:0;z-index:10;display:flex;justify-content:space-between;align-items:center}
.main{flex:1;padding:2rem;max-width:1400px;margin:0 auto;width:100%}
.sidebar{display:none;position:fixed;top:0;left:0;height:100vh;width:280px;background:var(--card);border-right:1px solid var(--line);padding:1rem;overflow-y:auto;z-index:20}
.sidebar.open{display:block}
.sidebar-overlay{display:none;position:fixed;top:0;left:0;right:0;bottom:0;background:rgba(0,0,0,0.5);z-index:15}
.sidebar-overlay.open{display:block}
@media (min-width: 1024px) {
  .app{display:grid;grid-template-columns:280px 1fr}
  .header{display:none}
  .main{padding:2rem;margin:0}
  .sidebar{display:block;position:sticky;top:0;height:100vh}
  .sidebar-overlay{display:none !important}
}
h1{font-size:2rem;margin:0 0 1rem;font-weight:700}
h2{font-size:1.25rem;margin:1.5rem 0 0.5rem;font-weight:600}
h3{font-size:1rem;margin:1rem 0 0.5rem;font-weight:600}
small, .muted{color:var(--muted)}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:1.5rem;margin-bottom:1rem}
.grid{display:grid;gap:1rem;grid-template-columns:repeat(auto-fit,minmax(200px,1fr))}
.kpi-grid{display:grid;gap:1rem;grid-template-columns:repeat(auto-fit,minmax(180px,1fr))}
.chart-grid{display:grid;gap:1rem;grid-template-columns:repeat(auto-fit,minmax(400px,1fr))}
@media (max-width: 768px) {
  .chart-grid{grid-template-columns:1fr}
  .main{padding:1rem}
}
.tbl{border-collapse:collapse;width:100%;font-size:0.875rem;margin-top:1rem}
.tbl th,.tbl td{padding:0.75rem 0.5rem;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
.tbl th{font-weight:600;background:var(--bg);position:sticky;top:0;z-index:5}
.tbl th.sticky,.tbl td.sticky{position:sticky;left:0;background:var(--card);z-index:6;min-width:120px}
.tbl-container{overflow-x:auto;margin-top:1rem;border:1px solid var(--line);border-radius:8px}
.kpi{background:var(--bg);border:1px solid var(--line);border-radius:12px;padding:1rem;text-align:center}
.kpi .lab{font-size:0.75rem;color:var(--muted);text-transform:uppercase;letter-spacing:0.05em;margin-bottom:0.25rem}
.kpi .val{font-size:1.5rem;font-weight:700;margin:0}
.pill{display:inline-block;padding:0.25rem 0.5rem;font-size:0.75rem;border-radius:999px;margin:0.125rem;background:var(--chip);white-space:nowrap}
.status{display:inline-flex;align-items:center;gap:0.5rem;font-weight:700;padding:0.5rem 1rem;border-radius:8px;background:var(--card);border:1px solid var(--line)}
.status::before{content:'';width:8px;height:8px;border-radius:50%}
.hero{background:linear-gradient(135deg, rgba(37,99,235,.1), rgba(37,99,235,0));border:1px solid var(--line);border-radius:12px;padding:1.5rem;margin-bottom:1.5rem}
.footer{color:var(--muted);font-size:0.75rem;margin-top:2rem;padding-top:1rem;border-top:1px solid var(--line)}
.chart{position:relative;min-height:300px;width:100%;overflow:hidden}
.tip{position:absolute;background:var(--card);border:1px solid var(--line);padding:0.5rem;border-radius:6px;font-size:0.75rem;pointer-events:none;opacity:0;z-index:100;box-shadow:0 4px 6px rgba(0,0,0,0.1)}
kbd{background:var(--chip);border:1px solid var(--line);border-radius:3px;padding:0.1em 0.3em;font-size:0.85em;font-family:monospace}
.menu-btn{display:block;background:none;border:1px solid var(--line);padding:0.5rem;border-radius:6px;cursor:pointer;color:var(--fg)}
@media (min-width: 1024px) {
  .menu-btn{display:none}
}
a{color:inherit;text-decoration:none}
a:hover{text-decoration:underline}
.sample-link{display:block;padding:0.5rem;border-radius:6px;margin:0.25rem 0;background:var(--bg);border:1px solid transparent}
.sample-link:hover{background:var(--card);border-color:var(--line)}
</style>
"""

JS = r"""
<script>
// ---------- State & helpers ----------
const DATA = JSON.parse(document.getElementById('payload').textContent || '{}');
const $ = (q, el=document)=>el.querySelector(q);
const $$ = (q, el=document)=>Array.from(el.querySelectorAll(q));
const fmt = n => (n==null||isNaN(n))? 'n/a' : Number(n).toLocaleString();
const pct = x => (x==null||isNaN(x))? 'n/a' : (100*Number(x)).toFixed(1)+'%';
const labelNA = s => (!s || String(s).trim()==='' || String(s).toUpperCase()==='NA') ? 'Unspecified' : String(s);

// ---------- Mobile sidebar toggle ----------
function toggleSidebar() {
  const sidebar = $('.sidebar');
  const overlay = $('.sidebar-overlay');
  sidebar.classList.toggle('open');
  overlay.classList.toggle('open');
}

function sampleColorMap(){
  const ids = [...new Set(DATA.rows.map(r=>r.sample_id))].sort();
  const m={}; for (let i=0;i<ids.length;i++){ m[ids[i]] = `hsl(${Math.round(360*(i/Math.max(1,ids.length)))},70%,50%)`; }
  return sid => m[sid] || '#3b82f6';
}
const colorOfSample = sampleColorMap();

const REGION_COLORS = {
  'pppolii':'rgb(243,132,0)','pppol':'rgb(243,132,0)',
  'divtx':'rgb(178,59,212)','ppdiv':'rgb(178,59,212)',
  'enhancers':'rgb(115,212,122)','enh':'rgb(115,212,122)',
  'genebody':'rgb(0,0,0)','gb':'rgb(0,0,0)',
  'cps':'rgb(103,200,249)','shortgenes':'rgb(253,218,13)','sg':'rgb(253,218,13)','tw':'rgb(255,54,98)'
};
function regionColor(label){
  const k=String(label||'').toLowerCase().replace(/\s+/g,'');
  if (REGION_COLORS[k]) return REGION_COLORS[k];
  if (k==='non-localized'||k==='nonlocalized') return 'rgb(160,170,180)';
  // fallback hash
  let h=0; for (let i=0;i<k.length;i++) h=(h*31+k.charCodeAt(i))|0;
  const P=['#3b82f6','#10b981','#f59e0b','#ef4444','#8b5cf6','#06b6d4','#84cc16','#ec4899','#22c55e','#f43f5e','#14b8a6','#a855f7'];
  return P[Math.abs(h)%P.length];
}

// ---------- Charts (SVG, exportable) ----------
function exportPNG(svg, name='chart.png'){
  const s=new XMLSerializer().serializeToString(svg),
        img=new Image(),
        url=URL.createObjectURL(new Blob([s],{type:'image/svg+xml'}));
  img.onload=()=>{ const r=svg.viewBox.baseVal, c=document.createElement('canvas');
    c.width=r&&r.width? r.width: svg.clientWidth*2; c.height=r&&r.height? r.height: svg.clientHeight*2;
    const ctx=c.getContext('2d'); ctx.fillStyle=getComputedStyle(document.documentElement).getPropertyValue('--bg')||'#fff';
    ctx.fillRect(0,0,c.width,c.height); ctx.drawImage(img,0,0,c.width,c.height);
    c.toBlob(b=>{const a=document.createElement('a'); a.href=URL.createObjectURL(b); a.download=name; a.click();});
  }; img.src=url;
}

function barChart(el, labels, values, colors, ylab, onClick){
  const containerWidth = el.parentElement.clientWidth || 600;
  const W = Math.min(Math.max(containerWidth - 40, 400), 1200);
  const H = 360;
  const L = 70;
  const B = 46;
  const T = 8;
  const gap = Math.max(8, W / 100);
  const n = labels.length || 1;
  const max = Math.max(...values, 1);
  const bw = Math.max(12, Math.floor((W - L - 16 - (n - 1) * gap) / n));
  
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.style.width = '100%';
  svg.style.height = 'auto';
  
  const tip = document.createElement('div');
  tip.className = 'tip';
  el.style.position = 'relative';
  el.appendChild(tip);

  // grid + axis ticks
  for(let i=0;i<=5;i++){ const y=T+(H-B-T)*i/5; const ln=document.createElementNS(svg.namespaceURI,'line');
    ln.setAttribute('x1',L); ln.setAttribute('x2',W-8); ln.setAttribute('y1',y); ln.setAttribute('y2',y);
    ln.setAttribute('stroke','currentColor'); ln.setAttribute('opacity','0.12'); svg.appendChild(ln);
    const t=document.createElementNS(svg.namespaceURI,'text'); t.setAttribute('x',14); t.setAttribute('y',y+4); t.setAttribute('font-size','12');
    t.setAttribute('fill','currentColor'); t.setAttribute('opacity','0.6'); t.textContent=Math.round(max*(1-i/5)).toLocaleString(); svg.appendChild(t);
  }

  const showEvery = n>24? Math.ceil(n/24): 1;
  labels.forEach((lab,i)=>{
    const v=values[i]||0, h=Math.round((v/max)*(H-B-T)), x=L+i*(bw+gap), y=(H-B)-h;
    const r=document.createElementNS(svg.namespaceURI,'rect'); r.setAttribute('x',x); r.setAttribute('y',y);
    r.setAttribute('width',bw); r.setAttribute('height',h); r.setAttribute('rx','4');
    r.setAttribute('fill', colors[i] || '#3b82f6'); r.setAttribute('opacity','0.95'); if (onClick) r.style.cursor='pointer';
    r.addEventListener('mousemove', (e)=>{ tip.style.left=(e.offsetX+12)+'px'; tip.style.top=(e.offsetY-10)+'px';
      tip.textContent=`${lab}: ${v.toLocaleString()}`; tip.style.opacity='1'; });
    r.addEventListener('mouseleave', ()=> tip.style.opacity='0');
    if (onClick) r.addEventListener('click', ()=> onClick(lab));
    svg.appendChild(r);
    if (i%showEvery===0){ const t=document.createElementNS(svg.namespaceURI,'text'); t.setAttribute('x',x+bw/2); t.setAttribute('y',H-12);
      t.setAttribute('text-anchor','middle'); t.setAttribute('font-size','12'); t.setAttribute('fill','currentColor'); t.setAttribute('opacity','0.8');
      t.textContent=lab; svg.appendChild(t); }
  });

  const yl=document.createElementNS(svg.namespaceURI,'text'); yl.setAttribute('x',16); yl.setAttribute('y',16);
  yl.setAttribute('font-size','12'); yl.setAttribute('fill','currentColor'); yl.setAttribute('opacity','0.8'); yl.textContent=ylab; svg.appendChild(yl);
  el.innerHTML=''; el.appendChild(svg);
}

function scatter(el, labels, xs, ys, colors, xl, yl, onClick){
  const containerWidth = el.parentElement.clientWidth || 600;
  const W = Math.min(Math.max(containerWidth - 40, 400), 1200);
  const H = 360;
  const L = 70;
  const B = 46;
  const T = 8;
  
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.style.width = '100%';
  svg.style.height = 'auto';
  
  const tip = document.createElement('div');
  tip.className = 'tip';
  el.style.position = 'relative';
  el.appendChild(tip);
  
  const maxX = Math.max(...xs.filter(v => v != null), 1);
  const maxY = Math.max(...ys.filter(v => v != null), 1);

  // axes
  const ax=document.createElementNS(svg.namespaceURI,'line'); ax.setAttribute('x1',L); ax.setAttribute('x2',W-8); ax.setAttribute('y1',H-B); ax.setAttribute('y2',H-B);
  ax.setAttribute('stroke','currentColor'); ax.setAttribute('opacity','0.2'); svg.appendChild(ax);
  const ay=document.createElementNS(svg.namespaceURI,'line'); ay.setAttribute('x1',L); ay.setAttribute('x2',L); ay.setAttribute('y1',T); ay.setAttribute('y2',H-B);
  ay.setAttribute('stroke','currentColor'); ay.setAttribute('opacity','0.2'); svg.appendChild(ay);

  labels.forEach((lab,i)=>{
    const x=L + (xs[i]/maxX)*(W-L-10), y=H-B - (ys[i]/maxY)*(H-B-T);
    const c=document.createElementNS(svg.namespaceURI,'circle'); c.setAttribute('cx',x); c.setAttribute('cy',y); c.setAttribute('r',5);
    c.setAttribute('fill', colors[i] || '#3b82f6'); c.setAttribute('opacity','0.7'); c.style.cursor='pointer';
    c.addEventListener('mousemove', (e)=>{ tip.style.left=(e.offsetX+12)+'px'; tip.style.top=(e.offsetY-10)+'px';
      tip.textContent=`${lab}: ${xl} ${(xs[i]??0).toFixed(2)}, ${yl} ${(ys[i]??0).toFixed(2)}`; tip.style.opacity='1'; });
    c.addEventListener('mouseleave', ()=> tip.style.opacity='0');
    if (onClick) c.addEventListener('click', ()=> onClick(lab));
    svg.appendChild(c);
  });
  const xt=document.createElementNS(svg.namespaceURI,'text'); xt.setAttribute('x',W/2); xt.setAttribute('y',H-16);
  xt.setAttribute('font-size','12'); xt.setAttribute('text-anchor','middle'); xt.setAttribute('fill','currentColor'); xt.textContent=xl; svg.appendChild(xt);
  const yt=document.createElementNS(svg.namespaceURI,'text'); yt.setAttribute('x',16); yt.setAttribute('y',H/2);
  yt.setAttribute('font-size','12'); yt.setAttribute('text-anchor','middle'); yt.setAttribute('transform',`rotate(-90 16,${H/2})`);
  yt.setAttribute('fill','currentColor'); yt.textContent=yl; svg.appendChild(yt);
  el.innerHTML=''; el.appendChild(svg);
}

function donut(el, labels, values, colorFor){
  const containerWidth = el.parentElement.clientWidth || 600;
  const W = Math.min(Math.max(containerWidth - 40, 400), 600);
  const H = 340;
  const R = Math.min(120, W / 4);
  const cx = W / 2;
  const cy = H / 2;
  const rIn = R * 0.5;
  
  const total = values.reduce((a,b) => a + b, 0) || 1;
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.style.width = '100%';
  svg.style.height = 'auto';
  
  const tip = document.createElement('div');
  tip.className = 'tip';
  el.style.position = 'relative';
  el.appendChild(tip);
  let a0=-Math.PI/2;
  values.forEach((v,i)=>{
    const a1=a0+(v/total)*2*Math.PI;
    const x0=cx+R*Math.cos(a0), y0=cy+R*Math.sin(a0), x1=cx+R*Math.cos(a1), y1=cy+R*Math.sin(a1), large=(a1-a0)>Math.PI?1:0;
    const path=document.createElementNS(svg.namespaceURI,'path');
    path.setAttribute('d',`M ${cx+rIn*Math.cos(a0)} ${cy+rIn*Math.sin(a0)} L ${x0} ${y0} A ${R} ${R} 0 ${large} 1 ${x1} ${y1} L ${cx+rIn*Math.cos(a1)} ${cy+rIn*Math.sin(a1)} A ${rIn} ${rIn} 0 ${large} 0 ${cx+rIn*Math.cos(a0)} ${cy+rIn*Math.sin(a0)} Z`);
    path.setAttribute('fill', colorFor(labels[i])); path.setAttribute('opacity','0.95');
    path.addEventListener('mousemove', (e)=>{ tip.style.left=(e.offsetX+12)+'px'; tip.style.top=(e.offsetY-10)+'px';
      tip.textContent=`${labels[i]}: ${values[i].toLocaleString()} (${((values[i]/total)*100).toFixed(1)}%)`; tip.style.opacity='1'; });
    path.addEventListener('mouseleave', ()=> tip.style.opacity='0');
    svg.appendChild(path); a0=a1;
  });
  el.innerHTML=''; el.appendChild(svg);
}

function build(){
  // Hero KPIs & badge
  const rows = DATA.rows || [];
  const dupVals = rows.map(r=>Number(r.duplicate_percent)).filter(v=>!Number.isNaN(v)).sort((a,b)=>a-b);
  const medianDup = dupVals.length ? (dupVals.length%2 ? dupVals[Math.floor(dupVals.length/2)] : (dupVals[dupVals.length/2-1]+dupVals[dupVals.length/2])/2) : null;
  const unlocVals = rows.map(r=>Number(r.unlocalized_fraction)).filter(v=>!Number.isNaN(v)).sort((a,b)=>a-b);
  const medianUnloc = unlocVals.length ? (unlocVals.length%2 ? unlocVals[Math.floor(unlocVals.length/2)] : (unlocVals[unlocVals.length/2-1]+unlocVals[unlocVals.length/2])/2) : null;

  const badge = (()=>{ const ok = (medianUnloc!=null && medianUnloc<=0.25) && (medianDup!=null && medianDup<=30);
    return { label: ok? 'PASS' : (medianUnloc!=null && medianUnloc<=0.4 ? 'WARN' : 'FAIL'), color: ok? 'var(--ok)' : (medianUnloc!=null && medianUnloc<=0.4 ? 'var(--warn)' : 'var(--fail)')};
  })();

  $('#badge').innerHTML = `<span class="status" style="--status:${badge.color}"><span style="width:10px;height:10px;border-radius:50%;background:${badge.color}"></span>${badge.label}</span>`;
  $('#kpi-n').textContent = rows.length;
  $('#kpi-reads').textContent = (rows.reduce((a,r)=>a+(Number(r.reads_total_functional)||0),0)).toLocaleString();
  $('#kpi-dup').textContent = (medianDup==null? 'n/a' : medianDup.toFixed(2)+'%');
  $('#kpi-reg').textContent = (DATA.region_keys||[]).length;
  $('#kpi-unloc').textContent = (medianUnloc==null? 'n/a' : (100*medianUnloc).toFixed(1)+'%');
  
  // Quality status summary
  const passCount = rows.filter(r => {
    const reads = r.input_reads || 0;
    const dup = r.duplicate_percent || 0;
    const unloc = r.unlocalized_fraction || 0;
    return reads >= 5e6 && dup < 15 && unloc < 0.25;
  }).length;
  const warnCount = rows.filter(r => {
    const reads = r.input_reads || 0;
    const dup = r.duplicate_percent || 0;
    const unloc = r.unlocalized_fraction || 0;
    return reads >= 2e6 && dup < 30 && unloc < 0.4 && !(reads >= 5e6 && dup < 15 && unloc < 0.25);
  }).length;
  const failCount = rows.length - passCount - warnCount;
  $('#kpi-quality').innerHTML = `<span style="color: #10b981">${passCount}P</span> <span style="color: #f59e0b">${warnCount}W</span> <span style="color: #ef4444">${failCount}F</span>`;

  // Charts
  const labels = rows.map(r=>r.sample_id);
  const colors = labels.map(colorOfSample);
  const reads  = rows.map(r=>Number(r.reads_total_functional)||0);
  const piMed  = rows.map(r=>Number(r.median_pausing_index)||0);
  const dens   = rows.map(r=>Number(r.median_density)||0);

  // Charts with enhanced information
  barChart($('#chart-reads'), labels, reads, colors, 'Reads in functional regions', sid => location.hash = '#sample/'+encodeURIComponent(sid));
  
  // Calculate correlation for scatter plot
  const validPairs = piMed.map((pi, i) => [pi, dens[i]]).filter(([pi, d]) => pi > 0 && d > 0);
  const correlation = validPairs.length > 1 ? 
    calculateCorrelation(validPairs.map(p => p[0]), validPairs.map(p => p[1])) : 0;
  
  scatter($('#chart-scatter'), labels, piMed, dens, colors, 'Median Pausing Index', 'Median Density', sid => location.hash = '#sample/'+encodeURIComponent(sid));
  
  // Add correlation info below scatter plot
  const scatterContainer = $('#chart-scatter').parentElement;
  if (scatterContainer && !scatterContainer.querySelector('.correlation-info')) {
    const corrInfo = document.createElement('div');
    corrInfo.className = 'correlation-info';
    corrInfo.style.cssText = 'margin-top: 0.5rem; text-align: center; font-size: 0.9em; color: var(--muted);';
    corrInfo.innerHTML = `<strong>Correlation (PI vs Density):</strong> r = ${correlation.toFixed(3)} ${Math.abs(correlation) > 0.7 ? '(Strong)' : Math.abs(correlation) > 0.3 ? '(Moderate)' : '(Weak)'}`;
    scatterContainer.appendChild(corrInfo);
  }
  
  const rk = DATA.region_keys||[]; donut($('#chart-pie'), rk, rk.map(k=>(DATA.region_totals||{})[k]||0), s=>regionColor(s));

  // Table
  const tb = $('#table-body'); tb.innerHTML='';
  rows.forEach(r=>{
    const tr=document.createElement('tr');
    const unlocCell = (r.unlocalized_fraction==null || isNaN(r.unlocalized_fraction)) ? 'n/a' : (100*Number(r.unlocalized_fraction)).toFixed(1)+'%';
    tr.innerHTML = `
      <td class="sticky">${r.sample_id}</td>
      <td>${labelNA(r.condition)}</td>
      <td>${labelNA(r.timepoint)}</td>
      <td>${labelNA(r.replicate)}</td>
      <td>${fmt(r.input_reads)}</td>
      <td>${fmt(r.dedup_reads)}</td>
      <td>${r.duplicate_percent!=null? (r.umi_deduplication_enabled? 'UMI: ' : '') + Number(r.duplicate_percent).toFixed(2)+'%':'n/a'}</td>
      <td>${fmt(r.divergent_regions)}</td>
      <td>${fmt(r.total_regions)}</td>
      <td>${fmt(r.reads_total_functional)}</td>
      <td>${r.median_pausing_index!=null? Number(r.median_pausing_index).toFixed(2):'n/a'}</td>
      <td>${r.median_density!=null? Number(r.median_density).toFixed(2):'n/a'}</td>
      <td>${r.cpm_factor!=null? Number(r.cpm_factor).toFixed(2):'n/a'}</td>
      <td>${r.crpmsi_factor!=null? Number(r.crpmsi_factor).toFixed(2):'n/a'}</td>
      <td>${unlocCell}</td>
      </td>`;
    tr.onclick = ()=> location.hash = '#sample/'+encodeURIComponent(r.sample_id);
    tb.appendChild(tr);
  });
  
  // Initialize filter options
  initializeFilters(rows);
  
  // Initialize table info
  const info = $('#table-info');
  if (info) info.textContent = `${rows.length} samples`;

  // Functional regions table
  const ftb = $('#functional-table-body'); ftb.innerHTML='';
  rows.forEach(r=>{
    const tr=document.createElement('tr');
    tr.innerHTML = `
      <td class="sticky">${r.sample_id}</td>
      <td>${labelNA(r.condition)}</td>
      <td>${fmt(r['func_Promoter'] || 0)}</td>
      <td>${fmt(r['func_Gene body'] || 0)}</td>
      <td>${fmt(r['func_CPS'] || 0)}</td>
      <td>${fmt(r['func_DivergentTx'] || 0)}</td>
      <td>${fmt(r['func_Enhancers'] || 0)}</td>
      <td>${fmt(r['func_Short genes'] || 0)}</td>
      <td>${fmt(r['func_Termination window'] || 0)}</td>
      <td>${fmt(r['func_Non-localized polymerase'] || 0)}</td>`;
    tr.onclick = ()=> location.hash = '#sample/'+encodeURIComponent(r.sample_id);
    ftb.appendChild(tr);
  });

  // Sidebar: samples list
  const sl = $('#sample-links'); sl.innerHTML='';
  rows.forEach(r=>{
    const a = document.createElement('a');
    a.href = '#sample/' + encodeURIComponent(r.sample_id);
    a.textContent = r.sample_id;
    a.className = 'sample-link';
    a.onclick = () => { if (window.innerWidth < 1024) toggleSidebar(); };
    sl.appendChild(a);
  });

  // Per-sample pages
  const pages = $('#sample-pages'); pages.innerHTML='';
  (DATA.samples||[]).forEach(s=>{
    const div=document.createElement('div'); div.id='sample/'+encodeURIComponent(s.sample_id); div.className='card';
    const totalReads = fmt(s.input_reads), funcReads = fmt(s.reads_total_functional);
    const unloc = (s.unlocalized_fraction==null || isNaN(s.unlocalized_fraction)) ? 'n/a' : (100*Number(s.unlocalized_fraction)).toFixed(1)+'%';
    div.innerHTML = `
      <h2>Sample: ${s.sample_id}</h2>
      <div class="kpi-grid">
        <div class="kpi"><div class="lab">Total Reads</div><div class="val">${totalReads}</div></div>
        <div class="kpi"><div class="lab">Functional Reads</div><div class="val">${funcReads}</div></div>
        <div class="kpi"><div class="lab">${s.umi_deduplication_enabled? 'UMI Dedup %' : 'Duplicate %'}</div><div class="val">${s.duplicate_percent!=null? Number(s.duplicate_percent).toFixed(2)+'%':'n/a'}</div></div>
        <div class="kpi"><div class="lab">Median PI</div><div class="val">${s.median_pausing_index!=null? Number(s.median_pausing_index).toFixed(2):'n/a'}</div></div>
        <div class="kpi"><div class="lab">Median Density</div><div class="val">${s.median_density!=null? Number(s.median_density).toFixed(2):'n/a'}</div></div>
        <div class="kpi"><div class="lab">Unlocalized</div><div class="val">${unloc}</div></div>
      </div>
      <div class="card">
        <h3>Region composition</h3>
        <div id="reg-${s.sample_id}" class="chart"></div>
      </div>`;
    pages.appendChild(div);
    const labels=(s.regions||[]).map(r=>r.region), values=(s.regions||[]).map(r=>Number(r.reads)||0),
          colors=labels.map(regionColor);
    barChart($('#reg-'+s.sample_id), labels, values, colors, 'Reads', null);
  });

  // Routing
  function show(hash){
    const h=hash||'#overview';
    $$('.section').forEach(s=> s.style.display='none');
    if (h.startsWith('#sample/')){
      $('#overview').style.display='none';
      $('#sample-section').style.display='block';
      const id = decodeURIComponent(h.split('/')[1]||'');
      const target = document.getElementById('sample/'+encodeURIComponent(id));
      if (target){ target.scrollIntoView({behavior:'smooth',block:'start'}); }
    } else {
      $('#sample-section').style.display='none';
      $('#overview').style.display='block';
    }
  }
  window.addEventListener('hashchange', ()=> show(location.hash));
  show(location.hash);
}

// ---------- Theme toggle ----------
function toggleTheme() {
  const currentTheme = document.documentElement.getAttribute('data-theme');
  const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
  document.documentElement.setAttribute('data-theme', newTheme);
  localStorage.setItem('theme', newTheme);
  
  // Update button icon
  const btn = document.querySelector('button[onclick="toggleTheme()"]');
  if (btn) btn.textContent = newTheme === 'dark' ? '☀️' : '🌙';
}

// Initialize theme on load
function initTheme() {
  const savedTheme = localStorage.getItem('theme') || 
    (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
  document.documentElement.setAttribute('data-theme', savedTheme);
  
  // Update button icon
  const btn = document.querySelector('button[onclick="toggleTheme()"]');
  if (btn) btn.textContent = savedTheme === 'dark' ? '☀️' : '🌙';
}

// Copy to clipboard functionality
function copyToClipboard(text) {
  navigator.clipboard.writeText(text).then(() => {
    // Show temporary feedback
    const btn = event.target;
    const originalText = btn.textContent;
    btn.textContent = '✅ Copied!';
    btn.style.background = 'var(--ok)';
    setTimeout(() => {
      btn.textContent = originalText;
      btn.style.background = 'var(--accent)';
    }, 2000);
  }).catch(() => {
    // Fallback for older browsers
    const textArea = document.createElement('textarea');
    textArea.value = text;
    document.body.appendChild(textArea);
    textArea.select();
    document.execCommand('copy');
    document.body.removeChild(textArea);
    
    const btn = event.target;
    const originalText = btn.textContent;
    btn.textContent = '✅ Copied!';
    setTimeout(() => btn.textContent = originalText, 2000);
  });
}

// Enhanced table filter functionality with statistics
function filterTable() {
  const searchTerm = document.getElementById('sample-search').value.toLowerCase();
  const conditionFilter = document.getElementById('condition-filter').value;
  const qualityFilter = document.getElementById('quality-filter').value;
  const rows = document.querySelectorAll('#table-body tr');
  
  let visibleCount = 0;
  const filteredData = [];
  
  rows.forEach(row => {
    const cells = row.querySelectorAll('td');
    const sampleId = cells[0]?.textContent || '';
    const condition = cells[1]?.textContent || '';
    const timepoint = cells[2]?.textContent || '';
    const replicate = cells[3]?.textContent || '';
    const inputReads = parseFloat(cells[4]?.textContent.replace(/,/g, '') || '0');
    const dupPercent = parseFloat(cells[6]?.textContent.replace(/[^0-9.]/g, '') || '0');
    const functionalReads = parseFloat(cells[9]?.textContent.replace(/,/g, '') || '0');
    const medianPI = parseFloat(cells[10]?.textContent || '0');
    const unlocPercent = parseFloat(cells[14]?.textContent.replace(/[^0-9.]/g, '') || '0');
    
    // Text search
    const textContent = [sampleId, condition, timepoint, replicate].join(' ').toLowerCase();
    const matchesText = !searchTerm || textContent.includes(searchTerm);
    
    // Condition filter
    const matchesCondition = !conditionFilter || condition === conditionFilter;
    
    // Quality filter
    let matchesQuality = true;
    if (qualityFilter) {
      const isPass = inputReads >= 5e6 && dupPercent < 15 && unlocPercent < 25;
      const isWarn = inputReads >= 2e6 && dupPercent < 30 && unlocPercent < 40;
      
      if (qualityFilter === 'pass') matchesQuality = isPass;
      else if (qualityFilter === 'warn') matchesQuality = isPass || isWarn;
      else if (qualityFilter === 'fail') matchesQuality = !isPass && !isWarn;
    }
    
    const isVisible = matchesText && matchesCondition && matchesQuality;
    row.style.display = isVisible ? '' : 'none';
    
    if (isVisible) {
      visibleCount++;
      filteredData.push({
        sampleId, condition, inputReads, functionalReads, medianPI, dupPercent, unlocPercent
      });
    }
  });
  
  updateTableStats(visibleCount, rows.length, filteredData);
}

function updateTableStats(visibleCount, totalCount, filteredData) {
  const info = document.getElementById('table-info');
  const stats = document.getElementById('filter-stats');
  
  if (info) {
    info.textContent = visibleCount === totalCount ? 
      `${totalCount} samples` : 
      `Showing ${visibleCount} of ${totalCount} samples`;
  }
  
  if (stats && filteredData.length > 0) {
    const avgReads = Math.round(filteredData.reduce((a, b) => a + b.functionalReads, 0) / filteredData.length);
    const medianPI = filteredData.map(d => d.medianPI).filter(v => v > 0).sort((a,b) => a-b);
    const medianPIValue = medianPI.length > 0 ? 
      (medianPI.length % 2 ? medianPI[Math.floor(medianPI.length/2)] : 
       (medianPI[medianPI.length/2-1] + medianPI[medianPI.length/2])/2) : 0;
    
    stats.innerHTML = `
      <strong>Selection Stats:</strong> 
      Avg Functional Reads: ${avgReads.toLocaleString()} | 
      Median PI: ${medianPIValue.toFixed(2)}
    `;
  } else if (stats) {
    stats.innerHTML = '';
  }
}

function clearFilters() {
  document.getElementById('sample-search').value = '';
  document.getElementById('condition-filter').value = '';
  document.getElementById('quality-filter').value = '';
  filterTable();
}

function exportFilteredCSV() {
  const rows = Array.from(document.querySelectorAll('#table-body tr'))
    .filter(row => row.style.display !== 'none');
  
  if (rows.length === 0) {
    alert('No data to export. Please adjust your filters.');
    return;
  }
  
  const headers = ['Sample', 'Condition', 'Timepoint', 'Replicate', 'Input Reads', 'Dedup Reads', 
                   'Duplicate %', 'Divergent Regions', 'Total Regions', 'Functional Reads', 
                   'Median PI', 'Median Density', 'CPM Factor', 'cRPMsi Factor', 'Unlocalized %'];
  
  let csv = headers.join(',') + '\\n';
  
  rows.forEach(row => {
    const cells = Array.from(row.querySelectorAll('td'));
    const values = cells.map(cell => {
      let text = cell.textContent.trim();
      // Escape commas and quotes in CSV
      if (text.includes(',') || text.includes('"')) {
        text = '"' + text.replace(/"/g, '""') + '"';
      }
      return text;
    });
    csv += values.join(',') + '\\n';
  });
  
  // Download CSV
  const blob = new Blob([csv], { type: 'text/csv' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `tracktx_filtered_samples_${new Date().toISOString().split('T')[0]}.csv`;
  a.click();
  URL.revokeObjectURL(url);
}

function initializeFilters(rows) {
  // Populate condition filter
  const conditionFilter = document.getElementById('condition-filter');
  if (conditionFilter) {
    const conditions = [...new Set(rows.map(r => r.condition))].filter(c => c && c !== 'NA').sort();
    conditions.forEach(condition => {
      const option = document.createElement('option');
      option.value = condition;
      option.textContent = condition;
      conditionFilter.appendChild(option);
    });
  }
}

function calculateCorrelation(xs, ys) {
  const n = xs.length;
  if (n < 2) return 0;
  
  const meanX = xs.reduce((a, b) => a + b, 0) / n;
  const meanY = ys.reduce((a, b) => a + b, 0) / n;
  
  let num = 0, denX = 0, denY = 0;
  for (let i = 0; i < n; i++) {
    const dx = xs[i] - meanX;
    const dy = ys[i] - meanY;
    num += dx * dy;
    denX += dx * dx;
    denY += dy * dy;
  }
  
  const den = Math.sqrt(denX * denY);
  return den === 0 ? 0 : num / den;
}

// Keyboard shortcuts
document.addEventListener('keydown', (e) => {
  // Ctrl+F or Cmd+F - focus search box
  if ((e.ctrlKey || e.metaKey) && e.key === 'f') {
    e.preventDefault();
    const searchBox = document.getElementById('sample-search');
    if (searchBox) {
      searchBox.focus();
      searchBox.select();
    }
  }
  
  // Escape - clear all filters
  if (e.key === 'Escape') {
    clearFilters();
  }
  
  // Ctrl+E or Cmd+E - export CSV
  if ((e.ctrlKey || e.metaKey) && e.key === 'e') {
    e.preventDefault();
    exportFilteredCSV();
  }
});

window.addEventListener('load', () => {
  initTheme();
  build();
});
</script>
"""

HTML = f"""<!doctype html>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TrackTx — Cohort Summary</title>
{CSS}
<div class="app">
  <header class="header">
    <h1>TrackTx Cohort Summary</h1>
    <div style="display: flex; gap: 0.5rem; align-items: center;">
      <button onclick="toggleTheme()" style="background: none; border: 1px solid var(--border); border-radius: 4px; padding: 0.5rem; cursor: pointer; font-size: 1.2em;" title="Toggle dark/light mode">🌙</button>
      <button class="menu-btn" onclick="toggleSidebar()">☰ Menu</button>
    </div>
  </header>
  
  <div class="sidebar-overlay" onclick="toggleSidebar()"></div>
  <aside class="sidebar">
    <div class="hero">
      <h2 style="margin:0 0 0.5rem">TrackTx</h2>
      <div id="badge"></div>
      <small class="muted">Single-file cohort report</small>
    </div>
    <h2>Samples</h2>
    <div id="sample-links"></div>
    <div class="footer">
      <div>TrackTx-NF v{args.pipeline_version}</div>
      <small class="muted">Run: {args.run_name} • Profile: {args.profile} • Duration: {args.duration}</small>
    </div>
  </aside>

  <main class="main">
    <section id="overview" class="section" style="display:block">
      <h1>Cohort Overview</h1>
      
      <!-- Run Metadata Card -->
      <div class="card">
        <h3>🔬 Analysis Summary</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 1rem; margin-top: 1rem;">
          <div>
            <strong>Pipeline:</strong> TrackTx-NF v{args.pipeline_version}<br>
            <strong>Profile:</strong> {args.profile}<br>
            <strong>Duration:</strong> {args.duration}
          </div>
          <div>
            <strong>Run Name:</strong> {args.run_name}<br>
            <strong>Generated:</strong> {datetime.datetime.now().strftime("%Y-%m-%d %H:%M %Z")}<br>
            <strong>Report Type:</strong> Interactive Cohort Summary
          </div>
        </div>
      </div>

      <div class="kpi-grid">
        <div class="kpi"><div class="lab">Total Samples</div><div id="kpi-n" class="val"></div></div>
        <div class="kpi"><div class="lab">Functional Reads</div><div id="kpi-reads" class="val"></div></div>
        <div class="kpi"><div class="lab">Median Duplicate %</div><div id="kpi-dup" class="val"></div></div>
        <div class="kpi"><div class="lab">Region Types</div><div id="kpi-reg" class="val"></div></div>
        <div class="kpi"><div class="lab">Median Unlocalized</div><div id="kpi-unloc" class="val"></div></div>
        <div class="kpi"><div class="lab">Quality Status</div><div id="kpi-quality" class="val"></div></div>
      </div>

      <h2>Quality Control Overview</h2>
      <div class="card">
        <h3>What do these metrics mean?</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1rem; margin-top: 1rem;">
          <div>
            <strong>Duplicate %:</strong> Percentage of reads removed as duplicates. UMI deduplication (shown as "UMI: X%") is more accurate than PCR duplicate detection.
          </div>
          <div>
            <strong>Unlocalized %:</strong> Percentage of reads not mapping to functional regions. Values >40% may indicate issues.
          </div>
          <div>
            <strong>Functional Reads:</strong> Reads mapping to biologically relevant regions (promoters, gene bodies, etc.).
          </div>
          <div>
            <strong>Pausing Index:</strong> Measures RNA polymerase pausing. Higher values indicate more transcriptional pausing.
          </div>
        </div>
      </div>

      <h2>Charts</h2>
      <div class="chart-grid">
        <div class="card"><div id="chart-reads" class="chart"></div></div>
        <div class="card"><div id="chart-scatter" class="chart"></div></div>
        <div class="card"><div id="chart-pie" class="chart"></div></div>
      </div>

      <h2>Samples Table</h2>
      <div class="card" style="margin: 1rem 0;">
        <h3>🔍 Advanced Filtering & Analysis</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1rem; margin-top: 1rem;">
          <div>
            <label style="display: block; font-weight: 600; margin-bottom: 0.25rem;">Text Search</label>
            <input type="text" id="sample-search" placeholder="🔍 Search samples, conditions, timepoints..." style="width: 100%; padding: 0.5rem; border: 1px solid var(--line); border-radius: 4px; background: var(--card); color: var(--fg);" oninput="filterTable()">
          </div>
          <div>
            <label style="display: block; font-weight: 600; margin-bottom: 0.25rem;">Condition</label>
            <select id="condition-filter" style="width: 100%; padding: 0.5rem; border: 1px solid var(--line); border-radius: 4px; background: var(--card); color: var(--fg);" onchange="filterTable()">
              <option value="">All Conditions</option>
            </select>
          </div>
          <div>
            <label style="display: block; font-weight: 600; margin-bottom: 0.25rem;">Quality Status</label>
            <select id="quality-filter" style="width: 100%; padding: 0.5rem; border: 1px solid var(--line); border-radius: 4px; background: var(--card); color: var(--fg);" onchange="filterTable()">
              <option value="">All Quality</option>
              <option value="pass">PASS Only</option>
              <option value="warn">WARN+ (PASS + WARN)</option>
              <option value="fail">FAIL Only</option>
            </select>
          </div>
          <div>
            <label style="display: block; font-weight: 600; margin-bottom: 0.25rem;">Actions</label>
            <div style="display: flex; gap: 0.5rem;">
              <button onclick="clearFilters()" style="flex: 1; padding: 0.5rem; border: 1px solid var(--line); border-radius: 4px; background: var(--card); color: var(--fg); cursor: pointer;">Clear All</button>
              <button onclick="exportFilteredCSV()" style="flex: 1; padding: 0.5rem; border: 1px solid var(--accent); border-radius: 4px; background: var(--accent); color: white; cursor: pointer;">📊 Export CSV</button>
            </div>
          </div>
        </div>
        <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 1rem;">
          <small id="table-info" class="muted"></small>
          <div id="filter-stats" style="font-size: 0.9em;"></div>
        </div>
        <div style="margin-top: 0.5rem; font-size: 0.85em; color: var(--muted);">
          💡 <strong>Tips:</strong> 
          <kbd>Ctrl+F</kbd> for quick search | 
          <kbd>Escape</kbd> to clear filters | 
          Click chart points to jump to samples
        </div>
      </div>
      <div class="tbl-container">
        <table class="tbl">
          <thead>
            <tr>
              <th class="sticky">Sample</th>
              <th>Condition</th>
              <th>Timepoint</th>
              <th>Replicate</th>
              <th>Input</th>
              <th>Dedup</th>
              <th>Dup %</th>
              <th>Divergent</th>
              <th>Total Regions</th>
              <th>Functional Reads</th>
              <th>Median PI</th>
              <th>Median Density</th>
              <th>CPM</th>
              <th>cRPMsi</th>
              <th>Unlocalized</th>
            </tr>
          </thead>
          <tbody id="table-body"></tbody>
        </table>
      </div>

      <h2>Functional Region Breakdown</h2>
      <div class="tbl-container">
        <table class="tbl">
          <thead>
            <tr>
              <th class="sticky">Sample</th>
              <th>Condition</th>
              <th>Promoter</th>
              <th>Gene Body</th>
              <th>CPS</th>
              <th>DivergentTx</th>
              <th>Enhancers</th>
              <th>Short Genes</th>
              <th>Termination</th>
              <th>Non-localized</th>
            </tr>
          </thead>
          <tbody id="functional-table-body"></tbody>
        </table>
      </div>

      <h2>📁 Results Navigation Guide</h2>
      <div class="card">
        <h3>🎯 Quick Access</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 1rem; margin-top: 1rem;">
          <div>
            <strong>📊 Interactive Reports</strong><br>
            <code>11_reports/cohort/global_summary.html</code> ← You are here!<br>
            <code>11_reports/samples/&lt;SAMPLE&gt;/&lt;SAMPLE&gt;.report.html</code><br>
            <small class="muted">Rich HTML reports with interactive visualizations</small>
          </div>
          <div>
            <strong>🧬 Genome Browser Tracks</strong><br>
            <code>03_genome_tracks/&lt;SAMPLE&gt;/3p/*.bw</code><br>
            <code>05_normalized_tracks/&lt;SAMPLE&gt;/*.bw</code><br>
            <small class="muted">BigWig files for IGV, UCSC Genome Browser</small>
          </div>
          <div>
            <strong>📋 Data Tables</strong><br>
            <code>global_summary.tsv</code> ← Cohort data<br>
            <code>11_reports/samples/&lt;SAMPLE&gt;/&lt;SAMPLE&gt;.report.tsv</code><br>
            <small class="muted">Machine-readable data for downstream analysis</small>
          </div>
          <div>
            <strong>🔬 Analysis Results</strong><br>
            <code>08_pol2_metrics/&lt;SAMPLE&gt;/</code> ← Pol II metrics<br>
            <code>06_divergent_tx/&lt;SAMPLE&gt;/</code> ← Divergent transcription<br>
            <small class="muted">Detailed analysis outputs and intermediate files</small>
          </div>
        </div>
      </div>

      <div class="card">
        <h3>📂 Complete Directory Structure</h3>
        <div style="font-family: monospace; font-size: 0.9em; line-height: 1.4; background: #f8f9fa; padding: 1rem; border-radius: 8px; overflow-x: auto;">
results/<br>
├── 📁 <strong>00_references/</strong> ← Genome references and indices<br>
├── 📁 <strong>01_trimmed_fastq/&lt;SAMPLE&gt;/</strong> ← Processed FASTQ files and QC<br>
├── 📁 <strong>02_alignments/&lt;SAMPLE&gt;/</strong> ← BAM files and alignment stats<br>
├── 📁 <strong>03_genome_tracks/&lt;SAMPLE&gt;/</strong> ← Raw 3′/5′ bedGraphs & BigWigs<br>
├── 📁 <strong>04_counts/&lt;SAMPLE&gt;/</strong> ← Read count matrices<br>
├── 📁 <strong>05_normalized_tracks/&lt;SAMPLE&gt;/</strong> ← CPM & spike-in normalized tracks<br>
├── 📁 <strong>06_divergent_tx/&lt;SAMPLE&gt;/</strong> ← Bidirectional transcription analysis<br>
├── 📁 <strong>07_functional_regions/&lt;SAMPLE&gt;/</strong> ← Region assignments & summaries<br>
├── 📁 <strong>08_pol2_metrics/&lt;SAMPLE&gt;/</strong> ← RNA Pol II density & pausing metrics<br>
├── 📁 <strong>09_pol2_aggregate/</strong> ← Cross-sample Pol II comparisons<br>
├── 📁 <strong>10_qc/&lt;SAMPLE&gt;/</strong> ← Quality control metrics & stats<br>
└── 📁 <strong>11_reports/</strong> ← Interactive HTML reports (start here!)<br>
    ├── 📁 <strong>cohort/</strong> ← This cohort summary<br>
    └── 📁 <strong>samples/&lt;SAMPLE&gt;/</strong> ← Individual sample reports
        </div>
      </div>

      <div class="card">
        <h3>🔄 Reproducibility</h3>
        <p><strong>To reproduce this analysis:</strong></p>
        <div style="position: relative;">
          <div style="font-family: monospace; background: #f1f5f9; padding: 1rem; border-radius: 8px; margin: 1rem 0; position: relative;">
            <button onclick="copyToClipboard(this.nextElementSibling.textContent)" style="position: absolute; top: 0.5rem; right: 0.5rem; background: var(--accent); color: white; border: none; border-radius: 4px; padding: 0.25rem 0.5rem; cursor: pointer; font-size: 0.8em;">📋 Copy</button>
            <div style="display: none;">nextflow run main.nf -entry TrackTx -profile {args.profile} -resume</div>
# Resume with same parameters<br>
nextflow run main.nf -entry TrackTx -profile {args.profile} -resume<br><br>
# Fresh run with monitoring<br>
./run_pipeline.sh -profile {args.profile} -with-report -with-timeline
          </div>
        </div>
        <p><small class="muted">💡 <strong>Tip:</strong> Use <code>./run_pipeline.sh --validate-only</code> to test your environment before running.</small></p>
      </div>
    </section>

    <section id="sample-section" class="section" style="display:none">
      <h1>Per-sample</h1>
      <div id="sample-pages" class="grid"></div>
    </section>
  </main>
</div>

<script type="application/json" id="payload">{DATA_JSON}</script>
{JS}
"""

with open(args.out_html, "w", encoding="utf-8") as fh:
    fh.write(HTML)

print("INFO: cohort report written.", file=sys.stderr)
print(f"[combine_py] done ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)
