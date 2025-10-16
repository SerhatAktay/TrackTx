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

    # Regions list: [{"region": name, "reads": float, "region_count": int,
    #                 "region_length_total_bp": float, "region_length_median_bp": float}, ...]
    regions = js.get("regions") or []
    func_totals: Dict[str, float] = {}
    region_counts: Dict[str, int] = {}
    region_len_totals: Dict[str, float] = {}
    region_len_medians: Dict[str, float] = {}
    for it in regions:
        try:
            k = str(it.get("region"))
            v = float(it.get("reads", 0) or 0.0)
            c = int(it.get("region_count", 0) or 0)
            lt = float(it.get("region_length_total_bp", 0) or 0.0)
            lm = float(it.get("region_length_median_bp", 0) or 0.0)
            if k:
                func_totals[k] = func_totals.get(k, 0.0) + v
                region_counts[k] = region_counts.get(k, 0) + c
                region_len_totals[k] = region_len_totals.get(k, 0.0) + lt
                # For medians across multiple entries of the same label in one sample,
                # take median of medians conservatively by averaging if multiple present
                if k in region_len_medians and region_len_medians[k] != 0:
                    region_len_medians[k] = (region_len_medians[k] + lm) / 2.0
                else:
                    region_len_medians[k] = lm
        except Exception:
            pass
    # Separate true functional reads from non-localized polymerase
    unlocalized_reads = sum(v for k, v in func_totals.items() if _UNLOC_RE.search(str(k) or ""))
    localized_reads = sum(v for k, v in func_totals.items() if not _UNLOC_RE.search(str(k) or ""))
    reads_total_functional = localized_reads  # Exclude non-localized from functional total
    
    # Use explicit value from JSON if present, otherwise compute
    # Unlocalized fraction = unlocalized / (localized + unlocalized)
    unloc = js.get("unlocalized_fraction")
    if unloc is None:
        total_with_unloc = localized_reads + unlocalized_reads
        unloc = unlocalized_reads / total_with_unloc if total_with_unloc > 0 else None
    unloc = _to_float(unloc)

    # Extract nested fields from metrics and qc objects
    metrics = js.get("metrics") or {}
    qc = js.get("qc") or {}
    
    # Determine duplicate percentage - prefer UMI deduplication if available
    dup_percent = None
    if qc.get("umi_deduplication_enabled", False) and qc.get("umi_deduplication_percent") is not None:
        dup_percent = qc.get("umi_deduplication_percent")
    else:
        # Accept multiple synonyms
        dup_percent = (
            qc.get("duplicate_percent")
            or qc.get("duplicate_perc_of_total")
            or js.get("duplicate_percent")
            or js.get("duplicate_perc_of_total")
        )
    
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
        density_source=metrics.get("density_source"),
        density_reason=metrics.get("density_reason"),
        input_reads=qc.get("total_reads_raw") or js.get("input_reads"),
        dedup_reads=qc.get("dedup_reads_mapq_ge") or js.get("dedup_reads"),
        duplicate_percent=dup_percent,
        umi_deduplication_enabled=qc.get("umi_deduplication_enabled", False),
        umi_deduplication_percent=qc.get("umi_deduplication_percent"),
        unlocalized_fraction=unloc,
        func_totals=func_totals,
        region_counts=region_counts,
        region_length_totals=region_len_totals,
        region_length_medians=region_len_medians,
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
    # Include optional fields for downstream diagnostics
    r["density_source"] = s.get("density_source")
    r["density_reason"] = s.get("density_reason")
    for k in region_keys:
        r[f"func_{k}"] = (s.get("func_totals") or {}).get(k, 0)
        r[f"count_{k}"] = (s.get("region_counts") or {}).get(k, 0)
        r[f"len_total_{k}"] = (s.get("region_length_totals") or {}).get(k, 0)
        r[f"len_median_{k}"] = (s.get("region_length_medians") or {}).get(k, 0)
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
    columns=core_cols
            + [f"func_{k}" for k in region_keys]
            + [f"count_{k}" for k in region_keys]
            + [f"len_total_{k}" for k in region_keys]
            + [f"len_median_{k}" for k in region_keys],
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
  --bg:#f8f9fa; --fg:#1a1a1a; --muted:#666; --line:#ddd; --card:#fff;
  --ok:#10b981; --warn:#f59e0b; --fail:#ef4444; --accent:#2563eb;
}
@media (prefers-color-scheme: dark){
  :root{--bg:#0f1419; --fg:#e5e7eb; --muted:#999; --line:#2d333b; --card:#1a1f28;}
}
*{box-sizing:border-box;margin:0;padding:0}
body{background:var(--bg);color:var(--fg);font-family:-apple-system,system-ui,sans-serif;line-height:1.6;padding:1rem}
.container{max-width:1400px;margin:0 auto}
h1{font-size:2rem;margin-bottom:1rem}
h2{font-size:1.4rem;margin:2rem 0 1rem;border-bottom:2px solid var(--line);padding-bottom:0.5rem}
h3{font-size:1.1rem;margin:1.5rem 0 0.75rem}
.header{background:var(--card);border-radius:12px;padding:1.5rem;margin-bottom:1.5rem;border:1px solid var(--line)}
.kpi-grid{display:grid;gap:1rem;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));margin:1rem 0}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:1rem;text-align:center}
.kpi-label{font-size:0.75rem;color:var(--muted);text-transform:uppercase;margin-bottom:0.25rem}
.kpi-value{font-size:1.8rem;font-weight:700}
.status{display:inline-flex;align-items:center;gap:0.5rem;padding:0.5rem 1rem;border-radius:6px;font-weight:600;background:var(--card);border:1px solid var(--line)}
.status-dot{width:10px;height:10px;border-radius:50%}
.chart{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:1.5rem;margin:1rem 0;min-height:350px;position:relative;overflow:visible}
.chart-grid{display:grid;gap:1.5rem;grid-template-columns:repeat(auto-fit,minmax(500px,1fr))}
.table-container{background:var(--card);border:1px solid var(--line);border-radius:8px;overflow:hidden;margin:1rem 0}
table{width:100%;border-collapse:collapse;font-size:0.9rem}
th,td{padding:0.75rem;text-align:left;border-bottom:1px solid var(--line)}
th{background:var(--bg);font-weight:600;position:sticky;top:0;z-index:10}
.sticky{position:sticky;left:0;background:var(--card);z-index:5}
.filters{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:1rem;margin:1rem 0;display:grid;gap:1rem;grid-template-columns:repeat(auto-fit,minmax(200px,1fr))}
.filters input,.filters select{padding:0.5rem;border:1px solid var(--line);border-radius:4px;background:var(--bg);color:var(--fg);width:100%}
.btn{padding:0.5rem 1rem;border:1px solid var(--line);border-radius:4px;background:var(--card);color:var(--fg);cursor:pointer;font-size:0.9rem}
.btn:hover{background:var(--accent);color:#fff;border-color:var(--accent)}
.btn-primary{background:var(--accent);color:#fff;border-color:var(--accent)}
.tabs{display:flex;gap:0.5rem;margin:1rem 0;flex-wrap:wrap}
.tab{padding:0.5rem 1rem;border-radius:6px;cursor:pointer;background:var(--bg);border:1px solid var(--line)}
.tab.active{background:var(--accent);color:#fff;border-color:var(--accent)}
.section{display:none}
.section.active{display:block}
.muted{color:var(--muted);font-size:0.9rem}
.tip{position:absolute;background:var(--card);border:2px solid var(--accent);padding:0.75rem;border-radius:6px;font-size:0.9rem;pointer-events:none;opacity:0;z-index:1000;box-shadow:0 6px 16px rgba(0,0,0,0.2);transition:opacity 0.15s ease;max-width:250px;line-height:1.4}
@media (max-width:768px){
  .chart-grid{grid-template-columns:1fr}
  .kpi-grid{grid-template-columns:repeat(2,1fr)}
  body{padding:0.5rem}
}
</style>
"""

JS = r"""
<script>
const DATA = JSON.parse(document.getElementById('payload').textContent || '{}');
const $ = q => document.querySelector(q);
const $$ = q => Array.from(document.querySelectorAll(q));
const fmt = n => (n==null||isNaN(n))? 'n/a' : Number(n).toLocaleString();
const pct = x => (x==null||isNaN(x))? 'n/a' : (100*Number(x)).toFixed(1)+'%';
const median = arr => {const s=arr.sort((a,b)=>a-b); return s.length%2? s[Math.floor(s.length/2)]:(s[s.length/2-1]+s[s.length/2])/2;};
const labelNA = v => (v==null||v===''||v==='NA')? 'n/a' : String(v);

// ---------- Tab Navigation ----------
function switchTab(tab) {
  // Remove active from all tabs
  $$('.tab').forEach(t => t.classList.remove('active'));
  
  // Add active to clicked tab
  const clickedTab = Array.from($$('.tab')).find(t => t.textContent.toLowerCase().includes(tab.toLowerCase()));
  if (clickedTab) clickedTab.classList.add('active');
  
  // Navigate based on tab
  if (tab === 'overview') {
    location.hash = '#overview';
  } else if (tab === 'samples') {
    location.hash = '#samples';
  } else if (tab === 'regions') {
    location.hash = '#regions';
  } else if (tab === 'quality') {
    location.hash = '#qc';
  }
}

// ---------- Mobile sidebar toggle ----------
function toggleSidebar() {
  const sidebar = $('.sidebar');
  const overlay = $('.sidebar-overlay');
  if (sidebar) sidebar.classList.toggle('open');
  if (overlay) overlay.classList.toggle('open');
}

function sampleColorMap(){
  const ids = [...new Set(DATA.rows.map(r=>r.sample_id))].sort();
  const m={}; for (let i=0;i<ids.length;i++){ m[ids[i]] = `hsl(${Math.round(360*(i/Math.max(1,ids.length)))},70%,50%)`; }
  return sid => m[sid] || '#3b82f6';
}
const colorOfSample = sampleColorMap();

// Canonical colors from functional_regions.py (RGB values)
const REGION_COLORS = {
  'promoter':'rgb(243,132,0)','activepromoter':'rgb(243,132,0)','pppolii':'rgb(243,132,0)','pppol':'rgb(243,132,0)',
  'divergenttx':'rgb(178,59,212)','divtx':'rgb(178,59,212)','divergent':'rgb(178,59,212)','ppdiv':'rgb(178,59,212)',
  'enhancers':'rgb(115,212,122)','enhancer':'rgb(115,212,122)','enh':'rgb(115,212,122)',
  'genebody':'rgb(0,0,0)','gene body':'rgb(0,0,0)','body':'rgb(0,0,0)','gb':'rgb(0,0,0)',
  'cps':'rgb(103,200,249)','cleavagepolyadenylation':'rgb(103,200,249)',
  'shortgenes':'rgb(253,218,13)','short genes':'rgb(253,218,13)','sg':'rgb(253,218,13)',
  'terminationwindow':'rgb(255,54,98)','termination window':'rgb(255,54,98)','termination':'rgb(255,54,98)','tw':'rgb(255,54,98)'
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
  const W = Math.min(Math.max(containerWidth - 16, 320), 1200);
  const H = 500;  // Increased height for better label space
  const L = 70;   // Increased left margin for y-axis labels
  const B = 150;  // Much larger bottom margin for long rotated labels
  const T = 40;   // Increased top margin for title
  const R = 20;   // Right margin
  const gap = Math.max(2, W / 150);
  const n = labels.length || 1;
  const dataMax = Math.max(...values, 1);
  const max = dataMax * 1.1;  // Extend 10% past max for better visibility
  const bw = Math.max(8, Math.floor((W - L - R - (n - 1) * gap) / n));
  
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.style.width = '100%';
  svg.style.height = 'auto';
  
  const tip = document.createElement('div');
  tip.className = 'tip';
  el.style.position = 'relative';
  el.appendChild(tip);

  // Title
  const title = document.createElementNS(svg.namespaceURI,'text');
  title.setAttribute('x', W/2); title.setAttribute('y', 20);
  title.setAttribute('text-anchor', 'middle'); title.setAttribute('font-size', '16');
  title.setAttribute('font-weight', 'bold'); title.setAttribute('fill', 'currentColor');
  title.textContent = ylab;
  svg.appendChild(title);

  // grid + axis ticks
  for(let i=0;i<=5;i++){ 
    const y = T + (H-B-T)*i/5;
    const ln = document.createElementNS(svg.namespaceURI,'line');
    ln.setAttribute('x1', L); ln.setAttribute('x2', W-R); 
    ln.setAttribute('y1', y); ln.setAttribute('y2', y);
    ln.setAttribute('stroke', 'currentColor'); ln.setAttribute('opacity', '0.15'); 
    ln.setAttribute('stroke-dasharray', '4,4');
    svg.appendChild(ln);
    
    const t = document.createElementNS(svg.namespaceURI,'text'); 
    t.setAttribute('x', L-8); t.setAttribute('y', y+5); 
    t.setAttribute('font-size', '13'); t.setAttribute('text-anchor', 'end');
    t.setAttribute('fill', 'currentColor'); t.setAttribute('opacity', '0.7'); 
    t.textContent = Math.round(max*(1-i/5)).toLocaleString(); 
    svg.appendChild(t);
  }

  // Smart label display: show fewer labels if too many bars, but ensure readability
  const showEvery = n > 30 ? Math.ceil(n/15) : n > 15 ? Math.ceil(n/10) : 1;
  
  labels.forEach((lab, i) => {
    const v = values[i] || 0;
    const h = Math.round((v/max) * (H-B-T));
    const x = L + i*(bw+gap);
    const y = (H-B) - h;
    
    const r = document.createElementNS(svg.namespaceURI, 'rect');
    r.setAttribute('x', x); r.setAttribute('y', y);
    r.setAttribute('width', bw); r.setAttribute('height', h); r.setAttribute('rx', '3');
    r.setAttribute('fill', colors[i] || '#3b82f6');
    r.setAttribute('opacity', '0.9');
    r.setAttribute('stroke', 'rgba(255,255,255,0.3)');
    r.setAttribute('stroke-width', '0.5');
    if (onClick) r.style.cursor = 'pointer';
    
    r.addEventListener('mousemove', (e) => {
      tip.style.left = (e.offsetX+12) + 'px';
      tip.style.top = (e.offsetY-10) + 'px';
      tip.textContent = `${lab}: ${v.toLocaleString()}`;
      tip.style.opacity = '1';
    });
    r.addEventListener('mouseleave', () => tip.style.opacity = '0');
    if (onClick) r.addEventListener('click', () => onClick(lab));
    svg.appendChild(r);
    
    // Add labels with 90-degree rotation for maximum readability
    if (i % showEvery === 0 || n <= 10) {
      const t = document.createElementNS(svg.namespaceURI, 'text');
      t.setAttribute('x', x + bw/2);
      t.setAttribute('y', H - B + 12);
      t.setAttribute('text-anchor', 'end');
      t.setAttribute('font-size', '10');
      t.setAttribute('fill', 'currentColor');
      t.setAttribute('opacity', '0.8');
      t.setAttribute('transform', `rotate(-90 ${x+bw/2} ${H-B+12})`);
      // Truncate to fit in 140px (B margin minus some padding)
      t.textContent = lab.length > 18 ? lab.substring(0, 16) + '..' : lab;
      svg.appendChild(t);
    }
  });

  el.innerHTML = '';
  el.appendChild(svg);
}

function scatter(el, labels, xs, ys, colors, xl, yl, onClick){
  const containerWidth = el.parentElement.clientWidth || 600;
  const W = Math.min(Math.max(containerWidth - 16, 320), 1200);
  const H = 420;
  const L = 80;
  const B = 70;
  const T = 50;
  const R = 20;
  
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.style.width = '100%';
  svg.style.height = 'auto';
  
  const tip = document.createElement('div');
  tip.className = 'tip';
  el.style.position = 'relative';
  el.appendChild(tip);
  
  const maxX = Math.max(...xs.filter(v => v != null && !isNaN(v)), 1);
  const maxY = Math.max(...ys.filter(v => v != null && !isNaN(v)), 1);

  // Title
  const title = document.createElementNS(svg.namespaceURI, 'text');
  title.setAttribute('x', W/2); title.setAttribute('y', 25);
  title.setAttribute('text-anchor', 'middle'); title.setAttribute('font-size', '16');
  title.setAttribute('font-weight', 'bold'); title.setAttribute('fill', 'currentColor');
  title.textContent = `${yl} vs ${xl}`;
  svg.appendChild(title);

  // Grid lines
  for(let i=0; i<=5; i++){
    const y = T + (H-B-T)*i/5;
    const ln = document.createElementNS(svg.namespaceURI, 'line');
    ln.setAttribute('x1', L); ln.setAttribute('x2', W-R);
    ln.setAttribute('y1', y); ln.setAttribute('y2', y);
    ln.setAttribute('stroke', 'currentColor'); ln.setAttribute('opacity', '0.1');
    ln.setAttribute('stroke-dasharray', '3,3');
    svg.appendChild(ln);
    
    const t = document.createElementNS(svg.namespaceURI, 'text');
    t.setAttribute('x', L-8); t.setAttribute('y', y+4);
    t.setAttribute('font-size', '11'); t.setAttribute('text-anchor', 'end');
    t.setAttribute('fill', 'currentColor'); t.setAttribute('opacity', '0.7');
    t.textContent = (maxY * (1-i/5)).toFixed(1);
    svg.appendChild(t);
  }
  
  for(let i=0; i<=5; i++){
    const x = L + (W-L-R)*i/5;
    const ln = document.createElementNS(svg.namespaceURI, 'line');
    ln.setAttribute('x1', x); ln.setAttribute('x2', x);
    ln.setAttribute('y1', T); ln.setAttribute('y2', H-B);
    ln.setAttribute('stroke', 'currentColor'); ln.setAttribute('opacity', '0.1');
    ln.setAttribute('stroke-dasharray', '3,3');
    svg.appendChild(ln);
    
    const t = document.createElementNS(svg.namespaceURI, 'text');
    t.setAttribute('x', x); t.setAttribute('y', H-B+18);
    t.setAttribute('font-size', '11'); t.setAttribute('text-anchor', 'middle');
    t.setAttribute('fill', 'currentColor'); t.setAttribute('opacity', '0.7');
    t.textContent = (maxX * i/5).toFixed(1);
    svg.appendChild(t);
  }

  // axes
  const ax = document.createElementNS(svg.namespaceURI, 'line');
  ax.setAttribute('x1', L); ax.setAttribute('x2', W-R);
  ax.setAttribute('y1', H-B); ax.setAttribute('y2', H-B);
  ax.setAttribute('stroke', 'currentColor'); ax.setAttribute('opacity', '0.3');
  ax.setAttribute('stroke-width', '2');
  svg.appendChild(ax);
  
  const ay = document.createElementNS(svg.namespaceURI, 'line');
  ay.setAttribute('x1', L); ay.setAttribute('x2', L);
  ay.setAttribute('y1', T); ay.setAttribute('y2', H-B);
  ay.setAttribute('stroke', 'currentColor'); ay.setAttribute('opacity', '0.3');
  ay.setAttribute('stroke-width', '2');
  svg.appendChild(ay);

  labels.forEach((lab, i) => {
    const xVal = xs[i] ?? 0;
    const yVal = ys[i] ?? 0;
    if (isNaN(xVal) || isNaN(yVal)) return;
    
    const x = L + (xVal/maxX) * (W-L-R);
    const y = H-B - (yVal/maxY) * (H-B-T);
    
    const c = document.createElementNS(svg.namespaceURI, 'circle');
    c.setAttribute('cx', x); c.setAttribute('cy', y); c.setAttribute('r', 6);
    c.setAttribute('fill', colors[i] || '#3b82f6');
    c.setAttribute('opacity', '0.7');
    c.setAttribute('stroke', 'white');
    c.setAttribute('stroke-width', '1.5');
    c.style.cursor = 'pointer';
    
    c.addEventListener('mousemove', (e) => {
      tip.style.left = (e.offsetX+12) + 'px';
      tip.style.top = (e.offsetY-10) + 'px';
      tip.innerHTML = `<strong>${lab}</strong><br>${xl}: ${xVal.toFixed(2)}<br>${yl}: ${yVal.toFixed(2)}`;
      tip.style.opacity = '1';
    });
    c.addEventListener('mouseleave', () => tip.style.opacity = '0');
    if (onClick) c.addEventListener('click', () => onClick(lab));
    svg.appendChild(c);
  });
  
  // Axis labels
  const xt = document.createElementNS(svg.namespaceURI, 'text');
  xt.setAttribute('x', W/2); xt.setAttribute('y', H-20);
  xt.setAttribute('font-size', '14'); xt.setAttribute('text-anchor', 'middle');
  xt.setAttribute('fill', 'currentColor'); xt.setAttribute('font-weight', 'bold');
  xt.textContent = xl;
  svg.appendChild(xt);
  
  const yt = document.createElementNS(svg.namespaceURI, 'text');
  yt.setAttribute('x', 20); yt.setAttribute('y', H/2);
  yt.setAttribute('font-size', '14'); yt.setAttribute('text-anchor', 'middle');
  yt.setAttribute('transform', `rotate(-90 20,${H/2})`);
  yt.setAttribute('fill', 'currentColor'); yt.setAttribute('font-weight', 'bold');
  yt.textContent = yl;
  svg.appendChild(yt);
  
  el.innerHTML = '';
  el.appendChild(svg);
}

function histogram(el, values, xlabel, ylabel, bins){
  if (!values || values.length === 0) {
    el.innerHTML = '<div class="muted" style="padding:2rem;text-align:center;">No data available</div>';
    return;
  }
  
  const containerWidth = el.parentElement.clientWidth || 600;
  const W = Math.min(Math.max(containerWidth - 16, 320), 1200);
  const H = 400;
  const L = 80;
  const B = 70;
  const T = 50;
  const R = 20;
  
  // Calculate histogram bins
  const min = Math.min(...values);
  const max = Math.max(...values);
  const binWidth = (max - min) / bins || 1;
  const binCounts = new Array(bins).fill(0);
  const binEdges = [];
  const binLabels = [];  // Store sample IDs per bin
  
  for (let i = 0; i <= bins; i++) {
    binEdges.push(min + i * binWidth);
  }
  
  // Track which samples fall in which bin
  for (let i = 0; i < bins; i++) {
    binLabels[i] = [];
  }
  
  values.forEach((v, idx) => {
    const binIdx = Math.min(Math.floor((v - min) / binWidth), bins - 1);
    if (binIdx >= 0 && binIdx < bins) {
      binCounts[binIdx]++;
      // Store the index for potential sample name lookup
    }
  });
  
  const dataMax = Math.max(...binCounts, 1);
  const maxCount = dataMax * 1.1;  // Extend 10% past max
  
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.style.width = '100%';
  svg.style.height = 'auto';
  
  const tip = document.createElement('div');
  tip.className = 'tip';
  el.style.position = 'relative';
  el.appendChild(tip);
  
  // Title
  const title = document.createElementNS(svg.namespaceURI, 'text');
  title.setAttribute('x', W/2); title.setAttribute('y', 25);
  title.setAttribute('text-anchor', 'middle'); title.setAttribute('font-size', '16');
  title.setAttribute('font-weight', 'bold'); title.setAttribute('fill', 'currentColor');
  title.textContent = `${xlabel} Distribution`;
  svg.appendChild(title);
  
  // Grid
  for(let i=0; i<=5; i++){ 
    const y = T + (H-B-T)*i/5;
    const ln = document.createElementNS(svg.namespaceURI,'line');
    ln.setAttribute('x1', L); ln.setAttribute('x2', W-R);
    ln.setAttribute('y1', y); ln.setAttribute('y2', y);
    ln.setAttribute('stroke', 'currentColor'); ln.setAttribute('opacity', '0.15');
    ln.setAttribute('stroke-dasharray', '4,4');
    svg.appendChild(ln);
    
    const t = document.createElementNS(svg.namespaceURI,'text');
    t.setAttribute('x', L-8); t.setAttribute('y', y+4);
    t.setAttribute('font-size', '11'); t.setAttribute('text-anchor', 'end');
    t.setAttribute('fill', 'currentColor'); t.setAttribute('opacity', '0.7');
    t.textContent = Math.round(maxCount*(1-i/5));
    svg.appendChild(t);
  }
  
  // Bars
  const barWidth = (W - L - R) / bins - 2;
  binCounts.forEach((count, i) => {
    const x = L + i * ((W - L - R) / bins);
    const h = (count / maxCount) * (H - B - T);
    const y = H - B - h;
    
    const rect = document.createElementNS(svg.namespaceURI, 'rect');
    rect.setAttribute('x', x);
    rect.setAttribute('y', y);
    rect.setAttribute('width', barWidth);
    rect.setAttribute('height', h);
    rect.setAttribute('fill', '#3b82f6');
    rect.setAttribute('opacity', '0.8');
    rect.setAttribute('rx', '2');
    
    rect.addEventListener('mousemove', (e) => {
      tip.style.left = (e.offsetX+12) + 'px';
      tip.style.top = (e.offsetY-10) + 'px';
      tip.innerHTML = `<strong>Range:</strong> ${binEdges[i].toFixed(1)} - ${binEdges[i+1].toFixed(1)}<br><strong>Count:</strong> ${count}`;
      tip.style.opacity = '1';
    });
    rect.addEventListener('mouseleave', () => tip.style.opacity = '0');
    
    svg.appendChild(rect);
  });
  
  // X-axis labels
  for(let i=0; i<=5; i++){
    const x = L + (W-L-R)*i/5;
    const val = min + (max-min)*i/5;
    const t = document.createElementNS(svg.namespaceURI,'text');
    t.setAttribute('x', x); t.setAttribute('y', H-B+20);
    t.setAttribute('font-size', '10'); t.setAttribute('text-anchor', 'middle');
    t.setAttribute('fill', 'currentColor'); t.setAttribute('opacity', '0.7');
    t.textContent = val.toFixed(1);
    svg.appendChild(t);
  }
  
  // Axis labels
  const xl = document.createElementNS(svg.namespaceURI, 'text');
  xl.setAttribute('x', W/2); xl.setAttribute('y', H-20);
  xl.setAttribute('font-size', '13'); xl.setAttribute('text-anchor', 'middle');
  xl.setAttribute('fill', 'currentColor'); xl.setAttribute('font-weight', 'bold');
  xl.textContent = xlabel;
  svg.appendChild(xl);
  
  const yl = document.createElementNS(svg.namespaceURI, 'text');
  yl.setAttribute('x', 20); yl.setAttribute('y', H/2);
  yl.setAttribute('font-size', '13'); yl.setAttribute('text-anchor', 'middle');
  yl.setAttribute('transform', `rotate(-90 20,${H/2})`);
  yl.setAttribute('fill', 'currentColor'); yl.setAttribute('font-weight', 'bold');
  yl.textContent = ylabel;
  svg.appendChild(yl);
  
  el.innerHTML = '';
  el.appendChild(svg);
}

function donut(el, labels, values, colorFor){
  const containerWidth = el.parentElement.clientWidth || 600;
  const W = Math.min(Math.max(containerWidth - 40, 500), 700);
  const chartH = 300;
  const legendH = Math.max(80, Math.ceil(labels.length / 3) * 25);
  const H = chartH + legendH + 40;
  const R = Math.min(120, W / 4);
  const cx = W / 2;
  const cy = chartH / 2 + 30;
  const rIn = R * 0.55;
  
  const total = values.reduce((a,b) => a + b, 0) || 1;
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.style.width = '100%';
  svg.style.height = 'auto';
  
  const tip = document.createElement('div');
  tip.className = 'tip';
  el.style.position = 'relative';
  el.appendChild(tip);
  
  // Title
  const title = document.createElementNS(svg.namespaceURI, 'text');
  title.setAttribute('x', W/2); title.setAttribute('y', 20);
  title.setAttribute('text-anchor', 'middle'); title.setAttribute('font-size', '16');
  title.setAttribute('font-weight', 'bold'); title.setAttribute('fill', 'currentColor');
  title.textContent = 'Region Composition';
  svg.appendChild(title);
  
  let a0 = -Math.PI/2;
  values.forEach((v, i) => {
    const a1 = a0 + (v/total) * 2 * Math.PI;
    const x0 = cx + R*Math.cos(a0), y0 = cy + R*Math.sin(a0);
    const x1 = cx + R*Math.cos(a1), y1 = cy + R*Math.sin(a1);
    const large = (a1-a0) > Math.PI ? 1 : 0;
    const path = document.createElementNS(svg.namespaceURI, 'path');
    path.setAttribute('d', `M ${cx+rIn*Math.cos(a0)} ${cy+rIn*Math.sin(a0)} L ${x0} ${y0} A ${R} ${R} 0 ${large} 1 ${x1} ${y1} L ${cx+rIn*Math.cos(a1)} ${cy+rIn*Math.sin(a1)} A ${rIn} ${rIn} 0 ${large} 0 ${cx+rIn*Math.cos(a0)} ${cy+rIn*Math.sin(a0)} Z`);
    path.setAttribute('fill', colorFor(labels[i]));
    path.setAttribute('opacity', '0.92');
    path.setAttribute('stroke', 'white');
    path.setAttribute('stroke-width', '2');
    path.style.cursor = 'pointer';
    
    path.addEventListener('mouseenter', function() {
      this.setAttribute('opacity', '1');
      this.setAttribute('transform', `translate(${5*Math.cos((a0+a1)/2)},${5*Math.sin((a0+a1)/2)})`);
    });
    path.addEventListener('mouseleave', function() {
      this.setAttribute('opacity', '0.92');
      this.setAttribute('transform', '');
    });
    
    path.addEventListener('mousemove', (e) => {
      tip.style.left = (e.offsetX+12) + 'px';
      tip.style.top = (e.offsetY-10) + 'px';
      tip.innerHTML = `<strong>${labels[i]}</strong><br>${values[i].toLocaleString()} reads<br>${((values[i]/total)*100).toFixed(1)}%`;
      tip.style.opacity = '1';
    });
    path.addEventListener('mouseleave', () => tip.style.opacity = '0');
    
    svg.appendChild(path);
    a0 = a1;
  });
  
  // Center label
  const centerText = document.createElementNS(svg.namespaceURI, 'text');
  centerText.setAttribute('x', cx); centerText.setAttribute('y', cy-5);
  centerText.setAttribute('text-anchor', 'middle'); centerText.setAttribute('font-size', '13');
  centerText.setAttribute('font-weight', 'bold'); centerText.setAttribute('fill', 'currentColor');
  centerText.textContent = 'Total';
  svg.appendChild(centerText);
  
  const centerVal = document.createElementNS(svg.namespaceURI, 'text');
  centerVal.setAttribute('x', cx); centerVal.setAttribute('y', cy+12);
  centerVal.setAttribute('text-anchor', 'middle'); centerVal.setAttribute('font-size', '15');
  centerVal.setAttribute('font-weight', 'bold'); centerVal.setAttribute('fill', 'currentColor');
  centerVal.textContent = total.toLocaleString();
  svg.appendChild(centerVal);
  
  // Legend below
  const legendY = chartH + 10;
  const itemsPerRow = 3;
  const itemWidth = W / itemsPerRow;
  labels.forEach((lab, i) => {
    const row = Math.floor(i / itemsPerRow);
    const col = i % itemsPerRow;
    const x = col * itemWidth + 15;
    const y = legendY + row * 25;
    
    const rect = document.createElementNS(svg.namespaceURI, 'rect');
    rect.setAttribute('x', x); rect.setAttribute('y', y);
    rect.setAttribute('width', 14); rect.setAttribute('height', 14);
    rect.setAttribute('rx', 2);
    rect.setAttribute('fill', colorFor(lab));
    svg.appendChild(rect);
    
    const text = document.createElementNS(svg.namespaceURI, 'text');
    text.setAttribute('x', x+20); text.setAttribute('y', y+11);
    text.setAttribute('font-size', '11'); text.setAttribute('fill', 'currentColor');
    text.textContent = `${lab} (${((values[i]/total)*100).toFixed(1)}%)`;
    svg.appendChild(text);
  });
  
  el.innerHTML = '';
  el.appendChild(svg);
}

function regionLegend(el, labels){
  if (!el) return;
  el.innerHTML = '<div style="background: var(--card); padding: 1rem; border-radius: 8px; border: 1px solid var(--line);"><h3 style="margin: 0 0 0.75rem; font-size: 1.1rem;">Region Color Legend</h3><div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap: 0.5rem;">' + 
    labels.map(l=>`<span class="pill" style="display:flex; align-items:center; gap:8px; padding:6px 10px; background: rgba(0,0,0,0.02); border-radius:6px; border:1px solid var(--line);"><span style="display:inline-block;width:16px;height:16px;border-radius:3px;background:${regionColor(l)};flex-shrink:0;"></span><span style="font-size:0.9rem;font-weight:500;">${l}</span></span>`).join('') +
    '</div></div>';
}

function toPercent(vals){
  const total = vals.reduce((a,b)=>a+b,0) || 1;
  return vals.map(v=> (100.0*v/total));
}

// Grouped bar chart for comparing two series across labels
function groupedBarChart(el, labels, series, colorsA, colorsB, ylab, legends=['A','B']){
  const containerWidth = el.parentElement.clientWidth || 600;
  const W = Math.min(Math.max(containerWidth - 40, 400), 1200);
  const H = 520, L = 80, B = 150, T = 50, R = 20;
  const n = labels.length || 1;
  const max = Math.max(...series[0], ...series[1], 1);
  const groupWidth = Math.max(24, Math.floor((W - L - R) / Math.max(1, n)) - 12);
  const bw = Math.max(10, Math.floor((groupWidth - 8) / 2));
  const gap = Math.max(8, Math.floor((W - L - R - n*groupWidth) / Math.max(1,n-1)));
  
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.style.width = '100%'; svg.style.height = 'auto';
  
  const tip = document.createElement('div'); tip.className = 'tip'; 
  el.style.position='relative'; el.appendChild(tip);
  
  // Title
  const title = document.createElementNS(svg.namespaceURI,'text');
  title.setAttribute('x', W/2); title.setAttribute('y', 25);
  title.setAttribute('text-anchor', 'middle'); title.setAttribute('font-size', '16');
  title.setAttribute('font-weight', 'bold'); title.setAttribute('fill', 'currentColor');
  title.textContent = ylab;
  svg.appendChild(title);
  
  // Grid
  for(let i=0;i<=5;i++){ 
    const y = T + (H-B-T)*i/5;
    const ln = document.createElementNS(svg.namespaceURI,'line');
    ln.setAttribute('x1', L); ln.setAttribute('x2', W-R);
    ln.setAttribute('y1', y); ln.setAttribute('y2', y);
    ln.setAttribute('stroke', 'currentColor'); ln.setAttribute('opacity', '0.15');
    ln.setAttribute('stroke-dasharray', '4,4');
    svg.appendChild(ln);
    
    const t = document.createElementNS(svg.namespaceURI,'text');
    t.setAttribute('x', L-8); t.setAttribute('y', y+4);
    t.setAttribute('font-size', '12'); t.setAttribute('text-anchor', 'end');
    t.setAttribute('fill', 'currentColor'); t.setAttribute('opacity', '0.7');
    t.textContent = Math.round(max*(1-i/5)).toLocaleString();
    svg.appendChild(t);
  }
  
  // Legend
  const leg1 = document.createElementNS(svg.namespaceURI, 'rect');
  leg1.setAttribute('x', W-R-150); leg1.setAttribute('y', T-25);
  leg1.setAttribute('width', 12); leg1.setAttribute('height', 12);
  leg1.setAttribute('fill', colorsA[0] || '#3b82f6');
  svg.appendChild(leg1);
  
  const txt1 = document.createElementNS(svg.namespaceURI, 'text');
  txt1.setAttribute('x', W-R-135); txt1.setAttribute('y', T-16);
  txt1.setAttribute('font-size', '11'); txt1.setAttribute('fill', 'currentColor');
  txt1.textContent = legends[0];
  svg.appendChild(txt1);
  
  const leg2 = document.createElementNS(svg.namespaceURI, 'rect');
  leg2.setAttribute('x', W-R-70); leg2.setAttribute('y', T-25);
  leg2.setAttribute('width', 12); leg2.setAttribute('height', 12);
  leg2.setAttribute('fill', colorsB[0] || '#10b981');
  svg.appendChild(leg2);
  
  const txt2 = document.createElementNS(svg.namespaceURI, 'text');
  txt2.setAttribute('x', W-R-55); txt2.setAttribute('y', T-16);
  txt2.setAttribute('font-size', '11'); txt2.setAttribute('fill', 'currentColor');
  txt2.textContent = legends[1];
  svg.appendChild(txt2);
  
  labels.forEach((lab, i) => {
    const x0 = L + i*(groupWidth+gap);
    const vA = series[0][i]||0, vB = series[1][i]||0;
    const hA = Math.round((vA/max)*(H-B-T));
    const hB = Math.round((vB/max)*(H-B-T));
    const yA = (H-B)-hA, yB = (H-B)-hB;
    
    const rA = document.createElementNS(svg.namespaceURI,'rect');
    rA.setAttribute('x', x0); rA.setAttribute('y', yA);
    rA.setAttribute('width', bw); rA.setAttribute('height', hA);
    rA.setAttribute('rx', '3'); rA.setAttribute('fill', colorsA[i]||'#3b82f6');
    rA.setAttribute('opacity', '0.9');
    rA.addEventListener('mousemove', (e) => {
      tip.style.left = (e.offsetX+12)+'px'; tip.style.top = (e.offsetY-10)+'px';
      tip.innerHTML = `<strong>${lab}</strong><br>${legends[0]}: ${vA.toLocaleString()}`;
      tip.style.opacity = '1';
    });
    rA.addEventListener('mouseleave', () => tip.style.opacity='0');
    svg.appendChild(rA);
    
    const rB = document.createElementNS(svg.namespaceURI,'rect');
    rB.setAttribute('x', x0+bw+4); rB.setAttribute('y', yB);
    rB.setAttribute('width', bw); rB.setAttribute('height', hB);
    rB.setAttribute('rx', '3'); rB.setAttribute('fill', colorsB[i]||'#10b981');
    rB.setAttribute('opacity', '0.9');
    rB.addEventListener('mousemove', (e) => {
      tip.style.left = (e.offsetX+12)+'px'; tip.style.top = (e.offsetY-10)+'px';
      tip.innerHTML = `<strong>${lab}</strong><br>${legends[1]}: ${vB.toLocaleString()}`;
      tip.style.opacity = '1';
    });
    rB.addEventListener('mouseleave', () => tip.style.opacity='0');
    svg.appendChild(rB);
    
    const t = document.createElementNS(svg.namespaceURI,'text');
    t.setAttribute('x', x0+groupWidth/2); t.setAttribute('y', H-B+12);
    t.setAttribute('text-anchor', 'end');
    t.setAttribute('font-size', '10');
    t.setAttribute('fill', 'currentColor'); t.setAttribute('opacity', '0.8');
    t.setAttribute('transform', `rotate(-90 ${x0+groupWidth/2} ${H-B+12})`);
    t.textContent = lab.length > 18 ? lab.substring(0, 16) + '..' : lab;
    svg.appendChild(t);
  });
  
  el.innerHTML = ''; el.appendChild(svg);
}

function renderConditionPies(rows){
  const rk = DATA.region_keys||[];
  const conds = [...new Set(rows.map(r=>r.condition).filter(Boolean))].sort();
  const wrap = document.getElementById('cond-pies'); if (!wrap) return;
  wrap.innerHTML = '';
  conds.forEach(cond=>{
    const subset = rows.filter(r=>r.condition===cond);
    const totals = rk.map(k=> subset.reduce((a,r)=> a + (Number(r['func_'+k])||0), 0));
    const card = document.createElement('div'); card.className='card';
    const title = document.createElement('h3'); title.textContent = `Condition: ${cond}`; title.style.margin = '0 0 8px'; card.appendChild(title);
    const div = document.createElement('div'); div.className='chart'; card.appendChild(div);
    wrap.appendChild(card);
    donut(div, rk, totals, s=>regionColor(s));
  });
}

function renderConditionCompare(rows){
  const rk = DATA.region_keys||[];
  const a = (document.getElementById('cmp-cond-a')||{}).value || '';
  const b = (document.getElementById('cmp-cond-b')||{}).value || '';
  const el = document.getElementById('chart-cond-compare'); if (!el) return;
  if (!a || !b || a===b){ el.innerHTML = '<div class="muted">Select two different conditions to compare.</div>'; return; }
  const rowsA = rows.filter(r=>r.condition===a);
  const rowsB = rows.filter(r=>r.condition===b);
  const valsA = rk.map(k=> rowsA.reduce((s,r)=> s + (Number(r['func_'+k])||0), 0));
  const valsB = rk.map(k=> rowsB.reduce((s,r)=> s + (Number(r['func_'+k])||0), 0));
  groupedBarChart(el, rk, [valsA, valsB], rk.map(regionColor), rk.map(regionColor), 'Reads by region', [a,b]);
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
  const unlocs = rows.map(r=>Number(r.unlocalized_fraction)||0);
  const dups   = rows.map(r=>Number(r.duplicate_percent)||0);

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
  regionLegend($('#region-legend'), rk);

  // New: Region-level bars across cohort
  const renderRegionCharts = () => {
    const regionReads = rk.map(k => (DATA.region_totals||{})[k] || 0);
    const pctToggle = document.getElementById('toggle-region-percent');
    const usePct = pctToggle && pctToggle.checked;
    const dispReads = usePct ? toPercent(regionReads) : regionReads;
    barChart($('#chart-reg-reads'), rk, dispReads, rk.map(regionColor), usePct? 'Reads by region (%)' : 'Total reads by region', null);
    
    const regionLenTotals = rk.map(k => rows.reduce((a,r)=> a + (Number(r['len_total_'+k])||0), 0));
    const dispLen = usePct ? toPercent(regionLenTotals) : regionLenTotals;
    barChart($('#chart-reg-lentot'), rk, dispLen, rk.map(regionColor), usePct? 'Total length by region (%)' : 'Total length (bp) by region', null);
  };
  
  renderRegionCharts();
  
  // Wire up the toggle
  const pctToggle = document.getElementById('toggle-region-percent');
  if (pctToggle) {
    pctToggle.addEventListener('change', renderRegionCharts);
  }

  // Sum counts across samples per region
  const regionCounts = rk.map(k => rows.reduce((a,r)=> a + (Number(r['count_'+k])||0), 0));
  barChart($('#chart-reg-counts'), rk, regionCounts, rk.map(regionColor), 'Total region count (cohort)', null);

  // (moved into renderRegionCharts function above)

  // Median of median lengths across samples per region
  const regionLenMedians = rk.map(k => {
    const vals = rows.map(r => Number(r['len_median_'+k])||0).filter(v=>v>0).sort((a,b)=>a-b);
    if (!vals.length) return 0;
    return vals.length%2? vals[Math.floor(vals.length/2)] : (vals[vals.length/2-1]+vals[vals.length/2])/2;
  });
  barChart($('#chart-reg-lenmed'), rk, regionLenMedians, rk.map(regionColor), 'Median length (bp) by region', null);

  // New scatters: Functional reads vs Unlocalized; Duplicate% vs Unlocalized
  scatter($('#chart-unloc-vs-reads'), labels, reads, unlocs.map(x=>x*100), colors, 'Reads (functional)', 'Unlocalized %', sid => location.hash = '#sample/'+encodeURIComponent(sid));
  scatter($('#chart-dup-vs-unloc'), labels, dups, unlocs.map(x=>x*100), colors, 'Duplicate %', 'Unlocalized %', sid => location.hash = '#sample/'+encodeURIComponent(sid));
  
  // Histograms for duplicate and unlocalized percentages
  histogram($('#chart-hist-dup'), dups.filter(v=>!isNaN(v)&&v>0), 'Duplicate %', 'Number of samples', 20);
  histogram($('#chart-hist-unloc'), unlocs.filter(v=>!isNaN(v)&&v>0).map(x=>x*100), 'Unlocalized %', 'Number of samples', 20);
  
  // CPM factor vs reads scatter
  const cpmFactors = rows.map(r=>Number(r.cpm_factor)||0).filter(v=>v>0);
  scatter($('#chart-cpm-vs-reads'), 
          rows.filter(r=>Number(r.cpm_factor)>0).map(r=>r.sample_id),
          cpmFactors,
          rows.filter(r=>Number(r.cpm_factor)>0).map(r=>Number(r.reads_total_functional)||0),
          rows.filter(r=>Number(r.cpm_factor)>0).map(r=>colorOfSample(r.sample_id)),
          'CPM Normalization Factor', 'Functional Reads',
          sid => location.hash = '#sample/'+encodeURIComponent(sid));

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
  populateChartFilters(rows);
  renderConditionPies(rows);
  const setCmpOptions = ()=>{
    const aSel = document.getElementById('cmp-cond-a');
    const bSel = document.getElementById('cmp-cond-b');
    if (!aSel || !bSel) return;
    const conds = [...new Set(rows.map(r=>r.condition).filter(Boolean))].sort();
    aSel.innerHTML = '<option value="">Select A</option>' + conds.map(c=>`<option value="${c}">${c}</option>`).join('');
    bSel.innerHTML = '<option value="">Select B</option>' + conds.map(c=>`<option value="${c}">${c}</option>`).join('');
  };
  setCmpOptions();
  const cmpBtn = document.getElementById('cmp-render'); if (cmpBtn) cmpBtn.onclick = ()=> renderConditionCompare(rows);
  
  // Initialize table info
  const info = $('#table-info');
  if (info) info.textContent = `${rows.length} samples`;

  // Missing inputs / n/a reasons
  const missingDiv = document.getElementById('missing-inputs');
  if (missingDiv){
    const issues = [];
    rows.forEach(r=>{
      if (r.median_density==null || isNaN(r.median_density)){
        issues.push(`• ${r.sample_id}: median density n/a${r.density_reason? ' — '+r.density_reason : ''}`);
      }
      if (r.duplicate_percent==null || isNaN(r.duplicate_percent)){
        issues.push(`• ${r.sample_id}: duplicate% n/a — missing in qc_pol2.json`);
      }
    });
    missingDiv.innerHTML = issues.length? issues.join('<br>') : 'No issues detected.';
  }
  const missingQC = document.getElementById('missing-inputs-qc');
  if (missingQC){
    missingQC.innerHTML = (document.getElementById('missing-inputs')||{}).innerHTML || '';
  }

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

  // Group visibility helpers
  const REGION_IDS = ['region-legend','chart-reg-reads','chart-reg-counts','chart-reg-lentot','chart-reg-lenmed','cond-pies','chart-cond-compare'];
  const QC_IDS = ['chart-unloc-vs-reads','chart-dup-vs-unloc','chart-hist-dup','chart-hist-unloc','chart-cpm-vs-reads','missing-inputs','missing-inputs-qc','chart-cond','chart-tp','chart-reset'];
  function setGroupVisible(ids, visible){
    ids.forEach(id=>{ const el=document.getElementById(id); if(!el) return; const card=el.closest('.card')||el; card.style.display = visible? '' : 'none'; });
  }

  // Routing
  function show(hash){
    const h = (hash || '#overview');
    $$('.section').forEach(s=> s.style.display='none');
    if (h.startsWith('#sample/')){
      $('#overview').style.display='none';
      $('#sample-section').style.display='block';
      const id = decodeURIComponent(h.split('/')[1]||'');
      const target = document.getElementById('sample/'+encodeURIComponent(id));
      if (target){ target.scrollIntoView({behavior:'smooth',block:'start'}); }
      // hide both groups on sample view
      setGroupVisible(REGION_IDS, false); setGroupVisible(QC_IDS, false);
    } else if (h === '#regions'){
      $('#overview').style.display='block';
      // show regions, hide qc
      setGroupVisible(REGION_IDS, true); setGroupVisible(QC_IDS, false);
      const el = document.getElementById('chart-reg-reads') || document.getElementById('region-legend') || document.getElementById('chart-pie');
      if (el) el.scrollIntoView({behavior:'smooth', block:'start'});
    } else if (h === '#qc'){
      $('#overview').style.display='block';
      // show qc, hide regions
      setGroupVisible(REGION_IDS, false); setGroupVisible(QC_IDS, true);
      const el = document.getElementById('chart-unloc-vs-reads') || document.getElementById('chart-hist-dup');
      if (el) el.scrollIntoView({behavior:'smooth', block:'start'});
    } else if (h === '#samples'){
      $('#overview').style.display='block';
      setGroupVisible(REGION_IDS, false); setGroupVisible(QC_IDS, false);
      const el = document.getElementById('table-body');
      if (el) el.scrollIntoView({behavior:'smooth', block:'start'});
    } else {
      $('#sample-section').style.display='none';
      $('#overview').style.display='block';
      // default: show a light overview (hide extra groups)
      setGroupVisible(REGION_IDS, false); setGroupVisible(QC_IDS, false);
    }
  }
  window.addEventListener('hashchange', ()=> show(location.hash));
  if (!location.hash) location.hash = '#overview';
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

function populateChartFilters(rows) {
  // Populate condition and timepoint filters for charts
  const condSel = document.getElementById('chart-cond');
  const tpSel = document.getElementById('chart-tp');
  
  if (condSel) {
    const conditions = [...new Set(rows.map(r => r.condition))].filter(c => c && c !== 'NA').sort();
    conditions.forEach(c => {
      const opt = document.createElement('option');
      opt.value = c;
      opt.textContent = c;
      condSel.appendChild(opt);
    });
  }
  
  if (tpSel) {
    const timepoints = [...new Set(rows.map(r => r.timepoint))].filter(t => t && t !== 'NA').sort();
    timepoints.forEach(t => {
      const opt = document.createElement('option');
      opt.value = t;
      opt.textContent = t;
      tpSel.appendChild(opt);
    });
  }
}

function filteredRowsForCharts(rows) {
  const condSel = document.getElementById('chart-cond');
  const tpSel = document.getElementById('chart-tp');
  
  let filtered = rows;
  if (condSel && condSel.value) {
    filtered = filtered.filter(r => r.condition === condSel.value);
  }
  if (tpSel && tpSel.value) {
    filtered = filtered.filter(r => r.timepoint === tpSel.value);
  }
  return filtered;
}

function renderCharts(rows) {
  // Re-render main charts with filtered data
  if (!rows || rows.length === 0) return;
  
  const labels = rows.map(r => r.sample_id);
  const colors = labels.map(colorOfSample);
  const reads = rows.map(r => Number(r.reads_total_functional) || 0);
  const piMed = rows.map(r => Number(r.median_pausing_index) || 0);
  const dens = rows.map(r => Number(r.median_density) || 0);
  
  barChart($('#chart-reads'), labels, reads, colors, 'Reads in functional regions', sid => location.hash = '#sample/' + encodeURIComponent(sid));
  scatter($('#chart-scatter'), labels, piMed, dens, colors, 'Median Pausing Index', 'Median Density', sid => location.hash = '#sample/' + encodeURIComponent(sid));
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
  // After build, wire up chart filters and initial render
  const rows = (DATA && DATA.rows) ? DATA.rows : [];
  if (rows && rows.length){
    // Delay to ensure DOM is ready
    setTimeout(()=>{
      if (typeof populateChartFilters==='function'){
        populateChartFilters(rows);
      }
      if (typeof renderCharts==='function'){
        const render = ()=> renderCharts(filteredRowsForCharts(rows));
        const condSel = document.getElementById('chart-cond');
        const tpSel   = document.getElementById('chart-tp');
        const rst     = document.getElementById('chart-reset');
        const pctTgl  = document.getElementById('toggle-region-percent');
        if (condSel) condSel.onchange = render;
        if (tpSel)   tpSel.onchange   = render;
        if (rst)     rst.onclick      = ()=>{ condSel && (condSel.value=''); tpSel && (tpSel.value=''); render(); };
        if (pctTgl)  pctTgl.onchange  = render;
        render();
      }
    }, 0);
  }
});
</script>
"""

HTML = f"""<!doctype html>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>TrackTx — Cohort Summary</title>
{CSS}
<div class="container">
  <div class="header">
    <h1>TrackTx Cohort Summary</h1>
    <div class="muted">{args.run_name} • v{args.pipeline_version} • {args.profile} • {args.duration}</div>
  </div>

  <div class="tabs">
    <div class="tab active" onclick="switchTab('overview')">Overview</div>
    <div class="tab" onclick="switchTab('samples')">Samples</div>
    <div class="tab" onclick="switchTab('regions')">Regions</div>
    <div class="tab" onclick="switchTab('quality')">Quality</div>
  </div>

  <div class="section active" id="overview">
      <!-- Hero Banner -->
      <div class="hero" style="background: linear-gradient(135deg, rgba(37,99,235,0.05) 0%, rgba(37,99,235,0.1) 50%, rgba(147,51,234,0.05) 100%); padding: 2.5rem; margin-bottom: 2rem; border: 2px solid var(--line);">
        <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 1.5rem;">
          <div>
            <h1 style="margin: 0 0 0.5rem; font-size: 2.5rem; background: linear-gradient(135deg, #2563eb, #9333ea); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;">Cohort Analysis Summary</h1>
            <p style="margin: 0; font-size: 1.1rem; color: var(--muted);">
              <strong>{args.run_name}</strong> • TrackTx v{args.pipeline_version} • {args.profile} profile
            </p>
            <p style="margin: 0.5rem 0 0; font-size: 0.9rem; color: var(--muted);">
              Generated {datetime.datetime.now().strftime("%B %d, %Y at %H:%M")} • Duration: {args.duration}
            </p>
          </div>
          <div id="badge"></div>
        </div>
      </div>

      <!-- Key Performance Indicators -->

      <div class="kpi-grid">
        <div class="kpi"><div class="lab">Total Samples</div><div id="kpi-n" class="val"></div></div>
        <div class="kpi"><div class="lab">Functional Reads</div><div id="kpi-reads" class="val"></div></div>
        <div class="kpi"><div class="lab">Median Duplicate %</div><div id="kpi-dup" class="val"></div></div>
        <div class="kpi"><div class="lab">Region Types</div><div id="kpi-reg" class="val"></div></div>
        <div class="kpi"><div class="lab">Median Unlocalized</div><div id="kpi-unloc" class="val"></div></div>
        <div class="kpi"><div class="lab">Quality Status</div><div id="kpi-quality" class="val"></div></div>
      </div>

      <!-- Section Divider -->
      <div style="border-top: 3px solid var(--line); margin: 3rem 0 2rem; padding-top: 1rem;">
        <h2 style="font-size: 1.75rem; margin-bottom: 0.5rem;">📊 Quality Control Overview</h2>
        <p class="muted" style="margin: 0;">Understand your data quality at a glance</p>
      </div>
      
      <div class="card">
        <h3>💡 What do these metrics mean?</h3>
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
          <div>
            <strong>Median Density:</strong> Median CPM/signal across functional regions. Source shown in the KPI. If n/a, hover the KPI to see the reason.
          </div>
        </div>
      </div>

      <div class="card">
        <h3>Missing Inputs / n/a Reasons</h3>
        <div id="missing-inputs" class="muted">Computed from sample metadata at load.</div>
      </div>

      <!-- Section Divider -->
      <div style="border-top: 3px solid var(--line); margin: 3rem 0 2rem; padding-top: 1rem;">
        <h2 id="regions" style="font-size: 1.75rem; margin-bottom: 0.5rem;">📈 Interactive Visualizations</h2>
        <p class="muted" style="margin: 0;">Explore data patterns through dynamic charts</p>
      </div>
      
      <div class="chart-grid">
        <div class="card"><div id="chart-reads" class="chart"></div></div>
        <div class="card"><div id="chart-scatter" class="chart"></div></div>
        <div class="card"><div id="chart-pie" class="chart"></div></div>
        <div class="card"><div id="region-legend"></div></div>
        <div class="card"><div id="chart-unloc-vs-reads" class="chart"></div></div>
        <div class="card"><div id="chart-dup-vs-unloc" class="chart"></div></div>
        <div class="card">
          <div style="display:flex; gap:8px; align-items:center; margin-bottom:6px;">
            <strong>Facet charts by:</strong>
            <select id="chart-cond" style="padding:4px 6px; border:1px solid var(--line); background: var(--card); color: var(--fg);">
              <option value="">All Conditions</option>
            </select>
            <select id="chart-tp" style="padding:4px 6px; border:1px solid var(--line); background: var(--card); color: var(--fg);">
              <option value="">All Timepoints</option>
            </select>
            <button id="chart-reset" style="margin-left:auto; padding:4px 8px; border:1px solid var(--line); background:var(--card); color:var(--fg); cursor:pointer;">Reset</button>
          </div>
          <div id="chart-hist-dup" class="chart"></div>
        </div>
        <div class="card"><div id="chart-hist-unloc" class="chart"></div></div>
        <div class="card"><div id="chart-cpm-vs-reads" class="chart"></div></div>
      </div>

      <!-- Section Divider -->
      <div style="border-top: 3px solid var(--line); margin: 3rem 0 2rem; padding-top: 1rem;">
        <h2 style="font-size: 1.75rem; margin-bottom: 0.5rem;">🧬 Functional Region Landscape</h2>
        <p class="muted" style="margin: 0;">Comprehensive analysis of genomic region assignment across the cohort</p>
      </div>
      
      <div class="chart-grid">
        <div class="card">
          <div style="display:flex; gap:10px; align-items:center; margin-bottom:6px;">
            <strong>Display:</strong>
            <label style="display:inline-flex;gap:6px;align-items:center;font-size:0.9em"><input id="toggle-region-percent" type="checkbox"> Percent composition</label>
          </div>
          <div id="chart-reg-reads" class="chart"></div>
        </div>
        <div class="card"><div id="chart-reg-counts" class="chart"></div></div>
        <div class="card"><div id="chart-reg-lentot" class="chart"></div></div>
        <div class="card"><div id="chart-reg-lenmed" class="chart"></div></div>
      </div>

      <h2>Per-Condition Composition</h2>
      <div id="cond-pies" class="chart-grid"></div>

      <h2>Compare Conditions</h2>
      <div class="card">
        <div style="display:flex; gap:8px; align-items:center; margin-bottom:6px;">
          <strong>Compare</strong>
          <select id="cmp-cond-a" style="padding:4px 6px; border:1px solid var(--line); background: var(--card); color: var(--fg);"></select>
          <strong>vs</strong>
          <select id="cmp-cond-b" style="padding:4px 6px; border:1px solid var(--line); background: var(--card); color: var(--fg);"></select>
          <button id="cmp-render" style="margin-left:auto; padding:4px 8px; border:1px solid var(--line); background:var(--card); color:var(--fg); cursor:pointer;">Render</button>
        </div>
        <div id="chart-cond-compare" class="chart"></div>
      </div>

      <!-- Section Divider -->
      <div style="border-top: 3px solid var(--line); margin: 3rem 0 2rem; padding-top: 1rem;">
        <h2 id="samples" style="font-size: 1.75rem; margin-bottom: 0.5rem;">📋 Sample Data Explorer</h2>
        <p class="muted" style="margin: 0;">Detailed metrics for all samples with advanced filtering</p>
      </div>
      
      <div class="card" style="margin: 1rem 0;">
        <h3>🔍 Filter & Search Samples</h3>
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

      <div class="card" style="margin: 2rem 0; background: linear-gradient(135deg, rgba(37,99,235,0.02) 0%, rgba(147,51,234,0.02) 100%);">
        <h3 style="font-size: 1.3rem; margin-bottom: 1rem;">🧬 Functional Region Assignment Details</h3>
        <div style="background: var(--bg); padding: 1rem; border-radius: 8px; margin-bottom: 1rem;">
          <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1rem; font-size: 0.9em;">
            <div>
              <strong style="color: rgb(243,132,0);">● Promoter</strong> — Active transcription start sites (DT sites overlapping gene promoters)
            </div>
            <div>
              <strong style="color: rgb(178,59,212);">● DivergentTx</strong> — Upstream divergent transcription regions (opposite strand from promoters)
            </div>
            <div>
              <strong style="color: rgb(115,212,122);">● Enhancers</strong> — Active enhancers (DT sites NOT near gene promoters)
            </div>
            <div>
              <strong style="color: rgb(0,0,0);">● Gene Body</strong> — Transcribed gene regions between promoter and CPS
            </div>
            <div>
              <strong style="color: rgb(103,200,249);">● CPS</strong> — Cleavage and polyadenylation sites (3' processing)
            </div>
            <div>
              <strong style="color: rgb(253,218,13);">● Short Genes</strong> — Genes ≤750bp (too short for region subdivision)
            </div>
            <div>
              <strong style="color: rgb(255,54,98);">● Termination</strong> — Termination windows downstream of genes
            </div>
            <div>
              <strong style="color: rgb(160,170,180);">● Non-localized</strong> — Reads not assigned to defined functional regions
            </div>
          </div>
        </div>
        
        <div class="tbl-container">
          <table class="tbl">
            <thead>
              <tr>
                <th class="sticky">Sample</th>
                <th>Condition</th>
                <th title="Active promoter sites (DT overlapping genes)">Promoter</th>
                <th title="Transcribed gene bodies">Gene Body</th>
                <th title="Cleavage/polyadenylation sites">CPS</th>
                <th title="Divergent transcription regions">DivergentTx</th>
                <th title="Active enhancers (DT not near genes)">Enhancers</th>
                <th title="Short genes (≤750bp)">Short Genes</th>
                <th title="Transcription termination windows">Termination</th>
                <th title="Reads outside functional regions">Non-localized</th>
              </tr>
            </thead>
            <tbody id="functional-table-body"></tbody>
          </table>
        </div>
      </div>

      <!-- Section Divider -->
      <div style="border-top: 3px solid var(--line); margin: 3rem 0 2rem; padding-top: 1rem;">
        <h2 style="font-size: 1.75rem; margin-bottom: 0.5rem;">📁 Results Navigation & Reproducibility</h2>
        <p class="muted" style="margin: 0;">Access your data and reproduce this analysis</p>
      </div>
      
      <div class="card">
        <h3>🎯 Quick Access to Output Files</h3>
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
