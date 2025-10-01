#!/usr/bin/env python3
# nf-monitor v11 — colorful TUI + log bootstrap + filesystem probe + per-task CPU/RSS
# stdlib-only, Python 3.8+, macOS/Linux
# Self-contained monitoring tool for TrackTx pipeline

import argparse, curses, io, json, os, re, signal, subprocess, sys, time
from dataclasses import dataclass, field
from typing import Dict, Optional, List, Tuple
from collections import deque

# Python version compatibility
try:
    from typing import Pattern
except ImportError:
    # For older Python versions, use Any
    from typing import Any as Pattern

# ─────────────────────────── CLI ───────────────────────────
def parse_args():
    ap = argparse.ArgumentParser(description="nf-monitor v11 — Nextflow TUI + live introspection")
    ap.add_argument("--log", default=".nextflow.log")
    ap.add_argument("--trace", default="results/trace/trace.txt")
    ap.add_argument("--work", default="")
    ap.add_argument("--refresh", type=float, default=0.8)
    ap.add_argument("--tail", type=int, default=20)
    ap.add_argument("--oneshot", action="store_true")
    ap.add_argument("--json", default="")
    ap.add_argument("--filter", default="")
    ap.add_argument("--no-alt", action="store_true")
    ap.add_argument("--simple", action="store_true", help="force simple TUI (no Rich)")
    ap.add_argument("--mini", action="store_true", help="alias for --simple")
    ap.add_argument("--all-logs", action="store_true", help="show live logs for all running tasks (stacked)")
    ap.add_argument("--log-rows", type=int, default=24, help="preferred height (rows) for live log panel in Rich UI")
    ap.add_argument("--resolve-hash", default="", help="print process info for a given task hash and exit")
    ap.add_argument("--from-start", action="store_true", help="ingest full .nextflow.log and trace history (default: tail current run only)")
    ap.add_argument("--retain-sec", type=int, default=0, help="seconds to retain terminal tasks in memory (0 = clear immediately)")
    ap.add_argument("--max-tasks", type=int, default=2000, help="approximate max tasks to keep (older terminal tasks pruned)")
    return ap.parse_args()

# ───────────────────────── Regexes ─────────────────────────
ANSI_RE   = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
RE_PREFIX = re.compile(r"(Submitted|Cached|Completed|Failed|Error executing|Killed|Terminated)\s+process\s*>\s*(?P<name>[^\s(]+)(?:\s*\((?P<tag>[^)]+)\))?(?:\s*\[(?P<id>[0-9a-f]{2}/[0-9a-f]+)\])?", re.I)
RE_SUFFIX = re.compile(r"process\s*>\s*(?P<name>[^\s(]+)(?:\s*\((?P<tag>[^)]+)\))?(?:\s*\[(?P<id>[0-9a-f]{2}/[0-9a-f]+)\])?.*?\b(Submitted|Cached|Completed|succeeded|Failed|FAILED|Error executing|Killed|Terminated|Retrying)\b", re.I)
RE_EXEC   = re.compile(r"(?:^|\s)executor\s*>\s*(?P<exec>[^\s]+)", re.I)
RE_RUN    = re.compile(r"(?:Workflow run name:|(?:^|\s)runName\s*[:=]\s*)(?P<name>[A-Za-z0-9_.:-]+)")
RE_SES    = re.compile(r"(?:^|\s)(?:session:|sessionId[ :=]+|Session id[ :=]+)(?P<sid>[A-Za-z0-9-]+)")
RE_WORK   = re.compile(r"(?:^|\s)(?:workDir|Working (?:dir|directory))\s*:\s*(?P<dir>\S+)")
RE_TASKDIR= re.compile(r".*/work/[0-9a-f]{2}/[0-9a-f]+$")

def norm_state(tok: str) -> str:
    t = tok.lower()
    if "retry" in t: return "RETRYING"
    if "submit" in t: return "SUBMITTED"
    if "cache"  in t: return "CACHED"
    if "succeed" in t or "completed" in t: return "COMPLETED"
    if "error executing" in t or "fail" in t: return "FAILED"
    if "kill" in t or "terminate" in t: return "KILLED"
    return ""

# ───────────────────────── Model ──────────────────────────
@dataclass
class Task:
    id: str                      # "aa/bbbbbbbb…" or provisional "_pre:name|tag"
    name: str = ""
    tag: str = ""
    state: str = ""              # SUBMITTED/RETRYING/RUNNING/COMPLETED/FAILED/CACHED/KILLED
    first_ts: float = field(default_factory=time.time)
    last_ts: float = field(default_factory=time.time)
    retries: int = 0
    workdir: str = ""            # resolved when known
    pid: Optional[int] = None    # from .command.pid if present
    cpu_pct: Optional[float] = None
    rss_mb: Optional[float] = None
    duration_ms: Optional[int] = None  # from trace duration
    stage: Optional[str] = None  # best-effort current stage parsed from user logs
    cpus: Optional[float] = None # requested/allocated cpus from trace
    cpu_hist: deque = field(default_factory=lambda: deque(maxlen=30))
    def label(self) -> str:
        return f"{self.name} ({self.tag})" if self.tag else self.name

@dataclass
class Meta:
    run_name: str = "?"
    session: str  = "?"
    executor: str = "?"
    work_root: str = ""

@dataclass
class World:
    tasks: Dict[str, Task] = field(default_factory=dict)
    meta: Meta = field(default_factory=Meta)
    recent_errors: deque = field(default_factory=lambda: deque(maxlen=6))
    start_ts: float = field(default_factory=time.time)
    last_log_pos: int = 0
    last_log_ino: Optional[int] = None
    trace_path: Optional[str] = None
    last_trace_pos: int = 0
    last_trace_ino: Optional[int] = None
    trace_header: Optional[List[str]] = None
    tail_n: int = 5
    filt: Optional[Pattern] = None
    retain_sec: int = 900
    max_tasks: int = 2000
    # system
    cpu_pct: int = 0
    mem_pct: int = 0
    load_1: str = "0.00"
    ncpu: int = 0
    # cumulative run stats (persist even if tasks are pruned)
    seen_ids: set = field(default_factory=set)
    completed_ids: set = field(default_factory=set)
    cum_seen: int = 0
    cum_done: int = 0
    cum_cached: int = 0
    cum_failed: int = 0
    cum_killed: int = 0

# ─────────────────────── FS & metrics ───────────────────────
def abspath(p: str) -> str: return p if os.path.isabs(p) else os.path.abspath(p)

def guess_work_root(log_path: str, cli: str) -> str:
    if cli: return abspath(cli)
    # from log
    try:
        with open(log_path, "rb") as fh:
            txt = ANSI_RE.sub("", fh.read().decode("utf-8","ignore"))
        m = RE_WORK.search(txt)
        if m: return abspath(m.group("dir"))
    except Exception: pass
    # heuristic: prefer ./work with hashed layout
    for candidate in ("work", "../work", "../../work"):
        c = abspath(candidate)
        if os.path.isdir(c):
            # looks like nf work?
            for root, dirs, files in os.walk(c):
                if RE_TASKDIR.match(root):
                    return c
            # fallback accept
            return c
    return abspath(os.environ.get("NXF_WORK") or "work")

def sys_metrics()->Tuple[int,int,str]:
    cpu=0; mem=0; load="0.00"
    try:
        if sys.platform.startswith("darwin"):
            ncpu=int(subprocess.check_output(["sysctl","-n","hw.ncpu"],text=True))
            vals=[float(x) for x in subprocess.check_output(["ps","-A","-o","%cpu="],text=True).split() if x.strip()]
            cpu=int(round(min(100.0, max(0.0, sum(vals)/max(1,ncpu)))))
        else:
            ncpu=int(subprocess.check_output(["bash","-lc","nproc 2>/dev/null || getconf _NPROCESSORS_ONLN"],text=True))
            vals=[float(x) for x in subprocess.check_output(["bash","-lc","ps -A -o %cpu= || true"],text=True).split() if x.strip()]
            cpu=int(round(min(100.0, max(0.0, sum(vals)/max(1,ncpu)))))
    except Exception: pass
    try:
        if sys.platform.startswith("darwin"):
            total=int(subprocess.check_output(["sysctl","-n","hw.memsize"],text=True))
            vm=subprocess.check_output(["vm_stat"],text=True)
            pages=0
            for ln in vm.splitlines():
                if any(k in ln for k in ("Pages free","Pages inactive","Pages speculative")):
                    pages += int(ln.split(":")[1].strip().strip("."))
            avail=pages*4096; used=max(0,total-avail)
            mem=int(round(used*100/max(1,total)))
        else:
            with open("/proc/meminfo") as fh:
                m={ln.split(":")[0]: int(ln.split(":")[1].strip().split()[0]) for ln in fh}
            total=m.get("MemTotal",0); avail=m.get("MemAvailable",0)
            mem=int(round((total-avail)*100/max(1,total)))
    except Exception: pass
    try:
        if sys.platform.startswith("darwin"):
            up=subprocess.check_output(["uptime"],text=True)
            for key in ("load averages:","load average:"):
                if key in up: load=up.split(key)[-1].split(",")[0].strip()
        else:
            with open("/proc/loadavg") as fh: load=fh.read().split()[0]
    except Exception: pass
    return cpu,mem,load

def du_h(path:str)->str:
    if not os.path.isdir(path): return "--"
    try: return subprocess.check_output(["du","-sh",path], text=True, stderr=subprocess.DEVNULL).split()[0]
    except Exception: return "--"

# ───────────────────── Log I/O + bootstrap ─────────────────────
def _filter_lines(s: str)->List[str]:
    return [l for l in s.splitlines() if not l.lstrip().startswith("DEBUG") and "Missed collect-file cache" not in l]

def read_all_lines(path:str)->List[str]:
    try:
        with open(path,"rb") as fh:
            return _filter_lines(ANSI_RE.sub("", fh.read().decode("utf-8","ignore")))
    except Exception:
        return []

def tail_lines(world:World, path:str)->List[str]:
    try: st=os.stat(path)
    except OSError: return []
    reopen = (world.last_log_ino is None or st.st_ino!=world.last_log_ino or st.st_size<world.last_log_pos)
    if reopen:
        f=open(path,"rb"); world.last_log_ino=st.st_ino; world.last_log_pos=0
    else:
        f=open(path,"rb"); f.seek(world.last_log_pos, io.SEEK_SET)
    chunk=f.read(); world.last_log_pos=f.tell(); f.close()
    if not chunk: return []
    return _filter_lines(ANSI_RE.sub("", chunk.decode("utf-8","ignore")))

# ───────────────────── Trace I/O (CSV incremental) ───────────────────
def _open_trace(world:World, path:str):
    try:
        st=os.stat(path)
    except OSError:
        return None
    f=open(path,"rb"); world.last_trace_ino=st.st_ino; world.last_trace_pos=0; world.trace_header=None
    return f

def _parse_trace_header(line:str)->Optional[List[str]]:
    # split CSV respecting quoted commas
    try:
        import csv
        return next(csv.reader([line]))
    except Exception:
        return None

def _parse_trace_row(header:List[str], line:str)->Optional[Dict[str,str]]:
    try:
        import csv
        vals=next(csv.reader([line]))
        if len(vals)!=len(header):
            return None
        return {header[i]: vals[i] for i in range(len(header))}
    except Exception:
        return None

def tail_trace(world:World)->List[Dict[str,str]]:
    p=world.trace_path
    if not p or not os.path.isfile(p):
        return []
    try:
        st=os.stat(p)
    except OSError:
        return []
    reopen = (world.last_trace_ino is None or st.st_ino!=world.last_trace_ino or st.st_size<world.last_trace_pos)
    if reopen:
        f=_open_trace(world,p)
        if not f: return []
    else:
        f=open(p,"rb"); f.seek(world.last_trace_pos, io.SEEK_SET)
    chunk=f.read(); world.last_trace_pos=f.tell(); f.close()
    if not chunk:
        return []
    text=chunk.decode("utf-8","ignore")
    lines=[ln for ln in text.splitlines() if ln.strip()]
    out=[]
    for ln in lines:
        if world.trace_header is None:
            hdr=_parse_trace_header(ln)
            if hdr: world.trace_header=hdr
            continue
        row=_parse_trace_row(world.trace_header, ln)
        if row: out.append(row)
    return out

def apply_trace_row(world:World, row:Dict[str,str], now:float):
    # common columns: task_id/hash, process, tag, status, workdir, %cpu, rss, start, complete, realtime
    name=(row.get("process") or row.get("name") or "").strip()
    tag=(row.get("tag") or "").strip()
    tid=(row.get("hash") or row.get("task_id") or "").strip()
    workdir=(row.get("workdir") or "").strip()
    status=(row.get("status") or "").strip().upper()
    t=_ensure(world, tid or _prekey(name,tag), name, tag, now, provisional=not bool(tid))
    if tid: t.id=tid
    if workdir: t.workdir=workdir
    if status: t.state=status
    # populate meta from trace when available
    try:
        rn=(row.get("runName") or row.get("run_name") or "").strip()
        if rn: world.meta.run_name = rn
    except Exception:
        pass
    try:
        sid=(row.get("session") or row.get("session_id") or row.get("sessionId") or "").strip()
        if sid: world.meta.session = sid
    except Exception:
        pass
    # resources
    def _parse_float(x:Optional[str])->Optional[float]:
        try:
            return float(x)
        except Exception:
            return None
    cpu=_parse_float(row.get("%cpu") or row.get("pcpu"))
    if cpu is not None: t.cpu_pct=cpu
    # cpus requested/allocated
    cpus_val = row.get("cpus") or row.get("cpu") or row.get("ncpus")
    cpus=_parse_float(cpus_val) if cpus_val is not None else None
    if cpus is not None:
        t.cpus = cpus
    rss=_parse_float(row.get("rss"))
    if rss is not None:
        # Nextflow trace rss is in MB by default
        t.rss_mb=rss
    # durations
    dur=row.get("duration") or ""
    if dur:
        # duration in ms
        try:
            t.duration_ms=int(dur)
        except Exception:
            pass
    elif row.get("realtime"):
        rt=row.get("realtime")
        # parse HH:MM:SS(.ms)
        try:
            parts=rt.split(":")
            if len(parts)==3:
                h=int(parts[0]); m=int(parts[1]); s=float(parts[2])
                t.duration_ms=int((h*3600+m*60+s)*1000)
        except Exception:
            pass

# ───────────────────── Parser / state machine ───────────────────
def _prekey(name:str, tag:str)->str: return f"_pre:{name}|{tag}"

def _ensure(world:World, key:str, name:str, tag:str, now:float, provisional:bool)->Task:
    t=world.tasks.get(key)
    if not t:
        t=Task(id=key, name=name or "", tag=tag or "", first_ts=now, last_ts=now)
        t.state = "SUBMITTED" if provisional else t.state
        world.tasks[key]=t
    else:
        if name and not t.name: t.name=name
        if tag and not t.tag:  t.tag=tag
        t.last_ts=now
    return t

def _apply(world:World, name:str, tag:str, tid:str, token:str, now:float):
    st = norm_state(token)
    if tid:
        t=_ensure(world, tid, name, tag, now, provisional=False)
        # cumulative tracking
        if tid not in world.seen_ids:
            world.seen_ids.add(tid); world.cum_seen += 1
        if st=="RETRYING": t.retries+=1
        if st: t.state=st
        pk=_prekey(name,tag)
        if pk in world.tasks:
            old=world.tasks.pop(pk)
            t.first_ts=min(t.first_ts, old.first_ts)
    else:
        pk=_prekey(name,tag)
        t=_ensure(world, pk, name, tag, now, provisional=True)
        if st: t.state=st

def update_meta(world:World, line:str):
    m=RE_EXEC.search(line); world.meta.executor = (m.group("exec") if m else world.meta.executor)
    m=RE_RUN.search(line);  world.meta.run_name = (m.group("name") if m else world.meta.run_name)
    m=RE_SES.search(line);  world.meta.session  = (m.group("sid")  if m else world.meta.session)
    m=RE_WORK.search(line); 
    if m:
        cand = abspath(m.group("dir"))
        # Normalize: if this looks like a task dir .../work/aa/hash, lift to .../work
        try:
            parts=cand.rstrip("/").split("/")
            if len(parts)>=3 and parts[-3].endswith("work") and re.fullmatch(r"[0-9a-f]{2}", parts[-2] or "") and re.fullmatch(r"[0-9a-f]+", parts[-1] or ""):
                cand = "/".join(parts[:-2])
        except Exception:
            pass
        world.meta.work_root = cand

def parse_line(world:World, line:str, now:float):
    m=RE_PREFIX.search(line)
    if m:
        _apply(world, (m.group("name") or "").strip(), (m.group("tag") or "").strip(), (m.group("id") or "").strip(), m.group(1) or "", now); return
    m=RE_SUFFIX.search(line)
    if m:
        tail=line[m.end()-40:].lower()
        token="Completed" if ("succeed" in tail or "completed" in tail) else tail
        _apply(world, (m.group("name") or "").strip(), (m.group("tag") or "").strip(), (m.group("id") or "").strip(), token, now)

def bootstrap(world:World, log_path:str):
    for ln in read_all_lines(log_path):
        update_meta(world, ln)
        parse_line(world, ln, time.time())
        if any(k in ln for k in ("ERROR","FAILED","Exception","Caused by:")) and "collect-file" not in ln:
            world.recent_errors.append(ln.strip())
    try: st=os.stat(log_path); world.last_log_ino=st.st_ino; world.last_log_pos=st.st_size
    except OSError: pass

# ───────────────────── Filesystem probe (authoritative running) ─────────────
def find_work_roots(seed: str) -> List[str]:
    roots=set()
    # primary
    if os.path.isdir(seed): roots.add(abspath(seed))
    # discover any nearby 'work' dirs with hashed layout
    for base in (os.getcwd(), os.path.dirname(os.getcwd())):
        cand=os.path.join(base,"work")
        if os.path.isdir(cand): roots.add(abspath(cand))
    # scan .nextflow* cache paths if present
    return list(sorted(roots))

def label_from_dir(d: str) -> Tuple[str,str]:
    """name, tag from .command.env/.command.begin/.command.log"""
    def get_env(k: str, f: str)->str:
        try:
            for ln in open(f,"r",encoding="utf-8",errors="ignore"):
                if ln.startswith(k+"="):
                    val=ln.split("=",1)[1].strip().strip('"')
                    return val
        except Exception: pass
        return ""
    env=os.path.join(d,".command.env")
    name=get_env("NXF_PROCESS", env) or get_env("NXF_TASK_NAME", env)
    tag =get_env("NXF_TASK_TAG", env)
    if not name:
        # .command.begin → 'process > NAME (TAG)'
        for src in (".command.begin",".command.log"):
            p=os.path.join(d,src)
            try:
                for ln in open(p,"r",encoding="utf-8",errors="ignore"):
                    if "process >" in ln:
                        nm=re.sub(r'.*process\s*>\s*','',ln).strip()
                        nm=re.sub(r'[\s(].*','',nm)
                        name=nm; m=re.search(r'\(([^)]+)\)', ln)
                        if m: tag=m.group(1)
                        break
            except Exception: pass
    return name or os.path.basename(d), tag or ""

def pid_from_dir(d: str) -> Optional[int]:
    for cand in (".command.pid",".nxf.pid"):  # try both
        p=os.path.join(d,cand)
        try:
            txt=open(p).read().strip()
            if txt.isdigit(): return int(txt)
        except Exception: pass
    return None

def cpus_from_dir(d: str) -> Optional[float]:
    """Best-effort read of allocated CPUs from .command.env (NXF_CPUS/task.cpus)."""
    try:
        env=os.path.join(d, ".command.env")
        val=""
        with open(env, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if ln.startswith("NXF_CPUS=") or ln.startswith("task.cpus="):
                    val=ln.split("=",1)[1].strip().strip('"')
                    break
        if val:
            try:
                return float(val)
            except Exception:
                pass
    except Exception:
        pass
    return None

def ps_for_pid(pid: int) -> Tuple[Optional[float], Optional[float]]:
    """Return (%cpu, rss_mb) for pid or (None,None)."""
    try:
        if sys.platform.startswith("darwin"):
            out=subprocess.check_output(["ps","-p",str(pid),"-o","%cpu=,rss="], text=True)
        else:
            out=subprocess.check_output(["ps","-p",str(pid),"-o","%cpu=,rss="], text=True)
        parts=[p for p in out.strip().split() if p]
        if len(parts)>=2:
            cpu=float(parts[0]); rss_kb=float(parts[1]); return cpu, rss_kb/1024.0
    except Exception: pass
    return None, None

def _ps_snapshot() -> List[Tuple[int, int, float, float, str]]:
    """Return a lightweight process snapshot [(pid, ppid, pcpu, rss_mb, cmdline)]"""
    out=[]
    try:
        if sys.platform.startswith("darwin"):
            txt = subprocess.check_output(["ps","-Ao","pid=,ppid=,pcpu=,rss=,command="], text=True)
        else:
            txt = subprocess.check_output(["bash","-lc","ps -Ao pid=,ppid=,pcpu=,rss=,command= || true"], text=True)
        for ln in txt.splitlines():
            try:
                parts=ln.strip().split(maxsplit=4)
                if len(parts)>=5:
                    pid=int(parts[0]); ppid=int(parts[1]); pcpu=float(parts[2]); rss_kb=float(parts[3]); cmd=parts[4]
                    out.append((pid, ppid, pcpu, rss_kb/1024.0, cmd))
            except Exception:
                continue
    except Exception:
        pass
    return out

def _docker_running_task_stats() -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    """Map nf-task-hash -> (cpu_pct, rss_mb) using docker ps/stats if available.
    Returns empty dict on any error.
    """
    try:
        # Build map: hash -> container ID
        ps_out = subprocess.check_output(
            ["bash","-lc","docker ps --format '{{.ID}} {{.Label ""nf-task-hash""}}'"],
            text=True, stderr=subprocess.DEVNULL
        )
    except Exception:
        return {}
    cid_by_hash = {}
    for ln in ps_out.splitlines():
        parts=[p for p in ln.strip().split() if p]
        if len(parts)>=2:
            cid, thash = parts[0], parts[1]
            if thash and thash != "<no value>":
                cid_by_hash[thash] = cid
    if not cid_by_hash:
        return {}
    # Query docker stats for these containers
    try:
        ids = " ".join(cid_by_hash.values())
        st_out = subprocess.check_output(
            ["bash","-lc", f"docker stats --no-stream --format '{{{{.Container}}}} {{{{.CPUPerc}}}} {{{{.MemUsage}}}}' {ids}"],
            text=True, stderr=subprocess.DEVNULL
        )
    except Exception:
        return {}
    # Build inverse map: cid -> hash
    hash_by_cid = {cid: h for h, cid in cid_by_hash.items()}
    stats = {}
    for ln in st_out.splitlines():
        parts=ln.strip().split()
        if not parts:
            continue
        cid = parts[0]
        cpu_str = parts[1] if len(parts)>1 else ""
        mem_str = parts[2] if len(parts)>2 else ""
        # Parse cpu like '12.34%'
        try:
            cpu = float(cpu_str.strip('%'))
        except Exception:
            cpu = None
        # Parse mem like '123.4MiB/..' or '1.2GiB/..'
        rss_mb = None
        try:
            m = mem_str.split('/')[0]
            if m.lower().endswith('gib') or m.lower().endswith('gb'):
                rss_mb = float(m[:-3]) * 1024.0
            elif m.lower().endswith('mib') or m.lower().endswith('mb'):
                rss_mb = float(m[:-3])
            elif m.lower().endswith('kib') or m.lower().endswith('kb'):
                rss_mb = float(m[:-3]) / 1024.0
        except Exception:
            pass
        h = hash_by_cid.get(cid)
        if h:
            stats[h] = (cpu, rss_mb)
    return stats

def running_dirs(work_roots: List[str]) -> List[str]:
    """Active task dirs: have .command.run/.command.sh and NO .exitcode."""
    out=[]
    for root in work_roots:
        # depth=2 layout: work/aa/abcdef... → scan two levels cheaply
        try:
            for aa in os.listdir(root):
                if len(aa)!=2: continue
                p1=os.path.join(root,aa)
                if not os.path.isdir(p1): continue
                for bb in os.listdir(p1):
                    d=os.path.join(p1,bb)
                    if not os.path.isdir(d): continue
                    if os.path.exists(os.path.join(d,".exitcode")): continue
                    if os.path.exists(os.path.join(d,".command.run")) or os.path.exists(os.path.join(d,".command.sh")):
                        out.append(d)
        except Exception: pass
    # sort by mtime (recent first)
    out.sort(key=lambda x: os.path.getmtime(x) if os.path.exists(x) else 0, reverse=True)
    return out

# ───────────────────── Introspection (payload/tail/outputs) ─────────────
WRAP_DROP = [
    re.compile(r"^\+{1,}"), re.compile(r"^set -"), re.compile(r"^trap "),
    re.compile(r"^ulimit "), re.compile(r"^exec > "), re.compile(r"^nxf_"),
    re.compile(r"^NXF_"), re.compile(r"^CAPSULE:"), re.compile(r"^Picked up _JAVA_OPTIONS"),
    re.compile(r"^Warning: .*illegal reflective access"), re.compile(r"^INFO +\(.*Nextflow.*\)"),
    re.compile(r"^ *\["), re.compile(r"read -t \d+ -r DONE")
]
def keep_line(s: str)->bool:
    s=ANSI_RE.sub("", s.rstrip("\r"))
    for rx in WRAP_DROP:
        if rx.search(s): return False
    if s.lstrip().startswith("echo "): return False
    return bool(s.strip())

def slurp_tail(path: str, n: int, scrub=True)->List[str]:
    try:
        with open(path,"rb") as fh:
            fh.seek(0, io.SEEK_END); size=fh.tell()
            block=4096; data=b""
            while len(data.splitlines())<=n+1 and fh.tell()>0:
                off=max(0, fh.tell()-block); fh.seek(off)
                data=fh.read(size-off)+data; fh.seek(off)
                if off==0: break; block*=2
        lines=[ANSI_RE.sub("", l) for l in data.decode("utf-8","ignore").splitlines()[-n:]]
        return [l for l in lines if keep_line(l)] if scrub else lines
    except Exception:
        return []

def payload_from(workdir:str)->str:
    for fn in (".command.sh",".command.run"):
        p=os.path.join(workdir,fn)
        try:
            last=""
            for ln in open(p,"r",encoding="utf-8",errors="ignore"):
                if keep_line(ln): last=ln.strip()
            if last: return last
        except Exception: pass
    return "<no payload detected>"

def best_log_for_tail(d:str)->Optional[str]:
    cand=[os.path.join(d, ".command.out"), os.path.join(d, ".command.err"), os.path.join(d, ".command.log")]
    # choose newest
    cand=[p for p in cand if os.path.exists(p)]
    if not cand: return None
    cand.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return cand[0]

def newest_user_log(d: str)->Optional[str]:
    try:
        newest=None; mt=-1
        for nm in os.listdir(d):
            if nm.startswith(".command.") or nm==".exitcode": continue
            if any(s in nm.lower() for s in ("stdout","stderr",".log",".err",".out")):
                p=os.path.join(d,nm); m=os.path.getmtime(p)
                if m>mt: newest, mt = p, m
        return newest
    except Exception: return None

def list_outputs(d: str, k=3)->List[Tuple[str,str]]:
    out=[]
    try:
        for nm in os.listdir(d):
            if nm.startswith("."): continue
            p=os.path.join(d,nm)
            if os.path.isfile(p):
                st=os.stat(p); out.append((nm, st.st_size, st.st_mtime))
        out.sort(key=lambda x:(-x[2], -x[1]))
        def human(n):
            for u in ("B","KB","MB","GB","TB","PB"):
                if n<1024: return f"{n:.1f}{u}"; n/=1024
            return f"{n:.1f}EB"
        return [(n, human(sz)) for n,sz,_ in out[:k]]
    except Exception: return []

# ───────────────────── Stage detection (best-effort) ───────────────
STEP_RE = re.compile(r"\bStep\s+\d+\s*/\s*\d+\s*:\s*(.+)")
def detect_stage_from_tail(lines: List[str]) -> Optional[str]:
    for ln in reversed(lines or []):
        m=STEP_RE.search(ln)
        if m:
            return m.group(0).strip()
    # also support INFO  [tracks] … style
    for ln in reversed(lines or []):
        if "INFO  Step" in ln:
            return ln.strip()
    return None

# ───────────────────── Classify (merge log + FS) ───────────────
def classify(world:World):
    # 1) start with log/trace-based states
    done=fail=cache=killed=0; run=[]; queued=[]
    for t in world.tasks.values():
        st=t.state or ""
        if st in ("COMPLETED","SUCCEEDED"): done+=1
        elif st=="FAILED":  fail+=1
        elif st=="CACHED":  cache+=1
        elif st=="KILLED":  killed+=1
        else: queued.append(t)

    # 2) overlay filesystem: authoritative running set
    roots=find_work_roots(world.meta.work_root or "work")
    active_dirs=running_dirs(roots)
    seen_ids=set()
    docker_stats = _docker_running_task_stats()
    for d in active_dirs:
        # id = aa/bbbbb…
        tid="/".join(d.rstrip("/").split("/")[-2:])
        seen_ids.add(tid)
        t = world.tasks.get(tid)
        if not t:
            # construct from dir (even if log never emitted id)
            name, tag = label_from_dir(d)
            t=Task(id=tid, name=name, tag=tag, state="RUNNING", first_ts=os.path.getmtime(d), last_ts=time.time())
            world.tasks[tid]=t
        t.workdir=d; t.state="RUNNING"
        # enrich with PID + ps
        pid=pid_from_dir(d)
        t.pid=pid
        if pid:
            cpu,rss=ps_for_pid(pid); t.cpu_pct=cpu; t.rss_mb=rss
        # learn cpus from env if missing
        if t.cpus is None:
            cc=cpus_from_dir(d)
            if cc is not None:
                t.cpus = cc
        # record cpu history (prefer docker ps override if present applied below)
        if t.cpu_pct is not None:
            try:
                t.cpu_hist.append(max(0.0, float(t.cpu_pct)))
            except Exception:
                pass
    # Fallback PID/CPU/RSS mapping: try to match processes to work dirs or hashes
    try:
        snap=_ps_snapshot()
        # build children index by ppid
        children={}  # ppid -> list[(pid, pcpu, rss_mb, cmd)]
        for pid, ppid, pcpu, rss_mb, cmd in snap:
            children.setdefault(ppid, []).append((pid, pcpu, rss_mb, cmd))
        for t in [tt for tt in world.tasks.values() if tt.state=="RUNNING" and (tt.pid is None or tt.cpu_pct is None or tt.rss_mb is None) ]:
            thash=(t.id or '').split('/')[-1]
            # match by workdir path fragment or hash in cmdline
            for pid, ppid, pcpu, rss_mb, cmd in snap:
                try:
                    hit=False
                    if thash and thash in cmd:
                        hit=True
                    elif t.workdir and t.workdir in cmd:
                        hit=True
                    if hit:
                        # if this looks like a wrapper (bash, nextflow, java), try to pick hottest child
                        base = cmd.split()[0] if cmd else ""
                        if any(k in base for k in ("/bash","bash","/sh","sh","nextflow","java")) and children.get(pid):
                            kids = children.get(pid) or []
                            kids.sort(key=lambda x: x[1], reverse=True)
                            kpid, kcpu, krss, kcmd = kids[0]
                            pid, pcpu, rss_mb, cmd = kpid, kcpu, krss, kcmd
                        t.pid = t.pid or pid
                        if t.cpu_pct is None: t.cpu_pct = pcpu
                        if t.rss_mb  is None: t.rss_mb  = rss_mb
                        if t.cpu_pct is not None:
                            try: t.cpu_hist.append(max(0.0,float(t.cpu_pct)))
                            except Exception: pass
                        break
                except Exception:
                    continue
    except Exception:
        pass
        # parse stage from latest user logs
        payload, tl, _ = insight(d, world.tail_n)
        st = detect_stage_from_tail(tl)
        if st: t.stage = st
        # docker-based enrichment by nf-task-hash label if present
        # task hash is last component without the aa/
        task_hash = tid.split('/')[-1]
        if task_hash in docker_stats:
            dcpu, drss = docker_stats[task_hash]
            if dcpu is not None: t.cpu_pct = dcpu
            if drss is not None: t.rss_mb = drss
            if t.cpu_pct is not None:
                try:
                    t.cpu_hist.append(max(0.0, float(t.cpu_pct)))
                except Exception:
                    pass

    # Any tasks previously marked RUNNING but not seen in active_dirs this tick
    # are no longer running. To avoid the running list filling up, drop them.
    try:
        for key, t in list(world.tasks.items()):
            if t.state == "RUNNING" and t.id and "/" in t.id:
                if t.id not in seen_ids:
                    # remove stale running entry; will be re-added if it runs again
                    del world.tasks[key]
    except Exception:
        pass

    # re-split after overlay
    now=time.time()
    run=[t for t in world.tasks.values() if t.state=="RUNNING"]
    queued=[t for t in world.tasks.values() if t.state not in ("RUNNING","COMPLETED","FAILED","CACHED","KILLED")]
    # prune terminal tasks per retention policy or immediately if retain_sec==0
    terminal=[(k,t) for k,t in world.tasks.items() if t.state in ("COMPLETED","FAILED","CACHED","KILLED")]
    # immediate prune mode
    if world.retain_sec == 0:
        for k,_ in terminal:
            try: del world.tasks[k]
            except Exception: pass
        terminal=[]
    else:
        # time-based prune first
        for k,t in terminal:
            if world.retain_sec>0 and (now - t.last_ts) > world.retain_sec:
                try: del world.tasks[k]
                except Exception: pass
    # cap-based prune (oldest terminal removed first)
    if len(world.tasks) > world.max_tasks:
        term=[(k,t) for k,t in world.tasks.items() if t.state in ("COMPLETED","FAILED","CACHED","KILLED")]
        term.sort(key=lambda kv: kv[1].last_ts)
        for k,_ in term:
            if len(world.tasks) <= world.max_tasks: break
            try: del world.tasks[k]
            except Exception: pass
    # best-effort name resolution for items lacking names
    for t in run + queued:
        if not (t.name or t.tag):
            ensure_task_label_from_fs(world, t)

    # filter
    if world.filt:
        run=[t for t in run if world.filt.search(t.label() or "")]
        queued=[t for t in queued if world.filt.search(t.label() or "")]

    run.sort(key=lambda x: x.first_ts); queued.sort(key=lambda x: x.first_ts)
    # update cumulative finished counters
    world.cum_done += done
    world.cum_cached += cache
    world.cum_failed += fail
    world.cum_killed += killed
    return run, queued, dict(done=done, failed=fail, cached=cache, killed=killed, running=len(run), queued=len(queued), cum_seen=world.cum_seen, cum_done=world.cum_done, cum_cached=world.cum_cached, cum_failed=world.cum_failed, cum_killed=world.cum_killed)

# ───────────────────── Introspection bundle ───────────────
def insight(d: str, tail_n: int)->Tuple[str, List[str], List[Tuple[str,str]]]:
    payload=payload_from(d)
    user=newest_user_log(d)
    src=user or best_log_for_tail(d)
    tl=slurp_tail(src, tail_n, scrub=True) if src else ["<no output yet>"]
    outs=list_outputs(d)
    return payload, tl, outs

# ───────────────────── UI (colorful) ───────────────
def sec2hms(t:int)->str: h=t//3600; m=(t%3600)//60; s=t%60; return f"{h:02d}:{m:02d}:{s:02d}"
def truncate(s:str,w:int)->str: return s if len(s)<=w else (s[:max(0,w-3)]+"...")
def looks_like_hash(s:str)->bool:
    return bool(re.fullmatch(r"[0-9a-f]{24,64}", s or ""))

def short_hash(h: Optional[str])->str:
    try:
        if not h: return "?"
        return h.replace("/","")[:12]
    except Exception:
        return "?"

def ensure_task_label_from_fs(w: World, t: Task):
    """Populate t.name/tag and workdir from filesystem if missing.
    If workdir is unknown but id looks like aa/bbbbb..., construct from work_root.
    """
    try:
        # determine candidate dir
        d = t.workdir
        if (not d) and t.id and "/" in t.id and (w.meta.work_root or ""):
            d = os.path.join(w.meta.work_root, t.id)
        if d and os.path.isdir(d):
            if not (t.name and t.tag):
                name, tag = label_from_dir(d)
                if name and not t.name: t.name = name
                if tag  and not t.tag:  t.tag  = tag
            if not t.workdir:
                t.workdir = d
    except Exception:
        pass

class TUI:
    def __init__(self, stdscr, w:World, a):
        self.s=stdscr; self.w=w; self.a=a; self.sel=0
        curses.start_color(); curses.use_default_colors()
        curses.init_pair(1, curses.COLOR_CYAN,  -1)  # title
        curses.init_pair(2, curses.COLOR_GREEN, -1)  # ok
        curses.init_pair(3, curses.COLOR_RED,   -1)  # fail
        curses.init_pair(4, curses.COLOR_YELLOW,-1)  # warn
        curses.init_pair(5, curses.COLOR_BLUE,  -1)  # section

    def draw(self):
        self.s.erase(); maxy,maxx=self.s.getmaxyx(); now=int(time.time()-self.w.start_ts)
        run,que,C=classify(self.w)
        # compute cores-in-use (sum cpus of running tasks where known)
        cores_in_use = 0.0
        for t in run:
            if t.cpus is not None:
                cores_in_use += max(0.0, t.cpus)
        # best-effort total cores on system
        try:
            if self.w.ncpu <= 0:
                if sys.platform.startswith("darwin"):
                    self.w.ncpu=int(subprocess.check_output(["sysctl","-n","hw.ncpu"],text=True))
                else:
                    self.w.ncpu=int(subprocess.check_output(["bash","-lc","nproc 2>/dev/null || getconf _NPROCESSORS_ONLN"],text=True))
        except Exception:
            pass
        total_cores = self.w.ncpu or 0
        # rough ETA: median duration per process * remaining non-terminal (run+queued) across processes
        # reuse proc_stats for medians
        # process summary (per process name)
        proc_stats={}
        for t in self.w.tasks.values():
            # prefer real process names; avoid bare hashes
            nm=t.name or "?"
            if looks_like_hash(nm) and t.tag:
                nm = f"{nm[:6]}… ({t.tag})"  # fallback readable
            ps=proc_stats.setdefault(nm, dict(done=0,fail=0,run=0,queued=0,durs=[]))
            st=t.state or ""
            if st in ("COMPLETED","SUCCEEDED"): ps["done"]+=1; 
            elif st=="FAILED": ps["fail"]+=1
            elif st=="RUNNING": ps["run"]+=1
            else: ps["queued"]+=1
            if t.duration_ms: ps["durs"].append(t.duration_ms)
        total=C["running"]+C["queued"]+C["done"]+C["failed"]+C["cached"]+C["killed"]
        pct=int(((C["done"]+C["cached"])*100/total)) if total else 0
        # header
        self._add(0,0, f" nf-monitor ", curses.color_pair(1)|curses.A_BOLD)
        self._add(0,12, f"{time.strftime('%H:%M:%S')}  (run:{self.w.meta.run_name}  session:{self.w.meta.session}  exec:{self.w.meta.executor})")
        self._add(1,0,  f" Root: {os.getcwd()}")
        # fallback cores-in-use using sum of per-task cpu_pct when cpus unknown
        try:
            if cores_in_use <= 0 and total_cores>0:
                approx = 0.0
                for t in run:
                    if t.cpu_pct is not None:
                        approx += max(0.0, t.cpu_pct/100.0)
                cores_in_use = min(float(total_cores), approx)
        except Exception:
            pass
        cores_disp = int(round(cores_in_use)) if cores_in_use>0 else 0
        cores_str = (f" | CORES:{cores_disp}/{total_cores}" if total_cores>0 else "")
        self._add(2,0,  f" Uptime {sec2hms(now)} | CPU:{self.w.cpu_pct}% | MEM:{self.w.mem_pct}% | LOAD1:{self.w.load_1}{cores_str}")
        # bar
        self._add(4,0,"┌─ Pipeline ─┐", curses.color_pair(5))
        # use cumulative progress when available
        total_so_far = C.get('cum_done',0)+C.get('cum_cached',0)+C.get('cum_failed',0)+C.get('cum_killed',0)+C['running']+C['queued']
        pct_cum = int(((C.get('cum_done',0)+C.get('cum_cached',0))*100/max(1,total_so_far))) if total_so_far else 0
        barw=max(10,maxx-12); fill=int(barw*pct_cum/100)
        self._add(5,0,"["+("#"*fill)+"."*(barw-fill)+f"] {pct_cum:3d}%")
        self._add(6,0,f"  🚀 total:{total}  ⚙ run:{C['running']}  ",0)
        self._add(6,28,f"✓ done:{C['done']}  ", curses.color_pair(2))
        self._add(6,42,f"✗ fail:{C['failed']}  ", curses.color_pair(3))
        self._add(6,58,f"cache:{C['cached']}", curses.color_pair(4))
        # start running list below header
        row=8
        # running list
        # continue below whatever row currently is
        self._add(row,0,"┌─ Running tasks ─┐", curses.color_pair(5)); row+=1
        if run:
            # show between 12 and 18 running tasks when space allows, reducing reserved space
            avail=min(18, max(12, maxy-row-6))
            start=max(0, min(self.sel, max(0, len(run)-avail)))
            end=min(len(run), start+avail)
            for i in range(start,end):
                t=run[i]
                age=sec2hms(int(time.time()-t.first_ts))
                # CPU string + ASCII bar
                cpu=(f"{t.cpu_pct:.0f}%" if t.cpu_pct is not None else "--")
                barw = 10
                val = int(min(100, max(0, int(t.cpu_pct))) if t.cpu_pct is not None else 0)
                filled = int(barw * val / 100)
                cpu_bar = ("[" + ("#"*filled) + ("."*(barw-filled)) + "]") if barw>0 else ""
                rss=(f"{t.rss_mb:.0f}MB" if t.rss_mb is not None else "--")
                mark="▶" if i==self.sel else " "
                label = t.label() or (t.id[:10] if t.id else "?")
                if looks_like_hash(label) and t.tag:
                    label = f"{t.name or 'process'} ({t.tag})"
                ih=short_hash(t.id or "")
                stage = (t.stage or "").replace("INFO  ","")
                meta = (f"  • {stage}" if stage else "")
                self._add(row,0, f" {mark} {truncate(label, maxx-70):<{maxx-70}}  id:{ih:<12}  PID:{str(t.pid or '-'):>6}  CPU:{cpu_bar} {cpu:>4}  RSS:{rss:>6}  {age:>8}{meta}")
                row+=1
            self._add(row,0,"└─┘"); row+=1
            # detail pane(s)
            if getattr(self.a, 'all_logs', False):
                self._add(row,0, "┌─ Live logs (all running) ─┐", curses.color_pair(5)); row+=1
                # stack logs for all currently visible running tasks until screen filled
                for i in range(start, end):
                    t=run[i]
                    payload, tl, outs = insight(t.workdir, self.w.tail_n)
                    hdr = t.label() or t.id
                    self._add(row,2, truncate(f"▶ {hdr} (id:{short_hash(t.id or '')}) • {payload}", maxx-4)); row+=1
                    for ln in tl:
                        if row>=maxy-4: break
                        self._add(row,2, truncate(ln, maxx-4)); row+=1
                    if outs and row<maxy-4:
                        self._add(row,2, "outputs: "+", ".join([f"{n} ({sz})" for n,sz in outs])); row+=1
                    if row>=maxy-3: break
                    # separator between tasks
                    if i<end-1 and row<maxy-3:
                        self._add(row,2, "─"*(maxx-4), curses.color_pair(5)); row+=1
                self._add(row,0, "└─┘"); row+=1
            else:
                t=run[self.sel if self.sel< len(run) else len(run)-1]
                payload, tl, outs = insight(t.workdir, self.w.tail_n)
                self._add(row,0, "┌─ Details ─┐", curses.color_pair(5)); row+=1
                hdr = t.label() or t.id
                self._add(row,2, truncate(f"▶ {hdr} (id:{short_hash(t.id or '')}) • {payload}", maxx-4)); row+=1
                for ln in tl: self._add(row,2, truncate(ln, maxx-4)); row+=1
                if outs:
                    self._add(row,2, "outputs: "+", ".join([f"{n} ({sz})" for n,sz in outs])); row+=1
                self._add(row,0, "└─┘"); row+=1
        else:
            self._add(row,2,"none"); row+=1; self._add(row,0,"└─┘"); row+=1
        # (Queued section removed)
        # errors
        self._add(row,0,"┌─ Recent errors ─┐", curses.color_pair(5)); row+=1
        if self.w.recent_errors:
            for e in list(self.w.recent_errors)[-5:]: self._add(row,2, "✗ "+truncate(e,maxx-4), curses.color_pair(3)); row+=1
        else:
            self._add(row,2,"✓ none", curses.color_pair(2)); row+=1
        self._add(row,0,"└─┘"); self.s.refresh()
        # footer with key hints (curses mode)
        try:
            hint = "j/k move  q quit"
            maxy,maxx=self.s.getmaxyx()
            self.s.addstr(maxy-1, 0, hint[:max(0,maxx-1)], curses.A_DIM)
        except curses.error:
            pass

    def _add(self,y,x,s,attr=0):
        try: self.s.addstr(y,x, s[:max(0,self.s.getmaxyx()[1]-1)], attr)
        except curses.error: pass

    def keys(self):
        self.s.nodelay(True); ch=self.s.getch()
        if ch in (ord('q'),27): raise SystemExit
        if ch in (curses.KEY_DOWN, ord('j')): self.sel += 1
        if ch in (curses.KEY_UP,   ord('k')): self.sel = max(0, self.sel-1)
        if ch in (ord('a'),):
            # toggle all-logs mode at runtime
            try:
                setattr(self.a, 'all_logs', not getattr(self.a, 'all_logs', False))
            except Exception:
                pass
        if ch in (ord('+'), ord('=')):
            # no separate log_rows in curses; increase tail lines instead
            try:
                self.w.tail_n = max(1, int(self.w.tail_n) + 2)
            except Exception:
                pass
        if ch == ord('-'):
            try:
                self.w.tail_n = max(1, int(self.w.tail_n) - 2)
            except Exception:
                pass
        if ch == ord(']'):
            try:
                self.w.tail_n = max(1, int(self.w.tail_n) + 2)
            except Exception:
                pass
        if ch == ord('['):
            try:
                self.w.tail_n = max(1, int(self.w.tail_n) - 2)
            except Exception:
                pass

# ───────────────────── JSON / oneshot ───────────────
def build_json(w:World)->dict:
    run,que,C=classify(w)
    return {
        "run_name": w.meta.run_name, "session": w.meta.session, "timestamp": int(time.time()),
        "executor": w.meta.executor, "work_dir": w.meta.work_root,
        "cpu_pct": w.cpu_pct, "mem_pct": w.mem_pct, "load_1m": w.load_1,
        "done": C["done"], "failed": C["failed"], "cached": C["cached"], "killed": C["killed"],
        "running": C["running"], "queued": C["queued"],
        "tasks_running": [dict(id=t.id, label=t.label(), epoch=int(t.first_ts), workdir=t.workdir, pid=t.pid, cpu_pct=t.cpu_pct, rss_mb=t.rss_mb) for t in run],
        "tasks_queued":  [dict(id=t.id, label=t.label(), epoch=int(t.first_ts), state=t.state or "SUBMITTED") for t in que],
        "recent_errors": list(w.recent_errors),
    }

def oneshot(w:World):
    run,que,C=classify(w)
    total=sum(C.values()); pct=int(((C["done"]+C["cached"])*100/total)) if total else 0
    print(f"nf-monitor  {time.strftime('%H:%M:%S')}  (run:{w.meta.run_name}  session:{w.meta.session}  exec:{w.meta.executor})")
    print(f" Root: {os.getcwd()}"); print(f" Uptime {sec2hms(int(time.time()-w.start_ts))} | CPU:{w.cpu_pct}% | MEM:{w.mem_pct}% | LOAD:{w.load_1} | WORK:{du_h(w.meta.work_root)}\n")
    bw=80; fill=int(bw*pct/100); print("┌─ Pipeline ─┐"); print("["+("#"*fill)+"."*(bw-fill)+f"] {pct:3d}%")
    print(f"  🚀 total:{total}  ⚙ run:{C['running']}  💤 queued:{C['queued']}  ✓ done:{C['done']}  ✗ fail:{C['failed']}  cache:{C['cached']}\n")
    print("┌─ Running ─┐")
    for t in run[:20]:
        print(f"  • {t.label():<60} PID:{str(t.pid or '-'):>6} CPU:{(str(int(t.cpu_pct))+'%') if t.cpu_pct is not None else '--':>4} RSS:{(str(int(t.rss_mb))+'MB') if t.rss_mb is not None else '--':>6}")
    if not run: print("  none")
    print("└─┘")
    if que:
        print("\n┌─ Queued ─┐")
        for t in que[:20]: print(f"  • {t.label()}")
        print("└─┘")
    print("\n┌─ Recent errors ─┐")
    if w.recent_errors:
        for e in list(w.recent_errors)[-5:]: print("  ✗", e)
    else: print("  ✓ none")
    print("└─┘")

# ───────────────────── Main ───────────────
def main():
    a=parse_args()
    w=World(); w.tail_n=max(1,a.tail); w.retain_sec = max(0, getattr(a,'retain_sec',0)); w.max_tasks = max(100, getattr(a,'max_tasks',2000))
    if a.filter: w.filt=re.compile(a.filter)
    log=abspath(a.log) if a.log else ""
    w.trace_path=abspath(a.trace) if a.trace else None
    if os.path.isfile(log):
        w.meta.work_root=guess_work_root(log, a.work)
        bootstrap(w, log)
    else:
        # Prefer trace as primary; do not require the Nextflow log
        if a.work:
            w.meta.work_root=abspath(a.work)
        else:
            # best-effort default; may be updated later when trace rows provide workdir
            w.meta.work_root=guess_work_root(".nextflow.log", a.work)
        print(f"Warning: log not found, running in trace-first mode (log={a.log})", file=sys.stderr)
    w.cpu_pct, w.mem_pct, w.load_1 = sys_metrics()

    if a.resolve_hash:
        # best-effort resolve using trace + filesystem
        # normalize: allow full "aa/bbbbb" or just hash
        key=a.resolve_hash.strip()
        wanted=key.split("/")[-1]
        # ingest initial streams so tasks map is populated
        if os.path.isfile(log):
            bootstrap(w, log)
        for row in tail_trace(w):
            apply_trace_row(w, row, time.time())
        # search known tasks
        hit=None
        for t in w.tasks.values():
            if t.id and t.id.split("/")[-1].startswith(wanted):
                hit=t; break
        if not hit:
            # scan filesystem running dirs
            for d in running_dirs(find_work_roots(w.meta.work_root or "work")):
                tid="/".join(d.rstrip("/").split("/")[-2:])
                if tid.split("/")[-1].startswith(wanted):
                    name, tag = label_from_dir(d)
                    hit=Task(id=tid, name=name, tag=tag, state="RUNNING", workdir=d)
                    break
        if hit:
            out=dict(
                id=hit.id,
                short=short_hash(hit.id),
                process=hit.name or "?",
                tag=hit.tag or "",
                workdir=hit.workdir or "",
                state=hit.state or "",
            )
            print(json.dumps(out, indent=2))
        else:
            print(json.dumps({"error":"not found","query":key}, indent=2))
        return

    if a.oneshot:
        if a.json:
            tmp=a.json+".tmp"; open(tmp,"w").write(json.dumps(build_json(w),indent=2)); os.replace(tmp,a.json)
        oneshot(w); return

    def _sig(_s,_f): raise SystemExit
    signal.signal(signal.SIGINT,_sig)

    def loop(stdscr):
        stdscr.nodelay(True); stdscr.keypad(True)
        ui=TUI(stdscr,w,a); next_tick=0.0
        while True:
            now=time.time()
            if now>=next_tick:
                # nextflow log stream (meta + provisional states + errors) — optional
                if log and os.path.isfile(log):
                    for ln in tail_lines(w, log):
                        update_meta(w, ln); parse_line(w, ln, time.time())
                        if any(k in ln for k in ("ERROR","FAILED","Exception","Caused by:")) and "collect-file" not in ln:
                            w.recent_errors.append(ln.strip())
                # trace stream (authoritative states + resources)
                for row in tail_trace(w):
                    apply_trace_row(w, row, time.time())
                    # opportunistically learn work_root from first workdir
                    if not w.meta.work_root and row.get("workdir"):
                        try:
                            wd=row.get("workdir").strip()
                            if wd:
                                # set root two levels up from task dir
                                parts=wd.rstrip("/").split("/")
                                if len(parts)>=3:
                                    w.meta.work_root="/".join(parts[:-2])
                        except Exception:
                            pass
                w.cpu_pct, w.mem_pct, w.load_1 = sys_metrics()
                if a.json:
                    tmp=a.json+".tmp"; open(tmp,"w").write(json.dumps(build_json(w),indent=2)); os.replace(tmp,a.json)
                next_tick = now + max(0.3, a.refresh)
            ui.draw(); ui.keys(); time.sleep(0.05)

    # Default to Rich UI when available unless --simple/--mini
    if not (a.simple or a.mini):
        try:
            from rich.live import Live
            from rich.table import Table
            from rich.panel import Panel
            from rich.progress import Progress, BarColumn, TextColumn
            from rich.layout import Layout
            from rich.align import Align
            from rich.text import Text
            import hashlib

            def _task_label(t: Task) -> str:
                lab = t.label() or (t.id or "?")
                if looks_like_hash(lab):
                    base = t.name or lab[:10]
                    return f"{base} ({t.tag})" if t.tag else base
                return lab

            def _name_style(name: str) -> str:
                # deterministic color per process name
                palette = [
                    "bright_cyan","bright_green","bright_yellow","bright_magenta",
                    "bright_blue","bright_red","cyan","green","yellow","magenta","blue","red"
                ]
                if not name:
                    return "white"
                h = int(hashlib.md5(name.encode("utf-8")).hexdigest(), 16)
                return palette[h % len(palette)]

            # remember focused task across renders and support j/k to switch
            focused = {"id": None, "idx": 0}
            stop_keys = {"stop": False}
            filter_state = {"proc": None}
            def _key_thread():
                try:
                    import sys, termios, tty, select
                    fd = sys.stdin.fileno()
                    old = termios.tcgetattr(fd)
                    tty.setcbreak(fd)
                    while not stop_keys["stop"]:
                        r,_,_ = select.select([sys.stdin], [], [], 0.1)
                        if r:
                            ch = sys.stdin.read(1)
                            if ch in ('q','\x03','\x1b'):
                                stop_keys["stop"] = True
                                raise SystemExit
                            if ch == 'j':
                                focused["idx"] += 1
                            if ch == 'k':
                                focused["idx"] = max(0, focused["idx"]-1)
                            if ch == 'f':
                                # cycle process filter among seen process names (plus All)
                                names = sorted({ (t.name or "?") for t in w.tasks.values() })
                                if filter_state["proc"] is None:
                                    filter_state["proc"] = names[0] if names else None
                                else:
                                    try:
                                        i = names.index(filter_state["proc"]) if filter_state["proc"] in names else -1
                                    except ValueError:
                                        i = -1
                                    nxt = (i+1) % (len(names)+1)
                                    filter_state["proc"] = (None if nxt==len(names) else names[nxt])
                            if ch == 'a':
                                try:
                                    setattr(a, 'all_logs', not getattr(a, 'all_logs', False))
                                except Exception:
                                    pass
                            if ch == 's':
                                # toggle sort: default -> cpu -> rss -> age
                                curr = filter_state.get('sort') or 'default'
                                order = ['default','cpu','rss','age']
                                try:
                                    idx = order.index(curr)
                                except ValueError:
                                    idx = 0
                                filter_state['sort'] = order[(idx+1) % len(order)]
                            if ch in ('+', '='):
                                try:
                                    curr = int(getattr(a, 'log_rows', 24))
                                    setattr(a, 'log_rows', min(60, curr + 2))
                                except Exception:
                                    pass
                            if ch == '-':
                                try:
                                    curr = int(getattr(a, 'log_rows', 24))
                                    setattr(a, 'log_rows', max(8, curr - 2))
                                except Exception:
                                    pass
                            if ch == ']':
                                try:
                                    w.tail_n = max(1, int(w.tail_n) + 2)
                                except Exception:
                                    pass
                            if ch == '[':
                                try:
                                    w.tail_n = max(1, int(w.tail_n) - 2)
                                except Exception:
                                    pass
                    termios.tcsetattr(fd, termios.TCSADRAIN, old)
                except Exception:
                    pass

            def render():
                run,que,C=classify(w)
                total=C["running"]+C["queued"]+C["done"]+C["failed"]+C["cached"]+C["killed"]
                pct=int(((C["done"]+C["cached"])*100/max(1,total))) if total else 0

                # header with run/session and system stats
                hdr = Text()
                # hide placeholders for run/session when unknown
                rn = (w.meta.run_name or "?")
                ss = (w.meta.session or "?")
                hdr.append(f"Run: {rn if rn!='?' else 'n/a'}  ")
                hdr.append(f"Session: {ss if ss!='?' else 'n/a'}  ")
                hdr.append(f"Exec: {w.meta.executor}\n")
                # compute cores-in-use and ETA similar to curses UI
                run,que,C=classify(w)
                cores_in_use = sum([max(0.0, t.cpus) for t in run if t.cpus is not None])
                # determine total cores
                if w.ncpu <= 0:
                    try:
                        if sys.platform.startswith("darwin"):
                            w.ncpu=int(subprocess.check_output(["sysctl","-n","hw.ncpu"],text=True))
                        else:
                            w.ncpu=int(subprocess.check_output(["bash","-lc","nproc 2>/dev/null || getconf _NPROCESSORS_ONLN"],text=True))
                    except Exception:
                        w.ncpu = 0
                # fallback cores-in-use using sum of cpu_pct when cpus missing
                if cores_in_use <= 0 and w.ncpu>0:
                    approx = 0.0
                    for t in run:
                        if t.cpu_pct is not None:
                            approx += max(0.0, t.cpu_pct/100.0)
                    cores_in_use = min(float(w.ncpu), approx)
                cores_txt = f"   CORES {int(round(cores_in_use))}/{w.ncpu}" if w.ncpu>0 else ""
                hdr.append(f"CPU {w.cpu_pct}%   MEM {w.mem_pct}%   LOAD1 {w.load_1}{cores_txt}")
                if filter_state["proc"]:
                    hdr.append(f"\nFilter: {filter_state['proc']}", style="yellow")
                hdr.append("\n[dim]Keys: j/k focus • f filter • s sort • +/- log size • ]/[ tail • a toggle all-logs • q quit" + (" • ALL LOGS" if getattr(a, 'all_logs', False) else "") + "[/dim]")
                header_panel = Panel(hdr, title="nf-monitor", border_style="cyan")

                # pipeline progress bar (cumulative)
                total_so_far = C.get('cum_done',0)+C.get('cum_cached',0)+C.get('cum_failed',0)+C.get('cum_killed',0)+C['running']+C['queued']
                pct_cum = int(((C.get('cum_done',0)+C.get('cum_cached',0))*100/max(1,total_so_far))) if total_so_far else 0
                prog=Progress(TextColumn("[bold]Pipeline"), BarColumn(), TextColumn(f" {pct_cum}%"), expand=True)
                prog.add_task("pipe", total=100, completed=pct_cum)

                # remove per-process summary; show placeholder
                proc_tbl = Table(title="", expand=True, show_edge=False)

                # running table
                rt_title = "Running (all logs)" if getattr(a, 'all_logs', False) else "Running (j/k to focus)"
                rt = Table(title=rt_title, expand=True)
                rt.add_column("Module", ratio=2); rt.add_column("Sample", ratio=2); rt.add_column("id"); rt.add_column("PID"); rt.add_column("CPU"); rt.add_column("CPUS"); rt.add_column("RSS"); rt.add_column("Age"); rt.add_column("Stage", ratio=2); rt.add_column("CPU hist")
                # ensure labels from FS if missing
                for t in run:
                    ensure_task_label_from_fs(w, t)
                run_filtered = [t for t in run if (not filter_state["proc"] or (t.name or "?") == filter_state["proc"]) ]
                # apply sort option
                sort_mode = filter_state.get('sort') or 'default'
                if sort_mode == 'cpu':
                    run_filtered.sort(key=lambda t: (-(t.cpu_pct if t.cpu_pct is not None else -1), -(t.rss_mb if t.rss_mb is not None else -1), -(time.time()-t.first_ts)))
                elif sort_mode == 'rss':
                    run_filtered.sort(key=lambda t: (-(t.rss_mb if t.rss_mb is not None else -1), -(t.cpu_pct if t.cpu_pct is not None else -1), -(time.time()-t.first_ts)))
                elif sort_mode == 'age':
                    run_filtered.sort(key=lambda t: (-(time.time()-t.first_ts)))
                for idx, t in enumerate(run_filtered[:18]):
                    cpu=(f"{t.cpu_pct:.0f}%" if t.cpu_pct is not None else "--")
                    # colorize CPU and RSS
                    cpu_style = ("green" if (t.cpu_pct or 0) < 50 else ("yellow" if (t.cpu_pct or 0) < 80 else "red")) if t.cpu_pct is not None else "dim"
                    cpus = (f"{t.cpus:.1f}" if t.cpus is not None else "--")
                    rss=(f"{t.rss_mb:.0f}MB" if t.rss_mb is not None else "--")
                    rss_style = ("green" if (t.rss_mb or 0) < 2048 else ("yellow" if (t.rss_mb or 0) < 8192 else "red")) if t.rss_mb is not None else "dim"
                    # split into module (process name) and sample/tag
                    module_name = (t.name or "?")
                    sample_tag  = (t.tag or "-")
                    lab = Text(module_name, style=_name_style(module_name))
                    if (not getattr(a, 'all_logs', False)) and idx == focused["idx"]:
                        lab.stylize("bold yellow")
                    stage = (t.stage or "")
                    # sparkline from cpu_hist
                    hist = list(t.cpu_hist)
                    spark = ""
                    if hist:
                        try:
                            # map to 0-7 levels using unicode blocks
                            blocks = "▁▂▃▄▅▆▇█"
                            mx = max(1.0, max(hist))
                            spark = "".join([blocks[min(7, int((v/mx)*7))] for v in hist[-20:]])
                        except Exception:
                            spark = ""
                    rt.add_row(lab, sample_tag, short_hash(t.id or ""), str(t.pid or '-'), Text(cpu, style=cpu_style), cpus, Text(rss, style=rss_style), sec2hms(int(time.time()-t.first_ts)), stage, Text(spark, style="cyan"))

                # live log tail panel(s)
                log_panel = Panel(Text("No running tasks", style="yellow"), title="Live log", border_style="blue")
                if run_filtered:
                    if getattr(a, 'all_logs', False):
                        # stack multiple logs into one panel
                        text = Text()
                        limit = min(len(run_filtered), 6)
                        for i in range(limit):
                            t = run_filtered[i]
                            payload, tail, outs = insight(t.workdir, w.tail_n)
                            text.append(f"{_task_label(t)} (id:{short_hash(t.id or '')})\n", style="bold")
                            if payload:
                                text.append(f"▶ {payload}\n", style="cyan")
                            for ln in tail:
                                text.append(ln+"\n")
                            if outs:
                                text.append("outputs: "+", ".join([f"{n} ({sz})" for n,sz in outs]))
                            if i<limit-1:
                                text.append("\n", style="dim")
                        log_panel = Panel(text, title="Live logs", border_style="blue")
                    else:
                        # single focused log
                        if focused["idx"] >= len(run_filtered):
                            focused["idx"] = len(run_filtered)-1
                        ft = run_filtered[max(0, focused["idx"])]
                        focused["id"] = ft.id
                        payload, tail, outs = insight(ft.workdir, w.tail_n)
                        text = Text()
                        text.append(f"{_task_label(ft)} (id:{short_hash(ft.id or '')})\n", style="bold")
                        if payload:
                            text.append(f"▶ {payload}\n", style="cyan")
                        for ln in tail:
                            text.append(ln+"\n")
                        if outs:
                            text.append("outputs: "+", ".join([f"{n} ({sz})" for n,sz in outs]))
                        log_panel = Panel(text, title="Live log", border_style="blue")

                # queued table
                # remove queued table
                qt = Table(title="", expand=True)

                # errors panel
                err_txt = Text("✓ none", style="green") if not w.recent_errors else Text("\n".join(list(w.recent_errors)[-5:]), style="red")
                err_panel = Panel(err_txt, title="Recent errors", border_style="red")

                # help panel
                help_txt = Text("j/k focus • f filter by process • q quit\nLogs: right panel shows live output\nColors: per-process name", style="dim")
                help_panel = Panel(help_txt, title="Help", border_style="white")

                # layout
                layout=Layout()
                layout.split_column(
                    Layout(header_panel, size=4),
                    Layout(Panel(prog), size=3),
                    Layout(name="body")
                )
                layout["body"].split_row(
                    Layout(name="right", ratio=1)
                )
                # Compute dynamic sizes: prioritize log panel height
                log_rows = max(8, int(getattr(a, 'log_rows', 24)))
                run_rows = 14
                layout["right"].split_column(
                    Layout(rt, size=run_rows),
                    Layout(log_panel, size=log_rows),
                    Layout(err_panel, size=6),
                    Layout(help_panel, size=5)
                )
                return layout

            import threading
            kt = threading.Thread(target=_key_thread, daemon=True)
            kt.start()
            with Live(render(), refresh_per_second=max(2,int(1.0/max(0.2,a.refresh)))) as live:
                next_tick=0.0
                while True:
                    now=time.time()
                    if now>=next_tick:
                        for ln in tail_lines(w, log):
                            update_meta(w, ln); parse_line(w, ln, time.time())
                            if any(k in ln for k in ("ERROR","FAILED","Exception","Caused by:")) and "collect-file" not in ln:
                                w.recent_errors.append(ln.strip())
                        for row in tail_trace(w):
                            apply_trace_row(w, row, time.time())
                        w.cpu_pct, w.mem_pct, w.load_1 = sys_metrics()
                        live.update(render())
                        next_tick = now + max(0.3, a.refresh)
                    time.sleep(0.05)
            stop_keys["stop"] = True
            return
        except Exception:
            pass

    # Curses fallback
    if a.no_alt: curses.wrapper(loop)
    else:        curses.wrapper(loop)

if __name__=="__main__":
    main()
