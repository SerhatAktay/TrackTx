#!/usr/bin/env python3
# nf-monitor — Nextflow pipeline monitoring
# Enhanced monitoring tool for Nextflow pipelines with:
# - Sparklines & trend visualization
# - Predictive ETA & bottleneck detection
# - Smart alerts & resource tracking
# - Export capabilities (JSON/HTML/CSV/MD)
# - Multiple interactive views
# stdlib-only core, Python 3.8+, macOS/Linux

import argparse, curses, io, json, os, re, signal, subprocess, sys, time
from dataclasses import dataclass, field
from typing import Dict, Optional, List, Set, Tuple
from collections import deque, defaultdict
from datetime import datetime

# Python version compatibility
try:
    from typing import Pattern
except ImportError:
    from typing import Any as Pattern

# ═══════════════════════════════════════════════════════════════
# ─── Constants & Configuration ─────────────────────────────────
# ═══════════════════════════════════════════════════════════════

SPARKLINE_BLOCKS = "▁▂▃▄▅▆▇█"
SPARKLINE_WIDTH = 20

# ═══════════════════════════════════════════════════════════════
# ─── CLI Arguments ──────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

def parse_args():
    ap = argparse.ArgumentParser(
        description="nf-monitor — Nextflow pipeline monitoring with trends & predictions",
        epilog="Interactive keys: j/k=nav t=trends p=process e=export h=help q=quit"
    )
    ap.add_argument("--log", default=".nextflow.log", help="Path to .nextflow.log")
    ap.add_argument("--trace", default="", help="Path to trace file (auto-detect if empty)")
    ap.add_argument("--work", default="", help="Work directory (auto-detected if not specified)")
    ap.add_argument("--refresh", type=float, default=0.4, help="Refresh interval in seconds")
    ap.add_argument("--tail", type=int, default=20, help="Lines to tail from task logs")
    ap.add_argument("--oneshot", action="store_true", help="Run once and print summary")
    ap.add_argument("--json", default="", help="Export JSON snapshot path")
    ap.add_argument("--html", default="", help="Export HTML dashboard path")
    ap.add_argument("--csv", default="", help="Export CSV process stats path")
    ap.add_argument("--filter", default="", help="Regex filter for task names")
    ap.add_argument("--no-alt", action="store_true", help="Disable alternate screen (curses)")
    ap.add_argument("--simple", action="store_true", help="Force simple curses TUI (no Rich)")
    ap.add_argument("--all-logs", action="store_true", help="Show stacked logs for all running tasks")
    ap.add_argument("--log-rows", type=int, default=24, help="Log panel height for Rich UI")
    ap.add_argument("--resolve-hash", default="", help="Print process info for task hash and exit")
    ap.add_argument("--from-start", action="store_true", help="Ingest full log history")
    ap.add_argument("--retain-sec", type=int, default=900, help="Seconds to retain completed tasks (0=immediate clear)")
    ap.add_argument("--max-tasks", type=int, default=2000, help="Max tasks to keep in memory")
    return ap.parse_args()

# ═══════════════════════════════════════════════════════════════
# ─── Regex Patterns ─────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
RE_PREFIX = re.compile(
    r"(?:\[(?P<id_pre>[0-9a-f]{2}/[0-9a-f]+)\]\s+)?(?P<status>Submitted|Cached|Completed|Failed|Error executing|Killed|Terminated)\s+process\s*>\s*"
    r"(?P<name>[^\s(]+)(?:\s*\((?P<tag>[^)]+)\))?(?:\s*\[(?P<id>[0-9a-f]{2}/[0-9a-f]+)\])?", re.I
)
RE_SUFFIX = re.compile(
    r"process\s*>\s*(?P<name>[^\s(]+)(?:\s*\((?P<tag>[^)]+)\))?"
    r"(?:\s*\[(?P<id>[0-9a-f]{2}/[0-9a-f]+)\])?.*?"
    r"\b(Submitted|Cached|Completed|succeeded|Failed|FAILED|Error executing|Killed|Terminated|Retrying)\b", re.I
)
RE_HANDLER = re.compile(
    r"TaskHandler\[id:\s*\d+;\s*name:\s*(?P<name>[^;(]+)(?:\s*\((?P<tag>[^)]+)\))?;\s*status:\s*(?P<status>[^;]+);\s*.*?workDir:\s*(?P<workdir>[^\]]+)\]", re.I
)
RE_EXEC = re.compile(r"(?:^|\s)executor\s*>\s*(?P<exec>[^\s]+)", re.I)
RE_RUN = re.compile(
    r"(?:"
    r"(?:Workflow\s+run\s+name|Run\s+name)\s*:\s*"   # "Run name: xxx" (NF 25.x stdout/log)
    r"|(?:^|\s)runName\s*[:=]\s*"                    # "runName: xxx" (older NF log)
    r"|Launching\s+`[^`]+`\s+\["                     # "Launching `main.nf` [run_name]" (NF 25.x)
    r")(?P<name>[A-Za-z0-9_-]+)"
)
RE_SES = re.compile(
    r"(?:"
    r"Session\s+UUID\s*[:=]+\s*"         # "Session UUID: uuid" (NF 25.x DEBUG log)
    r"|Session\s+ID\s*[:=]+\s*"          # "Session ID: uuid" (NF stdout)
    r"|(?:^|\s)sessionId\s*[:=]+\s*"     # "sessionId: uuid" (older NF)
    r"|Session\s+id\s*[:=]+\s*"          # "Session id: uuid"
    r"|Workflow\s+session\s*[:=]\s*"     # "Workflow session: uuid"
    r"|(?:^|\s)session:\s*"              # "session: uuid"
    r")(?P<sid>[A-Za-z0-9-]+)"
)
RE_WORK = re.compile(r"(?:^|\s)(?:Work-dir|Working (?:dir|directory))\s*:\s*(?P<dir>\S+)")
RE_WORK_HANDLER = re.compile(r"\bworkDir:\s*(?P<dir>[^\]]+)\]")
RE_TASKDIR = re.compile(r".*/work/[0-9a-f]{2}/[0-9a-f]+$")
RE_CONTAINER = re.compile(r"(?:^|\s)container\s*>\s*'?(?P<container>[^'\s]+)'?", re.I)

def norm_state(tok: str) -> str:
    """Normalize state token to standard format"""
    t = tok.lower()
    if "retry" in t: return "RETRYING"
    if "submit" in t or "new" in t: return "SUBMITTED"
    if "cache" in t: return "CACHED"
    if "succeed" in t or "completed" in t: return "COMPLETED"
    if "error executing" in t or "fail" in t: return "FAILED"
    if "kill" in t or "terminate" in t: return "KILLED"
    return ""

# ═══════════════════════════════════════════════════════════════
# ─── Enhanced Data Models ───────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

@dataclass
class TaskMetrics:
    """Resource metrics with history for sparklines"""
    cpu_pct: Optional[float] = None
    rss_mb: Optional[float] = None
    cpus: Optional[float] = None  # allocated CPUs
    cpu_hist: deque = field(default_factory=lambda: deque(maxlen=SPARKLINE_WIDTH))
    rss_hist: deque = field(default_factory=lambda: deque(maxlen=SPARKLINE_WIDTH))
    
    def update(self):
        """Add current values to history"""
        if self.cpu_pct is not None:
            self.cpu_hist.append(max(0.0, float(self.cpu_pct)))
        if self.rss_mb is not None:
            self.rss_hist.append(max(0.0, float(self.rss_mb)))
    
    def cpu_sparkline(self) -> str:
        """Generate CPU usage sparkline (fixed 0-100% scale)"""
        return sparkline_pct(list(self.cpu_hist))
    
    def rss_sparkline(self) -> str:
        """Generate memory usage sparkline"""
        return sparkline(list(self.rss_hist))
    
    def cpu_trend(self) -> str:
        """Detect CPU trend: rising/falling/stable"""
        if len(self.cpu_hist) < 5:
            return "stable"
        recent = list(self.cpu_hist)[-5:]
        avg_early = sum(recent[:2]) / 2
        avg_late = sum(recent[-2:]) / 2
        diff = avg_late - avg_early
        if abs(diff) < 5:
            return "stable"
        return "rising" if diff > 0 else "falling"

@dataclass
class Task:
    """Enhanced task representation"""
    id: str
    name: str = ""
    tag: str = ""
    state: str = ""  # SUBMITTED/RUNNING/COMPLETED/FAILED/CACHED/KILLED/RETRYING
    first_ts: float = field(default_factory=time.time)
    last_ts: float = field(default_factory=time.time)
    retries: int = 0
    workdir: str = ""
    pid: Optional[int] = None
    duration_ms: Optional[int] = None
    stage: Optional[str] = None  # current stage from logs
    metrics: TaskMetrics = field(default_factory=TaskMetrics)
    slow_flag: bool = False  # marked if slower than process average
    
    def label(self) -> str:
        return f"{self.name} ({self.tag})" if self.tag else self.name
    
    def age(self) -> float:
        """Age in seconds"""
        return time.time() - self.first_ts
    
    def runtime(self) -> Optional[float]:
        """Runtime in seconds if running"""
        if self.duration_ms:
            return self.duration_ms / 1000.0
        if self.state == "RUNNING":
            return time.time() - self.first_ts
        return None

@dataclass
class ProcessStats:
    """Per-process aggregate statistics"""
    name: str
    total: int = 0
    completed: int = 0
    failed: int = 0
    cached: int = 0
    running: int = 0
    retries: int = 0
    durations: List[float] = field(default_factory=list)
    
    @property
    def avg_duration(self) -> Optional[float]:
        """Average duration in seconds"""
        if not self.durations:
            return None
        return sum(self.durations) / len(self.durations)
    
    @property
    def median_duration(self) -> Optional[float]:
        """Median duration in seconds"""
        if not self.durations:
            return None
        sorted_durs = sorted(self.durations)
        mid = len(sorted_durs) // 2
        return sorted_durs[mid]
    
    @property
    def p95_duration(self) -> Optional[float]:
        """95th percentile duration"""
        if not self.durations:
            return None
        sorted_durs = sorted(self.durations)
        idx = int(len(sorted_durs) * 0.95)
        return sorted_durs[min(idx, len(sorted_durs)-1)]
    
    def update(self, task: Task):
        """Update stats from completed task"""
        if task.state in ("COMPLETED", "CACHED") and task.duration_ms:
            self.durations.append(task.duration_ms / 1000.0)
            # Keep only recent 100 durations to avoid memory bloat
            if len(self.durations) > 100:
                self.durations = self.durations[-100:]

@dataclass
class Alert:
    """System alert/warning"""
    level: str  # INFO/WARNING/ERROR
    message: str
    task_id: Optional[str] = None
    process: Optional[str] = None
    timestamp: float = field(default_factory=time.time)
    
    def age_str(self) -> str:
        age = time.time() - self.timestamp
        if age < 60:
            return f"{int(age)}s ago"
        elif age < 3600:
            return f"{int(age/60)}m ago"
        else:
            return f"{int(age/3600)}h ago"

@dataclass
class Meta:
    """Pipeline metadata"""
    run_name: str = "?"
    session: str = "?"
    executor: str = "?"
    work_root: str = ""
    # Container / environment metadata (best-effort, may be empty)
    container_engine: str = ""
    container_image: str = ""

@dataclass
class World:
    """Complete monitoring state"""
    tasks: Dict[str, Task] = field(default_factory=dict)
    proc_stats: Dict[str, ProcessStats] = field(default_factory=dict)
    meta: Meta = field(default_factory=Meta)
    alerts: deque = field(default_factory=lambda: deque(maxlen=50))
    recent_errors: deque = field(default_factory=lambda: deque(maxlen=6))
    start_ts: float = field(default_factory=time.time)
    
    # Log/trace tailing state
    last_log_pos: int = 0
    last_log_ino: Optional[int] = None
    trace_path: Optional[str] = None
    last_trace_pos: int = 0
    last_trace_ino: Optional[int] = None
    trace_header: Optional[List[str]] = None
    
    # Configuration
    tail_n: int = 20
    filt: Optional[Pattern] = None
    retain_sec: int = 900
    max_tasks: int = 2000
    
    # System metrics
    cpu_pct: int = 0
    mem_pct: int = 0
    load_1: str = "0.00"
    ncpu: int = 0
    
    # Cumulative stats (survive task pruning)
    seen_ids: set = field(default_factory=set)
    cum_seen: int = 0
    cum_done: int = 0
    cum_cached: int = 0
    cum_failed: int = 0
    cum_killed: int = 0

    # Per-process CPU allocation learned from completed trace rows.
    # Key: process short name (e.g. "align_reads_to_genome"), value: cpus float.
    # Used to back-fill running tasks of the same type before their trace row
    # is written (NF 25.x local executor never writes NXF_TASK_CPUS to disk).
    process_cpus: Dict[str, float] = field(default_factory=dict)
    # Operating mode: full (log+trace), log_only, trace_only, fs_only
    mode: str = "fs_only"
    
    def add_alert(self, level: str, message: str, task_id: Optional[str] = None, process: Optional[str] = None):
        """Add new alert"""
        self.alerts.append(Alert(level=level, message=message, task_id=task_id, process=process))

# ═══════════════════════════════════════════════════════════════
# ─── Visualization Helpers ──────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

def sparkline(values: List[float], width: int = SPARKLINE_WIDTH) -> str:
    """Generate ASCII sparkline from values (relative scale — max value = full bar)."""
    if not values:
        return ""
    max_val = max(values) if max(values) > 0 else 1.0
    normalized = [min(7, int((v / max_val) * 7)) for v in values]
    if len(normalized) < width:
        normalized = ([0] * (width - len(normalized))) + normalized
    else:
        normalized = normalized[-width:]
    return "".join([SPARKLINE_BLOCKS[v] for v in normalized])

def sparkline_pct(values: List[float], width: int = SPARKLINE_WIDTH) -> str:
    """Generate ASCII sparkline fixed to a 0–100% scale.

    Unlike sparkline(), the tallest block always represents 100%, so a bar at
    half-height genuinely means ~50% CPU — making different tasks comparable.
    Uses SPARKLINE_BLOCKS (starting with ▁) so even idle tasks show a thin
    bar rather than invisible spaces.
    """
    if not values:
        return ""
    # Scale 0-100% to levels 0-7.  Use ceiling so even 1% shows ▁ rather than
    # being rounded away to a space.
    def _level(v: float) -> int:
        v = min(100.0, max(0.0, v))
        if v == 0.0:
            return 0
        return min(7, max(0, int((v / 100.0) * 8) - 1 + 1))  # 1..100 → 0..7

    normalized = [_level(v) for v in values]
    if len(normalized) < width:
        normalized = ([0] * (width - len(normalized))) + normalized
    else:
        normalized = normalized[-width:]
    return "".join([SPARKLINE_BLOCKS[v] for v in normalized])

def colorize_cpu_pct(pct: float) -> Tuple[str, str]:
    """Return (color_name, style) for CPU percentage"""
    if pct < 50:
        return ("green", "dim")
    elif pct < 80:
        return ("yellow", "normal")
    else:
        return ("red", "bold")

def colorize_mem_mb(mb: float) -> Tuple[str, str]:
    """Return (color_name, style) for memory MB"""
    if mb < 2048:  # < 2GB
        return ("green", "dim")
    elif mb < 8192:  # < 8GB
        return ("yellow", "normal")
    else:
        return ("red", "bold")

def format_duration(seconds: float) -> str:
    """Format duration as human-readable string"""
    if seconds < 60:
        return f"{int(seconds)}s"
    elif seconds < 3600:
        return f"{int(seconds/60)}m {int(seconds%60)}s"
    else:
        h = int(seconds / 3600)
        m = int((seconds % 3600) / 60)
        return f"{h}h {m}m"

def format_eta(seconds: float) -> str:
    """Format ETA in human-readable form"""
    if seconds < 60:
        return f"~{int(seconds)}s"
    elif seconds < 3600:
        return f"~{int(seconds/60)}m"
    else:
        return f"~{int(seconds/3600)}h"

def sec2hms(t: int) -> str:
    """Convert seconds to HH:MM:SS"""
    h = t // 3600
    m = (t % 3600) // 60
    s = t % 60
    return f"{h:02d}:{m:02d}:{s:02d}"

def truncate(s: str, w: int) -> str:
    """Truncate string to width"""
    return s if len(s) <= w else (s[:max(0, w-3)] + "...")

def short_hash(h: Optional[str]) -> str:
    """Shorten task hash for display"""
    try:
        if not h:
            return "?"
        return h.replace("/", "")[:12]
    except Exception:
        return "?"

def looks_like_hash(s: str) -> bool:
    """Check if string looks like a hash"""
    return bool(re.fullmatch(r"[0-9a-f]{24,64}", s or ""))

# ═══════════════════════════════════════════════════════════════
# ─── Analytics & Predictions ────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

def predict_completion_eta(world: World) -> Optional[float]:
    """Estimate seconds until pipeline completion based on process stats"""
    running = [t for t in world.tasks.values() if t.state == "RUNNING"]
    if not running:
        return 0.0
    
    # Use process averages to estimate remaining time
    max_remaining = 0.0
    for task in running:
        stats = world.proc_stats.get(task.name)
        if stats and stats.avg_duration:
            runtime = task.runtime() or 0
            remaining = max(0, stats.avg_duration - runtime)
            max_remaining = max(max_remaining, remaining)
    
    return max_remaining if max_remaining > 0 else None

def identify_slow_tasks(world: World) -> List[Task]:
    """Find tasks running >1.5x their process average"""
    slow = []
    for task in world.tasks.values():
        if task.state != "RUNNING":
            continue
        stats = world.proc_stats.get(task.name)
        if not stats or not stats.avg_duration:
            continue
        runtime = task.runtime()
        if runtime and runtime > stats.avg_duration * 1.5:
            slow.append(task)
    return slow

def detect_bottlenecks(world: World) -> List[str]:
    """Identify processes that are bottlenecks (high failure rate or slow)"""
    bottlenecks = []
    for name, stats in world.proc_stats.items():
        # High failure rate
        if stats.total > 5 and stats.failed / stats.total > 0.2:
            bottlenecks.append(f"{name}: high failure rate ({stats.failed}/{stats.total})")
        
        # Excessive retries
        if stats.total > 0 and stats.retries > stats.total * 0.5:
            bottlenecks.append(f"{name}: many retries ({stats.retries})")
    
    return bottlenecks

def detect_resource_issues(task: Task, world: World) -> Optional[str]:
    """Detect resource-related issues for a task"""
    m = task.metrics
    
    # Low CPU efficiency
    if m.cpus and m.cpu_pct is not None:
        efficiency = (m.cpu_pct / (m.cpus * 100)) * 100
        if efficiency < 30:
            return f"Low CPU efficiency: {efficiency:.0f}%"
    
    # Memory leak? (rising trend)
    if m.rss_hist and len(m.rss_hist) >= 10:
        recent = list(m.rss_hist)[-10:]
        trend = (recent[-1] - recent[0]) / max(1, recent[0]) * 100
        if trend > 50:  # >50% increase
            return f"Memory rising: +{trend:.0f}%"
    
    return None

# ═══════════════════════════════════════════════════════════════
# ─── Filesystem & System Metrics ────────────────────────────────
# ═══════════════════════════════════════════════════════════════

def abspath(p: str) -> str:
    return p if os.path.isabs(p) else os.path.abspath(p)

def find_nextflow_log(seed: str) -> Optional[str]:
    """Search upwards from seed directory for a .nextflow.log file"""
    try:
        d = abspath(seed or ".")
        while True:
            cand = os.path.join(d, ".nextflow.log")
            if os.path.isfile(cand):
                return cand
            parent = os.path.dirname(d)
            if not parent or parent == d:
                break
            d = parent
    except Exception:
        pass
    return None

def guess_work_root(log_path: str, cli: str) -> str:
    """Guess work directory from log or filesystem"""
    if cli:
        return abspath(cli)
    # From log
    try:
        with open(log_path, "rb") as fh:
            txt = ANSI_RE.sub("", fh.read().decode("utf-8", "ignore"))
        m = RE_WORK.search(txt)
        if m:
            return abspath(m.group("dir"))
    except Exception:
        pass
    # Heuristic: prefer ./work with hashed layout
    for candidate in ("work", "../work", "../../work"):
        c = abspath(candidate)
        if os.path.isdir(c):
            for root, dirs, files in os.walk(c):
                if RE_TASKDIR.match(root):
                    return c
            return c
    return abspath(os.environ.get("NXF_WORK") or "work")

def sys_metrics() -> Tuple[int, int, str]:
    """Get system CPU%, memory%, and load average"""
    cpu = 0
    mem = 0
    load = "0.00"
    
    try:
        if sys.platform.startswith("darwin"):
            ncpu = int(subprocess.check_output(["sysctl", "-n", "hw.ncpu"], text=True))
            vals = [float(x) for x in subprocess.check_output(["ps", "-A", "-o", "%cpu="], text=True).split() if x.strip()]
            cpu = int(round(min(100.0, max(0.0, sum(vals) / max(1, ncpu)))))
        else:
            ncpu = int(subprocess.check_output(["bash", "-lc", "nproc 2>/dev/null || getconf _NPROCESSORS_ONLN"], text=True))
            vals = [float(x) for x in subprocess.check_output(["bash", "-lc", "ps -A -o %cpu= || true"], text=True).split() if x.strip()]
            cpu = int(round(min(100.0, max(0.0, sum(vals) / max(1, ncpu)))))
    except Exception:
        pass
    
    try:
        if sys.platform.startswith("darwin"):
            total = int(subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True))
            vm = subprocess.check_output(["vm_stat"], text=True)
            pages = 0
            for ln in vm.splitlines():
                if any(k in ln for k in ("Pages free", "Pages inactive", "Pages speculative")):
                    pages += int(ln.split(":")[1].strip().strip("."))
            avail = pages * 4096
            used = max(0, total - avail)
            mem = int(round(used * 100 / max(1, total)))
        else:
            with open("/proc/meminfo") as fh:
                m = {ln.split(":")[0]: int(ln.split(":")[1].strip().split()[0]) for ln in fh}
            total = m.get("MemTotal", 0)
            avail = m.get("MemAvailable", 0)
            mem = int(round((total - avail) * 100 / max(1, total)))
    except Exception:
        pass
    
    try:
        if sys.platform.startswith("darwin"):
            up = subprocess.check_output(["uptime"], text=True)
            for key in ("load averages:", "load average:"):
                if key in up:
                    load = up.split(key)[-1].split(",")[0].strip()
        else:
            with open("/proc/loadavg") as fh:
                load = fh.read().split()[0]
    except Exception:
        pass
    
    return cpu, mem, load

def du_h(path: str) -> str:
    """Get human-readable disk usage"""
    if not os.path.isdir(path):
        return "--"
    try:
        return subprocess.check_output(["du", "-sh", path], text=True, stderr=subprocess.DEVNULL).split()[0]
    except Exception:
        return "--"

# ═══════════════════════════════════════════════════════════════
# ─── Log & Trace I/O ────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

def _filter_lines(s: str) -> List[str]:
    """Filter out debug/noise lines"""
    return [l for l in s.splitlines() if not l.lstrip().startswith("DEBUG") and "Missed collect-file cache" not in l]

def read_all_lines(path: str) -> List[str]:
    """Read entire file with ANSI filtering"""
    try:
        with open(path, "rb") as fh:
            return _filter_lines(ANSI_RE.sub("", fh.read().decode("utf-8", "ignore")))
    except Exception:
        return []

def tail_lines(world: World, path: str) -> List[str]:
    """Tail new lines from log file"""
    try:
        st = os.stat(path)
    except OSError:
        return []
    
    reopen = (world.last_log_ino is None or st.st_ino != world.last_log_ino or st.st_size < world.last_log_pos)
    if reopen:
        f = open(path, "rb")
        world.last_log_ino = st.st_ino
        world.last_log_pos = 0
    else:
        f = open(path, "rb")
        f.seek(world.last_log_pos, io.SEEK_SET)
    
    chunk = f.read()
    world.last_log_pos = f.tell()
    f.close()
    
    if not chunk:
        return []
    return _filter_lines(ANSI_RE.sub("", chunk.decode("utf-8", "ignore")))

def _parse_trace_header(line: str) -> Optional[List[str]]:
    """Parse CSV header from trace file"""
    try:
        import csv
        return next(csv.reader([line]))
    except Exception:
        return None

def _parse_trace_row(header: List[str], line: str) -> Optional[Dict[str, str]]:
    """Parse CSV row from trace file"""
    try:
        import csv
        vals = next(csv.reader([line]))
        if len(vals) != len(header):
            return None
        return {header[i]: vals[i] for i in range(len(header))}
    except Exception:
        return None

def tail_trace(world: World) -> List[Dict[str, str]]:
    """Tail new rows from trace file"""
    p = world.trace_path
    if not p or not os.path.isfile(p):
        return []
    
    try:
        st = os.stat(p)
    except OSError:
        return []
    
    reopen = (world.last_trace_ino is None or st.st_ino != world.last_trace_ino or st.st_size < world.last_trace_pos)
    if reopen:
        f = open(p, "rb")
        world.last_trace_ino = st.st_ino
        world.last_trace_pos = 0
        world.trace_header = None
    else:
        f = open(p, "rb")
        f.seek(world.last_trace_pos, io.SEEK_SET)
    
    chunk = f.read()
    world.last_trace_pos = f.tell()
    f.close()
    
    if not chunk:
        return []
    
    text = chunk.decode("utf-8", "ignore")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    out = []
    
    for ln in lines:
        if world.trace_header is None:
            hdr = _parse_trace_header(ln)
            if hdr:
                world.trace_header = hdr
            continue
        row = _parse_trace_row(world.trace_header, ln)
        if row:
            out.append(row)
    
    return out

def apply_trace_row(world: World, row: Dict[str, str], now: float):
    """Apply trace row data to world state"""
    name = (row.get("process") or row.get("name") or "").strip()
    tag = (row.get("tag") or "").strip()
    tid = (row.get("hash") or row.get("task_id") or "").strip()
    workdir = (row.get("workdir") or "").strip()
    status = (row.get("status") or "").strip().upper()
    
    t = _ensure(world, tid or _prekey(name, tag), name, tag, now, provisional=not bool(tid))
    
    if tid:
        t.id = tid
    if workdir:
        t.workdir = workdir
    if status:
        t.state = status
    
    # Metadata
    try:
        rn = (row.get("runName") or row.get("run_name") or "").strip()
        if rn:
            world.meta.run_name = rn
    except Exception:
        pass
    
    try:
        sid = (row.get("session") or row.get("session_id") or row.get("sessionId") or "").strip()
        if sid:
            world.meta.session = sid
    except Exception:
        pass
    
    # Resources
    def _parse_float(x: Optional[str]) -> Optional[float]:
        try:
            return float(x)
        except Exception:
            return None
    
    cpu = _parse_float(row.get("%cpu") or row.get("pcpu"))
    if cpu is not None:
        t.metrics.cpu_pct = cpu
    
    cpus_val = row.get("cpus") or row.get("cpu") or row.get("ncpus")
    cpus = _parse_float(cpus_val) if cpus_val is not None else None
    if cpus is not None:
        t.metrics.cpus = cpus
        # Cache per process name so running tasks of the same type can be back-filled
        if name:
            short = name.split(":")[-1]   # strip workflow prefix (e.g. "TrackTx:")
            world.process_cpus[short] = cpus
            world.process_cpus[name]  = cpus  # also store full name
    
    rss = _parse_float(row.get("rss"))
    if rss is not None:
        # Nextflow trace rss is in bytes, convert to MB
        t.metrics.rss_mb = rss / 1024 / 1024
    
    # Duration
    dur = row.get("duration") or ""
    if dur:
        try:
            t.duration_ms = int(dur)
        except Exception:
            pass
    elif row.get("realtime"):
        rt = row.get("realtime")
        try:
            parts = rt.split(":")
            if len(parts) == 3:
                h = int(parts[0])
                m = int(parts[1])
                s = float(parts[2])
                t.duration_ms = int((h * 3600 + m * 60 + s) * 1000)
        except Exception:
            pass
    
    # If this is a terminal state and we have a duration, approximate the true start time
    if status in ("COMPLETED", "CACHED", "FAILED", "KILLED") and t.duration_ms:
        try:
            start_ts = now - (t.duration_ms / 1000.0)
            if start_ts < t.first_ts:
                t.first_ts = start_ts
        except Exception:
            pass

    # Update metrics history
    t.metrics.update()

# ═══════════════════════════════════════════════════════════════
# ─── Parser & State Machine ─────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

def _prekey(name: str, tag: str) -> str:
    """Generate provisional key for tasks without hash yet"""
    return f"_pre:{name}|{tag}"

def _ensure(world: World, key: str, name: str, tag: str, now: float, provisional: bool) -> Task:
    """Get or create task"""
    t = world.tasks.get(key)
    if not t:
        t = Task(id=key, name=name or "", tag=tag or "", first_ts=now, last_ts=now)
        t.state = "SUBMITTED" if provisional else t.state
        world.tasks[key] = t
    else:
        if name and not t.name:
            t.name = name
        if tag and not t.tag:
            t.tag = tag
        t.last_ts = now
    return t

def _apply(world: World, name: str, tag: str, tid: str, token: str, now: float):
    """Apply state change from log line"""
    st = norm_state(token)
    
    if tid:
        t = _ensure(world, tid, name, tag, now, provisional=False)
        
        # Cumulative tracking
        if tid not in world.seen_ids:
            world.seen_ids.add(tid)
            world.cum_seen += 1
        
        if st == "RETRYING":
            t.retries += 1
            # Update process stats
            if t.name:
                if t.name not in world.proc_stats:
                    world.proc_stats[t.name] = ProcessStats(name=t.name)
                world.proc_stats[t.name].retries += 1
        
        if st:
            old_state = t.state
            t.state = st
            
            # Alert on failures
            if st == "FAILED" and old_state != "FAILED":
                world.add_alert("ERROR", f"Task failed: {t.label()}", task_id=t.id, process=t.name)
            elif st == "KILLED" and old_state != "KILLED":
                world.add_alert("WARNING", f"Task killed: {t.label()}", task_id=t.id, process=t.name)
        
        # Merge provisional task if exists
        pk = _prekey(name, tag)
        if pk in world.tasks:
            old = world.tasks.pop(pk)
            t.first_ts = min(t.first_ts, old.first_ts)
    else:
        pk = _prekey(name, tag)
        t = _ensure(world, pk, name, tag, now, provisional=True)
        if st:
            t.state = st

def update_meta(world: World, line: str):
    """Extract metadata from log line"""
    m = RE_EXEC.search(line)
    if m:
        world.meta.executor = m.group("exec")
    
    m = RE_RUN.search(line)
    if m:
        world.meta.run_name = m.group("name")
    
    m = RE_SES.search(line)
    if m:
        world.meta.session = m.group("sid")
    # Fallback: NF session IDs are always UUIDs.  If a session-related log line
    # carries a UUID that the structured regex missed, grab it directly.
    elif (world.meta.session in ("", "?") and "session" in line.lower()):
        m_uuid = re.search(
            r'\b([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})\b',
            line, re.I
        )
        if m_uuid:
            world.meta.session = m_uuid.group(1)
    
    # Only match global Work-dir declaration (not TaskHandler workDir)
    m = RE_WORK.search(line)
    if m and "TaskHandler" not in line:
        cand = m.group("dir").rstrip("]").rstrip()
        cand = abspath(cand)
        # Normalize: if this looks like a task dir .../work/aa/hash, lift to .../work
        try:
            parts = cand.rstrip("/").split("/")
            if len(parts) >= 3 and parts[-3].endswith("work") and re.fullmatch(r"[0-9a-f]{2}", parts[-2] or "") and re.fullmatch(r"[0-9a-f]+", parts[-1] or ""):
                cand = "/".join(parts[:-2])
        except Exception:
            pass
        world.meta.work_root = cand

    # Container engine / image (best-effort heuristics)
    m = RE_CONTAINER.search(line)
    if m:
        world.meta.container_image = m.group("container")
    if not world.meta.container_engine:
        low = line.lower()
        # Only match on lines that definitively declare the engine — not any mention of the word.
        # "No docker container", "checking for docker", etc. must not trigger.
        if re.search(r"container\s*engine\s*[:=]\s*docker|executor\s*>\s*docker|docker\s*executor", low):
            world.meta.container_engine = "docker"
        elif re.search(r"container\s*engine\s*[:=]\s*singularity|singularity\s*executor|executor\s*>\s*singularity", low):
            world.meta.container_engine = "singularity"
        elif re.search(r"container\s*engine\s*[:=]\s*podman|podman\s*executor|executor\s*>\s*podman", low):
            world.meta.container_engine = "podman"
        elif re.search(r"activating\s+conda|conda\s*env(?:ironment)?\s*[:=]|executor\s*>\s*conda|conda.server", low):
            world.meta.container_engine = "conda"

def parse_line(world: World, line: str, now: float):
    """Parse log line for task events"""
    # Try TaskHandler first (rich info)
    m = RE_HANDLER.search(line)
    if m:
        wd = m.group("workdir").strip()
        # Extract hash from workdir if possible
        tid = ""
        if "/work/" in wd:
            parts = wd.split("/work/")[-1].split("/")
            if len(parts) >= 2:
                tid = f"{parts[0]}/{parts[1]}"
        
        name = (m.group("name") or "").strip()
        tag = (m.group("tag") or "").strip()
        status = m.group("status")
        
        _apply(world, name, tag, tid, status, now)
        
        # Also update workdir if we have the task
        if tid and tid in world.tasks:
             world.tasks[tid].workdir = wd
        return

    m = RE_PREFIX.search(line)
    if m:
        tid = (m.group("id") or m.group("id_pre") or "").strip()
        status = m.group("status")
        _apply(world, (m.group("name") or "").strip(), (m.group("tag") or "").strip(), tid, status, now)
        return
    
    m = RE_SUFFIX.search(line)
    if m:
        tail = line[m.end()-40:].lower()
        token = "Completed" if ("succeed" in tail or "completed" in tail) else tail
        _apply(world, (m.group("name") or "").strip(), (m.group("tag") or "").strip(), (m.group("id") or "").strip(), token, now)

def bootstrap(world: World, log_path: str):
    """Bootstrap world state from full log"""
    for ln in read_all_lines(log_path):
        update_meta(world, ln)
        parse_line(world, ln, time.time())
        
        if any(k in ln for k in ("ERROR", "FAILED", "Exception", "Caused by:")) and "collect-file" not in ln:
            world.recent_errors.append(ln.strip())
    
    try:
        st = os.stat(log_path)
        world.last_log_ino = st.st_ino
        world.last_log_pos = st.st_size
    except OSError:
        pass

# ═══════════════════════════════════════════════════════════════
# ─── Filesystem Probe (Running Tasks) ──────────────────────────
# ═══════════════════════════════════════════════════════════════

def find_work_roots(seed: str) -> List[str]:
    """Find work directories"""
    roots = set()
    if os.path.isdir(seed):
        roots.add(abspath(seed))
    
    for base in (os.getcwd(), os.path.dirname(os.getcwd())):
        cand = os.path.join(base, "work")
        if os.path.isdir(cand):
            roots.add(abspath(cand))
    
    return list(sorted(roots))

def label_from_dir(d: str) -> Tuple[str, str]:
    """Extract name, tag from .command.env/.command.begin"""
    def get_env(k: str, f: str) -> str:
        try:
            for ln in open(f, "r", encoding="utf-8", errors="ignore"):
                if ln.startswith(k + "="):
                    val = ln.split("=", 1)[1].strip().strip('"')
                    return val
        except Exception:
            pass
        return ""
    
    env = os.path.join(d, ".command.env")
    name = get_env("NXF_PROCESS", env) or get_env("NXF_TASK_NAME", env)
    tag = get_env("NXF_TASK_TAG", env)
    
    if not name:
        for src in (".command.begin", ".command.log"):
            p = os.path.join(d, src)
            try:
                for ln in open(p, "r", encoding="utf-8", errors="ignore"):
                    if "process >" in ln:
                        nm = re.sub(r'.*process\s*>\s*', '', ln).strip()
                        nm = re.sub(r'[\s(].*', '', nm)
                        name = nm
                        m = re.search(r'\(([^)]+)\)', ln)
                        if m:
                            tag = m.group(1)
                        break
            except Exception:
                pass
    
    return name or os.path.basename(d), tag or ""

def pid_from_dir(d: str) -> Optional[int]:
    """Read PID of the task running in work dir d.

    With scratch=false: .command.pid lives in d itself.
    With scratch=true:  Nextflow runs the task from a /tmp/nxf-* temp dir and
    writes .command.pid there.  We find the matching process via /proc by
    looking for a bash wrapper running .command.run whose cwd or script body
    references the NFS work dir hash.
    """
    # Fast path: .command.pid in the NFS work dir (scratch=false)
    for cand in (".command.pid", ".nxf.pid"):
        p = os.path.join(d, cand)
        try:
            txt = open(p).read().strip()
            if txt.isdigit():
                return int(txt)
        except Exception:
            pass

    # Slow path (Linux only): scan /proc for bash processes running .command.run
    # whose scratch cwd or script body references our work dir.
    if not os.path.isdir("/proc"):
        return None

    thash      = os.path.basename(d)          # long hex hash, e.g. "c1a2b3d4..."
    work_norm  = os.path.normpath(d)
    # Only bother if the hash looks like a real NF hash (>=8 hex chars)
    if len(thash) < 8 or not re.fullmatch(r"[0-9a-f]+", thash, re.I):
        return None

    try:
        for entry in os.scandir("/proc"):
            if not entry.name.isdigit():
                continue
            pid = int(entry.name)
            try:
                with open(f"/proc/{pid}/cmdline", "rb") as f:
                    cmdline = f.read(512).decode("utf-8", "replace").replace("\x00", " ")
                if ".command.run" not in cmdline:
                    continue

                # Read cwd symlink — this is the scratch temp dir
                try:
                    cwd = os.readlink(f"/proc/{pid}/cwd")
                except OSError:
                    continue

                # Case 1: cwd IS the NFS work dir (scratch=false, PID file missing)
                if os.path.normpath(cwd) == work_norm:
                    for pf in (".command.pid", ".nxf.pid"):
                        try:
                            txt = open(os.path.join(cwd, pf)).read().strip()
                            if txt.isdigit(): return int(txt)
                        except Exception:
                            pass
                    return pid

                # Case 2: cwd is a scratch temp dir — check its .command.run for
                # a reference to our work dir hash (NF writes copy-back paths)
                run_file = os.path.join(cwd, ".command.run")
                try:
                    with open(run_file, "r", errors="ignore") as f:
                        head = f.read(4096)
                    if thash in head or work_norm in head:
                        # Found the scratch dir for our task — get PID from there
                        for pf in (".command.pid", ".nxf.pid"):
                            try:
                                txt = open(os.path.join(cwd, pf)).read().strip()
                                if txt.isdigit(): return int(txt)
                            except Exception:
                                pass
                        return pid   # fall back to wrapper PID itself
                except Exception:
                    pass
            except Exception:
                continue
    except Exception:
        pass
    return None

def cpus_from_dir(d: str) -> Optional[float]:
    """Read allocated CPUs from task work dir files.

    NF 25.x local executor no longer writes NXF_TASK_CPUS as a static variable.
    The beforeScript in nextflow.config renders task.cpus into .command.sh as
    thread-count env vars (OMP_NUM_THREADS, BOWTIE2_THREADS, etc.), so we look
    for those.  Fallback patterns cover older NF versions.
    """
    # Primary: NF 25.x renders task.cpus into .command.sh via beforeScript
    # e.g.  export OMP_NUM_THREADS=4   export BOWTIE2_THREADS=4
    _thread_re = re.compile(
        r"^export\s+(?:"
        r"OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS"
        r"|BOWTIE2_THREADS|SAMTOOLS_THREADS|BEDTOOLS_THREADS"
        r"|STAR_THREADS|HISAT2_THREADS|MINIMAP2_THREADS"
        r"|PIGZ_THREADS|KALLISTO_THREADS|SALMON_THREADS"
        r")\s*=\s*[\"']?(?P<val>\d+)[\"']?",
        re.I
    )
    # Legacy: NXF_TASK_CPUS, nxf_cpus, etc. (older NF versions)
    _cpu_re = re.compile(
        r"^(?:export\s+)?(?:"
        r"NXF_CPUS|NXF_TASK_CPUS|NXF_TASK_PROCESS_CPUS"
        r"|nxf_cpus|nxf_task_cpus|nxf_num_cpus|task\.cpus"
        r")\s*=\s*[\"']?(?P<val>\d+(?:\.\d+)?)[\"']?",
        re.I
    )
    # Scan .command.sh first (beforeScript output), then fallbacks
    for fname in (".command.sh", ".command.env", ".command.run"):
        fpath = os.path.join(d, fname)
        try:
            with open(fpath, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    stripped = ln.strip()
                    m = _thread_re.match(stripped) or _cpu_re.match(stripped)
                    if m:
                        try:
                            return float(m.group("val"))
                        except (ValueError, TypeError):
                            pass
        except Exception:
            pass
    return None

def ps_for_pid(pid: int) -> Tuple[Optional[float], Optional[float]]:
    """Get %cpu, rss_mb for PID"""
    try:
        if sys.platform.startswith("darwin"):
            out = subprocess.check_output(["ps", "-p", str(pid), "-o", "%cpu=,rss="], text=True)
        else:
            out = subprocess.check_output(["ps", "-p", str(pid), "-o", "%cpu=,rss="], text=True)
        parts = [p for p in out.strip().split() if p]
        if len(parts) >= 2:
            cpu = float(parts[0])
            rss_kb = float(parts[1])
            return cpu, rss_kb / 1024.0
    except Exception:
        pass
    return None, None

def _ps_snapshot() -> List[Tuple[int, int, float, float, str]]:
    """Return a lightweight process snapshot [(pid, ppid, pcpu, rss_mb, cmdline)]"""
    out = []
    try:
        if sys.platform.startswith("darwin"):
            txt = subprocess.check_output(["ps", "-Ao", "pid=,ppid=,pcpu=,rss=,command="], text=True)
        else:
            txt = subprocess.check_output(["bash", "-lc", "ps -Ao pid=,ppid=,pcpu=,rss=,command= || true"], text=True)
        for ln in txt.splitlines():
            try:
                parts = ln.strip().split(maxsplit=4)
                if len(parts) >= 5:
                    pid = int(parts[0])
                    ppid = int(parts[1])
                    pcpu = float(parts[2])
                    rss_kb = float(parts[3])
                    cmd = parts[4]
                    out.append((pid, ppid, pcpu, rss_kb / 1024.0, cmd))
            except Exception:
                continue
    except Exception:
        pass
    return out

def _tree_cpu_rss(root_pid: int, snap: Optional[List] = None) -> Tuple[float, float]:
    """Sum CPU% and max RSS (MB) across root_pid and all its descendants.

    The PID stored in .command.pid is the bash wrapper that runs .command.run.
    The actual compute (bowtie2, samtools, python, …) runs as a child of that
    wrapper, so we must walk the tree to get meaningful CPU/RSS numbers.
    """
    if snap is None:
        snap = _ps_snapshot()
    by_ppid: Dict[int, List] = {}
    for p, pp, c, r, _ in snap:
        by_ppid.setdefault(pp, []).append((p, c, r))

    total_cpu = 0.0
    max_rss   = 0.0
    queue: List[int] = [root_pid]
    seen:  Set[int]  = set()
    while queue:
        pid = queue.pop(0)
        if pid in seen:
            continue
        seen.add(pid)
        for p, pp, c, r, _ in snap:
            if p == pid:
                total_cpu += c
                max_rss    = max(max_rss, r)
                break
        queue.extend(cp for cp, _, _ in by_ppid.get(pid, []))

    return total_cpu, max_rss

def docker_containers_map() -> Dict[str, str]:
    """Map work directories to Docker container IDs"""
    mapping = {}
    try:
        out = subprocess.check_output(
            ["docker", "ps", "-q"],
            text=True, stderr=subprocess.DEVNULL
        )
        container_ids = [cid.strip() for cid in out.strip().splitlines() if cid.strip()]

        for container_id in container_ids:
            try:
                cmd_out = subprocess.check_output(
                    ["docker", "inspect", container_id, "--format", "{{.Path}} {{join .Args \" \"}}"],
                    text=True, stderr=subprocess.DEVNULL
                )
                # Match work dir: .../work/xx/hash (with optional trailing / or .command.run)
                match = re.search(r'(/[^\s]+/work/[0-9a-f]{2}/[0-9a-f]+)', cmd_out)
                if match:
                    workdir = os.path.normpath(match.group(1).rstrip("/"))
                    mapping[workdir] = container_id
            except Exception:
                continue
    except Exception:
        pass
    return mapping

def docker_stats_for_container(container_id: str) -> Tuple[Optional[float], Optional[float]]:
    """Get %cpu, rss_mb for Docker container"""
    try:
        # Get stats without stream (single sample)
        out = subprocess.check_output(
            ["docker", "stats", "--no-stream", "--format", "{{.CPUPerc}}\t{{.MemUsage}}", container_id],
            text=True, stderr=subprocess.DEVNULL, timeout=2
        )
        parts = out.strip().split('\t')
        if len(parts) >= 2:
            # CPU like "12.34%"
            cpu_str = parts[0].strip().rstrip('%')
            cpu = float(cpu_str) if cpu_str else 0.0
            
            # MemUsage like "123.4MiB / 8GiB"
            mem_str = parts[1].split('/')[0].strip()
            rss_mb = 0.0
            if 'GiB' in mem_str or 'GB' in mem_str:
                rss_mb = float(mem_str.replace('GiB', '').replace('GB', '').strip()) * 1024
            elif 'MiB' in mem_str or 'MB' in mem_str:
                rss_mb = float(mem_str.replace('MiB', '').replace('MB', '').strip())
            elif 'KiB' in mem_str or 'KB' in mem_str:
                rss_mb = float(mem_str.replace('KiB', '').replace('KB', '').strip()) / 1024
            
            return cpu, rss_mb
    except Exception:
        pass
    return None, None

def running_dirs(work_roots: List[str]) -> List[str]:
    """Find active task directories (no .exitcode)"""
    out = []
    for root in work_roots:
        try:
            for aa in os.listdir(root):
                if len(aa) != 2:
                    continue
                p1 = os.path.join(root, aa)
                if not os.path.isdir(p1):
                    continue
                for bb in os.listdir(p1):
                    d = os.path.join(p1, bb)
                    if not os.path.isdir(d):
                        continue
                    if os.path.exists(os.path.join(d, ".exitcode")):
                        continue
                    if os.path.exists(os.path.join(d, ".command.run")) or os.path.exists(os.path.join(d, ".command.sh")):
                        out.append(d)
        except Exception:
            pass
    
    out.sort(key=lambda x: os.path.getmtime(x) if os.path.exists(x) else 0, reverse=True)
    return out

# ═══════════════════════════════════════════════════════════════
# ─── Introspection (Payload/Logs/Outputs) ──────────────────────
# ═══════════════════════════════════════════════════════════════

WRAP_DROP = [
    re.compile(r"^\+{1,}"), re.compile(r"^set -"), re.compile(r"^trap "),
    re.compile(r"^ulimit "), re.compile(r"^exec > "), re.compile(r"^nxf_"),
    re.compile(r"^NXF_"), re.compile(r"^CAPSULE:"), re.compile(r"^Picked up _JAVA_OPTIONS"),
    re.compile(r"^Warning: .*illegal reflective access"), re.compile(r"^INFO +\(.*Nextflow.*\)"),
    re.compile(r"read -t \d+ -r DONE"),
    re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=.*"),  # VAR=value or VAR=$(...) - skip variable assignments
]

def keep_line(s: str) -> bool:
    """Filter out wrapper noise from task logs"""
    s = ANSI_RE.sub("", s.rstrip("\r"))
    for rx in WRAP_DROP:
        if rx.search(s):
            return False
    if s.lstrip().startswith("echo "):
        return False
    return bool(s.strip())

def slurp_tail(path: str, n: int, scrub=True) -> List[str]:
    """Tail last n lines from file"""
    try:
        with open(path, "rb") as fh:
            fh.seek(0, io.SEEK_END)
            size = fh.tell()
            block = 4096
            data = b""
            while len(data.splitlines()) <= n + 1 and fh.tell() > 0:
                off = max(0, fh.tell() - block)
                fh.seek(off)
                data = fh.read(size - off) + data
                fh.seek(off)
                if off == 0:
                    break
                block *= 2
        lines = [ANSI_RE.sub("", l) for l in data.decode("utf-8", "ignore").splitlines()[-n:]]
        return [l for l in lines if keep_line(l)] if scrub else lines
    except Exception:
        return []

def payload_from(workdir: str) -> str:
    """Extract actual command from task directory (skip var assignments, prefer real commands)"""
    # Shell noise: variable assignments, control keywords, if/while/for conditions
    PAYLOAD_SKIP = re.compile(
        r"^(?:\s*#.*|"  # pure comments like '# Get final file sizes'"
        r"[A-Za-z_][A-Za-z0-9_]*=.*|\s*(?:fi|done|else|elif|then|esac)\s*$|"
        r"\s*if\s+\[\[|\s*if\s+\[|\s*elif\s+\[\[|\s*elif\s+\[|\s*while\s+\[\[|\s*while\s+\[|\s*for\s+\w+\s+in)"
    )
    for fn in (".command.sh", ".command.run"):
        p = os.path.join(workdir, fn)
        try:
            lines = []
            for ln in open(p, "r", encoding="utf-8", errors="ignore"):
                if keep_line(ln):
                    lines.append(ln.strip())
            # Prefer last line that looks like a real command
            for ln in reversed(lines):
                if ln and not PAYLOAD_SKIP.match(ln):
                    return ln
            if lines:
                return lines[-1]
        except Exception:
            pass
    return "<no payload detected>"

def best_log_for_tail(d: str) -> Optional[str]:
    """Find best log file to tail"""
    cand = [os.path.join(d, ".command.out"), os.path.join(d, ".command.err"), os.path.join(d, ".command.log")]
    cand = [p for p in cand if os.path.exists(p)]
    if not cand:
        return None
    cand.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return cand[0]

def newest_user_log(d: str) -> Optional[str]:
    """Find newest user-generated log file"""
    try:
        # Prefer known TrackTx module logs when present
        tracktx_candidates = [
            "preprocess_reads.log",
            "align_reads.log",
            "tracks.log",
            "quantify_reads_per_gene.log",
            "divergent.log",
            "qc.log",
            "functional_regions.log",
            "pol_metrics.log",
            "aggregate.log",
            "normalize_coverage_tracks.log",
            "combine.log",
        ]
        best = None
        mt = -1.0
        for nm in tracktx_candidates:
            p = os.path.join(d, nm)
            if os.path.isfile(p):
                try:
                    m = os.path.getmtime(p)
                except OSError:
                    continue
                if m > mt:
                    best, mt = p, m
        if best:
            return best

        # Generic heuristic: newest non-wrapper log-like file
        newest = None
        mt = -1.0
        for nm in os.listdir(d):
            if nm.startswith(".command.") or nm == ".exitcode":
                continue
            if any(s in nm.lower() for s in ("stdout", "stderr", ".log", ".err", ".out")):
                p = os.path.join(d, nm)
                try:
                    m = os.path.getmtime(p)
                except OSError:
                    continue
                if m > mt:
                    newest, mt = p, m
        return newest
    except Exception:
        return None

def list_outputs(d: str, k=3) -> List[Tuple[str, str]]:
    """List output files with sizes"""
    out = []
    try:
        for nm in os.listdir(d):
            if nm.startswith("."):
                continue
            p = os.path.join(d, nm)
            if os.path.isfile(p):
                st = os.stat(p)
                out.append((nm, st.st_size, st.st_mtime))
        out.sort(key=lambda x: (-x[2], -x[1]))
        
        def human(n):
            for u in ("B", "KB", "MB", "GB", "TB", "PB"):
                if n < 1024:
                    return f"{n:.1f}{u}"
                n /= 1024
            return f"{n:.1f}EB"
        
        return [(n, human(sz)) for n, sz, _ in out[:k]]
    except Exception:
        return []

# Stage detection from logs
STEP_RE = re.compile(r"\bStep\s+\d+\s*/\s*\d+\s*:\s*(.+)")
# TrackTx module format: "PREP | START |", "TRACKS | BIGWIG |", "DIVERGENT | CONFIG |", etc.
TRACKTX_STAGE_RE = re.compile(r"^([A-Za-z_]+\s*\|\s*[A-Za-z_]+\s*\|.*)$")

def detect_stage_from_tail(lines: List[str]) -> Optional[str]:
    """Extract current processing stage from logs"""
    for ln in reversed(lines or []):
        m = STEP_RE.search(ln)
        if m:
            return m.group(0).strip()
    for ln in reversed(lines or []):
        if "INFO  Step" in ln:
            return ln.strip()
    # TrackTx modules: PREP |, TRACKS |, DIVERGENT |, GTF |, ALIGN |, POL |, QC |, NORM |, etc.
    for ln in reversed(lines or []):
        s = ln.strip()
        if TRACKTX_STAGE_RE.match(s) and len(s) < 120:
            return s
    # Fallback: last non-empty line that looks like progress (contains |)
    for ln in reversed(lines or []):
        s = ln.strip()
        if s and " | " in s and len(s) < 100:
            return s
    return None

def extract_timestamp_from_log(lines: List[str]) -> Optional[str]:
    """Extract most recent timestamp from log (ts=YYYY-MM-DDTHH:MM:SSZ); prefer over script payload."""
    for ln in reversed(lines or []):
        m = re.search(r"ts=(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)", ln)
        if m:
            ts = m.group(1)
            if "COMPLETE" in ln:
                return f"Completed {ts}"
            if "START" in ln:
                return f"Started {ts}"
            return ts
    return None

def insight(d: str, tail_n: int) -> Tuple[str, List[str], List[Tuple[str, str]]]:
    """Bundle introspection: payload, log tail, outputs"""
    user = newest_user_log(d)
    src = user or best_log_for_tail(d)
    # Always try to tail something so the log pane isn't empty
    tl = slurp_tail(src, tail_n, scrub=True) if src else ["<no output yet>"]
    # If strict filtering removed everything:
    #  - For TrackTx module logs, fall back to a raw tail so we keep all PREP/ALIGN/TRACKS lines.
    #  - For Nextflow wrapper logs ('.command.*'), avoid dumping bash noise and instead show a
    #    clear placeholder so the UI isn't misleadingly empty.
    if src and not tl:
        base = os.path.basename(src)
        if base.startswith(".command."):
            tl = ["<waiting for module logs>"]
        else:
            tl = slurp_tail(src, tail_n, scrub=False)
    # Prefer actual timestamp or stage line from module logs over script commands
    ts = extract_timestamp_from_log(tl)
    stage = detect_stage_from_tail(tl)
    if ts:
        payload = ts
    elif stage:
        payload = stage
    else:
        base = os.path.basename(src) if src else ""
        if base.startswith(".command."):
            # For wrapper logs, prefer a meaningful tool/progress line from the tail
            # (e.g. bam_sort_core, bowtie2, samtools) and avoid echoing shell glue.
            meaningful = None
            for ln in reversed(tl):
                s = (ln or "").strip()
                if not s:
                    continue
                if s.startswith("#"):
                    continue
                if s.startswith("sed "):
                    continue
                meaningful = s
                break
            if meaningful:
                payload = meaningful
            else:
                # Nothing but wrapper noise so far
                payload = "<waiting for module logs>"
        else:
            payload = payload_from(d)
    outs = list_outputs(d)
    return payload, tl, outs

def ensure_task_label_from_fs(w: World, t: Task):
    """Populate task name/tag/workdir from filesystem if missing"""
    try:
        d = t.workdir
        if (not d) and t.id and "/" in t.id and (w.meta.work_root or ""):
            d = os.path.join(w.meta.work_root, t.id)
        if d and os.path.isdir(d):
            if not (t.name and t.tag):
                name, tag = label_from_dir(d)
                if name and not t.name:
                    t.name = name
                if tag and not t.tag:
                    t.tag = tag
            if not t.workdir:
                t.workdir = d
    except Exception:
        pass

# ═══════════════════════════════════════════════════════════════
# ─── Task Classification (Merge Log + FS) ──────────────────────
# ═══════════════════════════════════════════════════════════════

def classify(world: World):
    """Classify tasks and enrich with filesystem probe"""
    # Start with log/trace-based states
    done = fail = cache = killed = 0
    run = []
    queued = []
    
    for t in world.tasks.values():
        st = t.state or ""
        if st in ("COMPLETED", "SUCCEEDED"):
            done += 1
        elif st == "FAILED":
            fail += 1
        elif st == "CACHED":
            cache += 1
        elif st == "KILLED":
            killed += 1
        else:
            queued.append(t)
    
    # Overlay filesystem: authoritative running set
    roots = find_work_roots(world.meta.work_root or "work")
    active_dirs = running_dirs(roots)
    seen_ids = set()
    
    # Get Docker container mapping only when Docker is actually in use.
    # Calling `docker ps` on every classify() cycle is expensive and noisy for
    # conda/local runs where Docker is not involved.
    docker_map = docker_containers_map() if world.meta.container_engine == "docker" else {}
    
    for d in active_dirs:
        tid = "/".join(d.rstrip("/").split("/")[-2:])
        seen_ids.add(tid)
        
        # Try to find task by full hash or by short hash prefix
        t = world.tasks.get(tid)
        if not t:
            # Try matching by short hash prefix (log uses short, filesystem uses full)
            short_tid = tid[:9]  # e.g., "a6/124fb1" from "a6/124fb197de9c92cb7845294a2788ae"
            for key, task in world.tasks.items():
                if key.startswith(short_tid) or short_tid.startswith(key[:9]):
                    t = task
                    # Update task ID to full hash for consistency
                    t.id = tid
                    # Also update the dictionary key
                    world.tasks[tid] = t
                    if key != tid:
                        try:
                            del world.tasks[key]
                        except:
                            pass
                    break
        
        if not t:
            # Construct from dir
            name, tag = label_from_dir(d)
            t = Task(id=tid, name=name, tag=tag, state="RUNNING", first_ts=os.path.getmtime(d), last_ts=time.time())
            world.tasks[tid] = t
        
        t.workdir = d
        t.state = "RUNNING"
        
        # Try Docker stats first (for Docker executor)
        got_metrics = False
        d_norm = os.path.normpath(d.rstrip("/"))
        container_id = docker_map.get(d_norm) or docker_map.get(d)
        if container_id:
            cpu, rss = docker_stats_for_container(container_id)
            if cpu is not None:
                t.metrics.cpu_pct = cpu
                got_metrics = True
            if rss is not None:
                t.metrics.rss_mb = rss
                got_metrics = True

        # Fallback to PID-based monitoring (for local executor)
        if not got_metrics:
            pid = pid_from_dir(d)
            t.pid = pid

            if pid:
                # The PID in .command.pid is the bash wrapper (.command.run).
                # The actual tool (bowtie2, samtools, …) runs as a child.
                # Walk the whole process subtree so we capture real CPU/RSS.
                snap = _ps_snapshot()
                cpu, rss = _tree_cpu_rss(pid, snap)
                t.metrics.cpu_pct = cpu
                t.metrics.rss_mb  = rss

        # Fallback: try to match process by command line when PID file is absent
        if t.state == "RUNNING" and t.pid is None:
            snap = _ps_snapshot()
            # build children index by ppid
            children: Dict[int, List] = {}
            for pid, ppid, pcpu, rss_mb, cmd in snap:
                if ppid not in children: children[ppid] = []
                children[ppid].append((pid, pcpu, rss_mb, cmd))

            thash = (t.id or '').split('/')[-1]
            for pid, ppid, pcpu, rss_mb, cmd in snap:
                hit = False
                if thash and thash in cmd: hit = True
                elif t.workdir and t.workdir in cmd: hit = True

                if hit:
                    # if wrapper, try to find hottest child
                    if any(k in cmd for k in ("/bash", "bash", "/sh", "sh", "nextflow", "java")) and pid in children:
                        kids = children[pid]
                        kids.sort(key=lambda x: x[1], reverse=True)
                        if kids:
                            pid, pcpu, rss_mb, cmd = kids[0]

                    t.pid = pid
                    if t.metrics.cpu_pct is None: t.metrics.cpu_pct = pcpu
                    if t.metrics.rss_mb  is None: t.metrics.rss_mb  = rss_mb
                    break
        
        # Learn allocated CPUs if still unknown.
        # Priority: (1) task file scan, (2) per-process cache from completed trace rows
        if t.metrics.cpus is None:
            cc = cpus_from_dir(d)
            if cc is not None:
                t.metrics.cpus = cc
            else:
                # NF 25.x local executor does not write NXF_TASK_CPUS to any file.
                # Fall back to the per-process cache populated from completed trace rows.
                short = (t.name or "").split(":")[-1]
                cached = world.process_cpus.get(t.name) or world.process_cpus.get(short)
                if cached is not None:
                    t.metrics.cpus = cached
        
        # Update metrics history
        t.metrics.update()
        
        # Parse stage from latest logs
        payload, tl, _ = insight(d, world.tail_n)
        st = detect_stage_from_tail(tl)
        if st:
            t.stage = st
        
        # Check for resource issues
        issue = detect_resource_issues(t, world)
        if issue:
            world.add_alert("INFO", f"{t.label()}: {issue}", task_id=t.id, process=t.name)
    
    # Remove stale running tasks
    try:
        for key, t in list(world.tasks.items()):
            if t.state == "RUNNING" and t.id and "/" in t.id:
                if t.id not in seen_ids:
                    del world.tasks[key]
    except Exception:
        pass
    
    # Re-classify after overlay
    now = time.time()
    run = [t for t in world.tasks.values() if t.state == "RUNNING"]
    queued = [t for t in world.tasks.values() if t.state not in ("RUNNING", "COMPLETED", "FAILED", "CACHED", "KILLED")]
    
    # Prune terminal tasks
    terminal = [(k, t) for k, t in world.tasks.items() if t.state in ("COMPLETED", "FAILED", "CACHED", "KILLED")]
    
    if world.retain_sec == 0:
        for k, _ in terminal:
            try:
                del world.tasks[k]
            except Exception:
                pass
        terminal = []
    else:
        for k, t in terminal:
            if world.retain_sec > 0 and (now - t.last_ts) > world.retain_sec:
                try:
                    del world.tasks[k]
                except Exception:
                    pass
    
    # Cap-based prune
    if len(world.tasks) > world.max_tasks:
        term = [(k, t) for k, t in world.tasks.items() if t.state in ("COMPLETED", "FAILED", "CACHED", "KILLED")]
        term.sort(key=lambda kv: kv[1].last_ts)
        for k, _ in term:
            if len(world.tasks) <= world.max_tasks:
                break
            try:
                del world.tasks[k]
            except Exception:
                pass
    
    # Ensure labels from fs
    for t in run + queued:
        if not (t.name or t.tag):
            ensure_task_label_from_fs(world, t)
    
    # Update process stats
    for t in world.tasks.values():
        if not t.name:
            continue
        if t.name not in world.proc_stats:
            world.proc_stats[t.name] = ProcessStats(name=t.name)
        
        ps = world.proc_stats[t.name]
        ps.total = max(ps.total, 1)  # Ensure at least 1
        
        if t.state == "RUNNING":
            ps.running = len([x for x in world.tasks.values() if x.name == t.name and x.state == "RUNNING"])
        elif t.state == "COMPLETED":
            ps.completed = len([x for x in world.tasks.values() if x.name == t.name and x.state == "COMPLETED"])
            ps.update(t)
        elif t.state == "FAILED":
            ps.failed = len([x for x in world.tasks.values() if x.name == t.name and x.state == "FAILED"])
        elif t.state == "CACHED":
            ps.cached = len([x for x in world.tasks.values() if x.name == t.name and x.state == "CACHED"])
    
    # Mark slow tasks
    for t in run:
        ps = world.proc_stats.get(t.name)
        if ps and ps.avg_duration:
            runtime = t.runtime()
            if runtime and runtime > ps.avg_duration * 1.5:
                t.slow_flag = True
    
    # Filter
    if world.filt:
        run = [t for t in run if world.filt.search(t.label() or "")]
        queued = [t for t in queued if world.filt.search(t.label() or "")]
    
    run.sort(key=lambda x: x.first_ts)
    queued.sort(key=lambda x: x.first_ts)
    
    # Update cumulative counters
    world.cum_done = max(world.cum_done, done)
    world.cum_cached = max(world.cum_cached, cache)
    world.cum_failed = max(world.cum_failed, fail)
    world.cum_killed = max(world.cum_killed, killed)
    
    # Progress: use cumulative counters so pruned tasks still count.
    # cum_done/cum_cached are monotonically non-decreasing; cum_seen tracks
    # every distinct task hash seen in the log.
    completed_like = world.cum_done + world.cum_cached
    snapshot_total = len(run) + len(queued) + done + fail + cache + killed
    total = max(world.cum_seen, snapshot_total)
    progress_pct = int((completed_like * 100 / max(1, total))) if total else 0
    
    return run, queued, dict(
        done=done, failed=fail, cached=cache, killed=killed,
        running=len(run), queued=len(queued),
        cum_seen=world.cum_seen, cum_done=world.cum_done,
        cum_cached=world.cum_cached, cum_failed=world.cum_failed,
        cum_killed=world.cum_killed,
        total=total, progress_pct=progress_pct
    )

# ═══════════════════════════════════════════════════════════════
# ─── Export Functions ───────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

def export_json(world: World, path: str):
    """Export snapshot to JSON"""
    run, que, C = classify(world)
    
    data = {
        "timestamp": datetime.now().isoformat(),
        "run_name": world.meta.run_name,
        "session": world.meta.session,
        "executor": world.meta.executor,
        "work_dir": world.meta.work_root,
        "cpu_pct": world.cpu_pct,
        "mem_pct": world.mem_pct,
        "load_1m": world.load_1,
        "done": C["done"],
        "failed": C["failed"],
        "cached": C["cached"],
        "killed": C["killed"],
        "running": C["running"],
        "queued": C["queued"],
        "tasks_running": [
            {
                "id": t.id,
                "label": t.label(),
                "workdir": t.workdir,
                "pid": t.pid,
                "cpu_pct": t.metrics.cpu_pct,
                "rss_mb": t.metrics.rss_mb,
                "runtime": t.runtime(),
                "slow": t.slow_flag,
            }
            for t in run
        ],
        "process_stats": {
            name: {
                "total": ps.total,
                "completed": ps.completed,
                "failed": ps.failed,
                "cached": ps.cached,
                "avg_duration": ps.avg_duration,
                "median_duration": ps.median_duration,
                "p95_duration": ps.p95_duration,
            }
            for name, ps in world.proc_stats.items()
        },
        "alerts": [
            {
                "level": a.level,
                "message": a.message,
                "age": a.age_str(),
            }
            for a in list(world.alerts)[-20:]
        ],
        "recent_errors": list(world.recent_errors),
    }
    
    with open(path, "w") as f:
        json.dump(data, f, indent=2)

def export_html(world: World, path: str):
    """Export dashboard to HTML"""
    run, que, C = classify(world)
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>nfmon Dashboard - {world.meta.run_name}</title>
    <style>
        body {{ font-family: monospace; background: #1e1e1e; color: #d4d4d4; padding: 20px; }}
        h1 {{ color: #4ec9b0; }}
        h2 {{ color: #569cd6; margin-top: 30px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
        th {{ background: #2d2d30; padding: 8px; text-align: left; border: 1px solid #3e3e42; }}
        td {{ padding: 8px; border: 1px solid #3e3e42; }}
        .running {{ color: #4ec9b0; }}
        .failed {{ color: #f48771; }}
        .done {{ color: #608b4e; }}
        .slow {{ background: #3e2723; }}
        .alert-error {{ color: #f48771; }}
        .alert-warning {{ color: #dcdcaa; }}
    </style>
</head>
<body>
    <h1>nfmon Dashboard</h1>
    <p><strong>Pipeline:</strong> {world.meta.run_name}</p>
    <p><strong>Session:</strong> {world.meta.session}</p>
    <p><strong>Executor:</strong> {world.meta.executor}</p>
    <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    
    <h2>Pipeline Progress</h2>
    <table>
        <tr>
            <th>Total</th><th>Running</th><th>Completed</th><th>Failed</th><th>Cached</th>
        </tr>
        <tr>
            <td>{C['running'] + C['queued'] + C['done'] + C['failed'] + C['cached']}</td>
            <td class="running">{C['running']}</td>
            <td class="done">{C['done']}</td>
            <td class="failed">{C['failed']}</td>
            <td>{C['cached']}</td>
        </tr>
    </table>
    
    <h2>Running Tasks</h2>
    <table>
        <tr>
            <th>Process</th><th>Tag</th><th>CPU%</th><th>Memory</th><th>Runtime</th><th>Status</th>
        </tr>
"""
    
    for t in run[:50]:
        cpu = f"{t.metrics.cpu_pct:.1f}%" if t.metrics.cpu_pct else "-"
        mem = f"{t.metrics.rss_mb:.0f}MB" if t.metrics.rss_mb else "-"
        runtime = format_duration(t.runtime()) if t.runtime() else "-"
        slow_class = ' class="slow"' if t.slow_flag else ''
        
        html += f"""        <tr{slow_class}>
            <td>{t.name}</td>
            <td>{t.tag or '-'}</td>
            <td>{cpu}</td>
            <td>{mem}</td>
            <td>{runtime}</td>
            <td>{t.stage or '-'}</td>
        </tr>
"""
    
    html += """    </table>
    
    <h2>Process Statistics</h2>
    <table>
        <tr>
            <th>Process</th><th>Total</th><th>Done</th><th>Failed</th><th>Avg Duration</th>
        </tr>
"""
    
    for name, ps in sorted(world.proc_stats.items()):
        avg = format_duration(ps.avg_duration) if ps.avg_duration else "-"
        html += f"""        <tr>
            <td>{name}</td>
            <td>{ps.total}</td>
            <td>{ps.completed}</td>
            <td>{ps.failed}</td>
            <td>{avg}</td>
        </tr>
"""
    
    html += """    </table>
    
    <h2>Recent Alerts</h2>
    <ul>
"""
    
    for alert in list(world.alerts)[-20:]:
        level_class = f"alert-{alert.level.lower()}"
        html += f'        <li class="{level_class}">[{alert.level}] {alert.message} ({alert.age_str()})</li>\n'
    
    html += """    </ul>
</body>
</html>
"""
    
    with open(path, "w") as f:
        f.write(html)

def export_csv(world: World, path: str):
    """Export process stats to CSV"""
    import csv
    
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Process", "Total", "Completed", "Failed", "Cached", "Avg_Duration_Sec", "Median_Duration_Sec", "P95_Duration_Sec"])
        
        for name, ps in sorted(world.proc_stats.items()):
            writer.writerow([
                name,
                ps.total,
                ps.completed,
                ps.failed,
                ps.cached,
                f"{ps.avg_duration:.2f}" if ps.avg_duration else "",
                f"{ps.median_duration:.2f}" if ps.median_duration else "",
                f"{ps.p95_duration:.2f}" if ps.p95_duration else "",
            ])

# ═══════════════════════════════════════════════════════════════
# ─── Curses TUI (Simple Fallback) ──────────────────────────────
# ═══════════════════════════════════════════════════════════════

class TUI:
    """Curses-based TUI"""
    def __init__(self, stdscr, w: World, a):
        self.s = stdscr
        self.w = w
        self.a = a
        self.sel = 0
        self.show_trends = False
        
        curses.start_color()
        curses.use_default_colors()
        curses.init_pair(1, curses.COLOR_CYAN, -1)
        curses.init_pair(2, curses.COLOR_GREEN, -1)
        curses.init_pair(3, curses.COLOR_RED, -1)
        curses.init_pair(4, curses.COLOR_YELLOW, -1)
        curses.init_pair(5, curses.COLOR_BLUE, -1)
    
    def draw(self):
        self.s.erase()
        maxy, maxx = self.s.getmaxyx()
        now = int(time.time() - self.w.start_ts)
        run, que, C = classify(self.w)
        
        # Compute cores in use
        cores_in_use = sum([max(0.0, t.metrics.cpus or 0) for t in run])
        if self.w.ncpu <= 0:
            try:
                if sys.platform.startswith("darwin"):
                    self.w.ncpu = int(subprocess.check_output(["sysctl", "-n", "hw.ncpu"], text=True))
                else:
                    self.w.ncpu = int(subprocess.check_output(["bash", "-lc", "nproc 2>/dev/null || getconf _NPROCESSORS_ONLN"], text=True))
            except Exception:
                pass
        
        # Progress: snapshot-based (completed+cached vs all known tasks)
        pct = C.get("progress_pct", 0)
        total = C.get("total", C["running"] + C["queued"] + C["done"] + C["failed"] + C["cached"] + C["killed"])

        # Header (compact)
        rn = (self.w.meta.run_name or "").strip()
        run_label = rn[:18] if rn and rn != "?" else "nf-monitor"
        cores_str = f" C{int(cores_in_use)}/{self.w.ncpu}" if self.w.ncpu > 0 else ""
        mode = getattr(self.w, "mode", "")
        mode_str = f" [{mode}]" if mode else ""
        self._add(0, 0, f"{run_label} {time.strftime('%H:%M:%S')}{mode_str}", curses.color_pair(1) | curses.A_BOLD)
        self._add(1, 0, f"[{sec2hms(now)}] CPU:{self.w.cpu_pct}% MEM:{self.w.mem_pct}% LOAD:{self.w.load_1}{cores_str} | {C['running']} run {C['done']} done {C['failed']} fail")
        eta = predict_completion_eta(self.w)
        if eta:
            self._add(1, maxx - 12, format_eta(eta), curses.color_pair(4))
        # Progress bar
        barw = max(10, maxx - 12)
        fill = int(barw * pct / 100)
        self._add(2, 0, ("━" * fill) + ("─" * (barw - fill)) + f"  {pct:3d}%")
        # Layout: header 3 rows, rest for list + details
        header_h = 4
        avail_h = max(0, maxy - header_h)
        
        list_h = int(avail_h * 0.6)
        details_h = avail_h - list_h
        
        # Ensure minimums
        if details_h < 10:
            list_h = avail_h
            details_h = 0
        
        # Task list
        row = header_h
        self._add(row, 0, "Running:", curses.color_pair(5))
        row += 1
        
        list_content_h = max(0, list_h - 1)
        
        if run:
            start = max(0, min(self.sel, max(0, len(run) - list_content_h)))
            end = min(len(run), start + list_content_h)
            
            for i in range(start, end):
                t = run[i]
                age = sec2hms(int(t.age()))
                cpu = f"{t.metrics.cpu_pct:.0f}%" if t.metrics.cpu_pct is not None else "--"
                rss = f"{t.metrics.rss_mb:.0f}MB" if t.metrics.rss_mb is not None else "--"
                mark = "▶" if i == self.sel else " "
                label = t.label() or t.id[:10]
                
                # Sparkline if trends mode
                spark = ""
                if self.show_trends and t.metrics.cpu_hist:
                    spark = " " + t.metrics.cpu_sparkline()
                
                # Slow indicator
                slow = "🐌" if t.slow_flag else ""
                
                line = f"{mark} {truncate(label, 30):30} CPU:{cpu:>4} RSS:{rss:>7} {age:>8}{spark}{' ' + slow if slow else ''}"
                attr = curses.color_pair(4) if t.slow_flag else 0
                self._add(row, 0, truncate(line, maxx - 1), attr)
                row += 1
            
            for _ in range(end - start, list_content_h):
                self._add(row, 0, "")
                row += 1
        else:
            self._add(row, 0, "none")
            row += 1
        
        # ─── Draw Details Pane ───
        if details_h > 0:
            target_t = run[self.sel] if run and self.sel < len(run) else None
            # Reserve last row for keys hint
            details_h = maxy - (header_h + list_h) - 1
            if details_h > 0:
                self.draw_details_pane(target_t, header_h + list_h, maxy - 1, maxx)

        # Keys hint at bottom (compact but informative)
        if maxy > 0:
            hint = " j/k=up-down  t=trends  e=export  q=quit "
            self._add(maxy - 1, 0, hint[:maxx], curses.A_DIM)

        self.s.refresh()

    def draw_details_pane(self, t: Optional[Task], start_y: int, maxy: int, maxx: int):
        """Draw process details pane at bottom"""
        h = maxy - start_y
        w = maxx
        
        if h < 5 or w < 20:
            return
            
        try:
            win = curses.newwin(h, w, start_y, 0)
            # win.box() # Optional: box the whole pane
        except curses.error:
            return
            
        if not t:
            win.addstr(0, 0, "Details: (none)", curses.color_pair(5))
            win.overwrite(self.s)
            return

        win.addstr(0, 0, truncate(f"Details: {t.label()}", w - 1), curses.color_pair(5))
        if t.workdir and os.path.isdir(t.workdir):
            payload, logs, outs = insight(t.workdir, max(5, h - 3))
        else:
            payload, logs, outs = "Workdir not found", [], []

        inner_w = w - 2
        # Line 1: cmd | line 2+: logs (skip Outputs label to save space)
        win.addstr(1, 0, truncate(payload, inner_w), curses.A_DIM)
        avail_log_h = h - 2
        
        for i, ln in enumerate(logs):
            if i >= avail_log_h:
                break
            try:
                win.addstr(2 + i, 0, truncate(ln, inner_w), 0)
            except curses.error:
                pass
                
        win.overwrite(self.s)
    
    def _add(self, y, x, s, attr=0):
        try:
            self.s.addstr(y, x, s[:max(0, self.s.getmaxyx()[1] - x - 1)], attr)
        except curses.error:
            pass
    
    def keys(self):
        self.s.nodelay(True)
        ch = self.s.getch()
        
        if ch in (ord('q'), 27):
            raise SystemExit
        if ch in (curses.KEY_DOWN, ord('j')):
            self.sel += 1
        if ch in (curses.KEY_UP, ord('k')):
            self.sel = max(0, self.sel - 1)
        if ch == ord('t'):
            self.show_trends = not self.show_trends
        if ch == ord('e'):
            # Quick export
            ts = datetime.now().strftime('%Y%m%d_%H%M%S')
            export_json(self.w, f"nfmon_snapshot_{ts}.json")
            export_html(self.w, f"nfmon_dashboard_{ts}.html")

# ═══════════════════════════════════════════════════════════════
# ─── Oneshot Mode ───────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

def oneshot(w: World):
    """Print snapshot to terminal"""
    run, que, C = classify(w)
    total = C.get("total", sum([C[k] for k in ["running", "queued", "done", "failed", "cached", "killed"]]))
    pct = C.get("progress_pct", 0)

    container = getattr(w.meta, "container_engine", "") or "none"
    print(f"nf-monitor  {time.strftime('%H:%M:%S')}  (run:{w.meta.run_name}  session:{w.meta.session}  exec:{w.meta.executor}  container:{container})")
    print(f" Root: {os.getcwd()}")
    print(f" Uptime {sec2hms(int(time.time() - w.start_ts))} | CPU:{w.cpu_pct}% | MEM:{w.mem_pct}% | LOAD:{w.load_1}")
    print()
    bw = 80
    fill = int(bw * pct / 100)
    print("┌─ Pipeline ─┐")
    print(("━" * fill) + ("─" * (bw - fill)) + f"  {pct:3d}%")
    print(f"  total:{total}  run:{C['running']}  done:{C['done']}  fail:{C['failed']}  cache:{C['cached']}")
    
    # ETA
    eta = predict_completion_eta(w)
    if eta:
        print(f"  Estimated completion: {format_eta(eta)}")
    
    print("┌─ Running ─┐")
    for t in run[:20]:
        cpu = f"{int(t.metrics.cpu_pct)}%" if t.metrics.cpu_pct is not None else "--"
        rss = f"{int(t.metrics.rss_mb)}MB" if t.metrics.rss_mb is not None else "--"
        slow = " [slow]" if t.slow_flag else ""
        print(f"  {t.label():54} CPU:{cpu:>4} RSS:{rss:>6}{slow}")
    if not run:
        print("  none")
    print("└─┘")
    
    # Process stats
    if w.proc_stats:
        print("\n┌─ Process Statistics ─┐")
        for name, ps in sorted(w.proc_stats.items()):
            avg = format_duration(ps.avg_duration) if ps.avg_duration else "-"
            print(f"  {name:30}  done:{ps.completed:>3}  fail:{ps.failed:>2}  avg:{avg}")
        print("└─┘")
    
    # Alerts
    print("\n┌─ Recent Alerts ─┐")
    if w.alerts:
        for alert in list(w.alerts)[-5:]:
            print(f"  [{alert.level}] {alert.message}")
    else:
        print("  ✓ none")
    print("└─┘")

# ═══════════════════════════════════════════════════════════════
# ─── Main Entry Point ───────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════

def main():
    a = parse_args()
    w = World()
    w.tail_n = max(1, a.tail)
    w.retain_sec = max(0, a.retain_sec)
    w.max_tasks = max(100, a.max_tasks)
    
    if a.filter:
        w.filt = re.compile(a.filter)
    
    log = abspath(a.log) if a.log else ""
    # Trace: explicit path, or auto-detect (results/trace or results_*/trace, most recent)
    if a.trace:
        w.trace_path = abspath(a.trace)
    else:
        candidates = []
        try:
            for p in ["results/trace/trace.txt"] + [
                os.path.join(d, "trace", "trace.txt") for d in (os.listdir(".") or [])
                if d.startswith("results") and os.path.isdir(d)
            ]:
                if os.path.isfile(p):
                    candidates.append((p, os.path.getmtime(p)))
        except OSError:
            pass
        w.trace_path = abspath(max(candidates, key=lambda x: x[1])[0]) if candidates else None
    
    if os.path.isfile(log):
        w.meta.work_root = guess_work_root(log, a.work)
        bootstrap(w, log)
    else:
        # Try to infer work root first
        if a.work:
            w.meta.work_root = abspath(a.work)
        else:
            w.meta.work_root = guess_work_root(".nextflow.log", a.work)
        # Search upwards from work root (or cwd) for a .nextflow.log
        alt_log = find_nextflow_log(w.meta.work_root or os.getcwd())
        if alt_log and os.path.isfile(alt_log):
            log = alt_log
            w.meta.work_root = guess_work_root(log, a.work)
            bootstrap(w, log)
        else:
            print(f"Warning: log not found, running in trace-first mode (log={a.log})", file=sys.stderr)
    
    w.cpu_pct, w.mem_pct, w.load_1 = sys_metrics()

    # Determine operating mode based on available metadata
    if os.path.isfile(log) and w.trace_path:
        w.mode = "full"
    elif os.path.isfile(log):
        w.mode = "log_only"
    elif w.trace_path:
        w.mode = "trace_only"
    else:
        w.mode = "fs_only"
    
    # Resolve hash mode
    if a.resolve_hash:
        key = a.resolve_hash.strip()
        wanted = key.split("/")[-1]
        
        if os.path.isfile(log):
            bootstrap(w, log)
        
        for row in tail_trace(w):
            apply_trace_row(w, row, time.time())
        
        hit = None
        for t in w.tasks.values():
            if t.id and t.id.split("/")[-1].startswith(wanted):
                hit = t
                break
        
        if not hit:
            for d in running_dirs(find_work_roots(w.meta.work_root or "work")):
                tid = "/".join(d.rstrip("/").split("/")[-2:])
                if tid.split("/")[-1].startswith(wanted):
                    name, tag = label_from_dir(d)
                    hit = Task(id=tid, name=name, tag=tag, state="RUNNING", workdir=d)
                    break
        
        if hit:
            out = dict(
                id=hit.id,
                short=short_hash(hit.id),
                process=hit.name or "?",
                tag=hit.tag or "",
                workdir=hit.workdir or "",
                state=hit.state or "",
            )
            print(json.dumps(out, indent=2))
        else:
            print(json.dumps({"error": "not found", "query": key}, indent=2))
        return
    
    # Oneshot mode
    if a.oneshot:
        if a.json:
            export_json(w, a.json)
        if a.html:
            export_html(w, a.html)
        if a.csv:
            export_csv(w, a.csv)
        oneshot(w)
        return
    
    def _sig(_s, _f):
        raise SystemExit
    
    signal.signal(signal.SIGINT, _sig)
    
    def loop(stdscr):
        stdscr.nodelay(True)
        stdscr.keypad(True)
        ui = TUI(stdscr, w, a)
        next_tick = 0.0
        
        while True:
            now = time.time()
            if now >= next_tick:
                # Nextflow log stream
                if log and os.path.isfile(log):
                    for ln in tail_lines(w, log):
                        update_meta(w, ln)
                        parse_line(w, ln, time.time())
                        if any(k in ln for k in ("ERROR", "FAILED", "Exception", "Caused by:")) and "collect-file" not in ln:
                            w.recent_errors.append(ln.strip())
                
                # Trace stream
                for row in tail_trace(w):
                    apply_trace_row(w, row, time.time())
                    if not w.meta.work_root and row.get("workdir"):
                        try:
                            wd = row.get("workdir").strip()
                            if wd:
                                parts = wd.rstrip("/").split("/")
                                if len(parts) >= 3:
                                    w.meta.work_root = "/".join(parts[:-2])
                        except Exception:
                            pass
                
                w.cpu_pct, w.mem_pct, w.load_1 = sys_metrics()
                
                if a.json:
                    tmp = a.json + ".tmp"
                    export_json(w, tmp)
                    os.replace(tmp, a.json)
                
                next_tick = now + max(0.2, a.refresh)
            
            ui.draw()
            ui.keys()
            time.sleep(0.02)
    
    # Try Rich UI first (unless --simple)
    if not a.simple:
        try:
            from rich.live import Live
            from rich.table import Table
            from rich.panel import Panel
            from rich.layout import Layout
            from rich.align import Align
            from rich.text import Text
            from rich.rule import Rule
            import rich.box
            import hashlib

            def _task_label(t: Task) -> str:
                """Generate human-readable task label, avoiding cryptic hashes"""
                # Try to extract from filesystem if name is missing or looks like hash
                if not t.name or looks_like_hash(t.name):
                    ensure_task_label_from_fs(w, t)
                
                # If we have a proper name, use it
                if t.name and not looks_like_hash(t.name):
                    return f"{t.name} ({t.tag})" if t.tag else t.name
                
                # Fallback: try to extract from workdir path
                if t.workdir:
                    try:
                        # workdir is typically .../work/xx/yyyyyyyy
                        # Try to find .command.env for process name
                        name, tag = label_from_dir(t.workdir)
                        if name and not looks_like_hash(name):
                            t.name = name
                            t.tag = tag
                            return f"{name} ({tag})" if tag else name
                    except Exception:
                        pass
                
                # Last resort: show short hash with prefix
                lab = t.label() or (t.id or "?")
                if looks_like_hash(lab):
                    return f"task_{lab[:12]}"
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

            def _vstack(*items):
                """Stack renderables without Panel borders.
                First item (typically a Rule) spans full width; subsequent items
                get a 1-line top gap and 2-char left/right gutters."""
                outer = Table(show_header=False, show_edge=False,
                              box=None, padding=0, expand=True)
                outer.add_column(ratio=1)
                inner = Table(show_header=False, show_edge=False,
                              box=None, padding=(0, 2), expand=True)
                inner.add_column(ratio=1)
                for i, item in enumerate(items):
                    if i == 0:
                        outer.add_row(item)       # Rule: full width, no padding
                    else:
                        inner.add_row(Text(""))   # 1-line breathing gap
                        inner.add_row(item)       # content with 2-char side gutters
                outer.add_row(inner)
                return outer

            def mini_bar(pct: float, width: int = 6,
                         full_style: str = "green", empty_style: str = "dim") -> Text:
                """Return a Rich Text mini fill bar e.g. ▓▓▓▓░░"""
                filled = min(width, max(0, int(round(pct / 100.0 * width))))
                t = Text()
                t.append("▓" * filled,          style=full_style)
                t.append("░" * (width - filled), style=empty_style)
                return t

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
                        r,_,_ = select.select([sys.stdin], [], [], 0.05)
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
                            if ch == 'h':
                                new_val = not filter_state.get('show_help', False)
                                filter_state['show_help'] = new_val
                                if new_val:
                                    filter_state['help_opened_at'] = time.time()
                    termios.tcsetattr(fd, termios.TCSADRAIN, old)
                except Exception:
                    pass

            def render():
                import shutil
                try:
                    term_height = shutil.get_terminal_size().lines
                except Exception:
                    term_height = 40
                run,que,C=classify(w)   # single call — reused throughout render()

                # Auto-close help after 60 s of being open
                if filter_state.get('show_help') and \
                        time.time() - filter_state.get('help_opened_at', 0) > 60:
                    filter_state['show_help'] = False

                # ── Header (borderless) ─────────────────────────────────────
                show_help = filter_state.get('show_help', False)
                sort_str  = filter_state.get('sort', 'default')
                rn = (w.meta.run_name or "?")
                ss = (w.meta.session or "?")
                ss_short = ss[:8] if (ss != "?" and len(ss) > 8) else ss
                container = getattr(w.meta, "container_engine", "") or "none"

                hdr = Text()
                # Line 1: branded title + identity + uptime
                hdr.append(" ◈ nf-monitor ", style="bold white on dark_cyan")
                hdr.append("  run:", style="dim")
                hdr.append(f" {rn if rn!='?' else 'n/a'}", style="bold cyan")
                hdr.append("  sess:", style="dim")
                hdr.append(f" {ss_short if ss!='?' else 'n/a'}", style="")
                hdr.append(f"  {w.meta.executor}  {container}", style="dim")
                hdr.append(f"  ⏱ {sec2hms(int(time.time() - w.start_ts))}\n", style="dim")

                # Line 2: system resources
                cores_in_use = sum([max(0.0, t.metrics.cpus) for t in run if t.metrics.cpus is not None])
                if w.ncpu <= 0:
                    try:
                        if sys.platform.startswith("darwin"):
                            w.ncpu=int(subprocess.check_output(["sysctl","-n","hw.ncpu"],text=True))
                        else:
                            w.ncpu=int(subprocess.check_output(["bash","-lc","nproc 2>/dev/null || getconf _NPROCESSORS_ONLN"],text=True))
                    except Exception:
                        w.ncpu = 0
                if cores_in_use <= 0 and w.ncpu > 0:
                    approx = sum(max(0.0, t.metrics.cpu_pct/100.0) for t in run if t.metrics.cpu_pct is not None)
                    cores_in_use = min(float(w.ncpu), approx)
                cores_txt = f"  cores {int(round(cores_in_use))}/{w.ncpu}" if w.ncpu > 0 else ""
                hdr.append(f"sys  cpu {w.cpu_pct}%  mem {w.mem_pct}%  load {w.load_1}{cores_txt}\n\n", style="dim")

                # Line 3: progress bar + task badges (blank line above provided by \n\n)
                pct = C.get("progress_pct", 0)
                try:
                    term_w = shutil.get_terminal_size().columns
                except Exception:
                    term_w = 80
                bar_w = max(4, term_w - 27)
                fill = int(bar_w * pct / 100)
                hdr.append("━" * fill,           style="green")
                hdr.append("━" * (bar_w - fill), style="dim")
                hdr.append(f"  {pct}%  ", style="bold")
                n_run  = C.get("running", 0);  n_done = C.get("done", 0) + C.get("cached", 0)
                n_fail = C.get("failed",  0);  n_que  = C.get("queued", 0)
                if n_run:  hdr.append(f"▶{n_run} ",  style="bold cyan")
                if n_done: hdr.append(f"✓{n_done} ", style="green")
                if n_fail: hdr.append(f"✗{n_fail} ", style="bold red")
                if n_que:  hdr.append(f"⏳{n_que}",  style="dim")
                hdr.append("\n")  # trailing blank line before the Rule separator

                # _vstack puts Rule first (full-width) then hdr with side gutters
                # but header is special: hdr comes BEFORE the Rule, so we build it manually
                _hdr_outer = Table(show_header=False, show_edge=False,
                                   box=None, padding=0, expand=True)
                _hdr_outer.add_column(ratio=1)
                _hdr_inner = Table(show_header=False, show_edge=False,
                                   box=None, padding=(0, 2), expand=True)
                _hdr_inner.add_column(ratio=1)
                _hdr_inner.add_row(hdr)
                _hdr_outer.add_row(_hdr_inner)
                _hdr_outer.add_row(Rule(style="cyan"))
                header_grp = _hdr_outer

                # running table
                rt_title = "Running (all logs)" if getattr(a, 'all_logs', False) else "Running (j/k to focus)"
                rt = Table(title=None, expand=True, padding=(0, 1), show_edge=False,
                           box=rich.box.SIMPLE, header_style="bold dim", show_lines=False)
                rt.add_column("",       width=1,       no_wrap=True)   # focus/load dot
                rt.add_column("Module", max_width=30, no_wrap=True)   # capped so Sample gets room
                rt.add_column("Sample", ratio=1,      no_wrap=True)   # fills remaining flex space
                rt.add_column("CPU",    justify="right", width=16) # "▓▓▓▓░░ 45%/4c"
                rt.add_column("MEM",    justify="right", width=12) # "▓▓░░ 2.1G"
                rt.add_column("Age",    justify="right", width=7)
                rt.add_column("▲",      width=10)                  # CPU sparkline
                # ensure labels from FS if missing
                for t in run:
                    ensure_task_label_from_fs(w, t)
                run_filtered = [t for t in run if (not filter_state["proc"] or (t.name or "?") == filter_state["proc"]) ]
                # apply sort option
                sort_mode = filter_state.get('sort') or 'default'
                if sort_mode == 'cpu':
                    run_filtered.sort(key=lambda t: (-(t.metrics.cpu_pct if t.metrics.cpu_pct is not None else -1), -(t.metrics.rss_mb if t.metrics.rss_mb is not None else -1), -(time.time()-t.first_ts)))
                elif sort_mode == 'rss':
                    run_filtered.sort(key=lambda t: (-(t.metrics.rss_mb if t.metrics.rss_mb is not None else -1), -(t.metrics.cpu_pct if t.metrics.cpu_pct is not None else -1), -(time.time()-t.first_ts)))
                elif sort_mode == 'age':
                    run_filtered.sort(key=lambda t: (-(time.time()-t.first_ts)))
                max_rows = max(6, min(18, (term_height - 12) // 2))
                for idx, t in enumerate(run_filtered[:max_rows]):
                    # Ensure we have proper labels from filesystem
                    ensure_task_label_from_fs(w, t)
                    
                    # CPU — mini bar + number
                    cpu_pct = t.metrics.cpu_pct
                    cpus_n  = t.metrics.cpus
                    cpu_style = ("green" if (cpu_pct or 0) < 50 else ("yellow" if (cpu_pct or 0) < 80 else "red")) if cpu_pct is not None else "dim"
                    if cpu_pct is not None:
                        cpu_col = mini_bar(cpu_pct, 6, cpu_style)
                        suffix = f" {cpu_pct:.0f}%"
                        if cpus_n is not None: suffix += f"/{int(cpus_n)}c"
                        cpu_col.append(suffix, style=cpu_style)
                    elif cpus_n is not None:
                        cpu_col = Text(f"░░░░░░ ──/{int(cpus_n)}c", style="dim")
                    else:
                        cpu_col = Text("░░░░░░  ──", style="dim")

                    # MEM — mini bar + number (scale: 32 GB = 100%)
                    _rss = t.metrics.rss_mb
                    rss_style = ("green" if (_rss or 0) < 2048 else ("yellow" if (_rss or 0) < 8192 else "red")) if _rss is not None else "dim"
                    if _rss is not None:
                        rss_pct = min(100.0, _rss / 32768.0 * 100)
                        mem_col = mini_bar(rss_pct, 4, rss_style)
                        if _rss >= 1024:
                            mem_col.append(f" {_rss/1024:.1f}G", style=rss_style)
                        else:
                            mem_col.append(f" {_rss:.0f}M", style=rss_style)
                    else:
                        mem_col = Text("░░░░ ──", style="dim")

                    # Module name: strip workflow prefix ("TrackTx:align_reads_to_genome" → "align_reads_to_genome")
                    if not t.name or looks_like_hash(t.name):
                        if t.workdir:
                            name_temp, tag_temp = label_from_dir(t.workdir)
                            if name_temp and not looks_like_hash(name_temp):
                                t.name = name_temp; t.tag = tag_temp
                    module_name = t.name if (t.name and not looks_like_hash(t.name)) else f"task_{short_hash(t.id or '?')}"
                    if ':' in module_name:
                        module_name = module_name.split(':', 1)[-1]   # drop "TrackTx:" prefix
                    sample_tag = (t.tag or "─")
                    is_focused = (not getattr(a, 'all_logs', False)) and idx == focused["idx"]
                    if is_focused:
                        dot = Text("►", style="bold cyan")
                    elif cpu_pct is not None:
                        if cpu_pct > 80:   dot = Text("●", style="bold red")
                        elif cpu_pct > 50: dot = Text("●", style="yellow")
                        elif cpu_pct > 10: dot = Text("●", style="green")
                        else:              dot = Text("●", style="dim")
                    else:
                        dot = Text("·", style="dim")
                    lab = Text(module_name, style=_name_style(module_name))
                    if is_focused:
                        lab.stylize("bold underline")

                    hist  = list(t.metrics.cpu_hist)
                    spark = ""
                    if hist:
                        try:
                            spark = sparkline_pct(hist[-20:])
                        except Exception:
                            spark = ""
                    rt.add_row(dot, lab, sample_tag, cpu_col, mem_col,
                               sec2hms(int(time.time()-t.first_ts)), Text(spark, style="cyan"))

                # Wrap running table in a borderless group with a Rule title.
                # _vstack adds: 1-line gap + 2-char side gutters after the Rule.
                # Extra Text("") row adds a blank line below the table before the next section.
                _rt_title = f"▶ Running [{len(run_filtered)}]" if run_filtered else "▶ Running  (none)"
                if not getattr(a, 'all_logs', False): _rt_title += "  j/k navigate"
                rt_group = _vstack(Rule(title=_rt_title, style="cyan", align="left"), rt, Text(""))

                # live log tail panel(s) - dynamically sized based on available space
                log_panel = _vstack(Rule(title="Log", style="blue", align="left"),
                                    Text("No running tasks", style="yellow"))

                # help — comprehensive reference panel
                _secs_open = int(time.time() - filter_state.get('help_opened_at', time.time()))
                _secs_left = max(0, 60 - _secs_open)
                help_txt = Text()

                help_txt.append("KEYBOARD SHORTCUTS\n", style="bold cyan")
                help_txt.append("  j / ↓      ", style="bold yellow")
                help_txt.append("Move focus down to next running task\n")
                help_txt.append("  k / ↑      ", style="bold yellow")
                help_txt.append("Move focus up to previous running task\n")
                help_txt.append("  f          ", style="bold yellow")
                help_txt.append("Cycle process filter: All → each module type in turn → All\n")
                help_txt.append("  s          ", style="bold yellow")
                help_txt.append("Cycle sort order: default → by CPU% → by memory (RSS) → by age\n")
                help_txt.append("  a          ", style="bold yellow")
                help_txt.append("Toggle all-logs: focused task only ↔ all running tasks stacked\n")
                help_txt.append("  h          ", style="bold yellow")
                help_txt.append("Toggle this help panel (auto-closes after 60 s)\n")
                help_txt.append("  q / Esc    ", style="bold yellow")
                help_txt.append("Quit nf-monitor\n")
                help_txt.append("\n")

                help_txt.append("COLUMN REFERENCE\n", style="bold cyan")
                help_txt.append("  ●/►   ", style="bold yellow")
                help_txt.append("Focus & load indicator — see legend below\n")
                help_txt.append("  Module", style="bold yellow")
                help_txt.append("  Nextflow process name (workflow prefix stripped)\n")
                help_txt.append("  Sample", style="bold yellow")
                help_txt.append("  Task tag — usually the sample ID or input file name\n")
                help_txt.append("  CPU   ", style="bold yellow")
                help_txt.append("  ▓▓▓░░░ bar (6-wide) + actual%/allocated cores\n")
                help_txt.append("         ", style="bold yellow")
                help_txt.append("  bar is proportional to CPU%; ── means no data yet\n", style="dim")
                help_txt.append("  MEM   ", style="bold yellow")
                help_txt.append("  ▓▓░░ bar (4-wide, 32 GB=full) + current RSS\n")
                help_txt.append("         ", style="bold yellow")
                help_txt.append("  shows G (gigabytes) or M (megabytes)\n", style="dim")
                help_txt.append("  Age   ", style="bold yellow")
                help_txt.append("  Wall-clock time since the task was submitted\n")
                help_txt.append("  ▲     ", style="bold yellow")
                help_txt.append("  CPU% sparkline — last 20 samples, ▁=low █=high\n")
                help_txt.append("\n")

                help_txt.append("INDICATOR LEGEND\n", style="bold cyan")
                help_txt.append("  ► ", style="bold cyan")
                help_txt.append("cyan      focused row (log panel tracks this task)\n")
                help_txt.append("  ● ", style="green")
                help_txt.append("green     CPU 10–50%   moderate load\n")
                help_txt.append("  ● ", style="yellow")
                help_txt.append("yellow    CPU 50–80%   high load\n")
                help_txt.append("  ● ", style="bold red")
                help_txt.append("red       CPU > 80%    saturated / bottleneck\n")
                help_txt.append("  · ", style="dim")
                help_txt.append("dim       CPU unknown or task idle\n")
                help_txt.append("\n")

                help_txt.append("PROGRESS BAR BADGES\n", style="bold cyan")
                help_txt.append("  ▶N  ", style="bold cyan")
                help_txt.append("N tasks currently running\n")
                help_txt.append("  ✓N  ", style="green")
                help_txt.append("N tasks completed successfully (includes cached)\n")
                help_txt.append("  ✗N  ", style="bold red")
                help_txt.append("N tasks failed or errored\n")
                help_txt.append("  ⏳N ", style="dim")
                help_txt.append("N tasks queued (submitted but not yet started)\n")
                help_txt.append("\n")

                help_txt.append("LOG PANEL\n", style="bold cyan")
                help_txt.append("The log panel shows live tail of the focused task's ")
                help_txt.append(".command.log", style="bold")
                help_txt.append(" file.\n")
                help_txt.append("Summary line at top: process name, sample tag, elapsed, CPU, RSS.\n")
                help_txt.append("A ▶ line shows the most recent meaningful output (insight).\n")
                help_txt.append("Press ")
                help_txt.append("a", style="bold yellow")
                help_txt.append(" to switch to all-logs mode — stacks every running task.\n")
                help_txt.append("\n")

                help_txt.append("HOW CPU & MEMORY ARE MEASURED\n", style="bold cyan")
                help_txt.append("CPU and RSS are summed across the ")
                help_txt.append("full process subtree", style="bold")
                help_txt.append(" (the bash wrapper\n")
                help_txt.append("+ all child processes it spawned). This gives accurate numbers\n")
                help_txt.append("even for tools like STAR or BWA that fork many threads.\n")
                help_txt.append("CPU allocation (/Nc) is read from ")
                help_txt.append("OMP_NUM_THREADS", style="bold")
                help_txt.append(" in .command.sh\n")
                help_txt.append("(written by Nextflow's beforeScript) or from completed trace rows.\n")
                help_txt.append("\n")

                help_txt.append("TIPS\n", style="bold cyan")
                help_txt.append("  • ", style="dim")
                help_txt.append("Use --refresh 1.0 to reduce overhead on large cluster runs\n", style="dim")
                help_txt.append("  • ", style="dim")
                help_txt.append("Use --filter 'align' to pre-filter tasks by regex from startup\n", style="dim")
                help_txt.append("  • ", style="dim")
                help_txt.append("On scratch=true (NFS), PID lookup scans /proc — Linux only\n", style="dim")
                help_txt.append("  • ", style="dim")
                help_txt.append("--retain-sec 0 clears completed tasks immediately from the list\n", style="dim")
                help_txt.append("  • ", style="dim")
                help_txt.append("Run with --oneshot to print a summary snapshot and exit\n", style="dim")

                _help_rule_title = f"Help  ·  h to close  ·  auto-closes in {_secs_left}s"
                help_panel = _vstack(Rule(title=_help_rule_title, style="bright_black"), help_txt)

                # ── Footer bar: key hints + active state + error status ────────
                foot = Text()
                for key, desc in [("j/k","nav"), ("f","filter"), ("s","sort"),
                                   ("a","all-logs"), ("h","help"), ("q","quit")]:
                    foot.append(f" {key} ", style="bold yellow")
                    foot.append(f"{desc}  ", style="dim")
                if filter_state["proc"]: foot.append(f"│ proc:{filter_state['proc']}  ", style="yellow")
                if sort_str != 'default': foot.append(f"│ sort:{sort_str}  ", style="yellow")
                if getattr(a, 'all_logs', False): foot.append("│ ALL-LOGS  ", style="yellow")
                if show_help: foot.append("│ HELP  ", style="bright_black")
                foot.append("  │  ", style="dim")
                if w.recent_errors:
                    last_err = list(w.recent_errors)[-1].strip()
                    foot.append(f"✗ {last_err[:90]}", style="bold red")
                else:
                    foot.append("✓ no errors", style="green")
                footer_grp = _vstack(Rule(style="bright_black"), foot)

                # Estimate log lines from terminal (header=6, footer=3, log Rule+gap=2)
                available_log_lines = max(5, term_height - 11)
                
                if run_filtered:
                    if getattr(a, 'all_logs', False):
                        # Calculate how many processes and lines per process
                        num_processes = len(run_filtered)
                        
                        # Lines needed per process: title(1) + payload(1) + separator(1) = 3 base lines
                        # Plus variable log lines
                        base_lines_per_process = 3
                        
                        # How many processes can we show?
                        # Start by assuming we want at least 3 log lines per process
                        min_log_lines_per_process = 3
                        lines_per_process = base_lines_per_process + min_log_lines_per_process
                        
                        max_processes = max(1, available_log_lines // lines_per_process)
                        processes_to_show = min(num_processes, max_processes)
                        
                        # Now calculate actual lines per process
                        if processes_to_show > 0:
                            lines_per_process = available_log_lines // processes_to_show
                            log_lines_per_process = max(3, lines_per_process - base_lines_per_process)
                        else:
                            log_lines_per_process = 10
                        
                        # stack multiple logs into one panel
                        text = Text()
                        for i in range(processes_to_show):
                            t = run_filtered[i]
                            payload, tail, outs = insight(t.workdir, log_lines_per_process)
                            text.append(f"{_task_label(t)} (id:{short_hash(t.id or '')})\n", style="bold")
                            if payload:
                                text.append(f"▶ {payload}\n", style="cyan")
                            # Show only the calculated number of lines
                            for ln in tail[:log_lines_per_process]:
                                text.append(ln+"\n")
                            if outs:
                                text.append("outputs: "+", ".join([f"{n} ({sz})" for n,sz in outs[:2]]))
                            if i < processes_to_show - 1:
                                text.append("\n", style="dim")
                        log_panel = _vstack(
                            Rule(title=f"Log ({processes_to_show}/{num_processes})  ·  a to toggle", style="blue", align="left"),
                            text)
                    else:
                        # single focused log - use ALL available lines
                        if focused["idx"] >= len(run_filtered):
                            focused["idx"] = len(run_filtered)-1
                        ft = run_filtered[max(0, focused["idx"])]
                        focused["id"] = ft.id
                        
                        # Use all available lines for single process
                        # Account for: summary header(2) + payload(1) + outputs(1) = 4 overhead lines
                        single_log_lines = max(5, available_log_lines - 4)

                        payload, tail, outs = insight(ft.workdir, single_log_lines)
                        text = Text()
                        # Task summary mini-header inside log panel
                        elapsed_s = sec2hms(int(time.time() - ft.first_ts))
                        _cpu_s = ""
                        if ft.metrics.cpu_pct is not None and ft.metrics.cpus is not None:
                            _cpu_s = f"{ft.metrics.cpu_pct:.0f}%/{int(ft.metrics.cpus)}c"
                        elif ft.metrics.cpu_pct is not None:
                            _cpu_s = f"{ft.metrics.cpu_pct:.0f}%"
                        _rss_s = f"{ft.metrics.rss_mb:.0f} MB" if ft.metrics.rss_mb is not None else ""
                        _meta = "  ".join(x for x in [elapsed_s, _cpu_s, _rss_s] if x)
                        _mn = (ft.name or "?").split(":", 1)[-1]
                        text.append(f"● {_mn}", style="bold cyan")
                        if ft.tag: text.append(f"  ({ft.tag})", style="")
                        if _meta: text.append(f"  {_meta}", style="dim")
                        text.append("\n")
                        text.append("─" * 60 + "\n", style="dim")
                        if payload:
                            text.append(f"▶ {payload}\n", style="cyan")
                        for ln in tail[:single_log_lines]:
                            text.append(ln+"\n")
                        if outs:
                            text.append("outputs: "+", ".join([f"{n} ({sz})" for n,sz in outs]))
                        log_panel = _vstack(
                            Rule(title=f"Log — {_mn}  ({focused['idx']+1}/{len(run_filtered)})", style="blue", align="left"),
                            text)
                
                layout = Layout()
                layout.split_column(
                    Layout(header_grp, size=6),
                    Layout(name="body"),
                    Layout(footer_grp, size=3),
                )
                layout["body"].split_row(Layout(name="right", ratio=1))
                # Dynamic height for running table: only as tall as it needs to be.
                # SIMPLE box: 1 header line + 1 separator line + N data rows.
                # Plus 1 for the Rule title above = N + 3 total.
                _body_h = term_height - 6 - 2  # header_grp(6) + footer_grp(2)
                _rt_rows = len(run_filtered[:max_rows])
                # _vstack layout: Rule(1) + gap(1) + table header+sep(2) + rows(N) + trailing blank(1) = N+5
                _rt_height = max(5, min(_rt_rows + 5, _body_h - 8))  # always leave ≥8 for log
                # Stacked: running tasks (top, dynamic) | log (middle) | optional help (bottom)
                if show_help:
                    layout["right"].split_column(
                        Layout(rt_group,   size=_rt_height),
                        Layout(log_panel,  ratio=2, minimum_size=5),
                        Layout(help_panel, ratio=1, minimum_size=10),
                    )
                else:
                    layout["right"].split_column(
                        Layout(rt_group,  size=_rt_height),
                        Layout(log_panel, ratio=1, minimum_size=5),
                    )
                return layout

            import threading
            kt = threading.Thread(target=_key_thread, daemon=True)
            kt.start()
            with Live(render(), refresh_per_second=max(5, int(1.0 / max(0.15, a.refresh)))) as live:
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
                        next_tick = now + max(0.2, a.refresh)
                    time.sleep(0.02)
            stop_keys["stop"] = True
            return
        except Exception as e:
            if not a.simple and isinstance(e, ImportError):
                print("Tip: pip install rich  for enhanced UI", file=sys.stderr)

    # Curses fallback
    curses.wrapper(loop)

if __name__ == "__main__":
    main()

