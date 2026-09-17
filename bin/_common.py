#!/usr/bin/env python3
"""
Shared logging / __main__ / atomic-write helpers for TrackTx's bin/*.py
scripts.

Before this module existed, three scripts (calculate_pol_metrics.py,
combine_reports.py, compare_pol_metrics.py) each independently defined the
byte-for-byte identical log_info/log_warning/log_error/log_progress
functions and __main__ try/except block, and five other scripts each did
something different again (a quiet-flag timestamp logger, a one-line stderr
print, or no logging helper at all -- bare print()). This module gives every
script the same underlying implementation to call instead of re-defining it.

Import note: Nextflow puts `bin/` on PATH (so these scripts run directly as
commands) but NOT on PYTHONPATH, so a script can't just `import _common`
from an arbitrary working directory. Every script using this module starts
with:

    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from _common import make_logger, run_main, AtomicFileWriter

`os.path.abspath(__file__)` resolves to the script's real location
(bin/<script>.py) regardless of the Nextflow task's working directory, so
this works the same whether the script is invoked via PATH lookup or a
direct path.
"""
import os
import sys
import tempfile


# =============================================================================
# TRACKTX_THEME — canonical color tokens for every TrackTx-generated report
# =============================================================================
# Single source of truth for the palette shared by all five report/UI
# surfaces: the per-sample report (render_sample_report.py), the cohort
# report (combine_reports.py), the run landing page (inline Python in
# modules/15_combine_reports_into_cohort.nf), TrackTx_config_generator.html,
# and the per-sample plots page (also render_sample_report.py).
#
# Before this existed, each of those hand-rolled its own `:root{...}` colour
# block independently and drifted to four different palettes (confirmed by
# diffing them) even though three shared the same CSS variable NAMES
# (--bg/--fg/--accent/...). This module does not inject CSS into any of
# them (they're independently generated -- three different Python files plus
# one static HTML file, no shared build/import step) -- it's the value each
# of those `:root` blocks is kept manually in sync with. When changing the
# palette, edit here first, then update the four `:root` blocks (grep for
# "TRACKTX_THEME" in render_sample_report.py, combine_reports.py,
# modules/15_combine_reports_into_cohort.nf, and TrackTx_config_generator.html).
TRACKTX_THEME = {
    "dark": {
        "bg": "#0f1117", "surface": "#1a1d27", "surface2": "#232636",
        "border": "#2e3248", "accent": "#6c8ef5", "accent2": "#4ecdc4",
        "text": "#e8eaf0", "muted": "#7b82a0",
    },
    "light": {
        "bg": "#ffffff", "surface": "#f8fafc", "surface2": "#eef1f6",
        "border": "#e2e5ec", "accent": "#4a6cf0", "accent2": "#2bb3a8",
        "text": "#1a1d23", "muted": "#6b7280",
    },
}


def make_logger(prefix: str):
    """
    Returns (log_info, log_warning, log_error, log_progress) bound to
    `prefix`, e.g.:
        log_info, log_warning, log_error, log_progress = make_logger("POL_CALC")
    Output format: "[PREFIX] LEVEL | message" -- the convention already
    used by calculate_pol_metrics.py/combine_reports.py/compare_pol_metrics.py
    before this module existed; kept identical so downstream log-reading
    (Nextflow task logs, humans debugging a run) doesn't need to change.
    """
    def log_info(message: str):
        print(f"[{prefix}] INFO | {message}", flush=True)

    def log_warning(message: str):
        print(f"[{prefix}] WARNING | {message}", flush=True)

    def log_error(message: str):
        print(f"[{prefix}] ERROR | {message}", file=sys.stderr, flush=True)

    def log_progress(section: str, current: int, total: int):
        if total > 0:
            percent = (current / total) * 100
            print(f"[{prefix}] {section} | Progress: {current}/{total} ({percent:.1f}%)", flush=True)

    return log_info, log_warning, log_error, log_progress


def run_main(main_func, log_error) -> int:
    """
    Standard __main__ entry point, replacing the identical try/except block
    previously duplicated across three scripts (and adopted by the rest
    during this cleanup). Use as:

        if __name__ == "__main__":
            sys.exit(run_main(main, log_error))

    Handles KeyboardInterrupt (exit 130) and any other exception (prints
    traceback, exit 1) uniformly; SystemExit (e.g. from argparse or an
    explicit sys.exit inside main_func) passes through untouched so an
    intentional exit code isn't overwritten.
    """
    try:
        return main_func()
    except KeyboardInterrupt:
        log_error("Interrupted by user")
        return 130
    except SystemExit:
        raise
    except Exception as e:
        log_error(f"Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


class AtomicFileWriter:
    """
    Context manager for atomically writing a file: writes go to a temp file
    in the same directory, which is moved into place with os.replace() only
    on a clean exit. On any exception, the partial temp file is removed and
    the exception propagates -- the real output path is never left holding
    truncated content from a killed or crashed process.

    Matches the mktemp+mv pattern already used on the shell-script side of
    this pipeline (e.g. modules 01/04's genome/index caching); this brings
    the same guarantee to bin/*.py scripts' primary outputs, which
    previously wrote directly via open(path, 'w').

    Usage:
        with AtomicFileWriter(path) as f:
            f.write(...)
    """

    def __init__(self, path, mode: str = "w", encoding: str = "utf-8"):
        self.path = path
        self.mode = mode
        self.encoding = None if "b" in mode else encoding
        self._tmp_path = None
        self._fh = None

    def __enter__(self):
        d = os.path.dirname(os.path.abspath(self.path)) or "."
        fd, self._tmp_path = tempfile.mkstemp(dir=d, prefix=".tmp_", suffix=".partial")
        # mkstemp creates the temp file at mode 0600 (owner-only), which
        # os.replace() would otherwise carry over to the final path -- unlike
        # a plain open(path, "w"), which gets the umask-derived default
        # (typically 644). Match that default so switching to atomic writes
        # doesn't silently make pipeline outputs unreadable by other users.
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(fd, 0o666 & ~umask)
        self._fh = os.fdopen(fd, self.mode) if self.encoding is None else os.fdopen(fd, self.mode, encoding=self.encoding)
        return self._fh

    def __exit__(self, exc_type, exc, tb):
        self._fh.close()
        if exc_type is None:
            os.replace(self._tmp_path, self.path)
        else:
            try:
                os.remove(self._tmp_path)
            except OSError:
                pass
        return False


def _selftest():
    """Runs on every import via __main__ below -- exercises AtomicFileWriter's
    both paths (commit on success, discard on failure), since that's the one
    piece of real logic in this module worth a regression check."""
    import shutil

    d = tempfile.mkdtemp(prefix="_common_selftest_")
    try:
        target = os.path.join(d, "out.txt")

        with AtomicFileWriter(target) as f:
            f.write("hello")
        assert os.path.exists(target)
        with open(target) as f:
            assert f.read() == "hello"
        assert not any(n.startswith(".tmp_") for n in os.listdir(d))

        try:
            with AtomicFileWriter(target) as f:
                f.write("should not land")
                raise ValueError("boom")
        except ValueError:
            pass
        with open(target) as f:
            assert f.read() == "hello"  # unchanged -- failed write never committed
        assert not any(n.startswith(".tmp_") for n in os.listdir(d))  # temp cleaned up

        assert run_main(lambda: 0, lambda m: None) == 0
        # run_main prints a traceback for unexpected exceptions by design;
        # redirect stderr so a passing self-test doesn't look like a crash.
        import contextlib
        import io
        with contextlib.redirect_stderr(io.StringIO()):
            assert run_main(lambda: (_ for _ in ()).throw(ValueError("x")), lambda m: None) == 1
    finally:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    _selftest()
    print("_common.py self-test passed")
