#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Robust decoy/randbox/realbox Vina sweep runner.

Why this exists:
  The earlier `xargs -P 8 bash -c '{}'` path worked for *execution* but hit
  two problems under the Windows/Cygwin sandbox:
    1. A grep option-parsing bug (`grep -oE '--out ...'`) broke the skip
       logic, so already-complete tasks were re-run.
    2. Repeated Cygwin `fork()` failures (`dofork: child died ... 0xC000026B`)
       could stall the whole sweep.

This runner avoids both:
  * It spawns vina.exe DIRECTLY as a native Windows process via
    subprocess(list, shell=False) -> no Cygwin fork, no cmd.exe $VAR issue.
  * Skip logic is done in Python by checking the --out file size, so
    completed tasks are never re-run. The sweep is therefore resumable:
    just re-run the script; it skips anything already >= MIN_BYTES.

Usage:
  python run_decoy_py.py            # full sweep (resumes automatically)
  python run_decoy_py.py --cap 12   # smoke test: only first 12 pending tasks
  python run_decoy_py.py --workers 6 --cpu 4
"""
import os
import sys
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

# ---- paths (derived from this script's own location) --------------------
DC = os.path.dirname(os.path.abspath(__file__))                 # .../decoy_control
ROOT = os.path.dirname(os.path.dirname(DC))                     # .../方案二_...
REC_PQ = os.path.join(ROOT, "06_docking", "hpc_bundle_autodl", "06_docking", "rec_pdbqt")
LIG_PQ = os.path.join(ROOT, "06_docking", "hpc_bundle_autodl", "06_docking", "lig_pdbqt")
DEC_PQ = os.path.join(DC, "lig_decoy_pdbqt")
OUTDIR = DC

ENV = {
    "$REC_PQ": REC_PQ,
    "$LIG_PQ": LIG_PQ,
    "$DEC_PQ": DEC_PQ,
    "$OUTDIR": OUTDIR,
}

CMDS_FILE = os.path.join(DC, "vina_cmds_abs.txt")
MIN_BYTES = 5000          # a complete Vina --out pdbqt is > 5 KB
VINA_CPU = 4              # cap vina internal threads to limit total load
WORKERS = 8
TIMEOUT = 600             # seconds per docking job


def expand_token(tok: str) -> str:
    for k, v in ENV.items():
        if k in tok:
            tok = tok.replace(k, v)
    return tok


def parse_out_path(tokens):
    """Return the value of --out (expanded)."""
    for i, t in enumerate(tokens):
        if t == "--out" and i + 1 < len(tokens):
            return expand_token(tokens[i + 1])
        if t.startswith("--out="):
            return expand_token(t.split("=", 1)[1])
    return None


def build_cmd(line: str):
    raw = line.strip()
    if not raw:
        return None
    tokens = shlex.split(raw)            # safe split, no shell expansion
    tokens = [expand_token(t) for t in tokens]
    # cap vina internal threads (does not change the score, only speed)
    has_cpu = any(t == "--cpu" for t in tokens)
    if not has_cpu:
        tokens += ["--cpu", str(VINA_CPU)]
    # quieten Vina's verbose banner to keep things light
    has_verb = any(t == "--verbosity" for t in tokens)
    if not has_verb:
        tokens += ["--verbosity", "1"]
    return tokens


def worker(idx_cmd):
    idx, cmd = idx_cmd
    tokens = build_cmd(cmd)
    if tokens is None:
        return (idx, "skip-empty", None)
    out_path = parse_out_path(tokens)
    # skip if already complete
    if out_path and os.path.isfile(out_path) and os.path.getsize(out_path) >= MIN_BYTES:
        return (idx, "skip-done", out_path)
    try:
        proc = subprocess.run(
            tokens,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            timeout=TIMEOUT,
        )
        rc = proc.returncode
    except subprocess.TimeoutExpired:
        return (idx, "timeout", out_path)
    except Exception as e:  # pragma: no cover
        return (idx, f"error:{e}", out_path)
    # verify output
    if out_path and os.path.isfile(out_path) and os.path.getsize(out_path) >= MIN_BYTES:
        return (idx, "ok", out_path)
    # failed: surface stderr tail for debugging
    err = (proc.stderr or b"").decode("utf-8", "replace")[-400:]
    return (idx, f"fail-rc{rc}:{err}", out_path)


def main():
    global WORKERS
    cap = None
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--cap" and i + 1 < len(args):
            cap = int(args[i + 1]); i += 2
        elif a == "--workers" and i + 1 < len(args):
            WORKERS = int(args[i + 1]); i += 2
        elif a == "--cpu" and i + 1 < len(args):
            globals()["VINA_CPU"] = int(args[i + 1]); i += 2
        else:
            i += 1

    with open(CMDS_FILE, "r", encoding="utf-8") as f:
        lines = [l for l in f.read().splitlines() if l.strip()]

    # pre-filter: which are still pending (avoids spinning threads on done ones)
    pending = []
    already = 0
    for idx, line in enumerate(lines):
        tokens = build_cmd(line)
        out_path = parse_out_path(tokens) if tokens else None
        if out_path and os.path.isfile(out_path) and os.path.getsize(out_path) >= MIN_BYTES:
            already += 1
        else:
            pending.append((idx, line))

    if cap is not None:
        pending = pending[:cap]

    total = len(lines)
    print(f"[run_decoy_py] total_cmds={total} already_done={already} "
          f"pending_this_run={len(pending)} workers={WORKERS} cpu={VINA_CPU}", flush=True)

    if not pending:
        print("[run_decoy_py] NOTHING PENDING — sweep complete.", flush=True)
        return

    done = 0
    ok = 0
    bad = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(worker, pc) for pc in pending]
        for fut in as_completed(futs):
            idx, status, outp = fut.result()
            done += 1
            if status == "ok":
                ok += 1
            elif status.startswith("skip"):
                pass
            else:
                bad += 1
                print(f"  [#{idx}] {status}  out={outp}", flush=True)
            if done % 25 == 0 or done == len(pending):
                print(f"  progress {done}/{len(pending)}  ok={ok} bad={bad}", flush=True)

    print(f"[run_decoy_py] FINISHED this run: ok={ok} bad={bad} "
          f"(already_done={already} of {total})", flush=True)


if __name__ == "__main__":
    main()
