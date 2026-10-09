#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RESUME runner for the 630-task decoy/randbox/realbox Vina 1.2.7 sweep.

- Reads vina_cmds_abs.txt (self-contained vina commands using $REC_PQ/$LIG_PQ/$DEC_PQ/$OUTDIR).
- Skips any task whose --out file already exists AND is >= 5000 bytes (complete Vina output).
- Runs the remaining tasks with 8 parallel workers.
- Logs progress; prints DONE_ALL_DECOY_RUNS and final counts.
This is safe to launch repeatedly: completed outputs are never re-run, so it just
picks up where a stalled/killed run left off.
"""
import os, re, sys, time, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = r"D:/2026.9/极速交付9月会员日优惠套路/05_多组学+虚拟敲除药物发现/方案二_阿片耐受与痛觉过敏小胶质枢纽基因虚拟敲除"
DEC  = ROOT + "/06_docking/decoy_control"
VINA = ROOT + "/06_docking/vina_bin/vina.exe"
CMDS = DEC + "/vina_cmds_abs.txt"

# environment required by the command templates (substituted into the command
# strings below, because the runner executes under cmd.exe which does NOT expand
# bash-style $VAR references)
ENV = {
    "REC_PQ": ROOT + "/06_docking/hpc_bundle_autodl/06_docking/rec_pdbqt",
    "LIG_PQ": ROOT + "/06_docking/hpc_bundle_autodl/06_docking/lig_pdbqt",
    "DEC_PQ": DEC + "/lig_decoy_pdbqt",
    "OUTDIR": DEC,
}

MIN_DONE_BYTES = 5000  # a complete Vina output is ~14-24 KB
WORKERS = 8

def expand(cmd):
    for k, v in ENV.items():
        cmd = cmd.replace("$" + k, v)
    return cmd

def out_path(cmd):
    m = re.search(r"--out\s+(\S+)", cmd)
    return m.group(1) if m else None

def is_done(cmd):
    p = expand(out_path(cmd)) if out_path(cmd) else None
    if not p: return False
    try:
        return os.path.getsize(p) >= MIN_DONE_BYTES
    except OSError:
        return False

def run_one(cmd):
    try:
        r = subprocess.run(expand(cmd), shell=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, timeout=600)
        return r.returncode
    except Exception as e:
        return str(e)

def main():
    with open(CMDS, encoding="utf-8") as f:
        cmds = [ln.strip() for ln in f if ln.strip()]
    done = [c for c in cmds if is_done(c)]
    todo = [c for c in cmds if not is_done(c)]
    print("total=%d already_done=%d to_run=%d  (%s)" % (len(cmds), len(done), len(todo), time.strftime("%H:%M:%S")))
    if not todo:
        print("NOTHING TO RUN — sweep already complete.")
        return
    completed = 0; failed = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(run_one, c): c for c in todo}
        for fut in as_completed(futs):
            rc = fut.result()
            completed += 1
            if rc not in (0, None):
                failed += 1
            if completed % 20 == 0 or completed == len(todo):
                el = time.time() - t0
                print("progress %d/%d  failed=%d  elapsed=%.0fs  (%s)" % (
                    completed, len(todo), failed, el, time.strftime("%H:%M:%S")), flush=True)
    rb = len([1 for c in cmds if "/results_randbox/" in (out_path(c) or "") and is_done(c)])
    dc = len([1 for c in cmds if "/results_decoy/" in (out_path(c) or "") and is_done(c)])
    rl = len([1 for c in cmds if "/results_realbox/" in (out_path(c) or "") and is_done(c)])
    print("DONE_ALL_DECOY_RUNS (%s)  realbox=%d/30 randbox=%d/360 decoy=%d/240 failed=%d" % (
        time.strftime("%H:%M:%S"), rl, rb, dc, failed), flush=True)

if __name__ == "__main__":
    main()
