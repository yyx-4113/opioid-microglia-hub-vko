#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-process builder of the Vina resume list (NO forking -> sandbox-safe).

Reads vina_cmds_abs.txt, extracts each task's --out path, expands $OUTDIR,
and keeps only tasks whose output file is missing or < 5000 bytes.
Writes the result to vina_cmds_resume.txt. This replaces the previous
bash `while read` loop that forked 1260 external processes and got killed
by the sandbox's process-rate monitor.
"""
import os

DC = os.path.dirname(os.path.abspath(__file__))
OUTDIR = DC
SRC = os.path.join(DC, "vina_cmds_abs.txt")
DST = os.path.join(DC, "vina_cmds_resume.txt")
MIN_BYTES = 5000


def main():
    with open(SRC, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    pending = []
    done = 0
    for ln in lines:
        if not ln.strip():
            continue
        toks = ln.split()
        out = None
        for i, t in enumerate(toks):
            if t == "--out" and i + 1 < len(toks):
                out = toks[i + 1]
                break
            if t.startswith("--out="):
                out = t.split("=", 1)[1]
                break
        if out is None:
            pending.append(ln)
            continue
        out = out.replace("$OUTDIR", OUTDIR)
        if os.path.isfile(out) and os.path.getsize(out) >= MIN_BYTES:
            done += 1
        else:
            pending.append(ln)
    with open(DST, "w", encoding="utf-8") as f:
        if pending:
            f.write("\n".join(pending) + "\n")
    print(f"build_resume: total={len(lines)} done={done} pending={len(pending)} -> {DST}")


if __name__ == "__main__":
    main()
