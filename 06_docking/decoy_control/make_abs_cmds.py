#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rewrite vina_cmds.txt so each line calls the Vina exe by absolute path
(avoids Git Bash PATH translation issues with D:/ style paths)."""
import io, sys, os

ROOT = r"D:/2026.9/极速交付9月会员日优惠套路/05_多组学+虚拟敲除药物发现/方案二_阿片耐受与痛觉过敏小胶质枢纽基因虚拟敲除"
VINA = ROOT + "/06_docking/vina_bin/vina.exe"
SRC = ROOT + "/06_docking/decoy_control/vina_cmds.txt"
DST = ROOT + "/06_docking/decoy_control/vina_cmds_abs.txt"

with io.open(SRC, "r", encoding="utf-8") as f:
    lines = f.readlines()

out = []
n_fixed = 0
for ln in lines:
    s = ln.rstrip("\n")
    if s.startswith("vina "):
        out.append(VINA + s[4:])
        n_fixed += 1
    else:
        out.append(s)

with io.open(DST, "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")

print("total lines:", len(lines))
print("fixed (vina->abs):", n_fixed)
print("written:", DST)
print("sample:", out[0][:160])
