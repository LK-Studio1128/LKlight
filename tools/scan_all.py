# -*- coding: utf-8 -*-
"""扫描全部 pose：评分 vs 深clash，找 Pareto 最优解"""
import math
src = open('/tmp/clash_analysis/clash_analyze.py').read()
ns = {}
exec(src.split("def main")[0], ns)
parse_pdb, VDW, elem = ns['parse_pdb'], ns['VDW'], ns['elem']

rec_ref = parse_pdb('lightdock_receptor.pdb')
n_rec = len(rec_ref)

def analyze(pose_file):
    pose = parse_pdb(pose_file)
    rec, lig = pose[:n_rec], pose[n_rec:]
    # 配体网格 cell=2.4 (4.5/2)
    CELL = 2.4
    grid = {}
    for j, a in enumerate(lig):
        key = (int(a[4]//CELL), int(a[5]//CELL), int(a[6]//CELL))
        grid.setdefault(key, []).append(j)
    sev = 0; maxpen = 0.0; contacts = 0
    for i, a in enumerate(rec):
        xi, yi, zi = a[4], a[5], a[6]
        ri = VDW.get(elem(a[0]), 1.7)
        cx, cy, cz = int(xi//CELL), int(yi//CELL), int(zi//CELL)
        for gx in (cx-1, cx, cx+1):
            for gy in (cy-1, cy, cy+1):
                for gz in (cz-1, cz, cz+1):
                    for j in grid.get((gx,gy,gz), ()):
                        b = lig[j]
                        dx = xi-b[4]; dy = yi-b[5]; dz = zi-b[6]
                        d2 = dx*dx+dy*dy+dz*dz
                        if d2 >= 20.25: continue
                        d = math.sqrt(d2)
                        if d < 1e-6: continue
                        contacts += 1
                        rj = VDW.get(elem(b[0]), 1.7)
                        if d < 0.6*(ri+rj):
                            sev += 1
                            pen = ri+rj-d
                            if pen > maxpen: maxpen = pen
    return contacts, sev, maxpen

rows = []
for ln in open('ranking.list'):
    if ln.startswith('#'): continue
    p = ln.split()
    if len(p) >= 7:
        rows.append((float(p[3]), p[5]))  # luciferin 评分(p[3]更稳), pdb

# 全部 120
import os
results = []
for score, pdb in rows:
    if not os.path.exists(pdb): continue
    c, sev, mp = analyze(pdb)
    results.append((score, pdb, c, sev, mp))

print(f"{'评分(低=好)':>10} {'深clash':>6} {'最大穿透':>6}  pose")
for score, pdb, c, sev, mp in sorted(results, key=lambda r: -r[0]):
    flag = '  <-- 低clash候选' if sev <= 8 else ''
    print(f"{score:>10.1f} {sev:>6} {mp:>6.2f}  {pdb}{flag}")
