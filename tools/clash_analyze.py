# -*- coding: utf-8 -*-
"""
PPI docking pose quality analyzer
=================================
对 LightDock/LKlight 输出的 pose PDB 做量化 clash / 界面接触分析：
  1. 自动区分受体/配体原子（pose 前 N_rec 个原子 = 受体，其余 = 配体；
     N_rec 取引擎 setup 写出的 lightdock_rec.pdb 原子数，若存在）
  2. 界面接触：跨分子 原子间距 d < 4.5 Å 的对（残基对/原子对）
  3. vdW clash：用 AMBER 风格 vdW 半径，判断 d < 0.75*(ri+rj) 轻度 / <0.6 严重，
     并给出最大穿透深度
  4. 输出最严重的 8 对 clash，方便到 PyMOL/VMD 定位

用法: python3 clash_analyze.py pose.pdb [lightdock_rec.pdb]
"""
import sys, math
from collections import defaultdict

# AMBER/常用 vdW 半径（Å），元素符号 -> 半径
VDW = {'C':1.70,'N':1.55,'O':1.52,'H':1.10,'P':1.80,'S':1.80,
       'F':1.47,'CL':1.75,'BR':1.85,'I':1.98,'MG':1.18,'ZN':1.10,
       'CA':1.70,'NA':1.36,'K':2.02,'FE':1.40,'MN':1.40,'SE':1.90}

def elem(atom_name):
    """从 PDB 原子名推断元素（右对齐原子名: 第13列起4字符，去空格取首字母+可能的小写）"""
    a = atom_name.strip()
    if not a:
        return 'C'
    # 处理如 "CA"、"FE"、"CL" 双字母元素
    two = {'CL','BR','MG','ZN','CA','NA','FE','MN','SE','NI','CU','CO'}
    if a[:2].upper() in two:
        return a[:2].upper()
    if a[0].isdigit():
        a = a[1:]
    if not a:
        return 'C'
    return a[0].upper()

def parse_pdb(path):
    """返回 [(serial, name, resname, chain, resnum, x, y, z, element), ...]"""
    atoms = []
    for ln in open(path):
        if ln.startswith(('ATOM','HETATM')):
            try:
                x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
            except ValueError:
                continue
            atoms.append((ln[12:16].strip(), ln[17:20].strip(), ln[21:22],
                          int(ln[22:26]) if ln[22:26].strip() else 0,
                          x, y, z))
    return atoms

def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    pose_path = sys.argv[1]
    rec_ref = sys.argv[2] if len(sys.argv) > 2 else None

    pose = parse_pdb(pose_path)
    if not pose:
        print('pose 无原子'); return

    # 受体原子数：用 lightdock_rec.pdb 参照（引擎实际写入的坐标系文件）
    n_rec = None
    if rec_ref:
        try:
            n_rec = len(parse_pdb(rec_ref))
        except Exception:
            n_rec = None
    if n_rec is None or n_rec <= 0 or n_rec > len(pose):
        n_rec = None

    if n_rec is not None:
        rec = pose[:n_rec]; lig = pose[n_rec:]
        mode = f'按参照文件切分: 受体 {len(rec)} / 配体 {len(lig)}'
    else:
        # fallback：按链切分（多链时受体取最大链集合）
        chains = defaultdict(list)
        for a in pose:
            chains[a[2]].append(a)
        if len(chains) >= 2:
            sorted_chains = sorted(chains.items(), key=lambda kv: -len(kv[1]))
            rec = sorted_chains[0][1]; lig = sum((v for k,v in sorted_chains[1:]), [])
            mode = f'按链切分: 受体(链{sorted_chains[0][0]}) {len(rec)} / 配体 {len(lig)}'
        else:
            print('无法区分受体/配体'); return
    print(f'=== {pose_path} ===')
    print(f'归属: {mode}')
    print(f'受体链/残基: {sorted(set(a[2] for a in rec))} 残基数 {len(set((a[2],a[3]) for a in rec))}')
    print(f'配体链/残基: {sorted(set(a[2] for a in lig))} 残基数 {len(set((a[2],a[3]) for a in lig))}')

    # 原子坐标
    def coords(atoms):
        return [(a[4],a[5],a[6]) for a in atoms]
    rc = coords(rec); lc = coords(lig)

    contacts = []   # (d, ri, rj) 原子对索引
    clash_light = []
    clash_severe = []
    max_penetration = 0.0
    worst = []
    # 双重循环跨分子（可能较大，用朴素 O(N^2)，对 pose 数千原子 OK）
    for i in range(len(rc)):
        xi, yi, zi = rc[i]
        ei = elem(rec[i][0]); rvi = VDW.get(ei, 1.7)
        for j in range(len(lc)):
            dx = xi - lc[j][0]; dy = yi - lc[j][1]; dz = zi - lc[j][2]
            d2 = dx*dx + dy*dy + dz*dz
            if d2 >= 4.5*4.5:
                continue
            d = math.sqrt(d2)
            if d < 1e-6:
                continue
            contacts.append((d, i, j))
            ej = elem(lig[j][0]); rvj = VDW.get(ej, 1.7)
            sum_r = rvi + rvj
            if d < 0.60 * sum_r:
                clash_severe.append((d, i, j))
                pen = sum_r - d
                if pen > max_penetration:
                    max_penetration = pen
                worst.append((pen, d, i, j, rec[i], lig[j], ei, ej))
            elif d < 0.75 * sum_r:
                clash_light.append((d, i, j))
                pen = sum_r - d
                if pen > max_penetration:
                    max_penetration = pen
                worst.append((pen, d, i, j, rec[i], lig[j], ei, ej))

    # 界面残基对
    res_pairs = set()
    for d, i, j in contacts:
        res_pairs.add((rec[i][2], rec[i][3], lig[j][2], lig[j][3]))

    n_heavy_i = len(rec)
    n_heavy_j = len(lig)

    print(f'\n── 界面接触 ──')
    print(f'  接触原子对 (<4.5Å): {len(contacts)}')
    print(f'  接触残基对: {len(res_pairs)}')
    # 区分：含氢/重原子
    heavy_contacts = sum(1 for d,i,j in contacts
                         if elem(rec[i][0])!='H' and elem(lig[j][0])!='H')
    print(f'  其中 重原子-重原子 接触对: {heavy_contacts}')

    print(f'\n── vdW clash 检查 ──')
    tot_pairs = len(contacts)
    print(f'  轻度过近 (d < 0.75×vdW和): {len(clash_light)} 对')
    print(f'  严重穿透 (d < 0.60×vdW和): {len(clash_severe)} 对')
    print(f'  最大穿透深度: {max_penetration:.2f} Å'
          f' (vdW 和超出量, >0.5Å 需警惕, >1.0Å 基本可判为无效 pose)')
    if tot_pairs:
        print(f'  严重 clash 占界面接触比例: {100.0*len(clash_severe)/max(tot_pairs,1):.1f}%')

    print(f'\n── 最严重 clash 前 8 对 (到 PyMOL 用 resi 定位) ──')
    worst.sort(reverse=True)
    for pen, d, i, j, ra, la, ei, ej in worst[:8]:
        rn = f'{ra[1]}{ra[3]}{ra[2]}'  # e.g. GLY123A
        ln = f'{la[1]}{la[3]}{la[2]}'
        print(f'  受体 {rn:<12} {ra[0]:<5}({ei}) × 配体 {ln:<12} {la[0]:<5}({ej})'
              f' | d={d:.2f}Å  vdW和={VDW.get(ei,1.7)+VDW.get(ej,1.7):.2f}Å  穿透={pen:.2f}Å')

if __name__ == '__main__':
    main()
