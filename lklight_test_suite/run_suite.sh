#!/bin/bash
# LKlight 综合测试套件：评分矩阵 + 参数扫描 + 模块功能
set -u
ROOT=/Users/luoxiaowen/Desktop/LKDock/byi/LKlight
BIN=$ROOT/target/release/LKlight
OUT=/tmp/lktest
mkdir -p $OUT/{scores,swarms,params,modules}

METHODS="dfire fastdfire dfire2 dna ddna mj3h pydock cpydock sd vdw pisa sipper tobi"

# ---------- Part 1: 评分函数 × 案例矩阵（score 模块） ----------
CASES=(
  "1azp $ROOT/tests/1azp/1azp_receptor.pdb $ROOT/tests/1azp/1azp_ligand.pdb"
  "2oob $ROOT/tests/2oob/2oob_receptor.pdb $ROOT/tests/2oob/2oob_ligand.pdb"
)
echo "case,method,score,time_ms" > $OUT/scores/matrix.csv
for c in "${CASES[@]}"; do
  name=$(echo $c | cut -d' ' -f1); rec=$(echo $c | cut -d' ' -f2); lig=$(echo $c | cut -d' ' -f3)
  for m in $METHODS; do
    start=$(python3 -c 'import time;print(time.time())')
    score=$($BIN score "$rec" "$lig" $m 2>$OUT/scores/$name.$m.err | grep -oE '[-0-9.eE]+' | head -1)
    end=$(python3 -c 'import time;print(time.time())')
    tms=$(python3 -c "print(f'{($end-$start)*1000:.1f}')")
    [ -z "$score" ] && score="ERR"
    echo "$name,$m,$score,$tms" >> $OUT/scores/matrix.csv
  done
done
echo "== score matrix done =="

# ---------- Part 2: 参数扫描（1AZP pydock + dfire） ----------
# 维度: glowworms g (25/50/100/200/400), steps (10/25/50/100/200)
for m in pydock dfire; do
  for g in 25 50 100 200 400; do
    d=$OUT/params/${m}_g$g; rm -rf $d; mkdir -p $d; cd $d
    $BIN setup $ROOT/tests/1azp/1azp_receptor.pdb $ROOT/tests/1azp/1azp_ligand.pdb -s 1 -g $g >/dev/null 2>&1
    start=$(python3 -c 'import time;print(time.time())')
    $BIN run setup.json initial_positions_0.dat 100 $m > run.log 2>&1
    end=$(python3 -c 'import time;print(time.time())')
    echo "$m,g=$g,steps=100,$(echo "$end-$start"|bc -l)" >> $OUT/params/scan.csv
  done
  for st in 10 25 50 100 200; do
    d=$OUT/params/${m}_s$st; rm -rf $d; mkdir -p $d; cd $d
    $BIN setup $ROOT/tests/1azp/1azp_receptor.pdb $ROOT/tests/1azp/1azp_ligand.pdb -s 1 -g 200 >/dev/null 2>&1
    start=$(python3 -c 'import time;print(time.time())')
    $BIN run setup.json initial_positions_0.dat $st $m > run.log 2>&1
    end=$(python3 -c 'import time;print(time.time())')
    echo "$m,g=200,steps=$st,$(echo "$end-$start"|bc -l)" >> $OUT/params/scan.csv
    # 记录最终最优评分
    best=$(grep -oE 'gso_score.*' run.log | tail -1)
    echo "$m,$st,$best" >> $OUT/params/best_scores.csv
  done
done
echo "== param scan done =="

# ---------- Part 3: 模块功能测试（1AZP, pydock, 1 swarm 200g 100步） ----------
d=$OUT/modules/base; rm -rf $d; mkdir -p $d; cd $d
$BIN setup $ROOT/tests/1azp/1azp_receptor.pdb $ROOT/tests/1azp/1azp_ligand.pdb -s 2 -g 200 --anm > setup.log 2>&1
for sw in 0 1; do
  $BIN run setup.json initial_positions_$sw.dat 100 pydock > run_$sw.log 2>&1
done
$BIN rank 2 100 > rank.log 2>&1
$BIN cluster gso_0.out --cutoff 4.0 > cluster.log 2>&1
$BIN top rank_by_scoring.list 10 > top10.log 2>&1
$BIN gso_to_csv rank_by_scoring.list ranking.csv > csv.log 2>&1
$BIN diameter $ROOT/tests/1azp/1azp_receptor.pdb > diameter.log 2>&1
$BIN map_contacts $ROOT/tests/1azp/1azp_receptor.pdb $ROOT/tests/1azp/1azp_ligand.pdb gso_0.out > contacts.log 2>&1
$BIN trajectory $ROOT/tests/1azp/1azp_receptor.pdb $ROOT/tests/1azp/1azp_ligand.pdb 0 1 > traj.log 2>&1
$BIN generate $ROOT/tests/1azp/1azp_receptor.pdb $ROOT/tests/1azp/1azp_ligand.pdb gso_0.out 5 > gen.log 2>&1
$BIN reference_points $ROOT/tests/1azp/1azp_receptor.pdb --save > refpts.log 2>&1
$BIN pipeline $ROOT/tests/1azp/1azp_receptor.pdb $ROOT/tests/1azp/1azp_ligand.pdb pydock > pipeline.log 2>&1
echo "== modules done =="
ls -la $OUT/modules/base | head -40
