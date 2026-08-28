#!/bin/bash
# 场景测试: 4 个生物医学相关复合物 × 代表性评分函数
set -u
ROOT=/Users/luoxiaowen/Desktop/LKDock/byi/LKlight
BIN=$ROOT/target/release/LKlight
OUT=/tmp/lktest
mkdir -p $OUT/scenarios

echo "scenario,method,g=200,steps=100,time_s,best_gso_score,status" > $OUT/scenarios/results.csv

run_case () {
  local name=$1 rec=$2 lig=$3; shift 3
  local methods="$@"
  for m in $methods; do
    d=$OUT/scenarios/${name}_${m}; rm -rf $d; mkdir -p $d; cd $d
    # 蛋白-DNA 案例: 只在 dna/ddna/dfire2 上有意义, 但也测通用势
    if ! $BIN setup "$rec" "$lig" -s 1 -g 200 > setup.log 2>&1; then
      echo "$name,$m,-,-,-,SETUP_FAIL" >> $OUT/scenarios/results.csv; continue
    fi
    start=$(python3 -c 'import time;print(time.time())')
    if $BIN run setup.json initial_positions_0.dat 100 $m > run.log 2>&1; then
      end=$(python3 -c 'import time;print(time.time())')
      ts=$(python3 -c "print(f'{$end-$start:.2f}')")
      best=$(grep -oE 'gso_score of glowworm [0-9]+ = [-0-9.eE]+' run.log | awk '{print $NF}' | sort -g | head -1)
      echo "$name,$m,200,100,$ts,${best:-NA},OK" >> $OUT/scenarios/results.csv
    else
      end=$(python3 -c 'import time;print(time.time())')
      ts=$(python3 -c "print(f'{$end-$start:.2f}')")
      echo "$name,$m,200,100,$ts,-,RUN_FAIL" >> $OUT/scenarios/results.csv
    fi
  done
}

# 场景1: p53 肿瘤抑制核心域 - DNA (癌症位点研究, 蛋白-DNA)
run_case p53DNA /tmp/lktest/cases/1DIZ_p53.pdb /tmp/lktest/cases/1DIZ_DNA.pdb dna ddna dfire2 cpydock dfire
# 场景2: 抗体-溶菌酶 (Ab-Ag 经典, 1VFB)
run_case AbLyso /tmp/lktest/cases/1VFB_antibody.pdb /tmp/lktest/cases/1VFB_lysozyme.pdb pydock cpydock dfire sd vdw
# 场景3: 抗体-HIV 肽抗原 (1DQJ)
run_case AbHIVpep /tmp/lktest/cases/1DQJ_antibody.pdb /tmp/lktest/cases/1DQJ_peptide.pdb pydock cpydock dfire sd vdw
# 场景4: SARS-CoV-2 RBD - hACE2 (病毒-宿主, 6M0J)
run_case RBD_ACE2 /tmp/lktest/cases/6M0J_ACE2.pdb /tmp/lktest/cases/6M0J_RBD.pdb pydock cpydock dfire sd vdw

echo "== scenarios done =="
cat $OUT/scenarios/results.csv
