#!/usr/bin/env bash
# Assemble v1.2.0 release assets for LKlight (exact, packaged) / LKlight-grid (bare) / LKlight-GPU (bare).
set -e
BYI=/Users/luoxiaowen/Desktop/LKDock/byi
RL="/Users/luoxiaowen/Desktop/LKDock/LKlight论文/functest/release_v120"
WIN="/Users/luoxiaowen/Desktop/LKDock/LKlight论文/functest/engines_win5090_v120"
LIN=$RL/linux
V=v1.2.0

echo "== 1/3 exact packaged assets =="
rm -rf $RL/exact_out && mkdir -p $RL/exact_out
for p in linux-x86_64 mac-arm64 win-x64; do
  mkdir -p $RL/exact_out/LKlight-$V-$p
  for f in LICENSE CHANGELOG.md NOTICE README.md; do cp $BYI/LKlight/$f $RL/exact_out/LKlight-$V-$p/; done
done
cp $LIN/LKlight-linux-x64-exact   $RL/exact_out/LKlight-$V-linux-x86_64/LKlight
cp $BYI/LKlight/target/release/LKlight $RL/exact_out/LKlight-$V-mac-arm64/LKlight
cp $WIN/LKlight_exact_v1.2.0_win_x64_20260905.exe $RL/exact_out/LKlight-$V-win-x64/LKlight.exe
chmod +x $RL/exact_out/LKlight-$V-linux-x86_64/LKlight $RL/exact_out/LKlight-$V-mac-arm64/LKlight
( cd $RL/exact_out && rm -f ../LKlight-$V-*.tar.gz ../LKlight-$V-*.zip && tar czf ../LKlight-$V-linux-x86_64.tar.gz LKlight-$V-linux-x86_64 && zip -qr ../LKlight-$V-mac-arm64.zip LKlight-$V-mac-arm64 && zip -qr ../LKlight-$V-win-x64.zip LKlight-$V-win-x64 )
echo "exact archives: $(ls -la $RL/LKlight-$V-*.tar.gz $RL/LKlight-$V-*.zip | awk '{print $NF, $5}')"

echo "== 2/3 grid bare assets =="
rm -rf $RL/grid && mkdir -p $RL/grid
cp $LIN/LKlight-linux-x64-grid $RL/grid/LKlight-linux-x64
cp $BYI/LKlight-grid/target/release/LKlight $RL/grid/LKlight-mac-arm64
cp $WIN/LKlight_grid_v1.2.0_win_x64_20260905.exe $RL/grid/LKlight-win64.exe
chmod +x $RL/grid/*
ls -la $RL/grid | awk '{print $NF, $5}'

echo "== 3/3 gpu bare assets =="
rm -rf $RL/gpu && mkdir -p $RL/gpu
cp $LIN/LKlight-linux-x64-cuda $RL/gpu/LKlight-linux-cuda
cp $WIN/LKlight_gpu_v1.2.0_win_x64_sm120_fix_20260905.exe $RL/gpu/LKlight-win64-cuda.exe
chmod +x $RL/gpu/*
ls -la $RL/gpu | awk '{print $NF, $5}'

echo "== md5 summary =="
md5 $RL/LKlight-$V-* $RL/grid/* $RL/gpu/* | sed "s|$RL/||"
echo ASSEMBLE_DONE
