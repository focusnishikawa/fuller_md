#!/bin/bash
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025, Takeshi Nishikawa
#===========================================================================
#  Test_fuller.sh — フラーレン結晶 NPT-MD 動作検証スクリプト
#
#  使い方:
#    cd fuller_md && src/Test_fuller.sh
#
#  概要:
#    全ての実行モジュールを最小ステップ数で実行し、
#    正常終了するかを検証する。
#    実行はfuller_md/ディレクトリから行うこと
#    (FullereneLib/ への相対パスを利用するため)。
#
#---------------------------------------------------------------------------
#  選択可能なフラーレン (--fullerene= オプション):
#
#  ● C60-76 系列 (FullereneLib/C60-76/):
#    --fullerene=C60       C60 バッキーボール (Ih対称, 60原子) ※デフォルト
#    --fullerene=C70       C70 (D5h対称, 70原子)
#    --fullerene=C72       C72 (D6d対称, 72原子)
#    --fullerene=C74       C74 (D3h対称, 74原子)
#    --fullerene=C76:D2    C76 異性体 (D2対称, 76原子)
#    --fullerene=C76:Td    C76 異性体 (Td対称, 76原子)
#
#  ● C84 系列 (FullereneLib/C84/, 24異性体):
#    --fullerene=C84:1     C84 No.01 (D2対称)
#    --fullerene=C84:2     C84 No.02 (C2対称)
#    --fullerene=C84:3     C84 No.03 (Cs対称)
#    --fullerene=C84:4     C84 No.04 (D2d対称)
#    --fullerene=C84:5     C84 No.05 (D2対称)
#    --fullerene=C84:6     C84 No.06 (C2v対称)
#    --fullerene=C84:7     C84 No.07 (C2v対称)
#    --fullerene=C84:8     C84 No.08 (C2対称)
#    --fullerene=C84:9     C84 No.09 (C2対称)
#    --fullerene=C84:10    C84 No.10 (Cs対称)
#    --fullerene=C84:11    C84 No.11 (C2対称)
#    --fullerene=C84:12    C84 No.12 (C1対称)
#    --fullerene=C84:13    C84 No.13 (C2対称)
#    --fullerene=C84:14    C84 No.14 (Cs対称)
#    --fullerene=C84:15    C84 No.15 (Cs対称)
#    --fullerene=C84:16    C84 No.16 (Cs対称)
#    --fullerene=C84:17    C84 No.17 (C2v対称)
#    --fullerene=C84:18    C84 No.18 (C2v対称)
#    --fullerene=C84:19    C84 No.19 (D3d対称)
#    --fullerene=C84:20    C84 No.20 (Td対称)   ※最も安定な異性体の一つ
#    --fullerene=C84:21    C84 No.21 (D2対称)
#    --fullerene=C84:22    C84 No.22 (D2対称)   ※最も安定な異性体の一つ
#    --fullerene=C84:23    C84 No.23 (D2d対称)  ※最も安定な異性体の一つ
#    --fullerene=C84:24    C84 No.24 (D6h対称)
#
#    対称群を明示する書式: --fullerene=C84:20:Td
#
#---------------------------------------------------------------------------
#  全オプション一覧 (フル版 [3][4][5] 共通):
#
#    --help                  ヘルプ表示
#    --fullerene=<名前>      フラーレン種 (デフォルト: C60)
#    --crystal=<fcc|hcp|bcc> 結晶構造 (デフォルト: fcc)
#    --cell=<nc>             単位胞の繰り返し数 (デフォルト: 3)
#    --temp=<K>              目標温度 [K] (デフォルト: 298.0)
#    --pres=<GPa>            目標圧力 [GPa] (デフォルト: 0.0)
#    --step=<N>              本計算ステップ数 (デフォルト: 10000)
#    --dt=<fs>               時間刻み [fs] (LJ: 1.0, MMMD: 0.1, AIREBO: 0.5)
#    --init_scale=<s>        格子定数スケール因子 (デフォルト: 1.0)
#    --seed=<n>              乱数シード (デフォルト: 42)
#    --coldstart=<N>         低温(4K)ステップ数 (デフォルト: 0)
#    --warmup=<N>            昇温ステップ数 4K→T (デフォルト: 0)
#    --from=<step>           平均開始ステップ (デフォルト: 本計算の3/4地点)
#    --to=<step>             平均終了ステップ (デフォルト: nsteps)
#    --mon=<N>               モニタリング出力間隔 (デフォルト: 自動)
#    --warmup_mon=<mode>     昇温中の出力頻度 norm|freq|some (デフォルト: norm)
#    --ovito=<N>             OVITO XYZ出力間隔 (0=無効, デフォルト: 0)
#    --ofile=<filename>      OVITO出力ファイル名 (LJ版のみ, デフォルト: 自動)
#    --restart=<N>           リスタート保存間隔 (0=無効, デフォルト: 0)
#    --resfile=<path>        リスタートファイルから再開
#    --libdir=<path>         フラーレンライブラリ (デフォルト: FullereneLib)
#
#  分子力学版 [4] 追加オプション:
#    --ff_kb=<kcal/mol>      結合伸縮力定数 (デフォルト: 469.0)
#    --ff_kth=<kcal/mol>     角度曲げ力定数 (デフォルト: 63.0)
#    --ff_v2=<kcal/mol>      二面角力定数 (デフォルト: 14.5)
#    --ff_kimp=<kcal/mol>    不適切二面角力定数 (デフォルト: 15.0)
#
#---------------------------------------------------------------------------
#  全オプションを網羅した実行例:
#
#  ● LJ剛体版 (dt=1.0fs):
#    bin/fuller_LJ_npt_md_omp \
#      --fullerene=C60 --crystal=fcc --cell=3 \
#      --temp=298 --pres=0.0 --step=50000 --dt=1.0 \
#      --init_scale=1.0 --seed=42 \
#      --coldstart=2000 --warmup=3000 \
#      --from=40000 --to=50000 --mon=500 --warmup_mon=norm \
#      --ovito=100 --ofile=my_traj.xyz \
#      --restart=5000 --libdir=FullereneLib
#
#  ● LJ剛体版 リスタート再開:
#    bin/fuller_LJ_npt_md_omp --resfile=restart_LJ_omp_00025000.rst
#
#  ● 分子力学版 (dt=0.1fs):
#    bin/fuller_LJ_npt_mmmd_omp \
#      --fullerene=C70 --crystal=fcc --cell=3 \
#      --temp=500 --pres=1.0 --step=100000 --dt=0.1 \
#      --init_scale=1.0 --seed=123 \
#      --coldstart=5000 --warmup=5000 \
#      --from=80000 --to=100000 --mon=1000 --warmup_mon=freq \
#      --ovito=200 --restart=10000 --libdir=FullereneLib \
#      --ff_kb=469 --ff_kth=63 --ff_v2=14.5 --ff_kimp=15.0
#
#  ● 分子力学版 リスタート再開:
#    bin/fuller_LJ_npt_mmmd_omp --resfile=restart_mmmd_omp_00050000.rst
#
#  ● AIREBO版 (dt=0.5fs):
#    bin/fuller_airebo_npt_md_omp \
#      --fullerene=C84:20:Td --crystal=fcc --cell=3 \
#      --temp=300 --pres=0.0 --step=50000 --dt=0.5 \
#      --init_scale=1.0 --seed=99 \
#      --coldstart=3000 --warmup=2000 \
#      --from=40000 --to=50000 --mon=500 --warmup_mon=some \
#      --ovito=100 --restart=5000 --libdir=FullereneLib
#
#  ● AIREBO版 リスタート再開:
#    bin/fuller_airebo_npt_md_omp --resfile=restart_airebo_omp_00025000.rst
#
#  ● 異なるフラーレン種の実行例:
#    bin/fuller_LJ_npt_md_omp --fullerene=C70 --step=10000
#    bin/fuller_LJ_npt_md_omp --fullerene=C72 --step=10000
#    bin/fuller_LJ_npt_md_omp --fullerene=C74 --step=10000
#    bin/fuller_LJ_npt_md_omp --fullerene=C76:D2 --step=10000
#    bin/fuller_LJ_npt_md_omp --fullerene=C76:Td --step=10000
#    bin/fuller_LJ_npt_md_omp --fullerene=C84:20:Td --step=10000
#    bin/fuller_LJ_npt_md_omp --fullerene=C84:23:D2d --step=10000
#    bin/fuller_LJ_npt_md_omp --fullerene=C84:22:D2 --step=10000
#
#  ● 停止制御:
#    mkdir abort.md   # 即座に停止 (リスタート有効時は保存して終了)
#    mkdir stop.md    # 次のリスタートチェックポイントで停止
#
#===========================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BINDIR="$BASE_DIR/bin"

# テストはfuller_md/ディレクトリで実行する (FullereneLib参照のため)
cd "$BASE_DIR"

echo "================================================================"
echo "  Test_fuller.sh — 動作検証"
echo "================================================================"
echo "  Working dir : $(pwd)"
echo "  Bin dir     : ${BINDIR}/"
echo "================================================================"
echo ""

PASS=0
FAIL=0
SKIP=0

#--- テスト関数 ------------------------------------------------------------
run_test() {
    local name="$1"
    local exe="$2"
    shift 2
    local args=("$@")

    printf "  %-52s ... " "$name"
    if [ ! -x "$exe" ]; then
        echo "SKIP (not built)"
        SKIP=$((SKIP + 1))
        return
    fi
    if "$exe" "${args[@]}" > /tmp/test_fuller_out.txt 2>&1; then
        # 出力に "Done" が含まれるか確認 (正常終了の目印)
        if grep -q "Done" /tmp/test_fuller_out.txt; then
            echo "PASS"
            PASS=$((PASS + 1))
        else
            echo "FAIL (no 'Done' in output)"
            tail -3 /tmp/test_fuller_out.txt
            FAIL=$((FAIL + 1))
        fi
    else
        echo "FAIL (exit code $?)"
        tail -3 /tmp/test_fuller_out.txt
        FAIL=$((FAIL + 1))
    fi
}

#===========================================================================
#  [1] fuller_LJ_core_serial_pure — コア版 Serial専用
#      引数: [nc] のみ。リスタートなし、OVITOなし。
#===========================================================================
echo "[1/5] LJ Core Serial Pure (100 steps)"
run_test "fuller_LJ_core_serial_pure" \
    "$BINDIR/fuller_LJ_core_serial_pure"
echo ""

#===========================================================================
#  [2] fuller_LJ_core — コア版 Serial/OMP/ACC
#      引数: [nc] のみ。リスタートなし、OVITOなし。
#===========================================================================
echo "[2/5] LJ Core Serial/OMP (100 steps)"
run_test "fuller_LJ_core_serial" \
    "$BINDIR/fuller_LJ_core_serial"
run_test "fuller_LJ_core_omp" \
    "$BINDIR/fuller_LJ_core_omp"
run_test "fuller_LJ_core_gpu" \
    "$BINDIR/fuller_LJ_core_gpu"
echo ""

#===========================================================================
#  [3] fuller_LJ_npt_md — LJ剛体 フル版
#      最小ステップ (step=200) でリスタート・OVITO含む基本検証。
#===========================================================================
echo "[3/5] LJ Rigid-body Full (200 steps, with restart+ovito)"
run_test "fuller_LJ_npt_md_serial" \
    "$BINDIR/fuller_LJ_npt_md_serial" \
    --step=200 --mon=100
run_test "fuller_LJ_npt_md_omp" \
    "$BINDIR/fuller_LJ_npt_md_omp" \
    --step=200 --mon=100
run_test "fuller_LJ_npt_md_gpu" \
    "$BINDIR/fuller_LJ_npt_md_gpu" \
    --step=200 --mon=100
run_test "fuller_LJ_npt_md_serial (ovito)" \
    "$BINDIR/fuller_LJ_npt_md_serial" \
    --step=200 --mon=100 --ovito=100
run_test "fuller_LJ_npt_md_serial (restart)" \
    "$BINDIR/fuller_LJ_npt_md_serial" \
    --step=200 --mon=100 --restart=100
echo ""

#===========================================================================
#  [4] fuller_LJ_npt_mmmd — 分子力学 フル版
#      最小ステップ (step=200) でリスタート・OVITO含む基本検証。
#===========================================================================
echo "[4/5] Molecular Mechanics Full (200 steps, with restart+ovito)"
run_test "fuller_LJ_npt_mmmd_serial" \
    "$BINDIR/fuller_LJ_npt_mmmd_serial" \
    --step=200 --mon=100
run_test "fuller_LJ_npt_mmmd_omp" \
    "$BINDIR/fuller_LJ_npt_mmmd_omp" \
    --step=200 --mon=100
run_test "fuller_LJ_npt_mmmd_gpu" \
    "$BINDIR/fuller_LJ_npt_mmmd_gpu" \
    --step=200 --mon=100
run_test "fuller_LJ_npt_mmmd_serial (ovito)" \
    "$BINDIR/fuller_LJ_npt_mmmd_serial" \
    --step=200 --mon=100 --ovito=100
run_test "fuller_LJ_npt_mmmd_serial (restart)" \
    "$BINDIR/fuller_LJ_npt_mmmd_serial" \
    --step=200 --mon=100 --restart=100
echo ""

#===========================================================================
#  [5] fuller_airebo_npt_md — AIREBO フル版
#      最小ステップ (step=200) でリスタート・OVITO含む基本検証。
#===========================================================================
echo "[5/5] AIREBO Full (200 steps, with restart+ovito)"
run_test "fuller_airebo_npt_md_serial" \
    "$BINDIR/fuller_airebo_npt_md_serial" \
    --step=200 --mon=100
run_test "fuller_airebo_npt_md_omp" \
    "$BINDIR/fuller_airebo_npt_md_omp" \
    --step=200 --mon=100
run_test "fuller_airebo_npt_md_gpu" \
    "$BINDIR/fuller_airebo_npt_md_gpu" \
    --step=200 --mon=100
run_test "fuller_airebo_npt_md_serial (ovito)" \
    "$BINDIR/fuller_airebo_npt_md_serial" \
    --step=200 --mon=100 --ovito=100
run_test "fuller_airebo_npt_md_serial (restart)" \
    "$BINDIR/fuller_airebo_npt_md_serial" \
    --step=200 --mon=100 --restart=100
echo ""

#--- テスト生成ファイルの清掃 -----------------------------------------------
echo "Cleaning up test output files..."
rm -f ovito_traj_*.xyz restart_*.rst
echo ""

#--- 結果サマリ ------------------------------------------------------------
echo "================================================================"
echo "  Test complete:  PASS=$PASS  FAIL=$FAIL  SKIP=$SKIP"
echo "================================================================"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
exit 0
