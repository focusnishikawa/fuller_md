#!/bin/bash
#===========================================================================
#  Build_fuller.sh — フラーレン結晶 NPT-MD ビルドスクリプト
#
#  使い方:
#    ./Build_fuller.sh              # Serial + OpenMP 全ファイルビルド
#    ./Build_fuller.sh serial       # Serial のみ
#    ./Build_fuller.sh omp          # OpenMP のみ
#    ./Build_fuller.sh acc          # OpenACC GPU のみ
#    ./Build_fuller.sh all          # Serial + OpenMP + OpenACC 全部
#    ./Build_fuller.sh clean        # bin/ 内の実行ファイルを削除
#
#  ディレクトリ構造:
#    fuller_md/
#    ├── src/           ← このスクリプトとソースコード
#    ├── bin/           ← 実行モジュール出力先
#    └── FullereneLib/  ← フラーレン座標データ
#
#---------------------------------------------------------------------------
#  ビルド対象ソースファイル (5本):
#
#  [1] fuller_LJ_npt_md_core_serial.cpp         — LJ剛体 コア版 (Serial専用)
#  [2] fuller_LJ_npt_md_core_serial_omp_acc.cpp — LJ剛体 コア版 (Serial/OMP/ACC)
#  [3] fuller_LJ_npt_md_serial_omp_acc.cpp      — LJ剛体 フル版 (Serial/OMP/ACC)
#  [4] fuller_LJ_npt_mmmd_serial_omp_acc.cpp    — 分子力学 フル版 (Serial/OMP/ACC)
#  [5] fuller_airebo_npt_md_serial_omp_acc.cpp  — AIREBO フル版 (Serial/OMP/ACC)
#
#  コア版 [1][2]: パラメータ固定、引数は nc (セルサイズ) のみ
#                 リスタート機能なし、OVITO出力なし
#  フル版 [3][4][5]: 全ランタイムオプション対応
#                    リスタート保存/再開、OVITO XYZ出力対応
#
#---------------------------------------------------------------------------
#  実行例 — コア版 [1][2]:
#
#    bin/fuller_LJ_core_serial_pure       # 3x3x3 (N=108) デフォルト
#    bin/fuller_LJ_core_serial_pure 5     # 5x5x5 (N=500)
#    bin/fuller_LJ_core_omp 4             # OpenMP, 4x4x4 (N=256)
#
#---------------------------------------------------------------------------
#  実行例 — LJ剛体 フル版 [3]:
#
#    # 基本実行 (C60 FCC 3x3x3, 298K, 10000ステップ)
#    bin/fuller_LJ_npt_md_serial
#
#    # 温度・圧力・ステップ数を指定
#    bin/fuller_LJ_npt_md_omp --temp=500 --pres=1.0 --step=50000
#
#    # 低温開始 + 昇温 + 本計算
#    bin/fuller_LJ_npt_md_serial --coldstart=2000 --warmup=3000 --step=20000
#
#    # OVITO XYZ出力 (100ステップ毎にトラジェクトリ書き出し)
#    bin/fuller_LJ_npt_md_omp --step=10000 --ovito=100
#
#    # リスタート保存 (5000ステップ毎 + 最終ステップ)
#    bin/fuller_LJ_npt_md_serial --step=50000 --restart=5000
#
#    # リスタートファイルから再開
#    bin/fuller_LJ_npt_md_serial --resfile=restart_LJ_serial_00025000.rst
#
#    # OVITO + リスタートを同時使用
#    bin/fuller_LJ_npt_md_omp --step=100000 --ovito=500 --restart=10000
#
#    # 全オプション一覧
#    bin/fuller_LJ_npt_md_serial --help
#
#---------------------------------------------------------------------------
#  実行例 — 分子力学 フル版 [4]:
#
#    bin/fuller_LJ_npt_mmmd_serial --step=20000
#    bin/fuller_LJ_npt_mmmd_omp --temp=500 --step=100000 --ovito=200
#    bin/fuller_LJ_npt_mmmd_serial --step=100000 --restart=10000
#    bin/fuller_LJ_npt_mmmd_serial --resfile=restart_mmmd_serial_00050000.rst
#    bin/fuller_LJ_npt_mmmd_omp --ff_kb=500 --ff_kth=70 --step=20000
#    bin/fuller_LJ_npt_mmmd_serial --help
#
#---------------------------------------------------------------------------
#  実行例 — AIREBO フル版 [5]:
#
#    bin/fuller_airebo_npt_md_serial --step=10000
#    bin/fuller_airebo_npt_md_omp --temp=500 --step=50000 --ovito=100
#    bin/fuller_airebo_npt_md_serial --step=50000 --restart=5000
#    bin/fuller_airebo_npt_md_serial --resfile=restart_airebo_serial_00025000.rst
#    bin/fuller_airebo_npt_md_omp --fullerene=C84:20:Td --cell=4 --step=20000
#    bin/fuller_airebo_npt_md_serial --help
#
#---------------------------------------------------------------------------
#  停止制御 (フル版 [3][4][5] 共通):
#    実行中にカレントディレクトリに以下を作成:
#      mkdir abort.md   → 即座に停止 (リスタート有効時は保存して終了)
#      mkdir stop.md    → 次のリスタートチェックポイントで停止
#
#===========================================================================
set -e

# スクリプトのあるディレクトリを基準にパスを設定
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SRCDIR="$SCRIPT_DIR"
BINDIR="$BASE_DIR/bin"
CXXSTD="-std=c++17"
OPT="-O3"

#--- コンパイラ検出 --------------------------------------------------------
# macOS の g++ は clang のエイリアスで -fopenmp 非対応のため、
# Homebrew の実 GCC を探す
find_gxx() {
    if [ -n "$CXX" ]; then
        echo "$CXX"
        return
    fi
    for v in 15 14 13 12; do
        for p in /opt/homebrew/bin /usr/local/bin; do
            if [ -x "$p/g++-$v" ]; then
                echo "$p/g++-$v"
                return
            fi
        done
    done
    echo "g++"
}

GXX=$(find_gxx)
NVCXX="nvc++"

#--- ビルドモード判定 ------------------------------------------------------
MODE="${1:-serial_omp}"
case "$MODE" in
    serial)     DO_SERIAL=1; DO_OMP=0; DO_ACC=0 ;;
    omp)        DO_SERIAL=0; DO_OMP=1; DO_ACC=0 ;;
    acc|gpu)    DO_SERIAL=0; DO_OMP=0; DO_ACC=1 ;;
    all)        DO_SERIAL=1; DO_OMP=1; DO_ACC=1 ;;
    serial_omp) DO_SERIAL=1; DO_OMP=1; DO_ACC=0 ;;
    clean)
        echo "Cleaning ${BINDIR}/ ..."
        rm -f "${BINDIR}"/*
        echo "Done."
        exit 0
        ;;
    *)
        echo "Usage: $0 [serial|omp|acc|all|clean]"
        exit 1
        ;;
esac

mkdir -p "$BINDIR"

#--- コンパイラ確認 --------------------------------------------------------
echo "================================================================"
echo "  Build_fuller.sh"
echo "================================================================"
echo "  g++ compiler : $GXX"
if [ "$DO_ACC" -eq 1 ]; then
    if ! command -v "$NVCXX" &>/dev/null; then
        echo "  WARNING: $NVCXX not found — OpenACC builds will be skipped"
        DO_ACC=0
    else
        echo "  nvc++ compiler: $NVCXX"
    fi
fi
echo "  Build mode   : serial=$DO_SERIAL omp=$DO_OMP acc=$DO_ACC"
echo "  Source dir   : $SRCDIR/"
echo "  Output dir   : $BINDIR/"
echo "================================================================"
echo ""

OK=0
FAIL=0

#--- ビルド関数 ------------------------------------------------------------
build() {
    local compiler="$1"
    local flags="$2"
    local src="$3"
    local out="$4"

    printf "  %-50s ... " "$out"
    if $compiler $CXXSTD $OPT $flags -o "${BINDIR}/${out}" "${SRCDIR}/${src}" -lm 2>/tmp/build_fuller_err.txt; then
        echo "OK"
        OK=$((OK + 1))
    else
        echo "FAIL"
        cat /tmp/build_fuller_err.txt | head -5
        FAIL=$((FAIL + 1))
    fi
}

#===========================================================================
#  1) fuller_LJ_npt_md_core_serial.cpp  (Serial専用)
#===========================================================================
echo "[1/5] fuller_LJ_npt_md_core_serial.cpp (Serial only)"
if [ "$DO_SERIAL" -eq 1 ]; then
    build "$GXX" "" \
        "fuller_LJ_npt_md_core_serial.cpp" \
        "fuller_LJ_core_serial_pure"
fi
echo ""

#===========================================================================
#  2) fuller_LJ_npt_md_core_serial_omp_acc.cpp  (Serial/OMP/ACC)
#===========================================================================
echo "[2/5] fuller_LJ_npt_md_core_serial_omp_acc.cpp"
if [ "$DO_SERIAL" -eq 1 ]; then
    build "$GXX" "-Wno-unknown-pragmas" \
        "fuller_LJ_npt_md_core_serial_omp_acc.cpp" \
        "fuller_LJ_core_serial"
fi
if [ "$DO_OMP" -eq 1 ]; then
    build "$GXX" "-fopenmp -Wno-unknown-pragmas" \
        "fuller_LJ_npt_md_core_serial_omp_acc.cpp" \
        "fuller_LJ_core_omp"
fi
if [ "$DO_ACC" -eq 1 ]; then
    build "$NVCXX" "-acc -gpu=cc80 -Minfo=accel" \
        "fuller_LJ_npt_md_core_serial_omp_acc.cpp" \
        "fuller_LJ_core_gpu"
fi
echo ""

#===========================================================================
#  3) fuller_LJ_npt_md_serial_omp_acc.cpp  (Serial/OMP/ACC)
#===========================================================================
echo "[3/5] fuller_LJ_npt_md_serial_omp_acc.cpp"
if [ "$DO_SERIAL" -eq 1 ]; then
    build "$GXX" "-Wno-unknown-pragmas" \
        "fuller_LJ_npt_md_serial_omp_acc.cpp" \
        "fuller_LJ_npt_md_serial"
fi
if [ "$DO_OMP" -eq 1 ]; then
    build "$GXX" "-fopenmp -Wno-unknown-pragmas" \
        "fuller_LJ_npt_md_serial_omp_acc.cpp" \
        "fuller_LJ_npt_md_omp"
fi
if [ "$DO_ACC" -eq 1 ]; then
    build "$NVCXX" "-acc -gpu=cc80 -Minfo=accel" \
        "fuller_LJ_npt_md_serial_omp_acc.cpp" \
        "fuller_LJ_npt_md_gpu"
fi
echo ""

#===========================================================================
#  4) fuller_LJ_npt_mmmd_serial_omp_acc.cpp  (Serial/OMP/ACC)
#===========================================================================
echo "[4/5] fuller_LJ_npt_mmmd_serial_omp_acc.cpp"
if [ "$DO_SERIAL" -eq 1 ]; then
    build "$GXX" "-Wno-unknown-pragmas" \
        "fuller_LJ_npt_mmmd_serial_omp_acc.cpp" \
        "fuller_LJ_npt_mmmd_serial"
fi
if [ "$DO_OMP" -eq 1 ]; then
    build "$GXX" "-fopenmp -Wno-unknown-pragmas" \
        "fuller_LJ_npt_mmmd_serial_omp_acc.cpp" \
        "fuller_LJ_npt_mmmd_omp"
fi
if [ "$DO_ACC" -eq 1 ]; then
    build "$NVCXX" "-acc -gpu=cc80 -Minfo=accel" \
        "fuller_LJ_npt_mmmd_serial_omp_acc.cpp" \
        "fuller_LJ_npt_mmmd_gpu"
fi
echo ""

#===========================================================================
#  5) fuller_airebo_npt_md_serial_omp_acc.cpp  (Serial/OMP/ACC)
#===========================================================================
echo "[5/5] fuller_airebo_npt_md_serial_omp_acc.cpp"
if [ "$DO_SERIAL" -eq 1 ]; then
    build "$GXX" "-Wno-unknown-pragmas" \
        "fuller_airebo_npt_md_serial_omp_acc.cpp" \
        "fuller_airebo_npt_md_serial"
fi
if [ "$DO_OMP" -eq 1 ]; then
    build "$GXX" "-fopenmp -Wno-unknown-pragmas" \
        "fuller_airebo_npt_md_serial_omp_acc.cpp" \
        "fuller_airebo_npt_md_omp"
fi
if [ "$DO_ACC" -eq 1 ]; then
    build "$NVCXX" "-acc -gpu=cc80 -Minfo=accel" \
        "fuller_airebo_npt_md_serial_omp_acc.cpp" \
        "fuller_airebo_npt_md_gpu"
fi
echo ""

#--- 結果サマリ ------------------------------------------------------------
echo "================================================================"
echo "  Build complete:  OK=$OK  FAIL=$FAIL"
if [ "$OK" -gt 0 ]; then
    echo ""
    echo "  Executables in ${BINDIR}/:"
    ls -lh "${BINDIR}"/ 2>/dev/null | grep -v "^total" | awk '{printf "    %-40s %s\n", $NF, $5}'
fi
echo "================================================================"
echo ""
echo "実行例:"
echo ""
echo "  [コア版] パラメータ固定、引数は nc (セルサイズ) のみ"
echo "    ${BINDIR}/fuller_LJ_core_serial_pure         # Serial, デフォルト 3x3x3"
echo "    ${BINDIR}/fuller_LJ_core_serial_pure 5       # Serial, 5x5x5 (N=500)"
echo "    ${BINDIR}/fuller_LJ_core_omp 4               # OpenMP, 4x4x4 (N=256)"
echo ""
echo "  [LJ剛体 フル版] 全ランタイムオプション対応"
echo "    ${BINDIR}/fuller_LJ_npt_md_serial --help                          # ヘルプ表示"
echo "    ${BINDIR}/fuller_LJ_npt_md_omp --temp=500 --step=50000            # 温度指定"
echo "    ${BINDIR}/fuller_LJ_npt_md_serial --step=10000 --ovito=100        # OVITO出力"
echo "    ${BINDIR}/fuller_LJ_npt_md_serial --step=50000 --restart=5000     # リスタート保存"
echo "    ${BINDIR}/fuller_LJ_npt_md_serial --resfile=restart_*.rst         # リスタート再開"
echo ""
echo "  [分子力学 フル版]"
echo "    ${BINDIR}/fuller_LJ_npt_mmmd_serial --help                        # ヘルプ表示"
echo "    ${BINDIR}/fuller_LJ_npt_mmmd_omp --step=100000 --ovito=200        # OVITO出力"
echo "    ${BINDIR}/fuller_LJ_npt_mmmd_serial --step=100000 --restart=10000 # リスタート保存"
echo "    ${BINDIR}/fuller_LJ_npt_mmmd_serial --resfile=restart_*.rst       # リスタート再開"
echo ""
echo "  [AIREBO フル版]"
echo "    ${BINDIR}/fuller_airebo_npt_md_serial --help                      # ヘルプ表示"
echo "    ${BINDIR}/fuller_airebo_npt_md_omp --step=50000 --ovito=100       # OVITO出力"
echo "    ${BINDIR}/fuller_airebo_npt_md_serial --step=50000 --restart=5000 # リスタート保存"
echo "    ${BINDIR}/fuller_airebo_npt_md_serial --resfile=restart_*.rst     # リスタート再開"
echo ""
echo "  [停止制御] フル版 [3][4][5] 共通"
echo "    mkdir abort.md   # 即座に停止 (リスタート有効時は保存して終了)"
echo "    mkdir stop.md    # 次のリスタートチェックポイントで停止"
echo ""

exit $FAIL
