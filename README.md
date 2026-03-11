# fuller_md — Fullerene Crystal NPT Molecular Dynamics

C60/C70/C72/C74/C76/C84 フラーレン結晶の NPT 分子動力学シミュレーションコード。
3種の力場モデル（LJ剛体・分子力学・AIREBO）を実装し、
Serial / OpenMP / OpenACC GPU の3モードをコンパイル時スイッチで切替可能。

## ディレクトリ構造

```
fuller_md/
├── README.md              ← このファイル
├── FullereneLib/          ← フラーレン分子座標データ (.cc1)
│   ├── C60-76/            ←   C60(Ih), C70(D5h), C72(D6d), C74(D3h), C76(D2,Td)
│   └── C84/               ←   C84 No.01〜No.24 (24異性体)
├── src/                   ← ソースコード・スクリプト
│   ├── Build_fuller.sh    ←   ビルドスクリプト
│   ├── Test_fuller.sh     ←   動作検証スクリプト
│   ├── fuller_LJ_npt_md_core_serial.cpp          [1] LJ剛体コア版 (Serial専用)
│   ├── fuller_LJ_npt_md_core_serial_omp_acc.cpp   [2] LJ剛体コア版 (Serial/OMP/ACC)
│   ├── fuller_LJ_npt_md_serial_omp_acc.cpp        [3] LJ剛体フル版 (Serial/OMP/ACC)
│   ├── fuller_LJ_npt_mmmd_serial_omp_acc.cpp      [4] 分子力学フル版 (Serial/OMP/ACC)
│   └── fuller_airebo_npt_md_serial_omp_acc.cpp    [5] AIREBOフル版 (Serial/OMP/ACC)
└── bin/                   ← 実行モジュール出力先 (ビルド時に自動作成)
```

## ソースファイル一覧

### コア版 [1][2] — 学習・ベンチマーク用

パラメータはソースコード内に固定（T=300K, P=0GPa, dt=1fs, 1000ステップ）。
引数は `nc`（セルサイズ）のみ。リスタート機能・OVITO出力なし。

| # | ファイル | 説明 | 並列化 |
|---|---------|------|--------|
| 1 | `fuller_LJ_npt_md_core_serial.cpp` | LJ剛体 Serial専用。並列化ガイドコメント付き | Serial |
| 2 | `fuller_LJ_npt_md_core_serial_omp_acc.cpp` | LJ剛体 3モード統合。`#ifdef`による切替 | Serial/OMP/ACC |

### フル版 [3][4][5] — 本番計算用

全ランタイムオプション対応。リスタート保存/再開、OVITO XYZ出力、停止制御（abort.md/stop.md）を実装。

| # | ファイル | 力場 | dt既定値 |
|---|---------|------|----------|
| 3 | `fuller_LJ_npt_md_serial_omp_acc.cpp` | LJ 剛体分子間ポテンシャル | 1.0 fs |
| 4 | `fuller_LJ_npt_mmmd_serial_omp_acc.cpp` | 分子力学 (Bond+Angle+Dihedral+Improper+LJ+Coulomb) | 0.1 fs |
| 5 | `fuller_airebo_npt_md_serial_omp_acc.cpp` | AIREBO (REBO-II + LJ) | 0.5 fs |

## ビルド方法

### 前提条件

- **Serial/OpenMP**: GCC 12以上 (g++) または同等のC++17対応コンパイラ
- **OpenACC GPU**: NVIDIA HPC SDK (nvc++)
- macOSの場合: Homebrew GCC (`brew install gcc`) が必要（clang++は`-fopenmp`非対応）

### ビルドスクリプト

```bash
cd fuller_md

# Serial + OpenMP をビルド（デフォルト）
src/Build_fuller.sh

# Serial のみ
src/Build_fuller.sh serial

# OpenMP のみ
src/Build_fuller.sh omp

# OpenACC GPU のみ
src/Build_fuller.sh acc

# 全モード (Serial + OpenMP + OpenACC)
src/Build_fuller.sh all

# ビルド成果物を削除
src/Build_fuller.sh clean
```

スクリプトはmacOS上のHomebrew GCC (`g++-15`, `g++-14`, ...) を自動検出する。
環境変数 `CXX` でコンパイラを明示指定することも可能:

```bash
CXX=/opt/gcc/bin/g++ src/Build_fuller.sh
```

### 生成される実行モジュール

| 実行モジュール | ソース | モード |
|---------------|--------|--------|
| `fuller_LJ_core_serial_pure` | [1] | Serial |
| `fuller_LJ_core_serial` | [2] | Serial |
| `fuller_LJ_core_omp` | [2] | OpenMP |
| `fuller_LJ_core_gpu` | [2] | OpenACC |
| `fuller_LJ_npt_md_serial` | [3] | Serial |
| `fuller_LJ_npt_md_omp` | [3] | OpenMP |
| `fuller_LJ_npt_md_gpu` | [3] | OpenACC |
| `fuller_LJ_npt_mmmd_serial` | [4] | Serial |
| `fuller_LJ_npt_mmmd_omp` | [4] | OpenMP |
| `fuller_LJ_npt_mmmd_gpu` | [4] | OpenACC |
| `fuller_airebo_npt_md_serial` | [5] | Serial |
| `fuller_airebo_npt_md_omp` | [5] | OpenMP |
| `fuller_airebo_npt_md_gpu` | [5] | OpenACC |

## 動作検証

```bash
cd fuller_md
src/Test_fuller.sh
```

全実行モジュールを最小ステップ数で実行し、OVITO出力・リスタート保存を含む基本動作を検証する。
GPU版はビルドされていなければ自動的にSKIPされる。

## 実行方法

### コア版

```bash
cd fuller_md

# デフォルト (3x3x3, N=108分子, 1000ステップ)
bin/fuller_LJ_core_serial_pure

# セルサイズ指定 (5x5x5, N=500分子)
bin/fuller_LJ_core_omp 5
bin/fuller_LJ_core_omp --cell=5
```

### フル版 — 基本実行

```bash
cd fuller_md

# LJ剛体 (デフォルト: C60, FCC 3x3x3, 298K, 0GPa, 10000ステップ)
bin/fuller_LJ_npt_md_serial

# 温度・圧力・ステップ数を指定
bin/fuller_LJ_npt_md_omp --temp=500 --pres=1.0 --step=50000

# 低温開始(4K) + 昇温 + 本計算
bin/fuller_LJ_npt_md_serial --coldstart=2000 --warmup=3000 --step=20000
```

### フル版 — OVITO出力

指定間隔でXYZ形式のトラジェクトリファイルを出力する。OVITOで可視化可能。

```bash
# 100ステップ毎に出力
bin/fuller_LJ_npt_md_omp --step=10000 --ovito=100
# → ovito_traj_LJ_omp_*.xyz が生成される

# 分子力学版
bin/fuller_LJ_npt_mmmd_omp --step=20000 --ovito=200

# AIREBO版
bin/fuller_airebo_npt_md_omp --step=10000 --ovito=100
```

### フル版 — リスタート

指定間隔でリスタートファイルを保存する。最終ステップでも自動保存。

```bash
# 5000ステップ毎にリスタート保存
bin/fuller_LJ_npt_md_serial --step=50000 --restart=5000
# → restart_LJ_serial_*_00005000.rst, ..._00010000.rst, ... が生成

# リスタートファイルから再開
bin/fuller_LJ_npt_md_serial --resfile=restart_LJ_serial_00025000.rst

# OVITO + リスタートを同時使用
bin/fuller_LJ_npt_md_omp --step=100000 --ovito=500 --restart=10000
```

### フル版 — 停止制御

実行中にカレントディレクトリにディレクトリを作成して制御:

```bash
# 即座に停止（リスタート有効時は保存してから終了）
mkdir abort.md

# 次のリスタートチェックポイントで停止
mkdir stop.md
```

### フル版 — ヘルプ表示

```bash
bin/fuller_LJ_npt_md_serial --help
bin/fuller_LJ_npt_mmmd_serial --help
bin/fuller_airebo_npt_md_serial --help
```

## ランタイムオプション一覧（フル版共通）

| オプション | 説明 | デフォルト |
|-----------|------|-----------|
| `--help` | ヘルプ表示 | — |
| `--fullerene=<名前>` | フラーレン種 | C60 |
| `--crystal=<fcc\|hcp\|bcc>` | 結晶構造 | fcc |
| `--cell=<nc>` | 単位胞の繰り返し数 | 3 |
| `--temp=<K>` | 目標温度 [K] | 298.0 |
| `--pres=<GPa>` | 目標圧力 [GPa] | 0.0 |
| `--step=<N>` | 本計算ステップ数 | 10000 |
| `--dt=<fs>` | 時間刻み [fs] | 力場依存 |
| `--init_scale=<s>` | 格子定数スケール因子 | 1.0 |
| `--seed=<n>` | 乱数シード | 42 |
| `--coldstart=<N>` | 低温(4K)ステップ数 | 0 |
| `--warmup=<N>` | 昇温ステップ数 4K→T | 0 |
| `--from=<step>` | 平均開始ステップ | 自動(3/4地点) |
| `--to=<step>` | 平均終了ステップ | nsteps |
| `--mon=<N>` | モニタリング出力間隔 | 自動 |
| `--warmup_mon=<mode>` | 昇温中の出力頻度 (norm\|freq\|some) | norm |
| `--ovito=<N>` | OVITO XYZ出力間隔 (0=無効) | 0 |
| `--ofile=<filename>` | OVITO出力ファイル名 (LJ版のみ) | 自動生成 |
| `--restart=<N>` | リスタート保存間隔 (0=無効) | 0 |
| `--resfile=<path>` | リスタートファイルから再開 | — |
| `--libdir=<path>` | フラーレンライブラリ | FullereneLib |

### 分子力学版 [4] 追加オプション

| オプション | 説明 | デフォルト |
|-----------|------|-----------|
| `--ff_kb=<kcal/mol>` | 結合伸縮力定数 | 469.0 |
| `--ff_kth=<kcal/mol>` | 角度曲げ力定数 | 63.0 |
| `--ff_v2=<kcal/mol>` | 二面角力定数 | 14.5 |
| `--ff_kimp=<kcal/mol>` | 不適切二面角力定数 | 15.0 |

## 選択可能なフラーレン

`--fullerene=` オプションで指定。座標データは `FullereneLib/` に格納。

### C60-76 系列

| 指定値 | 原子数 | 対称群 |
|--------|--------|--------|
| `C60` (デフォルト) | 60 | Ih |
| `C70` | 70 | D5h |
| `C72` | 72 | D6d |
| `C74` | 74 | D3h |
| `C76:D2` | 76 | D2 |
| `C76:Td` | 76 | Td |

### C84 系列 (24異性体)

`--fullerene=C84:<番号>` または `--fullerene=C84:<番号>:<対称群>` で指定。

| 指定値 | 対称群 | 備考 |
|--------|--------|------|
| `C84:1` | D2 | |
| `C84:2` | C2 | |
| `C84:3` | Cs | |
| `C84:4` | D2d | |
| `C84:5` | D2 | |
| `C84:6` | C2v | |
| `C84:7` | C2v | |
| `C84:8` | C2 | |
| `C84:9` | C2 | |
| `C84:10` | Cs | |
| `C84:11` | C2 | |
| `C84:12` | C1 | |
| `C84:13` | C2 | |
| `C84:14` | Cs | |
| `C84:15` | Cs | |
| `C84:16` | Cs | |
| `C84:17` | C2v | |
| `C84:18` | C2v | |
| `C84:19` | D3d | |
| `C84:20` | Td | 安定異性体 |
| `C84:21` | D2 | |
| `C84:22` | D2 | 安定異性体 |
| `C84:23` | D2d | 安定異性体 |
| `C84:24` | D6h | |

## 物理モデル

- **アンサンブル**: NPT (等温等圧)
- **熱浴**: Nose-Hoover chain
- **圧力制御**: Parrinello-Rahman
- **時間積分**: Velocity-Verlet
- **周期境界条件**: 3次元 (triclinic cell)
- **近接リスト**: 対称フルリスト (Newton第3法則不使用、GPU競合回避)

### LJ剛体版 [1][2][3]

- 分子間: Lennard-Jones (C-C, sigma=3.4A, epsilon=2.63meV)
- 剛体回転: 四元数表現
- 自由度: 重心並進(3) + 回転(3) × N分子

### 分子力学版 [4]

- 分子内: Bond stretching + Angle bending + Dihedral + Improper
- 分子間: LJ + Coulomb (電荷がある場合)
- 全原子自由度

### AIREBO版 [5]

- REBO-II (Brenner 2002): 共有結合 (結合次数ポテンシャル)
- LJ (Stuart 2000): 分子間 van der Waals
- 全原子自由度

## 単位系

| 物理量 | 単位 |
|--------|------|
| 距離 | A (オングストローム) |
| 質量 | amu (原子質量単位) |
| エネルギー | eV (電子ボルト) |
| 時間 | fs (フェムト秒) |
| 温度 | K (ケルビン) |
| 圧力 | GPa (ギガパスカル) |

## 動作環境

- macOS (Homebrew GCC) — Serial/OpenMP
- Linux (GCC 12+) — Serial/OpenMP
- Linux (NVIDIA HPC SDK + GPU) — Serial/OpenMP/OpenACC
- 検証済み: FOCUS スパコン (GCC 12.2.0 + NVIDIA HPC SDK 24.5)

## ライセンス

本プロジェクトは [BSD 3-Clause License](LICENSE) の下で公開されています。

## 他言語・他バージョン

| リポジトリ | 言語 | 説明 |
|-----------|------|------|
| [fuller_md](https://github.com/focusnishikawa/fuller_md) | C++ (日本語) | このリポジトリ |
| [fuller_md_en](https://github.com/focusnishikawa/fuller_md_en) | C++ (English) | C++版 英語 |
| [fuller_md_Julia](https://github.com/focusnishikawa/fuller_md_Julia) | Julia (English) | Julia版 英語 |
| [fuller_md_Julia_ja](https://github.com/focusnishikawa/fuller_md_Julia_ja) | Julia (日本語) | Julia版 日本語 |
