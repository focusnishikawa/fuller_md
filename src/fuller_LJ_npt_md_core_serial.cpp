// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025, Takeshi Nishikawa
/*===========================================================================
  fuller_LJ_npt_md_core_serial.cpp
  C60フラーレン結晶 NPT分子動力学シミュレーション
  (LJ剛体モデル・コア版 — シングルスレッド Serial)

  コンパイル:
    g++ -std=c++17 -O3 -o fuller_LJ_core_serial fuller_LJ_npt_md_core_serial.cpp -lm

  実行時オプション:
    ./fuller_LJ_core_serial [nc]

    nc (整数, デフォルト: 3, 最大: 8)
      FCC単位胞の繰り返し数。分子数 N = 4*nc^3。
      nc=3 → N=108, nc=4 → N=256, nc=5 → N=500

    指定方法:
      ./fuller_LJ_core_serial 3         # 位置引数
      ./fuller_LJ_core_serial --cell=5  # キーワード引数

  実行例:
    # デフォルト (3x3x3, N=108分子, 1000ステップ)
    ./fuller_LJ_core_serial

    # 大きなシステム (5x5x5, N=500分子)
    ./fuller_LJ_core_serial 5

  固定パラメータ (ソースコード内で変更):
    温度 T       = 300 K
    圧力 Pe      = 0.0 GPa
    時間刻み dt   = 1.0 fs
    ステップ数    = 1000
    出力間隔      = 100 ステップ
    近接リスト更新 = 25 ステップ

  機能:
    - C60剛体分子のLJ分子間相互作用によるNPT-MDシミュレーション
    - Nose-Hoover熱浴 + Parrinello-Rahman圧力制御
    - 四元数による剛体回転
    - Velocity-Verlet時間積分
    - FCC結晶構造の初期配置を自動生成
    - リスタート機能なし、OVITO出力なし (コア版)

  並列化について:
    このコードはシングルスレッドで動作する。
    各計算ループは分子単位で独立しているため、以下の並列化が可能:

    ● OpenMP化:
      各ループの前に #pragma omp parallel for を追加。
      ビリアルWm9への加算は #pragma omp atomic update で保護。
      コンパイル: g++ -std=c++17 -O3 -fopenmp ...

    ● OpenACC GPU化:
      配列をGPUに配置: #pragma acc data copy(pos[0:N*3], ...)
      各ループ: #pragma acc parallel loop present(...)
      デバイス関数: #pragma acc routine seq
      ビリアル: #pragma acc atomic update
      コンパイル: nvc++ -std=c++17 -O3 -acc -gpu=cc80 ...

    ● GPU最適化方針:
      forces(): gang並列(分子) × vector並列(分子内原子ai)
      - 対称フルリスト使用: Newton第3法則不使用→競合回避
      - ai原子ループをvector化: 最大128スレッドのwarp並列
      - ビリアル9成分をスカラーreduction変数に展開
      データ常駐: #pragma acc data で配列をGPU上に保持

  単位系: A (距離), amu (質量), eV (エネルギー), fs (時間), K (温度), GPa (圧力)
===========================================================================*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <string>


/* ═══════════════ 物理定数・単位変換 ═══════════════ */
/*  ニュートンの運動方程式 F = m × a を A/amu/eV/fs 単位系で使うための変換係数。
    a [A/fs^2] = (F [eV/A]) / (m [amu]) × CONV                              */
constexpr double CONV       = 9.64853321e-3;   // eV*fs^2/(amu*A^2)
constexpr double kB         = 8.617333262e-5;   // ボルツマン定数 [eV/K]
constexpr double eV2GPa     = 160.21766208;      // eV/A^3 → GPa
constexpr double eV2kcalmol = 23.06054783;       // eV → kcal/mol

/* ═══════════════ LJポテンシャルパラメータ ═══════════════ */
/*  V(r) = 4 eps [ (sigma/r)^12 - (sigma/r)^6 ]
    RCUT: カットオフ距離 (3*sigma ≒ 10.29 A)。r > RCUT の原子ペアは計算をスキップ。
    VSHFT: V(RCUT) を引いてカットオフでエネルギーが不連続にならないようにする。     */
constexpr double sigma_LJ = 3.431;              // C-C LJ sigma [A]
constexpr double eps_LJ   = 2.635e-3;           // C-C LJ epsilon [eV]
constexpr double RCUT      = 3.0*sigma_LJ;      // カットオフ距離 ~10.29 A
constexpr double RCUT2     = RCUT*RCUT;
constexpr double sig2_LJ   = sigma_LJ*sigma_LJ;
constexpr double mC        = 12.011;             // 炭素原子量 [amu]

/* VSHFT: constexpr化 — sigma/RCUT = 1/3, (1/3)^6 = 1/729 */
constexpr double _sr_v  = 1.0/3.0;
constexpr double _sr2_v = _sr_v * _sr_v;
constexpr double _sr6_v = _sr2_v * _sr2_v * _sr2_v;
constexpr double VSHFT  = 4.0*eps_LJ*(_sr6_v*_sr6_v - _sr6_v);

/* ═══════════════ 分子パラメータ ═══════════════ */
constexpr int C60_NATOM = 60;                    // 1分子あたりの原子数
constexpr int MAX_NATOM = 84;                    // C84対応の最大原子数
constexpr double MC60   = C60_NATOM * mC;        // C60分子量 [amu]
constexpr double RC60   = 3.55;                  // C60半径 [A]
constexpr double RMCUT   = RCUT + 2*RC60 + 1.0;  // 分子ペアカットオフ ~18.4 A
constexpr double RMCUT2  = RMCUT*RMCUT;

/*  MAX_NEIGH: 近傍リストの最大近傍数
    FCC結晶 (a0=14.17A) での近傍数:
      1st shell: 12 at a0/sqrt(2) = 10.0 A
      2nd shell:  6 at a0         = 14.2 A
      3rd shell: 24 at a0*sqrt(3/2) = 17.4 A
      合計: ~42 (対称リストなので同数)
    スキン3.0Aの余裕を含め、80で十分安全。                                     */
constexpr int MAX_NEIGH = 80;

/*  VECTOR_LENGTH: OpenACC GPU版でのgang内vectorスレッド数 (参考値)
    128 = 4 warp (NVIDIA GPU warp=32スレッド)
    C60 (natom=60): ceil(60/128)=1反復, C84 (natom=84): ceil(84/128)=1反復    */
constexpr int VECTOR_LENGTH = 128;


/* ═══════════════ H行列フラット9成分操作 ═══════════════ */
/*  h[0..8] = { H00,H01,H02, H10,H11,H12, H20,H21,H22 }
    行i,列j → h[3*i+j]
    シミュレーションセルの形状を表す3x3行列。
    体積 V = |det(H)|, 逆行列 Hi で実座標→分率座標の変換を行う。             */
#define H_(h,i,j) ((h)[3*(i)+(j)])

/* OpenACC化: #pragma acc routine seq を関数の前に追加してデバイス関数にする */
static inline double mat_det9(const double* h){
    return H_(h,0,0)*(H_(h,1,1)*H_(h,2,2)-H_(h,1,2)*H_(h,2,1))
          -H_(h,0,1)*(H_(h,1,0)*H_(h,2,2)-H_(h,1,2)*H_(h,2,0))
          +H_(h,0,2)*(H_(h,1,0)*H_(h,2,1)-H_(h,1,1)*H_(h,2,0));
}

static inline double mat_tr9(const double* h){
    return H_(h,0,0)+H_(h,1,1)+H_(h,2,2);
}

static void mat_inv9(const double* h, double* hi){
    double d=mat_det9(h), id=1.0/d;
    hi[0]=id*(H_(h,1,1)*H_(h,2,2)-H_(h,1,2)*H_(h,2,1));
    hi[1]=id*(H_(h,0,2)*H_(h,2,1)-H_(h,0,1)*H_(h,2,2));
    hi[2]=id*(H_(h,0,1)*H_(h,1,2)-H_(h,0,2)*H_(h,1,1));
    hi[3]=id*(H_(h,1,2)*H_(h,2,0)-H_(h,1,0)*H_(h,2,2));
    hi[4]=id*(H_(h,0,0)*H_(h,2,2)-H_(h,0,2)*H_(h,2,0));
    hi[5]=id*(H_(h,0,2)*H_(h,1,0)-H_(h,0,0)*H_(h,1,2));
    hi[6]=id*(H_(h,1,0)*H_(h,2,1)-H_(h,1,1)*H_(h,2,0));
    hi[7]=id*(H_(h,0,1)*H_(h,2,0)-H_(h,0,0)*H_(h,2,1));
    hi[8]=id*(H_(h,0,0)*H_(h,1,1)-H_(h,0,1)*H_(h,1,0));
}


/* ═══════════════ 最小像規約 ═══════════════ */
/*  2粒子間の「最短」距離ベクトルを求める (周期境界条件)。
    分率座標で [-0.5, 0.5) に丸めてから実座標に戻す。
    forces() の最内側で呼ばれるホット関数。
    OpenACC化: #pragma acc routine seq でデバイス関数にする              */
static inline void mimg_flat(double &dx, double &dy, double &dz,
                             const double* hi, const double* h){
    double s0=hi[0]*dx+hi[1]*dy+hi[2]*dz;
    double s1=hi[3]*dx+hi[4]*dy+hi[5]*dz;
    double s2=hi[6]*dx+hi[7]*dy+hi[8]*dz;
    s0-=round(s0); s1-=round(s1); s2-=round(s2);
    dx=h[0]*s0+h[1]*s1+h[2]*s2;
    dy=h[3]*s0+h[4]*s1+h[5]*s2;
    dz=h[6]*s0+h[7]*s1+h[8]*s2;
}


/* ═══════════════ 四元数操作 ═══════════════ */
/*  四元数 q = [w,x,y,z] で3次元回転を表現。ジンバルロックなし、正規化だけで済む。
    - q2R_flat:      四元数→回転行列 (分子内原子の実空間座標の計算に使用)
    - qmul_flat:     四元数の積 (2つの回転を合成)
    - qnorm_flat:    正規化 |q|=1 (数値誤差の蓄積を修正)
    - omega2dq_flat: 角速度ω × 時間dt → 微小回転の四元数
    OpenACC化: 各関数の前に #pragma acc routine seq を追加               */

static inline void q2R_flat(const double* q, double* R){
    double w=q[0],x=q[1],y=q[2],z=q[3];
    R[0]=1-2*(y*y+z*z); R[1]=2*(x*y-w*z);   R[2]=2*(x*z+w*y);
    R[3]=2*(x*y+w*z);   R[4]=1-2*(x*x+z*z); R[5]=2*(y*z-w*x);
    R[6]=2*(x*z-w*y);   R[7]=2*(y*z+w*x);   R[8]=1-2*(x*x+y*y);
}

static inline void qmul_flat(const double* a, const double* b, double* out){
    out[0]=a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
    out[1]=a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2];
    out[2]=a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1];
    out[3]=a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0];
}

static inline void qnorm_flat(double* q){
    double n=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
    double inv=1.0/n;
    q[0]*=inv; q[1]*=inv; q[2]*=inv; q[3]*=inv;
}

static inline void omega2dq_flat(double wx, double wy, double wz,
                                 double dt, double* dq){
    double wm=sqrt(wx*wx+wy*wy+wz*wz);
    double th=wm*dt*0.5;
    if(th<1e-14){
        dq[0]=1.0; dq[1]=0.5*dt*wx; dq[2]=0.5*dt*wy; dq[3]=0.5*dt*wz;
    } else {
        double s=sin(th)/wm;
        dq[0]=cos(th); dq[1]=s*wx; dq[2]=s*wy; dq[3]=s*wz;
    }
}


/* ═══════════════ C60座標の生成 ═══════════════ */
/*  黄金比 phi = (1+sqrt(5))/2 を使い、3群で60頂点を生成:
      群1: (0, ±1, ±3*phi) の巡回置換 → 12頂点
      群2: (±2, ±(1+2*phi), ±phi) の巡回置換 → 24頂点
      群3: (±1, ±(2+phi), ±2*phi) の巡回置換 → 24頂点
    0.72倍スケールで実際のC60サイズ (半径≒3.55 A) に。
    I0 = 等方慣性モーメント: 剛体の回転方程式 tau = I * alpha で使用。       */
struct C60Data { double coords[60*3]; double I0, Mmol, Rmol; };

static C60Data generate_c60(){
    C60Data d; d.Mmol=MC60;
    double phi=(1.0+sqrt(5.0))/2.0;
    int n=0, cyc[3][3]={{0,1,2},{1,2,0},{2,0,1}};
    double tmp[60][3]={};
    for(int p=0;p<3;p++) for(int s2:{-1,1}) for(int s3:{-1,1}){
        tmp[n][cyc[p][1]]=s2; tmp[n][cyc[p][2]]=s3*3*phi; n++; }
    for(int p=0;p<3;p++) for(int s1:{-1,1}) for(int s2:{-1,1}) for(int s3:{-1,1}){
        tmp[n][cyc[p][0]]=s1*2; tmp[n][cyc[p][1]]=s2*(1+2*phi); tmp[n][cyc[p][2]]=s3*phi; n++; }
    for(int p=0;p<3;p++) for(int s1:{-1,1}) for(int s2:{-1,1}) for(int s3:{-1,1}){
        tmp[n][cyc[p][0]]=s1; tmp[n][cyc[p][1]]=s2*(2+phi); tmp[n][cyc[p][2]]=s3*2*phi; n++; }
    double cm[3]={0,0,0};
    for(int i=0;i<60;i++) for(int a=0;a<3;a++) cm[a]+=tmp[i][a];
    for(int a=0;a<3;a++) cm[a]/=60.0;
    for(int i=0;i<60;i++) for(int a=0;a<3;a++) tmp[i][a]=(tmp[i][a]-cm[a])*0.72;
    d.Rmol=0; double Isum=0;
    for(int i=0;i<60;i++){
        double r2=tmp[i][0]*tmp[i][0]+tmp[i][1]*tmp[i][1]+tmp[i][2]*tmp[i][2];
        double r=sqrt(r2); if(r>d.Rmol) d.Rmol=r;
        Isum+=mC*r2;
    }
    d.I0=Isum*2.0/3.0;
    for(int i=0;i<60;i++) for(int a=0;a<3;a++) d.coords[i*3+a]=tmp[i][a];
    return d;
}


/* ═══════════════ FCC結晶生成 ═══════════════ */
/*  FCC (面心立方格子) は C60 結晶の室温での安定構造。
    1単位格子に4分子: (0,0,0), (a/2,a/2,0), (a/2,0,a/2), (0,a/2,a/2)
    nc x nc x nc 格子 → 4 * nc^3 分子 (nc=3 → 108分子)                     */
static int make_fcc(double a, int nc, double* pos, double* h){
    double bas[4][3]={{0,0,0},{.5*a,.5*a,0},{.5*a,0,.5*a},{0,.5*a,.5*a}};
    int n=0;
    for(int ix=0;ix<nc;ix++) for(int iy=0;iy<nc;iy++) for(int iz=0;iz<nc;iz++)
        for(int b=0;b<4;b++){
            pos[n*3+0]=a*ix+bas[b][0];
            pos[n*3+1]=a*iy+bas[b][1];
            pos[n*3+2]=a*iz+bas[b][2];
            n++;
        }
    for(int i=0;i<9;i++) h[i]=0.0;
    h[0]=h[4]=h[8]=nc*a;
    return n;
}


/* ═══════════════ 近傍リスト構築 (対称フルリスト) ═══════════════ */
/*  対称リスト: i→j と j→i の両方を格納。
    forces() で各分子が自分への力のみ集計 → 書込み競合なし (GPU並列化に最適)。
    nl_count[i]: 分子iの近傍数
    nl_list[i*MAX_NEIGH+k]: 分子iのk番目の近傍インデックス                   */
static void nlist_build_sym(const double* pos, const double* h, const double* hi,
                            int N, double rmcut, int* nl_count, int* nl_list){
    double rc2=(rmcut+3.0)*(rmcut+3.0);
    for(int i=0;i<N;i++) nl_count[i]=0;

    for(int i=0;i<N;i++){
        for(int j=i+1;j<N;j++){
            double dx=pos[j*3]-pos[i*3];
            double dy=pos[j*3+1]-pos[i*3+1];
            double dz=pos[j*3+2]-pos[i*3+2];
            mimg_flat(dx,dy,dz,hi,h);
            double r2=dx*dx+dy*dy+dz*dz;
            if(r2<rc2){
                int ci=nl_count[i], cj=nl_count[j];
                if(ci<MAX_NEIGH){ nl_list[i*MAX_NEIGH+ci]=j; nl_count[i]++; }
                else { printf("WARNING: nl overflow mol %d (count=%d)\n",i,ci); }
                if(cj<MAX_NEIGH){ nl_list[j*MAX_NEIGH+cj]=i; nl_count[j]++; }
                else { printf("WARNING: nl overflow mol %d (count=%d)\n",j,cj); }
            }
        }
    }
}


/* ═══════════════ PBC適用 ═══════════════ */
/*  全分子をセル内 (分率座標 s ∈ [0,1)) に戻す。
    各分子は独立に処理可能。
    OpenMP化:  #pragma omp parallel for schedule(static)
    OpenACC化: #pragma acc parallel loop present(pos,h,hi)               */
static void apply_pbc(double* pos, const double* h, const double* hi, int N){
    for(int i=0;i<N;i++){
        double px=pos[i*3], py=pos[i*3+1], pz=pos[i*3+2];
        double s0=hi[0]*px+hi[1]*py+hi[2]*pz;
        double s1=hi[3]*px+hi[4]*py+hi[5]*pz;
        double s2=hi[6]*px+hi[7]*py+hi[8]*pz;
        s0-=floor(s0); s1-=floor(s1); s2-=floor(s2);
        pos[i*3  ]=h[0]*s0+h[1]*s1+h[2]*s2;
        pos[i*3+1]=h[3]*s0+h[4]*s1+h[5]*s2;
        pos[i*3+2]=h[6]*s0+h[7]*s1+h[8]*s2;
    }
}


/* ═══════════════ 力・トルク・ビリアル計算 (メインカーネル) ═══════════════ */
/*  ★ 全計算時間の90%以上を占める最重量計算 — 並列化の最優先ターゲット

    処理の流れ:
      1. lab座標計算: 四元数→回転行列→分子内原子の実空間座標(重心基準)
      2. メイン力ループ: 近傍リスト内の全ペア(i,j)について60×60原子ペアのLJ力を計算

    対称フルリスト方式:
      各分子iは自分への力のみ集計 → 書込み競合なし (GPU並列化に最適)
      Newton第3法則不使用 → 計算量2倍だが完全並列化可能
      エネルギー・ビリアルは0.5倍して重複カウントを補正

    並列化ガイド:
      OpenMP化:
        lab計算ループ:  #pragma omp parallel for schedule(static)
        メイン力ループ: #pragma omp parallel for schedule(dynamic,1) reduction(+:Ep)
        Wm9加算:        #pragma omp atomic update
      OpenACC化:
        lab計算ループ:  #pragma acc parallel loop gang vector_length(128) present(...)
        原子内ループ:   #pragma acc loop vector
        メイン力ループ: #pragma acc parallel loop gang vector_length(128) present(...) reduction(+:Ep)
        原子ペアループ: #pragma acc loop vector reduction(+:fi0,...,w22)
        Wm9加算:        #pragma acc atomic update                         */
static double forces(double* Fv, double* Tv, double* Wm9,
                     const double* pos, const double* qv,
                     const double* body, const double* h, const double* hi,
                     const int* nl_count, const int* nl_list,
                     int N, int natom, double rmcut2,
                     double* lab)
{
    /* --- ステップ1: lab座標計算 (四元数→回転行列→分子内原子の実空間座標) --- */
    for(int i=0;i<N;i++){
        double R[9];
        q2R_flat(&qv[i*4], R);
        for(int a=0;a<natom;a++){
            double bx=body[a*3], by=body[a*3+1], bz=body[a*3+2];
            int idx=i*natom*3+a*3;
            lab[idx  ]=R[0]*bx+R[1]*by+R[2]*bz;
            lab[idx+1]=R[3]*bx+R[4]*by+R[5]*bz;
            lab[idx+2]=R[6]*bx+R[7]*by+R[8]*bz;
        }
    }

    /* --- ゼロ初期化 --- */
    for(int i=0;i<N*3;i++){ Fv[i]=0.0; Tv[i]=0.0; }
    for(int i=0;i<9;i++) Wm9[i]=0.0;

    /* --- ステップ2: LJ力計算メインカーネル --- */
    double Ep=0.0;

    for(int i=0;i<N;i++){
        double fi0=0, fi1=0, fi2=0;   /* 分子iへの力の局所蓄積 */
        double ti0=0, ti1=0, ti2=0;   /* 分子iへのトルクの局所蓄積 */
        double my_Ep=0;
        /* ビリアル9成分をスカラーに展開 (GPU vector reduction互換のため配列不使用) */
        double w00=0,w01=0,w02=0, w10=0,w11=0,w12=0, w20=0,w21=0,w22=0;

        int nni=nl_count[i];
        for(int k=0;k<nni;k++){
            int j=nl_list[i*MAX_NEIGH+k];
            /* 分子ij間の最小像距離 */
            double dmx=pos[j*3]-pos[i*3];
            double dmy=pos[j*3+1]-pos[i*3+1];
            double dmz=pos[j*3+2]-pos[i*3+2];
            mimg_flat(dmx,dmy,dmz,hi,h);
            if(dmx*dmx+dmy*dmy+dmz*dmz > rmcut2) continue;

            /* 60×60 原子ペアループ (シングルスレッドでは逐次実行) */
            for(int ai=0;ai<natom;ai++){
                int ia=i*natom*3+ai*3;
                double rax=lab[ia], ray=lab[ia+1], raz=lab[ia+2];

                for(int bj=0;bj<natom;bj++){
                    int jb=j*natom*3+bj*3;
                    double rbx=lab[jb], rby=lab[jb+1], rbz=lab[jb+2];
                    double ddx=dmx+rbx-rax;
                    double ddy=dmy+rby-ray;
                    double ddz=dmz+rbz-raz;
                    double r2=ddx*ddx+ddy*ddy+ddz*ddz;

                    if(r2<RCUT2){
                        if(r2<0.25) r2=0.25;   /* NaN防止クランプ */
                        double ri2=1.0/r2;
                        double sr2=sig2_LJ*ri2;
                        double sr6=sr2*sr2*sr2;
                        double sr12=sr6*sr6;
                        double fm=24.0*eps_LJ*(2.0*sr12-sr6)*ri2;
                        double fx=fm*ddx, fy=fm*ddy, fz=fm*ddz;

                        /* 分子iへの力 (j側はjのループで計算される) */
                        fi0-=fx; fi1-=fy; fi2-=fz;
                        /* トルク = -ra × F */
                        ti0-=(ray*fz-raz*fy);
                        ti1-=(raz*fx-rax*fz);
                        ti2-=(rax*fy-ray*fx);
                        /* エネルギー (半分: 対称リストで(i,j)と(j,i)の両方で計算される) */
                        my_Ep+=0.5*(4.0*eps_LJ*(sr12-sr6)-VSHFT);
                        /* ビリアル (半分) */
                        w00+=0.5*ddx*fx; w01+=0.5*ddx*fy; w02+=0.5*ddx*fz;
                        w10+=0.5*ddy*fx; w11+=0.5*ddy*fy; w12+=0.5*ddy*fz;
                        w20+=0.5*ddz*fx; w21+=0.5*ddz*fy; w22+=0.5*ddz*fz;
                    }
                }
            }
        }

        /* 分子iへの力・トルクをグローバル配列に書込み (競合なし: 各iが自分だけ書く) */
        Fv[i*3]=fi0; Fv[i*3+1]=fi1; Fv[i*3+2]=fi2;
        Tv[i*3]=ti0; Tv[i*3+1]=ti1; Tv[i*3+2]=ti2;

        /* ビリアルは全分子からの寄与を集約
           並列化時は atomic update で保護する:
             OpenMP化:  #pragma omp atomic update
             OpenACC化: #pragma acc atomic update                        */
        Wm9[0]+=w00; Wm9[1]+=w01; Wm9[2]+=w02;
        Wm9[3]+=w10; Wm9[4]+=w11; Wm9[5]+=w12;
        Wm9[6]+=w20; Wm9[7]+=w21; Wm9[8]+=w22;

        Ep+=my_Ep;
    }
    return Ep;
}


/* ═══════════════ 運動エネルギー ═══════════════ */
/*  並列化ガイド:
      OpenMP化:  #pragma omp parallel for reduction(+:s)
      OpenACC化: #pragma acc parallel loop present(vel) reduction(+:s)   */
static double ke_trans(const double* vel, int N, double Mmol){
    double s=0;
    for(int i=0;i<N;i++)
        s+=vel[i*3]*vel[i*3]+vel[i*3+1]*vel[i*3+1]+vel[i*3+2]*vel[i*3+2];
    return 0.5*Mmol*s/CONV;
}

static double ke_rot(const double* omg, int N, double I0){
    double s=0;
    for(int i=0;i<N;i++)
        s+=omg[i*3]*omg[i*3]+omg[i*3+1]*omg[i*3+1]+omg[i*3+2]*omg[i*3+2];
    return 0.5*I0*s/CONV;
}

static inline double inst_T(double KE, int Nf){return 2*KE/(Nf*kB);}
static inline double inst_P(const double* W, double KEt, double V){
    return (2*KEt+W[0]+W[4]+W[8])/(3*V)*eV2GPa;
}


/* ═══════════════ NPT状態変数 ═══════════════ */
/*  Nose-Hoover サーモスタット + Parrinello-Rahman バロスタット:
      xi:    サーモスタット変数 (温度制御の摩擦係数)
      Q:     サーモスタット質量 (温度制御の応答速度)
      Vg[9]: バロスタット速度 (セル変形の速度)
      W:     バロスタット質量 (圧力制御の応答速度)
      Pe:    目標圧力 [GPa], Tt: 目標温度 [K], Nf: 自由度              */
struct NPTState { double xi,Q,Vg[9],W,Pe,Tt; int Nf; };

static NPTState make_npt(double T, double Pe, int N){
    int Nf=6*N-3;  /* 並進3自由度 + 回転3自由度 - 重心拘束3 */
    NPTState s;
    s.xi=0; s.Q=std::max(Nf*kB*T*100.0*100.0,1e-20);
    for(int i=0;i<9;i++) s.Vg[i]=0;
    s.W=std::max((Nf+9)*kB*T*1000.0*1000.0,1e-20);
    s.Pe=Pe; s.Tt=T; s.Nf=Nf;
    return s;
}


/* ═══════════════ NPT速度ベルレ1ステップ ═══════════════ */
/*  NPT (一定温度・一定圧力) アンサンブルでの時間発展。
    Velocity-Verlet法: 力の計算を1ステップに1回だけで済ませる。
      (A) サーモスタット前半 → (B) バロスタット前半 → (C) 速度前半更新
      → (D) 座標更新 → (E) セル更新 → (F) 分率→実座標 → (G) 四元数更新
      → (H) 力の再計算 → (I) 速度後半更新 → (J)(K) サーモ/バロ後半

    並列化ガイド:
      各ループ (C)(D)(F)(G)(I) は全分子独立:
        OpenMP化:  #pragma omp parallel for schedule(static)
        OpenACC化: #pragma acc parallel loop present(...)
      (E) のH行列更新はホスト側のスカラー操作 (並列化不要)
      OpenACC化ではH行列更新後に #pragma acc update device(h[0:9]) が必要 */
static std::pair<double,double>
step_npt(double* pos, double* vel, double* qv, double* omg,
         double* Fv, double* Tv, double* Wm9,
         double* h, double* hi,
         const double* body, double I0, double Mmol,
         int N, int natom, double rmcut2, double dt, NPTState& npt,
         const int* nl_count, const int* nl_list, double* lab)
{
    double hdt=0.5*dt;
    mat_inv9(h,hi);

    double V=fabs(mat_det9(h));
    double kt=ke_trans(vel,N,Mmol), kr=ke_rot(omg,N,I0), KE=kt+kr;

    /* (A) サーモスタット前半 */
    npt.xi+=hdt*(2*KE-npt.Nf*kB*npt.Tt)/npt.Q;
    npt.xi=std::clamp(npt.xi,-0.1,0.1);

    /* (B) バロスタット前半 */
    double dP=inst_P(Wm9,kt,V)-npt.Pe;
    for(int a=0;a<3;a++) npt.Vg[a*4]+=hdt*V*dP/(npt.W*eV2GPa);
    for(int a=0;a<3;a++) npt.Vg[a*4]=std::clamp(npt.Vg[a*4],-0.01,0.01);

    double eps_tr=npt.Vg[0]*hi[0]+npt.Vg[4]*hi[4]+npt.Vg[8]*hi[8];
    double sc_nh=exp(-hdt*npt.xi);
    double sc_pr=exp(-hdt*eps_tr/3.0);
    double sc_v=sc_nh*sc_pr;
    double cF=CONV/Mmol, cT=CONV/I0;

    /* (C) 速度の前半更新 */
    for(int i=0;i<N;i++){
        for(int a=0;a<3;a++){
            vel[i*3+a]=vel[i*3+a]*sc_v+hdt*Fv[i*3+a]*cF;
            omg[i*3+a]=omg[i*3+a]*sc_nh+hdt*Tv[i*3+a]*cT;
        }
    }

    /* (D) 座標更新 (分率座標で積分 + PBC) */
    for(int i=0;i<N;i++){
        double px=pos[i*3],py=pos[i*3+1],pz=pos[i*3+2];
        double vx=vel[i*3],vy=vel[i*3+1],vz=vel[i*3+2];
        double sx=hi[0]*px+hi[1]*py+hi[2]*pz;
        double sy=hi[3]*px+hi[4]*py+hi[5]*pz;
        double sz=hi[6]*px+hi[7]*py+hi[8]*pz;
        double vsx=hi[0]*vx+hi[1]*vy+hi[2]*vz;
        double vsy=hi[3]*vx+hi[4]*vy+hi[5]*vz;
        double vsz=hi[6]*vx+hi[7]*vy+hi[8]*vz;
        sx+=dt*vsx; sy+=dt*vsy; sz+=dt*vsz;
        sx-=floor(sx); sy-=floor(sy); sz-=floor(sz);
        pos[i*3]=sx; pos[i*3+1]=sy; pos[i*3+2]=sz;
    }

    /* (E) セルH行列の更新 */
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) h[a*3+b]+=dt*npt.Vg[a*3+b];

    /* (F) 分率座標→実座標 */
    for(int i=0;i<N;i++){
        double sx=pos[i*3],sy=pos[i*3+1],sz=pos[i*3+2];
        pos[i*3  ]=h[0]*sx+h[1]*sy+h[2]*sz;
        pos[i*3+1]=h[3]*sx+h[4]*sy+h[5]*sz;
        pos[i*3+2]=h[6]*sx+h[7]*sy+h[8]*sz;
    }

    /* (G) 四元数更新 */
    for(int i=0;i<N;i++){
        double dq[4],tmp[4];
        omega2dq_flat(omg[i*3],omg[i*3+1],omg[i*3+2],dt,dq);
        qmul_flat(&qv[i*4],dq,tmp);
        qv[i*4]=tmp[0]; qv[i*4+1]=tmp[1]; qv[i*4+2]=tmp[2]; qv[i*4+3]=tmp[3];
        qnorm_flat(&qv[i*4]);
    }

    /* (H) 力の再計算 */
    mat_inv9(h,hi);
    double Ep=forces(Fv,Tv,Wm9,pos,qv,body,h,hi,nl_count,nl_list,
                     N,natom,rmcut2,lab);

    /* (I) 速度の後半更新 */
    double eps_tr2=npt.Vg[0]*hi[0]+npt.Vg[4]*hi[4]+npt.Vg[8]*hi[8];
    double sc_v2=sc_nh*exp(-hdt*eps_tr2/3.0);
    for(int i=0;i<N;i++){
        for(int a=0;a<3;a++){
            vel[i*3+a]=(vel[i*3+a]+hdt*Fv[i*3+a]*cF)*sc_v2;
            omg[i*3+a]=(omg[i*3+a]+hdt*Tv[i*3+a]*cT)*sc_nh;
        }
    }

    /* (J)(K) サーモ/バロ後半更新 */
    kt=ke_trans(vel,N,Mmol); kr=ke_rot(omg,N,I0); KE=kt+kr;
    npt.xi+=hdt*(2*KE-npt.Nf*kB*npt.Tt)/npt.Q;
    npt.xi=std::clamp(npt.xi,-0.1,0.1);
    double V2=fabs(mat_det9(h));
    dP=inst_P(Wm9,kt,V2)-npt.Pe;
    for(int a=0;a<3;a++) npt.Vg[a*4]+=hdt*V2*dP/(npt.W*eV2GPa);
    for(int a=0;a<3;a++) npt.Vg[a*4]=std::clamp(npt.Vg[a*4],-0.01,0.01);

    return {Ep,KE};
}


/* ═══════════════ メインプログラム ═══════════════ */
int main(int argc, char** argv){
    /* --- パラメータ設定 --- */
    int nc=3;
    for(int i=1;i<argc;i++){
        std::string a=argv[i];
        if(a.substr(0,7)=="--cell=") nc=std::atoi(a.substr(7).c_str());
        else if(a[0]!='-') nc=atoi(argv[i]);
    }
    if(nc<1||nc>8){ printf("Error: nc must be 1-8 (got %d)\n",nc); return 1; }

    constexpr int nsteps=1000;
    constexpr int mon=100;
    constexpr int nlup=25;
    constexpr double T=300.0;
    constexpr double Pe=0.0;
    constexpr double dt=1.0;
    constexpr double a0=14.17;
    int avg_from=nsteps-nsteps/4;

    /* --- C60分子座標の生成 --- */
    C60Data c60=generate_c60();
    int natom=C60_NATOM;

    /* --- 配列確保 (生配列: OpenACC data region と互換性のある構造) --- */
    int Nmax=4*nc*nc*nc;
    double* pos      =new double[Nmax*3]();
    double* vel      =new double[Nmax*3]();
    double* omg      =new double[Nmax*3]();
    double* qv       =new double[Nmax*4]();
    double* Fv       =new double[Nmax*3]();
    double* Tv       =new double[Nmax*3]();
    double* lab      =new double[Nmax*natom*3]();
    double  h[9]={}, hi[9]={}, Wm9[9]={};
    double* body     =new double[natom*3];
    int*    nl_count =new int[Nmax]();
    int*    nl_list  =new int[Nmax*MAX_NEIGH]();

    for(int i=0;i<natom*3;i++) body[i]=c60.coords[i];

    int N=make_fcc(a0,nc,pos,h);
    mat_inv9(h,hi);

    /* --- バナー表示 --- */
    printf("================================================================\n");
    printf("  C60 LJ NPT-MD Core (Serial)\n");
    printf("================================================================\n");
    printf("  FCC cell        : %dx%dx%d  N=%d molecules\n", nc,nc,nc,N);
    printf("  Atoms/molecule  : %d\n", natom);
    printf("  a0=%.2f A  T=%.0f K  P=%.1f GPa  dt=%.1f fs  steps=%d\n",
           a0,T,Pe,dt,nsteps);
    printf("  MAX_NEIGH=%d  VECTOR_LENGTH=%d\n", MAX_NEIGH, VECTOR_LENGTH);

    double mem_lab = (double)N*natom*3*8/1024/1024;
    double mem_nl  = (double)N*MAX_NEIGH*4/1024/1024;
    double mem_total = (double)(N*3*8*5 + N*4*8 + N*natom*3*8 + natom*3*8
                                + N*4 + N*MAX_NEIGH*4 + 9*8*3) / 1024/1024;
    printf("  Memory          : lab=%.2f MB  nl=%.2f MB  total=%.2f MB\n",
           mem_lab, mem_nl, mem_total);
    printf("================================================================\n\n");

    /* --- 初期速度 (Maxwell-Boltzmann分布) --- */
    std::mt19937 rng(42);
    std::normal_distribution<double> gauss(0,1);
    double sv=sqrt(kB*T*CONV/c60.Mmol), sw=sqrt(kB*T*CONV/c60.I0);
    for(int i=0;i<N;i++){
        for(int a=0;a<3;a++){vel[i*3+a]=sv*gauss(rng); omg[i*3+a]=sw*gauss(rng);}
        for(int a=0;a<4;a++) qv[i*4+a]=gauss(rng);
        double n=sqrt(qv[i*4]*qv[i*4]+qv[i*4+1]*qv[i*4+1]
                     +qv[i*4+2]*qv[i*4+2]+qv[i*4+3]*qv[i*4+3]);
        for(int a=0;a<4;a++) qv[i*4+a]/=n;
    }
    /* 重心速度除去 */
    double vcm[3]={0,0,0};
    for(int i=0;i<N;i++) for(int a=0;a<3;a++) vcm[a]+=vel[i*3+a];
    for(int a=0;a<3;a++) vcm[a]/=N;
    for(int i=0;i<N;i++) for(int a=0;a<3;a++) vel[i*3+a]-=vcm[a];

    NPTState npt=make_npt(T,Pe,N);

    /* --- 初回の近傍リスト構築 --- */
    nlist_build_sym(pos,h,hi,N,RMCUT,nl_count,nl_list);

    /* PBC適用 + 初回力計算 */
    /*  OpenACC化: この前に #pragma acc data copy(...) copyin(...) create(...)
        でGPUメモリに配列を配置し、MDループ全体をデータ領域内で実行する      */
    apply_pbc(pos,h,hi,N);
    forces(Fv,Tv,Wm9,pos,qv,body,h,hi,nl_count,nl_list,
           N,natom,RMCUT2,lab);

    double sT=0,sP=0,sa=0,sEp=0; int nav=0;
    auto t0=std::chrono::steady_clock::now();
    printf("%8s %7s %9s %8s %10s %7s\n",
           "step","T[K]","P[GPa]","a[A]","Ecoh[eV]","t[s]");

    /* ═══ MDメインループ ═══ */
    for(int g=1;g<=nsteps;g++){

        /* 近傍リスト再構築 */
        /*  OpenACC化: 再構築前に #pragma acc update self(pos[0:N*3]) で
            座標をホストに転送し、構築後に
            #pragma acc update device(nl_count[0:N], nl_list[0:N*MAX_NEIGH])
            でデバイスに転送する                                             */
        if(g%nlup==0){
            mat_inv9(h,hi);
            nlist_build_sym(pos,h,hi,N,RMCUT,nl_count,nl_list);
        }

        /* 1ステップの時間発展 */
        auto [Ep,KE]=step_npt(pos,vel,qv,omg,Fv,Tv,Wm9,
                               h,hi,body,c60.I0,c60.Mmol,
                               N,natom,RMCUT2,dt,npt,
                               nl_count,nl_list,lab);

        /* 瞬時物理量 */
        double kt=ke_trans(vel,N,c60.Mmol);
        double V=fabs(mat_det9(h));
        double Tn=inst_T(KE,npt.Nf), Pn=inst_P(Wm9,kt,V);
        double Ec=Ep/N, an=h[0]/nc;
        if(g>=avg_from){sT+=Tn;sP+=Pn;sa+=an;sEp+=Ec;nav++;}

        /* モニタリング出力 */
        if(g%mon==0||g==nsteps){
            double el=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-t0).count();
            printf("%8d %7.1f %9.3f %8.3f %10.5f %7.0f\n",
                   g,Tn,Pn,an,Ec,el);
        }
    }

    if(nav>0) printf("Avg(%d): T=%.2f P=%.4f a=%.4f Ecoh=%.5f\n",
                      nav,sT/nav,sP/nav,sa/nav,sEp/nav);
    printf("Done %.1fs\n",std::chrono::duration<double>(
        std::chrono::steady_clock::now()-t0).count());

    /* --- 解放 --- */
    delete[] pos; delete[] vel; delete[] omg; delete[] qv;
    delete[] Fv;  delete[] Tv;  delete[] lab; delete[] body;
    delete[] nl_count; delete[] nl_list;
    return 0;
}
