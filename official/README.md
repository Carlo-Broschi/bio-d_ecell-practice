# 公式トラック — E-Cell4 公式チュートリアル/例題の移植

このディレクトリは **E-Cell4 公式ドキュメント**（https://ecell4.e-cell.org/ ）の\
**Examples（例題モデル）を、出典リンクつきで忠実に移植**したもの。`../notebooks/`（bio-x オリジナルの自作ノート）とは区別する。

各ノートは冒頭に**出典 URL** を明記し、公式コードを可能な限りそのまま使う。インストール版（ecell4 1.2.2 / NumPy 2.5）で\
走らせるための最小限の調整（重いソルバの差し替え・依存追加・doc 誤植の訂正・非互換モジュールの回避）は各ノートに明示する。

**公式 Examples 全 14 本を収録**（`example01`〜`example14` に対応）。

| ノート | 公式 | テーマ | 実行 |
|---|---|---|---|
| `o1_simple_equilibrium.ipynb` | example10 | 低レベル World/Simulator・binding/unbinding rule・ソルバ差し替え | ✅ |
| `o2_lotka_volterra.ipynb` | example07 | 捕食者-被食者。rate-law ODE + 素反応 Gillespie | ✅ |
| `o3_dual_phosphorylation.ipynb` | example03 | 二重リン酸化サイクル（MAPK 型）。デコレータ・酵素反応連結記法 | ✅ |
| `o4_attractors.ipynb` | example01 | カオス力学系 5 種（Lorenz 他）。汎用 ODE ソルバとしての E-Cell4 | ✅ |
| `o5_drosophila_clock.ipynb` | example02 | 概日時計（Goldbeter 1995）。負フィードバック振動 | ✅ |
| `o6_action_potential.ipynb` | example06 | Hodgkin-Huxley 神経活動電位 | ✅ |
| `o7_tyson_cell_cycle.ipynb` | example11 | Tyson1991 細胞周期振動 | ✅ |
| `o8_glycolysis_mca.ipynb` | example05 | 解糖系 + 代謝制御解析（MCA） | ⚠️ ODEのみ実行（MCAは `ecell4.mca` が NumPy<2 依存で非実行・コード掲載） |
| `o9_unit_system.ipynb` | example12 | 単位系（pint）。次元つきパラメータ | ✅（`pint` 追加要） |
| `o10_egfr_rulebased.ipynb` | example04 | ルールベース（BNGL 風）EGFR シグナル。site/state・expand | ✅ |
| `o11_minde_meso.ipynb` | example08 | MinDE 極間振動（メソスコピック空間） | ⛔ 重い空間（本環境で duration=30 が ~17分）・可視化インタラクティブ→**非自動実行・コード忠実掲載** |
| `o12_minde_spatiocyte.ipynb` | example09 | MinDE（Spatiocyte 単分子空間） | ⛔ 同上（さらに重い duration=240）→**非自動実行・コード忠実掲載** |
| `o13_gpcr.ipynb` | example13 | GPCR 多量体複合体モデル。N=0(TCM) を実行＋N=1..5 概説 | ✅（N=0 のみ。上位は概説） |
| `o14_sgfrd.ipynb` | example14 | sGFRD（曲面上の単分子反応拡散）の解説＋API 形 | ⚠️ 出典ページが大きくコード自動取得不可→解説＋API スケッチ（出典参照） |

## 実行と注記

- ✅ = `nbconvert --execute` で実行・出力あり。⚠️ = 一部のみ実行（理由を本文に明記）。⛔ = 重い空間系のため非自動実行（公式コードは忠実に掲載、手元で実行可）。
- 空間系（o11/o12）は `plotting.plot_world`/`plot_movie` などインタラクティブ可視化を使う。`uv run jupyter lab` から実行を。

## チュートリアル（`tutorials/`）— 公式 Tutorials 全 10 本

公式の how-to チュートリアルを実行可能な形で移植（すべて ✅ 実行済み）。

| ノート | 公式 | 内容 |
|---|---|---|
| `tut01_brief_tour` | tutorial01 | 最短ワークフロー（1 モデルを ode/gillespie/meso で） |
| `tut02_build_a_model` | tutorial02 | モデルの 4 通りの書き方（DSL/オブジェクト/ファクトリ/デコレータ） |
| `tut03_initial_condition` | tutorial03 | World・add_molecules・volume・bind_to |
| `tut04_run_a_simulation` | tutorial04 | Simulator・step/run・Factory パターン |
| `tut05_log_and_visualize` | tutorial05 | Observer（記録）と可視化 |
| `tut06_rate_law_odes` | tutorial06 | rate-law（MM・任意式）・Pyfunc 速度則 |
| `tut07_rule_based` | tutorial07 | サイト/状態・ワイルドカード・expand |
| `tut08_more_about_brief_tour` | tutorial08 | World/Simulator 分離・Real3 演算・ソルバ付替え |
| `tut09_spatial_gillespie` | tutorial09 | meso 反応拡散・初期分離の効果（小規模で実行） |
| `tut10_spatiocyte` | tutorial10 | Spatiocyte 単分子・軌跡・構造（小規模で実行） |

## テストモデル（`tests/`）— 公式 Test models 全 5 本

各ソルバの整合性検証用の最小モデル。ODE+Gillespie を実行、空間系（meso/spatiocyte/egfrd/bd）は忠実コードを参照掲載。

| ノート | 公式 | 内容 |
|---|---|---|
| `t1_birth_death` | Birth-Death | 生成/死滅のみ。定常=birth/death |
| `t2_homodimerization` | Homodimerization/Annihilation | A+A→∅。巨視的 kon vs 微視的 ka |
| `t3_reversible` | Reversible | A+B⇌C（非拡散律速→全ソルバ一致） |
| `t4_reversible_diffusion_limited` | Reversible (Diff-limited) | 拡散律速→ well-mixed と空間がずれる |
| `t5_msd` | MSD | 単分子の平均二乗変位 = 6Dt（空間・コードのみ） |

## 収録状況

公式が提供する **チュートリアル 10 / 例題 14 / テスト 5 = 全て**を収録した（"練習問題(演習)"は公式に存在しない）。\
空間・単分子で重いもの（o11/o12 MinDE、t5 MSD、各テストの spatiocyte/egfrd/bd 呼び出し）は\
公式コードを忠実に掲載しつつ自動実行は避け、理由と実行方法を各ノートに明記している。
