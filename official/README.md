# 公式トラック — E-Cell4 公式チュートリアル/例題の移植

このディレクトリは **E-Cell4 公式ドキュメント**（https://ecell4.e-cell.org/ ）の\
**Examples（例題モデル）を、出典リンクつきで忠実に移植**したもの。`../notebooks/`（bio-d オリジナルの自作ノート）とは区別する。

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

## 公式が用意している素材の全体像（参考）

公式は「解くべき練習問題」ではなく、**チュートリアル（how-to, 10 本）と例題モデル（Examples, 14 本）とテストモデル（5 本）**を提供している。\
本トラックは **Examples 14 本を全て**収録した。チュートリアル（Brief Tour / Build a Model / Rule-based / Spatial Gillespie / Spatiocyte 等）や\
テストモデル（Birth-Death, Homodimerization, Reversible, MSD）も、同じ要領で追加できる。
