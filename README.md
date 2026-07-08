# bio-d: E-Cell4 practice

[E-Cell4](https://ecell4.e-cell.org/) を使った細胞シミュレーションの練習リポジトリ。
決定論（ODE）と確率論（Gillespie）を中心に、反応ネットワークのモデル化を段階的に手を動かして学ぶ。

## 環境

- Python 3.14 / macOS arm64（`ecell4_base` の公式 wheel が対応）
- パッケージ管理は [`uv`](https://docs.astral.sh/uv/)

```bash
git clone https://github.com/Carlo-Broschi/bio-d_ecell-practice.git
cd bio-d_ecell-practice
uv sync                       # .venv を作成し依存を同期
uv run jupyter lab            # Notebook で練習する場合
```

スクリプトを直接実行する場合:

```bash
uv run python scripts/01_reversible_reaction.py   # outputs/ に PNG を保存
```

## 構成

| パス | 内容 |
|---|---|
| `notebooks/01_getting_started.ipynb` | 可逆反応 A ⇌ B を ODE と Gillespie で解く最初の練習 |
| `notebooks/02_srna_silencing.ipynb`  | Hfq 依存 sRNA による標的 mRNA 抑制の最小回路（`bio-a` と地続き） |
| `notebooks/03_translation_readout.ipynb` | 02 に翻訳段 `m→m+P` を追加。sRNA 制御のタンパク質出力への伝播とノイズ濾過 |
| `notebooks/04_catalytic_vs_stoichiometric.ipynb` | 化学量論 `m+s→∅` vs 触媒 `m+s→s`。応答が閾値スイッチ / MM 的滴定に分岐 |
| `notebooks/05_consumption_vs_recycling.ipynb` | 04 の批判的検算。速度を揃えると「触媒=graded」は崩れ、閾値の正体は sRNA 消費だと分かる |
| `notebooks/06_hfq_robustness.ipynb` | **論文再現①** Adamson & Lim 2011。Hfq を陽に入れ、duplex 形成が Hfq に対し釣鐘状（不足／隔離）＋解離が窓を広げる |
| `notebooks/07_hfq_setpoint.ipynb` | **論文再現②** Sagawa & Lim 2015。random-order 結合だと高 Hfq で逆に抑制（Hfq set-point）、compulsory-order なら単調飽和 |
| `notebooks/08_srna_competition.ipynb` | **論文再現③** Moon & Gottesman 2011。2 sRNA が Hfq を共有→競合者の過剰発現でレポーターsRNA抑制、Hfq追加で回復 |
| `notebooks/09_partnered_vs_unpartnered.ipynb` | **論文再現④** Hussein & Lim 2011。unpartnered な sRNA/mRNA 単独は Hfq を抱え妨害、partnered（共転写）なら duplex で Hfq を返し妨害小 |
| `notebooks/10_cotranscriptional.ipynb` | **論文再現⑤** Reyer 2021。実測パラメータで post-転写だけでは抑制が浅く、共転写制御(α_m→0.46α_m)を足すと実測の深い抑制に合う |
| `notebooks/11_active_exchange.ipynb` | **論文再現⑥** Park/Fei 2021。Hfqは常時ほぼmRNA占有＋強固な結合。能動交換(sRNAがmRNAを置換)で分単位の応答が可能に |
| `notebooks/12_active_cycling.ipynb` | **論文再現⑦** Wagner 2013。滞在RNAの見かけ解離速度が競合濃度で上がり飽和（能動循環の指紋、半減期>150分→~1.5分） |
| `notebooks/13_autoregulation.ipynb` | **論文再現⑧** Morita & Aiba 2019。Hfqが自分のmRNAを抑える負の自己制御→転写8×でも蛋白4×どまり(恒常性)＋応答3倍速 |
| `notebooks/14_fabmos.ipynb` | **論文再現⑨** Hung 2014 FABMOS（Gillespie）。平均一定でスイッチ頻度を変えると分布が bimodal↔unimodal、CVだけを制御 |
| `scripts/01_reversible_reaction.py`  | 01 のスクリプト版（ヘッドレスで図を保存・動作検証用） |
| `scripts/02_srna_silencing.py`       | 02 のスクリプト版（応答曲線 + Fano ノイズ図を保存） |
| `scripts/03_translation_readout.py`  | 03 のスクリプト版 |
| `scripts/04_catalytic_vs_stoichiometric.py` | 04 のスクリプト版 |
| `scripts/05_consumption_vs_recycling.py` | 05 のスクリプト版 |
| `outputs/` | 生成物（`.gitignore` 済み。再現可能なので追跡しない） |

### 03 / 04 のポイント

- **03 翻訳段**: `m→m+P`, `P→∅` を追加。タンパク質は `<P>=(k_p/b_p)<m>` で閾値応答を継承。長寿命タンパク質 (`b_p<b_m`) は mRNA ノイズを時間平均で濾過し `CV_P<CV_m` になるが、閾値近傍のノイズ増幅は減衰しつつ伝播する。
- **04 作用機構**: sRNA が消費される化学量論的共分解 (`m+s→∅`) は鋭い閾値スイッチ、sRNA が再利用される触媒的分解 (`m+s→s`) は `<m>=a_m/(b_m+k_cat·a_s/b_s)` の双曲線（MM 的）でゼロにならない。同じ「抑制」でもデジタル/アナログが分かれる。
- **05 検算**: 04 の「触媒=graded」は触媒側に弱い速度を与えた見かけだった。ペアリング率 `k_pair` を揃え分岐確率 `φ`（sRNA 再利用確率）でつなぐと、φ↑（触媒寄り）ほど silencing はむしろ**強く・低 a_s 側**に。閾値の正体は触媒/消費の別ではなく **sRNA が消費（滴定）されること**（+ 鋭さには強いペアリングも必要）。モデル比較は速度定数を揃えて再検算する、という作法の教材。

### 02 のポイント（研究との接点）

sRNA サイレンシングの最小回路（∅→m, ∅→s, m→∅, s→∅, m+s→∅）。sRNA 転写率 `a_s` を掃引すると、
`a_s > a_m` で標的 mRNA が閾値的に急落（*threshold-linear* 応答）し、その閾値近傍で Gillespie の
ゆらぎ（Fano factor）が Poisson の 1 から ~3 へ増幅する。Hfq / sRNA の**進化**を扱う `bio-a` に対する
「回路としての振る舞い」側の相棒。**掃引で毎回モデルを組み直すため、DSL ではなく `ReactionRule` の明示 API を使う**
（`reaction_rules()` はグローバルモデルに累積するので掃引に不向き）。

## E-Cell4 の最小 API メモ

```python
from ecell4 import get_model, reaction_rules, run_simulation

with reaction_rules():
    A > B | 0.25      # A -> B, 速度定数 0.25
    B > A | 1.0
model = get_model()

ret = run_simulation(10.0, y0={'A': 60, 'B': 0}, model=model,
                     solver='ode',        # 'ode' | 'gillespie' | 'meso' | ...
                     species_list=['A', 'B'])
ret.as_array()        # ndarray（各行 = [t, A, B]）
ret.as_dataframe()    # pandas DataFrame（要 pandas）
ret.plot()            # 組み込みプロット
```

- 可逆反応の平衡比は `[A]/[B] = kr/kf`。上の例なら `1.0/0.25 = 4`（60 分子 → A:B = 48:12）。
- `solver` を差し替えるだけで決定論／確率論を切り替えられるのが E-Cell4 の勘所。

## 参考

- E-Cell4 ドキュメント: https://ecell4.e-cell.org/
- プロジェクトページ: https://www.e-cell.org/projects/ecell4.html
