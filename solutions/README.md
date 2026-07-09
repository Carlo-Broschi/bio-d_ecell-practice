# 発展課題（自作） 解答（solutions）

各 notebook（`../notebooks/NN_*.ipynb`）末尾の「発展課題（自作）」を解いたもの。
計算問題はコード＋図＋数値で解答し、`[[bio-a_hfq-phylogenetics]]` 接続の考察問題は
文章で答えつつ、定量化は spin-off プロジェクト **bio-d** で行う旨を示している。

すべて `uv run jupyter nbconvert --execute` で実行検証済み。

| 解答 | 対象 | 主な内容 |
|---|---|---|
| `01_solutions.ipynb` | 01 | 平衡比 kr/kf、少数分子ノイズ、ode/gillespie/meso 比較、2分子反応 |
| `02_solutions.ipynb` | 02 | ペアリング k と閾値の鋭さ、翻訳のノイズ濾過、触媒的 sRNA の MM 応答 |
| `03_solutions.ipynb` | 03 | b_p とフィルタ効果、翻訳バーストと Fano、CV vs Fano の使い分け |
| `04_solutions.ipynb` | 04 | 触媒半減点 b_m·b_s/k_cat、両モードのノイズ、混在の中間応答 |
| `05_solutions.ipynb` | 05 | k_pair と鋭い閾値の復活、閾値位置＝a_m、φ の系統形質解釈 |
| `06_solutions.ipynb` | 06 | K5 と釣鐘、非対称転写、協同結合で窓拡大 |
| `07_solutions.ipynb` | 07 | ka と set-point 位置、mRNA-first も飽和、**mRNA set-point** |
| `08_solutions.ipynb` | 08 | 強い競合者(ka)、親和性ヒエラルキー(3 sRNA)、標的量と感受性 |
| `09_solutions.ipynb` | 09 | k5 遅延で partnered→unpartnered 化、ka と差、活性化型の考察 |
| `10_solutions.ipynb` | 10 | P_CO と抑制深さ、K_XS=0 で翻訳/分解分離、sRNA 動的誘導の遅れ |
| `11_solutions.ipynb` | 11 | k_ex と応答時間、強固結合で差拡大、**RNase E デコイ**最小モデル |
| `12_solutions.ipynb` | 12 | K_REAR と飽和値、good/poor competitor、cycling×set-point 考察 |
| `13_solutions.ipynb` | 13 | 協同で緩衝は強まらない（隔離型の限界）、外乱の整定、競合への組込み |
| `14_solutions.ipynb` | 14 | bias で平均・頻度でばらつき（独立）、f/D_DEG 境界、1細胞矩形波軌道 |
