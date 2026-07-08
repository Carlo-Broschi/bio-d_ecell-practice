"""10. sRNA は「転写後」だけでなく「共転写」でも標的を抑える（Reyer et al., Cell Reports 2021）。

出典: Reyer MA, et al. (2021) "Kinetic modeling reveals additional regulation at
co-transcriptional level by post-transcriptional sRNA regulators." Cell Rep 36:109764.

SgrS は「転写後制御 sRNA」の代表とされてきた（成熟 mRNA に結合し翻訳阻害＋共分解）。
だがこの論文は単一細胞イメージング＋速度論モデルで、SgrS が **転写途中の nascent mRNA にも効き**、
実効的な mRNA 生成そのものを下げている（共転写制御）ことを示した。

モデル（本notebookは論文の構造を代表パラメータで実装。比 α_ms/α_m≈0.46 は論文 WT 値）:
  種: s(sRNA, pre-induced で一定), m(成熟mRNA), ms(sRNA-mRNA複合体), p(タンパク質)
    ∅ -> m            (実効転写 α_m。共転写制御ありなら α_m→0.46·α_m)
    m -> ∅            (分解 β_m)
    m + s <=> ms      (結合 k_on / 解離 k_off)
    ms -> ∅           (共分解 β_ms)
    m -> m + p        (翻訳 k_x)
    ms -> ms + p      (結合 mRNA の翻訳 k_xs << k_x)
    p -> ∅            (希釈 β_p, タンパク質半減期は長い)

3条件を比較:  minus（sRNAなし）／ post（転写後のみ）／ co（転写後＋共転写）。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

S, M, MS, P = (Species(x) for x in ["s", "m", "ms", "p"])
# 代表パラメータ（/min）。構造と α_ms/α_m≈0.46 は論文、絶対値は代表値。
A_M, B_M, K_X, B_P = 60.0, 0.2, 2.0, 0.008
A_S, B_S = 15.0, 0.1            # sRNA を ~150 copies に保つ（pre-induced）
K_ON, K_OFF, B_MS, K_XS = 0.002, 0.05, 0.4, 0.2
P_CO = 0.46                     # 共転写制御での実効転写の低下（WT 実測比）
SP = ["s", "m", "ms", "p"]


def R(a, b, k):
    return ReactionRule(a, b, k)


def build(mode):
    a_m_eff = A_M * (P_CO if mode == "co" else 1.0)
    rules = [R([], [M], a_m_eff), R([M], [], B_M), R([M], [M, P], K_X), R([P], [], B_P)]
    if mode != "minus":
        rules += [R([], [S], A_S), R([S], [], B_S),
                  R([M, S], [MS], K_ON), R([MS], [M, S], K_OFF),
                  R([MS], [], B_MS), R([MS], [MS, P], K_XS)]
    mdl = NetworkModel()
    for r in rules:
        mdl.add_reaction_rule(r)
    return mdl


def run(mode, t_end=25.0, ndiv=100):
    y0 = {"s": 150.0} if mode != "minus" else {}
    ret = run_simulation(t_end, y0=y0, model=build(mode), solver="ode",
                         ndiv=ndiv, species_list=SP)
    return ret.as_array()


def main():
    styles = [("minus", "C7", "no sRNA (Δsgr)"),
              ("post", "C0", "post-transcriptional only"),
              ("co", "C3", "post + co-transcriptional")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    folds = {}
    base = run("minus")
    for mode, c, lab in styles:
        a = run(mode)
        t = a[:, 0]
        ax1.plot(t, a[:, 2], "-", color=c, lw=2, label=lab)       # mRNA (col index 2 = m)
        ax2.plot(t, a[:, 4], "-", color=c, lw=2, label=lab)       # protein (col 4 = p)
        folds[mode] = (base[-1][2] / a[-1][2], base[-1][4] / a[-1][4])
    ax1.set_xlabel("time after mRNA induction (min)")
    ax1.set_ylabel("target mRNA (copies)")
    ax1.set_title("target mRNA: co-transcriptional deepens repression")
    ax1.legend(fontsize=9)
    ax2.set_xlabel("time after mRNA induction (min)")
    ax2.set_ylabel("protein (a.u.)")
    ax2.set_title("protein output follows")
    ax2.legend(fontsize=9)

    # 抑制倍率の注記
    txt = ("repression fold (mRNA / protein):\n"
           f"  post-only:    {folds['post'][0]:.1f}x / {folds['post'][1]:.1f}x\n"
           f"  post + co:    {folds['co'][0]:.1f}x / {folds['co'][1]:.1f}x")
    ax1.text(0.03, 0.03, txt, transform=ax1.transAxes, fontsize=8,
             va="bottom", family="monospace",
             bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=.8))

    fig.tight_layout()
    out = "outputs/10_cotranscriptional.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for mode in ["minus", "post", "co"]:
        a = run(mode)
        print(f"{mode:6s}: mRNA(25)={a[-1][2]:6.1f}  protein(25)={a[-1][4]:8.0f}")


if __name__ == "__main__":
    main()
