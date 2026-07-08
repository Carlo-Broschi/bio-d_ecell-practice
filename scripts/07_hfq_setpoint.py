"""07. なぜ Hfq を増やすと逆に効かなくなるのか — random-order 結合の帰結
       （Sagawa, Shin, Hussein & Lim, Nucleic Acids Res 2015）。

出典: Sagawa S, Shin J-E, Hussein R, Lim HN (2015) "Paradoxical suppression of
small RNA activity at high Hfq concentrations due to random-order binding."
Nucleic Acids Res 43(17):8502-8515.

06 で見た「Hfq 釣鐘」の正体を、結合の順序という観点で切る。sRNA と mRNA は
どちらが先に Hfq へ乗ってもよい（random order, bi-uni 反応）。この論文の主張は:

    random-order だからこそ、高 Hfq で sRNA 活性がむしろ落ちる（paradoxical suppression）。
    もし compulsory-order（必ず sRNA が先に Hfq に乗る）なら、高 Hfq でも単調に飽和するだけで抑制は起きない。

つまり「Hfq set-point（最適 Hfq）」は結合順序が生む性質。

モデル: 06 と同じ骨格。random と compulsory の違いは1点だけ——
  random    : s+H<=>sH と m+H<=>mH の両方（mRNA も単独で Hfq に乗れる）
  compulsory: s+H<=>sH のみ（mRNA は sH に対してしか結合しない = 必ず sRNA が先）
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

H, S, M, SH, MH, T, D = (Species(x) for x in ["H", "s", "m", "sH", "mH", "T", "D"])
A_S = A_M = 10.0
BETA, KA, KD, K5 = 1.0, 1.0, 1.0, 10.0
SP = ["H", "s", "m", "sH", "mH", "T", "D"]


def R(a, b, k):
    return ReactionRule(a, b, k)


def build(h_tot, mode):
    """mode = 'random'（両順）または 'compulsory'（sRNA が必ず先）。"""
    common = [
        R([], [S], A_S), R([], [M], A_M), R([S], [], BETA), R([M], [], BETA),
        R([S, H], [SH], KA), R([SH], [S, H], KD),          # sRNA が Hfq に乗る
        R([SH, M], [T], KA), R([T], [SH, M], KD),          # sH + mRNA -> 三者
        R([T], [D, H], K5),
        R([SH], [H], BETA), R([T], [H], BETA), R([D], [], BETA),
    ]
    random_extra = [
        R([M, H], [MH], KA), R([MH], [M, H], KD),          # mRNA も単独で Hfq に乗れる
        R([MH, S], [T], KA), R([T], [MH, S], KD),          # mH + sRNA -> 三者
        R([MH], [H], BETA),
    ]
    mdl = NetworkModel()
    for r in common + (random_extra if mode == "random" else []):
        mdl.add_reaction_rule(r)
    return mdl


def steady(h_tot, mode):
    ret = run_simulation(400.0, y0={"H": h_tot}, model=build(h_tot, mode),
                         solver="ode", ndiv=1, species_list=SP)
    return dict(zip(SP, ret.as_array()[-1][1:]))


def pct_duplex(h_tot, mode):
    v = steady(h_tot, mode)
    tot_m = v["m"] + v["mH"] + v["T"] + v["D"]
    return 100 * v["D"] / tot_m if tot_m > 0 else 0.0


def main():
    rel = np.logspace(-2.5, 3.5, 43)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # (左) random（釣鐘・抑制あり） vs compulsory（単調飽和・抑制なし）
    for mode, c, lab in [("random", "C3", "random-order (both can bind Hfq first)"),
                         ("compulsory", "C0", "compulsory-order (sRNA binds first)")]:
        y = np.array([pct_duplex(r * 10, mode) for r in rel])
        ax1.plot(np.log10(rel), y, "-", color=c, lw=2, label=lab)
    ax1.axvline(0, ls=":", c="gray")
    ax1.annotate("Hfq set-point\n(random-order only)", xy=(0, 55), xytext=(-2.2, 44),
                 fontsize=8, color="C3",
                 arrowprops=dict(arrowstyle="->", color="C3", lw=1))
    ax1.set_xlabel("relative Hfq   log10( H_tot / target mRNA )")
    ax1.set_ylabel("% target mRNA in duplex")
    ax1.set_title("binding order decides whether high Hfq suppresses")
    ax1.legend(fontsize=8, loc="upper right")

    # (右) なぜ random だと落ちるか: sRNA の内訳（高 Hfq で sH に囚われ duplex 崩壊）
    parts, plab = ["s", "sH", "T", "D"], ["free s", "s.Hfq (singly)", "ternary", "duplex"]
    cols = ["#b8c4be", "#e0a24b", "#8c78e6", "#0e9e6e"]
    data = []
    for r in [0.3, 1.0, 100.0]:
        v = steady(r * 10, "random")
        tot_s = v["s"] + v["sH"] + v["T"] + v["D"]
        data.append([v[p] / tot_s for p in parts])
    data = np.array(data)
    bottom = np.zeros(3); x = np.arange(3)
    for j, (pl, col) in enumerate(zip(plab, cols)):
        ax2.bar(x, data[:, j], bottom=bottom, color=col, label=pl, width=.6)
        bottom += data[:, j]
    ax2.set_xticks(x)
    ax2.set_xticklabels(["low\n(relHfq=0.3)", "optimal\n(relHfq=1)", "high\n(relHfq=100)"], fontsize=9)
    ax2.set_ylabel("fraction of sRNA")
    ax2.set_title("random-order: high Hfq traps sRNA in singly-bound s.Hfq")
    ax2.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))

    fig.tight_layout()
    out = "outputs/07_hfq_setpoint.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print("saved:", out)
    for r in [0.1, 1, 3, 10, 100]:
        print(f"relHfq={r:6g}  random={pct_duplex(r*10,'random'):5.1f}  "
              f"compulsory={pct_duplex(r*10,'compulsory'):5.1f}")


if __name__ == "__main__":
    main()
