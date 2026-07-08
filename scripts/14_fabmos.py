"""14. スイッチの「頻度」で表現型のばらつきを操る — FABMOS（Hung et al., Nat Commun 2014）。

出典: Hung M, Chang E, Hussein R, et al. (2014) "Modulating the frequency and bias of
stochastic switching to control phenotypic variation." Nat Commun 5:4574.

シリーズ唯一の確率論（Gillespie）。細胞は遺伝子の OFF/ON 状態を確率的に行き来し、ON のときだけ
タンパク質を作る。FABMOS の発見:

    平均発現を保ったまま「スイッチ頻度」だけを変えると、集団の分布が
    bimodal（低頻度：OFF 群と ON 群にくっきり二峰）と unimodal（高頻度：中央に一峰）を行き来する。

直感: スイッチが遅い(低頻度)と、各細胞は OFF か ON に長く留まり、タンパク質が 0 か最大値に振り切れる → 二峰。
      スイッチが速い(高頻度)と、留まる暇なく状態が入れ替わり、タンパク質は時間平均されて中央に集まる → 一峰。

モデル（遺伝子1コピー: Goff + Gon = 1）:
    Goff -> Gon | k_on,   Gon -> Goff | k_off      （確率的スイッチ）
    Gon  -> Gon + P | s                             （ON でのみ産生）
    P -> ∅ | d                                       （分解）
bias = k_on/k_off を固定（平均を一定に保つ）。頻度 f は k_on,k_off を一律に倍する。
各細胞 = 1本の Gillespie 軌道。多数の細胞の終端 P でヒストグラムを作る。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

Goff, Gon, P = (Species(x) for x in ["Goff", "Gon", "P"])
S_PROD, D_DEG = 50.0, 1.0     # ON 定常 = s/d = 50、タンパク質寿命 = 1
SP = ["Goff", "Gon", "P"]


def R(a, b, k):
    return ReactionRule(a, b, k)


def build(f, bias=1.0):
    k_off = f
    k_on = f * bias               # bias = k_on/k_off（=1 で ON 割合 0.5、平均 25）
    mdl = NetworkModel()
    for x in [R([Goff], [Gon], k_on), R([Gon], [Goff], k_off),
              R([Gon], [Gon, P], S_PROD), R([P], [], D_DEG)]:
        mdl.add_reaction_rule(x)
    return mdl


def population(f, n=300, t_end=30.0):
    vals = np.empty(n)
    for i in range(n):
        y0 = {"Goff": 1} if i % 2 == 0 else {"Gon": 1}   # 初期状態の偏りを回避
        a = run_simulation(t_end, y0=y0, model=build(f), solver="gillespie",
                           ndiv=1, species_list=SP, rndseed=i)
        vals[i] = a.as_array()[-1][3]
    return vals


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # (左) 同じ平均・違う分布: 低頻度=bimodal / 高頻度=unimodal
    bins = np.linspace(0, 70, 30)
    for f, c, lab in [(0.1, "C3", "low frequency"), (10.0, "C0", "high frequency")]:
        v = population(f, n=400)
        ax1.hist(v, bins=bins, alpha=0.55, color=c, density=True,
                 label=f"{lab}  (mean={v.mean():.0f}, CV={v.std()/v.mean():.2f})")
    ax1.set_xlabel("protein per cell")
    ax1.set_ylabel("fraction of cells")
    ax1.set_title("same mean, different shape:\nlow freq = bimodal, high freq = unimodal")
    ax1.legend(fontsize=8.5)

    # (右) 頻度を掃引: CV(ばらつき) が bimodal->unimodal で下がる。平均は一定。
    f_grid = np.logspace(-1.3, 1.5, 11)
    means, cvs = [], []
    for f in f_grid:
        v = population(f, n=250)
        means.append(v.mean()); cvs.append(v.std() / v.mean())
    ax2.plot(np.log10(f_grid), cvs, "o-", color="C4", lw=2, label="CV (spread)")
    ax2.set_xlabel("switching frequency  log10(f)")
    ax2.set_ylabel("CV of protein", color="C4")
    ax2.tick_params(axis="y", labelcolor="C4")
    ax2.set_title("frequency is the knob: CV falls as switching speeds up")
    ax2b = ax2.twinx()
    ax2b.plot(np.log10(f_grid), means, "s--", color="gray", lw=1.5, label="mean")
    ax2b.set_ylabel("mean protein", color="gray")
    ax2b.set_ylim(0, 50)
    ax2b.tick_params(axis="y", labelcolor="gray")
    ax2.text(0.03, 0.5, "bimodal", transform=ax2.transAxes, color="C3", fontsize=9)
    ax2.text(0.7, 0.15, "unimodal", transform=ax2.transAxes, color="C0", fontsize=9)

    fig.tight_layout()
    out = "outputs/14_fabmos.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for f in [0.1, 1.0, 10.0]:
        v = population(f, n=250)
        print(f"f={f:5g}: mean={v.mean():5.1f}  CV={v.std()/v.mean():.2f}")


if __name__ == "__main__":
    main()
