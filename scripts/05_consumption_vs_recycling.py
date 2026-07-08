"""05. 閾値の正体は「sRNA の消費」— 化学量論 vs 触媒を1つの分岐確率でつなぐ。

04 では化学量論(m+s→∅)を threshold、触媒(m+s→s)を graded と対比した。
だがそれは触媒側にわざと弱い速度 (k_cat=0.1) を与えた見かけだった。
ここでは両機構を **同じペアリング率 k_pair** の下で、分岐確率 φ で連続的につなぐ:

    m + s -> ∅   | k_pair * (1 - φ)    (消費: sRNA も死ぬ)
    m + s -> s   | k_pair * φ          (再利用: sRNA は触媒として残る)

mRNA は φ によらず常に k_pair*[m][s] で殺される。φ は「sRNA が消費されるか」だけを変える。

得られる正直な描像:
  - φ=0 (完全消費): sRNA は mRNA と 1:1 で滴定され、s は a_s < a_m の間ほぼ 0 に
    抑えられる → a_s = a_m に固定した閾値。
  - φ→1 (触媒): s = a_s/b_s と無制限に蓄積し、少ない a_s でも持続的に mRNA を壊す
    → silencing はむしろ強まり、閾値は低 a_s 側へ動く。「触媒=なだらか」ではない。

結論: 閾値スイッチを生むのは触媒/消費の別ではなく **sRNA が消費（滴定）されること**。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

M, S = Species("m"), Species("s")
A_M, B_M, B_S, K_PAIR = 10.0, 1.0, 1.0, 5.0


def build(a_s, phi):
    mdl = NetworkModel()
    rules = [
        ReactionRule([], [M], A_M),
        ReactionRule([], [S], a_s),
        ReactionRule([M], [], B_M),
        ReactionRule([S], [], B_S),
    ]
    if (1 - phi) > 0:
        rules.append(ReactionRule([M, S], [], K_PAIR * (1 - phi)))    # 消費
    if phi > 0:
        rules.append(ReactionRule([M, S], [S], K_PAIR * phi))         # 再利用
    for rr in rules:
        mdl.add_reaction_rule(rr)
    return mdl


def steady(a_s, phi):
    ret = run_simulation(200.0, y0={"m": 0, "s": 0}, model=build(a_s, phi),
                         solver="ode", ndiv=2000, species_list=["m", "s"])
    m, s = ret.as_array()[-1][1:3]
    return m, s


def main():
    a_grid = np.linspace(0, 20, 41)
    phis = [0.0, 0.25, 0.5, 0.75, 1.0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # (左) φ を上げるほど silencing は強く・閾値は低 a_s 側へ
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(phis)))
    for phi, c in zip(phis, cmap):
        m = np.array([steady(a, phi)[0] for a in a_grid])
        label = f"φ={phi:.2f}" + (" (consume)" if phi == 0 else " (recycle)" if phi == 1 else "")
        ax1.plot(a_grid, m / A_M, "-", color=c, label=label)
    ax1.axvline(A_M, ls=":", c="gray")
    ax1.set_xlabel("sRNA transcription rate  a_s")
    ax1.set_ylabel("relative mRNA  <m>/<m>(0)")
    ax1.set_title("more recycling (φ↑) → stronger, lower-threshold silencing")
    ax1.legend()

    # (右) なぜか: 消費は s を閾値まで抑え込み、再利用は s を蓄積させる
    for phi, c, ls in [(0.0, "C0", "-"), (1.0, "C3", "-")]:
        s = np.array([steady(a, phi)[1] for a in a_grid])
        tag = "φ=0 consume" if phi == 0 else "φ=1 recycle"
        ax2.plot(a_grid, s, ls, color=c, label=f"sRNA level, {tag}")
    ax2.plot(a_grid, a_grid / B_S, ":", color="gray", label="a_s/b_s (uncapped)")
    ax2.axvline(A_M, ls=":", c="gray")
    ax2.set_xlabel("sRNA transcription rate  a_s")
    ax2.set_ylabel("steady-state sRNA  <s>")
    ax2.set_title("consumption titrates sRNA away below the threshold")
    ax2.legend()

    fig.tight_layout()
    out = "outputs/05_consumption_vs_recycling.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for a in [0, 5, 10, 15, 20]:
        row = [steady(float(a), p)[0] / A_M for p in phis]
        print(f"a_s={a:4.0f}  m/m0 by φ{phis} = " + " ".join(f"{x:4.2f}" for x in row))


if __name__ == "__main__":
    main()
