"""Hfq 依存 sRNA による標的 mRNA 抑制の最小回路（E-Cell4）。

細菌の転写後制御: sRNA が Hfq を介して標的 mRNA と塩基対合し、
両者が「共分解（co-degradation）」で一緒に消える。最小モデルは 5 反応:

    ∅ -> m        (a_m : mRNA 転写)
    ∅ -> s        (a_s : sRNA 転写)     <- これを掃引する制御ノブ
    m -> ∅        (b_m : mRNA 分解)
    s -> ∅        (b_s : sRNA 分解)
    m + s -> ∅     (k   : Hfq 依存ペアリングによる共分解)

決定論(ODE)では、a_s が a_m を超えると mRNA が閾値的に急落する
（Levine et al. 2007 の "threshold-linear" 応答）。
確率論(Gillespie)では、その閾値近傍で mRNA のゆらぎ（CV）が最大になる
= sRNA 制御は閾値付近で最もノイジー、という生物学的に重要な帰結が見える。

bio-a (Hfq 系統解析) と地続きの題材: 進化が保存してきた Hfq/sRNA 制御が
「回路としてどう振る舞うか」を数式で回す。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

M, S = Species("m"), Species("s")
A_M, B_M, B_S, K = 10.0, 1.0, 1.0, 100.0  # 転写/mRNA分解/sRNA分解/ペアリング


def build_model(a_s):
    mdl = NetworkModel()
    rules = [
        ReactionRule([], [M], A_M),
        ReactionRule([], [S], a_s),
        ReactionRule([M], [], B_M),
        ReactionRule([S], [], B_S),
        ReactionRule([M, S], [], K),
    ]
    for rr in rules:
        mdl.add_reaction_rule(rr)
    return mdl


def ode_steady(a_s):
    ret = run_simulation(80.0, y0={"m": 0, "s": 0}, model=build_model(a_s),
                         solver="ode", ndiv=800, species_list=["m", "s"])
    return ret.as_array()[-1][1]  # 定常 mRNA


def gillespie_stats(a_s, t_end=400.0, ndiv=4000, burn=0.25, seed=0):
    """1 本の長い軌道の後半から mRNA の平均・標準偏差を得る。"""
    ret = run_simulation(t_end, y0={"m": 0, "s": 0}, model=build_model(a_s),
                         solver="gillespie", ndiv=ndiv,
                         species_list=["m", "s"], rndseed=seed)
    arr = ret.as_array()
    tail = arr[int(len(arr) * burn):, 1]  # burn-in 後の mRNA
    return tail.mean(), tail.std()


def main():
    a_s_grid = np.array([0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 20], dtype=float)

    ode = np.array([ode_steady(a) for a in a_s_grid])
    g_mean, g_std = np.array([gillespie_stats(a) for a in a_s_grid]).T

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # (左) 応答曲線: ODE 閾値 + Gillespie 平均±SD
    ax1.plot(a_s_grid, ode, "-", color="C0", label="ODE (deterministic)")
    ax1.errorbar(a_s_grid, g_mean, yerr=g_std, fmt="o", color="C1", capsize=3,
                 label="Gillespie mean ± SD")
    ax1.axvline(A_M, ls=":", c="gray")
    ax1.text(A_M + 0.2, ax1.get_ylim()[1] * 0.88, "a_s = a_m\n(threshold)", fontsize=9)
    ax1.set_xlabel("sRNA transcription rate  a_s")
    ax1.set_ylabel("steady-state mRNA  <m>")
    ax1.set_title("threshold-linear silencing")
    ax1.legend()

    # (右) ノイズ増幅: Fano factor = var/mean。Poisson なら 1、超過分が増幅。
    # 平均が小さすぎる領域は推定が不安定なので描かない (nan)。
    fano = np.where(g_mean > 0.3, g_std ** 2 / np.where(g_mean > 0, g_mean, 1), np.nan)
    ax2.axhline(1.0, ls="--", c="gray", lw=1, label="Poisson (Fano = 1)")
    ax2.plot(a_s_grid, fano, "o-", color="C3")
    ax2.axvline(A_M, ls=":", c="gray")
    ax2.set_xlabel("sRNA transcription rate  a_s")
    ax2.set_ylabel("noise  Fano = var/mean  of mRNA")
    ax2.set_title("noise amplified near / above threshold")
    ax2.legend()

    fig.tight_layout()
    out = "outputs/02_srna_silencing.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for a, o, gm, gs in zip(a_s_grid, ode, g_mean, g_std):
        print(f"a_s={a:4.0f}  ODE={o:6.2f}  Gillespie={gm:6.2f}±{gs:4.2f}")


if __name__ == "__main__":
    main()
